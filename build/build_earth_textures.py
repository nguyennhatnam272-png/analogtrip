#!/usr/bin/env python3
"""Build the Earth textures used by the Analog Atlas 3D globe.

Inputs (in ./source, see README.md for where each comes from):
  earth-topography-10arcmin.nc  NOAA ETOPO1 ice-surface elevation, 10 arc-minute grid (public domain)
  earth-blue-marble.jpg         NASA Blue Marble surface colour, 4096 x 2048 (public domain)
  earth-water.png               water mask that matches the Blue Marble image, 1600 x 800

Outputs (in ./earth, published next to the page):
  earth-height.png   2048 x 1024, height encoded in 16 bits: R*256 + G = metres + 11000
  earth-water.png    2048 x 1024, 8-bit water coverage (255 = water), smooth coastlines
  earth-color.jpg    4096 x 2048, Blue Marble land with ocean re-coloured from real bathymetry
  earth-clouds.jpg   2048 x 1024, cloud opacity, procedural and shaped by the page's pressure model

The page derives normals, displacement, the Terrain view and the Air pressure view
from earth-height.png in the browser, using one exaggeration constant
(EARTH_RELIEF_EXAGGERATION in the page script).

Run:  python3 build_earth_textures.py      (needs numpy, Pillow, netCDF4)
"""
import pathlib
import numpy as np
import netCDF4
from PIL import Image

HERE = pathlib.Path(__file__).resolve().parent
SRC = HERE / 'source'
OUT = HERE / 'earth'
W, H = 2048, 1024          # height, water and cloud rasters (the source grid is 2160 x 1080)
CW, CH = 4096, 2048        # colour texture, matches the Blue Marble image
HEIGHT_OFFSET = 11000      # metres added before 16-bit encoding


def cell_centres(w, h):
    lat = 90 - (np.arange(h) + 0.5) * 180 / h
    lon = -180 + (np.arange(w) + 0.5) * 360 / w
    return lat, lon


def resample_dem(w, h):
    """Bilinear resample of the grid-registered ETOPO grid to cell centres, wrapping at 180 degrees."""
    ds = netCDF4.Dataset(SRC / 'earth-topography-10arcmin.nc')
    z = np.asarray(ds.variables['topography'][:], dtype=np.float32)   # (1081, 2161), lat -90..90, lon -180..180
    glat = np.asarray(ds.variables['latitude'][:])
    glon = np.asarray(ds.variables['longitude'][:])
    z = z[:, :-1]                                # drop the duplicate +180 column so columns wrap cleanly
    step = glon[1] - glon[0]
    lat, lon = cell_centres(w, h)
    fy = (lat - glat[0]) / step
    fx = (lon - glon[0]) / step
    y0 = np.clip(np.floor(fy).astype(int), 0, z.shape[0] - 2)
    ty = (fy - y0)[:, None]
    x0 = np.floor(fx).astype(int)
    tx = (fx - x0)[None, :]
    nx = z.shape[1]
    x0m, x1m = x0 % nx, (x0 + 1) % nx
    r0 = z[y0][:, x0m] * (1 - tx) + z[y0][:, x1m] * tx
    r1 = z[y0 + 1][:, x0m] * (1 - tx) + z[y0 + 1][:, x1m] * tx
    return (r0 * (1 - ty) + r1 * ty).astype(np.float32)


def smoothstep(e0, e1, x):
    t = np.clip((x - e0) / (e1 - e0), 0, 1)
    return t * t * (3 - 2 * t)


# ---------- noise for clouds ----------
_rng = np.random.default_rng(7)
_perm = np.concatenate([_rng.permutation(256)] * 2).astype(np.int32)
_grads = np.array([[1, 1, 0], [-1, 1, 0], [1, -1, 0], [-1, -1, 0], [1, 0, 1], [-1, 0, 1], [1, 0, -1], [-1, 0, -1],
                   [0, 1, 1], [0, -1, 1], [0, 1, -1], [0, -1, -1]], dtype=np.float32)


def perlin(x, y, z):
    xi, yi, zi = np.floor(x).astype(np.int32), np.floor(y).astype(np.int32), np.floor(z).astype(np.int32)
    xf, yf, zf = x - xi, y - yi, z - zi
    xi, yi, zi = xi & 255, yi & 255, zi & 255
    u, v, w = [t * t * t * (t * (t * 6 - 15) + 10) for t in (xf, yf, zf)]

    def g(ix, iy, iz, dx, dy, dz):
        gv = _grads[_perm[_perm[_perm[ix] + iy] + iz] % 12]
        return gv[..., 0] * dx + gv[..., 1] * dy + gv[..., 2] * dz

    x1 = xi + 1
    y1 = yi + 1
    z1 = zi + 1
    n000 = g(xi, yi, zi, xf, yf, zf); n100 = g(x1, yi, zi, xf - 1, yf, zf)
    n010 = g(xi, y1, zi, xf, yf - 1, zf); n110 = g(x1, y1, zi, xf - 1, yf - 1, zf)
    n001 = g(xi, yi, z1, xf, yf, zf - 1); n101 = g(x1, yi, z1, xf - 1, yf, zf - 1)
    n011 = g(xi, y1, z1, xf, yf - 1, zf - 1); n111 = g(x1, y1, z1, xf - 1, yf - 1, zf - 1)
    a = n000 + u * (n100 - n000); b = n010 + u * (n110 - n010)
    c = n001 + u * (n101 - n001); d = n011 + u * (n111 - n011)
    e = a + v * (b - a); f = c + v * (d - c)
    return e + w * (f - e)


def fbm(p, octaves, lac=2.0, gain=0.5):
    total, amp, norm = 0, 0.5, 0
    x, y, z = p
    for _ in range(octaves):
        total = total + amp * perlin(x, y, z)
        norm += amp
        x, y, z = x * lac, y * lac, z * lac
        amp *= gain
    return total / norm


# Same typical highs and lows as the page's pressure model (lat, lon, radius deg, anomaly hPa)
SLP_CENTERS = [(36, -30, 12, 6), (34, -145, 14, 6), (-28, -12, 12, 4), (-32, 75, 14, 4), (-33, -100, 14, 4),
               (50, 100, 15, 8), (62, -28, 10, -8), (52, -175, 10, -8)]


def slp(lat, lon):
    alat = np.abs(lat)
    p = (1012 + 6 * np.exp(-((alat - 32) / 10) ** 2) - 12 * np.exp(-((lat + 62) / 7) ** 2)
         - 4 * np.exp(-((lat - 60) / 9) ** 2) - 3 * np.exp(-(lat / 8) ** 2))
    c = np.cos(np.radians(lat))
    for la, lo, r, a in SLP_CENTERS:
        dn = (((lon - lo) + 540) % 360 - 180) * c
        p = p + a * np.exp(-((lat - la) ** 2 + dn ** 2) / r ** 2)
    return p


def build_clouds(lat, lon):
    la, lo = np.radians(lat)[:, None], np.radians(lon)[None, :]
    p = np.stack([np.cos(la) * np.cos(lo), np.sin(la) * np.ones_like(lo), np.cos(la) * np.sin(lo)]).astype(np.float32)
    # Spiral storms at the usual low-pressure centres: twist space around each centre
    for clat, clon, rad, turns in [(60, -30, 11, 0.42), (52, -170, 11, 0.4), (-55, 20, 9, -0.38), (-58, 130, 9, -0.36),
                                   (-52, -100, 9, -0.32), (45, 160, 8, 0.3)]:
        c = np.array([np.cos(np.radians(clat)) * np.cos(np.radians(clon)), np.sin(np.radians(clat)),
                      np.cos(np.radians(clat)) * np.sin(np.radians(clon))], dtype=np.float32)[:, None, None]
        d = np.degrees(np.arccos(np.clip((p * c).sum(0), -1, 1)))
        th = turns * 2 * np.pi * np.exp(-(d / rad) ** 2)
        cosT, sinT = np.cos(th), np.sin(th)
        cxp = np.cross(c, p, axis=0)
        cdp = (c * p).sum(0)
        p = p * cosT + cxp * sinT + c * cdp * (1 - cosT)
    # Domain warp so cloud bands bend and fold
    q = np.stack([fbm(p * 1.7 + o, 4) for o in (np.float32(1.3), np.float32(7.1), np.float32(-4.2))])
    pw = p + 0.38 * q
    # Stretch along east-west so mid-latitude clouds form streaks along fronts
    stretch = 1 + 0.9 * np.exp(-((np.abs(lat)[:, None] - 50) / 15) ** 2)
    big = fbm((pw[0] * 2.3, pw[1] * 2.3 * stretch, pw[2] * 2.3), 6)
    detail = fbm(pw * 9.0, 3)
    field = big + 0.22 * detail
    alat = np.abs(lat)[:, None]
    lat2 = lat[:, None]
    lon2 = lon[None, :]
    # Cloudier under low pressure, clearer under the subtropical highs
    cover = (0.33 + 0.18 * np.exp(-(lat2 / 6) ** 2) + 0.20 * np.exp(-((alat - 58) / 12) ** 2)
             - 0.22 * np.exp(-((alat - 22) / 9) ** 2) + (1013 - slp(lat2, lon2)) / 14 * 0.10)
    th = (0.5 - cover) * 0.9
    alpha = smoothstep(th - 0.04, th + 0.17, field) * (0.78 + 0.22 * smoothstep(-0.3, 0.3, detail))
    alpha = alpha * 0.96
    alpha *= 1 - 0.6 * smoothstep(70, 85, -lat2)           # interior Antarctica is mostly clear
    return np.clip(alpha * 255, 0, 255).astype(np.uint8)


def main():
    OUT.mkdir(exist_ok=True)
    lat, lon = cell_centres(W, H)

    # 1. Height
    h = resample_dem(W, H)
    u16 = np.clip(np.round(h + HEIGHT_OFFSET), 0, 65535).astype(np.uint16)
    rgb = np.zeros((H, W, 3), np.uint8)
    rgb[..., 0] = u16 >> 8
    rgb[..., 1] = u16 & 255
    Image.fromarray(rgb, 'RGB').save(OUT / 'earth-height.png', optimize=True)
    print('height', h.min(), h.max())

    # 2. Water coverage, smooth at the coast (matches the Blue Marble coastline)
    water_src = Image.open(SRC / 'earth-water.png').convert('L')
    water = water_src.resize((W, H), Image.BILINEAR)
    water.save(OUT / 'earth-water.png', optimize=True)

    # 3. Colour: Blue Marble land, ocean shaded by real depth (shelves lighter, abyss darker)
    bm = np.asarray(Image.open(SRC / 'earth-blue-marble.jpg').convert('RGB'), dtype=np.float32)
    wc = np.asarray(water_src.resize((CW, CH), Image.BILINEAR), dtype=np.float32) / 255
    hc = np.asarray(Image.fromarray(h, 'F').resize((CW, CH), Image.BILINEAR), dtype=np.float32)
    depth = np.minimum(hc, -5)
    stops = np.array([0, -150, -400, -1500, -3500, -6000, -9000], np.float32)
    cols = np.array([[34, 92, 124], [28, 80, 114], [20, 58, 96], [13, 38, 78], [9, 27, 63], [6, 19, 50], [4, 14, 40]], np.float32)
    idx = np.clip(np.searchsorted(-stops, -depth) - 1, 0, len(stops) - 2)
    t = np.clip((depth - stops[idx]) / (stops[idx + 1] - stops[idx]), 0, 1)[..., None]
    ramp = cols[idx] * (1 - t) + cols[idx + 1] * t
    lum = bm.mean(axis=2)
    oceanness = 1 - smoothstep(120, 175, lum)          # keep sea ice and ice shelves from the image
    a = (wc * oceanness)[..., None]
    out = bm * (1 - a) + (ramp * 0.78 + bm * 0.22) * a
    Image.fromarray(np.clip(out, 0, 255).astype(np.uint8), 'RGB').save(OUT / 'earth-color.jpg', quality=86, optimize=True, progressive=True)

    # 4. Clouds
    clouds = build_clouds(lat, lon)
    Image.fromarray(clouds, 'L').save(OUT / 'earth-clouds.jpg', quality=82, optimize=True)

    for f in sorted(OUT.iterdir()):
        print(f.name, f.stat().st_size // 1024, 'KB')


if __name__ == '__main__':
    main()
