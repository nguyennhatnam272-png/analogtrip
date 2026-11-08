#!/usr/bin/env python3
"""Build Earth's elevation grid (DEM) for the hover readout and the 3D relief.

Source: "Terrain Tiles" on AWS Open Data (s3://elevation-tiles-prod, Terrarium PNG
encoding, zoom 5, about 4.9 km per pixel at the equator). It merges SRTM, GMTED2010,
ETOPO1 bathymetry and other national datasets. Attribution is required; see README.
The tiles stop at 85.05 deg N/S (Web Mercator), so the polar caps come from the
NOAA ETOPO1 10 arc-minute grid already in source/.

Outputs (in earth/):
  earth-dem-0.1.b64.txt  Global grid, 0.1 deg (3600 x 1800 cells), rows north to south, cell centres at
                     lon -180+0.05+0.1i, lat 90-0.05-0.1j. Values in 5 m units (deep ocean rounded
                     to 25 m), MED-predicted, low bytes then high bytes, zlib-deflated, then base64
                     (static hosts serve .txt everywhere; they usually gzip it in transit).
  dem/{i}_{j}.b64.txt    0.05 deg detail, 30 x 30 deg bundles (i = 0..11 from 180W, j = 0..5 from 90N),
                     602 x 602 cells including a 1-cell border on every side, same encoding.
                     True metres, no exaggeration; the page samples with bilinear interpolation.
  earth-height.png  2048 x 1024 display texture (R*256+G = metres + 11000), area-averaged
                    from the same grid, used only to draw the relief.

Run:  python3 build_earth_dem.py        (needs numpy, Pillow, netCDF4, scipy)
      Download the tiles first:  python3 build_earth_dem.py --download
"""
import base64, pathlib, sys, zlib, urllib.request, concurrent.futures as cf
import numpy as np, netCDF4
from PIL import Image

HERE = pathlib.Path(__file__).resolve().parent
TILES = HERE / 'source' / 'terrarium-z5'
Z = 5
RES = 0.05


def download():
    TILES.mkdir(parents=True, exist_ok=True)
    def get(xy):
        x, y = xy
        f = TILES / f'{x}_{y}.png'
        if not f.exists():
            f.write_bytes(urllib.request.urlopen(f'https://elevation-tiles-prod.s3.amazonaws.com/terrarium/{Z}/{x}/{y}.png', timeout=60).read())
    with cf.ThreadPoolExecutor(16) as ex:
        list(ex.map(get, [(x, y) for x in range(2 ** Z) for y in range(2 ** Z)]))


def mosaic():
    n = 2 ** Z
    m = np.zeros((n * 256, n * 256), np.float32)
    for x in range(n):
        for y in range(n):
            a = np.asarray(Image.open(TILES / f'{x}_{y}.png').convert('RGB'), dtype=np.float32)
            m[y * 256:(y + 1) * 256, x * 256:(x + 1) * 256] = a[..., 0] * 256 + a[..., 1] + a[..., 2] / 256 - 32768
    return m


def bilinear(img, fx, fy, wrap_x=True):
    h, w = img.shape
    x0 = np.floor(fx).astype(np.int64); y0 = np.clip(np.floor(fy).astype(np.int64), 0, h - 2)
    tx = fx - x0; ty = np.clip(fy - y0, 0, 1)
    x1 = x0 + 1
    if wrap_x:
        x0 %= w; x1 %= w
    a = img[y0, x0] * (1 - tx) + img[y0, x1] * tx
    b = img[y0 + 1, x0] * (1 - tx) + img[y0 + 1, x1] * tx
    return a * (1 - ty) + b * ty


def etopo_polar(lat, lon):
    ds = netCDF4.Dataset(HERE / 'source' / 'earth-topography-10arcmin.nc')
    z = np.asarray(ds.variables['topography'][:], dtype=np.float32)[:, :-1]
    LAT, LON = np.meshgrid(lat, lon, indexing='ij')
    return bilinear(z, (LON + 180) * 6, (LAT + 90) * 6)


def build():
    m = mosaic()
    N = m.shape[0]
    W, H = int(360 / RES), int(180 / RES)
    lat = 90 - (np.arange(H) + 0.5) * RES
    lon = -180 + (np.arange(W) + 0.5) * RES
    dem = np.zeros((H, W), np.float32)
    fx = (lon + 180) / 360 * N - 0.5
    for j, la in enumerate(lat):
        if abs(la) < 85.0:
            fy = (1 - np.log(np.tan(np.pi / 4 + np.radians(la) / 2)) / np.pi) / 2 * N - 0.5
            dem[j] = bilinear(m, fx, np.full_like(fx, fy))
    # The tiles hold bedrock under the Greenland and Antarctic ice sheets; use the ETOPO1 ice surface there
    ice_rows = (lat < -59) | (lat > 58)
    ice = etopo_polar(lat[ice_rows], lon)
    LATI, LONI = np.meshgrid(lat[ice_rows], lon, indexing='ij')
    sheet = (LATI < -59) | ((LATI > 58) & (LONI > -75) & (LONI < -10))
    d = dem[ice_rows]
    use = sheet & (ice - d > 30)
    d[use] = ice[use]
    dem[ice_rows] = d
    polar = np.abs(lat) >= 84.0
    pol = etopo_polar(lat[polar], lon)
    t = np.clip((np.abs(lat[polar]) - 84.0) / 1.0, 0, 1)[:, None]          # blend over 84-85 deg
    dem[polar] = dem[polar] * (1 - t) + pol * t
    # ETOPO1's 10-arcmin grid comes from a spherical-harmonic model that degrades near the poles
    # (2,742 m at the South Pole against a surveyed 2,835 m). Inside 85 deg S, over land, use ERA5's
    # orography instead (2,789 m at the pole); it is also the terrain the pressure data sits on.
    era = HERE / 'source' / 'era5-climatology.npz'
    if era.exists():
        d = np.load(era)
        ez, el = np.roll(d['z'], 720, axis=1), np.roll(d['lsm'], 720, axis=1)   # node grid, lon -180..179.75
        rows = np.where(lat < -84.5)[0]
        LATP, LONP = np.meshgrid(lat[rows], lon, indexing='ij')
        fx, fy = (LONP + 180) / 0.25, (90 - LATP) / 0.25
        zz, ll = bilinear(ez, fx, fy), bilinear(el, fx, fy)
        t = np.clip((-LATP - 84.5) / 1.0, 0, 1) * (ll > 0.5)
        dem[rows] = dem[rows] * (1 - t) + zz * t
    return lat, lon, dem


def encode(g):
    """Quantise (5 m on land, 25 m below -200 m), MED-predict, split bytes, deflate (zlib)."""
    q = np.where(g < -200, np.round(g / 25) * 25, np.round(g / 5) * 5)
    q = np.clip(np.round(q / 5), -2200, 1800).astype(np.int32)          # stored in 5 m units
    a = np.zeros_like(q); b = np.zeros_like(q); c = np.zeros_like(q)
    a[:, 1:] = q[:, :-1]; b[1:, :] = q[:-1, :]; c[1:, 1:] = q[:-1, :-1]
    mx, mn = np.maximum(a, b), np.minimum(a, b)
    r = (q - np.where(c >= mx, mn, np.where(c <= mn, mx, a + b - c))).astype('<i2')
    split = r.view(np.uint8).reshape(-1, 2).T.copy()                      # all low bytes, then all high bytes
    return zlib.compress(split.tobytes(), 9)


def write(dem):
    out = HERE / 'earth'
    (out / 'dem').mkdir(parents=True, exist_ok=True)
    H, W = dem.shape                                                        # 3600 x 7200 at 0.05 deg
    # 1. global 0.1 deg grid (3600 x 1800), loaded the first time Terrain or Air pressure is opened
    g10 = dem.reshape(H // 2, 2, W // 2, 2).mean(axis=(1, 3))
    (out / 'earth-dem-0.1.b64.txt').write_bytes(base64.b64encode(encode(g10)))
    # 2. 0.05 deg detail in 30 x 30 deg bundles with a 1-cell border (602 x 602), loaded on hover
    total = 0
    for bj in range(6):
        for bi in range(12):
            rows = np.clip(np.arange(bj * 600 - 1, bj * 600 + 601), 0, H - 1)
            cols = np.arange(bi * 600 - 1, bi * 600 + 601) % W
            data = encode(dem[np.ix_(rows, cols)])
            (out / 'dem' / f'{bi}_{bj}.b64.txt').write_bytes(base64.b64encode(data))
            total += len(data)
    # 3. display texture for the relief, area-averaged to 2048 x 1024
    from scipy.ndimage import zoom
    small = zoom(dem, (1024 / H, 2048 / W), order=1, grid_mode=True, mode='grid-wrap')
    u16 = np.clip(np.round(small + 11000), 0, 65535).astype(np.uint16)
    rgb = np.zeros(small.shape + (3,), np.uint8); rgb[..., 0] = u16 >> 8; rgb[..., 1] = u16 & 255
    Image.fromarray(rgb, 'RGB').save(out / 'earth-height.png', optimize=True)
    print('earth-dem-0.1.b64.txt', round((out / 'earth-dem-0.1.b64.txt').stat().st_size / 1e6, 2), 'MB;',
          '72 detail bundles', round(total / 1e6, 2), 'MB;',
          'earth-height.png', round((out / 'earth-height.png').stat().st_size / 1e6, 2), 'MB')


if __name__ == '__main__':
    if '--download' in sys.argv:
        download()
    lat, lon, dem = build()
    np.save(HERE / 'source' / 'earth-dem-005.npy', dem)
    write(dem)
