#!/usr/bin/env python3
"""Build the Moon and Mars textures used by the Analogtrix 3D globes.

Inputs (in ./source, see README.md):
  moon-topography-10arcmin.nc  LOLA-based lunar topography, 10 arc-minute grid (Wieczorek 2015, CC BY)
  moon-albedo-1024.jpg         Moon colour map from the three.js examples (MIT licence)
  mars-elevation-1deg.csv      MOLA-based Mars elevation, 1-degree grid (rows north to south, columns -180..180)
  mars-viking-color.jpg        USGS Viking colour mosaic, 1024 x 512 (public domain)

Outputs:
  moon/moon-height.png  2048 x 1024   R*256 + G = metres + 11000
  moon/moon-color.jpg   2048 x 1024
  mars/mars-height.png  1536 x 768    same encoding; MOLA shape plus procedural small craters (the 1-degree grid has none)
  mars/mars-color.jpg   2048 x 1024

Run:  python3 build_moon_mars_textures.py      (needs numpy, Pillow, netCDF4)
"""
import pathlib
import numpy as np
import netCDF4
from PIL import Image, ImageFilter

HERE = pathlib.Path(__file__).resolve().parent
SRC = HERE / 'source'
OFFSET = 11000


def save_height(h, path):
    u16 = np.clip(np.round(h + OFFSET), 0, 65535).astype(np.uint16)
    rgb = np.zeros(h.shape + (3,), np.uint8)
    rgb[..., 0] = u16 >> 8
    rgb[..., 1] = u16 & 255
    path.parent.mkdir(exist_ok=True)
    Image.fromarray(rgb, 'RGB').save(path, optimize=True)


def upscale_colour(src, w, h, path):
    im = Image.open(src).convert('RGB').resize((w, h), Image.BICUBIC)
    im = im.filter(ImageFilter.UnsharpMask(radius=2, percent=60, threshold=2))
    path.parent.mkdir(exist_ok=True)
    im.save(path, quality=86, optimize=True, progressive=True)


def moon():
    W, H = 2048, 1024
    ds = netCDF4.Dataset(SRC / 'moon-topography-10arcmin.nc')
    z = np.asarray(ds.variables['topography'][:], dtype=np.float32)   # lat 90..-90, lon 0..360
    glat = np.asarray(ds.variables['latitude'][:]); glon = np.asarray(ds.variables['longitude'][:])
    z = z[:, :-1]
    step = abs(glat[1] - glat[0])
    lat = 90 - (np.arange(H) + 0.5) * 180 / H
    lon = -180 + (np.arange(W) + 0.5) * 360 / W
    fy = (glat[0] - lat) / step
    fx = (lon % 360 - glon[0]) / step
    y0 = np.clip(np.floor(fy).astype(int), 0, z.shape[0] - 2); ty = (fy - y0)[:, None]
    x0 = np.floor(fx).astype(int); tx = (fx - x0)[None, :]
    nx = z.shape[1]; a, b = x0 % nx, (x0 + 1) % nx
    r0 = z[y0][:, a] * (1 - tx) + z[y0][:, b] * tx
    r1 = z[y0 + 1][:, a] * (1 - tx) + z[y0 + 1][:, b] * tx
    h = r0 * (1 - ty) + r1 * ty
    save_height(h, HERE / 'moon' / 'moon-height.png')
    moon_albedo(HERE / 'moon' / 'moon-color.jpg')
    print('moon height', h.min(), h.max())


def moon_albedo(path, W=2048, H=1024):
    """Albedo only: keep the large-scale brightness of the source map (maria, highlands) and blur away
    its baked-in crater shadows, then add the bright ray systems of young craters at their real positions.
    Relief shading comes from the LOLA heights in the page's lighting, not from this map."""
    from scipy.ndimage import gaussian_filter, median_filter
    t0 = np.asarray(Image.open(SRC / 'moon-albedo-1024.jpg').convert('L'), dtype=np.float32) / 255
    t0 = median_filter(t0, size=2, mode='nearest')
    low = gaussian_filter(t0, 4.0, mode=('nearest', 'wrap'))
    # keep a third of the fine detail: bright young craters and ray streaks stay, baked-in crater shadows fade
    mixed = low + 0.35 * (t0 - low)
    a = np.asarray(Image.fromarray(mixed, 'F').resize((W, H), Image.BICUBIC), dtype=np.float32)
    lo, hi = np.percentile(a, 2), np.percentile(a, 99.5)
    a = np.clip((a - lo) / (hi - lo), 0, 1)
    a = 0.32 + 0.6 * a                                   # maria about 0.6 of highland brightness, like the real Moon
    alat = np.abs(90 - (np.arange(H)[:, None] + 0.5) * 180 / H)
    t = np.clip((alat - 76) / 12, 0, 1); t = t * t * (3 - 2 * t)
    a = a * (1 - t) + np.median(a) * t                  # the source map smears at the poles; fade to typical highland brightness
    lat = np.radians(90 - (np.arange(H) + 0.5) * 180 / H)[:, None]
    lon = np.radians(-180 + (np.arange(W) + 0.5) * 360 / W)[None, :]
    rng = np.random.default_rng(5)
    # name, lat, lon, crater radius km, ray length km, number of rays, strength
    for _, clat, clon, rkm, length, n, k in [('Tycho', -43.3, -11.2, 43, 650, 70, 0.14), ('Copernicus', 9.6, -20.1, 47, 360, 50, 0.12),
                                             ('Kepler', 8.1, -38.0, 16, 240, 24, 0.2), ('Aristarchus', 23.7, -47.4, 20, 180, 18, 0.26),
                                             ('Proclus', 16.1, 46.9, 14, 200, 16, 0.14), ('Giordano Bruno', 35.9, 102.8, 11, 220, 18, 0.18)]:
        la0, lo0 = np.radians(clat), np.radians(clon)
        cosd = np.sin(la0) * np.sin(lat) + np.cos(la0) * np.cos(lat) * np.cos(lon - lo0)
        d = np.arccos(np.clip(cosd, -1, 1)) * 1737.4                   # km from the crater
        brg = np.arctan2(np.sin(lon - lo0) * np.cos(lat), np.cos(la0) * np.sin(lat) - np.sin(la0) * np.cos(lat) * np.cos(lon - lo0))
        rays = np.zeros_like(d)
        for _ in range(n):
            th = rng.uniform(-np.pi, np.pi); L = length * rng.uniform(0.3, 1.0); w = rng.uniform(12, 40)
            dth = (brg - th + np.pi) % (2 * np.pi) - np.pi
            across = np.abs(dth) * np.maximum(d, 1)
            rays = np.maximum(rays, rng.uniform(0.5, 1.0) * np.exp(-(across / w) ** 2) * np.exp(-d / L) * (d > rkm * 0.9) * (d < L * 2.5))
        halo = np.exp(-(np.maximum(d - rkm, 0) / (rkm * 1.6)) ** 2)
        a = a + k * np.maximum(rays * 0.75, halo * 0.9) * (1 - a * 0.5)
    img = np.clip(a * 235, 0, 255)
    mare = np.clip((0.62 - a) / 0.18, 0, 1)                                        # dark plains
    la0, lo0 = np.radians(8.5), np.radians(31.4)                                    # Mare Tranquillitatis
    dT = np.degrees(np.arccos(np.clip(np.sin(la0) * np.sin(lat) + np.cos(la0) * np.cos(lat) * np.cos(lon - lo0), -1, 1)))
    blue = mare * np.exp(-(dT / 16) ** 2)
    warm = mare * (1 - np.exp(-(dT / 16) ** 2))
    r = img * (1.0 - 0.035 * blue + 0.02 * warm); g = img * (0.985 - 0.01 * blue); b = img * (0.955 + 0.045 * blue - 0.02 * warm)
    rgb = np.clip(np.stack([r, g, b], axis=-1), 0, 255).astype(np.uint8)
    path.parent.mkdir(exist_ok=True)
    Image.fromarray(rgb, 'RGB').save(path, quality=88, optimize=True, progressive=True)


def mars():
    W, H = 1536, 768
    g = np.loadtxt(SRC / 'mars-elevation-1deg.csv', delimiter=',').astype(np.float32)   # (180, 360)
    pad = 4
    gp = np.concatenate([g[:, -pad:], g, g[:, :pad]], axis=1)          # wrap east-west before smooth resampling
    scale = W / 360
    big = Image.fromarray(gp, 'F').resize((int(round(gp.shape[1] * scale)), H), Image.BICUBIC)
    h = np.asarray(big, dtype=np.float32)[:, int(round(pad * scale)):int(round(pad * scale)) + W].copy()

    # Procedural small craters: the 1-degree grid is too coarse to show any
    rng = np.random.default_rng(11)
    lat_c = 90 - (np.arange(H) + 0.5) * 180 / H
    km_per_px_y = 2 * np.pi * 3389.5 / 2 / H
    count = 0
    for _ in range(9000):
        la = np.degrees(np.arcsin(rng.uniform(-1, 1)))
        lo = rng.uniform(-180, 180)
        if la > 25 and rng.random() < 0.75:          # the northern plains are younger and smoother
            continue
        if la > 80 or la < -82:
            continue
        D = 4 * (1 - rng.random()) ** (-1 / 1.6)       # km, power-law sizes
        if D > 160:
            continue
        depth = min(0.2 * D, 1.2 + 0.015 * D) * 1000
        rim = 0.25 * depth
        r_px = D / 2 / km_per_px_y
        cy = (90 - la) / 180 * H; cx = (lo + 180) / 360 * W
        cos_l = max(0.15, np.cos(np.radians(la)))
        ry, rx = r_px * 1.6, r_px * 1.6 / cos_l
        y0, y1 = int(max(0, cy - ry - 1)), int(min(H, cy + ry + 2))
        xs = np.arange(int(cx - rx - 1), int(cx + rx + 2))
        if y1 <= y0 or len(xs) < 2:
            continue
        yy = (np.arange(y0, y1)[:, None] + 0.5 - cy) / r_px
        xx = ((xs[None, :] + 0.5 - cx) * cos_l) / r_px
        r = np.sqrt(xx * xx + yy * yy)
        prof = np.where(r < 1, -depth * (1 - r ** 2) + rim * r ** 4, rim * np.exp(-((r - 1) / 0.28) ** 2))
        prof = np.where(r < 1.6, prof, 0)
        h[y0:y1][:, xs % W] += prof.astype(np.float32)
        count += 1
    save_height(h, HERE / 'mars' / 'mars-height.png')
    mars_colour(HERE / 'mars' / 'mars-color.jpg')
    print('mars height', h.min(), h.max(), 'craters', count)


def mars_colour(path, W=2048, H=1024):
    """Viking colour mosaic with its smeared polar regions replaced by the real cap shapes:
    the north residual cap with its spiral troughs and Chasma Boreale, and the small south residual cap,
    which sits about 3 degrees off the pole toward 315 degrees east."""
    from scipy.ndimage import gaussian_filter
    im = Image.open(SRC / 'mars-viking-color.jpg').convert('RGB').resize((W, H), Image.BICUBIC)
    im = im.filter(ImageFilter.UnsharpMask(radius=1.5, percent=45, threshold=2))
    a = np.asarray(im, dtype=np.float32) / 255
    lat = 90 - (np.arange(H) + 0.5) * 180 / H
    lon = -180 + (np.arange(W) + 0.5) * 360 / W
    LAT, LON = np.meshgrid(lat, lon, indexing='ij')
    def ss(e0, e1, x):
        t = np.clip((x - e0) / (e1 - e0), 0, 1); return t * t * (3 - 2 * t)
    rng = np.random.default_rng(3)
    noise = gaussian_filter(rng.standard_normal((H, W)).astype(np.float32), (3, 6), mode=('nearest', 'wrap'))
    noise /= noise.std()
    fine = gaussian_filter(rng.standard_normal((H, W)).astype(np.float32), 1.2, mode=('nearest', 'wrap')); fine /= fine.std()
    def polar_noise(sign, sigma, seed):
        """Smooth noise laid out around a pole, so it has no radial streaks from the map projection."""
        g = gaussian_filter(np.random.default_rng(seed).standard_normal((512, 512)).astype(np.float32), sigma)
        g /= g.std()
        colat = 90 - sign * LAT
        u = (colat * np.cos(np.radians(LON)) / 40 + 0.5) * 511
        v = (colat * np.sin(np.radians(LON)) / 40 + 0.5) * 511
        return g[np.clip(v, 0, 511).astype(int), np.clip(u, 0, 511).astype(int)]
    def blend(base, col, w):
        w = np.clip(w, 0, 1)[..., None]
        return base * (1 - w) + col * w
    # 1. the source mosaic is smeared above about 75 degrees: replace it with polar layered-terrain colours
    nN, nS = polar_noise(1, 6, 21), polar_noise(-1, 6, 22)
    nP = np.where(LAT > 0, nN, nS)
    lum = a.mean(axis=2)
    dk = ss(0.42, 0.22, lum)[..., None]                                  # Viking tints the dark regions blue; real ones are grey-brown
    a = a * (1 + dk * np.array([0.16, 0.05, -0.3], np.float32))
    ring_n = np.median(a[(LAT > 60) & (LAT < 68)], axis=0); ring_s = np.median(a[(LAT < -60) & (LAT > -68)], axis=0)
    ground = np.where((LAT > 0)[..., None], ring_n, ring_s) * (1.04 + 0.07 * nP)[..., None]   # the real colour just outside the smear
    a = blend(a, ground, ss(68, 76, np.abs(LAT) + nP * 1.2))
    ice = np.array([0.95, 0.93, 0.9], np.float32)
    # 2. north residual cap, about 1,000 km across, cut by spiral troughs that wind out from the pole, and Chasma Boreale
    colat = 90 - LAT
    th = np.radians(LON)
    fN = polar_noise(1, 2.5, 23)
    cap_n = ss(9.2, 7.6, colat + nN * 0.7)
    arms = np.sin(5 * th + colat * 1.15 + fN * 0.9)
    troughs = ss(0.66, 0.95, arms) * ss(1.2, 3.0, colat) * ss(-1.4, 0.4, fN)
    chasma = np.exp(-((((LON + 35 + 540) % 360) - 180) / (6 + colat * 0.6)) ** 2) * ss(9.5, 4.0, colat) * ss(1.5, 3.5, colat)
    trough_col = np.array([0.78, 0.64, 0.52], np.float32)
    capc = blend(np.broadcast_to(ice * (0.95 + 0.03 * fN)[..., None], a.shape), trough_col, troughs * 0.55)
    capc = blend(capc, ground, chasma * 0.95)
    dunes = ss(11.5, 9.5, colat) * ss(8.0, 9.8, colat) * ss(-0.2, 1.0, nN)          # dark dune collar (Olympia Undae)
    a = a * (1 - 0.32 * dunes[..., None])
    a = blend(a, capc, cap_n)
    # 3. south residual cap: small and offset about 3 degrees from the pole toward 315 degrees east
    la0, lo0 = np.radians(-87.0), np.radians(-45.0)
    d_s = np.degrees(np.arccos(np.clip(np.sin(la0) * np.sin(np.radians(LAT)) + np.cos(la0) * np.cos(np.radians(LAT)) * np.cos(np.radians(LON) - lo0), -1, 1)))
    a = blend(a, ice * (0.95 + 0.03 * nS)[..., None], ss(4.6, 3.2, d_s + nS * 0.5))
    a = blend(a, ice, ss(13, 6, d_s + nS) * 0.16)                                   # faint seasonal frost around it
    path.parent.mkdir(exist_ok=True)
    Image.fromarray(np.clip(a * 255, 0, 255).astype(np.uint8), 'RGB').save(path, quality=88, optimize=True, progressive=True)


if __name__ == '__main__':
    moon()
    mars()
    for d in ('moon', 'mars'):
        for f in sorted((HERE / d).iterdir()):
            print(f.name, f.stat().st_size // 1024, 'KB')
