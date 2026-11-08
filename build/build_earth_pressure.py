#!/usr/bin/env python3
"""Turn the ERA5 climatology into the compact pressure file the page loads.

Input:  source/era5-climatology.npz  (from fetch_era5_pressure.py)
Output: earth/earth-pressure.b64.txt base64 of zlib-deflated data: three grids of 1440 x 721 (0.25 deg,
                                     node-registered: lon -180..179.75, lat 90..-90):
                                       1. p_sl  Int16 LE, 0.1 hPa   pressure reduced to sea level
                                       2. T_sl  Int16 LE, 0.1 K     temperature reduced to sea level
                                       3. lsm   Uint8, 0-100        ERA5 land fraction

How the page gets the surface pressure at a point with ground height h (m, from the DEM):
    T_sl, p_sl = bilinear(grids, lat, lon)
    p = p_sl * exp(-g * h / (Rd * (T_sl - GAMMA * h / 2)))
Ground below sea level counts as land only where ERA5's land fraction is at least 0.5 (the Dead
Sea is 0.60, a coastal shelf about 0.37); any other point with h < 0 is sea, at h = 0.
The reduction here uses the same formula from ERA5's own orography z, so at h = z the page gets
back exactly ERA5's annual-mean surface pressure; where the finer DEM differs from ERA5's
smoothed terrain (valleys, the Dead Sea, ice-sheet domes) the hypsometric equation adjusts it.
"""
import base64, pathlib, zlib
import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
G, RD, GAMMA = 9.80665, 287.05, 0.0065


def fields():
    d = np.load(HERE / 'source' / 'era5-climatology.npz')
    sp, t2m, z = d['sp'] / 100.0, d['t2m'], d['z']          # hPa, K, m
    T_sl = t2m + GAMMA * z
    p_sl = sp * np.exp(G * z / (RD * (t2m + GAMMA * z / 2)))
    # ERA5 longitudes run 0..359.75; shift to -180..179.75
    roll = lambda a: np.roll(a, 720, axis=1)
    return roll(p_sl), roll(T_sl), roll(sp), roll(z), roll(t2m), roll(d['lsm'])


def surface_pressure(p_sl, T_sl, h):
    return p_sl * np.exp(-G * h / (RD * (T_sl - GAMMA * h / 2)))


def main():
    p_sl, T_sl, sp, z, t2m, lsm = fields()
    out = np.round(p_sl * 10).astype('<i2').tobytes() + np.round(T_sl * 10).astype('<i2').tobytes() + np.round(lsm * 100).astype(np.uint8).tobytes()
    (HERE / 'earth' / 'earth-pressure.b64.txt').write_bytes(base64.b64encode(zlib.compress(out, 9)))
    print('earth-pressure.b64.txt', round((HERE / 'earth' / 'earth-pressure.b64.txt').stat().st_size / 1e6, 2), 'MB;',
          'p_sl range', round(float(p_sl.min()), 1), round(float(p_sl.max()), 1))


if __name__ == '__main__':
    main()
