#!/usr/bin/env python3
"""Download an annual-mean surface-pressure climatology (1991-2020) from ERA5.

Source: ECMWF ERA5 reanalysis, copy hosted by NSF NCAR on AWS Open Data
(s3://nsf-ncar-era5, dataset ds633.0). Licence: Copernicus licence; credit
"Contains modified Copernicus Climate Change Service information (1991-2020)".

ERA5 stores hourly fields in ~1 GB monthly NetCDF files, compressed in blocks of
27 hours. Instead of downloading whole files, this script reads one 27-hour block
per sampled month over HTTP byte ranges and averages it (which also averages out
the daily pressure tide).

Sampling: 30 years x 4 months for surface pressure (each calendar month 10 times,
spread over the whole 1991-2020 period); 30 samples for 2 m temperature.
Synoptic noise in the mean is about 1 hPa.

Output: source/era5-climatology.npz with float32 grids (721 x 1440, 0.25 deg,
lat 90..-90, lon 0..359.75):  sp (Pa), t2m (K), z (orography, m), lsm (0-1).

Needs: numpy, h5py, fsspec, aiohttp.   Run:  python3 fetch_era5_pressure.py
"""
import concurrent.futures as cf, pathlib, sys, time
import numpy as np, fsspec, h5py

BASE = 'https://nsf-ncar-era5.s3.amazonaws.com/'
OUT = pathlib.Path(__file__).resolve().parent / 'source' / 'era5-climatology.npz'
FS = fsspec.filesystem('http', block_size=2 ** 22)


def month_file(var, code, y, m):
    nd = [31, 29 if y % 4 == 0 else 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31][m - 1]
    return f'{BASE}e5.oper.an.sfc/{y}{m:02d}/e5.oper.an.sfc.{code}_{var}.ll025sc.{y}{m:02d}0100_{y}{m:02d}{nd}23.nc'


def block_mean(url, key):
    for attempt in range(5):
        try:
            with FS.open(url, 'rb', cache_type='blockcache') as f:
                d = h5py.File(f, 'r')[key]
                return d[0:d.chunks[0]].astype(np.float64).mean(axis=0)
        except Exception as e:  # network hiccup: retry
            err = e
            time.sleep(2 ** attempt)
    raise RuntimeError(f'{url}: {err}')


def mean_of(jobs, key, label):
    acc, n, t0 = None, 0, time.time()
    with cf.ThreadPoolExecutor(6) as ex:
        for a in ex.map(lambda u: block_mean(u, key), jobs):
            acc = a if acc is None else acc + a
            n += 1
            print(f'{label} {n}/{len(jobs)}  {time.time() - t0:.0f}s', flush=True)
    return (acc / n).astype(np.float32)


def invariant(code, var, key):
    url = f'{BASE}e5.oper.invariant/197901/e5.oper.invariant.{code}_{var}.ll025sc.1979010100_1979010100.nc'
    with FS.open(url, 'rb') as f:
        return np.asarray(h5py.File(f, 'r')[key][0], dtype=np.float32)


def main():
    years = range(1991, 2021)
    sp_jobs = [month_file('sp', '128_134', y, m) for i, y in enumerate(years) for m in range(1, 13) if (m - i) % 3 == 0]
    t_jobs = [month_file('2t', '128_167', y, (i * 5) % 12 + 1) for i, y in enumerate(years)]
    z = invariant('128_129', 'z', 'Z') / 9.80665
    lsm = invariant('128_172', 'lsm', 'LSM')
    sp = mean_of(sp_jobs, 'SP', 'surface pressure')
    t2m = mean_of(t_jobs, 'VAR_2T', '2 m temperature')
    np.savez_compressed(OUT, sp=sp, t2m=t2m, z=z, lsm=lsm, n_sp=len(sp_jobs), n_t=len(t_jobs))
    print('wrote', OUT)


if __name__ == '__main__':
    main()
