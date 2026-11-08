# Globe textures for Analogtrix (Earth, Moon, Mars)

`build_earth_textures.py` turns real elevation data and satellite imagery into the
textures the page loads from `earth/`. Run it from this folder:

```bash
python3 build_earth_textures.py
```

It needs Python 3 with `numpy`, `Pillow` and `netCDF4`.

## Sources and licences

| File in `source/` | What it is | Where it came from | Licence |
|---|---|---|---|
| `earth-topography-10arcmin.nc` | NOAA ETOPO1 ice-surface elevation and bathymetry, 10 arc-minute grid (2161 × 1081), from a degree-720 spherical-harmonic model, so detail is about 25–30 km | Fatiando a Terra data release `fatiando-data/earth-topography-10arcmin` v1 (GitHub release asset). Original: NOAA NCEI, doi:10.7289/V5C8276M | Public domain |
| `earth-blue-marble.jpg` | NASA Blue Marble surface colour, 4096 × 2048 | `vasturiano/three-globe`, `example/img/earth-blue-marble.jpg`. Original: NASA Visible Earth | Public domain (NASA imagery) |
| `earth-water.png` | Water mask that matches the Blue Marble coastline, 1600 × 800 | `vasturiano/three-globe`, `example/img/earth-water.png` | Derived from NASA imagery |

Download commands:

```bash
curl -L -o source/earth-topography-10arcmin.nc https://github.com/fatiando-data/earth-topography-10arcmin/releases/download/v1/earth-topography-10arcmin.nc
curl -L -o source/earth-blue-marble.jpg https://raw.githubusercontent.com/vasturiano/three-globe/master/example/img/earth-blue-marble.jpg
curl -L -o source/earth-water.png https://raw.githubusercontent.com/vasturiano/three-globe/master/example/img/earth-water.png
```

## Outputs (`earth/`)

| File | Size | Use |
|---|---|---|
| `earth-height.png` | 2048 × 1024 | Height in 16 bits: `R × 256 + G − 11000` metres. The page builds displacement, lighting normals, the Terrain view and the Atmospheric pressure view from it. |
| `earth-water.png` | 2048 × 1024 | Smooth water coverage (255 = water) for coastlines, ocean shine and the data views |
| `earth-color.jpg` | 4096 × 2048 | Blue Marble land, with the ocean re-coloured from real depth so continental shelves look lighter and the abyss darker |
| `earth-clouds.jpg` | 2048 × 1024 | Procedural cloud opacity: domain-warped noise with spiral storms at the usual low-pressure centres, denser under low pressure and thinner under the subtropical highs |

## Limits

- The 10 arc-minute model smooths the highest peaks: Everest reads about 6 km instead of 8.8 km. The site cards in the page use real heights.
- A finer grid (ETOPO 2022 at 30 or 60 arc-seconds) would sharpen ridges and valleys. The build script works with any equirectangular NetCDF grid with `latitude`, `longitude` and `topography` variables.
- The clouds are procedural, not a satellite cloud image, because no clearly licensed cloud texture was available.
- Vertical exaggeration is set once in the page script as `EARTH_RELIEF_EXAGGERATION`.

## Moon and Mars

`build_moon_mars_textures.py` builds `moon/` and `mars/` the same way:

```bash
python3 build_moon_mars_textures.py
```

| File in `source/` | What it is | Where it came from | Licence |
|---|---|---|---|
| `moon-topography-10arcmin.nc` | Lunar topography (LOLA), 10 arc-minute grid | Fatiando a Terra data release `fatiando-data/moon-topography-10arcmin` v1. Model: Wieczorek (2015), doi:10.5281/zenodo.3870924 | CC BY (credit Wieczorek 2015) |
| `moon-albedo-1024.jpg` | Moon colour map, 1024 × 512 | `mrdoob/three.js`, `examples/textures/planets/moon_1024.jpg` | MIT (three.js) |
| `mars-elevation-1deg.csv` | Mars elevation (MOLA), 1-degree grid, 180 × 360 | `AI4EPS/EPS88_PyEarth`, `data/mars_elevation.csv` | NASA MOLA data, public domain |
| `mars-viking-color.jpg` | USGS Viking colour mosaic, 1024 × 512 | `DaveSteadman/KoreGlobe`, `Resources/mars_viking_clrmosaic_global_1024.jpg` | USGS / NASA, public domain |

The Moon colour is an albedo map: the source image is split into large-scale brightness (maria,
highlands) and fine detail, and only a third of the fine detail is kept, so bright young craters stay
while most baked-in crater shadows fade; the eastern maria get a faint blue tint and the western maria a faint warm one, the polar smear is faded out, and bright ray systems are added around
Tycho, Copernicus, Kepler, Aristarchus, Proclus and Giordano Bruno at their real positions. All relief
shading on the Moon comes from the LOLA heights, lit in the page with a Lommel–Seeliger model; the
Moon displacement uses 2× exaggeration (`RELIEF.moon.exag`) and its normals use true scale.

Limits: the Mars grid is coarse (about 60 km), so the script adds procedural small craters
(fewer on the younger northern plains). The real LROC WAC colour mosaic would give sharper albedo
detail than the blurred source map.
Higher-resolution NASA maps (LROC colour, LOLA and MOLA DEMs) are stored in Git LFS or on
NASA/USGS servers that were not reachable when these files were built; the scripts accept them
if you replace the source files.

## Earth heights and atmospheric pressure (hover readout, Terrain and Atmospheric pressure views)

Three scripts, run in this order:

```bash
python3 build_earth_dem.py --download    # AWS Terrain Tiles -> earth/earth-dem-0.1.b64.txt, earth/dem/*.b64.txt, earth/earth-height.png
python3 fetch_era5_pressure.py           # ERA5 1991-2020 samples over HTTP range requests -> source/era5-climatology.npz (~3 min)
python3 build_earth_pressure.py          # -> earth/earth-pressure.b64.txt
python3 build_earth_dem.py               # run again once the ERA5 file exists: it uses ERA5 orography inside 85 deg S
```

Needs numpy, Pillow, netCDF4, scipy, h5py, fsspec and aiohttp.

| Data | Source | Licence / credit |
|---|---|---|
| Heights (hover readout, relief) | Terrain Tiles on AWS Open Data (`s3://elevation-tiles-prod`, Terrarium PNG, zoom 5, about 4.9 km), built from SRTM, GMTED2010, ETOPO1 bathymetry and national DEMs | Open data; attribution required: see https://github.com/tilezen/joerd/blob/master/docs/attribution.md |
| Greenland and Antarctic ice surface | NOAA ETOPO1 (ice surface), 10 arc-minute grid from Fatiando a Terra; the tiles hold bedrock under the ice sheets | Public domain |
| Polar cap inside 85 deg S | ERA5 orography (the ETOPO1 10-arcmin grid degrades near the poles) | As ERA5 below |
| Surface pressure, 2 m temperature, orography, land-sea mask | ECMWF ERA5 reanalysis, NSF NCAR copy on AWS (`s3://nsf-ncar-era5`, ds633.0), 0.25 deg, 1991-2020: 120 sampled days for pressure (each month 10 times), 30 for temperature, each a 27-hour mean | Copernicus licence. Contains modified Copernicus Climate Change Service information |

How the page uses them:

- **Height:** bilinear interpolation of the true DEM, in metres with no exaggeration, rounded to 10 m.
  The global 0.1 deg grid (11 km) loads the first time Terrain or Atmospheric pressure is opened; a 0.05 deg
  (5.5 km) bundle for the 30 x 30 deg square under the cursor loads on hover. The 3D relief uses the
  same heights, exaggerated by `EARTH_RELIEF_EXAGGERATION`, but readouts never read the display texture.
- **Pressure:** ERA5 annual-mean surface pressure is reduced to sea level with ERA5's own terrain and
  temperature, then brought back down to the true ground height with the hypsometric equation:
  `p = p_sl * exp(-g h / (Rd (T_sl - 0.0065 h / 2)))`. At ERA5's own terrain this returns ERA5 exactly;
  where the finer DEM differs (valleys, the Dead Sea) the equation adjusts it. "% of sea level" is p / 1013.25.
- **Picking:** one function, `screenToLatLon()`, serves the 3D globe and the 2D map. On the globe it marches
  the mouse ray through the displaced relief, using the same per-vertex displacement the GPU draws, so the
  point under the cursor and the readout always agree. Inverses: `latLonToUV`, `latLonToWorld`, `latLonToScreen`.

### Tests

Open the page with `#refcheck` for an overlay table, or run in the browser console:

```js
await analogtrixReferenceTest()   // 13 heights and 9 pressures against reference values
analogtrixCoordTest()             // lat/lon <-> UV and lat/lon -> ray pick -> lat/lon round trips
```

Known gap: the reference "Mid North Atlantic 30N 40W, -4,000 to -5,500 m" reads about -3,400 m. That point
is on the flank of the Mid-Atlantic Ridge; ETOPO1 independently gives -3,300 to -3,900 m there.
