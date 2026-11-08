# Analogtrix — The gate to the space (Cửa ngõ vươn ra vũ trụ) — source code

Earth analogs for permanent Moon and Mars bases: an interactive site with 3D globes,
2D maps, place details and a dictionary of sites, in English and Vietnamese.

## Folders

- `site/` — the website
  - `index.html` — all the page code: HTML, CSS (in `<style>`) and JavaScript (in the last `<script>`).
    Place data (Moon, Mars and Earth sites, fun facts, terrain and climate numbers) is near the top of that script.
    Vietnamese text is in `I18N_VI` (interface), `SITE_VI` (site text by id) and `BODIES_VI`, just above the data.
  - `earth/`, `moon/`, `mars/` — map images the globes and 2D maps load.
- `build/` — how the map images were made
  - `build_earth_textures.py` — Earth colour, height, water and cloud maps
  - `build_moon_mars_textures.py` — Moon and Mars colour and height maps
  - `build_earth_dem.py` — true Earth heights (AWS Terrain Tiles) for the hover readout and relief
  - `fetch_era5_pressure.py`, `build_earth_pressure.py` — ERA5 1991–2020 annual-mean surface pressure
  - `make_single_file.py` — packs `site/` into one self-contained HTML file
  - `source/` — the raw data (NOAA ETOPO1, NASA Blue Marble, LOLA, MOLA, USGS Viking)
  - `README.md` — data sources, licences and limits

## Run the site

Serve the `site/` folder with any web server (opening index.html directly blocks the map images):

    cd site
    python3 -m http.server 8000

Then open http://localhost:8000.

## Rebuild the map images

Needs Python 3 with numpy, Pillow, netCDF4 and scipy:

    cd build
    python3 build_earth_textures.py
    python3 build_moon_mars_textures.py

They write into `build/earth`, `build/moon`, `build/mars`; copy those folders into `site/`.

## Make the single shareable file

    python3 build/make_single_file.py
