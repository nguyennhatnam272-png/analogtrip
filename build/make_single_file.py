#!/usr/bin/env python3
"""Make one self-contained HTML file from site/: embeds the map images and three.js.

Run from the project folder:  python3 build/make_single_file.py
Needs internet once, to download three.js r128 (MIT licence) from cdnjs.
Writes analogtrix.html next to the site/ folder.
"""
import base64, json, pathlib, urllib.request
root = pathlib.Path(__file__).resolve().parent.parent
site = root / 'site'
s = (site / 'index.html').read_text()
url = 'https://cdnjs.cloudflare.com/ajax/libs/three.js/r128/three.min.js'
three = urllib.request.urlopen(url).read().decode()
tag = '<script src="' + url + '"></script>'
s = s.replace(tag, '<script>/* three.js r128, MIT licence */\n' + three + '\n</script>')
mime = {'.jpg': 'image/jpeg', '.png': 'image/png'}
embed = {f'{d}/{f.name}': 'data:' + mime[f.suffix] + ';base64,' + base64.b64encode(f.read_bytes()).decode()
         for d in ('earth', 'moon', 'mars') for f in sorted((site / d).iterdir()) if f.suffix in mime}
# the global height grid and the pressure data (already base64 text); the 0.05 deg detail bundles are left out (about 14 MB)
for k in ('earth/earth-dem-0.1.b64.txt', 'earth/earth-pressure.b64.txt'):
    embed[k] = 'data:application/octet-stream;base64,' + (site / k).read_text().strip()
i = s.index('<script>/* three.js r128')
s = s[:i] + '<script>window.ANALOGTRIX_IMAGES = ' + json.dumps(embed) + ';</script>\n' + s[i:]
s = s.replace("  function loadImage(src) {\n", "  function loadImage(src) {\n    if (window.ANALOGTRIX_IMAGES && window.ANALOGTRIX_IMAGES[src]) src = window.ANALOGTRIX_IMAGES[src];\n", 1)
out = root / 'analogtrix.html'
out.write_text(s)
print('wrote', out, round(out.stat().st_size / 1e6, 1), 'MB')
