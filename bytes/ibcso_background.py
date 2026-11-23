"""ibcso_background.py -- COG of the IBCSO v1.0 shaded background image.

IBCSO_background_hq.zip (raad mirror of hs.pangaea.de/Maps/bathy/IBCSO_v1/,
2013): an RGB render of IBCSO v1.0 at 541.57 m in polar stereographic, true
scale 65 deg S on WGS 84 (EPSG:9354), with its own world file and CRS. The
COG keeps the pixels and georeferencing as published; compared with the v1
ice surface, the ocean edge lies within about 1 km (two pixels) of the grid.

Usage:
  python bytes/ibcso_background.py --tif RAW/unz/IBCSO_background_hq/ibcso_background_hq.tif \
    --fetch RAW/fetch.tsv --out OUT
"""
import argparse
import hashlib
import json
import os
import sys

import numpy as np
from osgeo import gdal

gdal.UseExceptions()


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--tif", required=True)
    p.add_argument("--fetch", required=True)
    p.add_argument("--out", required=True)
    a = p.parse_args()
    os.makedirs(a.out, exist_ok=True)
    with open(a.fetch) as f:
        fetched = {r[0]: r[1:] for r in (l.rstrip("\n").split("\t") for l in f)}
    h = sha256(a.tif)
    src = gdal.Open(a.tif)
    out = os.path.join(a.out, "ibcso-v1-background-hq.tif")
    gdal.Translate(out, src, format="COG", metadataOptions=[
        "SOURCE_DOI=10.1594/PANGAEA.805734", "SOURCE_FILE=ibcso_background_hq.tif",
        "SOURCE_ZIP=IBCSO_background_hq.zip", f"SOURCE_SHA256={h}"], creationOptions=[
        "COMPRESS=ZSTD", "LEVEL=17", "PREDICTOR=YES", "BLOCKSIZE=512",
        "OVERVIEW_RESAMPLING=AVERAGE", "NUM_THREADS=ALL_CPUS"])
    got = gdal.Open(out)
    if not np.array_equal(got.ReadAsArray(), src.ReadAsArray()) or got.GetGeoTransform() != src.GetGeoTransform():
        sys.exit("ibcso_background: the COG differs from the source")
    doc = {"source": {"file": "ibcso_background_hq.tif", "sha256": h, "zip": "IBCSO_background_hq.zip",
                      "zip_sha256": fetched["IBCSO_background_hq.zip"][1],
                      "origin": "raad mirror of hs.pangaea.de/Maps/bathy/IBCSO_v1/"},
           "license": "CC-BY-3.0", "cog": {"file": os.path.basename(out), "bytes": os.path.getsize(out),
                                         "sha256": sha256(out), "pixels_different": 0},
           "georeferencing": "as published (EPSG:9354 equivalent); about 1 km from the v1 ice surface "
                             "ocean edge", "gdal": gdal.__version__}
    with open(os.path.join(a.out, "provenance.json"), "w") as f:
        json.dump(doc, f, indent=1)
    print(f"background: {os.path.getsize(out) / 1e6:.1f} MB, identical to the source")


if __name__ == "__main__":
    main()
