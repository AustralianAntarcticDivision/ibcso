"""ibcso_v1_cog.py -- Cloud Optimized GeoTIFFs of the IBCSO v1.0 grids.

IBCSO v1.0 (Arndt et al. 2013, doi:10.1594/PANGAEA.805734, CC BY 3.0), from
the raad mirror of hs.pangaea.de/Maps/bathy/IBCSO_v1/: three GMT grids at
500 m in polar stereographic, true scale 71 deg S (EPSG:3031), 13335 x 13341,
Float32. The grids carry no CRS; EPSG:3031 is stated by the file names (PS71)
and the GeoTIFF of the ice surface, IBCSO_v1_is.tif, in the same release.

  bed          ibcso_v1_bed.grd   GMT native netCDF (1D z, x_range, dimension)
  ice_surface  ibcso_v1_is.grd    GMT native netCDF
  sid          ibcso_v1_sid.grd   COARDS netCDF (source identifier)

The native grids are read as rows from the top. Their nodata is +inf, kept
as stored so the pixels are bit-identical; sid's nodata is -3.4028235e+38.
Each COG is compared with its grid (NaN equal to NaN, inf to inf); any
differing pixel is a failure. ice_surface is also compared with
IBCSO_v1_is.tif: the nodata mask must match exactly, and the values agree
within 0.18 m without being identical (two encodings of the same surface in
the release); the maximum difference is recorded in provenance.json.

Settings are those of the Bedmap3 COGs (raad.virtualize py/bedmap3_cog.py): ZSTD level 17, no predictor, 512 x
512 tiles, SPARSE_OK=FALSE; overviews AVERAGE, or NEAREST for sid.

Usage:
  python bytes/ibcso_v1_cog.py --raw RAW/unz --fetch RAW/fetch.tsv --out OUT
"""
import argparse
import hashlib
import json
import os
import sys
import time

import netCDF4
import numpy as np
from osgeo import gdal, osr

gdal.UseExceptions()

DOI = "10.1594/PANGAEA.805734"
GT = (-3333500.0, 500.0, 0.0, 3333500.0, 0.0, -500.0)
GRIDS = {  # array: grid file, zip, overview resampling, long name
    "bed": ("ibcso_v1_bed.grd", "IBCSO_v1_bed_PS71_500m_grd.zip", "AVERAGE", "Bed elevation"),
    "ice_surface": ("ibcso_v1_is.grd", "IBCSO_v1_is_PS71_500m_grd.zip", "AVERAGE", "Ice surface elevation"),
    "sid": ("ibcso_v1_sid.grd", "IBCSO_v1_sid_PS71_500m_grd.zip", "NEAREST", "Source identifier"),
}


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def read_native(path):
    """GMT native netCDF grid as a (ny, nx) array, rows from the top, with
    its node-registered x_range and y_range checked against GT."""
    d = netCDF4.Dataset(path)
    d.set_auto_maskandscale(False)
    nx, ny = (int(v) for v in d["dimension"][:])
    xr, yr = d["x_range"][:], d["y_range"][:]
    if int(d["z"].getncattr("node_offset")) != 0 or \
            abs(xr[0] - (GT[0] + 250)) > 1e-3 or abs(yr[1] - (GT[3] - 250)) > 1e-3:
        sys.exit(f"ibcso_v1_cog: {path}: x_range {xr}, y_range {yr} do not match the grid")
    return np.asarray(d["z"][:]).reshape(ny, nx)


def same(a, b):
    return np.array_equal(a, b, equal_nan=True)


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--raw", required=True, help="directory of the unzipped grids")
    p.add_argument("--fetch", required=True, help="fetch.tsv of the zips")
    p.add_argument("--out", required=True)
    a = p.parse_args()
    os.makedirs(a.out, exist_ok=True)
    with open(a.fetch) as f:
        fetched = {r[0]: r[1:] for r in (l.rstrip("\n").split("\t") for l in f)}
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(3031)
    files = {}
    for name, (grid, zipname, resampling, long_name) in GRIDS.items():
        t0 = time.time()
        src = os.path.join(a.raw, grid)
        h = sha256(src)
        if grid == "ibcso_v1_sid.grd":
            s = gdal.Open(src)
            # the netCDF driver reports the origin with floating noise
            # (-3333500.0000000014); within 1 micrometre it is GT
            if not np.allclose(s.GetGeoTransform(), GT, rtol=0, atol=1e-6):
                sys.exit(f"ibcso_v1_cog: {grid} geotransform {s.GetGeoTransform()}")
            z = s.ReadAsArray()
            nodata = float(np.float32(-3.4028234663852886e+38))
        else:
            z = read_native(src)
            nodata = float("inf")
        mem = gdal.GetDriverByName("MEM").Create("", z.shape[1], z.shape[0], 1, gdal.GDT_Float32)
        mem.SetGeoTransform(GT)
        mem.SetSpatialRef(srs)
        b = mem.GetRasterBand(1)
        b.WriteArray(z)
        b.SetNoDataValue(nodata)
        b.SetDescription(name)
        mem.SetMetadata({"SOURCE_DOI": DOI, "SOURCE_FILE": grid, "SOURCE_ZIP": zipname,
                         "SOURCE_SHA256": h, "IBCSO_RELEASE": "v1", "long_name": long_name,
                         "SRS_SOURCE": "EPSG:3031 from the file names (PS71) and IBCSO_v1_is.tif"})
        out = os.path.join(a.out, f"{name}.tif")
        gdal.Translate(out, mem, format="COG", creationOptions=[
            "COMPRESS=ZSTD", "LEVEL=17", "PREDICTOR=NO", "BLOCKSIZE=512", "SPARSE_OK=FALSE",
            f"OVERVIEW_RESAMPLING={resampling}", "NUM_THREADS=ALL_CPUS"])
        got = gdal.Open(out)
        v = got.ReadAsArray()
        if v.dtype != np.float32 or not same(v, z.astype(np.float32)) or got.GetGeoTransform() != GT:
            sys.exit(f"ibcso_v1_cog: {name}: the COG differs from {grid}")
        entry = {"file": f"{name}.tif", "source": grid, "source_zip": zipname, "source_sha256": h,
                 "zip_sha256": fetched[zipname][1], "bytes": os.path.getsize(out), "sha256": sha256(out),
                 "pixels_compared": int(v.size), "pixels_different": 0, "nodata": str(nodata),
                 "overview_resampling": resampling}
        if name == "ice_surface":
            t = gdal.Open(os.path.join(a.raw, "ibcso_v1_is.tif")).ReadAsArray()
            tv, zv = t > -1e38, np.isfinite(z)
            if not np.array_equal(tv, zv):
                sys.exit("ibcso_v1_cog: ice_surface nodata mask differs from IBCSO_v1_is.tif")
            d = float(np.max(np.abs(t[tv] - z[zv])))
            entry["compared_with_geotiff"] = {"file": "ibcso_v1_is.tif", "mask_identical": True,
                                              "max_abs_difference_m": d,
                                              "cells_identical": int((t[tv] == z[zv]).sum()),
                                              "cells_valid": int(tv.sum())}
        files[name] = entry
        print(f"{name}: {entry['bytes'] / 1e6:.1f} MB, identical to {grid}, {time.time() - t0:.0f} s",
              flush=True)
    doc = {"release": "v1", "source": {"doi": f"https://doi.org/{DOI}",
                                       "origin": "raad mirror of hs.pangaea.de/Maps/bathy/IBCSO_v1/"},
           "license": "CC-BY-3.0", "gdal": gdal.__version__, "cogs": files,
           "fetched": {k: {"bytes": int(v[0]), "sha256": v[1], "url": v[2]} for k, v in fetched.items()}}
    with open(os.path.join(a.out, "provenance.json"), "w") as f:
        json.dump(doc, f, indent=1)


if __name__ == "__main__":
    main()
