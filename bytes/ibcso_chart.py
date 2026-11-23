"""ibcso_chart.py -- georeferenced COG of an IBCSO digital chart PDF.

The IBCSO chart PDFs carry no georeferencing. Each page is drawn in polar
stereographic, true scale at 65 deg S on WGS 84 (EPSG:9354), with the South
Pole at the centre of a disc whose rim is a parallel: v2 at 1:11,500,000 to
50 deg S (doi:10.1594/PANGAEA.937574), v1 at 1:7,000,000 to 60 deg S
(doi:10.1594/PANGAEA.805735). This script renders the page at --dpi, takes
the pixel size from the printed scale (--scale), locates the pole as the
centre of a circle fitted to the edge of the ocean disc, and checks the rim
stroke radius against the projected radius of --rim-lat: a difference of
more than 0.1% stops the script.

Outputs, in OUT:
  <name>.tif          full page as a COG (JPEG, YCbCr), margins included
  <name>.gdalg.json   GDAL streamed algorithm: the COG at --url, clipped to
                      neatline-<lat>S.geojson at --url, with an alpha band
                      (GDAL >= 3.11)
  neatline-<lat>S.geojson  the rim parallel as a polygon in EPSG:9354
  <name>.json         STAC item (roles visual, overview)
  <name>.provenance.json  source checksum, fit and georeferencing

Usage:
  python bytes/ibcso_chart.py --pdf RAW/IBSCO_v2_rev1_digital_chart.pdf \
    --name ibcso-v2-rev1-chart --title "..." --datetime 2024-08-20T12:14:01Z \
    --url https://.../v2/chart --out OUT [--dpi 300]
"""
import argparse
import hashlib
import json
import os
import sys

import numpy as np
from osgeo import gdal, osr

gdal.UseExceptions()
EPSG = 9354


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def fit_circle(P):
    A = np.column_stack([P[:, 0], P[:, 1], np.ones(len(P))])
    D, E, F = np.linalg.lstsq(A, -(P[:, 0] ** 2 + P[:, 1] ** 2), rcond=None)[0]
    cx, cy = -D / 2, -E / 2
    return cx, cy, np.sqrt(cx * cx + cy * cy - F)


def disc_edge(rgb, guess, rmin, rmax, n=1440):
    """Outermost ocean (blue) pixel along n rays from guess, pixel-corner
    coordinates; a robust circle fit to those points."""
    r, b = rgb[0].astype(np.int16), rgb[2].astype(np.int16)
    blue = (b - r > 40) & (b > 60)
    H, W = blue.shape
    pts = []
    for th in np.linspace(0, 2 * np.pi, n, endpoint=False):
        rr = np.arange(rmin, rmax)
        x = (guess[0] + rr * np.cos(th)).astype(int)
        y = (guess[1] + rr * np.sin(th)).astype(int)
        ok = (x >= 0) & (x < W) & (y >= 0) & (y < H)
        hit = np.flatnonzero(blue[y[ok], x[ok]])
        if len(hit):
            pts.append((x[ok][hit[-1]] + 0.5, y[ok][hit[-1]] + 0.5))
    P = np.array(pts, float)
    for _ in range(10):
        cx, cy, R = fit_circle(P)
        res = np.hypot(P[:, 0] - cx, P[:, 1] - cy) - R
        keep = np.abs(res) < max(3 * 1.4826 * np.median(np.abs(res)), 1.0)
        if keep.all():
            break
        P = P[keep]
    return cx, cy, R, float(res[keep].std()), len(P), blue


def rim_stroke(rgb, blue, cx, cy, R, n=2880):
    """Radius of the dark rim stroke just outside the disc edge (median over
    rays), pixel-corner coordinates."""
    s = rgb.astype(np.int16).sum(0)
    out = []
    for th in np.linspace(0, 2 * np.pi, n, endpoint=False):
        rr = np.arange(R - 6, R + 12, 0.25)
        x = (cx + rr * np.cos(th)).astype(int)
        y = (cy + rr * np.sin(th)).astype(int)
        v = s[y, x]
        w = np.clip(765 - v, 0, None).astype(float)
        w[(v > 700) | blue[y, x]] = 0
        if w.sum() > 0 and not blue[y[-20:], x[-20:]].any():
            out.append((rr * w).sum() / w.sum())
    out = np.array(out)
    med = np.median(out)
    good = np.abs(out - med) < 2
    return float(np.median(out[good])), float(out[good].std()), int(good.sum())


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--pdf", required=True)
    p.add_argument("--name", required=True)
    p.add_argument("--title", required=True)
    p.add_argument("--datetime", required=True, help="RFC 3339 date of the chart")
    p.add_argument("--url", required=True, help="public URL of OUT")
    p.add_argument("--out", required=True)
    p.add_argument("--dpi", type=int, default=300)
    p.add_argument("--scale", type=float, default=11_500_000, help="printed scale denominator")
    p.add_argument("--rim-lat", type=float, default=-50.0, help="latitude of the disc rim")
    p.add_argument("--doi", default="10.1594/PANGAEA.937574")
    p.add_argument("--source-href", help="URL of the PDF (default: PANGAEA v2 files)")
    p.add_argument("--license", default="CC-BY-4.0")
    a = p.parse_args()
    SCALE, LAT_RIM = a.scale, a.rim_lat
    neat_name = f"neatline-{abs(LAT_RIM):g}S"
    source_href = a.source_href or f"https://download.pangaea.de/dataset/937574/files/{os.path.basename(a.pdf)}"
    os.makedirs(a.out, exist_ok=True)

    gdal.SetConfigOption("GDAL_PDF_DPI", str(a.dpi))
    src = gdal.Open(a.pdf)
    rgb = src.ReadAsArray()[:3]
    H, W = rgb.shape[1:]
    px = 0.0254 / a.dpi * SCALE  # metres per pixel at the printed scale
    k = a.dpi / 150
    cx, cy, R, rms, nfit, blue = disc_edge(rgb, (W / 2, H * 0.42), int(1500 * k), int(3000 * k))
    rim, rim_sd, nrim = rim_stroke(rgb, blue, cx, cy, R)

    srs = osr.SpatialReference()
    srs.ImportFromEPSG(EPSG)
    geo = srs.CloneGeogCS()
    for s in (srs, geo):
        s.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    to_xy = osr.CoordinateTransformation(geo, srs)
    r50 = abs(to_xy.TransformPoint(0.0, LAT_RIM)[1])
    err = rim * px / r50 - 1
    print(f"{a.name}: {W} x {H} px at {a.dpi} dpi, {px:.4f} m/px; pole ({cx:.2f}, {cy:.2f}) px, "
          f"disc fit rms {rms:.2f} px ({nfit} rays); rim {rim:.2f} px (sd {rim_sd:.2f}, {nrim} rays), "
          f"{rim * px:.0f} m against {r50:.0f} m for {neat_name[9:]} ({100 * err:+.3f}%)", flush=True)
    if abs(err) > 0.001:
        sys.exit(f"ibcso_chart: the rim is {100 * err:+.3f}% off {neat_name[9:]} at 1:{SCALE:g}; not georeferenced")

    gt = (-cx * px, px, 0.0, cy * px, 0.0, -px)
    h = sha256(a.pdf)
    tif = os.path.join(a.out, f"{a.name}.tif")
    mem = gdal.Translate("", src, format="MEM", bandList=[1, 2, 3])
    mem.SetGeoTransform(gt)
    mem.SetSpatialRef(srs)
    mem.SetMetadata({"SOURCE_DOI": a.doi, "SOURCE_FILE": os.path.basename(a.pdf),
                     "SOURCE_SHA256": h, "TITLE": a.title,
                     "GEOREFERENCING": f"1:{SCALE:g} at {a.dpi} dpi; pole fitted to the {neat_name[9:]} disc "
                                       f"({cx:.2f}, {cy:.2f}) px; rim {100 * err:+.3f}% from {neat_name[9:]} "
                                       "(py/ibcso_chart.py)"})
    for i, c in enumerate((gdal.GCI_RedBand, gdal.GCI_GreenBand, gdal.GCI_BlueBand), 1):
        mem.GetRasterBand(i).SetColorInterpretation(c)
    gdal.Translate(tif, mem, format="COG", creationOptions=[
        "COMPRESS=JPEG", "QUALITY=90", "BLOCKSIZE=512", "OVERVIEW_RESAMPLING=AVERAGE",
        "NUM_THREADS=ALL_CPUS"])

    # the rim parallel as a polygon, 1440 vertices, closed, counter-clockwise
    lon = np.linspace(-180, 180, 1441)
    ring = [list(to_xy.TransformPoint(float(l), LAT_RIM)[:2]) for l in lon]
    ring[-1] = ring[0]
    if sum(x0 * y1 - x1 * y0 for (x0, y0), (x1, y1) in zip(ring[:-1], ring[1:])) < 0:
        ring = ring[::-1]
    ring = [[round(x, 2), round(y, 2)] for x, y in ring]
    neat = {"type": "FeatureCollection", "name": neat_name,
            "crs": {"type": "name", "properties": {"name": f"urn:ogc:def:crs:EPSG::{EPSG}"}},
            "features": [{"type": "Feature", "properties": {"latitude": LAT_RIM},
                          "geometry": {"type": "Polygon", "coordinates": [ring]}}]}
    with open(os.path.join(a.out, f"{neat_name}.geojson"), "w") as f:
        json.dump(neat, f)

    # streamed clip of the published COG to the disc
    url = f"/vsicurl/{a.url.rstrip('/')}/{a.name}.tif"
    neat_url = f"/vsicurl/{a.url.rstrip('/')}/{neat_name}.geojson"
    gdalg = {"type": "gdal_streamed_alg",
             # --like a published GeoJSON: an inline WKT geometry is split on
             # its negative coordinates when GDAL parses the command line
             "command_line": f"gdal raster pipeline read --input {url} ! clip --like "
                             f"{neat_url} --add-alpha",
             "gdal_version": "3130300"}
    with open(os.path.join(a.out, f"{a.name}.gdalg.json"), "w") as f:
        json.dump(gdalg, f, indent=1)

    # STAC item
    to_ll = osr.CoordinateTransformation(srs, geo)
    lats = [to_ll.TransformPoint(x, y)[1] for x, y in
            [(gt[0], gt[3]), (gt[0] + W * px, gt[3]), (gt[0], gt[3] - H * px), (gt[0] + W * px, gt[3] - H * px)]]
    item = {"type": "Feature", "stac_version": "1.1.0",
            "stac_extensions": ["https://stac-extensions.github.io/projection/v2.0.0/schema.json"],
            "id": a.name, "bbox": [-180, -90, 180, round(max(lats), 6)],
            "geometry": {"type": "Polygon", "coordinates": [[[-180, -90], [180, -90], [180, max(lats)],
                                                             [-180, max(lats)], [-180, -90]]]},
            "properties": {"title": a.title, "datetime": a.datetime, "license": a.license,
                           "proj:code": f"EPSG:{EPSG}", "proj:shape": [H, W],
                           "proj:transform": [gt[1], gt[2], gt[0], gt[4], gt[5], gt[3]],
                           "description": "Digital chart page, georeferenced from its printed scale "
                                          f"and the {abs(LAT_RIM):g} deg S rim; for display, not analysis."},
            "links": [{"rel": "cite-as", "href": f"https://doi.org/{a.doi}"}],
            "assets": {
                "visual": {"href": f"{a.url.rstrip('/')}/{a.name}.tif",
                           "type": "image/tiff; application=geotiff; profile=cloud-optimized",
                           "roles": ["visual", "overview"], "title": "Full page"},
                "visual-clipped": {"href": f"{a.url.rstrip('/')}/{a.name}.gdalg.json",
                                   "type": "application/json", "roles": ["visual"],
                                   "title": f"Chart clipped to {abs(LAT_RIM):g} deg S, alpha band (GDAL streamed algorithm)"},
                "neatline": {"href": f"{a.url.rstrip('/')}/{neat_name}.geojson",
                             "type": "application/geo+json", "roles": ["metadata"],
                             "title": f"{abs(LAT_RIM):g} deg S parallel, EPSG:9354"},
                "source": {"href": source_href,
                           "type": "application/pdf", "roles": ["source"], "file:checksum_sha256": h}}}
    with open(os.path.join(a.out, f"{a.name}.json"), "w") as f:
        json.dump(item, f, indent=1)
    prov = {"source": {"file": os.path.basename(a.pdf), "sha256": h, "doi": a.doi, "href": source_href},
            "dpi": a.dpi, "metres_per_pixel": px, "geotransform": gt,
            "pole_px": [cx, cy], "disc_fit_rms_px": rms, "disc_fit_rays": nfit,
            "rim_px": rim, "rim_sd_px": rim_sd, "rim_rays": nrim, "rim_vs_parallel": err, "rim_lat": LAT_RIM, "scale": SCALE,
            "cog": {"bytes": os.path.getsize(tif), "sha256": sha256(tif)}, "gdal": gdal.__version__}
    with open(os.path.join(a.out, f"{a.name}.provenance.json"), "w") as f:
        json.dump(prov, f, indent=1)
    print(f"{a.name}: {os.path.getsize(tif) / 1e6:.1f} MB", flush=True)


if __name__ == "__main__":
    main()
