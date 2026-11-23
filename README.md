# IBCSO on source.coop

Product repository for [`ausantarctic/ibcso`](https://source.coop/ausantarctic/ibcso):
Cloud Optimized GeoTIFF copies of the International Bathymetric Chart of the
Southern Ocean grids, display products made from the release charts, and the
virtual multidimensional indexes over them.

## Published layout

```
ibcso/
  README.md
  v1/cog/  v1/chart/  v1/background/      IBCSO v1.0 (2013)
  v2/cog/  v2/chart/                      IBCSO v2 (2022)
  annual-2024/cog/  annual-2025/cog/      archived annual releases
  virtual/                                STAC catalog and the indexes
```

Each `cog/` directory holds one COG per grid and a `provenance.json` with the
checksum of every source file and the result of a pixel comparison with it.

## Byte makers

Scripts that write the published bytes. Each needs GDAL (Python bindings) and
numpy; `ibcso_v1_cog.py` also needs netCDF4. Inputs are the release files as
fetched, with a `fetch.tsv` of their URL, size and SHA-256.

| script | writes | from |
|---|---|---|
| `bytes/ibcso_v1_cog.py` | `v1/cog/{bed,ice_surface,sid}.tif`, `provenance.json` | the v1.0 GMT grids (doi:10.1594/PANGAEA.805734), EPSG:3031 |
| `bytes/ibcso_background.py` | `v1/background/ibcso-v1-background-hq.tif`, `provenance.json` | `IBCSO_background_hq.zip` of the v1.0 release |
| `bytes/ibcso_chart.py` | `<release>/chart/<name>.tif`, `.gdalg.json`, `neatline-<lat>S.geojson`, a STAC item, `.provenance.json` | a release chart PDF, georeferenced from its printed scale and the fitted pole (EPSG:9354) |

The charts and the background are display products: they sit beside the
COGs and are not indexed.

```sh
python bytes/ibcso_v1_cog.py --raw RAW/unz --fetch RAW/fetch.tsv --out v1/cog
python bytes/ibcso_background.py --tif RAW/unz/IBCSO_background_hq/ibcso_background_hq.tif \
  --fetch RAW/fetch.tsv --out v1/background
python bytes/ibcso_chart.py --pdf RAW/IBSCO_v2_digital_chart.pdf --name ibcso-v2-chart \
  --title "..." --datetime 2022-05-20T12:03:02Z \
  --url https://data.source.coop/ausantarctic/ibcso/v2/chart --out v2/chart
```

The published `v1/`, `v1/background/` and `v2/chart/` objects of 2026-11-03
were written by these scripts.

## State of the move

This repository is the IBCSO product repo of the per-product repo map
(product files, dataset config, byte makers and runner here; source lists
from the [aad-filelist](https://github.com/mdsumner/aad-filelist) corpus;
engine from [raad.virtualize](https://github.com/mdsumner/raad.virtualize)
at a tag). The move is staged. Still in raad.virtualize, until it has a tag
to pin:

| piece | path in raad.virtualize |
|---|---|
| product files | `products/ibcso/` (`product.json`, `datasets.tsv`, `README.md`) |
| dataset config | `config/ibcso.sh` |
| COGs of v2 and the annual releases | `py/ibcso_cog.py` |
| runner | `hpc/vm/ibcso-publish.sh` |

## Citation and terms

IBCSO products are not for navigation or any other purpose involving safety
at sea. See the published README for the citation of each release.
