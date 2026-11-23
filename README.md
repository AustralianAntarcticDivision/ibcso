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
