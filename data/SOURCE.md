# Data sources

## TESS light curve

- File: `tess2019357164649-s0020-0000000086396382-0165-s_lc.fits`
- Archive: Mikulski Archive for Space Telescopes (MAST), TESS SPOC light-curve product
- TESS sector: 20
- TIC target ID: 86396382
- MAST observation ID: 27242856
- MAST data URI: `mast:TESS/product/tess2019357164649-s0020-0000000086396382-0165-s_lc.fits`
- Exact download URL: <https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:TESS%2Fproduct%2Ftess2019357164649-s0020-0000000086396382-0165-s_lc.fits>
- Collection DOI: [10.17909/t9-nmc8-f686](https://doi.org/10.17909/t9-nmc8-f686) (TESS 2-minute light curves, all sectors; sector 20 used here)
- Retrieved: 2026-08-15
- SHA-256: `497b2e27ffc8f6a2f1c76b388d7a5f4943ecd4e48b89c3882e7d014be91857b7`

The FITS file is stored unmodified. The analysis reads `TIME`, `PDCSAP_FLUX`,
`PDCSAP_FLUX_ERR`, and `QUALITY`. PDCSAP flux is the SPOC light curve with common
instrumental trends removed and aperture/crowding corrections applied; this does
not make it free of residual stellar or instrumental systematics.

## System parameters

- File: `system_parameters.csv`
- Service: NASA Exoplanet Archive TAP, `pscomppars` table
- Exact query: <https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query=select+pl_name%2Chostname%2Cra%2Cdec%2Cpl_orbper%2Cpl_tranmid%2Cpl_trandur%2Cpl_rade%2Cpl_bmasse%2Cpl_eqt%2Cpl_orbsmax%2Csy_dist%2Csy_tmag%2Cst_teff%2Cst_rad%2Cst_mass%2Cdisc_year%2Cdiscoverymethod%2Cdisc_refname%2Cdisc_pubdate%2Cdisc_facility+from+pscomppars+where+pl_name%3D%27WASP-12+b%27&format=csv>
- Retrieved: 2026-08-15

The saved row is the input actually used by `scripts/analyze_transit.py`; the
analysis does not query a changing live service at run time.


## Additional TESS sectors for robustness analysis

All are unmodified standard-cadence SPOC light curves from the same [MAST TESS collection](https://doi.org/10.17909/t9-nmc8-f686).

- Sector 20: `tess2019357164649-s0020-0000000086396382-0165-s_lc.fits` (1,926,720 bytes)
  - MAST URI: `mast:TESS/product/tess2019357164649-s0020-0000000086396382-0165-s_lc.fits`
  - SHA-256: `497b2e27ffc8f6a2f1c76b388d7a5f4943ecd4e48b89c3882e7d014be91857b7`
- Sector 45: `tess2021310001228-s0045-0000000086396382-0216-s_lc.fits` (1,840,320 bytes)
  - MAST URI: `mast:TESS/product/tess2021310001228-s0045-0000000086396382-0216-s_lc.fits`
  - SHA-256: `1ac919f1a4b472b7d6c1911901e4f967e222cec8559204c58c3b091d995f5d06`
- Sector 43: `tess2021258175143-s0043-0000000086396382-0214-s_lc.fits` (1,811,520 bytes)
  - MAST URI: `mast:TESS/product/tess2021258175143-s0043-0000000086396382-0214-s_lc.fits`
  - SHA-256: `e3d5c6f76884014ec1dec1eb60f04eaf734330a7abd8eef99141f49df03b2493`
