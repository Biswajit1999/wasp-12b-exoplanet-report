# Changelog

## 2.1.0 — 2026-10-04

- Added a new scientifically informed artistic concept of WASP-12 b and its late-F host.
- Labelled the image in the website and README as conceptual, non-observational and not to scale.
- Documented which visual cues come from the orbital, tidal-deformation and phase-curve literature and which details remain interpretive.

## 2.0.0 — 2026-10-04

- Completed the homogeneous 120-second TESS inventory with Sectors 44, 71 and 72, bringing the independent photometric check to six sectors and 119 supported event timings.
- Added the verbatim 379-row Leonardi et al. (2024) VizieR transit-timing compilation and fail-closed SHA-256 validation for every scientific input.
- Reproduced Ṗ = −31.29 ms yr⁻¹ and ΔBIC = 1687 for quadratic versus linear ephemerides; exposed both the ±0.76 ms yr⁻¹ formal error and the ±1.09 ms yr⁻¹ scatter-scaled error.
- Added leave-one-literature-source-out refits; all major-source deletions retain Ṗ between −31.93 and −30.35 ms yr⁻¹.
- Added a bounded opposing-sign sinusoidal timing proxy on the Yee et al. transit-and-occultation compilation. The quadratic decay model is preferred by ΔBIC = 20.16.
- Rebuilt the hosted site as a time-domain scientific dossier and removed synthetic artwork, animation, generic dashboard cards and portfolio link clutter.
- Added 12 regression/integrity/accessibility tests and pinned CI actions to immutable commit SHAs.

## 1.1.0 — 2026-08-16

- Added the initial three-sector TESS and 158-event timing reproduction.

## 1.0.0 — 2026-08-15

- Initial timing-adjusted TESS transit report.
