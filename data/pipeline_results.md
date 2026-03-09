# Pipeline Processing Results: 2018–2025

Full pipeline (scan, calibrate, QC, astrometry, photometry, PSF) run across
all 8 years of OHP M2 observing data. Years 2024–2025 were processed
previously; 2018–2023 were batch-processed in March 2026.

## Data Inventory

| Year | T080 Raw | T120 Raw | T080 Reduced | T120 Reduced | Notes |
|------|----------|----------|-------------|-------------|-------|
| 2018 | 469      | 144      | 452         | 34          | T080: no bias, flat-only reduction |
| 2019 | 895      | 87       | 1,231       | 87          | T080: 5 nights of transit photometry |
| 2020 | —        | 586      | —           | 583         | T120 only |
| 2021 | 148      | 564      | 0           | 564         | T080: all wrong shape (skipped) |
| 2022 | 90       | 117      | 0           | 0           | All 1x1 binned (see below) |
| 2023 | 118      | 923      | 9           | 504         | T080: only r-band flats survived QC |
| 2024 | 83       | 763      | 83          | 763         | |
| 2025 | 114      | 359      | 114         | 359         | |

**Total reduced frames: ~4,872** across all years and telescopes.

## Calibration Notes

- **2018 T080**: No bias frames in archive. Dome flats available for Clear and
  r bands only; 4/8 Clear and 4/7 r flats survived QC (low-signal rejects).
  Science frames flat-corrected without bias subtraction.
- **2019 T080**: 120 bias frames, 210 flats (Clear/g/r). Largest single-year
  dataset — 1,231 science frames across 5 nights of exoplanet transit
  monitoring (HAT-P-32b, WASP-33b, WASP-104b, HAT-P-10, HAT-P-36).
- **2019 T120**: French directory names in archive ("Deuxième nuit",
  "Première nuit") cause SExtractor failures on 36/87 frames due to accented
  characters in paths.
- **2020 T120**: Good calibration (30 bias, ~145 flats across 6 filters).
  Non-standard night directory names (2020-03-02-M2AAIS format).
- **2021 T080**: 148 raw frames but all skipped — shape mismatch with expected
  T080 dimensions. No science output.
- **2022**: Complete loss. All science frames taken at **1x1 binning**
  (2048x2048) instead of expected 2x2 (1024x1024). Master bias/flats were
  built but no science frames could be reduced. Would require separate 1x1
  calibration pipeline configuration.
- **2023 T080**: Only 9 science frames reduced. g-band flats rejected (low
  signal ~6100 ADU); only r-band flats survived. 6 frames got bias-only
  reduction (no matching flat).
- **2023 T120**: V-band flats all saturated (~60,000–61,000 ADU). Many r-band
  "flats" were actually dark frames (~180 ADU). 504/509 frames reduced
  despite challenging calibration data.

## Astrometry

T120 WCS solutions via `astrometry.net` `solve-field` (blind plate solving)
with Gaia DR3 cross-match validation.  T080 still uses the original pipeline
Gaia cross-matching solver (pending reprocessing — needs smaller-scale index
files for its 12'×8' FOV).

### T120 (solve-field)

| Year | Solved | Rate  | Median RMS | Notes |
|------|--------|-------|-----------|-------|
| 2018 | 33/34  | 97%   | 0.37"     | 1 failure (IC10) |
| 2019 | 87/87  | 100%  | 0.36"     | |
| 2020 | 583/583 | 100% | 0.35"     | |
| 2021 | 549/564 | 97%  | 0.42"     | 15 failures |
| 2023 | 500/504 | 99%  | 0.38"     | |
| 2024 | 753/763 | 99%  | 0.35"     | |
| 2025 | 343/359 | 95%  | 0.45"     | 15 M49/NGC5813 failures (sparse field) |

**Total: 2,848/2,894 solved (98.4%), typical RMS 0.3–0.5".**

### T080 (legacy pipeline solver — treat with caution)

| Year | Solved | Rate  | Median RMS | Notes |
|------|--------|-------|-----------|-------|
| 2018 | 452/452 | 100% | 1.9"     | Single field (HAT-P-32) |
| 2019 | 1,012/1,231 | 82% | 2.1" | 13 failures, 243 HIGH_RMS |
| 2023 | 0/9     | 0%   | —        | All skipped (no WCS input) |
| 2024 | 36/83   | 43%  | 2.5"     | |
| 2025 | 112/114 | 98%  | 2.0"     | |

T080 astrometry has known issues (RMS ~2" vs <0.5" for T120 solve-field).
Pending reprocessing with solve-field once index-4203 tiles are complete.

## Photometric Zero Points

Derived from PS1 DR2 cross-matching.  Values are 3σ-clipped medians;
scatter σ reflects frame-to-frame variation (weather, airmass,
flat-fielding).  No colour terms applied — B and V are calibrated against
PS1 g-band.

Convention: `mag = ZP - 2.5 * log10(counts / exptime)`

### T120 (solve-field astrometry)

| Year | Band | ZP (mag) | σ (mag) | N frames | Raw filter names |
|------|------|----------|---------|----------|-----------------|
| 2018 | B    | 21.23    | 0.14    | 5        | B_Cousins |
| 2018 | R    | 20.76    | 1.79    | 26       | R_Cousins |
| 2019 | B    | 23.14    | 0.02    | 9        | B Cousins |
| 2019 | R    | 24.19    | 0.02    | 24       | R Cousins, R_Cousins |
| 2020 | V    | 23.23    | 2.08    | 10       | v_Gunn |
| 2020 | g    | 22.76    | 0.08    | 70       | g_Gunn |
| 2020 | R    | 23.45    | 0.03    | 335      | r_Gunn |
| 2020 | i    | 23.41    | 0.03    | 12       | i_Gunn |
| 2021 | B    | 22.47    | 0.25    | 74       | B_Cousins |
| 2021 | V    | 23.76    | 0.20    | 74       | V_Cousins |
| 2021 | g    | 22.21    | 0.24    | 59       | g_Gunn |
| 2021 | R    | 22.93    | 0.47    | 264      | r_Gunn |
| 2021 | i    | 23.10    | 0.09    | 28       | i_Gunn |
| 2023 | B    | 23.41    | 0.06    | 60       | B_cousins |
| 2023 | V    | 24.46    | 0.03    | 31       | V_cousins |
| 2023 | g    | 23.25    | 0.05    | 43       | g |
| 2023 | R    | 23.72    | 0.10    | 182      | r |
| 2024 | V    | 24.30    | 0.01    | 5        | V |
| 2024 | g    | 23.10    | 0.05    | 127      | G |
| 2024 | R    | 24.25    | 0.05    | 331      | R |
| 2025 | B    | 23.15    | 0.49    | 47       | B |
| 2025 | V    | 24.12    | 0.30    | 45       | V |
| 2025 | g    | 22.53    | 0.03    | 23       | G |
| 2025 | R    | 24.10    | 0.07    | 144      | R |

**Total: 2,349/2,535 T120 frames calibrated (93%).**  Uncalibrated frames are
narrowband (Hα, SII, OIII) which have no PS1 equivalent.

### T080 (legacy astrometry — treat with caution)

| Year | Band | ZP (mag) | σ (mag) | N frames |
|------|------|----------|---------|----------|
| 2019 | g    | 21.90    | 0.64    | 135      |
| 2019 | R    | 22.05    | 0.87    | 138      |
| 2019 | i    | 23.08    | 0.17    | 2        |
| 2023 | g    | 22.32    | 0.93    | 3        |
| 2023 | R    | 24.56    | 0.97    | 3        |
| 2023 | i    | 22.19    | 0.32    | 3        |
| 2024 | g    | 21.51    | 2.52    | 4        |
| 2024 | R    | 19.98    | 0.79    | 4        |
| 2024 | i    | 21.18    | 1.03    | 7        |
| 2025 | g    | 20.66    | 0.85    | 37       |
| 2025 | R    | 21.35    | 1.19    | 35       |
| 2025 | i    | 21.20    | 1.25    | 24       |

T080 ZPs have large scatter partly due to unreliable WCS (pending
solve-field reprocessing).

### Notes

- **2018 T120 is anomalous**: R-band ZP=20.8 with σ=1.8, ~3 mag fainter than
  other years.  Only 34 frames; likely a different gain setting.
- **2022**: No reduced frames (all 1×1 binned).
- Year-to-year T120 R-band ZP variation (~1 mag) likely reflects real
  instrumental changes (mirror recoating, CCD replacement) plus the mix of
  photometric and non-photometric conditions.

## PSF Models

Per-frame PSF models built via SExtractor + PSFex (degree-2 spatial variation).

| Year | Telescope | Models Built | Median FWHM | Median Ellipticity |
|------|-----------|-------------|-------------|-------------------|
| 2018 | T080      | 451/452     | 4.2"        | 0.081             |
| 2018 | T120      | 32/34       | 2.5"        | 0.067             |
| 2019 | T080      | 552/1,231   | 0.9"        | 0.021             |
| 2019 | T120      | 51/87       | 3.1"        | 0.067             |
| 2020 | T120      | 581/583     | 3.8"        | 0.048             |
| 2021 | T120      | 523/564     | 2.2"        | 0.067             |
| 2023 | T080      | 6/9         | —           | —                 |
| 2023 | T120      | 491/504     | 3.4"        | 0.043             |
| 2024 | T080      | 68/83       | 0.9"        | 0.062             |
| 2024 | T120      | 717/763     | 2.9"        | 0.069             |
| 2025 | T080      | 99/114      | 3.7"        | —                 |
| 2025 | T120      | 359/359     | 3.2"        | —                 |

The low 2019 T080 FWHM (0.9") and build rate (45%) reflect short-exposure
transit frames with few bright, isolated stars — PSFex needs >= 20 clean
sources for a reliable model. The 2019 T120 failures (36/87) are caused by
SExtractor's inability to handle accented characters in file paths.

## Known Issues

1. **2022 binning mismatch** — All science frames at 1×1 binning. Requires
   adding a 1×1 expected shape to config.py or separate processing.
2. **French directory names (2019 T120)** — "Deuxième nuit", "Première nuit",
   "Troisième nuit" break SExtractor. Fix: rename directories or quote paths
   in the SExtractor subprocess call.
3. **2021 T080** — Raw frames exist but were all skipped during reduction
   (shape mismatch). Needs investigation.
4. **T080 astrometry** — Still using legacy Gaia cross-matching solver with
   ~2" RMS. Needs reprocessing with solve-field once index-4203 tiles
   (5.5–8' features) are downloaded for all healpix positions.
5. **No colour terms** — B and V ZPs are calibrated against PS1 g-band
   without colour correction. SCAMP refinement planned as next step.
