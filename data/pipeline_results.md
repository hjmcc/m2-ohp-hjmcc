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

WCS solutions via Gaia DR3 cross-matching with two-pass solver.

| Year | Telescope | Solved | Rate | Median RMS | Notes |
|------|-----------|--------|------|-----------|-------|
| 2018 | T080      | 452/452 | 100% | 1.9"     | Single field (HAT-P-32) |
| 2018 | T120      | 34/34   | 100% | 2.7"     | |
| 2019 | T080      | 1,012/1,231 | 82% | 2.1" | 13 failures, 243 HIGH_RMS |
| 2019 | T120      | 85/87   | 98%  | 2.9"     | |
| 2020 | T120      | 30/583  | 5%   | 3.8"     | Most objects unresolvable |
| 2021 | T120      | 447/564 | 79%  | 5.0"     | 117 failures |
| 2023 | T080      | 0/9     | 0%   | —        | All skipped (no WCS input) |
| 2023 | T120      | 17/504  | 3%   | 5.7"     | 191 failures, 296 skipped |
| 2024 | T080      | 36/83   | 43%  | 2.5"     | |
| 2024 | T120      | 0/763   | 0%   | 12.3"    | All failed |
| 2025 | T080      | 112/114 | 98%  | 2.0"     | |
| 2025 | T120      | 308/359 | 86%  | 3.9"     | |

**Low solve rates in 2020, 2023, 2024 T120** are caused by student-chosen
object names that Simbad cannot resolve to sky coordinates (e.g. "McCracken",
"saMar30", "SCAX_new", "2023_EO"). Without an initial pointing estimate the
solver cannot query the correct Gaia field.

## Photometric Zero Points

Derived from PS1 DR2 cross-matching. Only frames with PS1-compatible filters
(g/r/i) AND successful astrometric solutions produce zero points.

Convention: `mag = ZP - 2.5 * log10(counts / exptime)`

### T080

| Year | Filter | ZP (mag) | Scatter (MAD) | N frames |
|------|--------|----------|---------------|----------|
| 2019 | g'     | 21.9     | 0.4           | 136      |
| 2019 | r'     | 22.0     | 0.6           | 140      |
| 2019 | i'     | 23.1     | 0.2           | 2        |
| 2023 | g'     | 22.3     | 0.8           | 3        |
| 2023 | r'     | 24.6     | 0.7           | 3        |
| 2023 | i'     | 22.2     | 0.1           | 3        |
| 2024 | g'     | 21.5     | 2.9           | 4        |
| 2024 | r'     | 20.6     | 3.2           | 6        |
| 2024 | i'     | 21.2     | 1.1           | 7        |
| 2025 | g'     | 20.6     | 1.0           | 38       |
| 2025 | r'     | 21.4     | 1.2           | 35       |
| 2025 | i'     | 21.2     | 1.3           | 24       |

### T120

| Year | Filter | ZP (mag) | Scatter (MAD) | N frames |
|------|--------|----------|---------------|----------|
| 2023 | g      | 23.2     | 1.1           | 41       |
| 2023 | r      | 23.7     | 1.1           | 91       |
| 2024 | V      | 25.1     | 0.2           | 5        |
| 2024 | R      | 24.2     | 1.4           | 198      |
| 2024 | g      | 23.9     | 1.2           | 93       |
| 2025 | B      | 25.3     | 1.7           | 37       |
| 2025 | V      | 25.9     | 1.4           | 40       |
| 2025 | R      | 25.0     | 1.5           | 112      |
| 2025 | g      | 24.7     | 1.9           | 24       |

### Years with no zero points

- **2018**: T080 used Clear filter (no PS1 match); T120 used B/R/V/Halpha.
- **2020**: T120 used B/R/V/g/Halpha/i, but only 30/583 frames had WCS
  solutions and those were in non-PS1 bands.
- **2021**: T120 used B/R/V/g/Halpha/i, but no photometry matches despite
  having WCS for 447 frames.
- **2022**: No reduced frames at all.

The ~1–2 mag year-to-year variation is consistent with non-photometric
conditions. The 2023 T080 r' ZP = 24.6 is anomalous (3 frames only). These
values are suitable for order-of-magnitude flux estimates only.

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

1. **2022 binning mismatch** — All science frames at 1x1 binning. Requires
   adding a 1x1 expected shape to config.py or separate processing.
2. **French directory names (2019 T120)** — "Deuxième nuit", "Première nuit",
   "Troisième nuit" break SExtractor. Fix: rename directories or quote paths
   in the SExtractor subprocess call.
3. **Unresolvable object names** — Student-chosen names in 2020, 2023, 2024
   cannot be resolved by Simbad. Fix: add a manual coordinate lookup table
   for known student targets.
4. **2021 T080** — Raw frames exist but were all skipped during reduction
   (shape mismatch). Needs investigation.
5. **2024 T120 astrometry** — 0% solve rate despite 763 reduced frames.
   Median RMS 12.3" suggests systematic pointing error. Needs investigation.
