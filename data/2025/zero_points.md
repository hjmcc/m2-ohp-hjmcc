# Photometric Zero Points — OHP M2

Zero points derived from PS1 DR2 cross-matching of extracted sources in all
reduced science frames.  Astrometric solutions from `astrometry.net`
`solve-field`; Gaia DR3 validation gives typical RMS < 0.5″.

Convention: `mag = ZP - 2.5 * log10(counts / exptime)`

where `counts` is the total aperture flux in ADU and `exptime` is in seconds.
Values are sigma-clipped (3σ) medians.  The scatter σ reflects frame-to-frame
variation (weather, airmass, flat-fielding); it is **not** the uncertainty on
the median.

## T120 (1.2 m Newton)

### Per-Year Zero Points

| Filter | PS1 band | 2018 | 2019 | 2020 | 2021 | 2023 | 2024 | 2025 |
|--------|----------|------|------|------|------|------|------|------|
| B      | gmag     | 21.23 (5) | 23.14 (9) | — | 22.47 (74) | 23.41 (60) | — | 23.15 (47) |
| V      | gmag     | — | — | 23.23 (10) | 23.76 (74) | 24.46 (31) | 24.30 (5) | 24.12 (45) |
| g      | gmag     | — | — | 22.76 (70) | 22.21 (59) | 23.25 (43) | 23.10 (127) | 22.53 (23) |
| R      | rmag     | 20.76 (26) | 24.19 (24) | 23.45 (335) | 22.93 (264) | 23.72 (182) | 24.25 (331) | 24.10 (144) |
| i      | imag     | — | — | 23.41 (12) | 23.10 (28) | — | — | — |

Number of frames used in parentheses.

### Per-Year Scatter (σ, mag)

| Filter | 2018 | 2019 | 2020 | 2021 | 2023 | 2024 | 2025 |
|--------|------|------|------|------|------|------|------|
| B      | 0.14 | 0.02 | — | 0.25 | 0.06 | — | 0.49 |
| V      | — | — | 2.08 | 0.20 | 0.03 | 0.01 | 0.30 |
| g      | — | — | 0.08 | 0.24 | 0.05 | 0.05 | 0.03 |
| R      | 1.79 | 0.02 | 0.03 | 0.47 | 0.10 | 0.05 | 0.07 |
| i      | — | — | 0.03 | 0.09 | — | — | — |

### Filter Names by Year

Different filter designations have been used at the T120 across runs.  All
variants within a row are mapped to the same PS1 band.

| Band | 2018 | 2019 | 2020 | 2021 | 2023 | 2024 | 2025 |
|------|------|------|------|------|------|------|------|
| B | B_Cousins | B Cousins | — | B_Cousins | B_cousins | — | B |
| V | — | — | v_Gunn | V_Cousins | V_cousins | V | V |
| g | — | — | g_Gunn | g_Gunn | g | G | G |
| R | R_Cousins | R Cousins, R_Cousins | r_Gunn | r_Gunn | r | R | R |
| i | — | — | i_Gunn | i_Gunn | — | — | — |


## T080 (80 cm Cassegrain)

T080 astrometry not yet re-processed with `solve-field` (needs smaller-scale
index files for the 12′×8′ FOV).  These ZPs use the original pipeline
astrometry which has known issues; treat with caution.

| Filter | PS1 band | 2019 | 2023 | 2024 | 2025 |
|--------|----------|------|------|------|------|
| g (g') | gmag     | 21.90 (135) | 22.32 (3) | 21.51 (4) | 20.66 (37) |
| R (r') | rmag     | 22.05 (138) | 24.56 (3) | 19.98 (4) | 21.35 (35) |
| i (i') | imag     | 23.08 (2) | 22.19 (3) | 21.18 (7) | 21.20 (24) |

| Filter | σ 2019 | σ 2023 | σ 2024 | σ 2025 |
|--------|--------|--------|--------|--------|
| g      | 0.64   | 0.93   | 2.52   | 0.85   |
| R      | 0.87   | 0.97   | 0.79   | 1.19   |
| i      | 0.17   | 0.32   | 1.03   | 1.25   |

## Notes

- **2018 T120 is anomalous**: the R-band ZP=20.8 with σ=1.8 is ~3 mag fainter
  than other years.  Only 34 frames total; likely a different gain setting or
  pre-processing issue.
- **2022 has no data**: all frames were 1×1 binned, mismatching the expected
  2×2 binning, so calibration was skipped.
- **No colour terms** are applied.  B and V are calibrated against PS1 g-band
  without correction for the colour difference; this contributes to scatter in
  those bands.
- **T120 year-to-year variation** of ~0.5–1 mag in R-band likely reflects real
  instrumental changes (mirror recoating, CCD changes) plus the mix of
  photometric and non-photometric conditions.
- **T080 ZPs are unreliable** pending astrometry reprocessing.  The large
  scatter (σ ~ 1 mag) partly reflects incorrect cross-matches from poor WCS.
- Narrowband filters (Hα, SII, OIII, 6648_MD) have no PS1 equivalent and are
  excluded.
- **Astrometry success rate**: 98.8% across 2535 T120 frames (2018–2025),
  with sub-arcsecond RMS (typically 0.3–0.7″).
