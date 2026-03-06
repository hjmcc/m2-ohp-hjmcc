# Approximate Photometric Zero Points — OHP M2

Zero points derived from PS1 DR2 cross-matching of all science frames.
Conditions were **not photometric** on most nights, so these are rough
estimates only (typical scatter ~1–2 mag frame-to-frame).

Convention: `mag = ZP - 2.5 * log10(counts / exptime)`

where `counts` is the total aperture flux in ADU and `exptime` is in seconds.

## 2025 Run

### T080 (80 cm Cassegrain)

| Filter | ZP (mag) | Scatter (mag) | N frames |
|--------|----------|---------------|----------|
| g'     | 20.6     | 1.0           | 38       |
| r'     | 21.4     | 1.2           | 35       |
| i'     | 21.2     | 1.3           | 24       |

### T120 (1.2 m Newton)

| Filter | ZP (mag) | Scatter (mag) | N frames |
|--------|----------|---------------|----------|
| B      | 25.3     | 1.7           | 37       |
| V      | 25.9     | 1.4           | 40       |
| R      | 25.0     | 1.5           | 112      |
| g      | 24.7     | 1.9           | 24       |

## 2023 Run

### T080 (80 cm Cassegrain)

| Filter | ZP (mag) | Scatter (mag) | N frames |
|--------|----------|---------------|----------|
| g'     | 22.3     | 0.8           | 3        |
| r'     | 24.6     | 0.7           | 3        |
| i'     | 22.2     | 0.1           | 3        |

### T120 (1.2 m Newton)

| Filter | ZP (mag) | Scatter (mag) | N frames |
|--------|----------|---------------|----------|
| g      | 23.2     | 1.1           | 41       |
| r      | 23.7     | 1.1           | 91       |

## 2024 Run

### T080 (80 cm Cassegrain)

| Filter | ZP (mag) | Scatter (mag) | N frames |
|--------|----------|---------------|----------|
| g'     | 21.5     | 2.9           | 4        |
| r'     | 20.6     | 3.2           | 6        |
| i'     | 21.2     | 1.1           | 7        |

### T120 (1.2 m Newton)

| Filter | ZP (mag) | Scatter (mag) | N frames |
|--------|----------|---------------|----------|
| V      | 25.1     | 0.2           | 5        |
| R      | 24.2     | 1.4           | 198      |
| g      | 23.9     | 1.2           | 93       |

## 2019 Run

### T080 (80 cm Cassegrain)

| Filter | ZP (mag) | Scatter (mag) | N frames |
|--------|----------|---------------|----------|
| g'     | 21.9     | 0.4           | 136      |
| r'     | 22.0     | 0.6           | 140      |
| i'     | 23.1     | 0.2           | 2        |

## Year-to-Year Comparison

### T080

| Filter | 2019   | 2023   | 2024   | 2025   |
|--------|--------|--------|--------|--------|
| g'     | 21.9   | 22.3   | 21.5   | 20.6   |
| r'     | 22.0   | 24.6   | 20.6   | 21.4   |
| i'     | 23.1   | 22.2   | 21.2   | 21.2   |

### T120

| Filter | 2023   | 2024   | 2025   |
|--------|--------|--------|--------|
| g      | 23.2   | 23.9   | 24.7   |
| r      | 23.7   | —      | —      |
| R      | —      | 24.2   | 25.0   |
| V      | —      | 25.1   | 25.9   |
| B      | —      | —      | 25.3   |

## Notes

- Years 2018, 2020, 2021 produced no photometric zero points: 2018 T080 used
  Clear filter (no PS1 match), T120 used B/R/V/Halpha; 2020 and 2021 T120
  used B/R/V/Halpha/OIII with few astrometric solutions.
- 2022 had no reduced frames (all 1x1 binned, mismatched expected 2x2 shape).
- The 2023 T080 r' ZP=24.6 is anomalous (only 3 frames) and likely unreliable.
- The ~1–2 mag year-to-year variation is consistent with the large per-frame
  scatter and non-photometric conditions.
- ZP values are medians across all calibrated frames. Use these for
  order-of-magnitude flux estimates; do not rely on them for precision photometry.
