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

## Year-to-Year Comparison

| Telescope | Filter | 2024 ZP | 2025 ZP | Diff |
|-----------|--------|---------|---------|------|
| T080      | g'     | 21.5    | 20.6    | −0.9 |
| T080      | r'     | 20.6    | 21.4    | +0.8 |
| T080      | i'     | 21.2    | 21.2    | 0.0  |
| T120      | V      | 25.1    | 25.9    | +0.8 |
| T120      | R      | 24.2    | 25.0    | +0.8 |
| T120      | g      | 23.9    | 24.7    | +0.8 |

The ~1 mag year-to-year variation is consistent with the large per-frame scatter
and non-photometric conditions. The T080 2024 results are based on very few
frames (4–7) and should be treated with extra caution. T120 shows a systematic
~0.8 mag offset between years, likely due to different atmospheric conditions.

ZP values are medians across all calibrated frames. Use these for
order-of-magnitude flux estimates; do not rely on them for precision photometry.
