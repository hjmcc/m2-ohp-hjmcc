#!/usr/bin/env bash
# process_all.sh — Copy archive data and run full pipeline for 2018-2023
# Usage: bash scripts/process_all.sh [phase]
#   phase: setup, scan, astrometry, photometry, psf, all (default: all)
set -euo pipefail

ARCHIVE="/Users/hjmcc/Archive/ohp"
DATA_ROOT="/Users/hjmcc/Work/data"
PIPELINE_DIR="$(cd "$(dirname "$0")/.." && pwd)"
LOGDIR="$PIPELINE_DIR/scripts/logs"
mkdir -p "$LOGDIR"

phase="${1:-all}"

log() { echo "[$(date '+%H:%M:%S')] $*"; }

# ── Phase 1: Copy raw data into pipeline layout ─────────────────────────
setup_data() {
    log "=== Phase 1: Setting up raw data ==="

    # --- 2018: old format OHP_T{NNN}_{DATE} ---
    for tel in T080 T120; do
        dst="$DATA_ROOT/2018/$tel"
        if [ -d "$dst" ]; then
            log "  2018/$tel already exists, skipping"
            continue
        fi
        mkdir -p "$dst"
        for src in "$ARCHIVE/201801OHP_M2/DATA/OHP_${tel}_"*/; do
            night=$(basename "$src")
            log "  Copying $night → 2018/$tel/"
            cp -R "$src" "$dst/$night"
        done
    done

    # --- 2019: old format OHP_T{NNN}_{DATE} ---
    for tel in T080 T120; do
        dst="$DATA_ROOT/2019/$tel"
        if [ -d "$dst" ]; then
            log "  2019/$tel already exists, skipping"
            continue
        fi
        mkdir -p "$dst"
        for src in "$ARCHIVE/201901OHP_M2/DATA/OHP_${tel}_"*/; do
            night=$(basename "$src")
            log "  Copying $night → 2019/$tel/"
            cp -R "$src" "$dst/$night"
        done
    done

    # --- 2020: T120 only, non-standard night names + CALIB ---
    dst="$DATA_ROOT/2020/T120"
    if [ -d "$dst" ]; then
        log "  2020/T120 already exists, skipping"
    else
        mkdir -p "$dst"
        src="$ARCHIVE/202003OHP_M2/DATA/T120"
        # Copy night dirs (2020-03-*) and CALIB, skip calibration/ dir and stray files
        for d in "$src"/2020-03-*/; do
            night=$(basename "$d")
            log "  Copying $night → 2020/T120/"
            cp -R "$d" "$dst/$night"
        done
        if [ -d "$src/CALIB" ]; then
            log "  Copying CALIB → 2020/T120/"
            cp -R "$src/CALIB" "$dst/CALIB"
        fi
    fi

    # --- 2021: T080 and T120, skip processed/ and T152 ---
    for tel in T080 T120; do
        dst="$DATA_ROOT/2021/$tel"
        if [ -d "$dst" ]; then
            log "  2021/$tel already exists, skipping"
            continue
        fi
        src="$ARCHIVE/202103OHP_M2/DATA/$tel"
        if [ ! -d "$src" ]; then
            log "  2021/$tel not in archive, skipping"
            continue
        fi
        mkdir -p "$dst"
        for d in "$src"/20*/; do
            night=$(basename "$d")
            log "  Copying $night → 2021/$tel/"
            cp -R "$d" "$dst/$night"
        done
    done

    # --- 2022-2023: clean format, straight copy ---
    for year in 2022 2023; do
        archive_prefix="${year}03OHP_M2"
        for tel in T080 T120; do
            dst="$DATA_ROOT/$year/$tel"
            if [ -d "$dst" ]; then
                log "  $year/$tel already exists, skipping"
                continue
            fi
            src="$ARCHIVE/${archive_prefix}/DATA/$tel"
            if [ ! -d "$src" ]; then
                log "  $year/$tel not in archive, skipping"
                continue
            fi
            mkdir -p "$dst"
            for d in "$src"/20*/; do
                night=$(basename "$d")
                log "  Copying $night → $year/$tel/"
                cp -R "$d" "$dst/$night"
            done
        done
    done

    log "=== Phase 1 complete ==="
}

# ── Phase 2: Scan + Calibrate + QC (no network) ────────────────────────
run_scan() {
    log "=== Phase 2: Scan + Calibrate + QC ==="
    cd "$PIPELINE_DIR"
    for year in 2018 2019 2020 2021 2022 2023; do
        log "  Processing $year: scan"
        python -m pipeline scan --year "$year" 2>&1 | tee "$LOGDIR/${year}_scan.log"
        log "  Processing $year: calibrate"
        python -m pipeline calibrate --year "$year" 2>&1 | tee "$LOGDIR/${year}_calibrate.log"
        log "  Processing $year: qc"
        python -m pipeline qc --year "$year" 2>&1 | tee "$LOGDIR/${year}_qc.log"
    done
    log "=== Phase 2 complete ==="
}

# ── Phase 3: Astrometry (Gaia queries) ─────────────────────────────────
run_astrometry() {
    log "=== Phase 3: Astrometry ==="
    cd "$PIPELINE_DIR"
    for year in 2018 2019 2020 2021 2022 2023; do
        log "  Processing $year: astrometry"
        python -m pipeline astrometry --year "$year" 2>&1 | tee "$LOGDIR/${year}_astrometry.log"
        log "  Pausing 30s for rate limiting..."
        sleep 30
    done
    log "=== Phase 3 complete ==="
}

# ── Phase 4: Photometry (PS1 queries) ──────────────────────────────────
run_photometry() {
    log "=== Phase 4: Photometry ==="
    cd "$PIPELINE_DIR"
    for year in 2018 2019 2020 2021 2022 2023; do
        log "  Processing $year: photometry"
        python -m pipeline photometry --year "$year" 2>&1 | tee "$LOGDIR/${year}_photometry.log"
        log "  Pausing 30s for rate limiting..."
        sleep 30
    done
    log "=== Phase 4 complete ==="
}

# ── Phase 5: PSF models (local, includes 2024) ─────────────────────────
run_psf() {
    log "=== Phase 5: PSF ==="
    cd "$PIPELINE_DIR"
    for year in 2018 2019 2020 2021 2022 2023 2024; do
        log "  Processing $year: psf"
        python -m pipeline psf --year "$year" 2>&1 | tee "$LOGDIR/${year}_psf.log"
    done
    log "=== Phase 5 complete ==="
}

# ── Dispatch ────────────────────────────────────────────────────────────
case "$phase" in
    setup)      setup_data ;;
    scan)       run_scan ;;
    astrometry) run_astrometry ;;
    photometry) run_photometry ;;
    psf)        run_psf ;;
    all)
        setup_data
        run_scan
        run_astrometry
        run_photometry
        run_psf
        log "=== All phases complete ==="
        log "Run verification:"
        log "  for year in 2018 2019 2020 2021 2022 2023 2024 2025; do python -m pipeline status --year \$year; done"
        ;;
    *)
        echo "Usage: $0 [setup|scan|astrometry|photometry|psf|all]"
        exit 1
        ;;
esac
