"""CLI entry point: python -m pipeline <command> [options]."""

import argparse
import json
import sys
from pathlib import Path

from . import __version__, config
from .utils import setup_logging, load_archive, output_dir


def cmd_qc(args):
    from .quality import run_science_qc

    telescopes = [args.telescope] if args.telescope else list(config.TELESCOPES)

    for tel in telescopes:
        if tel not in config.TELESCOPES:
            print(f"Unknown telescope: {tel}")
            sys.exit(1)
        run_science_qc(args.year, tel, force=args.force)


def cmd_astrometry(args):
    from .astrometry import run_astrometry

    telescopes = [args.telescope] if args.telescope else list(config.TELESCOPES)

    for tel in telescopes:
        if tel not in config.TELESCOPES:
            print(f"Unknown telescope: {tel}")
            sys.exit(1)
        run_astrometry(args.year, tel, force=args.force)


def cmd_photometry(args):
    from .photometry import run_photometry

    telescopes = [args.telescope] if args.telescope else list(config.TELESCOPES)

    for tel in telescopes:
        if tel not in config.TELESCOPES:
            print(f"Unknown telescope: {tel}")
            sys.exit(1)
        run_photometry(args.year, tel, force=args.force)


def cmd_psf(args):
    from .psf import run_psf

    telescopes = [args.telescope] if args.telescope else list(config.TELESCOPES)

    for tel in telescopes:
        if tel not in config.TELESCOPES:
            print(f"Unknown telescope: {tel}")
            sys.exit(1)
        run_psf(args.year, tel, force=args.force)


def _load_archive_for(year, telescope):
    """Load archive JSON — prefer local archive.json, fall back to global."""
    if telescope:
        local = config.DATA_ROOT / year / telescope / "archive.json"
        if local.exists():
            with open(local) as f:
                return json.load(f)
    # Fall back to global archive
    return load_archive()


def cmd_scan(args):
    from .scanner import run_scan

    telescopes = [args.telescope] if args.telescope else list(config.TELESCOPES)

    for tel in telescopes:
        if tel not in config.TELESCOPES:
            print(f"Unknown telescope: {tel}")
            sys.exit(1)
        result = run_scan(args.year, tel)
        if result:
            print(f"  Wrote {result}")


def cmd_calibrate(args):
    from .calibration import run_calibration

    telescopes = [args.telescope] if args.telescope else list(config.TELESCOPES)

    for tel in telescopes:
        if tel not in config.TELESCOPES:
            print(f"Unknown telescope: {tel}")
            sys.exit(1)
        archive = _load_archive_for(args.year, tel)
        run_calibration(archive, args.year, tel, force=args.force)


def cmd_status(args):
    import statistics

    telescopes = [args.telescope] if args.telescope else list(config.TELESCOPES)

    for tel in telescopes:
        base = config.DATA_ROOT / args.year / tel
        print(f"\n{'=' * 50}")
        print(f"  {args.year} / {tel}")
        print(f"{'=' * 50}")

        for subdir in ["master_bias", "master_flat", "reduced",
                        "astrometry", "photometry", "psf", "qc"]:
            d = base / subdir
            if d.exists():
                fits_files = list(d.rglob("*.fits")) + list(d.rglob("*.fit"))
                json_files = list(d.rglob("*.json"))
                parts = []
                if fits_files:
                    parts.append(f"{len(fits_files)} FITS")
                if json_files:
                    parts.append(f"{len(json_files)} JSON")
                print(f"  {subdir + '/':20s} {', '.join(parts) if parts else 'empty'}")
            else:
                print(f"  {subdir + '/':20s} —")

        # Show calibration QC summary
        qc_file = base / "qc" / "qc_report.json"
        if qc_file.exists():
            with open(qc_file) as f:
                qc = json.load(f)
            print(f"\n  Calibration:")
            print(f"    Bias: {qc['bias']['n_accepted']}/{qc['bias']['n_input']} "
                  f"accepted (qual={qc['bias']['quality']})")
            for filt, fq in qc.get("flats", {}).items():
                print(f"    Flat [{filt}]: {fq['n_accepted']}/{fq['n_input']} accepted")
            sci = qc.get("science", {})
            print(f"    Science: {sci.get('n_reduced', 0)}/{sci.get('n_total', 0)} "
                  f"reduced ({sci.get('n_bias_only', 0)} bias-only, "
                  f"{sci.get('n_skipped', 0)} skipped)")

        # Show frame-level QC summary
        frame_stats_file = base / "qc" / "frame_stats.json"
        if frame_stats_file.exists():
            with open(frame_stats_file) as f:
                frame_stats = json.load(f)
            n_total = len(frame_stats)
            fwhm_vals = [s["fwhm_arcsec"] for s in frame_stats
                         if s.get("fwhm_arcsec") is not None]
            flags = {}
            for s in frame_stats:
                for flag in s.get("flag", "").split(","):
                    flag = flag.strip()
                    if flag and flag != "OK":
                        flags[flag] = flags.get(flag, 0) + 1
            n_ok = sum(1 for s in frame_stats if s.get("flag") == "OK")
            print(f"\n  Frame QC ({n_total} frames):")
            if fwhm_vals:
                med = statistics.median(fwhm_vals)
                print(f"    Median FWHM: {med:.2f}\"")
            print(f"    OK: {n_ok}, flagged: {n_total - n_ok}")
            if flags:
                flag_str = ", ".join(f"{k}={v}" for k, v in
                                     sorted(flags.items(), key=lambda x: -x[1]))
                print(f"    Flags: {flag_str}")

        # Show astrometry summary
        astrom_file = base / "astrometry" / "astrom_report.json"
        if astrom_file.exists():
            with open(astrom_file) as f:
                astrom = json.load(f)
            n_total = len(astrom)
            n_ok = sum(1 for r in astrom if r.get("status") == "OK")
            n_warn = sum(1 for r in astrom if r.get("status") == "HIGH_RMS")
            n_fail = sum(1 for r in astrom if r.get("status") == "FAILED")
            rms_vals = [r["rms_arcsec"] for r in astrom
                        if r.get("rms_arcsec") is not None]
            print(f"\n  Astrometry ({n_total} frames):")
            print(f"    Solved: {n_ok + n_warn}/{n_total} "
                  f"({100 * (n_ok + n_warn) / n_total:.0f}%)" if n_total else
                  f"    Solved: 0/0")
            if n_warn:
                print(f"    HIGH_RMS warnings: {n_warn}")
            if n_fail:
                print(f"    Failed: {n_fail}")
            if rms_vals:
                med_rms = statistics.median(rms_vals)
                print(f"    Median RMS: {med_rms:.3f}\"")

        # Show photometry summary
        phot_file = base / "photometry" / "phot_report.json"
        if phot_file.exists():
            with open(phot_file) as f:
                phot = json.load(f)
            n_total = len(phot)
            n_ok = sum(1 for r in phot if r.get("status") == "OK")
            zp_by_filt = {}
            for r in phot:
                if r.get("status") == "OK" and r.get("zp") is not None:
                    filt = r.get("filter", "?")
                    zp_by_filt.setdefault(filt, []).append(r["zp"])
            print(f"\n  Photometry ({n_total} frames):")
            print(f"    Calibrated: {n_ok}/{n_total}")
            for filt in sorted(zp_by_filt):
                zps = zp_by_filt[filt]
                med_zp = statistics.median(zps)
                print(f"    ZP [{filt}]: {med_zp:.2f} (n={len(zps)})")

        # Show PSF summary
        psf_file = base / "psf" / "psf_report.json"
        if psf_file.exists():
            with open(psf_file) as f:
                psf = json.load(f)
            n_total = len(psf)
            n_ok = sum(1 for r in psf if r.get("status") == "OK")
            fwhms = [r["fwhm_arcsec"] for r in psf
                     if r.get("status") == "OK" and r.get("fwhm_arcsec") is not None]
            print(f"\n  PSF ({n_total} frames):")
            print(f"    Models: {n_ok}/{n_total}")
            if fwhms:
                print(f"    Median FWHM: {statistics.median(fwhms):.2f}\"")
    print()


def cmd_inventory(args):
    from .stacking import build_target_inventory, print_inventory

    telescopes = [args.telescope] if args.telescope else list(config.TELESCOPES)
    years = [int(y) for y in args.year.split(",")] if args.year else None
    inventory = build_target_inventory(years=years, telescopes=telescopes)
    print_inventory(inventory)


def cmd_stack(args):
    from .stacking import run_stacking

    years = [int(y) for y in args.year.split(",")] if args.year else None
    run_stacking(
        target=args.target,
        years=years,
        telescope=args.telescope,
        force=args.force,
    )


def main():
    parser = argparse.ArgumentParser(
        prog="pipeline",
        description=f"OHP M2 Imaging Reduction Pipeline v{__version__}",
    )
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Debug-level logging")
    parser.add_argument("--data-root", type=Path, metavar="DIR",
                        help="Data directory (default: $OHP_DATA_ROOT or "
                             "/Users/hjmcc/Work/data)")
    sub = parser.add_subparsers(dest="command", required=True)

    # ── scan ───────────────────────────────────────────────────────────
    p_scan = sub.add_parser("scan",
                            help="Scan raw FITS directory and build archive JSON")
    p_scan.add_argument("--year", required=True, help="Observing year (e.g. 2025)")
    p_scan.add_argument("--telescope", choices=list(config.TELESCOPES),
                        help="Scan single telescope (default: all)")

    # ── calibrate / reduce ─────────────────────────────────────────────
    p_cal = sub.add_parser("calibrate", help="Build master cals and reduce science")
    p_cal.add_argument("--year", required=True, help="Observing year (e.g. 2025)")
    p_cal.add_argument("--telescope", choices=list(config.TELESCOPES),
                        help="Process single telescope (default: all)")
    p_cal.add_argument("--force", action="store_true",
                        help="Rebuild even if outputs exist")

    p_red = sub.add_parser("reduce", help="Alias for calibrate")
    p_red.add_argument("--year", required=True, help="Observing year (e.g. 2025)")
    p_red.add_argument("--telescope", choices=list(config.TELESCOPES),
                        help="Process single telescope (default: all)")
    p_red.add_argument("--force", action="store_true",
                        help="Rebuild even if outputs exist")

    # ── status ────────────────────────────────────────────────────────
    p_st = sub.add_parser("status", help="Show pipeline output status")
    p_st.add_argument("--year", required=True, help="Observing year")
    p_st.add_argument("--telescope", choices=list(config.TELESCOPES),
                        help="Show single telescope")

    # ── qc ─────────────────────────────────────────────────────────────
    p_qc = sub.add_parser("qc", help="Post-calibration QC stats on reduced frames")
    p_qc.add_argument("--year", required=True, help="Observing year (e.g. 2025)")
    p_qc.add_argument("--telescope", choices=list(config.TELESCOPES),
                        help="Process single telescope (default: all)")
    p_qc.add_argument("--force", action="store_true",
                        help="Recompute even if outputs exist")

    # ── astrometry ──────────────────────────────────────────────────
    p_ast = sub.add_parser("astrometry", help="Astrometric calibration via Gaia DR3")
    p_ast.add_argument("--year", required=True, help="Observing year (e.g. 2025)")
    p_ast.add_argument("--telescope", choices=list(config.TELESCOPES),
                        help="Process single telescope (default: all)")
    p_ast.add_argument("--force", action="store_true",
                        help="Re-solve even if WCS already present")

    # ── photometry ────────────────────────────────────────────────────
    p_phot = sub.add_parser("photometry",
                            help="Photometric calibration via PS1 cross-match")
    p_phot.add_argument("--year", required=True, help="Observing year (e.g. 2025)")
    p_phot.add_argument("--telescope", choices=list(config.TELESCOPES),
                        help="Process single telescope (default: all)")
    p_phot.add_argument("--force", action="store_true",
                        help="Re-calibrate even if PHOTZP already present")

    # ── psf ───────────────────────────────────────────────────────────
    p_psf = sub.add_parser("psf",
                           help="PSF modeling via SExtractor + PSFex")
    p_psf.add_argument("--year", required=True, help="Observing year (e.g. 2025)")
    p_psf.add_argument("--telescope", choices=list(config.TELESCOPES),
                        help="Process single telescope (default: all)")
    p_psf.add_argument("--force", action="store_true",
                        help="Reprocess even if PSF models exist")

    # ── inventory ────────────────────────────────────────────────────
    p_inv = sub.add_parser("inventory",
                           help="List stackable targets across all years")
    p_inv.add_argument("--year", help="Comma-separated years (default: all)")
    p_inv.add_argument("--telescope", choices=list(config.TELESCOPES))

    # ── stack ────────────────────────────────────────────────────────
    p_stack = sub.add_parser("stack",
                             help="Coadd images with SWarp (align + flux-scale)")
    p_stack.add_argument("--target", help="Stack single target (canonical name)")
    p_stack.add_argument("--year", help="Comma-separated years (default: all)")
    p_stack.add_argument("--telescope", choices=list(config.TELESCOPES))
    p_stack.add_argument("--force", action="store_true",
                         help="Re-stack even if output exists")

    # ── all ───────────────────────────────────────────────────────────
    p_all = sub.add_parser("all", help="Run full pipeline")
    p_all.add_argument("--year", required=True)
    p_all.add_argument("--telescope", choices=list(config.TELESCOPES))
    p_all.add_argument("--force", action="store_true")

    args = parser.parse_args()
    if args.data_root:
        config.DATA_ROOT = args.data_root
    setup_logging(args.verbose)

    if args.command == "scan":
        cmd_scan(args)
    elif args.command in ("calibrate", "reduce"):
        cmd_calibrate(args)
    elif args.command == "qc":
        cmd_qc(args)
    elif args.command == "status":
        cmd_status(args)
    elif args.command == "all":
        cmd_calibrate(args)
        cmd_qc(args)
        cmd_astrometry(args)
        cmd_photometry(args)
        cmd_psf(args)
    elif args.command == "astrometry":
        cmd_astrometry(args)
    elif args.command == "photometry":
        cmd_photometry(args)
    elif args.command == "psf":
        cmd_psf(args)
    elif args.command == "inventory":
        cmd_inventory(args)
    elif args.command == "stack":
        cmd_stack(args)


if __name__ == "__main__":
    main()
