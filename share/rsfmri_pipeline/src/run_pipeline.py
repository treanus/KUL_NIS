#!/usr/bin/env python3
"""
run_pipeline.py — Top-level orchestrator

Runs the full hierarchical rsfMRI pipeline for a given condition profile:

  Step 0  mri_synthseg   (optional, run separately if not done)
  Step 1  Seed-based analysis (SBA) → per-subject Fisher-Z maps
  Step 2  RSN-level FC matrix      → per-subject 10×10 matrix
  Step 3  Masked ICA               → per-subject ICA components

Usage:
  python run_pipeline.py --profile MDD
  python run_pipeline.py --profile MDD --subjects PT01 PT02
  python run_pipeline.py --profile MDD --hv-only
  python run_pipeline.py --profile MDD --steps 1 3
  python run_pipeline.py --list-profiles
  python run_pipeline.py --list-seeds MDD

Examples:
  # Full pipeline for all subjects, MDD profile
  python run_pipeline.py --profile MDD

  # Only SBA + FC matrix (no ICA) for healthy volunteers
  python run_pipeline.py --profile Full --hv-only --steps 1 2

  # Re-run ICA for a single patient
  python run_pipeline.py --profile OCD --subjects PT01 --steps 3 --force
"""

import argparse
import logging
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from config import CONDITION_PROFILES, SEED_CATALOG, get_profile

from step1_sba        import main as run_step1
from step2_rsn_fc     import main as run_step2
from step3_masked_ica import main as run_step3
from step4_compare    import main as run_step4

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

STEP_LABELS = {
    1: "Seed-Based Analysis (SBA)",
    2: "RSN-level FC Matrix",
    3: "Masked ICA",
    4: "Normative Comparison (z-scores)",
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def list_profiles() -> None:
    print("\nAvailable profiles:\n")
    for name, cfg in sorted(CONDITION_PROFILES.items()):
        desc  = cfg.get("description", "").strip().splitlines()[0][:70]
        seeds = cfg.get("seeds", [])
        print(f"  {name:<25}  {desc}")
        print(f"  {'':25}  Seeds ({len(seeds)}): {', '.join(seeds)}")
        print()


def list_seeds(profile_name: str) -> None:
    profile = get_profile(profile_name)
    seeds   = profile["seeds"]
    print(f"\nSeeds for profile '{profile_name}':\n")
    for sn in seeds:
        cfg  = SEED_CATALOG[sn]
        hemi = cfg.get("hemi", "both")
        print(f"  {sn:<30}  atlas={cfg['atlas']}  hemi={hemi}")
    print()


# ---------------------------------------------------------------------------
# Argument forwarding
# ---------------------------------------------------------------------------

def build_step_args(args: argparse.Namespace, step: int) -> list[str]:
    """Build sys.argv equivalent to pass into a step's main()."""
    argv = []
    # Step 1 requires --profile; step 2 accepts it optionally (enables the
    # seed<->RSN coupling matrix); step 3 has no --profile argument at all
    # (passing it would make step3's own argparse fail) and doesn't need one —
    # it always runs across every Smith-10 RSN regardless of profile.
    if step in (1, 2, 4):
        argv += ["--profile", args.profile]
    if args.subjects:
        # Step 4 uses --patients for patient IDs; steps 1-3 use --subjects
        flag = "--patients" if step == 4 else "--subjects"
        argv += [flag] + args.subjects
    if step != 4 and args.hv_only:
        argv.append("--hv-only")
    if args.force:
        argv.append("--force")
    if step == 1 and args.save_runs:
        argv.append("--save-runs")
    return argv


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--profile", metavar="NAME",
                   help="Condition profile (see --list-profiles)")
    p.add_argument("--subjects", nargs="+", metavar="ID",
                   help="Subject IDs to process (default: all)")
    p.add_argument("--hv-only",  action="store_true",
                   help="Process healthy volunteers only (HV01–HV09)")
    p.add_argument("--steps", nargs="+", type=int,
                   choices=[1, 2, 3, 4], default=[1, 2, 3, 4], metavar="{1,2,3,4}",
                   help="Which steps to run (default: 1 2 3)")
    p.add_argument("--force",     action="store_true",
                   help="Recompute even if outputs already exist")
    p.add_argument("--save-runs", action="store_true",
                   help="(Step 1) Also save per-run SBA maps")
    p.add_argument("--list-profiles", action="store_true",
                   help="Print available profiles and exit")
    p.add_argument("--list-seeds",    metavar="PROFILE",
                   help="Print seeds for a profile and exit")
    return p.parse_args()


def main() -> None:
    args = parse_args()

    if args.list_profiles:
        list_profiles()
        return

    if args.list_seeds:
        list_seeds(args.list_seeds)
        return

    if not args.profile:
        log.error("--profile is required. Use --list-profiles to see options.")
        sys.exit(1)

    profile = get_profile(args.profile)
    log.info("=" * 60)
    log.info("Pipeline start")
    log.info("  Profile  : %s", args.profile)
    log.info("  Steps    : %s", args.steps)
    log.info("  Subjects : %s", args.subjects or ("HV only" if args.hv_only else "all"))
    log.info("=" * 60)

    step_fns = {1: run_step1, 2: run_step2, 3: run_step3, 4: run_step4}

    for step in sorted(args.steps):
        log.info("")
        log.info(">>> Step %d — %s", step, STEP_LABELS[step])
        sys.argv = [f"step{step}"] + build_step_args(args, step)
        try:
            step_fns[step]()
        except SystemExit as exc:
            if exc.code not in (0, None):
                log.error("Step %d exited with code %s", step, exc.code)
                sys.exit(exc.code)

    log.info("")
    log.info("=" * 60)
    log.info("Pipeline complete.")


if __name__ == "__main__":
    main()
