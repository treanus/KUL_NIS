#!/usr/bin/env python3
#
# Stefan Sunaert, 2025-08-24
# Compute FAT1w weighted MR image from T1w and FA images using MRtrix3
# M&M taken from Goedemans et al., Imaging Neurosci 2024. doi: 10.1162/imag_a_00139
#
import argparse
import subprocess
import os
import sys

def print_usage_and_exit(parser):
    """Print usage information and exit."""
    print("\nUsage: Compute FAT1w image from T1w and FA images using MRtrix3\n")
    parser.print_help()
    sys.exit(1)

def main():
    parser = argparse.ArgumentParser(
        description="Compute FAT1w from T1w and FA images using MRtrix3",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Using BIDS participant ID (without 'sub-'):
  KUL_FAT1.py -p 001
  
  # Using explicit file paths:
  KUL_FAT1.py -t1 path/to/T1w.nii.gz -fa path/to/fa.nii.gz -o path/to/output/FAT1w.nii.gz
  
  # With smoothing (default 2.0mm):
  KUL_FAT1.py -p 001 -s
  
  # With custom smoothing:
  KUL_FAT1.py -p 001 -s 3.0
"""
    )

    parser.add_argument('-p', '--participant', type=str, help="Participant ID")
    parser.add_argument('-t1', type=str, help="Path to T1w image")
    parser.add_argument('-fa', type=str, help="Path to FA image")
    parser.add_argument('-o', '--output', type=str, help="Path to output FAT1w image")
    parser.add_argument('-s', '--smooth', nargs='?', const=2.0, type=float, default=None,
                        help="Smooth FA with FWHM in mm (default 2.0 if flag provided without value)")

    args = parser.parse_args()

    # Show help if no arguments or -h/--help is provided
    if len(sys.argv) == 1 or args == argparse.Namespace(participant=None, t1=None, fa=None, output=None, smooth=None):
        print_usage_and_exit(parser)

    if args.participant:
        participant = args.participant
        t1_path = f"./RESULTS/sub-{participant}/Anat/T1w.nii.gz"
        fa_path = f"./dwiprep/sub-{participant}/sub-{participant}/qa/fa_reg2T1w.nii.gz"
        out_path = f"./BIDS/derivatives/KUL_compute/sub-{participant}/KUL_FAT1/FAT1w.nii.gz"
    elif all([args.t1, args.fa, args.output]):
        t1_path = args.t1
        fa_path = args.fa
        out_path = args.output
    else:
        parser.error("Either provide -p or all of -t1 -fa -o")

    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    # Temporary files
    base_dir = os.path.dirname(out_path)
    fa_regrid = os.path.join(base_dir, "fa_regrid.nii.gz")
    fa_processed = fa_regrid

    # Regrid FA to T1w
    subprocess.run(["mrgrid", "-force", fa_path, "regrid", "-template", t1_path, fa_regrid], check=True)

    if args.smooth is not None:
        fa_smooth = os.path.join(base_dir, "fa_smooth.nii.gz")
        subprocess.run(["mrfilter", "-force", fa_regrid, "smooth", "-fwhm", str(args.smooth), fa_smooth], check=True)
        fa_processed = fa_smooth

    # Compute FAT1w = sqrt(FA) * T1w
    subprocess.run(["mrcalc", "-force", fa_processed, "-sqrt", t1_path, "-mult", out_path], check=True)

if __name__ == "__main__":
    main()