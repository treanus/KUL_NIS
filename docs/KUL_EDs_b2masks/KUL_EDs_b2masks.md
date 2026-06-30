# KUL_EDs_b2masks.py

Swiss-knife distance tool for binary NIfTI masks. Computes minimum Euclidean distance, Hausdorff distance, 95th-percentile Hausdorff distance, and average symmetric surface distance (ASSD) between one or more pairs of masks, and writes a per-pair report (PNG snapshot, HTML report, CSV).

## Modes

- **One-to-one**: `-a A.nii.gz -b B.nii.gz` — compare a single pair.
- **One-to-many**: `-a A.nii.gz -b B1.nii.gz B2.nii.gz ...` — compare one source mask against several targets (parallelizable with `-n`).
- **All-pairs**: `-a A.nii.gz B.nii.gz C.nii.gz` (no `-b`), or `-l masks.txt` (one path per line) — compare every combination of ≥2 masks.

## Multi-class masks

`--classes 1 2 --class-labels core edema` extracts the listed integer label values from mask A separately and reports distances per class, instead of treating the mask as a single binary volume.

## Usage

```
KUL_EDs_b2masks.py -a A.nii.gz -b B.nii.gz -o out
KUL_EDs_b2masks.py -a A.nii.gz -b B1.nii.gz B2.nii.gz -o out -n 4
KUL_EDs_b2masks.py -a A.nii.gz B.nii.gz C.nii.gz -o out
KUL_EDs_b2masks.py -l masks.txt -o out -n -1 --csv
KUL_EDs_b2masks.py -a seg.nii.gz -b tract.nii.gz -o out --classes 1 2 --class-labels core edema
```

## Options

| Option | Description |
|---|---|
| `-a, --mask-a MASK [MASK ...]` | Source mask(s). With `-b`: one-to-one/one-to-many. Without `-b` (≥2 masks): all-pairs. |
| `-l, --list FILE` | Text file of mask paths (one per line) for all-pairs mode. Mutually exclusive with `-a`. |
| `-b, --mask-b MASK [MASK ...]` | Target mask(s). Omit for all-pairs mode. |
| `-o, --out` | Output prefix/path (default: `KUL_EDs`). Results go to `<out>_output/`. |
| `--classes VAL [VAL ...]` | Integer label values to extract from mask A separately. |
| `--class-labels LABEL [LABEL ...]` | Human-readable names for `--classes` (must match count). |
| `--metric {euclidean,hausdorff,mean_surface,all}` | Metric to report (default `all`; all are always computed internally). |
| `--surface` | Restrict Hausdorff/ASSD to 1-voxel erosion edge voxels instead of all mask voxels via EDT lookup. |
| `--no-maps` | Skip writing float distance-map NIfTIs (faster; snapshot is still generated). |
| `--bg-image NII` | Anatomical reference (T1w, FA, ...) used as snapshot background. Without `--volume-render`: orthographic 3-panel `mrview` snapshot (needs `xvfb-run`). |
| `--volume-render` | VTK off-screen glass-brain 3-D snapshot (superior/lateral/posterior views) instead of `mrview`. Requires `--bg-image`. No display/`xvfb-run` needed. |
| `--csv` | Write a `_metrics.csv` summary alongside the per-pair output. |
| `-n, --workers N` | Parallel workers for multi-pair comparisons (`-1` = all CPUs; default `1`). |
| `-v, --verbose` | Verbose logging. |

## Output

For each mask pair, a subdirectory is written under `<out>_output/` containing:
- a text/metrics summary
- a PNG snapshot (`mrview` 3-panel or VTK glass-brain render, depending on `--bg-image`/`--volume-render`)
- an HTML report

With `--csv`, an aggregate `_metrics.csv` is also written across all computed pairs, and a summary table is printed to stdout.

## Performance notes

In one-to-many mode, mask A's distance-transform (EDT) is computed once per worker process (via `_init_worker`) and cached, so the per-pair cost only covers mask B's EDT and the comparison itself.

## Dependencies

Python 3, `numpy`, `scipy` (distance transforms / morphology), `nibabel` (NIfTI I/O). `--bg-image` snapshots additionally need `mrview` (MRtrix3) with `xvfb-run`, or VTK for `--volume-render`.
