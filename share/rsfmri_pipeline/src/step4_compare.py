#!/usr/bin/env python3
# Run with:  python3 step4_compare.py [options]
import matplotlib
matplotlib.use("Agg")

"""
Step 4 — Normative Comparison

Two z-score methods (select with --method):

  simple (default):
    z[voxel] = (patient_FC[voxel] − mean(HV_FC)[voxel]) / SD(HV_FC)[voxel]
    p-values via Crawford-Howell modified t-test:
      t = z · √(n_HV / (n_HV + 1)),   df = n_HV − 1
    Same test can be applied to per-HV LOO z-maps for QC and threshold calibration.
    No covariate correction; interpretable and appropriate for small n.

  glm:
    OLS GLM on HV FC maps: FC ~ 1 + age [+ sex] [+ n_runs]
    Prediction-interval z-score: z = (y_obs − ŷ) / √(σ²(1 + leverage))
    This is the covariate-adjusted generalisation of Crawford-Howell;
    distributes as t(n_HV − n_params).

Both methods apply:
  • FDR correction (Benjamini–Hochberg, q < 0.05)
  • LOO-calibrated empirical threshold (95th percentile of pooled |LOO z|)

Runs in two modes (independently or together):

  Seed mode (--profile required):
    Normative comparison on anatomical-seed SBA maps from step 1.

  RSN mode (--rsn flag):
    Normative comparison on RSN whole-brain FC maps from step 2.

Outputs
-------
  Seed mode:
    analysis/sba/normative/
      seed-<name>_HV_GLM_mean.nii.gz        HV mean (simple) or GLM-predicted mean (glm)
      seed-<name>_HV_residual_SD.nii.gz     HV SD (simple) or GLM residual SD (glm)
    analysis/sba/sub-<PT>/
      sub-<PT>_seed-<name>_desc-zscore_statmap.nii.gz
      sub-<PT>_seed-<name>_desc-pval_statmap.nii.gz
      sub-<PT>_seed-<name>_desc-fdrsig_statmap.nii.gz   (z at FDR q<0.05, 0 elsewhere)
      sub-<PT>_seed-<name>_desc-loosig_statmap.nii.gz   (z at |z|>LOO thr, 0 elsewhere)
    analysis/comparison/
      sub-<PT>_profile-<name>_roi_zscores.csv
      normative_QC.csv
      normative_covariates.csv  (glm mode only)

  RSN mode — same structure under analysis/rsn_fc/

Usage
-----
  python step4_compare.py --profile MDD --patients PT01 PT02
  python step4_compare.py --profile MDD --method glm --patients PT01 PT02
  python step4_compare.py --rsn --patients PT01 PT02
  python step4_compare.py --profile MDD --rsn --patients PT01 PT02 --force
  python step4_compare.py --profile MDD --roi-only
"""

import argparse
import logging
import sys
from collections import namedtuple
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
from config import (
    ANALYSIS_DIR, BASE_DIR, HV_SUBJECTS, PIPELINE_DIR, PT_SUBJECTS,
    RSN_FC_DIR, SBA_DIR, get_profile, get_seed, seed_subdir,
)
from utils import (
    get_auditory_group_mask, get_language_group_mask,
    get_language_network_mask, get_seed_mask, get_yeo17_subject_masks,
    resample_to, YEO17_NETWORK_NAMES,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

NORM_DIR     = SBA_DIR    / "normative"
RSN_NORM_DIR = RSN_FC_DIR / "normative"
COMP_DIR     = ANALYSIS_DIR / "comparison"
BIDS_DIR     = BASE_DIR / "BIDS"

for _d in [NORM_DIR, RSN_NORM_DIR, COMP_DIR]:
    _d.mkdir(parents=True, exist_ok=True)

MIN_SD             = 1e-4
DEFAULT_COVARIATES = ["age", "sex", "n_runs"]

_SEX_LUT = {
    "F": 0.5, "FEMALE": 0.5, "1": 0.5,
    "M": -0.5, "MALE": -0.5, "0": -0.5,
}

# Resampled Yeo17 atlas (produced once by step0 step 4a) — used as group-level
# ROI masks in RSN normative comparison.
_YEO17_RS_PATH = (
    PIPELINE_DIR
    / "Yeo_JNeurophysiol11_MNI152"
    / "Yeo2011_17Networks_MNI152_rs2mm.nii.gz"
)

_ALL_RSN_NAMES = YEO17_NETWORK_NAMES + ["Language", "Auditory"]


# ─────────────────────────────────────────────────────────────────────────────
# NormCtx — bundles path callables for one analysis item (seed or RSN)
# ─────────────────────────────────────────────────────────────────────────────

NormCtx = namedtuple("NormCtx", [
    "label",         # str: item name ("sgACC", "default_mode", …)
    "type_str",      # "seed" or "rsn"
    "hv_sba",        # Callable[str] → Path: HV SBA map
    "pt_sba",        # Callable[str] → Path: patient SBA map
    "zscore_out",    # Callable[str] → Path
    "pval_out",      # Callable[str] → Path
    "fdrsig_out",    # Callable[str] → Path
    "loosig_out",    # Callable[str] → Path
    "glm_mean_out",  # Callable[]    → Path  (HV mean for simple; GLM mean for glm)
    "resid_sd_out",  # Callable[]    → Path  (HV SD for simple; residual SD for glm)
])


# ─────────────────────────────────────────────────────────────────────────────
# Participant metadata
# ─────────────────────────────────────────────────────────────────────────────

def _count_runs(subject: str) -> int:
    fls = list((BIDS_DIR / f"sub-{subject}").rglob("*_bold.nii.gz"))
    return max(len(fls), 1)


def load_metadata(subjects: list[str]) -> pd.DataFrame:
    tsv = pd.read_csv(BIDS_DIR / "participants.tsv", sep="\t")
    tsv["subject"] = tsv["participant_id"].str.replace("sub-", "", regex=False)
    tsv = tsv.set_index("subject")

    rows: dict[str, dict] = {}
    for sid in subjects:
        row: dict = {"age": float(tsv.loc[sid, "age"])}
        if "sex" in tsv.columns:
            raw = str(tsv.loc[sid, "sex"]).strip().upper()
            if raw in _SEX_LUT:
                row["sex"] = _SEX_LUT[raw]
            elif raw not in {"NAN", ""}:
                try:
                    row["sex"] = float(raw)
                except ValueError:
                    log.warning("  Unknown sex value '%s' for %s", raw, sid)
        row["n_runs"] = _count_runs(sid)
        rows[sid] = row

    df = pd.DataFrame.from_dict(rows, orient="index")
    log.info("Metadata loaded — columns: %s", list(df.columns))
    return df


def select_covariates(
    requested: list[str],
    meta: pd.DataFrame,
    hv_ids: list[str],
) -> list[str]:
    hv = meta.loc[hv_ids]
    usable: list[str] = []
    for cov in requested:
        if cov not in meta.columns:
            log.warning("  Covariate '%s' not in participants.tsv — skipped", cov)
            continue
        if hv[cov].isna().any():
            log.warning("  Covariate '%s' has NaN in HV rows — skipped", cov)
            continue
        if hv[cov].std() < 1e-8:
            log.warning("  Covariate '%s' constant across HVs — skipped", cov)
            continue
        usable.append(cov)
    if usable:
        df_resid = len(hv_ids) - 1 - len(usable)
        log.info("  Active covariates: %s  |  residual df = %d", usable, df_resid)
        if df_resid < 3:
            log.warning("  Only %d residual df — estimates will be very uncertain", df_resid)
    else:
        log.warning("  No usable covariates — intercept-only model")
    return usable


# ─────────────────────────────────────────────────────────────────────────────
# NormativeGLM  (used by --method glm only)
# ─────────────────────────────────────────────────────────────────────────────

class NormativeGLM:
    """
    OLS normative model: Y ~ 1 + cov1 + cov2 + ...
    z-score uses the prediction interval: z = (y − ŷ) / √[σ²(1 + h)]
    Pass Y as (n_hv, n_vox) for voxelwise use, (n_hv,) for scalar.
    """

    def __init__(self, cov_names: list[str], age_mean: float = 0.0):
        self.cov_names = list(cov_names)
        self.age_mean  = age_mean
        self.beta_    = None
        self.sigma2_  = None
        self.XtXinv_  = None
        self.df_      = None
        self.hv_ids_  = None
        self._X_fit   = None
        self._Y_fit   = None

    def _row(self, sid: str, meta: pd.DataFrame) -> np.ndarray:
        r = [1.0]
        for c in self.cov_names:
            v = float(meta.loc[sid, c])
            if c == "age":
                v -= self.age_mean
            r.append(v)
        return np.array(r, dtype=np.float64)

    def _design(self, ids: list[str], meta: pd.DataFrame) -> np.ndarray:
        return np.stack([self._row(s, meta) for s in ids])

    def fit(self, hv_ids: list[str], Y: np.ndarray, meta: pd.DataFrame) -> "NormativeGLM":
        X    = self._design(hv_ids, meta)
        n, p = X.shape
        df   = n - p
        if df <= 0:
            raise RuntimeError(
                f"df exhausted: {n} HVs, {p} params ({['intercept'] + self.cov_names})"
            )
        XtXinv = np.linalg.pinv(X.T @ X)
        beta   = XtXinv @ X.T @ Y
        resid  = Y - X @ beta
        sigma2 = (
            (resid ** 2).sum(axis=0) / df if Y.ndim > 1
            else float((resid ** 2).sum() / df)
        )
        self.beta_   = beta
        self.sigma2_ = sigma2
        self.XtXinv_ = XtXinv
        self.df_     = df
        self.hv_ids_ = list(hv_ids)
        self._X_fit  = X
        self._Y_fit  = Y
        return self

    def predict(self, sid: str, meta: pd.DataFrame) -> np.ndarray:
        return self._row(sid, meta) @ self.beta_

    def _leverage(self, sid: str, meta: pd.DataFrame) -> float:
        x = self._row(sid, meta)
        return float(x @ self.XtXinv_ @ x)

    def z_score(self, y_obs: np.ndarray, sid: str, meta: pd.DataFrame) -> np.ndarray:
        pred_sd = np.sqrt(
            np.maximum(self.sigma2_, MIN_SD ** 2) * (1.0 + self._leverage(sid, meta))
        )
        return (y_obs - self.predict(sid, meta)) / pred_sd

    def t_critical(self, alpha: float = 0.05) -> float:
        return float(stats.t.ppf(1.0 - alpha / 2, df=self.df_))

    def mean_at_hv_centroid(self, meta: pd.DataFrame) -> np.ndarray:
        return self._X_fit.mean(axis=0) @ self.beta_

    def extrapolation_flags(self, sid: str, meta: pd.DataFrame) -> dict[str, bool]:
        flags: dict[str, bool] = {}
        for c in self.cov_names:
            hv_vals = [float(meta.loc[hv, c]) for hv in self.hv_ids_]
            pt_val  = float(meta.loc[sid, c])
            flags[c] = pt_val < min(hv_vals) or pt_val > max(hv_vals)
        return flags

    def loo_analysis(self, meta: pd.DataFrame, alpha: float = 0.05
                     ) -> tuple[pd.DataFrame, float]:
        """Single LOO pass: QC table + empirical threshold."""
        records:   list[dict]        = []
        all_abs_z: list[np.ndarray]  = []

        for i, hv in enumerate(self.hv_ids_):
            loo_ids = [s for j, s in enumerate(self.hv_ids_) if j != i]
            Y_loo   = (np.delete(self._Y_fit, i, axis=0)
                       if self._Y_fit.ndim > 1
                       else np.delete(self._Y_fit, i))
            z_s = np.nan
            try:
                m   = NormativeGLM(self.cov_names, self.age_mean).fit(loo_ids, Y_loo, meta)
                z   = m.z_score(self._Y_fit[i], hv, meta)
                z_s = float(np.nanmean(z))
                all_abs_z.append(np.abs(np.asarray(z).ravel()))
            except Exception as exc:
                log.warning("  LOO failed for %s: %s", hv, exc)
            records.append({
                "subject":      hv,
                "loo_z_global": round(z_s, 3) if not np.isnan(z_s) else np.nan,
                "outlier_flag": abs(z_s) > 2.0 if not np.isnan(z_s) else False,
            })

        qc_df = pd.DataFrame(records)

        if all_abs_z:
            pooled  = np.concatenate(all_abs_z)
            loo_thr = float(np.percentile(pooled, (1.0 - alpha) * 100.0))
            log.info(
                "  LOO empirical threshold (α=%.2f): %.3f  (parametric t_crit: %.3f)",
                alpha, loo_thr, self.t_critical(alpha),
            )
        else:
            log.warning("  LOO: all folds failed — using parametric t_critical")
            loo_thr = self.t_critical(alpha)

        return qc_df, loo_thr


# ─────────────────────────────────────────────────────────────────────────────
# BH-FDR
# ─────────────────────────────────────────────────────────────────────────────

def _network_extent_mm3(fc_img: nib.Nifti1Image,
                         z_thresh: float,
                         thresholded: bool,
                         vox_vol: float) -> float:
    """Volume (mm³) of suprathreshold voxels in a Fisher-Z FC map."""
    data = fc_img.get_fdata(dtype=np.float32)
    mask = (data != 0) if thresholded else (data >= z_thresh)
    return float(mask.sum()) * vox_vol


def _sig_volume_mm3(sig_img: nib.Nifti1Image, vox_vol: float) -> float:
    """Volume (mm³) of nonzero voxels in an FDR/LOO significance map."""
    return float((sig_img.get_fdata(dtype=np.float32) != 0).sum()) * vox_vol


def _bh_fdr(p_vals: np.ndarray, alpha: float = 0.05) -> np.ndarray:
    """Benjamini-Hochberg FDR. Returns q-values; significant where q <= alpha."""
    n     = len(p_vals)
    order = np.argsort(p_vals)
    rank  = np.empty(n, dtype=np.intp)
    rank[order] = np.arange(1, n + 1)
    q = np.minimum(p_vals * n / rank, 1.0)
    q_sorted = np.minimum.accumulate(q[order][::-1])[::-1]
    q_out    = np.empty(n)
    q_out[order] = q_sorted
    return q_out


# ─────────────────────────────────────────────────────────────────────────────
# Paths — seed mode
# ─────────────────────────────────────────────────────────────────────────────

def _seed_dir(sub: str, seed: str = "") -> Path:
    # seed_subdir() routes Lausanne-based seeds one level down, matching
    # step1_sba.py's writer. seed="" keeps the flat subject root for callers
    # that want the directory itself rather than a particular seed's file.
    d = SBA_DIR / f"sub-{sub}" / (seed_subdir(seed) if seed else "")
    d.mkdir(parents=True, exist_ok=True)
    return d

def sba_path(sub: str, seed: str, thresholded: bool = False) -> Path:
    desc = "sbaThresh" if thresholded else "sba"
    return _seed_dir(sub, seed) / f"sub-{sub}_seed-{seed}_desc-{desc}_statmap.nii.gz"

def zscore_path(sub: str, seed: str, thresholded: bool = False) -> Path:
    desc = "zscoreThresh" if thresholded else "zscore"
    return _seed_dir(sub, seed) / f"sub-{sub}_seed-{seed}_desc-{desc}_statmap.nii.gz"

def pval_path(sub: str, seed: str, thresholded: bool = False) -> Path:
    desc = "pvalThresh" if thresholded else "pval"
    return _seed_dir(sub, seed) / f"sub-{sub}_seed-{seed}_desc-{desc}_statmap.nii.gz"

def fdrsig_path(sub: str, seed: str, thresholded: bool = False) -> Path:
    desc = "fdrsigThresh" if thresholded else "fdrsig"
    return _seed_dir(sub, seed) / f"sub-{sub}_seed-{seed}_desc-{desc}_statmap.nii.gz"

def loosig_path(sub: str, seed: str, thresholded: bool = False) -> Path:
    desc = "loosigThresh" if thresholded else "loosig"
    return _seed_dir(sub, seed) / f"sub-{sub}_seed-{seed}_desc-{desc}_statmap.nii.gz"

def glm_mean_path(seed: str, thresholded: bool = False) -> Path:
    tag = "_thresh" if thresholded else ""
    return NORM_DIR / f"seed-{seed}_HV_GLM_mean{tag}.nii.gz"

def resid_sd_path(seed: str, thresholded: bool = False) -> Path:
    tag = "_thresh" if thresholded else ""
    return NORM_DIR / f"seed-{seed}_HV_residual_SD{tag}.nii.gz"


# ─────────────────────────────────────────────────────────────────────────────
# Paths — RSN mode
# ─────────────────────────────────────────────────────────────────────────────

def _rsn_dir(sub: str) -> Path:
    d = RSN_FC_DIR / f"sub-{sub}"
    d.mkdir(parents=True, exist_ok=True)
    return d

def rsn_sba_path(sub: str, rsn: str, thresholded: bool = False) -> Path:
    desc = "sbaThresh" if thresholded else "sba"
    return _rsn_dir(sub) / f"sub-{sub}_rsn-{rsn}_desc-{desc}_statmap.nii.gz"

def rsn_zscore_path(sub: str, rsn: str, thresholded: bool = False) -> Path:
    desc = "zscoreThresh" if thresholded else "zscore"
    return _rsn_dir(sub) / f"sub-{sub}_rsn-{rsn}_desc-{desc}_statmap.nii.gz"

def rsn_pval_path(sub: str, rsn: str, thresholded: bool = False) -> Path:
    desc = "pvalThresh" if thresholded else "pval"
    return _rsn_dir(sub) / f"sub-{sub}_rsn-{rsn}_desc-{desc}_statmap.nii.gz"

def rsn_fdrsig_path(sub: str, rsn: str, thresholded: bool = False) -> Path:
    desc = "fdrsigThresh" if thresholded else "fdrsig"
    return _rsn_dir(sub) / f"sub-{sub}_rsn-{rsn}_desc-{desc}_statmap.nii.gz"

def rsn_loosig_path(sub: str, rsn: str, thresholded: bool = False) -> Path:
    desc = "loosigThresh" if thresholded else "loosig"
    return _rsn_dir(sub) / f"sub-{sub}_rsn-{rsn}_desc-{desc}_statmap.nii.gz"

def rsn_glm_mean_path(rsn: str, thresholded: bool = False) -> Path:
    tag = "_thresh" if thresholded else ""
    return RSN_NORM_DIR / f"rsn-{rsn}_HV_GLM_mean{tag}.nii.gz"

def rsn_resid_sd_path(rsn: str, thresholded: bool = False) -> Path:
    tag = "_thresh" if thresholded else ""
    return RSN_NORM_DIR / f"rsn-{rsn}_HV_residual_SD{tag}.nii.gz"


# ─────────────────────────────────────────────────────────────────────────────
# Shared helper — load all HV maps into a (n_hv, n_vox) matrix
# ─────────────────────────────────────────────────────────────────────────────

def _load_hv_matrix(
    ctx:           NormCtx,
    positive_only: bool  = False,
) -> tuple[list[str], np.ndarray, "nib.Nifti1Image | None", np.ndarray]:
    """
    Load HV SBA maps, resample to common space, stack into float64 matrix.
    Returns (hv_ids, Y, ref_img, brain_mask) where:
      Y is (n_hv, n_vox) — columns outside the group brain mask are zeroed.
      brain_mask is a 1-D bool array (n_vox): True where ≥1 HV had a nonzero
      value in the *original* (pre-resample) map.  Used to suppress interpolation
      bleed at the brain boundary in all downstream output maps.
    Returns ([], empty, None, empty) if fewer than 1 map found.
    """
    hv_imgs: list[nib.Nifti1Image] = []
    hv_ids:  list[str] = []

    for hv in HV_SUBJECTS:
        p = ctx.hv_sba(hv)
        if p.exists():
            hv_imgs.append(nib.load(str(p)))
            hv_ids.append(hv)

    if not hv_imgs:
        return [], np.empty((0,)), None, np.empty((0,), dtype=bool)

    ref    = hv_imgs[0]
    stacks = [ref.get_fdata(dtype=np.float32).ravel()]
    for img in hv_imgs[1:]:
        stacks.append(
            resample_to(img, ref, interpolation="linear")
            .get_fdata(dtype=np.float32).ravel()
        )
    Y = np.stack(stacks, axis=0).astype(np.float64)   # (n_hv, n_vox)

    # Group brain mask: voxels nonzero in at least one HV (before interpolation bleed
    # could add spurious values, so we use the union of original nonzero patterns).
    # Linear resampling spreads signal into adjacent voxels at the brain boundary;
    # zeroing those columns prevents floating voxels outside the brain in all outputs.
    brain_mask = (Y != 0).any(axis=0)                 # (n_vox,)  bool
    Y[:, ~brain_mask] = 0.0

    if positive_only:
        Y = np.maximum(Y, 0.0)

    return hv_ids, Y, ref, brain_mask


# ─────────────────────────────────────────────────────────────────────────────
# Simple normative model
# ─────────────────────────────────────────────────────────────────────────────

def build_simple_model(
    ctx:           NormCtx,
    force:         bool  = False,
    positive_only: bool  = False,
) -> tuple[
    "np.ndarray | None",
    "np.ndarray | None",
    int,
    "nib.Nifti1Image | None",
    float,
    "pd.DataFrame | None",
]:
    """
    Voxelwise HV mean and SD, plus LOO-calibrated threshold.
    Returns (hv_mean, hv_sd, n_hv, ref_img, loo_thr, qc_df).
    hv_mean and hv_sd are 1-D float64 arrays of length n_vox.
    Returns (None, …) on failure.
    """
    hv_ids, Y, ref, brain_mask = _load_hv_matrix(
        ctx, positive_only=positive_only
    )

    if len(hv_ids) < 2:
        log.error("  %s '%s': need ≥2 HV maps, found %d — skipping",
                  ctx.type_str, ctx.label, len(hv_ids))
        return None, None, 0, None, float("nan"), None

    missing = [hv for hv in HV_SUBJECTS if hv not in hv_ids]
    if missing:
        log.warning("  %s '%s': no SBA map for %s", ctx.type_str, ctx.label, missing)

    n_hv  = len(hv_ids)
    shape = ref.shape
    mu    = Y.mean(axis=0)
    sd    = np.maximum(Y.std(axis=0, ddof=1), MIN_SD)

    log.info("  %s '%s': simple model on %d HVs", ctx.type_str, ctx.label, n_hv)

    if not ctx.glm_mean_out().exists() or force:
        mu_save = np.where(brain_mask, mu, 0.0)
        sd_save = np.where(brain_mask, sd, 0.0)
        nib.save(nib.Nifti1Image(mu_save.reshape(shape).astype(np.float32), ref.affine),
                 str(ctx.glm_mean_out()))
        nib.save(nib.Nifti1Image(sd_save.reshape(shape).astype(np.float32), ref.affine),
                 str(ctx.resid_sd_out()))
        log.info("  Saved HV mean + SD maps")

    # LOO calibration
    log.info("  LOO calibration…")
    records:   list[dict]       = []
    all_abs_z: list[np.ndarray] = []

    for i, hv in enumerate(hv_ids):
        loo_mask   = np.ones(n_hv, dtype=bool)
        loo_mask[i] = False
        Y_loo  = Y[loo_mask]
        mu_loo = Y_loo.mean(axis=0)
        sd_loo = np.maximum(Y_loo.std(axis=0, ddof=1), MIN_SD)
        z_i    = (Y[i] - mu_loo) / sd_loo
        z_glob = float(np.nanmean(z_i))
        all_abs_z.append(np.abs(z_i))
        records.append({
            "subject":      hv,
            "loo_z_global": round(z_glob, 3),
            "outlier_flag": abs(z_glob) > 2.0,
        })

    qc_df   = pd.DataFrame(records)
    pooled  = np.concatenate(all_abs_z)
    loo_thr = float(np.percentile(pooled, 95.0))
    log.info("  LOO empirical threshold (simple, α=0.05): %.3f", loo_thr)

    outliers = qc_df.loc[qc_df["outlier_flag"], "subject"].tolist()
    if outliers:
        log.warning("  LOO |z| > 2 for: %s", outliers)

    return mu, sd, n_hv, ref, loo_thr, qc_df, brain_mask


def compute_simple_zscore_map(
    patient:       str,
    ctx:           NormCtx,
    hv_mean:       np.ndarray,
    hv_sd:         np.ndarray,
    n_hv:          int,
    brain_mask:    np.ndarray,
    ref_img:       "nib.Nifti1Image",
    loo_thr:       float,
    positive_only: bool  = False,
    fc_threshold:  float = 0.0,
    force:         bool  = False,
) -> "nib.Nifti1Image | None":
    """
    Simple normative z-score + Crawford-Howell p-values + FDR + LOO maps.
    Saves: zscore, pval, fdrsig, loosig.
    """
    out_z   = ctx.zscore_out(patient)
    out_p   = ctx.pval_out(patient)
    out_fdr = ctx.fdrsig_out(patient)
    out_loo = ctx.loosig_out(patient)

    if all(f.exists() for f in [out_z, out_p, out_fdr, out_loo]) and not force:
        log.info("  [%s | %s] Maps exist — loading", patient, ctx.label)
        return nib.load(str(out_z))

    pt_p = ctx.pt_sba(patient)
    if not pt_p.exists():
        log.error("  [%s | %s] No SBA map — run step1/step2 first", patient, ctx.label)
        return None

    y_obs = (
        resample_to(nib.load(str(pt_p)), ref_img, interpolation="linear")
        .get_fdata(dtype=np.float32).ravel().astype(np.float64)
    )
    if positive_only:
        y_obs = np.maximum(y_obs, 0.0)
    shape = ref_img.shape

    # Crawford-Howell: t = z * sqrt(n / (n+1)),  df = n-1
    z_vec = (y_obs - hv_mean) / hv_sd
    z_vec = np.where(brain_mask, z_vec, 0.0)   # suppress interpolation bleed outside brain
    # Post-hoc FC threshold: zero z-scores at voxels where the HV group mean
    # is below threshold.  Applied after statistics so the normative distribution
    # is estimated from unmodified data; only the output is masked.
    if fc_threshold > 0.0:
        fc_mask = np.abs(hv_mean) >= fc_threshold
        z_vec = np.where(fc_mask, z_vec, 0.0)
    t_vec = z_vec * np.sqrt(n_hv / (n_hv + 1.0))
    p_vec = 2.0 * stats.t.sf(np.abs(t_vec), df=n_hv - 1)
    p_vec = np.where(brain_mask, p_vec, 1.0)
    if fc_threshold > 0.0:
        p_vec = np.where(fc_mask, p_vec, 1.0)

    z_img = nib.Nifti1Image(z_vec.reshape(shape).astype(np.float32), ref_img.affine)
    nib.save(z_img, str(out_z))
    nib.save(nib.Nifti1Image(p_vec.reshape(shape).astype(np.float32), ref_img.affine),
             str(out_p))

    q_vals  = _bh_fdr(p_vec)
    fdr_sig = (q_vals <= 0.05).reshape(shape)
    fdr_map = np.where(fdr_sig, z_vec.reshape(shape), 0.0).astype(np.float32)
    nib.save(nib.Nifti1Image(fdr_map, ref_img.affine), str(out_fdr))

    loo_sig = (np.abs(z_vec) > loo_thr).reshape(shape)
    loo_map = np.where(loo_sig, z_vec.reshape(shape), 0.0).astype(np.float32)
    nib.save(nib.Nifti1Image(loo_map, ref_img.affine), str(out_loo))

    log.info(
        "  [%s | %s] z-map (simple)  FDR sig: %d vox  LOO sig: %d vox  (LOO thr=%.2f)",
        patient, ctx.label, int(fdr_sig.sum()), int(loo_sig.sum()), loo_thr,
    )
    return z_img


def roi_summary_simple(
    patient:      str,
    ctx:          NormCtx,
    roi_mask_img: "nib.Nifti1Image",
    hv_mean:      np.ndarray,
    hv_sd:        np.ndarray,
    n_hv:         int,
    ref_img:      "nib.Nifti1Image",
) -> list[dict]:
    """ROI-level Crawford-Howell test for simple normative model."""
    mask = resample_to(roi_mask_img, ref_img).get_fdata(dtype=np.float32).ravel() > 0
    if not mask.any():
        log.warning("  ROI '%s': empty mask after resampling", ctx.label)
        return []

    hv_roi: list[float] = []
    hv_ok:  list[str]   = []
    for hv in HV_SUBJECTS:
        p = ctx.hv_sba(hv)
        if p.exists():
            img = resample_to(nib.load(str(p)), ref_img, interpolation="linear")
            hv_roi.append(float(img.get_fdata(dtype=np.float32).ravel()[mask].mean()))
            hv_ok.append(hv)

    if len(hv_roi) < 2:
        return []

    pt_p = ctx.pt_sba(patient)
    if not pt_p.exists():
        return []
    pt_fc = float(
        resample_to(nib.load(str(pt_p)), ref_img, interpolation="linear")
        .get_fdata(dtype=np.float32).ravel()[mask].mean()
    )

    hv_arr   = np.array(hv_roi)
    n        = len(hv_arr)
    roi_mean = float(hv_arr.mean())
    roi_sd   = max(float(hv_arr.std(ddof=1)), MIN_SD)

    z      = (pt_fc - roi_mean) / roi_sd
    t_CH   = z * np.sqrt(n / (n + 1.0))
    p_val  = float(2.0 * stats.t.sf(abs(t_CH), df=n - 1))
    t_crit = float(stats.t.ppf(0.975, df=n - 1))

    return [{
        "patient":         patient,
        ctx.type_str:      ctx.label,
        "patient_mean_FC": round(pt_fc, 4),
        "HV_mean_FC":      round(roi_mean, 4),
        "ROI_z_score":     round(z, 3),
        "ROI_t_CH":        round(float(t_CH), 3),
        "p_value":         round(p_val, 4),
        "significant_p05": p_val < 0.05,
        "t_critical":      round(t_crit, 2),
        "df":              n - 1,
        "method":          "simple",
        "n_HV":            n,
        "interpretation":  _interpret_z(z, t_crit),
    }]


# ─────────────────────────────────────────────────────────────────────────────
# GLM normative model
# ─────────────────────────────────────────────────────────────────────────────

def build_normative_model(
    ctx:           NormCtx,
    cov_names:     list[str],
    meta:          pd.DataFrame,
    force:         bool  = False,
    positive_only: bool  = False,
) -> tuple["NormativeGLM | None", "nib.Nifti1Image | None", float, "pd.DataFrame | None", np.ndarray]:
    """Load HV maps, fit voxelwise GLM, run LOO calibration."""
    hv_ids, Y, ref, brain_mask = _load_hv_matrix(
        ctx, positive_only=positive_only
    )

    if len(hv_ids) < 2:
        log.error("  %s '%s': need ≥2 HV maps, found %d — skipping",
                  ctx.type_str, ctx.label, len(hv_ids))
        return None, None, float("nan"), None, np.empty((0,), dtype=bool)

    missing = [hv for hv in HV_SUBJECTS if hv not in hv_ids]
    if missing:
        log.warning("  %s '%s': no SBA map for %s", ctx.type_str, ctx.label, missing)

    log.info("  %s '%s': GLM on %d HVs  covariates=%s",
             ctx.type_str, ctx.label, len(hv_ids), cov_names or ["intercept only"])

    shape    = ref.shape
    age_mean = float(meta.loc[hv_ids, "age"].mean()) if "age" in cov_names else 0.0
    glm      = NormativeGLM(cov_names, age_mean)
    try:
        glm.fit(hv_ids, Y, meta)
    except RuntimeError as exc:
        log.error("  GLM failed for '%s': %s", ctx.label, exc)
        return None, None, float("nan"), None, np.empty((0,), dtype=bool)

    if not ctx.glm_mean_out().exists() or force:
        mean_map = np.where(brain_mask, glm.mean_at_hv_centroid(meta), 0.0).reshape(shape).astype(np.float32)
        sd_map   = np.where(brain_mask, np.sqrt(np.maximum(glm.sigma2_, MIN_SD ** 2)), 0.0).reshape(shape).astype(np.float32)
        nib.save(nib.Nifti1Image(mean_map, ref.affine), str(ctx.glm_mean_out()))
        nib.save(nib.Nifti1Image(sd_map,   ref.affine), str(ctx.resid_sd_out()))
        log.info("  Saved GLM mean + residual SD")

    log.info("  LOO calibration…")
    qc_df, loo_thr = glm.loo_analysis(meta)
    outliers = qc_df.loc[qc_df["outlier_flag"], "subject"].tolist()
    if outliers:
        log.warning("  LOO |z| > 2 for: %s", outliers)

    return glm, ref, loo_thr, qc_df, brain_mask


def compute_zscore_map(
    patient:       str,
    ctx:           NormCtx,
    glm:           NormativeGLM,
    ref_img:       "nib.Nifti1Image",
    meta:          pd.DataFrame,
    loo_thr:       float,
    brain_mask:    np.ndarray,
    positive_only: bool  = False,
    fc_threshold:  float = 0.0,
    force:         bool  = False,
) -> "nib.Nifti1Image | None":
    """GLM prediction-interval z-score + p-values + FDR + LOO significance maps."""
    out_z   = ctx.zscore_out(patient)
    out_p   = ctx.pval_out(patient)
    out_fdr = ctx.fdrsig_out(patient)
    out_loo = ctx.loosig_out(patient)

    if all(f.exists() for f in [out_z, out_p, out_fdr, out_loo]) and not force:
        log.info("  [%s | %s] Maps exist — loading", patient, ctx.label)
        return nib.load(str(out_z))

    pt_p = ctx.pt_sba(patient)
    if not pt_p.exists():
        log.error("  [%s | %s] No SBA map — run step1/step2 first", patient, ctx.label)
        return None

    y_obs = (
        resample_to(nib.load(str(pt_p)), ref_img, interpolation="linear")
        .get_fdata(dtype=np.float32).ravel().astype(np.float64)
    )
    if positive_only:
        y_obs = np.maximum(y_obs, 0.0)
    shape = ref_img.shape

    z_vec = glm.z_score(y_obs, patient, meta)
    z_vec = np.where(brain_mask, z_vec, 0.0)   # suppress interpolation bleed outside brain
    # Post-hoc FC threshold: zero z-scores at voxels where the HV group predicted
    # mean is below threshold.  Applied after statistics so the normative model is
    # estimated from unmodified data; only the output is masked.
    if fc_threshold > 0.0:
        hv_mean_pred = glm.mean_at_hv_centroid(meta).ravel()
        fc_mask = np.abs(hv_mean_pred) >= fc_threshold
        z_vec = np.where(fc_mask, z_vec, 0.0)
    p_vec = 2.0 * stats.t.sf(np.abs(z_vec), df=glm.df_)
    p_vec = np.where(brain_mask, p_vec, 1.0)
    if fc_threshold > 0.0:
        p_vec = np.where(fc_mask, p_vec, 1.0)

    z_img = nib.Nifti1Image(z_vec.reshape(shape).astype(np.float32), ref_img.affine)
    nib.save(z_img, str(out_z))
    nib.save(nib.Nifti1Image(p_vec.reshape(shape).astype(np.float32), ref_img.affine),
             str(out_p))

    q_vals  = _bh_fdr(p_vec)
    fdr_sig = (q_vals <= 0.05).reshape(shape)
    fdr_map = np.where(fdr_sig, z_vec.reshape(shape), 0.0).astype(np.float32)
    nib.save(nib.Nifti1Image(fdr_map, ref_img.affine), str(out_fdr))

    loo_sig = (np.abs(z_vec) > loo_thr).reshape(shape)
    loo_map = np.where(loo_sig, z_vec.reshape(shape), 0.0).astype(np.float32)
    nib.save(nib.Nifti1Image(loo_map, ref_img.affine), str(out_loo))

    ext     = glm.extrapolation_flags(patient, meta)
    flagged = [c for c, v in ext.items() if v]
    if flagged:
        log.warning("  [%s | %s] Extrapolation: %s outside HV range",
                    patient, ctx.label, flagged)

    log.info(
        "  [%s | %s] z-map (GLM)  FDR sig: %d vox  LOO sig: %d vox  (LOO thr=%.2f)",
        patient, ctx.label, int(fdr_sig.sum()), int(loo_sig.sum()), loo_thr,
    )
    return z_img


def roi_summary(
    patient:      str,
    ctx:          NormCtx,
    roi_mask_img: "nib.Nifti1Image",
    glm:          NormativeGLM,
    ref_img:      "nib.Nifti1Image",
    meta:         pd.DataFrame,
) -> list[dict]:
    """ROI-level GLM normative score — scalar GLM refit on ROI-averaged FC."""
    mask = resample_to(roi_mask_img, ref_img).get_fdata(dtype=np.float32).ravel() > 0
    if not mask.any():
        log.warning("  ROI '%s': empty mask after resampling", ctx.label)
        return []

    hv_ids_ok: list[str]   = []
    hv_fc:     list[float] = []
    for hv in glm.hv_ids_:
        p = ctx.hv_sba(hv)
        if p.exists():
            img = resample_to(nib.load(str(p)), ref_img, interpolation="linear")
            hv_fc.append(float(img.get_fdata(dtype=np.float32).ravel()[mask].mean()))
            hv_ids_ok.append(hv)

    if len(hv_fc) < 2:
        return []

    pt_p = ctx.pt_sba(patient)
    if not pt_p.exists():
        return []
    pt_fc = float(
        resample_to(nib.load(str(pt_p)), ref_img, interpolation="linear")
        .get_fdata(dtype=np.float32).ravel()[mask].mean()
    )

    Y_roi    = np.array(hv_fc, dtype=np.float64)
    age_mean = float(meta.loc[hv_ids_ok, "age"].mean()) if "age" in glm.cov_names else 0.0
    roi_glm  = NormativeGLM(glm.cov_names, age_mean)
    try:
        roi_glm.fit(hv_ids_ok, Y_roi, meta)
    except RuntimeError as exc:
        log.warning("  ROI GLM for '%s': %s", ctx.label, exc)
        return []

    roi_z  = float(roi_glm.z_score(pt_fc, patient, meta))
    t_crit = roi_glm.t_critical()
    p_val  = float(2.0 * stats.t.sf(abs(roi_z), df=roi_glm.df_))
    centroid_fc = float(roi_glm.mean_at_hv_centroid(meta))
    flags   = roi_glm.extrapolation_flags(patient, meta)
    ext_note = ", ".join(c for c, v in flags.items() if v) or "none"

    return [{
        "patient":          patient,
        ctx.type_str:       ctx.label,
        "patient_mean_FC":  round(pt_fc, 4),
        "HV_centroid_FC":   round(centroid_fc, 4),
        "ROI_z_score":      round(roi_z, 3),
        "p_value":          round(p_val, 4),
        "significant_p05":  p_val < 0.05,
        "t_critical":       round(t_crit, 2),
        "df":               roi_glm.df_,
        "covariates_used":  ", ".join(glm.cov_names) or "intercept only",
        "method":           "glm",
        "extrapolation":    ext_note,
        "n_HV":             len(hv_ids_ok),
        "interpretation":   _interpret_z(roi_z, t_crit),
    }]


def _interpret_z(z: float, t_crit: float) -> str:
    az   = abs(z)
    dir_ = "increased" if z > 0 else "decreased"
    if az < 1.0:
        return "within normal range"
    if az < t_crit:
        return f"mildly {dir_} connectivity (sub-threshold)"
    if az < t_crit * 1.5:
        return f"moderately {dir_} connectivity"
    return f"markedly {dir_} connectivity"


# ─────────────────────────────────────────────────────────────────────────────
# Coupling comparison
# ─────────────────────────────────────────────────────────────────────────────

def _coupling_raw_path(sub: str, profile: str) -> Path:
    return RSN_FC_DIR / f"sub-{sub}" / f"sub-{sub}_profile-{profile}_seed_rsn_coupling.csv"

def _coupling_zscore_path(sub: str, profile: str) -> Path:
    return COMP_DIR / f"sub-{sub}_profile-{profile}_seed_rsn_coupling_zscore.csv"

def _coupling_sig_path(sub: str, profile: str) -> Path:
    return COMP_DIR / f"sub-{sub}_profile-{profile}_seed_rsn_coupling_sig.csv"


def run_coupling_comparison(
    patients:     list[str],
    profile_name: str,
    force:        bool = False,
) -> None:
    """
    Element-wise normative comparison of seed–RSN coupling matrices.
    Loads HV matrices produced by step2 --profile, applies Crawford-Howell
    z-scores per element, calibrates an LOO empirical threshold (95th percentile
    of pooled |LOO z| across all folds and all matrix elements), saves z-score
    and LOO-significance CSVs.
    """
    hv_mats: dict[str, pd.DataFrame] = {}
    for hv in HV_SUBJECTS:
        p = _coupling_raw_path(hv, profile_name)
        if p.exists():
            hv_mats[hv] = pd.read_csv(str(p), index_col=0)

    if len(hv_mats) < 2:
        log.warning("Coupling: fewer than 2 HV matrices for profile '%s' — skipping", profile_name)
        return

    hv_ids     = list(hv_mats.keys())
    ref_df     = hv_mats[hv_ids[0]]
    seed_names = list(ref_df.index)
    rsn_names  = list(ref_df.columns)
    n_hv       = len(hv_ids)

    hv_stack = np.stack([hv_mats[hv].reindex(index=seed_names, columns=rsn_names).values
                         for hv in hv_ids], axis=0).astype(np.float64)  # (n_hv, n_s, n_r)
    hv_mu = np.nanmean(hv_stack, axis=0)
    hv_sd = np.maximum(np.nanstd(hv_stack, axis=0, ddof=1), 1e-4)

    log.info("Coupling normative comparison: %d HVs  |  %d seeds × %d RSNs",
             n_hv, len(seed_names), len(rsn_names))

    # LOO calibration: pool |z| across all HV folds and all matrix elements
    all_abs_z: list[np.ndarray] = []
    for i in range(n_hv):
        loo_mask   = np.ones(n_hv, dtype=bool)
        loo_mask[i] = False
        Y_loo  = hv_stack[loo_mask]
        mu_loo = np.nanmean(Y_loo, axis=0)
        sd_loo = np.maximum(np.nanstd(Y_loo, axis=0, ddof=1), 1e-4)
        z_loo  = (hv_stack[i] - mu_loo) / sd_loo
        all_abs_z.append(np.abs(z_loo).ravel())

    if all_abs_z:
        pooled  = np.concatenate(all_abs_z)
        loo_thr = float(np.percentile(pooled[np.isfinite(pooled)], 95.0))
        log.info("  Coupling LOO threshold (95th pctile): %.3f", loo_thr)
    else:
        loo_thr = 1.96
        log.warning("  Coupling LOO calibration failed — using z = 1.96")

    for patient in patients:
        out_z   = _coupling_zscore_path(patient, profile_name)
        out_sig = _coupling_sig_path(patient, profile_name)
        if out_z.exists() and out_sig.exists() and not force:
            log.info("  Coupling [%s]: exists — skipping", patient)
            continue

        pt_p = _coupling_raw_path(patient, profile_name)
        if not pt_p.exists():
            log.warning("  Coupling [%s]: no matrix found — run step2 --profile %s first",
                        patient, profile_name)
            continue

        pt_vals = (pd.read_csv(str(pt_p), index_col=0)
                   .reindex(index=seed_names, columns=rsn_names).values.astype(np.float64))

        z_mat   = (pt_vals - hv_mu) / hv_sd
        sig_mat = np.abs(z_mat) > loo_thr

        pd.DataFrame(z_mat,   index=seed_names, columns=rsn_names).to_csv(str(out_z),   float_format="%.4f")
        pd.DataFrame(sig_mat, index=seed_names, columns=rsn_names).to_csv(str(out_sig))

        log.info("  Coupling [%s]: %d/%d elements LOO-significant (|z| > %.2f)",
                 patient, int(sig_mat.sum()), int(sig_mat.size), loo_thr)


# ─────────────────────────────────────────────────────────────────────────────
# Generic normative loop
# ─────────────────────────────────────────────────────────────────────────────

def _run_normative_loop(
    ctx_list:      list[tuple[NormCtx, "nib.Nifti1Image | None"]],
    patients:      list[str],
    cov_names:     list[str],
    meta:          pd.DataFrame,
    roi_only:      bool,
    force:         bool,
    method:        str,
    positive_only: bool  = False,
    fc_threshold:  float = 0.0,
    r_thresh:      float = 0.25,
    thresholded:   bool  = False,
) -> tuple[list[pd.DataFrame], list[dict]]:
    qc_rows:  list[pd.DataFrame] = []
    roi_rows: list[dict]         = []

    z_thresh_fc = float(np.arctanh(np.clip(r_thresh, 0.0, 0.9999)))

    for ctx, roi_mask_img in ctx_list:
        log.info("")
        log.info("=== %s: %s  [method=%s] ===", ctx.type_str.upper(), ctx.label, method)

        if method == "simple":
            hv_mean, hv_sd, n_hv, ref_img, loo_thr, qc_df, brain_mask = build_simple_model(
                ctx, force=force, positive_only=positive_only,
            )
            if hv_mean is None:
                continue
        else:
            glm, ref_img, loo_thr, qc_df, brain_mask = build_normative_model(
                ctx, cov_names, meta, force=force, positive_only=positive_only,
            )
            if glm is None:
                continue

        qc_df["item"]   = ctx.label
        qc_df["type"]   = ctx.type_str
        qc_df["method"] = method
        qc_rows.append(qc_df)

        # HV network extents for volume Crawford-Howell (computed once per ctx)
        vox_vol = float(np.prod(np.abs(np.diag(ref_img.affine)[:3])))
        hv_extents: list[float] = []
        for hv in HV_SUBJECTS:
            p = ctx.hv_sba(hv)
            if p.exists():
                hv_extents.append(
                    _network_extent_mm3(nib.load(str(p)), z_thresh_fc, thresholded, vox_vol)
                )

        for patient in patients:
            log.info("  Patient: %s", patient)

            if not roi_only:
                if method == "simple":
                    compute_simple_zscore_map(
                        patient, ctx, hv_mean, hv_sd, n_hv, brain_mask, ref_img, loo_thr,
                        positive_only=positive_only, fc_threshold=fc_threshold, force=force,
                    )
                else:
                    compute_zscore_map(
                        patient, ctx, glm, ref_img, meta, loo_thr, brain_mask,
                        positive_only=positive_only, fc_threshold=fc_threshold, force=force,
                    )

            # ── Volume metrics ───────────────────────────────────────────────
            pt_fc_p  = ctx.pt_sba(patient)
            pt_ext   = (_network_extent_mm3(nib.load(str(pt_fc_p)), z_thresh_fc, thresholded, vox_vol)
                        if pt_fc_p.exists() else np.nan)

            if len(hv_extents) >= 2 and not np.isnan(pt_ext):
                ext_arr = np.array(hv_extents)
                n_e  = len(ext_arr)
                mu_e = float(ext_arr.mean())
                sd_e = max(float(ext_arr.std(ddof=1)), 1.0)
                z_e  = (pt_ext - mu_e) / sd_e
                p_e  = float(2.0 * stats.t.sf(abs(z_e * np.sqrt(n_e / (n_e + 1.0))), df=n_e - 1))
            else:
                mu_e = sd_e = z_e = p_e = np.nan

            fdr_p = ctx.fdrsig_out(patient)
            loo_p = ctx.loosig_out(patient)
            fdr_vol = _sig_volume_mm3(nib.load(str(fdr_p)), vox_vol) if fdr_p.exists() else np.nan
            loo_vol = _sig_volume_mm3(nib.load(str(loo_p)), vox_vol) if loo_p.exists() else np.nan

            vol_extras = {
                "patient_network_mm3": round(pt_ext)  if not np.isnan(pt_ext) else None,
                "HV_mean_network_mm3": round(mu_e)    if not np.isnan(mu_e)   else None,
                "network_extent_z":    round(z_e, 3)  if not np.isnan(z_e)    else None,
                "network_extent_p":    round(p_e, 4)  if not np.isnan(p_e)    else None,
                "network_extent_sig":  (p_e < 0.05)   if not np.isnan(p_e)    else None,
                "fdr_sig_mm3":         round(fdr_vol) if not np.isnan(fdr_vol) else None,
                "loo_sig_mm3":         round(loo_vol) if not np.isnan(loo_vol) else None,
            }
            # ────────────────────────────────────────────────────────────────

            if roi_mask_img is not None:
                if method == "simple":
                    rows = roi_summary_simple(
                        patient, ctx, roi_mask_img, hv_mean, hv_sd, n_hv, ref_img
                    )
                else:
                    rows = roi_summary(patient, ctx, roi_mask_img, glm, ref_img, meta)
                for row in rows:
                    row.update(vol_extras)
                roi_rows.extend(rows)
            elif any(v is not None for v in vol_extras.values()):
                roi_rows.append({
                    "patient":    patient,
                    ctx.type_str: ctx.label,
                    "method":     method,
                    **vol_extras,
                })

    return qc_rows, roi_rows


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--profile",    default=None,
                   help="Condition profile — required for seed mode")
    p.add_argument("--rsn",        action="store_true",
                   help="Run RSN normative comparison (uses step2 FC maps)")
    p.add_argument("--method",     choices=["simple", "glm"], default="simple",
                   help="z-score method: 'simple' = (pt−HV_mean)/HV_SD + Crawford-Howell "
                        "t-test (default); 'glm' = prediction-interval GLM with covariates")
    p.add_argument("--patients",   nargs="+", default=PT_SUBJECTS, metavar="ID")
    p.add_argument("--covariates", nargs="*", default=DEFAULT_COVARIATES, metavar="COV",
                   help="Covariates for --method glm only (default: age sex n_runs)")
    p.add_argument("--positive-only", action="store_true",
                   help="Constrain analysis to positive FC values only (negative FC set to 0)")
    p.add_argument("--fc-threshold", type=float, default=0.0, metavar="Z",
                   help="Minimum |Fisher-Z| to include in analysis; sub-threshold voxels set to 0 "
                        "(e.g. 0.1 ≈ r=0.10, 0.2 ≈ r=0.20; default: no threshold)")
    p.add_argument("--roi-only",   action="store_true",
                   help="Skip voxelwise maps; output ROI summary only")
    p.add_argument("--r-thresh",    type=float, default=0.25, metavar="R",
                   help="r threshold used to measure network extent volume (default: 0.25); "
                        "ignored when --thresholded (maps already binarised at this threshold)")
    p.add_argument("--thresholded", action="store_true",
                   help="Use thresholded+cluster-filtered FC maps from step1/step2 "
                        "(desc-sbaThresh) as input; outputs are also tagged *Thresh")
    p.add_argument("--force",      action="store_true",
                   help="Recompute and overwrite existing outputs")
    return p.parse_args()


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main() -> None:
    args = parse_args()

    if not args.profile and not args.rsn:
        log.error("Provide --profile, --rsn, or both.")
        sys.exit(1)

    patients  = list(args.patients)
    all_ids   = HV_SUBJECTS + patients
    meta      = load_metadata(all_ids)
    cov_names = (
        select_covariates(args.covariates or [], meta, HV_SUBJECTS)
        if args.method == "glm" else []
    )

    thr = args.thresholded

    if args.method == "simple":
        log.info("Method     : simple  (Crawford-Howell t-test, no covariate correction)")
    else:
        log.info("Method     : glm  (prediction-interval, covariates=%s)",
                 cov_names or ["intercept only"])
    log.info("Patients   : %s", patients)
    log.info("Thresholded: %s", thr)

    # Save covariate summary for GLM mode
    if args.method == "glm" and cov_names:
        cov_rows = []
        for cov in cov_names:
            hv_vals = meta.loc[HV_SUBJECTS, cov].values
            row = {
                "covariate": cov,
                "HV_mean":   round(float(hv_vals.mean()), 3),
                "HV_SD":     round(float(hv_vals.std()),  3),
                "HV_range":  f"{hv_vals.min():.2f}–{hv_vals.max():.2f}",
            }
            for pt in patients:
                row[f"{pt}_value"] = round(float(meta.loc[pt, cov]), 2)
            cov_rows.append(row)
        pd.DataFrame(cov_rows).to_csv(COMP_DIR / "normative_covariates.csv", index=False)

    all_qc_rows:  list[pd.DataFrame] = []
    all_roi_rows: list[dict]         = []

    # ── Seed mode ─────────────────────────────────────────────────────────────
    if args.profile:
        profile = get_profile(args.profile)
        seeds   = profile["seeds"]
        log.info("Profile : %s | Seeds: %s", args.profile, seeds)

        seed_ctx_list: list[tuple[NormCtx, nib.Nifti1Image | None]] = []
        for seed_name in seeds:
            seed_cfg = get_seed(seed_name)
            ctx = NormCtx(
                label        = seed_name,
                type_str     = "seed",
                hv_sba       = lambda hv, s=seed_name: sba_path(hv, s, thr),
                pt_sba       = lambda pt, s=seed_name: sba_path(pt, s, thr),
                zscore_out   = lambda pt, s=seed_name: zscore_path(pt, s, thr),
                pval_out     = lambda pt, s=seed_name: pval_path(pt, s, thr),
                fdrsig_out   = lambda pt, s=seed_name: fdrsig_path(pt, s, thr),
                loosig_out   = lambda pt, s=seed_name: loosig_path(pt, s, thr),
                glm_mean_out = lambda s=seed_name: glm_mean_path(s, thr),
                resid_sd_out = lambda s=seed_name: resid_sd_path(s, thr),
            )
            # Try atlas-based mask; fall back to first HV for subject-specific seeds
            roi_img = None
            for _subj in [None] + HV_SUBJECTS:
                try:
                    roi_img = get_seed_mask(seed_name, seed_cfg, subject=_subj)
                    break
                except Exception:
                    continue
            if roi_img is None:
                log.warning("  Seed '%s': could not build ROI mask", seed_name)
            seed_ctx_list.append((ctx, roi_img))

        qc, roi = _run_normative_loop(
            seed_ctx_list, patients, cov_names, meta, args.roi_only, args.force, args.method,
            positive_only=args.positive_only, fc_threshold=args.fc_threshold,
            r_thresh=args.r_thresh, thresholded=thr,
        )
        all_qc_rows.extend(qc)
        all_roi_rows.extend(roi)

        # ── Coupling comparison ───────────────────────────────────────────────
        run_coupling_comparison(patients, profile_name=args.profile, force=args.force)

    # ── RSN mode ──────────────────────────────────────────────────────────────
    if args.rsn:
        log.info("RSN normative comparison")

        yeo_img = nib.load(str(_YEO17_RS_PATH)) if _YEO17_RS_PATH.exists() else None
        if yeo_img is None:
            log.warning("Yeo17 group atlas not found at %s — RSN ROI masks disabled",
                        _YEO17_RS_PATH)

        available_rsns = [
            rsn for rsn in _ALL_RSN_NAMES
            if any(rsn_sba_path(hv, rsn, thr).exists() for hv in HV_SUBJECTS)
        ]
        if not available_rsns:
            log.error("No RSN SBA maps found — run step2 first")
        else:
            log.info("RSNs with HV data (%d): %s", len(available_rsns), available_rsns)

            rsn_ctx_list: list[tuple[NormCtx, nib.Nifti1Image | None]] = []
            for rsn_name in available_rsns:
                ctx = NormCtx(
                    label        = rsn_name,
                    type_str     = "rsn",
                    hv_sba       = lambda hv, r=rsn_name: rsn_sba_path(hv, r, thr),
                    pt_sba       = lambda pt, r=rsn_name: rsn_sba_path(pt, r, thr),
                    zscore_out   = lambda pt, r=rsn_name: rsn_zscore_path(pt, r, thr),
                    pval_out     = lambda pt, r=rsn_name: rsn_pval_path(pt, r, thr),
                    fdrsig_out   = lambda pt, r=rsn_name: rsn_fdrsig_path(pt, r, thr),
                    loosig_out   = lambda pt, r=rsn_name: rsn_loosig_path(pt, r, thr),
                    glm_mean_out = lambda r=rsn_name: rsn_glm_mean_path(r, thr),
                    resid_sd_out = lambda r=rsn_name: rsn_resid_sd_path(r, thr),
                )
                roi_mask = _get_rsn_roi_mask(rsn_name, yeo_img)
                rsn_ctx_list.append((ctx, roi_mask))

            qc, roi = _run_normative_loop(
                rsn_ctx_list, patients, cov_names, meta, args.roi_only, args.force, args.method,
                positive_only=args.positive_only, fc_threshold=args.fc_threshold,
                r_thresh=args.r_thresh, thresholded=thr,
            )
            all_qc_rows.extend(qc)
            all_roi_rows.extend(roi)

    # ── Save combined outputs ─────────────────────────────────────────────────
    if all_qc_rows:
        qc_path = COMP_DIR / "normative_QC.csv"
        pd.concat(all_qc_rows, ignore_index=True).to_csv(str(qc_path), index=False)
        log.info("LOO QC → %s", qc_path)

    if all_roi_rows:
        for patient in patients:
            pt_rows = [r for r in all_roi_rows if r["patient"] == patient]

            seed_rows = [r for r in pt_rows if r.get("seed")]
            if seed_rows and args.profile:
                out = COMP_DIR / f"sub-{patient}_profile-{args.profile}_roi_zscores.csv"
                pd.DataFrame(seed_rows).to_csv(str(out), index=False)
                log.info("Seed ROI z-scores → %s", out.name)
                _print_summary(patient, seed_rows, "seed")

            rsn_rows = [r for r in pt_rows if r.get("rsn")]
            if rsn_rows:
                out = COMP_DIR / f"sub-{patient}_rsn_roi_zscores.csv"
                pd.DataFrame(rsn_rows).to_csv(str(out), index=False)
                log.info("RSN ROI z-scores → %s", out.name)
                _print_summary(patient, rsn_rows, "rsn")

    log.info("")
    log.info("=" * 60)
    log.info("Normative comparison complete.")


def _get_rsn_roi_mask(rsn_name: str,
                      yeo_img: "nib.Nifti1Image | None") -> "nib.Nifti1Image | None":
    """
    Return a group-level binary ROI mask for one RSN.
    - Yeo17 networks  → extracted from resampled group Yeo17 atlas
    - Language        → Harvard-Oxford bilateral Broca/Wernicke mask
    - Auditory        → Smith RSN10 component 6 (z≥3), no subject GM constraint
    """
    if rsn_name == "Language":
        return get_language_group_mask()
    if rsn_name == "Auditory":
        return get_auditory_group_mask()
    if yeo_img is None:
        return None
    try:
        label_idx = YEO17_NETWORK_NAMES.index(rsn_name) + 1
    except ValueError:
        return None
    data = np.round(yeo_img.get_fdata()).astype(np.int32)
    mask = (data == label_idx).astype(np.uint8)
    return nib.Nifti1Image(mask, yeo_img.affine) if mask.any() else None


def _print_summary(patient: str, rows: list[dict], mode: str) -> None:
    if not rows:
        return
    # Volume-only rows (no ROI mask) lack ROI_z_score — only print rows that have it
    scored_rows = [r for r in rows if "ROI_z_score" in r]
    if not scored_rows:
        return
    t_crit = scored_rows[0].get("t_critical", float("nan"))
    df_val = scored_rows[0].get("df", "?")
    meth   = scored_rows[0].get("method", "?")
    key    = "seed" if mode == "seed" else "rsn"
    log.info("")
    log.info("  %s [%s | %s]   t_crit=%.2f  df=%s",
             patient, mode.upper(), meth, t_crit, df_val)
    log.info("  %-30s  %7s  %7s  %-6s  %s", "Label", "Z", "p", "Sig?", "Interpretation")
    log.info("  " + "─" * 72)
    for r in sorted(scored_rows, key=lambda x: abs(x["ROI_z_score"]), reverse=True):
        sig   = "yes *" if r.get("significant_p05") else "no"
        p_val = r.get("p_value", float("nan"))
        log.info("  %-30s  %+7.2f  %7.4f  %-6s  %s",
                 r.get(key, "?"), r["ROI_z_score"], p_val, sig, r.get("interpretation", ""))
    extrap = [r.get(key, "?") for r in rows
              if r.get("extrapolation") is not None and r.get("extrapolation") != "none"]
    if extrap:
        log.info("  Extrapolation flagged for: %s", extrap)


if __name__ == "__main__":
    main()
