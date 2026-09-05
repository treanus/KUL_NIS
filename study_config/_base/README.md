# `study_config/_base/`

The default config set. `KUL_clinical_fmridti.sh -s` copies **`_base/` first**,
then copies the `-t`-specific directory over it, so a type directory needs to
contain only the files it genuinely changes.

```
_base/                     -> copied for every type
clinical_fmri_dmri/        -> types 1-4  (empty: base is already the fmri/dmri set)
clinical_dmri_dbs_drt/     -> type 5
clinical_dmri_dbs_hdp/     -> type 6
DTI_ALPS_proc/             -> type 7
```

Overriding is at **file** granularity, not line granularity: a type directory
that ships `run_fmriprep.txt` replaces the base file wholesale.

Before this layout, all four type directories carried a full copy of every
config, and they had drifted -- `run_freesurfer.txt` existed in three versions
whose only real difference was a hardcoded `freesurfer_ncpu`. Per-step `*_ncpu`
values are now written by `KUL_clinical_fmridti.sh` from its own `-n`, so they
are no longer a reason for a file to differ per type.

One file per sub-command, matching the flags of the script it drives:

| file | drives |
|---|---|
| `sequences.txt` | `KUL_dcm2bids.sh` |
| `run_fmriprep.txt`, `run_freesurfer.txt`, `run_dwiprep.txt` | `KUL_preproc_all.sh` |
| `run_dwiprep_lore_sd.txt` | `KUL_preproc_all.sh` (lore-sd FOD variant) |
| `run_vbg.txt` | `KUL_VBG.sh` |
| `run_fwt.txt` + `tracks_list.txt` | `KUL_FWT_make_VOIs.sh`, `KUL_FWT_make_TCKs.sh` |
| `run_fmri_glm.txt` | `KUL_fmriproc_spm_new.sh` / `KUL_fmriproc_nilearn_new.sh` |
| `run_rsfmri_networks.txt` | `KUL_run_rsfMRI_networks.sh` |
| `run_dsc.txt` | `KUL_dsc_perfusion.sh` |
| `run_multiparc.txt` | `KUL_FS_multiparc.sh` |
| `run_karawun.txt` | `KUL_karawun_prepare.sh` |
