# Conda → Miniforge Migration Guide

Per-env versions and install instructions. Full pinned dependency lists are in
`/home/aradwa0/conda_migration_exports/<env>_full.yml` — use those with
`mamba env create -f <env>_full.yml` for exact reproducibility, then follow the
editable-install fixup below where noted.

---

## scilpynew
- Python 3.12.13
- `scilpy` 2.3.0 — editable install, dev commit `b2bf4ac9` (176 commits past tag `2.2.2`), from `github.com/scilus/scilpy.git`

```bash
mamba env create -n scilpy -f /home/aradwa0/conda_migration_exports/scilpynew_full.yml
mamba activate scilpy
pip uninstall scilpy -y
git clone https://github.com/scilus/scilpy.git
cd scilpy && git checkout b2bf4ac95ab3dfbb622dfdb586988123ef88475e
pip install -e . --no-deps
```

---

## deepFCD
- Python 3.8.20 — **nothing else is currently installed**; this env was never finished
- Needed per upstream README (none present yet): `keras==2.2.4`, `theano==1.0.4`, `antspyx==0.4.2`, `antspynet==0.2.3`, `torch==1.8.2` (LTS), `h5py==2.10.0`, `pygpu==0.7.6`

```bash
mamba create -n deepFCD python=3.8 -y
mamba activate deepFCD
git clone --recurse-submodules -j2 https://github.com/NOEL-MNI/deepFCD.git
pip install keras==2.2.4 theano==1.0.4 antspyx==0.4.2 antspynet==0.2.3 h5py==2.10.0 pygpu==0.7.6
pip install torch==1.8.2+cu111 -f https://download.pytorch.org/whl/lts/1.8/torch_lts.html
```
Verify against the current `NOEL-MNI/deepFCD` README before running — these are old, largely unmaintained libraries (Theano was discontinued in 2017) and install steps may need adjusting for your CUDA version.

---

## hd-bet-env
- Python 3.10.16
- `hd-bet` 2.0.1 — editable install, commit `678e44d` (4 commits past tag `v2.0.1`, docs-only diff), from `github.com/MIC-DKFZ/HD-BET`
- `torch` 2.6.0, `nnunetv2` 2.5.2, plus ~110 more pip deps (see `hd-bet-env_full.yml`)
- **GPU note (only when NOT recreating from the full yml below):** if you instead
  do a fresh `mamba create -n hd-bet-env python=3.10` and then install `torch` and
  HD-BET as two separate `pip install` commands (no `--no-deps`, as
  `setup_environment.sh`'s `env-hdbet` section used to), that's a confirmed trap —
  HD-BET/nnunetv2's unconstrained `torch` requirement gets re-resolved against plain
  PyPI on the second command and silently replaces a working CUDA build with a
  newer, incompatible one (e.g. a driver capped at CUDA 12.2 ends up with an
  unusable cu130-tagged torch). Fix: install both together in one resolver pass,
  pinned explicitly so it can't drift —
  `pip install --index-url https://download.pytorch.org/whl/cu121 --extra-index-url https://pypi.org/simple "torch==2.5.1+cu121" -e .`
  (no `--no-deps` there — it needs to resolve nnunetv2's other ~110 deps too since
  there's no full yml pre-populating them; 2.5.1+cu121 was the latest version
  actually published on the cu121 index as of this writing). The recreation-from-yml
  path below doesn't hit this, since `--no-deps` there skips re-resolving `torch`
  entirely — the full yml already pins a working build.

```bash
mamba env create -n hd-bet-env -f /home/aradwa0/conda_migration_exports/hd-bet-env_full.yml
mamba activate hd-bet-env
pip uninstall hd-bet -y
git clone https://github.com/MIC-DKFZ/HD-BET.git
cd HD-BET && git checkout 678e44d546a84de0f2a7fc245f176b82b7d912fd
pip install -e . --no-deps
```

---

## KarawunDev — **the default since Aug 2026**
- Python 3.14
- `karawun` — editable install of the KUL fork, `github.com/Rad-dude/karawun.git`,
  branch `kul-extended-palette`, commit `72bbc955`
- `dcm2bids` >= 3.1

```bash
mamba create -n KarawunDev python=3.14 -y
mamba activate KarawunDev
pip install 'dcm2bids>=3.1'
git clone -b kul-extended-palette https://github.com/Rad-dude/karawun.git
cd karawun && git checkout 72bbc9558c139678a423b2ee43fbb134d656328c
pip install -e .
```

Despite the name this was never only for developing karawun. It is the default
because the conda-forge alternative is wrong in two independent ways — see
KarawunEnv below. The fork keeps upstream's packaging unchanged (`setup.py`,
`pyproject.toml`, `setup.cfg`, `versioneer.py`), so the editable install is the
same as it always was; only the remote and branch differ.

Previously: `treanus/karawun.git` @ `80ea5cf`, then upstream
`DevelopmentalImagingMCRI/karawun` @ `9ecbf8e6` (tag `v0.2.6.0`).

---

## KarawunEnv — legacy, `--karawun-stock` only
- Python 3.8.17
- `karawun` 0.2.5.4 — regular conda-forge package, **no git clone needed**

```bash
mamba create -n KarawunEnv python=3.8 karawun=0.2.5.4 -c conda-forge -y
```

No longer the default. Two independent problems, neither fixable from
conda-forge (there is no newer release there):

1. **Corrupts Enhanced/multi-frame DICOM output.** Given a Philips (or any)
   Enhanced donor, leftover `NumberOfFrames` /
   `PerFrameFunctionalGroupsSequence` tags survive into what should be
   single-frame slices. Confirmed on real patient data.
2. **31-entry colour palette.** `lookup_cie()` clamps every label index above
   30 to the last entry. `KUL_karawun_prepare.sh` assigns 1-41 to known tracts,
   42-49 to auto tracts, 50 to the lesion and 51-63 to fMRI activations — so
   with stock karawun every activation renders in one colour, indistinguishable
   from each other and colliding with the tracts.

---

## resseg
- Python 3.8.17
- `resseg` 0.3.7, `antspyx` 0.4.2 — regular PyPI packages, **no git clone needed**
  (originally pinned to antspyx 0.3.8, since removed from PyPI; 0.4.2 is the oldest still available)

```bash
mamba create -n resseg python=3.8 -y
mamba activate resseg
pip install resseg==0.3.7 antspyx==0.4.2
```
Use `resseg_full.yml` instead if you want the exact 2023-era pinned versions of matplotlib/torch/etc. rather than latest-compatible.

---

## FastSurfer
- Not a conda env — FastSurfer's native install upstream moved to `uv` + a
  per-checkout `.venv`. `env/environment_gpu.yml` (what an older version of
  this guide/script assumed) no longer exists in the current repo; there's
  only an empty `env/fastsurfer.yml` placeholder. Confirmed via their current
  `doc/overview/INSTALL.md`.
- `uv` itself installed as a plain pip package into the isolated miniforge3
  base env (`mamba/pip install uv`) — self-contained, no separate download or
  PATH entry needed.
- `requirements.txt` deliberately excludes nvidia-* CUDA runtime packages and
  torch itself from its pinned list (per its own header comment) — it must be
  resolved with `uv pip compile --torch-backend <backend> requirements.txt |
  uv pip sync --torch-backend <backend> --python .venv/bin/python -`, **not**
  a direct `uv pip sync requirements.txt`. Confirmed: `sync` alone only
  installs packages literally listed in the file, silently skipping the CUDA
  runtime libs and leaving `torch` unable to import
  (`libcudart.so.11.0: cannot open shared object file`).
- Backend: `cu118` — `requirements.txt` currently pins `torch==2.7.1`, which
  has no `+cu121` build published but does have `+cu118`; this machine's
  driver caps at CUDA 12.2, well within `cu118`'s compatibility range. Re-check
  both (`nvidia-smi`, `pip index versions torch --index-url
  https://download.pytorch.org/whl/<backend>`) if `requirements.txt`'s pin
  changes later.
- **Wrapper shim required**: `run_fastsurfer.sh` calls plain `python3` with no
  self-activation of its own venv, and it's invoked both by full path
  (`$FASTSURFER_HOME/run_fastsurfer.sh` in `KUL_anat_segment_tumor.sh`,
  `KUL_T1T2FLAIRMTR_ratio.sh`) and by bare name via `PATH` (`KUL_FS_multiparc.sh`,
  `KUL_clinical_fmridti.sh`, `KUL_VBG.sh` — 7 call sites total across both
  repos). Renamed the real script to `run_fastsurfer.real.sh` and put a small
  wrapper at the original `run_fastsurfer.sh` path that sources `.venv/bin/activate`
  then `exec`s the real script — covers every call site with zero edits to
  KUL_NIS/KUL_VBG source, and the activation is scoped to that one process
  (replaced by `exec`) so it never leaks into the calling shell.

```bash
git clone https://github.com/Deep-MI/FastSurfer.git
cd FastSurfer
uv venv --python 3.12
uv pip compile --torch-backend cu118 requirements.txt --python .venv/bin/python \
  | uv pip sync --torch-backend cu118 --python .venv/bin/python -
mv run_fastsurfer.sh run_fastsurfer.real.sh
# then write the wrapper shim at run_fastsurfer.sh — see setup_environment.sh's
# section_env_fastsurfer for the exact contents
```

---

## rsfmri_env
- Python 3.10, conda-forge: `nilearn`, `nibabel`, `numpy`, `scipy`, `pandas`, `matplotlib`, `pyyaml`
- Not present in any export — this env was previously undocumented/missing entirely.
  `KUL_run_rsfMRI_networks.sh` (`KUL_NIS/share/rsfmri_pipeline`) takes the env name via
  its own `-c <conda_env>` flag (no hardcoded default), the same pattern scilpy uses.
- Installed as a single conda-forge resolver pass (not piecemeal pip installs) so numpy/
  pandas/scipy/matplotlib end up as mutually ABI-compatible builds — a colleague's
  similarly-purposed personal env (`rob_nilearn`, pip-installed incrementally over time)
  hit a numpy 2.x / `bottleneck` ABI mismatch that breaks `pandas` on import; this
  avoids that by construction.

```bash
mamba create -n rsfmri_env -c conda-forge python=3.10 nilearn nibabel numpy scipy pandas matplotlib pyyaml -y
```

---

## lore_sd
- Python 3.10, `lore_sd` 1.0 — editable install from `github.com/SiebeLeysen/LoRE-SD.git`, `mrtrix_module` branch
- Despite the branch name and older docs suggesting an MRtrix3 external-module
  build (`./build -external <path>` against an MRtrix3 source tree), LoRE-SD is
  actually just a regular pip package (`setup.py` + pybind11/nlopt/cmake) —
  confirmed working via the sequence below. No MRtrix3 source tree needed.
- Used by `KUL_dwiprep -D run_dwiprep_lore_sd.txt`.
- Validated on the shared/legacy layout at `/mnt/DATA1/aradwa0/local_KUL_NIS/LoRE-SD_upstream`
  (editable install pointing outside any SOFTWARE_ROOT) — works fine there, but for
  a clean install put the clone under `$SOFTWARE_ROOT/src/LoRE-SD` instead, matching
  every other repo here (editable installs point straight at the source, so keep it
  inside SOFTWARE_ROOT so nothing breaks if that location ever moves).

```bash
mamba create -n lore_sd python=3.10 -y
mamba activate lore_sd
mamba install -c conda-forge nlopt cmake compilers -y
pip install pybind11 numpy setuptools wheel
git clone --branch mrtrix_module https://github.com/SiebeLeysen/LoRE-SD.git
cd LoRE-SD
pip install -e . --no-build-isolation
```
