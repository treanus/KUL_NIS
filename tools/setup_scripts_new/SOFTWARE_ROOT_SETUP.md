# Self-contained software root — setup guide

Goal: get everything (miniforge, conda envs, source repos, neuroimaging tools,
caches, tmp) off the shared `/usr/local/KUL_apps` and `/usr/local/anaconda3`
mess and onto your own drive, under one root:

```
/mnt/DATA1/aradwa0/software/
```

This is a **documentation-only** guide — nothing below has been executed. It's
written so you can run it yourself at your own pace, and so the exact
versions/commits/branches currently in use (many of which are deliberately
*not* the latest upstream) are recorded before they get lost.

---

## 1. Directory layout

```
/mnt/DATA1/aradwa0/software/
├── miniforge3/              # conda/mamba install; envs live in miniforge3/envs/
├── cache/
│   ├── pip/                 # PIP_CACHE_DIR
│   ├── conda_pkgs/          # CONDA_PKGS_DIRS
│   └── xdg/                 # XDG_CACHE_HOME
├── tmp/                     # TMPDIR, APPTAINER_TMPDIR
├── src/                     # all git-cloned source code, one dir per repo
│   ├── KUL_NIS/
│   ├── KUL_VBG/
│   ├── KUL_FWT/
│   ├── LoRE-SD/
│   ├── scilpy/
│   ├── HD-BET/
│   ├── karawun/
│   ├── mrtrix3/             # source + build tree
│   └── ANTs/                # source tree (build output goes to ANTs/../ANTs_install)
├── FSL/                     # FSLDIR — self-contained FSL install
├── freesurfer/              # FREESURFER_HOME — move your existing 8.2.0 copy here
└── ANTs_install/            # ANTSPATH — ANTs compiled binaries (separate from src/ANTs)
```

---

## 2. Environment variables (`~/.bashrc`)

Replace the current scattered exports (pointing at `/usr/local/fsl`,
`/usr/local/KUL_apps/...`, `/mnt/DATA1/aradwa0/FS8.2.0/...`, etc.) with:

```bash
export SOFTWARE_ROOT=/mnt/DATA1/aradwa0/software

# --- Conda / miniforge ---
export PATH="$SOFTWARE_ROOT/miniforge3/bin:$PATH"
source "$SOFTWARE_ROOT/miniforge3/etc/profile.d/conda.sh"
export CONDA_PKGS_DIRS="$SOFTWARE_ROOT/cache/conda_pkgs"

# --- Caches / tmp (already partly in place today, just re-rooted) ---
export XDG_CACHE_HOME="$SOFTWARE_ROOT/cache/xdg"
export PIP_CACHE_DIR="$SOFTWARE_ROOT/cache/pip"
export TMPDIR="$SOFTWARE_ROOT/tmp"
export APPTAINER_TMPDIR="$SOFTWARE_ROOT/tmp/apptainer"

# --- FSL ---
export FSLDIR="$SOFTWARE_ROOT/FSL"
source "$FSLDIR/etc/fslconf/fsl.sh"
export PATH="$FSLDIR/share/fsl/bin:$PATH"

# --- FreeSurfer ---
export FREESURFER_HOME="$SOFTWARE_ROOT/freesurfer"
source "$FREESURFER_HOME/SetUpFreeSurfer.sh"
export SUBJECTS_DIR="$FREESURFER_HOME/subjects"
export FS_LICENSE="$FREESURFER_HOME/license.txt"

# --- ANTs ---
export ANTSPATH="$SOFTWARE_ROOT/ANTs_install/bin/"
export PATH="$ANTSPATH:$PATH"

# --- mrtrix3 ---
export PATH="$SOFTWARE_ROOT/src/mrtrix3/bin:$PATH"
export PYTHONPATH="$SOFTWARE_ROOT/src/mrtrix3/lib:$PYTHONPATH"
```

---

## 3. Miniforge

```bash
mkdir -p /mnt/DATA1/aradwa0/software
bash Miniforge3-Linux-x86_64.sh -b -p /mnt/DATA1/aradwa0/software/miniforge3
```
`.condarc` (conda-forge only, matches what you already have configured):
```yaml
channels:
  - conda-forge
mirrored_channels:
  conda-forge:
    - https://conda.anaconda.org/conda-forge
    - https://prefix.dev/conda-forge
```

---

## 4. Source repos — exact pinned commits currently in use

**None of these are "latest main/master"** — each is pinned to a specific
commit that is known-good for this pipeline. Cloning `main`/`dev` fresh will
likely give you a *different and possibly incompatible* version.

| Repo | Remote | Branch | Commit | Notes |
|---|---|---|---|---|
| KUL_NIS | `treanus/KUL_NIS.git` | `KUL_NIS_v2.0_20260701` | `ba4009d6` | Not `main` — this dev branch |
| KUL_VBG | `KUL-Radneuron/KUL_VBG.git` | `KUL_VBG_2.0` | `b25ac249` | v2.0 rewrite branch |
| KUL_FWT | `KUL-Radneuron/KUL_FWT.git` | `KUL_FWT_v2.0` | `c4accb8c` | v2.0 rewrite branch |
| LoRE-SD | `SiebeLeysen/LoRE-SD.git` | `mrtrix_module` | tracks branch tip (was `4d259c7e` at last check, now past that — not commit-pinned like the rows below) | **Not `main`**. Despite the branch name, it's installed as a regular editable pip package (`setup.py` + pybind11/nlopt/cmake in its own `lore_sd` conda env), not an MRtrix3 external-module build — see `MIGRATION_GUIDE.md`'s `lore_sd` entry |
| scilpy | `scilus/scilpy.git` | (detached, dev) | `b2bf4ac9` | 176 commits past tag `2.2.2` — pip installed editable |
| HD-BET | `MIC-DKFZ/HD-BET` | `master` | `678e44d5` | 4 commits past tag `v2.0.1` (docs only, harmless) |
| karawun | `treanus/karawun.git` | (detached, dev) | `80ea5cf0` | 4 commits past tag `v0.2.5.4` — unreleased bugfix |

```bash
cd /mnt/DATA1/aradwa0/software/src

git clone https://github.com/treanus/KUL_NIS.git && cd KUL_NIS && git checkout KUL_NIS_v2.0_20260701 && cd ..
git clone https://github.com/KUL-Radneuron/KUL_VBG.git && cd KUL_VBG && git checkout KUL_VBG_2.0 && cd ..
git clone https://github.com/KUL-Radneuron/KUL_FWT.git && cd KUL_FWT && git checkout KUL_FWT_v2.0 && cd ..
git clone https://github.com/SiebeLeysen/LoRE-SD.git && cd LoRE-SD && git checkout mrtrix_module && cd ..
git clone https://github.com/scilus/scilpy.git && cd scilpy && git checkout b2bf4ac95ab3dfbb622dfdb586988123ef88475e && cd ..
git clone https://github.com/MIC-DKFZ/HD-BET.git && cd HD-BET && git checkout 678e44d546a84de0f2a7fc245f176b82b7d912fd && cd ..
git clone https://github.com/treanus/karawun.git && cd karawun && git checkout 80ea5cf06910c5c58732542e80822f0e47f49dc0 && cd ..
```

For the conda envs that use these editable (scilpy, HD-BET, karawun-dev), see
`MIGRATION_GUIDE.md` for the exact `mamba env create` + `pip install -e .
--no-deps` sequence — just point the clone destination at
`$SOFTWARE_ROOT/src/...` instead of an arbitrary path.

---

## 5. Neuroimaging tools — document only for now (still using the shared installs today)

### mrtrix3 — ⚠️ pin carefully, latest `dev` is currently buggy

Currently running: **`dev` branch, commit `5a3a8bf6` (`3.0.4-537-g5a3a8bf6`),
built 2023-07-14.** This is a deliberately old, verified-good checkout — do
**not** just `git pull` the `dev` branch to get something newer; per your own
notes the current tip of `dev` has known regressions.

Note: `KUL_NIS/README.md` documents the requirement as `3.0.4-543-g86eb1ea8`
— 6 commits newer than what's actually installed at
`/usr/local/KUL_apps/mrtrix3`. Worth reconciling which one is really
"good" before pinning the new install — they're close but not identical.

```bash
cd /mnt/DATA1/aradwa0/software/src
git clone https://github.com/MRtrix3/mrtrix3.git
cd mrtrix3
git checkout 5a3a8bf6143cd0462f9ba8a5fbcda208970339aa   # known-good dev build in use today
./configure
./build
```

### ANTs

Currently running: **`master` branch, commit `40ee2d22` (`v2.4.4.post20`),
built 2023-06-29.**

```bash
cd /mnt/DATA1/aradwa0/software/src
git clone https://github.com/ANTsX/ANTs.git
cd ANTs
git checkout 40ee2d221aecbd4694283426cf67e34262acf1da
mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=/mnt/DATA1/aradwa0/software/ANTs_install ..
make -j$(nproc)
cd ANTS-build && make install
```

### FSL

Currently running: **6.0.6.5** (shared install at `/usr/local/fsl`). FSL is
distributed via its own licensed installer, not conda/git:

```bash
python fslinstaller.py -d /mnt/DATA1/aradwa0/software/FSL -V 6.0.6.5
```
(Check `https://fsl.fmrib.ox.ac.uk/fsldownloads` for the current installer —
syntax has changed across FSL releases.)

### FreeSurfer

You already have a personal copy — no reinstall needed, just relocate it:

```bash
mkdir -p /mnt/DATA1/aradwa0/software
cp -a /mnt/DATA1/aradwa0/FS8.2.0/usr/local/freesurfer/8.2.0 /mnt/DATA1/aradwa0/software/freesurfer
```
Version: `freesurfer-linux-ubuntu22_x86_64-8.2.0-20260314-d932c45`.

---

## 6. dcm2bids — version note

KUL_NIS requires **dcm2bids ≥ 3.0** (v3 schema); a v2 fallback script
(`KUL_dcm2bids_v2bkup.sh`) exists for the old `custom_entities` /
`sidecar_changes` config schema, but v3 is the maintained path. Your envs
currently disagree on exact patch version — KarawunDev has `3.1.1`,
hd-bet-env has `3.2.0`. Any recent 3.x satisfies KUL_NIS's requirement; pick
one (e.g. latest 3.x) and use it consistently across envs that need it,
since config file schema only changes at the v2→v3 boundary, not between 3.x
patches.

---

## 7. Order of operations (suggested)

1. Install miniforge3 under the new root (§3).
2. Recreate the 6 conda envs per `MIGRATION_GUIDE.md`, using `$SOFTWARE_ROOT/src/...`
   as the clone destination for editable installs.
3. Clone the KUL_NIS / KUL_VBG / KUL_FWT / LoRE-SD repos at their pinned
   branches (§4).
4. Update `~/.bashrc` to the new env-var block (§2) — do this **last**, once
   everything above exists, so your shell doesn't break mid-migration.
5. FSL / ANTs / mrtrix3 / FreeSurfer: leave pointing at the shared installs
   until you're ready to build the self-contained copies in §5; nothing
   depends on doing that immediately.
