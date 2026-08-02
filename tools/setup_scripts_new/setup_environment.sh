#!/usr/bin/env bash
#
# setup_environment.sh — set up a fresh Ubuntu/Linux Mint box to run
# KUL_NIS + KUL_VBG + KUL_FWT.
#
# Sources this script is built from (all in this same directory):
#   - KUL_NIS/README.md          (authoritative "Requirements" table)
#   - KUL_VBG/README.md          ("Updated Dependencies" section)
#   - KUL_FWT/README.md          ("Prerequisites" section)
#   - SOFTWARE_ROOT_SETUP.md     (exact pinned commits/versions currently in use)
#   - MIGRATION_GUIDE.md         (per-conda-env install recipes)
#
# This script is idempotent where practical (safe to re-run; skips what's
# already present) but it DOES modify the system: it installs apt packages,
# Docker, and writes a block to ~/.bashrc. Read section toggles below before
# running on a machine you care about. Everything under $SOFTWARE_ROOT is
# self-contained and easy to delete if you want to start over; the apt/Docker
# install and the ~/.bashrc edit are the two things that touch the system
# outside that directory.
#
# Usage:
#   ./setup_environment.sh                 # run every enabled section, in order
#   ./setup_environment.sh --list          # show section names and exit
#   ./setup_environment.sh --only mrtrix3,ants   # run just these sections
#   ./setup_environment.sh --skip docker,fsl     # run everything except these
#   ./setup_environment.sh --dry-run       # print what would run, do nothing
#
# Re-run safely after fixing a failure partway through — completed sections
# are cheap to skip (they check for existing output first).

set -euo pipefail

# ── Configuration — edit these before running ────────────────────────────────

SOFTWARE_ROOT="${SOFTWARE_ROOT:-/mnt/DATA1/$(whoami)/software}"
NCPU="${NCPU:-$(nproc)}"

# Defend against inherited conda/mamba state. Real bug, found in practice on a
# machine with a pre-existing micromamba install: MAMBA_ROOT_PREFIX (exported by
# micromamba's own 'shell init' block in ~/.bashrc) silently redirects
# `mamba create -n <name>` to THAT root instead of $SOFTWARE_ROOT/miniforge3 --
# even when the miniforge3 mamba binary is invoked by its full path. The env
# then gets created in the wrong place, and a later `conda run -n <name>` (using
# miniforge3's own conda, which correctly looks under $SOFTWARE_ROOT) fails with
# EnvironmentLocationNotFound -- which, combined with `set -e`, silently aborted
# the whole script before any later env-creation sections ever ran. A leftover
# activated env from an unrelated conda install (CONDA_PREFIX/VIRTUAL_ENV/compiler
# wrapper vars all pointing elsewhere) can cause similar surprises. Pin the one
# variable directly responsible for the confirmed bug; if you hit further
# contamination from your own shell's conda/mamba setup, the robust fix is
# running this script from a fresh, non-activated shell.
export MAMBA_ROOT_PREFIX="$SOFTWARE_ROOT/miniforge3"

# Pinned versions/commits — see SOFTWARE_ROOT_SETUP.md for how these were
# determined. These are the versions this pipeline is actually validated
# against, not necessarily the latest upstream.
MRTRIX3_COMMIT="86eb1ea8"                 # KUL_NIS/README.md's documented requirement
                                           # (3.0.4-543-g86eb1ea8). SOFTWARE_ROOT_SETUP.md
                                           # notes the machine it audited was actually on
                                           # 5a3a8bf6 (6 commits older) — if this commit
                                           # fails to build cleanly, fall back to that one.
ANTS_COMMIT="40ee2d22"                    # v2.4.4.post20
FSL_VERSION="6.0.6.5"
SCILPY_COMMIT="b2bf4ac95ab3dfbb622dfdb586988123ef88475e"   # 176 commits past tag 2.2.2
HDBET_COMMIT="678e44d546a84de0f2a7fc245f176b82b7d912fd"    # 4 commits past tag v2.0.1
KARAWUN_COMMIT="80ea5cf06910c5c58732542e80822f0e47f49dc0"  # 4 commits past tag v0.2.5.4
FMRIPREP_VERSION="25.1.4"

KUL_NIS_BRANCH="KUL_NIS_v2.0_20260701"
KUL_VBG_BRANCH="KUL_VBG_2.0"
KUL_FWT_BRANCH="KUL_FWT_v2.0"

# Section toggles — set any of these to 0 to skip by default, or use
# --only/--skip on the command line instead of editing this file.
DO_APT=1
DO_DOCKER=1
DO_MINIFORGE=1
DO_ENV_SCILPY=1          # required — KUL_FWT's own filtering/RecoBundles step
DO_ENV_HDBET=1           # brain extraction, used by KUL_dwiprep/KUL_anat_register
DO_ENV_RESSEG=1          # resection-cavity segmentation
DO_ENV_KARAWUN=1         # Brainlab export (uses the simpler conda-forge KarawunEnv;
                         # see --karawun-dev below for the editable-install variant)
DO_ENV_FASTSURFER=1      # needs a GPU to be useful; safe to leave on, just slow/CPU-only without one
DO_ENV_RSFMRI=1          # KUL_NIS/share/rsfmri_pipeline (KUL_run_rsfMRI_networks.sh -c <env>)
DO_REPOS=1               # clone KUL_NIS/KUL_VBG/KUL_FWT at pinned branches
DO_ENV_LORE_SD=1         # LoRE-SD (KUL_dwiprep -D run_dwiprep_lore_sd.txt) — regular pip package
DO_MRTRIX3=1
DO_ANTS=1
DO_FSL=1
DO_FREESURFER_CHECK=1    # can't auto-install (license-gated); just verifies/guides
DO_DOCKER_IMAGES=1       # pulls fmriprep/mriqc/synb0-disco/MSBP images
DO_BASHRC=1
DO_VERIFY=1              # final read-only health check + summary table; safe to
                         # run on its own any time, doesn't modify anything

USE_KARAWUN_DEV=0        # 1 = editable git install (KarawunDev) instead of the plain
                         # conda-forge package (KarawunEnv) — only needed if you're
                         # developing karawun itself.
INCLUDE_HDGLIOAUTO=0     # off by default: GPU-heavy, narrower use (tumor auto-seg),
                         # and its install process is the least standardized of
                         # everything here — read HD-GLIO-AUTO's own README before
                         # enabling this section.

DRY_RUN=0

# ── End configuration ─────────────────────────────────────────────────────────

SCRIPT_SECTIONS="apt docker miniforge env-scilpy env-hdbet env-resseg env-karawun env-fastsurfer env-rsfmri env-lore-sd repos mrtrix3 ants fsl freesurfer-check docker-images bashrc verify"

# ── Helpers ────────────────────────────────────────────────────────────────────

c_blue=$'\033[1;34m'; c_green=$'\033[1;32m'; c_yellow=$'\033[1;33m'; c_red=$'\033[1;31m'; c_reset=$'\033[0m'

log()   { echo "${c_blue}==>${c_reset} $*"; }
ok()    { echo "${c_green}  ✓${c_reset} $*"; }
warn()  { echo "${c_yellow}  !${c_reset} $*"; }
fail()  { echo "${c_red}  ✗${c_reset} $*" >&2; }
run()   { if [ "$DRY_RUN" -eq 1 ]; then echo "    [dry-run] $*"; else eval "$@"; fi; }

have()  { command -v "$1" >/dev/null 2>&1; }

section_enabled() {
    local name="$1"
    if [ -n "${ONLY_SECTIONS:-}" ]; then
        [[ ",${ONLY_SECTIONS}," == *",${name},"* ]]
    elif [ -n "${SKIP_SECTIONS:-}" ]; then
        [[ ",${SKIP_SECTIONS}," != *",${name},"* ]]
    else
        return 0
    fi
}

# ── Argument parsing ───────────────────────────────────────────────────────────

while [ $# -gt 0 ]; do
    case "$1" in
        --list)
            echo "$SCRIPT_SECTIONS" | tr ' ' '\n'
            exit 0
            ;;
        --only)      ONLY_SECTIONS="$2"; shift 2 ;;
        --skip)      SKIP_SECTIONS="$2"; shift 2 ;;
        --dry-run)   DRY_RUN=1; shift ;;
        --karawun-dev) USE_KARAWUN_DEV=1; shift ;;
        --with-hdglioauto) INCLUDE_HDGLIOAUTO=1; shift ;;
        -h|--help)
            grep '^#' "$0" | sed 's/^# \{0,1\}//'
            exit 0
            ;;
        *) fail "Unknown argument: $1"; exit 1 ;;
    esac
done

log "Target software root: $SOFTWARE_ROOT"
log "CPUs available for builds: $NCPU"
[ "$DRY_RUN" -eq 1 ] && warn "DRY RUN — no commands will actually be executed"

if [ "$DRY_RUN" -eq 1 ]; then
    echo "    [dry-run] would create the $SOFTWARE_ROOT directory skeleton (src, cache, tmp, FSL, freesurfer, ANTs_install)"
else
    mkdir -p "$SOFTWARE_ROOT"/{src,cache/pip,cache/conda_pkgs,cache/xdg,tmp/apptainer,FSL,freesurfer,ANTs_install}
fi

# ── 1. System packages ────────────────────────────────────────────────────────

section_apt() {
    log "Installing system packages (apt)"
    run "sudo apt-get update"
    run "sudo DEBIAN_FRONTEND=noninteractive apt-get install -y \
        build-essential cmake ninja-build git curl wget ca-certificates gnupg unzip \
        python3 python3-dev python3-pip python3-venv \
        libeigen3-dev zlib1g-dev libfftw3-dev libtiff5-dev libpng-dev \
        libqt5opengl5-dev libqt5svg5-dev libgl1-mesa-dev \
        xvfb \
        p7zip-full \
        dcmtk"
    ok "apt packages installed (mrtrix3/ANTs build deps, xvfb for headless mrview screenshots \
per KUL_NIS's README, p7zip for -B encrypted backups, dcmtk for send_2_orthanc)"
}

# ── 2. Docker ──────────────────────────────────────────────────────────────────

section_docker() {
    if have docker; then
        ok "docker already installed ($(docker --version 2>/dev/null || echo 'version unknown'))"
    else
        log "Installing Docker (docker.io from the Ubuntu/Mint apt repo)"
        run "sudo apt-get install -y docker.io"
        run "sudo usermod -aG docker \$(whoami)"
        warn "Added $(whoami) to the docker group — log out/in (or 'newgrp docker') \
before docker commands work without sudo."
    fi
    if have nvidia-smi; then
        ok "nvidia-smi found — for GPU-enabled containers (fastsurfer_gpu, HD-BET, \
hd-glio-auto) also install the NVIDIA Container Toolkit: \
https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/latest/install-guide.html \
(hardware/driver-specific, not automated by this script)"
    else
        warn "No nvidia-smi found — GPU-accelerated tools (FastSurfer, HD-BET, hd-glio-auto) \
will fall back to CPU, which is much slower but still works."
    fi
}

# ── 3. Miniforge ───────────────────────────────────────────────────────────────

section_miniforge() {
    if [ -x "$SOFTWARE_ROOT/miniforge3/bin/conda" ]; then
        ok "miniforge already installed at $SOFTWARE_ROOT/miniforge3"
        return
    fi
    log "Installing miniforge3"
    local installer="$SOFTWARE_ROOT/tmp/Miniforge3-Linux-x86_64.sh"
    run "curl -fsSL -o '$installer' https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
    run "bash '$installer' -b -p '$SOFTWARE_ROOT/miniforge3'"
    if [ "$DRY_RUN" -eq 1 ]; then
        echo "    [dry-run] would write $SOFTWARE_ROOT/.condarc"
    else
        cat > "$SOFTWARE_ROOT/.condarc" <<'EOF'
channels:
  - conda-forge
mirrored_channels:
  conda-forge:
    - https://conda.anaconda.org/conda-forge
    - https://prefix.dev/conda-forge
EOF
    fi
    run "'$SOFTWARE_ROOT/miniforge3/bin/conda' config --file '$SOFTWARE_ROOT/.condarc' --set always_yes true 2>/dev/null || true"
    ok "miniforge installed, .condarc written"
}

conda_bin() { echo "$SOFTWARE_ROOT/miniforge3/bin/conda"; }
mamba_bin() { echo "$SOFTWARE_ROOT/miniforge3/bin/mamba"; }
# Full path to a binary inside one of *our* envs. Deliberately NOT using
# `conda run -n <name>`/`mamba run -n <name>` for this: confirmed in practice
# that name-based env resolution can silently resolve to a completely different
# environment of the same name under an unrelated conda/mamba install elsewhere
# on the machine (e.g. a shared install's own "scilpy" env), even when the
# conda/mamba binary itself is invoked by full path under $SOFTWARE_ROOT — some
# inherited environment variable (CONDA_PREFIX, MAMBA_ROOT_PREFIX, etc.) takes
# precedence over the invoked binary's own location. Calling the target env's
# binary directly by its own full, unambiguous path sidesteps that entirely.
env_bin() { echo "$SOFTWARE_ROOT/miniforge3/envs/$1/bin/$2"; }

env_exists() {
    "$(conda_bin)" env list 2>/dev/null | awk '{print $1}' | grep -qx "$1"
}

# ── 4a. scilpy env (required — KUL_FWT's filtering/RecoBundles step) ──────────

section_env_scilpy() {
    if env_exists scilpy; then
        ok "conda env 'scilpy' already exists"
        return
    fi
    log "Creating 'scilpy' env (python 3.12, scilpy $SCILPY_COMMIT — pinned per \
KUL_NIS/README.md: v2.3.0)"
    run "$(mamba_bin) create -n scilpy python=3.12 -y"
    run "'$(env_bin scilpy pip)' install scilpy==2.3.0"
    run "'$(env_bin scilpy pip)' uninstall scilpy -y"
    local dest="$SOFTWARE_ROOT/src/scilpy"
    [ -d "$dest" ] || run "git clone https://github.com/scilus/scilpy.git '$dest'"
    run "cd '$dest' && git checkout $SCILPY_COMMIT"
    run "cd '$dest' && '$(env_bin scilpy pip)' install -e . --no-deps"
    ok "scilpy env ready. Point KUL_FWT/KUL_clinical_fmridti.sh at it with -f scilpy"
}

# ── 4b. HD-BET env (brain extraction) ─────────────────────────────────────────

section_env_hdbet() {
    if env_exists hd-bet-env; then
        ok "conda env 'hd-bet-env' already exists"
        return
    fi
    log "Creating 'hd-bet-env' (python 3.10, torch 2.6, HD-BET $HDBET_COMMIT)"
    run "$(mamba_bin) create -n hd-bet-env python=3.10 -y"
    local dest="$SOFTWARE_ROOT/src/HD-BET"
    [ -d "$dest" ] || run "git clone https://github.com/MIC-DKFZ/HD-BET.git '$dest'"
    run "cd '$dest' && git checkout $HDBET_COMMIT"
    # Install HD-BET and torch together in one resolver pass, with the CUDA (or
    # CPU) wheel index as primary and PyPI as extra. Doing this as two separate
    # pip calls (torch first, then 'pip install -e .' for HD-BET/nnunetv2) is a
    # confirmed bug: the second call's unconstrained transitive "torch"
    # requirement gets re-resolved against plain PyPI (no index-url given) and
    # silently replaces the working CUDA build with a newer, plain-PyPI one
    # that doesn't match this machine's CUDA driver -- e.g. a driver capped at
    # CUDA 12.2 ends up with a cu130-tagged torch that can't see the GPU at all.
    if have nvidia-smi; then
        # Pinned explicitly (not just "torch") so pip's resolver can't silently
        # drift to a newer, plain-PyPI (non-cu121) build that the -e install
        # would otherwise consider "already satisfies HD_BET's unconstrained
        # torch requirement" and skip re-resolving. 2.5.1+cu121 is the latest
        # version actually published on the cu121 wheel index as of this
        # writing; bump it if a newer +cu121 build exists and your driver
        # supports it (check with: pip index versions torch --index-url
        # https://download.pytorch.org/whl/cu121).
        run "'$(env_bin hd-bet-env pip)' install --index-url https://download.pytorch.org/whl/cu121 --extra-index-url https://pypi.org/simple 'torch==2.5.1+cu121' -e '$dest'"
    else
        warn "No nvidia-smi — installing CPU-only torch for hd-bet-env (much slower, but avoids a multi-GB CUDA wheel download you can't use)"
        run "'$(env_bin hd-bet-env pip)' install --index-url https://download.pytorch.org/whl/cpu --extra-index-url https://pypi.org/simple -e '$dest'"
    fi
    ok "hd-bet-env ready (used by KUL_dwiprep -m 1, KUL_anat_register)"
}

# ── 4c. resseg env (resection-cavity segmentation) ────────────────────────────

section_env_resseg() {
    if env_exists resseg; then
        ok "conda env 'resseg' already exists"
        return
    fi
    log "Creating 'resseg' env (python 3.8, resseg 0.3.7, antspyx 0.4.2)"
    run "$(mamba_bin) create -n resseg python=3.8 -y"
    # Pin simpleitk to the last version with a cp38 manylinux wheel before installing
    # resseg/antspyx -- torchio (a resseg dependency) allows any simpleitk!=2.0.*,!=2.1.1.1
    # and pip's resolver otherwise picks the latest (2.5.5+), which has no cp38 wheel and
    # fails building from source (conda toolchain doesn't see system crypt.h).
    run "'$(env_bin resseg pip)' install simpleitk==2.4.1"
    run "'$(env_bin resseg pip)' install resseg==0.3.7 antspyx==0.4.2"
    ok "resseg env ready (used by KUL_anat_segment_tumor.sh)"
}

# ── 4d. Karawun env (Brainlab export) ─────────────────────────────────────────

section_env_karawun() {
    if [ "$USE_KARAWUN_DEV" -eq 1 ]; then
        if env_exists KarawunDev; then
            ok "conda env 'KarawunDev' already exists"
            return
        fi
        log "Creating 'KarawunDev' (editable install, karawun $KARAWUN_COMMIT + dcm2bids 3.1.1)"
        run "$(mamba_bin) create -n KarawunDev python=3.10 -y"
        run "'$(env_bin KarawunDev pip)' install dcm2bids>=3.1"
        local dest="$SOFTWARE_ROOT/src/karawun"
        [ -d "$dest" ] || run "git clone https://github.com/treanus/karawun.git '$dest'"
        run "cd '$dest' && git checkout $KARAWUN_COMMIT"
        run "cd '$dest' && '$(env_bin KarawunDev pip)' install -e . --no-deps"
        ok "KarawunDev ready"
    else
        if env_exists KarawunEnv; then
            ok "conda env 'KarawunEnv' already exists"
            return
        fi
        log "Creating 'KarawunEnv' (plain conda-forge package, no git clone needed)"
        run "$(mamba_bin) create -n KarawunEnv python=3.8 karawun=0.2.5.4 -c conda-forge -y"
        ok "KarawunEnv ready (used by KUL_karawun_prepare.sh / KUL_karawun2brainlab.sh)"
    fi
}

# ── 4e. FastSurfer env (faster recon-all alternative; needs FreeSurfer 8 for compat) ──

section_env_fastsurfer() {
    local dest="$SOFTWARE_ROOT/src/FastSurfer"
    # FastSurfer's native-install method upstream is uv + a per-checkout .venv,
    # not a conda env (confirmed against the current repo: env/environment_gpu.yml
    # no longer exists, only an empty env/fastsurfer.yml placeholder; their own
    # docs now say 'uv venv' + 'uv pip sync requirements.txt --torch-backend
    # <backend>'). This section installs uv into the isolated miniforge3 base env
    # (self-contained — just a pip package, no separate download/PATH entry) and
    # uses it to build FastSurfer's own venv the same way.
    if [ -d "$dest/.venv" ]; then
        ok "FastSurfer already cloned with a .venv at $dest"
        return
    fi
    log "Cloning FastSurfer and building its uv venv"
    [ -d "$dest" ] || run "git clone https://github.com/Deep-MI/FastSurfer.git '$dest'"
    if [ ! -x "$SOFTWARE_ROOT/miniforge3/bin/uv" ]; then
        run "'$SOFTWARE_ROOT/miniforge3/bin/pip' install uv"
    fi
    local uv="$SOFTWARE_ROOT/miniforge3/bin/uv"
    run "cd '$dest' && '$uv' venv --python 3.12"
    # torch-backend: pick the newest CUDA backend that (a) actually has a build
    # for whatever version requirements.txt currently pins and (b) your driver
    # supports (check with: nvidia-smi, then pip index versions torch
    # --index-url https://download.pytorch.org/whl/<backend>). cu118 is the
    # verified-working choice as of this writing (driver capped at CUDA 12.2;
    # requirements.txt pins torch==2.7.1, which has no cu121 build but does
    # have cu118). Re-check both if this section starts failing later.
    local torch_backend="cu118"
    have nvidia-smi || torch_backend="cpu"
    # Must be compile-then-sync, NOT a direct 'uv pip sync requirements.txt' —
    # confirmed bug: requirements.txt deliberately excludes nvidia-* CUDA
    # runtime packages from its pinned list (per its own header comment) and
    # expects the compile step to resolve them in transitively; 'sync' alone
    # only installs packages literally listed in the file, silently skipping
    # the runtime libs and leaving torch unable to import
    # ("libcudart.so.11.0: cannot open shared object file").
    run "cd '$dest' && '$uv' pip compile --torch-backend $torch_backend requirements.txt --python .venv/bin/python | '$uv' pip sync --torch-backend $torch_backend --python .venv/bin/python -"
    # Wrapper shim: run_fastsurfer.sh calls plain python3 with no self-activation
    # of its own venv, and it's invoked both by full path and by bare name (via
    # PATH) across different KUL_NIS/KUL_VBG scripts. Renaming the real script
    # and dropping a wrapper at the original filename covers both invocation
    # styles without touching any KUL_NIS/KUL_VBG source.
    if [ -f "$dest/run_fastsurfer.sh" ] && [ ! -f "$dest/run_fastsurfer.real.sh" ]; then
        run "mv '$dest/run_fastsurfer.sh' '$dest/run_fastsurfer.real.sh'"
    fi
    if [ "$DRY_RUN" -eq 1 ]; then
        echo "    [dry-run] would write the run_fastsurfer.sh venv-activation wrapper shim"
    else
        cat > "$dest/run_fastsurfer.sh" <<'WRAPPER_EOF'
#!/usr/bin/env bash
# Wrapper shim: activates FastSurfer's uv-managed .venv (run_fastsurfer.real.sh
# calls plain python3 with no self-activation), then hands off. Scoped to this
# process via exec -- never leaks back to the calling shell.
set -euo pipefail
if [[ -z "${BASH_SOURCE[0]}" ]]; then THIS_SCRIPT="$0"
else THIS_SCRIPT="${BASH_SOURCE[0]}"
fi
FASTSURFER_HOME="$(cd "$(dirname "$THIS_SCRIPT")" &>/dev/null && pwd)"
source "$FASTSURFER_HOME/.venv/bin/activate"
exec "$FASTSURFER_HOME/run_fastsurfer.real.sh" "$@"
WRAPPER_EOF
        chmod +x "$dest/run_fastsurfer.sh"
    fi
    warn "Set FASTSURFER_HOME=$dest in your shell — KUL_anat_segment_tumor.sh calls \
\$FASTSURFER_HOME/run_fastsurfer.sh directly (and KUL_FS_multiparc.sh/KUL_clinical_fmridti.sh/KUL_VBG.sh \
call it by bare name via PATH); the wrapper shim handles venv activation either way. \
KUL_NIS's README requires FastSurfer v3+ for FreeSurfer 8 compatibility."
    ok "FastSurfer section done (uv venv + wrapper shim — verify torch imports and reports cuda available if you have a GPU)"
}

# ── 4f. rsfMRI pipeline env (KUL_NIS/share/rsfmri_pipeline) ───────────────────

section_env_rsfmri() {
    if env_exists rsfmri_env; then
        ok "conda env 'rsfmri_env' already exists"
        return
    fi
    log "Creating 'rsfmri_env' (nilearn/nibabel/numpy/scipy/pandas/matplotlib/pyyaml — \
KUL_run_rsfMRI_networks.sh's -c flag names this env; no hardcoded version pins upstream, \
so this creates everything fresh together in one resolver pass to avoid ABI drift \
between packages installed piecemeal over time)"
    run "$(mamba_bin) create -n rsfmri_env -c conda-forge python=3.10 nilearn nibabel numpy scipy pandas matplotlib pyyaml -y"
    ok "rsfmri_env ready. Pass it to KUL_run_rsfMRI_networks.sh with -c rsfmri_env"
}

# ── 5. Sibling repos, pinned branches ──────────────────────────────────────────

section_repos() {
    log "Cloning KUL_NIS / KUL_VBG / KUL_FWT at their pinned branches"
    local src="$SOFTWARE_ROOT/src"

    clone_at() {
        local url="$1" dir="$2" branch="$3"
        if [ -d "$src/$dir" ]; then
            ok "$dir already cloned at $src/$dir"
        else
            run "git clone '$url' '$src/$dir'"
            run "cd '$src/$dir' && git checkout '$branch'"
            ok "$dir @ $branch"
        fi
    }

    clone_at https://github.com/treanus/KUL_NIS.git       KUL_NIS    "$KUL_NIS_BRANCH"
    clone_at https://github.com/KUL-Radneuron/KUL_VBG.git KUL_VBG    "$KUL_VBG_BRANCH"
    clone_at https://github.com/KUL-Radneuron/KUL_FWT.git KUL_FWT    "$KUL_FWT_BRANCH"

    warn "Add these to PATH (see the ~/.bashrc block this script writes): \
$src/KUL_NIS, $src/KUL_VBG, $src/KUL_FWT, and KUL_NIS/KUL_DTI_ALPS. \
KUL_FWT's own scripts must resolve on PATH before running KUL_clinical_fmridti.sh — \
it no longer auto-prepends a sibling ../KUL_FWT folder."
}

# ── 5b. LoRE-SD (DWI decomposition; regular pip package, NOT an MRtrix3 module build) ──
# Despite living on a branch named 'mrtrix_module' and mentioning MRtrix3 module
# builds in older docs, LoRE-SD is just installed as a normal editable pip
# package (setup.py + pybind11/nlopt/cmake) — confirmed working via a validated
# install sequence. No mrtrix3 source tree needed.

section_env_lore_sd() {
    if env_exists lore_sd; then
        ok "conda env 'lore_sd' already exists"
        return
    fi
    local dest="$SOFTWARE_ROOT/src/LoRE-SD"
    [ -d "$dest" ] || run "git clone --branch mrtrix_module https://github.com/SiebeLeysen/LoRE-SD.git '$dest'"
    log "Creating 'lore_sd' env (python 3.10, nlopt/cmake/compilers via conda-forge, editable pip install)"
    run "$(mamba_bin) create -n lore_sd python=3.10 -y"
    run "$(mamba_bin) install -n lore_sd -c conda-forge nlopt cmake compilers -y"
    run "'$(env_bin lore_sd pip)' install pybind11 numpy setuptools wheel"
    run "cd '$dest' && '$(env_bin lore_sd pip)' install -e . --no-build-isolation"
    ok "lore_sd env ready — lore_dwi2decomposition / lore_decomposition2contrast available in it \
(used by KUL_dwiprep -D run_dwiprep_lore_sd.txt)"
}

# ── 6. MRtrix3 (build from source) ────────────────────────────────────────────

section_mrtrix3() {
    if have mrconvert; then
        ok "mrconvert already on PATH ($(mrconvert -version 2>&1 | head -1))"
        return
    fi
    local dest="$SOFTWARE_ROOT/src/mrtrix3"
    log "Building MRtrix3 @ $MRTRIX3_COMMIT (this takes a while)"
    [ -d "$dest" ] || run "git clone https://github.com/MRtrix3/mrtrix3.git '$dest'"
    run "cd '$dest' && git checkout '$MRTRIX3_COMMIT'"
    run "cd '$dest' && ./configure"
    run "cd '$dest' && NUMBER_OF_PROCESSORS=$NCPU ./build"
    warn "If this commit fails to build, SOFTWARE_ROOT_SETUP.md flags an older \
known-good fallback: commit 5a3a8bf6 (3.0.4-537-g5a3a8bf6). Re-run with \
MRTRIX3_COMMIT=5a3a8bf6 ./setup_environment.sh --only mrtrix3 if so."
    ok "mrtrix3 built at $dest (add $dest/bin to PATH — done automatically in the bashrc section)"
}

# ── 7. ANTs (build from source) ───────────────────────────────────────────────

section_ants() {
    if have antsRegistrationSyN.sh; then
        ok "ANTs already on PATH"
        return
    fi
    local dest="$SOFTWARE_ROOT/src/ANTs"
    log "Building ANTs @ $ANTS_COMMIT (v2.4.4.post20) — this takes a long while"
    [ -d "$dest" ] || run "git clone https://github.com/ANTsX/ANTs.git '$dest'"
    run "cd '$dest' && git checkout '$ANTS_COMMIT'"
    run "mkdir -p '$dest/build'"
    run "cd '$dest/build' && cmake -DCMAKE_INSTALL_PREFIX='$SOFTWARE_ROOT/ANTs_install' .."
    run "cd '$dest/build' && make -j$NCPU"
    run "cd '$dest/build/ANTS-build' && make install"
    ok "ANTs installed to $SOFTWARE_ROOT/ANTs_install"
}

# ── 8. FSL ─────────────────────────────────────────────────────────────────────

section_fsl() {
    if [ -x "$SOFTWARE_ROOT/FSL/bin/flirt" ] || have flirt; then
        ok "FSL already installed"
        return
    fi
    log "Installing FSL $FSL_VERSION via the official installer"
    local installer="$SOFTWARE_ROOT/tmp/fslinstaller.py"
    run "curl -fsSL -o '$installer' https://fsl.fmrib.ox.ac.uk/fsldownloads/fslinstaller.py"
    warn "FSL installer syntax changes across releases (SOFTWARE_ROOT_SETUP.md flags this) — \
if the -V flag below is rejected, check 'python3 $installer --help' for the current syntax."
    run "python3 '$installer' -d '$SOFTWARE_ROOT/FSL' -V '$FSL_VERSION'"
    ok "FSL installed to $SOFTWARE_ROOT/FSL"
}

# ── 9. FreeSurfer (license-gated — can't be auto-downloaded) ──────────────────

section_freesurfer_check() {
    log "FreeSurfer 8.2.0 can't be fetched by this script — it requires a personal \
license registered at https://surfer.nmr.mgh.harvard.edu/registration.html"
    if [ -f "$SOFTWARE_ROOT/freesurfer/SetUpFreeSurfer.sh" ]; then
        ok "FreeSurfer already present at $SOFTWARE_ROOT/freesurfer"
    else
        warn "Nothing found at $SOFTWARE_ROOT/freesurfer. To finish this step:"
        echo "      1. Register and download FreeSurfer 8.2.0 (or your own already-licensed copy)"
        echo "         from https://surfer.nmr.mgh.harvard.edu/fswiki/rel8download"
        echo "      2. Extract it to: $SOFTWARE_ROOT/freesurfer"
        echo "      3. Place your license.txt at: $SOFTWARE_ROOT/freesurfer/license.txt"
        echo "      (If you already have a licensed copy elsewhere on this machine, just"
        echo "       'cp -a /path/to/existing/freesurfer $SOFTWARE_ROOT/freesurfer' instead.)"
    fi
    if [ -f "$SOFTWARE_ROOT/freesurfer/license.txt" ]; then
        ok "license.txt present"
    else
        warn "No license.txt found yet — KUL_VBG.sh, KUL_FS_multiparc.sh, MSBP (via Docker), \
and recon-all itself will all refuse to run without one."
    fi
}

# ── 10. Docker images used by the pipeline ────────────────────────────────────

section_docker_images() {
    if ! have docker; then
        warn "docker not found — skipping image pulls (run the docker section first)"
        return
    fi
    log "Pulling Docker images used by the clinical pipeline"
    run "docker pull nipreps/fmriprep:$FMRIPREP_VERSION"
    run "docker pull poldracklab/mriqc:latest"
    run "docker pull leonyichencai/synb0-disco:v3.0"
    run "docker pull sebastientourbier/multiscalebrainparcellator:v1.1.1"
    if [ "$INCLUDE_HDGLIOAUTO" -eq 1 ]; then
        warn "hd-glio-auto has no single canonical image name pinned in this codebase — \
check github.com/NeuroAI-HD/HD-GLIO-AUTO's current README for the image/tag to pull, \
or use its local-install path (python run.py) referenced in KUL_anat_segment_tumor.sh instead."
    fi
    ok "Docker images pulled. FreeSurfer license.txt must be bind-mounted at \
/opt/freesurfer/license.txt (fmriprep) or /usr/local/freesurfer/license.txt (MSBP) at run time — \
see each script's existing 'docker run' invocation for the exact mount."
}

# ── 11. ~/.bashrc environment block ───────────────────────────────────────────

section_bashrc() {
    local marker="# >>> KUL_NIS/KUL_VBG/KUL_FWT environment (managed block) >>>"
    local marker_end="# <<< KUL_NIS/KUL_VBG/KUL_FWT environment (managed block) <<<"
    if grep -qF "$marker" "$HOME/.bashrc" 2>/dev/null; then
        ok "~/.bashrc already has the managed environment block — not touching it. \
Edit it by hand, or delete the block between the marker lines and re-run this section."
        return
    fi
    log "Appending environment block to ~/.bashrc"
    if [ "$DRY_RUN" -eq 1 ]; then
        echo "    [dry-run] would append the managed environment block to $HOME/.bashrc"
        return
    fi
    local src="$SOFTWARE_ROOT/src"
    cat >> "$HOME/.bashrc" <<EOF

$marker
export SOFTWARE_ROOT="$SOFTWARE_ROOT"

# --- Conda / miniforge ---
export PATH="\$SOFTWARE_ROOT/miniforge3/bin:\$PATH"
source "\$SOFTWARE_ROOT/miniforge3/etc/profile.d/conda.sh"
# conda.sh alone only wires up 'conda activate/deactivate'. mamba has its own
# shell integration ('mamba shell hook') -- without this, 'mamba activate' fails
# with "mamba is running as a subprocess and can't modify the parent shell."
# Sourcing both means either 'conda activate <env>' or 'mamba activate <env>'
# works interchangeably (mamba is faster for create/install; activation itself
# is the same underlying mechanism either way).
if command -v mamba >/dev/null 2>&1; then
    eval "\$(mamba shell hook --shell bash)"
fi
export CONDA_PKGS_DIRS="\$SOFTWARE_ROOT/cache/conda_pkgs"

# --- Caches / tmp ---
export XDG_CACHE_HOME="\$SOFTWARE_ROOT/cache/xdg"
export PIP_CACHE_DIR="\$SOFTWARE_ROOT/cache/pip"
export TMPDIR="\$SOFTWARE_ROOT/tmp"
export APPTAINER_TMPDIR="\$SOFTWARE_ROOT/tmp/apptainer"

# --- FSL ---
export FSLDIR="\$SOFTWARE_ROOT/FSL"
[ -f "\$FSLDIR/etc/fslconf/fsl.sh" ] && source "\$FSLDIR/etc/fslconf/fsl.sh"
export PATH="\$FSLDIR/share/fsl/bin:\$PATH"

# --- FreeSurfer ---
export FREESURFER_HOME="\$SOFTWARE_ROOT/freesurfer"
[ -f "\$FREESURFER_HOME/SetUpFreeSurfer.sh" ] && source "\$FREESURFER_HOME/SetUpFreeSurfer.sh"
export SUBJECTS_DIR="\$FREESURFER_HOME/subjects"
export FS_LICENSE="\$FREESURFER_HOME/license.txt"

# --- FastSurfer ---
export FASTSURFER_HOME="\$SOFTWARE_ROOT/src/FastSurfer"

# --- ANTs ---
export ANTSPATH="\$SOFTWARE_ROOT/ANTs_install/bin/"
export PATH="\$ANTSPATH:\$PATH"

# --- mrtrix3 ---
export PATH="\$SOFTWARE_ROOT/src/mrtrix3/bin:\$PATH"
export PYTHONPATH="\$SOFTWARE_ROOT/src/mrtrix3/lib:\$PYTHONPATH"

# --- KUL_NIS / KUL_VBG / KUL_FWT / KUL_DTI_ALPS on PATH ---
export PATH="$src/KUL_NIS:$src/KUL_VBG:$src/KUL_FWT:$src/KUL_NIS/KUL_DTI_ALPS:\$PATH"
$marker_end
EOF
    ok "~/.bashrc updated — run 'source ~/.bashrc' or open a new shell to pick it up"
}

# ── 12. Verify — read-only health check + summary table ──────────────────────
# Doesn't install or modify anything. Safe to run any time, including on its
# own (--only verify) to sanity-check a machine you didn't set up with this
# script, or to see what's left after a partial/interrupted run.

VERIFY_OK=0
VERIFY_WARN=0
VERIFY_FAIL=0

vprint() { printf "  [%-4s] %-24s %s\n" "$1" "$2" "$3"; }
vok()    { vprint "OK"   "$1" "$2"; VERIFY_OK=$((VERIFY_OK + 1)); }
vwarn()  { vprint "WARN" "$1" "$2"; VERIFY_WARN=$((VERIFY_WARN + 1)); }
vfail()  { vprint "FAIL" "$1" "$2"; VERIFY_FAIL=$((VERIFY_FAIL + 1)); }

section_verify() {
    log "Verifying the environment (read-only — nothing here modifies anything)"
    VERIFY_OK=0; VERIFY_WARN=0; VERIFY_FAIL=0

    echo "  -- core imaging tools (checked on PATH, wherever they actually come from) --"
    if have mrconvert; then
        vok "MRtrix3" "$(mrconvert -version 2>&1 | head -1 | sed -e 's/^== *//' -e 's/ *==$//')"
    else
        vfail "MRtrix3" "mrconvert not found on PATH"
    fi

    if have antsRegistrationSyN.sh; then
        vok "ANTs" "antsRegistrationSyN.sh resolves"
    else
        vfail "ANTs" "antsRegistrationSyN.sh not found on PATH"
    fi

    if have flirt; then
        # 'flirt -version' only reports FLIRT's own major.minor (e.g. "6.0"), not
        # FSL's actual patch version -- real behavioural differences exist between
        # e.g. 6.0.5 and 6.0.6 downstream in this pipeline, so pull the precise
        # version FSL itself records instead.
        _fsl_dir="${FSLDIR:-$SOFTWARE_ROOT/FSL}"
        if [ ! -f "$_fsl_dir/etc/fslversion" ]; then
            _flirt_bin=$(command -v flirt)
            _fsl_dir="${_flirt_bin%/bin/flirt}"
            _fsl_dir="${_fsl_dir%/share/fsl}"
        fi
        if [ -f "$_fsl_dir/etc/fslversion" ]; then
            vok "FSL" "$(cat "$_fsl_dir/etc/fslversion") (exact, from etc/fslversion)"
        else
            vwarn "FSL" "$(flirt -version 2>&1 | head -1) — exact patch version not found; only major.minor available"
        fi
    else
        vfail "FSL" "flirt not found on PATH"
    fi

    if have recon-all; then
        local _fs_lic="${FS_LICENSE:-${FREESURFER_HOME:-$SOFTWARE_ROOT/freesurfer}/license.txt}"
        if [ -f "$_fs_lic" ]; then
            vok "FreeSurfer" "recon-all resolves, license.txt present"
        else
            vwarn "FreeSurfer" "recon-all resolves, but no license.txt found at $_fs_lic"
        fi
    else
        vfail "FreeSurfer" "recon-all not found on PATH (license-gated install — see freesurfer-check section)"
    fi

    echo
    echo "  -- Docker --"
    if have docker; then
        if docker info >/dev/null 2>&1; then
            vok "Docker" "installed, usable without sudo"
        else
            vwarn "Docker" "installed, but not usable yet (log out/in after group add?)"
        fi
    else
        vfail "Docker" "not installed"
    fi

    if have docker; then
        for img in "nipreps/fmriprep:$FMRIPREP_VERSION" "poldracklab/mriqc:latest" \
                   "leonyichencai/synb0-disco:v3.0" "sebastientourbier/multiscalebrainparcellator:v1.1.1"; do
            if docker image inspect "$img" >/dev/null 2>&1; then
                vok "docker image" "$img"
            else
                vwarn "docker image" "$img — not pulled"
            fi
        done
    fi

    echo
    echo "  -- conda environments --"
    if [ -x "$(conda_bin)" ]; then
        vok "miniforge" "$(conda_bin)"

        _verify_env_bin() {
            # $1=env name, $2=label, $3=binary inside envs/$1/bin, $4... = args to smoke-test it
            local env="$1" label="$2" bin="$3"; shift 3
            local path="$SOFTWARE_ROOT/miniforge3/envs/$env/bin/$bin"
            if [ ! -d "$SOFTWARE_ROOT/miniforge3/envs/$env" ]; then
                vfail "conda: $env" "environment does not exist"
            elif [ ! -x "$path" ]; then
                vwarn "conda: $env" "environment exists, but $bin not found in it"
            elif "$path" "$@" >/dev/null 2>&1; then
                vok "conda: $env" "$label"
            else
                vwarn "conda: $env" "environment exists, but '$bin $*' failed"
            fi
        }

        _verify_env_bin scilpy "scil_tractogram_filter_by_roi resolves" scil_tractogram_filter_by_roi --help
        _verify_env_bin hd-bet-env "hd-bet resolves" hd-bet --help
        _verify_env_bin resseg "resseg importable" python -c "import resseg, ants"
        if [ "$USE_KARAWUN_DEV" -eq 1 ]; then
            _verify_env_bin KarawunDev "karawun importable" python -c "import karawun"
        else
            _verify_env_bin KarawunEnv "karawun importable" python -c "import karawun"
        fi
        _verify_env_bin rsfmri_env "nilearn importable" python -c "import nilearn, nibabel, pandas, matplotlib, yaml"
        _verify_env_bin lore_sd "lore_dwi2decomposition resolves" lore_dwi2decomposition --help
    else
        vfail "miniforge" "not installed at $SOFTWARE_ROOT/miniforge3"
    fi

    echo
    echo "  -- FastSurfer (uv venv, not a conda env) --"
    local fsdest="$SOFTWARE_ROOT/src/FastSurfer"
    if [ ! -x "$fsdest/.venv/bin/python" ]; then
        vfail "FastSurfer .venv" "not found at $fsdest/.venv"
    elif "$fsdest/.venv/bin/python" -c "import torch" >/dev/null 2>&1; then
        vok "FastSurfer .venv" "torch importable"
    else
        vwarn "FastSurfer .venv" "exists, but torch import failed"
    fi
    if [ -x "$fsdest/run_fastsurfer.sh" ] && [ -f "$fsdest/run_fastsurfer.real.sh" ]; then
        vok "FastSurfer wrapper" "run_fastsurfer.sh shim + run_fastsurfer.real.sh both present"
    else
        vwarn "FastSurfer wrapper" "run_fastsurfer.sh shim or run_fastsurfer.real.sh missing at $fsdest"
    fi

    echo
    echo "  -- sibling repos --"
    _verify_repo() {
        local name="$1" expect_branch="$2"
        local dir="$SOFTWARE_ROOT/src/$name"
        if [ ! -d "$dir/.git" ]; then
            vfail "repo: $name" "not cloned"
            return
        fi
        local actual_branch
        actual_branch=$(git -C "$dir" branch --show-current 2>/dev/null)
        if [ "$actual_branch" = "$expect_branch" ]; then
            vok "repo: $name" "@ $actual_branch"
        else
            vwarn "repo: $name" "on '$actual_branch', expected '$expect_branch' (detached HEAD is normal for a pinned-commit checkout — verify the commit by hand)"
        fi
    }
    _verify_repo KUL_NIS "$KUL_NIS_BRANCH"
    _verify_repo KUL_VBG "$KUL_VBG_BRANCH"
    _verify_repo KUL_FWT "$KUL_FWT_BRANCH"
    if [ -d "$SOFTWARE_ROOT/src/LoRE-SD/.git" ]; then
        vok "repo: LoRE-SD" "cloned @ $(git -C "$SOFTWARE_ROOT/src/LoRE-SD" branch --show-current 2>/dev/null) (see 'conda: lore_sd' above for the actual install check)"
    else
        vfail "repo: LoRE-SD" "not cloned"
    fi

    echo
    echo "  -- shell integration --"
    if grep -qF "# >>> KUL_NIS/KUL_VBG/KUL_FWT environment (managed block) >>>" "$HOME/.bashrc" 2>/dev/null; then
        vok "~/.bashrc" "managed environment block present"
    else
        vwarn "~/.bashrc" "managed block not found — run the bashrc section (or add it by hand if you deliberately skipped it, e.g. on a machine with pre-existing FSL/ANTs/FreeSurfer/mrtrix3)"
    fi

    echo
    echo "  ────────────────────────────────────────────────────────────"
    printf "  %s OK   %s WARN   %s FAIL\n" "$VERIFY_OK" "$VERIFY_WARN" "$VERIFY_FAIL"
    if [ "$VERIFY_FAIL" -gt 0 ]; then
        fail "One or more required components are missing — see the FAIL lines above."
    elif [ "$VERIFY_WARN" -gt 0 ]; then
        warn "Everything required is present; some secondary/optional pieces need attention — see the WARN lines above."
    else
        ok "Everything checked out clean."
    fi
}

# ── Main ───────────────────────────────────────────────────────────────────────

maybe_run_section() {
    local name="$1" flag_var="$2" func="$3"
    if ! section_enabled "$name"; then
        return
    fi
    if [ "${!flag_var}" -eq 1 ]; then
        echo
        "$func"
    else
        warn "Skipping '$name' (disabled by config — set $flag_var=1 to enable)"
    fi
}

maybe_run_section apt              DO_APT               section_apt
maybe_run_section docker           DO_DOCKER             section_docker
maybe_run_section miniforge        DO_MINIFORGE          section_miniforge
maybe_run_section env-scilpy       DO_ENV_SCILPY         section_env_scilpy
maybe_run_section env-hdbet        DO_ENV_HDBET          section_env_hdbet
maybe_run_section env-resseg       DO_ENV_RESSEG          section_env_resseg
maybe_run_section env-karawun      DO_ENV_KARAWUN         section_env_karawun
maybe_run_section env-fastsurfer   DO_ENV_FASTSURFER      section_env_fastsurfer
maybe_run_section env-rsfmri       DO_ENV_RSFMRI          section_env_rsfmri
maybe_run_section env-lore-sd      DO_ENV_LORE_SD         section_env_lore_sd
maybe_run_section repos            DO_REPOS               section_repos
maybe_run_section mrtrix3          DO_MRTRIX3             section_mrtrix3
maybe_run_section ants             DO_ANTS                section_ants
maybe_run_section fsl              DO_FSL                 section_fsl
maybe_run_section freesurfer-check DO_FREESURFER_CHECK    section_freesurfer_check
maybe_run_section docker-images    DO_DOCKER_IMAGES       section_docker_images
maybe_run_section bashrc           DO_BASHRC              section_bashrc
maybe_run_section verify           DO_VERIFY              section_verify

echo
log "Done. Manual steps this script cannot do for you:"
echo "   - FreeSurfer license.txt (registration-gated, see the freesurfer-check section output above)"
echo "   - SPM12 + MATLAB, for KUL_fmriproc_spm's task-fMRI GLM (commercial software; nilearn"
echo "     via KUL_fmriproc_nilearn_new.sh is the license-free alternative already in the pipeline)"
echo "   - NVIDIA driver + Container Toolkit, if you want GPU acceleration for FastSurfer/HD-BET/hd-glio-auto"
echo "   - hd-glio-auto itself, if INCLUDE_HDGLIOAUTO=1 wasn't enough — check its own README"
echo "   - 'source ~/.bashrc' (or open a new terminal) to pick up the new PATH/env vars"
