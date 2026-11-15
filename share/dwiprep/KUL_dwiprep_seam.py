#!/usr/bin/env python3
"""KUL_dwiprep_seam.py -- screening gate + correction of the MB=2 band-seam dark slice.

Called by KUL_dwiprep.sh (option -z) on dwi/degibbs.mif, i.e. AFTER dwidenoise + mrdegibbs and BEFORE
dwifslpreproc / shard-recon (topup/eddy resample the data, after which the seam no longer sits on one slice).

Why this order: the dark slice is a smooth multiplicative gain on raw slice N//2 that does NOT scale the noise
(measured on two MB=2 test datasets: signal 0.74-0.79 of the neighbours, noise 1.01). Dividing by the gain before dwidenoise would feed MP-PCA a slice
with ~1.3x the noise of its neighbours; after dwidenoise the noise is mostly gone and only the small residual is rescaled.

What it does, per phase-encoding group (volumes sharing the same pe-table row; a reverse-PE series is its own group):
  1. Screen: dip = median parenchyma signal of the seam slice / mean of its two neighbours, on the mean of ALL volumes of
     the group (b0 + DWI; robust against a single noisy or corrupted b0). The seam slice is searched in {N//2-1, N//2}
     (the array may be flipped relative to the scanner slice order; the dark slice is the one with the lower ratio).
       dip <  GATE_CORRECT (0.95)         -> correct
       GATE_CORRECT <= dip < GATE_FLAG (0.98) -> correct + flag for review
       dip >= GATE_FLAG                    -> pass through (group untouched)
  2. Correct: smooth in-plane gain (normalised-convolution, sigma 6 voxels) estimated on the group's mean DWI parenchyma
     (CSF blended in with K_CSF = 1.0), divided out of the seam slice of every volume in the group.
  The whole step is skipped unless MultibandAccelerationFactor == 2 in the header (the seam position N/2 is a property of MB=2).
  b0 slice alternation (seen in a single-b0 acquisition) is only DETECTED and reported (b0_alt > 0.08); it is not corrected.

usage: KUL_dwiprep_seam.py <in.mif> <out.mif> <qc_dir> [--force-mb2] [--dry]
  out.mif is written only if at least one group was corrected.  stdout last line: SEAM_RESULT=corrected|pass|skipped
  qc_dir gets seam_gate.json and seam_gate.png.  Needs mrtrix3 on PATH, numpy, scipy, nibabel, matplotlib (optional).
Same estimator as the validated KUL_reanalysis code/seam_qc/seam_correct.py and seam_gate.py (46-subject screen 13 Nov 2026).
"""
import sys, os, json, shutil, subprocess, tempfile
import numpy as np
import nibabel as nib
from scipy import ndimage as ndi

GATE_CORRECT, GATE_FLAG, ALT_THR, SIGMA, K_CSF = 0.95, 0.98, 0.08, 6.0, 1.0


def run(cmd):
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def mrinfo(mif, *a):
    r = subprocess.run(['mrinfo', mif, *a], capture_output=True, text=True)
    return r.stdout.strip() if r.returncode == 0 else ''


def nc(x, w, s=SIGMA, wmin=0.05):
    num = ndi.gaussian_filter(x * w, s); den = ndi.gaussian_filter(w.astype(float), s)
    return np.where(den > wmin, num / np.maximum(den, 1e-9), np.nan)


def masks(b0, dw, have_dwi):
    thr = np.percentile(b0[b0 > 0], 99) * 0.15
    m = ndi.binary_fill_holes(ndi.binary_opening(b0 > thr, iterations=2))
    if not have_dwi:
        return m, np.zeros_like(b0), m, np.zeros_like(m)
    adc = np.where(m, -np.log(np.clip(dw, 1e-3, None) / np.clip(b0, 1e-3, None)) / 1000., 0)
    return m, adc, m & (adc > 0.4e-3) & (adc < 1.0e-3), m & (adc > 2.0e-3)


def dip_at(vol, sel, k):
    p = [np.median(vol[..., j][sel[..., j]]) if sel[..., j].sum() > 50 else np.nan for j in (k - 1, k, k + 1)]
    return float(p[1] / (0.5 * (p[0] + p[2])))


def dip_profile(vol, sel):
    nz = vol.shape[2]; out = np.full(nz, np.nan)
    for k in range(1, nz - 1):
        out[k] = dip_at(vol, sel, k)
    return out


def seam_gain(dwv, b0v, m, adc, k):
    mk = m[..., k] & m[..., k - 1] & m[..., k + 1]
    den_d = 0.5 * (dwv[..., k - 1] + dwv[..., k + 1])
    rd = np.where(mk & (den_d > 0), dwv[..., k] / np.maximum(den_d, 1e-6), 0)
    a = ndi.gaussian_filter(adc[..., k], 0.7)
    wc = np.clip((a - 1.0e-3) / 1.0e-3, 0, 1) * mk
    gp = nc(rd, mk * (1 - wc))
    sel = mk & (wc == 0)
    glob = np.median(rd[sel]) if sel.sum() > 100 else 1.0
    gp = np.where(np.isnan(gp), glob, gp)
    g = (1 - wc) * gp + wc * np.minimum(gp * K_CSF, 1.0)
    g = np.where(mk, g, gp)
    return np.clip(g, 0.5, 1.3).astype(np.float32)


def alternation(b0, dw, m):
    nz = b0.shape[2]
    zs = [z for z in range(nz) if m[..., z].sum() > 500]
    if len(zs) < 8:
        return float('nan')
    r = np.array([np.median(b0[..., z][m[..., z]] / np.maximum(dw[..., z][m[..., z]], 1)) for z in zs])
    a = r[1:-1] / (0.5 * (r[:-2] + r[2:])) - 1
    sgn = np.where(np.arange(len(a)) % 2 == 0, 1, -1)
    return float(np.mean(np.abs(ndi.gaussian_filter1d(a * sgn, 3))))


def main():
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    in_mif, out_mif, qcd = args[:3]
    force_mb2, dry = '--force-mb2' in sys.argv, '--dry' in sys.argv
    os.makedirs(qcd, exist_ok=True)
    qc = dict(input=in_mif, gate=dict(correct=GATE_CORRECT, flag=GATE_FLAG), groups=[], flags=[])

    def finish(status, why=None):
        qc['status'] = status
        if why: qc['why'] = why
        json.dump(qc, open(f'{qcd}/seam_gate.json', 'w'), indent=1)
        print(f'SEAM_RESULT={status}' + (f' ({why})' if why else ''))
        sys.exit(0)

    mb = mrinfo(in_mif, '-property', 'MultibandAccelerationFactor')
    qc['MultibandAccelerationFactor'] = mb or None
    if not force_mb2 and mb.split()[:1] != ['2']:
        finish('skipped', f'MultibandAccelerationFactor={mb or "absent"}, the seam correction applies to MB=2 data only')
    size = [int(x) for x in mrinfo(in_mif, '-size').split()]
    if len(size) != 4 or size[3] < 4:
        finish('skipped', f'not a 4D dwi series (size {size})')
    nz = size[2]; h = nz // 2

    with tempfile.TemporaryDirectory(prefix='seam_', dir=os.environ.get('MRTRIX_TMPFILE_DIR') or None) as tmp:
        run(['mrinfo', in_mif, '-export_grad_mrtrix', f'{tmp}/g.b'])
        gb = np.array([[float(v) for v in l.split()] for l in open(f'{tmp}/g.b') if l.strip() and not l.startswith('#')])
        bv = gb[:, 3]
        pet = mrinfo(in_mif, '-petable')
        rows = [tuple(l.split()) for l in pet.splitlines() if l.strip()] if pet else [('0',)] * size[3]
        if len(rows) != size[3]:
            rows = [('0',)] * size[3]
        groups = {}
        for i, r in enumerate(rows): groups.setdefault(r, []).append(i)
        run(['mrconvert', in_mif, f'{tmp}/in.nii.gz', '-strides', '1,2,3,4', '-force', '-quiet'])
        img = nib.load(f'{tmp}/in.nii.gz'); d = np.asarray(img.dataobj, dtype=np.float32)
        gain4 = np.ones(d.shape, np.float32)
        any_corr = False
        profiles = []
        for gi, (key, idx) in enumerate(sorted(groups.items(), key=lambda kv: -len(kv[1]))):
            idx = np.array(idx); bg = bv[idx]
            i0, idw = idx[bg < 50], idx[bg >= 50]
            rec = dict(pe_row=' '.join(key), n_vol=int(len(idx)), n_b0=int(len(i0)), n_dwi=int(len(idw)))
            if len(idx) < 2:
                rec['decision'] = 'skipped_too_few_volumes'; qc['groups'].append(rec); continue
            have_dwi = len(idw) > 0
            allm = d[..., idx].mean(-1)
            b0 = d[..., i0].mean(-1) if len(i0) else (d[..., idw].mean(-1))
            dw = d[..., idw].mean(-1) if have_dwi else b0
            m, adc, par, csf = masks(b0, dw, have_dwi and len(i0) > 0)
            sel = par if par.sum() > 1000 else m
            cand = {k: dip_at(allm, sel, k) for k in (h - 1, h)}
            k = min(cand, key=cand.get)
            dip = cand[k]
            rec.update(seam_slice=int(k), dip_all=round(dip, 4), dip_candidates={str(a): round(b, 4) for a, b in cand.items()})
            if len(i0) and have_dwi:
                a_amp = alternation(b0, dw, m)
                rec['b0_alt'] = round(a_amp, 4)
                if a_amp > ALT_THR:
                    qc['flags'].append(f'group {gi}: b0 slice alternation {a_amp:.3f} > {ALT_THR} (corrupted single b0? not corrected here)')
            if not have_dwi or len(i0) == 0:
                qc['flags'].append(f'group {gi}: {"no DWI" if not have_dwi else "no b0"} volumes, gain estimated from the mean of all volumes only')
            if dip >= GATE_FLAG:
                rec['decision'] = 'pass'
            else:
                rec['decision'] = 'correct' if dip < GATE_CORRECT else 'correct_flag'
                if rec['decision'] == 'correct_flag':
                    qc['flags'].append(f'group {gi}: dip {dip:.3f} in the grey zone {GATE_CORRECT}-{GATE_FLAG}')
                if cand[h if k == h - 1 else h - 1] < 0.97:
                    qc['flags'].append(f'group {gi}: dips at BOTH {h-1} and {h}, check the slice order')
                G = seam_gain(dw, b0, m, adc, k)
                gain4[:, :, k, idx] = (1.0 / G)[..., None]
                rec['gain'] = dict(p10=round(float(np.percentile(G[m[..., k]], 10)), 3),
                                   p50=round(float(np.percentile(G[m[..., k]], 50)), 3),
                                   p90=round(float(np.percentile(G[m[..., k]], 90)), 3))
                after = allm.copy(); after[..., k] /= G
                rec['dip_all_after'] = round(dip_at(after, sel, k), 4)
                any_corr = True
                profiles.append((gi, dip_profile(allm, sel), dip_profile(after, sel), G, k))
            if rec['decision'] == 'pass':
                profiles.append((gi, dip_profile(allm, sel), None, None, k))
            qc['groups'].append(rec)

        if not any_corr:
            _fig(qcd, profiles, h)
            finish('pass', 'no group has a dip below %.2f' % GATE_FLAG)
        if dry:
            _fig(qcd, profiles, h)
            finish('corrected', 'dry run: nothing written')
        nib.save(nib.Nifti1Image(gain4, img.affine, img.header), f'{tmp}/gain4.nii.gz')
        run(['mrcalc', in_mif, f'{tmp}/gain4.nii.gz', '-mult', out_mif, '-datatype', 'float32', '-force', '-quiet'])
        _fig(qcd, profiles, h)
        finish('corrected')


def _fig(qcd, profiles, h):
    try:
        import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    except Exception:
        return
    if not profiles: return
    fig, ax = plt.subplots(1, 2, figsize=(11, 3.6), gridspec_kw=dict(width_ratios=[2, 1]))
    for gi, before, after, G, k in profiles:
        z = np.arange(len(before))
        ax[0].plot(z, before, '-o', ms=3, label=f'group {gi} before')
        if after is not None: ax[0].plot(z, after, '-o', ms=3, label=f'group {gi} after')
    ax[0].axhline(1, color='gray', lw=.5); ax[0].axhline(GATE_CORRECT, color='r', lw=.5, ls=':')
    ax[0].set_xlabel('slice'); ax[0].set_ylabel('slice / mean of neighbours (parenchyma)'); ax[0].legend(fontsize=7)
    ax[0].set_title('seam screening: slice-to-neighbour ratio, mean of all volumes')
    g = next((p[3] for p in profiles if p[3] is not None), None)
    if g is not None:
        im = ax[1].imshow(g.T, origin='lower', cmap='viridis', vmin=0.6, vmax=1.0); plt.colorbar(im, ax=ax[1], fraction=.046)
        ax[1].set_title('estimated gain, seam slice'); ax[1].axis('off')
    else:
        ax[1].axis('off'); ax[1].text(.5, .5, 'no correction applied', ha='center')
    fig.tight_layout(); fig.savefig(f'{qcd}/seam_gate.png', dpi=110); plt.close(fig)


if __name__ == '__main__':
    main()
