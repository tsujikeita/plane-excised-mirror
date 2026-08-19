# -*- coding: utf-8 -*-
"""Polarization benchmarks of Sec. 4 (verified 2026-08-15).
Computes, in the cosmic-variance limit: H_corr (TE-correlation-only) and
H_geo (E inherits the T suppression 0.10 at l=2-4) expectations for the
fixed-axis statistic S+_E, and detection probabilities.
Requires: healpy, camb, and the CMBanom repository (masks + SMICA map):
  git clone --depth 1 https://github.com/LauraHerold/CMBanom.git
Outputs: benchmark numbers to stdout and figures/a3_fig5_polbench.{png,pdf}.
"""
import numpy as np, healpy as hp, camb, sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))
import phase2_core as p2, plane_mirror as pm
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

pars = camb.CAMBparams()
pars.set_cosmology(H0=67.36, ombh2=0.02237, omch2=0.1200, tau=0.0544)  # verify vs Planck 2018
pars.InitPower.set_params(As=np.exp(3.044)*1e-10, ns=0.9649)
pars.set_for_lmax(64, lens_potential_accuracy=1)
cls = camb.get_results(pars).get_cmb_power_spectra(pars, CMB_unit='muK',
                                                   raw_cl=True)['lensed_scalar']
LM = 32
CTT, CEE, CTE = cls[:LM+1,0], cls[:LM+1,1], cls[:LM+1,3]
src = hp.read_map('CMBanom/data/real/map_smica_nside_128.fits')
if np.std(src) < 5e-3: src *= 1e6
bl_src = hp.gauss_beam(np.radians(1.0), lmax=LM) * p2.pixwin_pad(128, LM)
aT = hp.almxfl(hp.map2alm(src, lmax=LM), 1.0/np.maximum(bl_src, 1e-12))
aEp = hp.almxfl(aT.copy(), np.where(CTT>0, CTE/np.maximum(CTT,1e-30), 0.0))
CEc = np.maximum(CEE - np.where(CTT>0, CTE**2/np.maximum(CTT,1e-30), 0.0), 0.0)
bl = p2.transfer(16, 'planck', LM)

for K in ['common', 'full']:
    mask = p2.make_mask(16, K)
    ms = p2.MirrorStat(16, mask)
    fa = pm.FixedAxisMirror(ms, 1134)                    # frozen axis
    es = lambda a: fa.S_plus(hp.alm2map(hp.almxfl(a, bl), 16), mondip=False)
    Sn, Sc, Sg = [], [], []
    w = np.ones(LM+1); w[2:5] = np.sqrt(0.10)            # H_geo suppression
    for s in range(1000):
        np.random.seed(50_000+s); Sn.append(es(hp.synalm(CEE, lmax=LM)))
        np.random.seed(60_000+s); Sc.append(es(aEp + hp.synalm(CEc, lmax=LM)))
        np.random.seed(70_000+s); Sg.append(es(hp.almxfl(hp.synalm(CEE, lmax=LM), w)))
    Sn, Sc, Sg = map(np.array, (Sn, Sc, Sg))
    for lab, S in [('H_corr', Sc), ('H_geo', Sg)]:
        z = (np.median(Sn)-np.median(S))/np.std(Sn)
        det1 = (S < np.percentile(Sn, 1)).mean()
        print(f'{K:6s} {lab}: z={z:.2f}  median p={max((Sn<np.median(S)).mean(),1e-3):.3f}  P(p<0.01)={det1:.2f}')
    if K == 'common':
        fig, ax = plt.subplots(figsize=(7,4))
        bins = np.linspace(0, 0.06, 50)
        ax.hist(Sn, bins=bins, alpha=0.5, label='null (Gaussian LCDM E)', color='gray')
        ax.hist(Sc, bins=bins, alpha=0.5, label='H_corr: via TE only', color='C0')
        ax.hist(Sg, bins=bins, alpha=0.5, label='H_geo: E suppressed like T', color='C3')
        ax.set_xlabel('S+_E at the frozen axis (CV limit)')
        ax.set_ylabel('realizations'); ax.legend(fontsize=8)
        os.makedirs('figures', exist_ok=True)
        fig.savefig('figures/a3_fig5_polbench.png', dpi=140, bbox_inches='tight')
        fig.savefig('figures/a3_fig5_polbench.pdf', bbox_inches='tight')
