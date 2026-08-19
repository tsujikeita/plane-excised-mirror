# -*- coding: utf-8 -*-
"""Regenerate paper figures 1-4 from data/ (fig5: scripts/forecast_polarization.py)."""
import numpy as np, pandas as pd, healpy as hp, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))
import plane_mirror as pm
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
D = os.path.join(os.path.dirname(__file__), '..', 'data')
F = os.path.join(os.path.dirname(__file__), '..', 'figures'); os.makedirs(F, exist_ok=True)
def dual(fig, name):
    fig.savefig(f'{F}/{name}.png', dpi=140, bbox_inches='tight')
    fig.savefig(f'{F}/{name}.pdf', bbox_inches='tight')

ax_df = pd.read_csv(f'{D}/a1_axes.csv')
fig, ax = plt.subplots(figsize=(6,5))
mk = {'full':'o','common':'s','ext':'^'}
for K in mk:
    for n, col in [(16,'C0'),(32,'C3')]:
        s = ax_df[(ax_df['mask']==K)&(ax_df.nside==n)]
        ll = np.where(s.l>180, s.l-360, s.l)
        lp = np.where(s.b<0, ll-180, ll); bp = np.abs(s.b)
        ax.scatter(lp, bp, marker=mk[K], s=40, facecolors='none', edgecolors=col,
                   label=f'N{n}/{K}')
ax.set_xlabel('l [deg]'); ax.set_ylabel('|b| [deg]'); ax.legend(fontsize=7, ncol=2)
ax.set_title('Best anti-symmetry axes: 10 maps x 6 configurations')
dual(fig, 'a3_fig1_axes'); plt.close()

bd = pd.read_csv(f'{D}/a1_bands.csv')
fig, ax = plt.subplots(figsize=(7,4.5))
order = ['S2_4','S5_8','S9_16','S17_32','S33_64']; xx = np.arange(5)
for mp_, g in bd.groupby('map'):
    y = [g[g.band==b_].ratio.iloc[0] for b_ in order]
    st = '-' if mp_.startswith('PR3') else ('--' if mp_.startswith('PR4') else ':')
    ax.plot(xx, y, st, marker='o', ms=3, alpha=0.7, label=mp_)
ax.axhline(1, color='gray', lw=0.7)
ax.set_xticks(xx); ax.set_xticklabels(['2-4','5-8','9-16','17-32','33-64'])
ax.set_xlabel('multipole band'); ax.set_ylabel('symmetric power / null median')
ax.legend(fontsize=6, ncol=2)
dual(fig, 'a3_fig2_bands'); plt.close()

z = np.load(f'{D}/a1_exclusion.npz')
plt.figure(figsize=(7,4.5))
hp.mollview(z['PR3_SMICA'], title='Exclusion scan (SMICA)', unit='Delta ln S+',
            cmap='RdBu_r', min=-0.1, max=0.1, hold=True)
hp.graticule(alpha=0.3)
dual(plt.gcf(), 'a3_fig3_exclusion'); plt.close()

fb = pd.read_csv(f'{D}/a1_floorbreak.csv')
fig, ax = plt.subplots(figsize=(7.5,4))
for i, (cid, g) in enumerate(fb.groupby('config')):
    ax.errorbar(np.arange(len(g))+i*0.18-0.27, g.p, fmt='o', ms=4, alpha=0.85,
                label=cid.replace('_mdON_harm','').replace('_Splanck',''))
ax.axhline(1e-3, color='gray', ls=':', lw=0.8)
ax.set_yscale('log'); ax.set_ylim(8e-5, 3e-3)
maps = list(fb[fb.config==fb.config.unique()[0]]['map'])
ax.set_xticks(range(len(maps)))
ax.set_xticklabels([m.replace('PR3_','').replace('PR4_','PR4 ').replace('Nofi_','')
                    for m in maps], rotation=45, fontsize=7)
ax.set_ylabel('p(min S+)'); ax.legend(fontsize=8)
dual(fig, 'a3_fig4_floorbreak'); plt.close()
print('figures 1-4 regenerated')
