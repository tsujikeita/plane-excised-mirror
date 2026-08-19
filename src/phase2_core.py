# -*- coding: utf-8 -*-
"""Phase 2 較正基盤コア（v0.1）
- 処理構成マニフェスト（計画書v1.0 §3準拠）
- CRNマスター実現（synalm, seed 0..999, lmax=128）
- 統計①②③④の null 分布計算（チェックポイント/再開対応）
"""
import os, json, time
import numpy as np, healpy as hp

# ---------- 基本設定 ----------
LMAX_MASTER = 128
NSIM_FULL = 1000
FID_CL_FILE = 'CMBanom/data/real/COM_PowerSpect_CMB-base-plikHM-TTTEEE-lowl-lowE-lensing-minimum-theory_R3.01.txt'
COMMON_MASK_128 = 'CMBanom/data/masks/com_mask_cutoff_0.9_nside_128.fits'

def load_fid_cl(lmax=LMAX_MASTER):
    dat = np.loadtxt(FID_CL_FILE, skiprows=1)
    ll = np.arange(lmax + 1); cl = np.zeros(lmax + 1)
    n = min(lmax - 1, dat.shape[0])
    cl[2:2 + n] = dat[:n, 1] * 2 * np.pi / (ll[2:2 + n] * (ll[2:2 + n] + 1))
    return cl

# ---------- 伝達関数・マスク ----------
PIXWIN_CACHE = 'pixwin_cache'
def pixwin_pad(nside, lmax):
    fn = os.path.join(PIXWIN_CACHE, f'pixel_window_n{nside:04d}.fits')
    if os.path.exists(fn):
        from astropy.io import fits as _f
        with _f.open(fn) as h:
            pw = np.asarray(h[1].data['TEMPERATURE']).ravel()
    else:
        pw = hp.pixwin(nside)   # Colabでは通常経路（初回のみDL）
    return np.pad(pw, (0, max(0, lmax + 1 - len(pw))), mode='edge')[:lmax + 1]

def transfer(nside, smooth, lmax=LMAX_MASTER):
    """構成の伝達関数 b_ℓ p_ℓ。smooth ∈ {'planck','none','fix5deg'}"""
    pw = pixwin_pad(nside, lmax)
    if smooth == 'planck':
        fwhm_arcmin = 640.0 * 16.0 / nside
        return hp.gauss_beam(np.radians(fwhm_arcmin / 60.), lmax=lmax) * pw
    if smooth == 'fix5deg':
        return hp.gauss_beam(np.radians(5.0), lmax=lmax) * pw
    if smooth == 'none':
        return pw.copy()
    raise ValueError(smooth)

def dilate_mask(bad, nside, deg):
    """マスク（bad=True）を約deg度拡張（近傍膨張の反復）"""
    pixsize_deg = np.degrees(hp.nside2resol(nside))
    n_iter = max(1, int(np.ceil(deg / pixsize_deg)))
    bad = bad.copy()
    for _ in range(n_iter):
        idx = np.where(bad)[0]
        nb = hp.get_all_neighbours(nside, idx)
        bad[nb[nb >= 0]] = True
    return bad

def _dilate_n(bad, nside, n_iter):
    bad = bad.copy()
    for _ in range(n_iter):
        idx = np.where(bad)[0]
        nb = hp.get_all_neighbours(nside, idx)
        bad[nb[nb >= 0]] = True
    return bad

def _erode_n(bad, nside, n_iter):
    return ~_dilate_n(~bad, nside, n_iter)

def make_mask(nside, level):
    """level ∈ {'full','common','ext'} → bool（True=使用画素）
    ext = Planck 2015 XVI Table 12方式：共通マスクの拡散（銀河）成分のみを5°拡張し
          点源穴は拡張しない（補助マスク規定）。128で構成してから縮退。"""
    npix = hp.nside2npix(nside)
    if level == 'full':
        return np.ones(npix, bool)
    m128 = hp.read_map(COMMON_MASK_128)
    if level == 'common':
        return hp.ud_grade(m128, nside) >= 0.9
    if level == 'ext':
        bad128 = m128 < 0.5
        pix_deg = np.degrees(hp.nside2resol(128))          # ≈0.46°
        n_open = 2                                          # 開演算で点源(≲1°)を除去
        diffuse = _dilate_n(_erode_n(bad128, 128, n_open), 128, n_open)
        n5 = int(np.ceil(5.0 / pix_deg))                    # 5°膨張
        bad_ext = bad128 | _dilate_n(diffuse, 128, n5)
        return hp.ud_grade((~bad_ext).astype(float), nside) >= 0.9
    raise ValueError(level)

# ---------- CRNマスター実現 ----------
def gen_master_alm(outfile, nsim=NSIM_FULL, lmax=LMAX_MASTER):
    if os.path.exists(outfile):
        return np.load(outfile)['alms']
    cl = load_fid_cl(lmax)
    alms = np.empty((nsim, hp.Alm.getsize(lmax)), dtype=np.complex128)
    for s in range(nsim):
        np.random.seed(s)
        alms[s] = hp.synalm(cl, lmax=lmax)
    np.savez_compressed(outfile, alms=alms, lmax=lmax, nsim=nsim,
                        note='CRN master: fiducial PR3 bestfit, seeds 0..nsim-1')
    return alms

def config_maps(alms, nside, smooth, route='harmonic', lmax=LMAX_MASTER):
    """マスター実現→構成マップ群 (nsim, npix)"""
    bl = transfer(nside, smooth, lmax)
    npix = hp.nside2npix(nside)
    out = np.empty((alms.shape[0], npix))
    if route == 'harmonic':
        for s in range(alms.shape[0]):
            out[s] = hp.alm2map(hp.almxfl(alms[s], bl), nside)
    elif route == 'udgrade':   # 高解像度で実体化→画素平均（quadrature軸）
        nhi = min(4 * nside, 64)
        bl_hi = transfer(nside, smooth, lmax) / pixwin_pad(nside, lmax) * pixwin_pad(nhi, lmax)
        for s in range(alms.shape[0]):
            out[s] = hp.ud_grade(hp.alm2map(hp.almxfl(alms[s], bl_hi), nhi), nside)
    else:
        raise ValueError(route)
    return out

# ---------- 統計②：鏡映パリティ（マスク対応） ----------
class MirrorStat:
    def __init__(self, nside, mask):
        npix = hp.nside2npix(nside)
        vecs = np.array(hp.pix2vec(nside, np.arange(npix))).T
        R = np.empty((npix, npix), dtype=np.int32)
        for i in range(npix):
            n = vecs[i]
            refl = vecs - 2.0 * np.outer(vecs @ n, n)
            R[i] = hp.vec2pix(nside, refl[:, 0], refl[:, 1], refl[:, 2])
        self.R = R
        self.valid = mask[None, :] & mask[R]          # 双方が有効な対のみ
        self.cnt = np.maximum(self.valid.sum(axis=1), 1)
        self.mask = mask

    def min_S(self, mp, chunk=1024, mondip=True):
        if mondip:
            mp = np.asarray(hp.remove_dipole(hp.ma(np.where(self.mask, mp, hp.UNSEEN))))
        mp = np.where(self.mask, mp, 0.0).astype(np.float32)
        nd = self.R.shape[0]
        Sp = np.empty(nd, np.float32); Sm = np.empty(nd, np.float32)
        for a in range(0, nd, chunk):
            Tr = mp[self.R[a:a + chunk]]
            V = self.valid[a:a + chunk]
            Sp[a:a + chunk] = np.sum(V * (0.5 * (mp[None, :] + Tr)) ** 2, axis=1) / self.cnt[a:a + chunk]
            Sm[a:a + chunk] = np.sum(V * (0.5 * (mp[None, :] - Tr)) ** 2, axis=1) / self.cnt[a:a + chunk]
        return float(Sp.min()), float(Sm.min())

# ---------- 統計①：R/D（QML / MASTER / naive fsky） ----------
def C_pm(cl_from2, lmax):
    ell = np.arange(2, lmax + 1)
    Dl = ell * (ell + 1.) / (2 * np.pi) * cl_from2[:lmax - 1]
    ev = (ell % 2 == 0)
    return Dl[ev].sum() / (lmax - 1), Dl[~ev].sum() / (lmax - 1)

def RD_traj(cl_from2, lmaxes):
    out = np.empty((len(lmaxes), 2))
    for i, L in enumerate(lmaxes):
        p, m = C_pm(cl_from2, L)
        out[i] = (p / m, p - m)
    return out

def lmax_est_of(nside):
    return int(min(40, 2.5 * nside))

# ---------- 統計③：多重極ベクトル（polyMV非依存・性質テスト済） ----------
from scipy.special import gammaln as _gln

def multipole_vectors(alm, lmax, ell):
    a = np.zeros(2*ell+1, dtype=complex)
    for m in range(0, ell+1):
        v = alm[hp.Alm.getidx(lmax, ell, m)]
        a[ell+m] = v
        if m: a[ell-m] = (-1)**m * np.conj(v)
    k = np.arange(2*ell+1)
    logC = _gln(2*ell+1) - _gln(k+1) - _gln(2*ell-k+1)
    c = np.sqrt(np.exp(logC)) * a
    roots = np.roots(c[::-1])
    th = 2*np.arctan(np.abs(roots)); ph = np.angle(roots) + np.pi
    v = np.stack([np.sin(th)*np.cos(ph), np.sin(th)*np.sin(ph), np.cos(th)], 1)
    v = np.where(v[:, 2:3] >= 0, v, -v)
    keep = []
    for i in range(len(v)):
        if not any(np.dot(v[i], v[j]) > 0.999 for j in keep): keep.append(i)
    return v[keep[:ell]]

def stat3_SQO(maps, mask, lmax_alm=8):
    """素朴カットスカイalm経路（第1弾で検証済みの規約）でS_QO"""
    vals = np.empty(maps.shape[0])
    fmask = mask.astype(float)
    for s in range(maps.shape[0]):
        alm = hp.map2alm(maps[s]*fmask, lmax=lmax_alm, iter=3)
        v2 = multipole_vectors(alm, lmax_alm, 2)
        v3 = multipole_vectors(alm, lmax_alm, 3)
        w2 = np.cross(v2[0], v2[1])
        w3 = [np.cross(v3[i], v3[j]) for i in range(3) for j in range(i+1, 3)]
        vals[s] = np.mean([abs(np.dot(w2, w)) for w in w3])
    return vals

# ---------- 実行系（チェックポイント） ----------
def null_dir(base, cfg_id):
    d = os.path.join(base, cfg_id); os.makedirs(d, exist_ok=True); return d

def true_varlvmap(lvmaps, lvmask, mean_lvmap):
    """正しい逆分散重み用の規格化分散（リポジトリ版get_varlvmapは定数(N-1)²/Nになるバグ）"""
    return np.where(lvmask == 1.,
                    np.mean((lvmaps - mean_lvmap) ** 2, axis=0) / np.maximum(mean_lvmap, 1e-30) ** 2,
                    1.)

def run_stat2(base, cfg_id, maps, nside, mask, chunk_ckpt=100, mondip=True):
    """②のnull分布（部分保存・再開対応）"""
    d = null_dir(base, cfg_id); fn = os.path.join(d, 'stat2.npz')
    done = 0; minSp = []; minSm = []
    if os.path.exists(fn):
        z = np.load(fn)
        if z['complete']: return
        minSp, minSm, done = list(z['minSp']), list(z['minSm']), int(z['done'])
    ms = MirrorStat(nside, mask)
    for s in range(done, maps.shape[0]):
        sp, sm = ms.min_S(maps[s], mondip=mondip)
        minSp.append(sp); minSm.append(sm)
        if (s + 1) % chunk_ckpt == 0 or s == maps.shape[0] - 1:
            np.savez(fn, minSp=minSp, minSm=minSm, done=s + 1,
                     complete=(s == maps.shape[0] - 1))
    return np.array(minSp), np.array(minSm)
