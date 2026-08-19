# -*- coding: utf-8 -*-
"""plane_mirror.py — テーマA: 鏡映反対称性の高速バッチ評価
MirrorStat（phase2_core）と厳密同一の定義で，方向走査をマップ束一括化。
検証済み（2026-08-15）: MirrorStatとの max相対差 1.2e-6（float32水準）。
計時: N16 33ms/マップ（10^5=55分）, N32 0.9s/マップ（10^4=2.5時間）。
"""
import numpy as np
import healpy as hp


class MirrorBatch:
    def __init__(self, ms):
        self.R, self.valid, self.cnt, self.mask = ms.R, ms.valid, ms.cnt, ms.mask

    def min_S_batch(self, maps, mondip=True, return_argmin=False):
        B = maps.shape[0]
        T = np.empty((maps.shape[1], B), np.float32)
        for b in range(B):
            m = maps[b]
            if mondip:
                m = np.asarray(hp.remove_dipole(hp.ma(np.where(self.mask, m, hp.UNSEEN))))
            T[:, b] = np.where(self.mask, m, 0.0)
        nd = self.R.shape[0]
        minSp = np.full(B, np.inf, np.float32); minSm = np.full(B, np.inf, np.float32)
        argp = np.zeros(B, np.int32); argm = np.zeros(B, np.int32)
        for d in range(nd):
            v = self.valid[d]
            if not v.any():
                continue
            Tr = T[self.R[d]]
            Tv, Trv = T[v], Tr[v]
            Sp = np.einsum('ib,ib->b', 0.5 * (Tv + Trv), 0.5 * (Tv + Trv)) / self.cnt[d]
            Sm = np.einsum('ib,ib->b', 0.5 * (Tv - Trv), 0.5 * (Tv - Trv)) / self.cnt[d]
            better = Sp < minSp
            if return_argmin:
                argp = np.where(better, d, argp)
            minSp = np.where(better, Sp, minSp)
            better = Sm < minSm
            if return_argmin:
                argm = np.where(better, d, argm)
            minSm = np.where(better, Sm, minSm)
        if return_argmin:
            return minSp, minSm, argp, argm
        return minSp, minSm


def axis_lb(nside, d):
    """方向画素番号 → 銀経緯 (l, b) [deg]"""
    th, ph = hp.pix2ang(nside, int(d))
    return float(np.degrees(ph)), float(90.0 - np.degrees(th))


class FixedAxisMirror:
    """凍結軸 d* での対称成分解析（O(npix)/マップ）"""
    def __init__(self, ms, d_star):
        self.r = ms.R[d_star]
        self.v = ms.valid[d_star]
        self.cnt = ms.cnt[d_star]
        self.mask = ms.mask

    def S_plus(self, mp, mondip=True):
        if mondip:
            mp = np.asarray(hp.remove_dipole(hp.ma(np.where(self.mask, mp, hp.UNSEEN))))
        T = np.where(self.mask, mp, 0.0)
        s = 0.5 * (T + T[self.r])
        return float(np.sum(self.v * s * s) / self.cnt)

    def S_plus_batch(self, maps, mondip=True):
        return np.array([self.S_plus(maps[b], mondip) for b in range(maps.shape[0])])

    def band_decompose(self, alm128, nside, transfer_fl, bands, mondip=True):
        """ソースalm(lmax128)を帯域分解し，帯域別S⁺と全帯域和・交差項を返す"""
        import numpy as _np
        LMAX = hp.Alm.getlmax(len(alm128))
        ell = _np.arange(LMAX + 1)
        out = {}
        s_parts = []
        for (l0, l1) in bands:
            w = ((ell >= l0) & (ell <= l1)).astype(float) * transfer_fl
            mb = hp.alm2map(hp.almxfl(alm128.copy(), w), nside)
            if mondip and l0 <= 1:
                pass
            T = _np.where(self.mask, mb, 0.0)
            if mondip:
                T = _np.where(self.mask,
                              _np.asarray(hp.remove_dipole(hp.ma(_np.where(self.mask, mb, hp.UNSEEN)))), 0.0)
            s = 0.5 * (T + T[self.r])
            s_parts.append(s)
            out[f'S{l0}_{l1}'] = float(_np.sum(self.v * s * s) / self.cnt)
        stot = _np.sum(s_parts, axis=0)
        out['S_sum_bands'] = float(_np.sum(self.v * stot * stot) / self.cnt)
        return out

    def exclusion_scan(self, mp, nside_scan=8, radius_deg=15.0, mondip=True):
        """半径radius_degの円盤を各走査位置で（鏡映相手も対称に）除外した際の
        ln S⁺ の変化 Δ(q) を返す（正＝除外で対称性が増える＝その領域が反対称の担い手）"""
        if mondip:
            mp = np.asarray(hp.remove_dipole(hp.ma(np.where(self.mask, mp, hp.UNSEEN))))
        T = np.where(self.mask, mp, 0.0)
        s = 0.5 * (T + T[self.r])
        base_num = np.sum(self.v * s * s)
        base = base_num / self.cnt
        npix_scan = hp.nside2npix(nside_scan)
        nside_map = hp.npix2nside(len(T))
        delta = np.zeros(npix_scan)
        for q in range(npix_scan):
            vec = hp.pix2vec(nside_scan, q)
            disc = hp.query_disc(nside_map, vec, np.radians(radius_deg))
            ex = np.zeros(len(T), bool); ex[disc] = True
            ex = ex | ex[self.r]                       # 対称除外
            v2 = self.v & ~ex
            c2 = max(v2.sum(), 1)
            S2 = np.sum(v2 * s * s) / c2
            delta[q] = np.log(S2 / base) if S2 > 0 else 0.0
        return delta, base


def with_mask(ms, mask):
    """R表を共有して別マスクのMirrorStat相当を作る（N32のR再構築20秒を節約）"""
    obj = type(ms).__new__(type(ms))
    obj.R = ms.R
    obj.mask = mask
    obj.valid = mask[None, :] & mask[ms.R]
    obj.cnt = np.maximum(obj.valid.sum(axis=1), 1)
    return obj


def scan_S(ms, mp, mondip=True):
    """1マップの全方向S±(n̂)地形を返す（軸の縮退・地形幅の解析用）"""
    if mondip:
        mp = np.asarray(hp.remove_dipole(hp.ma(np.where(ms.mask, mp, hp.UNSEEN))))
    T = np.where(ms.mask, mp, 0.0).astype(np.float32)
    nd = ms.R.shape[0]
    Sp = np.empty(nd, np.float32); Sm = np.empty(nd, np.float32)
    for d in range(nd):
        v = ms.valid[d]
        Tr = T[ms.R[d]]
        Sp[d] = np.sum(v * (0.5 * (T + Tr)) ** 2) / ms.cnt[d]
        Sm[d] = np.sum(v * (0.5 * (T - Tr)) ** 2) / ms.cnt[d]
    return Sp, Sm


def axis_sep_deg(nside, d1, d2):
    """2軸の分離角[deg]（鏡映面法線は±同一視 → 90°超は補角）"""
    v1 = np.array(hp.pix2vec(nside, int(d1)))
    v2 = np.array(hp.pix2vec(nside, int(d2)))
    ang = np.degrees(np.arccos(np.clip(abs(v1 @ v2), -1, 1)))
    return float(ang)
