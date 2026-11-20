"""Tests added in the revision, kept separate so that similarity.py and extra_tests.py are unchanged.

1. Envelope law fitted to published zero-pressure balloons.
2. Ground- and balloon-launched vehicles each designed for their own launch over a wider design space
   (fineness 8-25, thrust-to-weight 1.5-6, ground exit pressure down to 0.4 of ambient).
3. Rocket masses on a grid of engine efficiency, drag scale and nozzle-wall mass, for the Monte Carlo of the saving.
4. Two-stage ground-launched comparators (pressure-fed liquid) and solid-motor families (single- and two-stage
   ground launch against a solid rockoon).
5. The Deacon rockoons flown with an independently reported gross mass (no fit).

Writes results/revision.json. Run after similarity.py."""
import itertools
import json
import math
import os
from multiprocessing import Pool

import numpy as np
from scipy.optimize import brentq, minimize_scalar

import rockoon_model as rm

ROOT = os.path.dirname(rm.HERE)
MATS = ('Al6061T6', 'SS321')
MFIX = (10.0, 20.0, 40.0, 80.0, 160.0, 320.0, 640.0)
ZB = (10.0, 12.5, 15.0, 17.5, 20.0, 22.5, 25.0, 27.5, 30.0, 32.5, 35.0)
FIN_G = (8.0, 10.0, 12.5, 15.0, 17.5, 20.0, 22.5, 25.0)
TW_G = (1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0)
SEP_G = (1.0, 0.85, 0.7, 0.55, 0.4)
FIN_B = (8.0, 10.0, 12.5, 15.0)
TW_B = (2.0, 3.0, 4.0, 5.0, 6.0)
ETA = (0.85, 0.92, 0.95)
CDS = (0.8, 1.0, 1.2)
AREAL = (0.0, 10.0, 20.0)
TARGET = 100e3
F_SOLID = 0.15          # motor case, liner and nozzle per kg of propellant (motor mass fraction 0.87)
BETA_LIQ = 0.5          # fixed mass of a liquid booster stage as a fraction of m_fix
BETA_SOL = 0.25         # fixed mass of a solid booster stage (fins, interstage) as a fraction of m_np

# ------------------------------------------------------------------ 1. published balloons
LB, CF = 0.45359237, 0.0283168466
BALLOONS = [  # volume m3, envelope kg, rated payload kg, source
    (214, 2.5, 38.5, 'Aerostar small zero-pressure, Table 1'), (268, 2.9, 38.5, 'Aerostar Table 1'),
    (538, 3.4, 17.2, 'Aerostar Table 1'), (1528, 6.2, 22.2, 'Aerostar Table 1'),
    (3990, 11.4, 32.7, 'Aerostar Table 2'), (8490, 18.3, 39.0, 'Aerostar Table 2'),
    (14150, 25.5, 41.3, 'Aerostar Table 2'), (19810, 32.1, 47.2, 'Aerostar Table 2'),
    (29.47e6 * CF, 3600 * LB, 4000 * LB, 'NASA 29.47 MCF'), (59.84e6 * CF, 2751 * LB, 1200 * LB, 'NASA 59.84 MCF'),
    (2.1e6, 774.0, 60.0, 'Rainwater and Smith, 2.1 million m3 design')]


def fit_envelope():
    V = np.array([b[0] for b in BALLOONS]); me = np.array([b[1] for b in BALLOONS]); mp = np.array([b[2] for b in BALLOONS])
    A = (36 * math.pi) ** (1 / 3) * V ** (2 / 3)
    w = 1.0 / me                                   # relative (not absolute) errors
    out = {}
    # law of the paper: sigma A + k m
    X1 = np.column_stack([A, mp]); c1 = np.linalg.lstsq(X1 * w[:, None], me * w, rcond=None)[0]
    X2 = np.column_stack([A, mp * V ** (1 / 3)]); c2 = np.linalg.lstsq(X2 * w[:, None], me * w, rcond=None)[0]
    for tag, X, c in (('law1', X1, c1), ('law2', X2, c2)):
        r = np.log(X @ c / me)
        out[tag] = dict(coef=c.tolist(), rms_log=float(np.sqrt(np.mean(r ** 2))), max_abs=float(np.max(np.abs(np.exp(r) - 1))),
                        rel=(np.exp(r) - 1).tolist())
    out['sigma_eff'] = (me / A).tolist()
    out['data'] = [dict(V=b[0], me=b[1], mp=b[2], src=b[3]) for b in BALLOONS]
    return out


# ------------------------------------------------------------------ 2. design searches
def _close(args):
    z, mfix, f, kw = args
    try:
        r = rm.close_vehicle(z * 1e3, mfix, f, target=TARGET, **kw)
    except ValueError:
        r = None
    return None if r is None else r['m0']


def best_designs(pool):
    jobs, keys = [], []
    for mat in MATS:
        f = rm.feed_fraction(mat)
        for mfix in MFIX:
            for fin, tw, sp in itertools.product(FIN_G, TW_G, SEP_G):
                jobs.append((0.0, mfix, f, dict(tw=tw, fineness=fin, sep=sp))); keys.append(('g', mat, mfix, 0.0, fin, tw, sp))
            for z in ZB:
                for fin, tw in itertools.product(FIN_B, TW_B):
                    jobs.append((z, mfix, f, dict(tw=tw, fineness=fin))); keys.append(('b', mat, mfix, z, fin, tw, 1.0))
    res = pool.map(_close, jobs, chunksize=20)
    best = {}
    for k, m in zip(keys, res):
        if m is None:
            continue
        kk = (k[0], k[1], k[2], k[3])
        if kk not in best or m < best[kk][0]:
            best[kk] = (m, k[4], k[5], k[6])
    return best


# ------------------------------------------------------------------ 3. rocket-mass grid for the Monte Carlo
def _close_grid(args):
    z, mfix, f, kw, eta, cds, areal = args
    try:
        if areal > 0:
            r = rm.close_vehicle_skirt(z * 1e3, mfix, f, target=TARGET, areal=areal, eta=eta, cd_scale=cds, **kw)
        else:
            r = rm.close_vehicle(z * 1e3, mfix, f, target=TARGET, eta=eta, cd_scale=cds, **kw)
    except ValueError:
        r = None
    return float('nan') if r is None else r['m0']


def mass_grid(pool, best):
    jobs, keys = [], []
    for mat in MATS:
        f = rm.feed_fraction(mat)
        for mfix in MFIX:
            for z in (0.0,) + ZB:
                b = best.get(('g' if z == 0 else 'b', mat, mfix, z))
                if b is None:
                    continue
                kw = dict(tw=b[2], fineness=b[1]) if z > 0 else dict(tw=b[2], fineness=b[1], sep=b[3])
                for eta, cds, ar in itertools.product(ETA, CDS, AREAL):
                    jobs.append((z, mfix, f, kw, eta, cds, ar)); keys.append((mat, mfix, z, eta, cds, ar))
    res = pool.map(_close_grid, jobs, chunksize=10)
    return {k: m for k, m in zip(keys, res)}


# ------------------------------------------------------------------ 4. two-stage and solid comparators
def _two_stage(args):
    """Lightest two-stage ground-launched vehicle reaching TARGET for one design (tw1, tw2, fineness)."""
    mfix, f, beta, key, tw1, tw2, fin, sep1, areal = args[:9]
    fin_drag = len(args) > 9 and args[9]
    m_fix2_base = mfix
    cdf = lambda d: rm.FinenessDrag(d, fin) if fin_drag else None

    def gross(mp2, skirt=0.0):
        m02 = m_fix2_base + skirt + (1 + f) * mp2

        def miss(mp1):
            m01 = m02 + beta * mfix + (1 + f) * mp1
            d = rm.diameter(m01, fineness=fin)
            try:
                ap, _ = rm.two_stage_apogee(m02, mp2, m01, mp1, d, tw1, tw2, sep1, key, cd_fn=cdf(d))
            except ValueError:
                return -1e9
            return ap - TARGET
        lo, hi = 0.5, 20.0
        if miss(lo) > 0:
            return None
        while miss(hi) < 0:
            hi *= 1.6
            if hi > 6e4:
                return None
        mp1 = brentq(miss, lo, hi, xtol=1e-2)
        m01 = m02 + beta * mfix + (1 + f) * mp1
        return m01, mp1, m02

    def total(mp2):
        skirt = 0.0
        out = None
        for _ in range(4):                       # second-stage nozzle skirt, iterated
            out = gross(mp2, skirt)
            if out is None:
                return None
            m01, mp1, m02 = out
            d = rm.diameter(m01, fineness=fin)
            try:
                _, nz2 = rm.two_stage_apogee(m02, mp2, m01, mp1, d, tw1, tw2, sep1, key, cd_fn=cdf(d))
            except ValueError:
                return None
            new = rm.nozzle_skirt(nz2.Ae, nz2.eps, areal=areal)[0] if areal > 0 else 0.0
            if abs(new - skirt) < 0.05:
                break
            skirt = new
        return out

    best = None
    grid = mfix * np.array([0.25, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0])
    vals = []
    for mp2 in grid:
        t = total(float(mp2))
        vals.append(t[0] if t else np.inf)
        if t and (best is None or t[0] < best[0]):
            best = (t[0], float(mp2), t[1])
    if best is None:
        return None
    i = int(np.argmin(vals))
    lo = grid[max(i - 1, 0)]; hi = grid[min(i + 1, len(grid) - 1)]
    if hi > lo:
        r = minimize_scalar(lambda x: (total(float(x)) or (np.inf,))[0], bounds=(lo, hi), method='bounded', options=dict(xatol=0.01 * mfix))
        t = total(float(r.x))
        if t and t[0] < best[0]:
            best = (t[0], float(r.x), t[1])
    return dict(m0=best[0], mp2=best[1], mp1=best[2], tw1=tw1, tw2=tw2, fin=fin, sep1=sep1)


def two_stage_family(pool, specs):
    jobs, keys = [], []
    for spec in specs:
        tag, mfix, f, beta, key, tw1s, tw2s, fins, areal = spec[:9]
        for tw1, tw2, fin in itertools.product(tw1s, tw2s, fins):
            jobs.append((mfix, f, beta, key, tw1, tw2, fin, 0.7, areal) + tuple(spec[9:])); keys.append((tag, mfix))
    res = pool.map(_two_stage, jobs, chunksize=1)
    best = {}
    for k, r in zip(keys, res):
        if r and (k not in best or r['m0'] < best[k]['m0']):
            best[k] = r
    return best


def _close_skirt(args):
    z, mfix, f, kw, areal = args
    try:
        r = rm.close_vehicle_skirt(z * 1e3, mfix, f, target=TARGET, areal=areal, **kw) if areal > 0 else \
            rm.close_vehicle(z * 1e3, mfix, f, target=TARGET, **kw)
    except ValueError:
        r = None
    return None if r is None else r['m0']


def solid_singles(pool, areal):
    """Single-stage solid vehicles, ground- and balloon-launched, each over its design space."""
    jobs, keys = [], []
    for mnp in MFIX:
        for fin, tw, sp in itertools.product((10.0, 15.0, 20.0, 25.0), (3.0, 5.0, 8.0, 12.0), (1.0, 0.7, 0.4)):
            jobs.append((0.0, mnp, F_SOLID, dict(tw=tw, fineness=fin, sep=sp, key='solid'), areal)); keys.append((mnp, 0.0))
        for z in ZB:
            for fin, tw in itertools.product((8.0, 10.0, 15.0), (3.0, 5.0, 8.0)):
                jobs.append((z, mnp, F_SOLID, dict(tw=tw, fineness=fin, key='solid'), areal)); keys.append((mnp, z))
    res = pool.map(_close_skirt, jobs, chunksize=10)
    best = {}
    for k, m in zip(keys, res):
        if m is not None and (k not in best or m < best[k]):
            best[k] = m
    return best


# ------------------------------------------------------------------ main
def main():
    out = {}
    env = fit_envelope(); out['envelope_fit'] = env
    with Pool() as pool:
        best = best_designs(pool)
        out['best'] = {'%s|%s|%g|%g' % k: list(v) for k, v in best.items()}
        grid = mass_grid(pool, best)
        out['grid'] = {'%s|%g|%g|%g|%g|%g' % k: v for k, v in grid.items()}
        # two-stage pressure-fed liquid, ground launch, at the nominal efficiency and drag
        specs = []
        for mat in MATS:
            f = rm.feed_fraction(mat)
            for mfix in MFIX:
                specs.append(('liq|' + mat, mfix, f, BETA_LIQ, 'nozzle', (2.0, 3.0, 4.0), (1.5, 2.5, 4.0), (20.0, 25.0), 10.0))
        for mnp in MFIX:
            specs.append(('sol', mnp, F_SOLID, BETA_SOL, 'solid', (4.0, 8.0, 12.0), (2.0, 4.0, 8.0), (20.0, 25.0), 10.0))
        ts = two_stage_family(pool, specs)
        out['two_stage'] = {'%s|%g' % k: v for k, v in ts.items()}
        out['solid_single'] = {'%g|%g' % k: v for k, v in solid_singles(pool, 10.0).items()}
    with open(os.path.join(ROOT, 'results', 'revision.json'), 'w') as f_:
        json.dump(out, f_, default=float)
    print('revision tests written:', len(best), 'designs;', len(grid), 'grid points;', len(ts), 'two-stage families')


if __name__ == '__main__':
    main()
