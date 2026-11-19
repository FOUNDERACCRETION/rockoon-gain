"""Final set of tests on the mass-optimal float altitude and on whether a rockoon saves mass.

1. The envelope rule: at the optimum the envelope mass equals (3/2) H_rho |d m_s / dz|; a
   closed-form approximation for z*.
2. Comparators designed for their own launch over a wider search (fineness 8-25, thrust-to-weight
   ratio 1.5-6, ground exit pressure down to 0.4 of ambient), two-stage pressure-fed vehicles, and
   single- and two-stage solid-motor vehicles.
3. The rockoon with the hardware only it needs (divergent-section mass, a heavier gondola for the
   launcher and attitude control, load-scaled envelope).
4. A Monte Carlo over these penalties and over efficiency and drag, giving the probability that the
   rockoon launches less mass.

Writes results/final.json. Kept separate from similarity.py and extra_tests.py so that their
results are unchanged.
"""
import collections
import itertools
import json
import math
import os
from multiprocessing import Pool

import numpy as np
from scipy.optimize import brentq, curve_fit

import rockoon_model as rm

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MATS = ('Al6061T6', 'SS321')
MFIX = (10.0, 20.0, 40.0, 80.0, 160.0, 320.0, 640.0)
ZB = (10.0, 12.5, 15.0, 17.5, 20.0, 22.5, 25.0, 27.5, 30.0, 32.5, 35.0)
TARGET = 100e3
ZMC = (12.5, 15.0, 17.5, 20.0, 22.5, 25.0, 27.5, 30.0)
SOLID_F = 0.15          # case, liner and nozzle per kg of propellant for the generic solid motor
GAMMA_H2 = rm.M_H2 / (rm.M_AIR - rm.M_H2)


# ------------------------------------------------------------------ helpers
def f_of(mat):
    return SOLID_F if mat == 'solid' else rm.feed_fraction(mat)


def key_of(mat):
    return 'solid' if mat == 'solid' else 'nozzle'


def close1(z, mfix, mat, areal=0.0, **kw):
    """Single-stage closure, optionally with the divergent-section mass added (iterated)."""
    kw = dict(kw, key=key_of(mat))
    try:
        if areal > 0:
            return rm.close_vehicle_skirt(z, mfix, f_of(mat), target=TARGET, areal=areal, **kw)
        return rm.close_vehicle(z, mfix, f_of(mat), target=TARGET, **kw)
    except ValueError:
        return None


def two_stage_apogee(m01, mp1, m02, mp2, d, tw1, tw2, key, sep, areal):
    r1 = rm.vehicle_run(m01, mp1, 0.0, d=d, tw=tw1, key=key, sep=sep)
    r2 = rm.vehicle_run(m02, mp2, r1['h_b'], d=d, tw=tw2, key=key, v0=r1['v_b'])
    return r2['apogee'], r1, r2


def close2(mfix, mat, s, tw1, tw2, fin, sep, areal=10.0):
    """Two-stage closure: second stage carries mfix, the booster 0.5 mfix (0.25 mfix for solids)."""
    f = f_of(mat); key = key_of(mat); boost = (0.25 if mat == 'solid' else 0.5) * mfix
    sk = [0.0, 0.0]

    def masses(mp):
        mp2 = s * mp; mp1 = mp - mp2
        m02 = mfix + sk[1] + (1 + f) * mp2
        m01 = m02 + boost + sk[0] + (1 + f) * mp1
        return m01, mp1, m02, mp2

    def miss(mp):
        m01, mp1, m02, mp2 = masses(mp)
        return two_stage_apogee(m01, mp1, m02, mp2, rm.diameter(m01, fineness=fin), tw1, tw2, key, sep, areal)[0] - TARGET
    try:
        for _ in range(3):
            lo, hi = 4.0, 40.0
            if miss(lo) > 0:
                return None
            while miss(hi) < 0:
                hi *= 1.6
                if hi > 2e5:
                    return None
            mp = brentq(miss, lo, hi, xtol=1e-2)
            m01, mp1, m02, mp2 = masses(mp)
            _, r1, r2 = two_stage_apogee(m01, mp1, m02, mp2, rm.diameter(m01, fineness=fin), tw1, tw2, key, sep, areal)
            new = [rm.nozzle_skirt(r1['Ae'], r1['eps'], areal=areal)[0], rm.nozzle_skirt(r2['Ae'], r2['eps'], areal=areal)[0]] if areal > 0 else [0.0, 0.0]
            if max(abs(a - b) for a, b in zip(new, sk)) < 0.02:
                break
            sk = new
        m01, mp1, m02, mp2 = masses(mp)
        return dict(m0=m01, mp=mp, s=s, tw1=tw1, tw2=tw2, fin=fin, sep=sep, h_stage=r1['h_b'], v_stage=r1['v_b'])
    except ValueError:
        return None


# ------------------------------------------------------------------ parallel jobs
def job_ground1(a):
    mat, mfix, fin, tw, sep, areal = a
    r = close1(0.0, mfix, mat, areal=areal, tw=tw, fineness=fin, sep=sep)
    return (a, None if r is None else r['m0'])


def job_balloon(a):
    mat, mfix, z, fin, tw, areal, eta, cd = a
    r = close1(z * 1e3, mfix, mat, areal=areal, tw=tw, fineness=fin, eta=eta, cd_scale=cd)
    return (a, None if r is None else r['m0'])


def job_ground_mc(a):
    mat, mfix, fin, tw, sep, areal, eta, cd = a
    r = close1(0.0, mfix, mat, areal=areal, tw=tw, fineness=fin, sep=sep, eta=eta, cd_scale=cd)
    return (a, None if r is None else r['m0'])


def job_two(a):
    mat, mfix, s, tw1, tw2, fin, sep = a
    return (a, close2(mfix, mat, s, tw1, tw2, fin, sep))


# ------------------------------------------------------------------ balloon
def launched(m_rocket, z_km, gf=0.10, sigma=0.030, k_load=0.0):
    ms = (1.0 + gf) * m_rocket
    V, me, mg = rm.balloon(ms, z_km * 1e3, sigma=sigma, k_load=k_load)
    return ms + me + 1.1 * mg, V, me


def optimum(zs, tot):
    zs = np.asarray(zs, float); tot = np.asarray(tot, float)
    i = int(np.argmin(tot))
    if 0 < i < len(zs) - 1:
        c = np.polyfit(zs[i - 1:i + 2], tot[i - 1:i + 2], 2)
        zz = -c[1] / (2 * c[0]); return float(zz), float(np.polyval(c, zz))
    return float(zs[i]), float(tot[i])


def h_rho(z):
    """Local density scale height (m) at altitude z (m)."""
    e = 50.0
    return -rm.atmosphere(z)[2] / ((rm.atmosphere(z + e)[2] - rm.atmosphere(z - e)[2]) / (2 * e))


# ------------------------------------------------------------------ 1. envelope rule
def envelope_rule(system):
    fam = collections.defaultdict(list)
    for r in system:
        if r['target'] == TARGET and r['sigma'] == 0.030 and r['gf'] == 0.10 and r['m0'] is not None:
            fam[(r['mat'], r['mfix'])].append(r)
    out = []
    lift_z = np.linspace(15e3, 25e3, 21)
    dr = np.array([rm.atmosphere(z)[1] * (rm.M_AIR - rm.M_H2) / (rm.RSTAR * rm.atmosphere(z)[0]) for z in lift_z])
    cH = np.polyfit(lift_z / 1e3, np.log(dr), 1); H_lift = -1.0 / cH[0]; drho0 = math.exp(cH[1])
    for (mat, mfix), v in sorted(fam.items()):
        v = sorted(v, key=lambda r: r['z'])
        z = np.array([r['z'] for r in v]); mr = np.array([r['m0'] for r in v]); tot = np.array([r['total'] for r in v])
        zfull, _ = optimum(z[1:], tot[1:])
        sel = (z >= 10) & (z <= 35)
        p, _ = curve_fit(lambda x, a, b, L: a + b * np.exp(-x / L), z[sel], mr[sel], p0=(mr[sel].min() * 0.9, mr[sel].max(), 8.0), maxfev=20000)
        minf, B, L = p
        fit_err = float(np.max(np.abs((minf + B * np.exp(-z[sel] / L)) / mr[sel] - 1)))
        mfit = lambda x: minf + B * math.exp(-x / L)

        def gap(zk, slope=None):
            ms = 1.10 * mfit(zk)
            _, me, _ = rm.balloon(ms, zk * 1e3)
            dms = 1.10 * B / L * math.exp(-zk / L) / 1e3          # kg per m
            return me - 1.5 * h_rho(zk * 1e3) * dms
        z_rule = brentq(gap, 8.0, 38.0)
        # the two sides at the full-model optimum, from the fitted rocket curve
        _, me_star, _ = rm.balloon(1.10 * mfit(zfull), zfull * 1e3)
        dms_fd = 1.10 * B / L * math.exp(-zfull / L) / 1e3
        Hr = h_rho(zfull * 1e3)
        Leff = 1.0 / (1.0 / L + 2.0 / (3.0 * H_lift))
        mr_s = mfit(zfull)
        arg = 3 * H_lift * 1.10 * B / (2 * L * 0.030 * (36 * math.pi) ** (1 / 3) * (1.10 * mr_s / drho0) ** (2 / 3))
        z_formula = Leff * math.log(arg)
        out.append(dict(mat=mat, mfix=mfix, L=L, minf=minf, B=B, fit_err=fit_err, Hrho=Hr / 1e3, Leff=Leff, z_model=zfull,
                        z_rule=z_rule, z_formula=z_formula, me=me_star, rhs=1.5 * Hr * abs(dms_fd), ground=float(mr[0])))
    return dict(families=out, H_lift=H_lift, drho0=drho0)


# ------------------------------------------------------------------ main
def main():
    S = json.load(open(os.path.join(ROOT, 'results', 'similarity.json')))
    R = {}
    R['rule'] = envelope_rule(S['system'])

    fams = [(m, x) for m in MATS + ('solid',) for x in MFIX]
    FIN_G = (8.0, 10.0, 12.5, 15.0, 20.0, 25.0); TW_G = (1.5, 2.0, 2.5, 3.0, 4.0, 6.0); SEP_G = (1.0, 0.7, 0.55, 0.4)
    FIN_B = (8.0, 10.0, 12.5, 15.0); TW_B = (5.0, 6.0)
    with Pool() as pool:
        # 2a. single-stage ground vehicles designed for ground launch, without and with the nozzle wall
        jobs = [(m, x, fi, tw, sp, ar) for (m, x) in fams for fi in FIN_G for tw in TW_G for sp in (SEP_G if m != 'solid' else (1.0, 0.7, 0.4)) for ar in (0.0, 10.0)]
        g1 = collections.defaultdict(list)
        for a, m0 in pool.map(job_ground1, jobs, chunksize=4):
            if m0:
                g1[(a[0], a[1], a[5])].append((m0, a[2], a[3], a[4]))
        best_g1 = {k: min(v) for k, v in g1.items()}
        # 2b. two-stage ground vehicles (central nozzle wall of 10 kg/m^2)
        jobs = [(m, x, s, t1, t2, fi, sp) for (m, x) in fams for s in (0.35, 0.45, 0.55) for t1 in (2.0, 3.0, 4.0) for t2 in (1.5, 2.5, 4.0)
                for fi in (15.0, 25.0) for sp in ((0.7,) if m != 'solid' else (1.0,))]
        g2 = collections.defaultdict(list)
        for a, r in pool.map(job_two, jobs, chunksize=2):
            if r:
                g2[(a[0], a[1])].append(r)
        best_g2 = {k: min(v, key=lambda r: r['m0']) for k, v in g2.items()}
        # 2c/3. balloon vehicles designed for their own launch, without and with the nozzle wall
        jobs = [(m, x, z, fi, tw, ar, 0.92, 1.0) for (m, x) in fams for z in ZB for fi in FIN_B for tw in TW_B for ar in (0.0, 10.0)]
        b1 = collections.defaultdict(list)
        for a, m0 in pool.map(job_balloon, jobs, chunksize=4):
            if m0:
                b1[(a[0], a[1], a[2], a[5])].append((m0, a[3], a[4]))
        best_b1 = {k: min(v) for k, v in b1.items()}
        # 4. Monte Carlo surrogate: closures on a grid of efficiency, drag scale and nozzle wall
        ETA = (0.85, 0.90, 0.95); CD = (0.8, 1.0, 1.2); AR = (5.0, 20.0)
        jobs_b = []
        for (m, x) in fams:
            if m == 'solid':
                continue
            for z in ZMC:
                _, fi, tw = best_b1[(m, x, z, 10.0)]
                jobs_b += [(m, x, z, fi, tw, ar, e, c) for e in ETA for c in CD for ar in AR]
        mc_b = dict((a, m0) for a, m0 in pool.map(job_balloon, jobs_b, chunksize=4))
        jobs_g = []
        for (m, x) in fams:
            if m == 'solid':
                continue
            _, fi, tw, sp = best_g1[(m, x, 10.0)]
            jobs_g += [(m, x, fi, tw, sp, ar, e, c) for e in ETA for c in CD for ar in AR]
        mc_g = dict((a, m0) for a, m0 in pool.map(job_ground_mc, jobs_g, chunksize=4))

    R['ground1'] = {'%s/%g/%g' % k: list(v) for k, v in best_g1.items()}
    R['ground2'] = {'%s/%g' % k: v for k, v in best_g2.items()}
    R['balloon1'] = {'%s/%g/%g/%g' % k: list(v) for k, v in best_b1.items()}

    # 2/3. comparisons at the optimum
    def rockoon(m, x, areal, gf, k_load, sigma=0.030):
        tot = [launched(best_b1[(m, x, z, areal)][0], z, gf=gf, sigma=sigma, k_load=k_load)[0] for z in ZB]
        return optimum(ZB, tot)
    rows = []
    for (m, x) in fams:
        z_n, t_n = rockoon(m, x, 0.0, 0.10, 0.0)
        z_c, t_c = rockoon(m, x, 10.0, 0.20, 0.10)
        g1n = best_g1[(m, x, 0.0)][0]; g1c = best_g1[(m, x, 10.0)][0]; g2c = best_g2[(m, x)]['m0'] if (m, x) in best_g2 else float('inf')
        rows.append(dict(mat=m, mfix=x, ground1_plain=g1n, ground1=g1c, ground2=g2c, ground_best=min(g1c, g2c),
                         ground1_params=best_g1[(m, x, 10.0)][1:], ground2_params={k: best_g2[(m, x)][k] for k in ('s', 'tw1', 'tw2', 'fin')} if (m, x) in best_g2 else None,
                         z_plain=z_n, save_plain=1 - t_n / g1n, z_central=z_c, launched_central=t_c, save_central=1 - t_c / min(g1c, g2c),
                         save_central_vs1=1 - t_c / g1c))
    R['compare'] = rows

    # cost weighting with central penalties and fair comparators
    cost = {}
    for r_c in (0.3, 1.0, 3.0, 10.0):
        res = []
        for (m, x) in fams:
            c_ = []
            for z in ZB:
                mr = best_b1[(m, x, z, 10.0)][0]
                tot, V, me = launched(mr, z, gf=0.20, k_load=0.10)
                c_.append(mr + r_c * (tot - mr))
            zz, cc = optimum(ZB, c_)
            g = min(best_g1[(m, x, 10.0)][0], best_g2[(m, x)]['m0'] if (m, x) in best_g2 else 1e99)
            res.append(dict(mat=m, mfix=x, zstar=zz, dearer=bool(cc > g)))
        cost[str(r_c)] = res
    R['cost'] = cost

    # 4. Monte Carlo over the penalties (single-stage ground comparator at the same efficiency and drag)
    rng = np.random.default_rng(20261121)
    NMC = 4000
    mc = []
    from scipy.interpolate import RegularGridInterpolator
    for (m, x) in fams:
        if m == 'solid':
            continue
        def grid_b(z):
            _, fi, tw = best_b1[(m, x, z, 10.0)]
            arr = np.array([[[mc_b[(m, x, z, fi, tw, ar, e, c)] for ar in AR] for c in CD] for e in ETA])
            return RegularGridInterpolator((ETA, CD, AR), np.log(arr))
        _, fi, tw, sp = best_g1[(m, x, 10.0)]
        gG = RegularGridInterpolator((ETA, CD, AR), np.log(np.array([[[mc_g[(m, x, fi, tw, sp, ar, e, c)] for ar in AR] for c in CD] for e in ETA])))
        zg = ZMC; gB = [grid_b(z) for z in zg]
        e = rng.uniform(0.85, 0.95, NMC); c = rng.uniform(0.8, 1.2, NMC); ar = rng.uniform(5.0, 20.0, NMC)
        gf = rng.uniform(0.10, 0.30, NMC); sg = rng.uniform(0.020, 0.040, NMC); kl = rng.uniform(0.0, 0.20, NMC)
        save = np.empty(NMC); zst = np.empty(NMC)
        for j in range(NMC):
            pt = (e[j], c[j], ar[j])
            tot = [launched(math.exp(float(gB[i](pt))), zg[i], gf=gf[j], sigma=sg[j], k_load=kl[j])[0] for i in range(len(zg))]
            zz, tt = optimum(zg, tot)
            save[j] = 1 - tt / math.exp(float(gG(pt))); zst[j] = zz
        mc.append(dict(mat=m, mfix=x, p_save=float(np.mean(save > 0)), save_q=[float(q) for q in np.percentile(save, [5, 25, 50, 75, 95])],
                       z_q=[float(q) for q in np.percentile(zst, [5, 25, 50, 75, 95])]))
    R['mc'] = mc
    with open(os.path.join(ROOT, 'results', 'final.json'), 'w') as f:
        json.dump(R, f, default=float)
    print('final tests written:', len(rows), 'families compared;', len(mc), 'Monte Carlo families')


if __name__ == '__main__':
    main()
