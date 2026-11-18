"""Two further tests of the hypotheses, kept separate from similarity.py so that its results are unchanged.

1. H3 on several vehicles: total-effect Sobol index of the drag scale factor against ignition altitude
   for four vehicles with different Pi_D(0), to see whether the curves collapse when plotted against
   Pi_D(z0) rather than altitude.
2. H1 at a fixed float altitude: the 800-vehicle ensemble of similarity.py, redrawn with the same
   generator but every vehicle ignited at 30 km (and at 15 km), so that the float altitude cannot
   spread the shares.

Writes results/extra.json.
"""
import json
import math
import os
from multiprocessing import Pool

import numpy as np

import rockoon_model as rm
from similarity import one_vehicle, CD_REF, P0

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
XI = 151.27 / 301.27
VEHICLES = [('30 kg', 30.0, 10.0), ('301 kg', 301.27, 10.0), ('3000 kg', 3000.0, 10.0), ('301 kg slender', 301.27, 20.0)]
ZS = [0, 2.5, 5, 7.5, 10, 12.5, 15, 20, 25, 30, 40]
NS = 512


def sob_point(a):
    m_nom, fin, z0, cd, eta, tw, fdry = a
    mp = XI * m_nom
    m0 = mp + (m_nom - mp) * fdry
    try:
        return rm.vehicle_run(m0, mp, z0, d=rm.diameter(m_nom, fineness=fin), tw=tw, eta=eta, cd_scale=cd)['apogee']
    except ValueError:
        return float('nan')


def main():
    rng = np.random.default_rng(20261120)

    def draw(n):
        return np.column_stack([rng.uniform(0.8, 1.2, n), rng.uniform(0.85, 0.95, n), rng.uniform(3.0, 5.0, n),
                                rng.normal(1.0, 0.05, n).clip(0.85, 1.15)])
    A_, B_ = draw(NS), draw(NS)
    mats = [A_, B_] + [np.column_stack([B_[:, j] if j == i else A_[:, j] for j in range(4)]) for i in range(4)]
    jobs = [(m, fin, z * 1e3, *row) for _, m, fin in VEHICLES for z in ZS for M in mats for row in M.tolist()]

    # the ensemble of similarity.py, same generator and draws, with the float altitude replaced
    rng_e = np.random.default_rng(20261119)
    N = 800
    smp = np.column_stack([np.exp(rng_e.uniform(math.log(20), math.log(5000), N)), rng_e.uniform(8, 20, N), rng_e.uniform(500, 900, N),
                           rng_e.uniform(0.40, 0.60, N), rng_e.uniform(2.0, 6.0, N), rng_e.uniform(5e3, 40e3, N),
                           rng_e.integers(0, 4, N), rng_e.uniform(0.8, 1.2, N)])
    fixed = {}
    with Pool() as pool:
        y = np.array(pool.map(sob_point, jobs, chunksize=200))
        for zf in (15e3, 30e3):
            s2 = smp.copy(); s2[:, 5] = zf
            fixed[str(int(zf / 1e3))] = [r for r in pool.map(one_vehicle, s2.tolist(), chunksize=10) if r]
    y = y.reshape(len(VEHICLES), len(ZS), 6, NS)
    sob = []
    for v, (name, m, fin) in enumerate(VEHICLES):
        d = rm.diameter(m, fineness=fin)
        pi0 = CD_REF * math.pi * d * d / 4 * P0 / (m * rm.G0)
        for k, z in enumerate(ZS):
            fA, fB, fC = y[v, k, 0], y[v, k, 1], y[v, k, 2]
            ok = np.isfinite(fA) & np.isfinite(fB) & np.isfinite(fC)
            var = np.var(np.concatenate([fA[ok], fB[ok]]))
            st = float(np.mean((fA[ok] - fC[ok]) ** 2) / (2.0 * var))
            bs = np.random.default_rng(7).integers(0, ok.sum(), (300, ok.sum()))
            a_, c_ = fA[ok], fC[ok]
            ci = np.percentile([np.mean((a_[b] - c_[b]) ** 2) / (2.0 * var) for b in bs], [2.5, 97.5])
            sob.append(dict(vehicle=name, m0=m, fineness=fin, d=d, z=z, pi0=pi0, piz=pi0 * rm.atmosphere(z * 1e3)[1] / P0,
                            ST_cd=st, ci=[float(ci[0]), float(ci[1])], n=int(ok.sum())))
    os.makedirs(os.path.join(ROOT, 'results'), exist_ok=True)
    with open(os.path.join(ROOT, 'results', 'extra.json'), 'w') as f:
        json.dump(dict(sobol_vehicles=sob, fixed_float=fixed), f, default=float)
    print('extra tests written:', len(sob), 'Sobol points;', {k: len(v) for k, v in fixed.items()}, 'fixed-float vehicles')


if __name__ == '__main__':
    main()
