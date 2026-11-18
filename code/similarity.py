"""Similarity study: random vehicle ensemble, optimum float altitude, and Sobol indices.
Writes results/similarity.json. Run after rockoon_model.py has built thermo_table.json."""
import json, math, os
from multiprocessing import Pool
import numpy as np
import rockoon_model as rm

ROOT = os.path.dirname(rm.HERE)
ENGINES = ['nozzle_1MPa', 'nozzle', 'nozzle_4MPa', 'solid']
CD_REF = 0.40
P0 = 101325.0


def one_vehicle(a):
    m0, lam, rhob, xi, psi, z0, eng, cds = a
    d = (4.0 * m0 / (math.pi * lam * rhob)) ** (1.0 / 3.0)
    key = ENGINES[int(eng)]
    try:
        s = rm.shapley(m0, xi * m0, z0, d=d, tw=psi, cd_scale=cds, key=key)
        g = rm.vehicle_run(m0, xi * m0, 0.0, d=d, tw=psi, cd_scale=cds, key=key)
        b = rm.vehicle_run(m0, xi * m0, z0, d=d, tw=psi, cd_scale=cds, key=key)
    except ValueError:
        return None
    A = math.pi * d * d / 4.0
    tot = s['balloon'] - s['ground']
    pi0 = cds * CD_REF * A * P0 / (m0 * rm.G0)
    return dict(m0=m0, lam=lam, rhob=rhob, xi=xi, psi=psi, z0=z0, eng=int(eng), cds=cds, d=d,
                pi0=pi0, piz=pi0 * rm.atmosphere(z0)[1] / P0, pc=rm.table()[key]['pc'],
                ground=s['ground'], balloon=s['balloon'], f_head=s['head'] / tot, f_drag=s['drag'] / tot, f_noz=s['nozzle'] / tot,
                ld_g=g['loss_d'], dv_g=g['dv_thrust'], lg_g=g['loss_g'], ld_b=b['loss_d'], dv_b=b['dv_thrust'],
                isp_g=g['isp_ign'], isp_b=b['isp_ign'], eps_g=g['eps'], eps_b=b['eps'],
                pin=rm.atmosphere(0.0)[1] * g['Ae'] / g['F0'], pic=rm.table()[key]['pc'] * A / g['F0'])


def system_point(a):
    mat, mfix, z, target, sigma, gf = a
    f = rm.feed_fraction(mat)
    r = rm.close_vehicle(z * 1e3, mfix, f, target=target)
    if r is None:
        return dict(mat=mat, mfix=mfix, z=z, target=target, sigma=sigma, gf=gf, m0=None)
    out = dict(mat=mat, mfix=mfix, z=z, target=target, sigma=sigma, gf=gf, m0=r['m0'], mp=r['mp'], d=r['d'])
    if z > 0:
        V, me, mg = rm.balloon((1.0 + gf) * r['m0'], z * 1e3, sigma=sigma)
        out.update(V=V, env=me, gas=mg, total=(1.0 + gf) * r['m0'] + me + 1.1 * mg)
    else:
        out.update(V=0.0, env=0.0, gas=0.0, total=r['m0'])
    return out


def sobol_point(a):
    z0, cd, eta, tw, fdry = a
    M0, MP = 301.27, 151.27
    m0 = MP + (M0 - MP) * fdry
    return rm.vehicle_run(m0, MP, z0, d=0.40, tw=tw, eta=eta, cd_scale=cd)['apogee']


def main():
    rng = np.random.default_rng(20261119)
    N = 800
    smp = np.column_stack([np.exp(rng.uniform(math.log(20), math.log(5000), N)), rng.uniform(8, 20, N), rng.uniform(500, 900, N),
                           rng.uniform(0.40, 0.60, N), rng.uniform(2.0, 6.0, N), rng.uniform(5e3, 40e3, N),
                           rng.integers(0, 4, N), rng.uniform(0.8, 1.2, N)])
    grid = [(mat, mfix, z, tgt) for mat in ('Al6061T6', 'SS321') for mfix in (10.0, 20.0, 40.0, 80.0, 160.0, 320.0, 640.0)
            for z in (0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5, 40) for tgt in (100e3,)]
    grid += [('Al6061T6', mfix, z, tgt) for mfix in (20.0, 80.0, 320.0) for z in (0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5, 40)
             for tgt in (60e3, 150e3)]
    grid = [g_ + (0.030, 0.10) for g_ in grid]
    ZG = (0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5, 40)
    grid += [('Al6061T6', mfix, z, 100e3, sg, gf) for mfix in (20.0, 80.0, 320.0) for z in ZG
             for sg, gf in ((0.020, 0.10), (0.040, 0.10), (0.030, 0.05), (0.030, 0.20))]
    # Sobol: Jansen estimators, 4 inputs
    NS = 1024

    def draw(n):
        return np.column_stack([rng.uniform(0.8, 1.2, n), rng.uniform(0.85, 0.95, n), rng.uniform(3.0, 5.0, n),
                                rng.normal(1.0, 0.05, n).clip(0.85, 1.15)])
    A_, B_ = draw(NS), draw(NS)
    zs = [0, 2.5, 5, 7.5, 10, 15, 20, 25, 30, 35, 40]
    jobs = []
    for z in zs:
        for M in [A_, B_] + [np.column_stack([B_[:, j] if j == i else A_[:, j] for j in range(4)]) for i in range(4)]:
            jobs += [(z * 1e3, *row) for row in M.tolist()]
    with Pool() as pool:
        ens = [r for r in pool.map(one_vehicle, smp.tolist(), chunksize=10) if r]
        sysm = pool.map(system_point, grid, chunksize=4)
        y = np.array(pool.map(sobol_point, jobs, chunksize=200))
    y = y.reshape(len(zs), 6, NS)
    sob = []
    for k, z in enumerate(zs):
        fA, fB = y[k, 0], y[k, 1]
        var = np.var(np.concatenate([fA, fB]))
        S = [float(1.0 - np.mean((fB - y[k, 2 + i]) ** 2) / (2.0 * var)) for i in range(4)]
        ST = [float(np.mean((fA - y[k, 2 + i]) ** 2) / (2.0 * var)) for i in range(4)]
        sob.append(dict(z=z, S=S, ST=ST, mean=float(np.mean(fA)) / 1e3, sd=float(np.std(fA)) / 1e3,
                        piz=CD_REF * math.pi * 0.04 * rm.atmosphere(z * 1e3)[1] / (301.27 * rm.G0)))
    with open(os.path.join(ROOT, 'results', 'similarity.json'), 'w') as f:
        json.dump(dict(ensemble=ens, system=sysm, sobol=sob), f, default=float)
    print(len(ens), 'vehicles;', len(sysm), 'closure points;', len(jobs), 'Sobol runs')


if __name__ == '__main__':
    main()
