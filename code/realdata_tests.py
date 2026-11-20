"""Tests with published dry-mass data and with drag that grows with fineness. Run after revision_tests.py.

1. Solid motors with the inert mass per kilogram of propellant of flown sounding-rocket motors
   (Black Brant V and Super Loki) in place of the generic 0.15.
2. A pressure-fed liquid family with the tank-and-pressurant mass per kilogram of propellant of a flown
   pressure-fed stage (Aerobee 150A sustainer) in place of the ideal-sphere estimates.
3. Ground- and balloon-launched vehicles designed again with skin friction that grows with fineness.

Writes results/realdata.json."""
import itertools
import json
import os
from multiprocessing import Pool

import rockoon_model as rm
from revision_tests import (MATS, MFIX, ZB, FIN_G, TW_G, SEP_G, FIN_B, TW_B, TARGET, BETA_LIQ, BETA_SOL, two_stage_family)

ROOT = os.path.dirname(rm.HERE)
LB = 0.45359237
# inert mass (case, liner, nozzle, igniter) per kg of propellant
SOLID_REAL = {
    'BlackBrantV': (2803.0 - 2243.0) / 2243.0,          # loaded motor 2803 lb, propellant 2243 lb
    'SuperLoki': (28.32 - 16.87 - 6.12) / 16.87,        # launch 28.32 kg, propellant 16.87 kg, dart 6.12 kg (NASA CR-61238)
}
# Aerobee 150A sustainer: tank assembly with its helium tank 74.156 kg, helium 5.15 lb, propellant 758.2 + 303.3 lb
AEROBEE_PROP = (758.2 + 303.3) * LB
F_AEROBEE = (74.156 + 5.15 * LB) / AEROBEE_PROP
AREAL = 10.0


def _close(args):
    z, mfix, f, kw, areal = args
    try:
        r = rm.close_vehicle_skirt(z * 1e3, mfix, f, target=TARGET, areal=areal, **kw) if areal > 0 else \
            rm.close_vehicle(z * 1e3, mfix, f, target=TARGET, **kw)
    except ValueError:
        r = None
    return None if r is None else r['m0']


def solid_family(pool, f):
    """Single-stage solid vehicles with inert fraction f, over the design space of revision_tests.solid_singles."""
    jobs, keys = [], []
    for mnp in MFIX:
        for fin, tw, sp in itertools.product((10.0, 15.0, 20.0, 25.0), (3.0, 5.0, 8.0, 12.0), (1.0, 0.7, 0.4)):
            jobs.append((0.0, mnp, f, dict(tw=tw, fineness=fin, sep=sp, key='solid'), AREAL)); keys.append((mnp, 0.0))
        for z in ZB:
            for fin, tw in itertools.product((8.0, 10.0, 15.0), (3.0, 5.0, 8.0)):
                jobs.append((z, mnp, f, dict(tw=tw, fineness=fin, key='solid'), AREAL)); keys.append((mnp, z))
    best = {}
    for k, m in zip(keys, pool.map(_close, jobs, chunksize=10)):
        if m is not None and (k not in best or m < best[k]):
            best[k] = m
    return best


def liquid_family(pool, tag_f, extra=None):
    """Pressure-fed liquid families: best designs without the nozzle wall, then closed again with it
    (as in revision_tests.best_designs and mass_grid). tag_f maps a label to a feed fraction."""
    extra = extra or {}
    jobs, keys = [], []
    for tag, f in tag_f.items():
        for mfix in MFIX:
            for fin, tw, sp in itertools.product(FIN_G, TW_G, SEP_G):
                jobs.append((0.0, mfix, f, dict(tw=tw, fineness=fin, sep=sp, **extra), 0.0)); keys.append((tag, mfix, 0.0, fin, tw, sp))
            for z in ZB:
                for fin, tw in itertools.product(FIN_B, TW_B):
                    jobs.append((z, mfix, f, dict(tw=tw, fineness=fin, **extra), 0.0)); keys.append((tag, mfix, z, fin, tw, 1.0))
    best = {}
    for k, m in zip(keys, pool.map(_close, jobs, chunksize=20)):
        if m is None:
            continue
        kk = k[:3]
        if kk not in best or m < best[kk][0]:
            best[kk] = (m, k[3], k[4], k[5])
    jobs, keys = [], []
    for (tag, mfix, z), b in best.items():
        kw = dict(tw=b[2], fineness=b[1], **extra)
        if z == 0:
            kw['sep'] = b[3]
        jobs.append((z, mfix, tag_f[tag], kw, AREAL)); keys.append((tag, mfix, z))
    closed = dict(zip(keys, pool.map(_close, jobs, chunksize=4)))
    return best, closed


def main():
    out = {'f_solid_real': SOLID_REAL, 'f_aerobee': F_AEROBEE}
    with Pool() as pool:
        for name, f in SOLID_REAL.items():
            single = solid_family(pool, f)
            two = two_stage_family(pool, [('sol', mnp, f, BETA_SOL, 'solid', (4.0, 8.0, 12.0), (2.0, 4.0, 8.0), (20.0, 25.0), AREAL) for mnp in MFIX])
            out['solid_' + name] = dict(single={'%g|%g' % k: v for k, v in single.items()}, two={'%g' % k[1]: v for k, v in two.items()})
            print('solid', name, 'done', flush=True)
        best, closed = liquid_family(pool, {'aerobee': F_AEROBEE})
        two = two_stage_family(pool, [('liq', mfix, F_AEROBEE, BETA_LIQ, 'nozzle', (2.0, 3.0, 4.0), (1.5, 2.5, 4.0), (20.0, 25.0), AREAL) for mfix in MFIX])
        out['liquid_aerobee'] = dict(best={'%g|%g' % k[1:]: list(v) for k, v in best.items()}, closed={'%g|%g' % k[1:]: v for k, v in closed.items()},
                                     two={'%g' % k[1]: v for k, v in two.items()})
        print('liquid, flown feed fraction, done', flush=True)
        tag_f = {mat: rm.feed_fraction(mat) for mat in MATS}
        best, closed = liquid_family(pool, tag_f, extra=dict(fin_drag=True))
        two = two_stage_family(pool, [('liq|' + mat, mfix, tag_f[mat], BETA_LIQ, 'nozzle', (2.0, 3.0, 4.0), (1.5, 2.5, 4.0), (15.0, 20.0, 25.0), AREAL, True)
                                      for mat in MATS for mfix in MFIX])
        out['fin_drag'] = dict(best={'%s|%g|%g' % k: list(v) for k, v in best.items()}, closed={'%s|%g|%g' % k: v for k, v in closed.items()},
                               two={'%s|%g' % (k[0].split('|')[1], k[1]): v for k, v in two.items()})
        print('fineness-dependent drag done', flush=True)
    with open(os.path.join(ROOT, 'results', 'realdata.json'), 'w') as f_:
        json.dump(out, f_, default=float)
    print('real-data tests written')


if __name__ == '__main__':
    main()
