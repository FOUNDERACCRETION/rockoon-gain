"""Numbers and tables for the tests with published dry-mass data, for drag that grows with fineness, and for the
Aerobee 150 check (called from run_analysis.py after revision_analysis).

Reads results/realdata.json (realdata_tests.py) and results/aerobee_check.json (aerobee_check.py)."""
import json
import os

import numpy as np

import rockoon_model as rm
from revision_analysis import CENTRAL, optimum, thousands, total
from revision_tests import MATS, MFIX, ZB

ROOT = os.path.dirname(rm.HERE)


def rockoon(mr):
    """Optimum float altitude and launched mass of a rockoon with the central penalties; mr maps altitude to rocket mass."""
    zs = [z for z in ZB if z in mr]
    return optimum(zs, [total(mr[z], z, CENTRAL['g'], CENTRAL['env'])[0] for z in zs])


def rng(put, tag, vals, fmt='{:.0f}', scale=1.0):
    put(tag + 'Lo', scale * min(vals), fmt); put(tag + 'Hi', scale * max(vals), fmt)


def main(put, RES, TEX):
    D = json.load(open(os.path.join(ROOT, 'results', 'realdata.json')))
    rev = RES['revision']; out = {}
    comp = {(r['mat'], r['mfix']): r for r in rev['comparators']}
    # ---------------------------------------------------------- 1. solid motors with flown inert fractions
    allsave = []
    for name, T in (('BlackBrantV', 'Bb'), ('SuperLoki', 'Sl')):
        d = D['solid_' + name]
        single = {tuple(float(x) for x in k.split('|')): v for k, v in d['single'].items()}
        rows = []
        for mnp in MFIX:
            one = single[(mnp, 0.0)]; two = d['two'].get('%g' % mnp, {}).get('m0', float('nan'))
            zst, tst = rockoon({z: single[(mnp, z)] for z in ZB if (mnp, z) in single})
            rows.append(dict(mfix=mnp, one=one, two=two, rockoon=tst, zstar=zst, save=1 - tst / np.nanmin([one, two])))
        put('rdSol' + T + 'F', D['f_solid_real'][name], '{:.2f}')
        rng(put, 'rdSol' + T + 'Save', [r['save'] for r in rows], scale=100); rng(put, 'rdSol' + T + 'Z', [r['zstar'] for r in rows])
        rng(put, 'rdSol' + T + 'Two', [r['two'] / r['one'] for r in rows], '{:.2f}')
        TEX['rdSol' + T + 'TwoLighter'] = str(sum(r['two'] < r['one'] for r in rows))
        allsave += [r['save'] for r in rows]; out['solid_' + name] = rows
    rng(put, 'rdSolSave', allsave, scale=100)
    rng(put, 'rdSolF', list(D['f_solid_real'].values()), '{:.2f}')
    # ---------------------------------------------------------- 2. liquid family with the feed fraction of a flown stage
    d = D['liquid_aerobee']
    closed = {tuple(float(x) for x in k.split('|')): v for k, v in d['closed'].items()}
    best = {tuple(float(x) for x in k.split('|')): v for k, v in d['best'].items()}
    rows = []
    for mfix in MFIX:
        one = closed[(mfix, 0.0)]; two = d['two'].get('%g' % mfix, {}).get('m0', float('nan'))
        zst, tst = rockoon({z: closed[(mfix, z)] for z in ZB if closed.get((mfix, z)) is not None})
        rows.append(dict(mfix=mfix, one=one, two=two, rockoon=tst, zstar=zst, save=1 - tst / np.nanmin([one, two]), fin=best[(mfix, 0.0)][1]))
    put('rdAbF', D['f_aerobee'], '{:.2f}')
    put('rdAbRatioAl', rm.feed_fraction('Al6061T6') / D['f_aerobee'], '{:.1f}'); put('rdAbRatioSs', rm.feed_fraction('SS321') / D['f_aerobee'], '{:.1f}')
    rng(put, 'rdAbSave', [r['save'] for r in rows], scale=100); rng(put, 'rdAbZ', [r['zstar'] for r in rows])
    rng(put, 'rdAbTwo', [r['two'] / r['one'] for r in rows], '{:.2f}')
    TEX['rdAbNeg'] = str(sum(r['save'] < 0 for r in rows))
    out['liquid_aerobee'] = rows
    # ---------------------------------------------------------- 3. skin friction that grows with fineness
    d = D['fin_drag']
    key3 = lambda k: (k.split('|')[0], float(k.split('|')[1]), float(k.split('|')[2]))
    closed = {key3(k): v for k, v in d['closed'].items()}; best = {key3(k): v for k, v in d['best'].items()}
    rows = []
    for mat in MATS:
        for mfix in MFIX:
            one = closed[(mat, mfix, 0.0)]; two = d['two'].get('%s|%g' % (mat, mfix), {}).get('m0', float('nan'))
            zst, tst = rockoon({z: closed[(mat, mfix, z)] for z in ZB if closed.get((mat, mfix, z)) is not None})
            c = comp[(mat, mfix)]
            rows.append(dict(mat=mat, mfix=mfix, one=one, two=two, rockoon=tst, zstar=zst, save=1 - tst / np.nanmin([one, two]),
                             fin=best[(mat, mfix, 0.0)][1], heavier=one / c['one'] - 1, rk_heavier=tst / c['rockoon'] - 1, save_generic=c['save'],
                             two_fin=d['two'].get('%s|%g' % (mat, mfix), {}).get('fin', float('nan'))))
    rng(put, 'rdFdFin', [r['fin'] for r in rows]); TEX['rdFdAtBound'] = str(sum(r['fin'] >= 25.0 for r in rows))
    rng(put, 'rdFdHeavy', [r['heavier'] for r in rows], scale=100); rng(put, 'rdFdRkHeavy', [r['rk_heavier'] for r in rows], '{:.1f}', scale=100)
    rng(put, 'rdFdZ', [r['zstar'] for r in rows])
    for mat, T in (('Al6061T6', 'Al'), ('SS321', 'Ss')):
        rng(put, 'rdFdSave' + T, [r['save'] for r in rows if r['mat'] == mat], scale=100)
        rng(put, 'rdFdGain' + T, [r['save'] - r['save_generic'] for r in rows if r['mat'] == mat], '{:.0f}', scale=100)
    TEX['rdFdAlNeg'] = str(sum(r['save'] < 0 for r in rows if r['mat'] == 'Al6061T6'))
    out['fin_drag'] = rows
    for k in list(TEX):
        if k.startswith('rd') and isinstance(TEX[k], str) and TEX[k].startswith('-'):
            TEX[k] = '$-$' + TEX[k][1:]
    # ---------------------------------------------------------- 4. Aerobee 150
    A = json.load(open(os.path.join(ROOT, 'results', 'aerobee_check.json')))
    th = A['thermo']
    put('abTc', th['Tc'], '{:.0f}'); put('abCstar', th['cstar'], '{:.0f}'); put('abIspSlIdeal', th['isp_sl_ideal'], '{:.1f}')
    put('abIspVacIdeal', th['isp_vac_ideal'], '{:.1f}'); put('abIspSlPub', A['isp_sl_published'], '{:.1f}')
    put('abEtaSl', A['eta_sl'], '{:.2f}'); put('abEtaVac', A['eta_vac'], '{:.2f}'); put('abPc', 324.0 * 6894.757 / 1e6, '{:.1f}')
    put('abFeed', A['feed_fraction'], '{:.2f}'); put('abTankRatio', A['tank_pressure_ratio'], '{:.2f}')
    nom = [r for r in A['rows'] if r['nominal']]; flt = [r for r in A['rows'] if not r['nominal']]
    err = lambda r: 100 * (r['model_peak'] / r['peak'] - 1); verr = lambda r: 100 * (r['model_vb'] / r['vb'] - 1)
    put('abNomA', err(nom[0]), '{:+.1f}'); put('abNomB', err(nom[1]), '{:+.1f}')
    put('abNomVbA', verr(nom[0]), '{:+.1f}'); put('abNomVbB', verr(nom[1]), '{:+.1f}')
    put('abNomMax', max(abs(err(r)) for r in nom), '{:.0f}')
    rng(put, 'abFlt', [err(r) for r in flt]); put('abFltMean', np.mean([err(r) for r in flt]), '{:.0f}')
    put('abFltN', len(flt), '{:d}')
    lines = []
    for r in A['rows']:
        vb = ('%.0f & %.0f & $%+.1f$' % (r['vb'], r['model_vb'], verr(r))) if r['vb'] is not None else '-- & %.0f & --' % r['model_vb']
        lines.append('%s & %.0f & %s & %.0f & %.0f & $%+.1f$ & %s \\\\' % (r['label'], r['payload'], thousands(r['z0']) if r['z0'] else '0', r['peak'] / 1e3,
                                                                    r['model_peak'] / 1e3, err(r), vb))
        if r is nom[-1]:
            lines.append('\\midrule')
    with open(os.path.join(ROOT, 'paper', 'tab_aerobee.tex'), 'w') as f_:
        f_.write('\n'.join(lines) + '\n')
    out['aerobee'] = A
    RES['realdata'] = out
    return out
