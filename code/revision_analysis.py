"""Numbers, tables and figure data for the tests added in the revision (called from run_analysis.py).

Reads results/revision.json (revision_tests.py) and results/similarity.json (similarity.py)."""
import collections
import json
import math
import os

import numpy as np
from scipy.optimize import brentq, curve_fit

import rockoon_model as rm
from revision_tests import MATS, MFIX, ZB, ETA, CDS, AREAL, F_SOLID, BETA_LIQ, BETA_SOL

ROOT = os.path.dirname(rm.HERE)
GAS = 1.1                     # free-lift factor on the hydrogen
NOM = dict(eta=0.92, cds=1.0)
CENTRAL = dict(areal=10.0, g=0.20, env=('sigma', 0.030, 0.10))
WORST = dict(areal=20.0, g=0.30, env=('sigma', 0.030, 0.20))
NONE_ = dict(areal=0.0, g=0.10, env=('sigma', 0.030, 0.0))


def thousands(x):
    return '{:,.0f}'.format(x).replace(',', '\\,')


class Grid:
    """Rocket gross mass on the (efficiency, drag scale, nozzle wall) grid, interpolated multilinearly."""

    def __init__(self, R2):
        self.d = {}
        for k, v in R2['grid'].items():
            mat, mfix, z, eta, cds, ar = k.split('|')
            self.d[(mat, float(mfix), float(z), float(eta), float(cds), float(ar))] = v

    def __call__(self, mat, mfix, z, eta=0.92, cds=1.0, areal=0.0):
        def lerp(vals, grid, x):
            return float(np.interp(x, grid, vals))
        out = []
        for e in ETA:
            row = []
            for c in CDS:
                row.append(lerp([self.d[(mat, mfix, z, e, c, a)] for a in AREAL], AREAL, areal))
            out.append(lerp(row, CDS, cds))
        return lerp(out, ETA, eta)


def total(mr, z, g, env):
    ms = (1.0 + g) * mr
    if env[0] == 'sigma':
        V, me, mg = rm.balloon(ms, z * 1e3, sigma=env[1], k_load=env[2])
    else:
        V, me, mg = rm.balloon_law2(ms, z * 1e3, env[1], env[2], mult=env[3])
    return ms + me + GAS * mg, me, V


def optimum(z, tot):
    z = np.asarray(z); tot = np.asarray(tot)
    i = int(np.argmin(tot))
    if 0 < i < len(z) - 1:
        c = np.polyfit(z[i - 1:i + 2], tot[i - 1:i + 2], 2)
        zs = -c[1] / (2 * c[0]); return float(zs), float(np.polyval(c, zs))
    return float(z[i]), float(tot[i])


def family_case(G, mat, mfix, areal, g, env, eta=0.92, cds=1.0, cost_r=None):
    mg_ = G(mat, mfix, 0.0, eta, cds, areal)
    tots = []
    for z in ZB:
        mr = G(mat, mfix, z, eta, cds, areal)
        t, me, V = total(mr, z, g, env)
        if cost_r is not None:
            t = mr + cost_r * (t - mr)
        tots.append(t)
    zs, ts = optimum(ZB, tots)
    return dict(zstar=zs, total=ts, ground=mg_, save=1.0 - ts / mg_)


def scenario(G, sc, **kw):
    return {(mat, mfix): family_case(G, mat, mfix, sc['areal'], sc['g'], sc['env'], **kw) for mat in MATS for mfix in MFIX}


def rng_put(put, tag, res, key, scale=1.0, fmt='{:.0f}', mats=MATS):
    vals = [scale * v[key] for (m, _), v in res.items() if m in mats]
    put(tag + 'Lo', min(vals), fmt); put(tag + 'Hi', max(vals), fmt)


def analytic(S, RES, put):
    """Closed-form optimum: exact first-order condition and exponential approximation."""
    fam = collections.defaultdict(list)
    for r in S['system']:
        if r['target'] == 100e3 and r['sigma'] == 0.030 and r['gf'] == 0.10 and r['m0'] is not None:
            fam[(r['mat'], r['mfix'])].append(r)
    opt = {(o['mat'], o['mfix']): o for o in RES['optimum'] if o['target'] == 100.0}

    def Hrho(z):
        dz = 50.0
        return -2 * dz / math.log(rm.atmosphere(z + dz)[2] / rm.atmosphere(z - dz)[2])
    f = lambda z, a, b, L: a + b * np.exp(-z / L)
    rows = []
    for k, v in sorted(fam.items()):
        v = sorted(v, key=lambda r: r['z'])
        z = np.array([r['z'] for r in v]); m = np.array([r['m0'] for r in v])
        sel = (z >= 10) & (z <= 35)
        p, _ = curve_fit(f, z[sel], m[sel], p0=[m[sel][-1], m[sel][0], 10.0], maxfev=20000)
        a, B, L = p
        resid = float(np.max(np.abs(f(z[sel], *p) / m[sel] - 1)))
        zm = opt[k]['zstar']
        H = Hrho(zm * 1e3) / 1e3
        Leff = 1.0 / (1.0 / L + 2.0 / (3.0 * H))
        # exact first-order condition with the fitted rocket curve: (1+g)|m_r'| = (2/3H) m_e
        cond = lambda zz: 1.1 * B / L * math.exp(-zz / L) - 2.0 / (3.0 * Hrho(zz * 1e3) / 1e3) * rm.balloon(1.1 * f(zz, *p), zz * 1e3)[1]
        zr = brentq(cond, 8.0, 45.0)
        me = rm.balloon(1.1 * f(zm, *p), zm * 1e3)[1]
        rule = 1.5 * H * 1.1 * B / L * math.exp(-zm / L)
        rows.append(dict(mat=k[0], mfix=k[1], L=L, B=B, minf=a, fit_resid=resid, H=H, Leff=Leff, z_model=zm, z_rule=zr,
                         me=me, rule=rule, ground=opt[k]['ground']))
    put('anLLo', min(r['L'] for r in rows), '{:.1f}'); put('anLHi', max(r['L'] for r in rows), '{:.1f}')
    put('anHLo', min(r['H'] for r in rows), '{:.1f}'); put('anHHi', max(r['H'] for r in rows), '{:.1f}')
    put('anLeffLo', min(r['Leff'] for r in rows), '{:.1f}'); put('anLeffHi', max(r['Leff'] for r in rows), '{:.1f}')
    put('anFitResid', 100 * max(r['fit_resid'] for r in rows), '{:.1f}')
    put('anZdiff', max(abs(r['z_rule'] - r['z_model']) for r in rows), '{:.1f}')
    put('anRuleLo', 100 * min(r['me'] / r['rule'] - 1 for r in rows), '{:+.0f}'); put('anRuleHi', 100 * max(r['me'] / r['rule'] - 1 for r in rows), '{:+.0f}')
    # predicted sensitivities against the full model
    b80 = [r for r in rows if r['mat'] == 'Al6061T6' and r['mfix'] == 80.0][0]
    put('anLeffBase', b80['Leff'], '{:.1f}')
    put('anSigPred', b80['Leff'] * math.log(2.0), '{:.1f}')
    var = {(o['sigma'], o['gf'], o['mfix']): o['zstar'] for o in RES['optimum_var']}
    put('anSigModel', var[(0.020, 0.10, 80.0)] - var[(0.040, 0.10, 80.0)], '{:.1f}')
    put('anGonPred', b80['Leff'] / 3.0 * math.log(1.20 / 1.05), '{:.2f}')
    put('anGonModel', var[(0.030, 0.20, 80.0)] - var[(0.030, 0.05, 80.0)], '{:.2f}')
    al = [r for r in rows if r['mat'] == 'Al6061T6']
    sl_rule = np.polyfit(np.log10([r['ground'] for r in al]), [r['z_rule'] for r in al], 1)[0]
    put('anSlopeRule', sl_rule, '{:.1f}')
    # mass exponent of the benefit coefficient B e^{-z*/L}, against ground mass
    return rows


def main(put, RES, TEX):
    R2 = json.load(open(os.path.join(ROOT, 'results', 'revision.json')))
    S = json.load(open(os.path.join(ROOT, 'results', 'similarity.json')))
    G = Grid(R2)
    best = {}
    for k, v in R2['best'].items():
        t, mat, mfix, z = k.split('|')
        best[(t, mat, float(mfix), float(z))] = v
    out = {}
    # ---------------------------------------------------------- envelope law fitted to published balloons
    ef = R2['envelope_fit']
    a_fit, b_fit = ef['law2']['coef']
    put('envN', len(ef['data']), '{:d}')
    put('envA', 1000 * a_fit, '{:.1f}'); put('envB', 1000 * b_fit, '{:.2f}')
    put('envRms', 100 * ef['law2']['rms_log'], '{:.0f}'); put('envMax', 100 * ef['law2']['max_abs'], '{:.0f}')
    put('envOneRms', 100 * ef['law1']['rms_log'], '{:.0f}'); put('envOneMax', 100 * ef['law1']['max_abs'], '{:.0f}')
    put('envSigLo', 1000 * min(ef['sigma_eff']), '{:.0f}'); put('envSigHi', 1000 * max(ef['sigma_eff']), '{:.0f}')
    vols = [d['V'] for d in ef['data']]; pays = [d['mp'] for d in ef['data']]
    TEX['envVLo'] = thousands(min(vols)); TEX['envVHi'] = '%.1f' % (max(vols) / 1e6)
    put('envPayLo', min(pays), '{:.0f}'); TEX['envPayHi'] = thousands(max(pays))
    FIT = ('law2', a_fit, b_fit, 1.0)
    # the fitted law applied to the baseline balloon at 30 km
    V_, me_, mg_ = rm.balloon_law2(301.27 + 30.0, 30e3, a_fit, b_fit)
    put('envBaseMe', me_, '{:.0f}')
    out['envelope_fit'] = ef
    # ---------------------------------------------------------- designs over the wider space
    g_designs = [best[('g', mat, mfix, 0.0)] for mat in MATS for mfix in MFIX]
    put('dsGndFinLo', min(d[1] for d in g_designs), '{:.0f}'); put('dsGndFinHi', max(d[1] for d in g_designs), '{:.0f}')
    put('dsGndTwLo', min(d[2] for d in g_designs), '{:.1f}'); put('dsGndTwHi', max(d[2] for d in g_designs), '{:.1f}')
    put('dsGndSepLo', min(d[3] for d in g_designs), '{:.2f}'); put('dsGndSepHi', max(d[3] for d in g_designs), '{:.2f}')
    TEX['dsGndAtBound'] = str(sum(d[1] >= 25.0 for d in g_designs))
    same = {(o['mat'], o['mfix']): o['ground'] for o in RES['optimum'] if o['target'] == 100.0}
    red = [1 - best[('g', mat, mfix, 0.0)][0] / same[(mat, mfix)] for mat in MATS for mfix in MFIX]
    put('dsRedLo', 100 * min(red), '{:.0f}'); put('dsRedHi', 100 * max(red), '{:.0f}')
    old = RES['fair_ground']
    red2 = [1 - best[('g', mat, mfix, 0.0)][0] / old['%s/%g' % (mat, mfix)][0] for mat in MATS for mfix in MFIX]
    put('dsRedOldHi', 100 * max(red2), '{:.0f}')
    # how much a body slimmer than the bound of 25 would still save (three families)
    import itertools
    slim = []
    for mat, mfix in (('Al6061T6', 80.0), ('Al6061T6', 640.0), ('SS321', 80.0)):
        f = rm.feed_fraction(mat); row = {}
        for fin in (25.0, 30.0, 35.0):
            ms = []
            for tw, sp in itertools.product((1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0), (0.85, 0.7, 0.55, 0.4)):
                try:
                    r = rm.close_vehicle(0.0, mfix, f, tw=tw, fineness=fin, sep=sp)
                except ValueError:
                    r = None
                if r:
                    ms.append(r['m0'])
            row[fin] = min(ms)
        slim.append((1 - row[30.0] / row[25.0], 1 - row[35.0] / row[25.0]))
    put('slimThirtyLo', 100 * min(x[0] for x in slim), '{:.0f}'); put('slimThirtyHi', 100 * max(x[0] for x in slim), '{:.0f}')
    put('slimThirtyFiveLo', 100 * min(x[1] for x in slim), '{:.0f}'); put('slimThirtyFiveHi', 100 * max(x[1] for x in slim), '{:.0f}')
    # ---------------------------------------------------------- deterministic scenarios
    sc = {'none': scenario(G, NONE_), 'central': scenario(G, CENTRAL), 'worst': scenario(G, WORST),
          'fitted': scenario(G, dict(areal=10.0, g=0.20, env=FIT)), 'fittedNone': scenario(G, dict(areal=0.0, g=0.10, env=FIT))}
    for tag, res in sc.items():
        T = tag.title().replace('none', 'None')
        rng_put(put, 'sc' + T + 'Z', res, 'zstar', fmt='{:.0f}')
        rng_put(put, 'sc' + T + 'Al', res, 'save', 100, mats=('Al6061T6',))
        rng_put(put, 'sc' + T + 'Ss', res, 'save', 100, mats=('SS321',))
        rng_put(put, 'sc' + T + 'Save', res, 'save', 100)
        TEX['sc' + T + 'AlNeg'] = str(sum(v['save'] < 0 for (m, _), v in res.items() if m == 'Al6061T6'))
        out['scenario_' + tag] = {'%s|%g' % k: v for k, v in res.items()}
    # ---------------------------------------------------------- cost weighting with both vehicles designed for their own launch
    for r_c, tag in ((0.3, 'A'), (1.0, 'B'), (3.0, 'C'), (10.0, 'D')):
        res = {(mat, mfix): family_case(G, mat, mfix, CENTRAL['areal'], CENTRAL['g'], CENTRAL['env'], cost_r=r_c) for mat in MATS for mfix in MFIX}
        rng_put(put, 'fc' + tag + 'Z', res, 'zstar', fmt='{:.0f}')
        rng_put(put, 'fc' + tag + 'Save', res, 'save', 100)
        TEX['fc' + tag + 'Neg'] = str(sum(v['save'] < 0 for v in res.values()))
        out['cost_' + tag] = {'%s|%g' % k: v for k, v in res.items()}
    # ---------------------------------------------------------- Monte Carlo of the saving
    rng = np.random.default_rng(20261121)
    N = 4000
    s_fit = ef['law2']['rms_log']
    smp = dict(eta=rng.uniform(0.85, 0.95, N), cds=rng.uniform(0.8, 1.2, N), areal=rng.uniform(5.0, 20.0, N),
               g=rng.uniform(0.10, 0.30, N), law=rng.random(N) < 0.5, mult=np.exp(rng.normal(0.0, s_fit, N)),
               sig=rng.uniform(0.020, 0.040, N), k=rng.uniform(0.0, 0.20, N))
    mc = {}
    for mat in MATS:
        for mfix in MFIX:
            sv, zz = np.empty(N), np.empty(N)
            for i in range(N):
                env = ('law2', a_fit, b_fit, smp['mult'][i]) if smp['law'][i] else ('sigma', smp['sig'][i], smp['k'][i])
                r = family_case(G, mat, mfix, smp['areal'][i], smp['g'][i], env, eta=smp['eta'][i], cds=smp['cds'][i])
                sv[i], zz[i] = r['save'], r['zstar']
            mc[(mat, mfix)] = dict(p=float(np.mean(sv > 0)), save_med=float(np.median(sv)), save_lo=float(np.percentile(sv, 5)),
                                   save_hi=float(np.percentile(sv, 95)), z_med=float(np.median(zz)), z_lo=float(np.percentile(zz, 5)),
                                   z_hi=float(np.percentile(zz, 95)), save=sv.tolist()[:1000], z=zz.tolist()[:1000])
    out['mc_saving'] = {'%s|%g' % k: v for k, v in mc.items()}
    TEX['mcsN'] = thousands(N)
    for mat, tag in (('Al6061T6', 'Al'), ('SS321', 'Ss')):
        v = [mc[(mat, m)] for m in MFIX]
        put('mcs' + tag + 'PLo', 100 * min(x['p'] for x in v), '{:.1f}' if tag == 'Ss' else '{:.0f}'); put('mcs' + tag + 'PHi', 100 * max(x['p'] for x in v), '{:.0f}')
        put('mcs' + tag + 'MedLo', 100 * min(x['save_med'] for x in v), '{:.0f}'); put('mcs' + tag + 'MedHi', 100 * max(x['save_med'] for x in v), '{:.0f}')
    put('mcsZLo', min(x['z_lo'] for x in mc.values()), '{:.0f}'); put('mcsZHi', max(x['z_hi'] for x in mc.values()), '{:.0f}')
    put('mcsZWidth', max(x['z_hi'] - x['z_lo'] for x in mc.values()), '{:.1f}')
    put('mcsSaveWidth', 100 * max(x['save_hi'] - x['save_lo'] for x in mc.values()), '{:.0f}')
    b = mc[('Al6061T6', 80.0)]
    put('mcsBaseP', 100 * b['p'], '{:.0f}'); put('mcsBaseMed', 100 * b['save_med'], '{:.0f}')
    put('mcsBaseLo', 100 * b['save_lo'], '{:+.0f}'); put('mcsBaseHi', 100 * b['save_hi'], '{:.0f}')
    put('mcsBaseZLo', b['z_lo'], '{:.0f}'); put('mcsBaseZHi', b['z_hi'], '{:.0f}')
    # ---------------------------------------------------------- two-stage and solid comparators
    ts = {}
    for k, v in R2['two_stage'].items():
        tag, mfix = k.rsplit('|', 1); ts[(tag, float(mfix))] = v
    rows = []
    for mat in MATS:
        for mfix in MFIX:
            one = G(mat, mfix, 0.0, areal=CENTRAL['areal'])
            two = ts.get(('liq|' + mat, mfix), {}).get('m0', float('nan'))
            rk = sc['central'][(mat, mfix)]
            bestg = np.nanmin([one, two])
            rows.append(dict(kind='liq', mat=mat, mfix=mfix, same=same[(mat, mfix)], one=one, two=two, rockoon=rk['total'], zstar=rk['zstar'],
                             save=1 - rk['total'] / bestg, p=mc[(mat, mfix)]['p']))
    ss = {}
    for k, v in R2['solid_single'].items():
        mnp, z = k.split('|'); ss[(float(mnp), float(z))] = v
    for mnp in MFIX:
        one = ss.get((mnp, 0.0), float('nan'))
        two = ts.get(('sol', mnp), {}).get('m0', float('nan'))
        tots = [total(ss[(mnp, z)], z, CENTRAL['g'], CENTRAL['env'])[0] for z in ZB if (mnp, z) in ss]
        zs_ = [z for z in ZB if (mnp, z) in ss]
        zst, tst = optimum(zs_, tots)
        bestg = np.nanmin([one, two])
        rows.append(dict(kind='sol', mat='solid', mfix=mnp, same=float('nan'), one=one, two=two, rockoon=tst, zstar=zst, save=1 - tst / bestg, p=None))
    out['comparators'] = rows
    liq = [r for r in rows if r['kind'] == 'liq']; sol = [r for r in rows if r['kind'] == 'sol']
    tl = [r['two'] / r['one'] for r in liq if np.isfinite(r['two'])]
    put('tsLiqRatioLo', min(tl), '{:.2f}'); put('tsLiqRatioHi', max(tl), '{:.2f}')
    TEX['tsLiqLighter'] = str(sum(x < 1 for x in tl))
    tsol = [r['two'] / r['one'] for r in sol if np.isfinite(r['two'])]
    put('tsSolRatioLo', min(tsol), '{:.2f}'); put('tsSolRatioHi', max(tsol), '{:.2f}')
    TEX['tsSolLighter'] = str(sum(x < 1 for x in tsol))
    put('solZLo', min(r['zstar'] for r in sol), '{:.0f}'); put('solZHi', max(r['zstar'] for r in sol), '{:.0f}')
    put('solSaveLo', 100 * min(r['save'] for r in sol), '{:.0f}'); put('solSaveHi', 100 * max(r['save'] for r in sol), '{:.0f}')
    put('liqBestSaveAlLo', 100 * min(r['save'] for r in liq if r['mat'] == 'Al6061T6'), '{:.0f}')
    put('liqBestSaveAlHi', 100 * max(r['save'] for r in liq if r['mat'] == 'Al6061T6'), '{:.0f}')
    put('liqBestSaveSsLo', 100 * min(r['save'] for r in liq if r['mat'] == 'SS321'), '{:.0f}')
    put('liqBestSaveSsHi', 100 * max(r['save'] for r in liq if r['mat'] == 'SS321'), '{:.0f}')
    put('fSolid', F_SOLID, '{:.2f}'); put('betaLiq', BETA_LIQ, '{:.1f}'); put('betaSol', BETA_SOL, '{:.2f}')
    # table of comparators
    lines = []
    for r in rows:
        lab = {'Al6061T6': 'Al', 'SS321': 'steel', 'solid': 'solid'}[r['mat']]
        same_ = thousands(r['same']) if np.isfinite(r['same']) else '--'
        two_ = thousands(r['two']) if np.isfinite(r['two']) else '--'
        p_ = '%.0f' % (100 * r['p']) if r['p'] is not None else '--'
        lines.append('%s & %.0f & %s & %s & %s & %s & %.1f & $%+.0f$ & %s \\\\' % (lab, r['mfix'], same_, thousands(r['one']), two_, thousands(r['rockoon']),
                                                                              r['zstar'], 100 * r['save'], p_))
        if r['mat'] == 'Al6061T6' and r['mfix'] == 640.0 or r['mat'] == 'SS321' and r['mfix'] == 640.0:
            lines.append('\\midrule')
    with open(os.path.join(ROOT, 'paper', 'tab_comparators.tex'), 'w') as f_:
        f_.write('\n'.join(lines) + '\n')
    # ---------------------------------------------------------- analytic optimum
    an = analytic(S, RES, put)
    out['analytic'] = an
    with open(os.path.join(ROOT, 'paper', 'tab_analytic.tex'), 'w') as f_:
        f_.write('\n'.join('%s & %.0f & %.1f & %.1f & %.2f & %.1f & %.1f & %.0f & %.0f \\\\' % (
            'Al' if r['mat'].startswith('Al') else 'steel', r['mfix'], r['L'], r['H'], r['Leff'], r['z_model'], r['z_rule'], r['me'], r['rule'])
            for r in an) + '\n')
    # typeset negative values in running text with a minus sign
    for k in list(TEX):
        if k.startswith(('sc', 'fc', 'mcs', 'liq', 'sol', 'ds', 'slim', 'ts')) and isinstance(TEX[k], str) and TEX[k].startswith('-'):
            TEX[k] = '$-$' + TEX[k][1:]
    RES['revision'] = out
    return out
