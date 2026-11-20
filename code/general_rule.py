"""The envelope rule for any envelope law (called from run_analysis.py after revision_analysis).

The float altitude that minimises J = m_s + mu m_e, with the envelope m_e = E(V, m_s) and the buoyancy balance
V drho(z) = lift (m_s + m_e), satisfies exactly

    [1 + mu beta eps + (mu - 1) alpha eps / (1 + eps)] / Lambda = mu eps alpha / H_rho,

where eps = m_e/m_s, Lambda = m_s/|dm_s/dz|, alpha = dlnE/dlnV and beta = dlnE/dln m_s (mass-weighted over the
terms of E), and H_rho is the density scale height. mu = 1 is the launched mass. This module solves the
condition for several envelope laws on the closed rocket masses of similarity.py, compares it with direct
minimisation on the altitude grid, and tests the one-step estimate that follows from it."""
import collections
import json
import math
import os

import numpy as np
from scipy.optimize import brentq, minimize_scalar

import rockoon_model as rm
from revision_analysis import GAS, optimum
from revision_tests import MATS, MFIX, ZB

ROOT = os.path.dirname(rm.HERE)
K36 = (36.0 * math.pi) ** (1.0 / 3.0)
GON = 0.10                       # gondola fraction of Sec. 5.2
FT = 0.3048
# Latex sounding balloons, Kaymont/Totex HAB series: mass (g), burst diameter (ft)
LATEX = [(350, 13.5), (600, 19.8), (800, 23.0), (1000, 25.8), (1200, 28.3), (1500, 31.0), (2000, 34.6), (3000, 42.7),
         (4000, 50.59)]
CLUSTER_G = 3000                 # balloon used for the cluster law
LIFT_LATEX = 1.10                # a rising latex balloon carries its 10% free lift to the firing altitude
K_LOAD = 0.10                    # load-proportional envelope term of the central case
TABLE = ('area', 'load', 'fitted', 'cluster')
# Far Side balloon (1957): volume ft3, envelope lb, suspended load lb
FARSIDE = (3.75e6, 1500.0, 2400.0)


def drho(z_km):
    T, p, _ = rm.atmosphere(z_km * 1e3)
    return p * (rm.M_AIR - rm.M_H2) / (rm.RSTAR * T)


def Hrho(z_km, dz=50.0):
    z = z_km * 1e3
    return -2 * dz / math.log(rm.atmosphere(z + dz)[2] / rm.atmosphere(z - dz)[2]) / 1e3


class Law:
    """m_e = sum_i c_i V^a_i m_s^b_i, on a balloon whose volume gives `lift` times the weight it carries."""

    def __init__(self, tag, label, terms, lift=1.0):
        self.tag, self.label, self.terms, self.lift = tag, label, terms, lift

    def parts(self, V, ms):
        return np.array([c * V ** a * ms ** b for c, a, b in self.terms])

    def at(self, ms, z_km):
        """Volume and envelope mass at altitude z_km, or None if the balloon cannot lift the load there."""
        d = drho(z_km) / self.lift
        f = lambda V: V * d - ms - self.parts(V, ms).sum()
        lo, hi = 1e-6, 10.0
        for _ in range(60):
            if f(hi) > 0:
                break
            hi *= 4.0
        else:
            return None
        V = brentq(f, lo, hi, xtol=1e-10, rtol=1e-13)
        return V, float(self.parts(V, ms).sum())

    def exponents(self, V, ms):
        p = self.parts(V, ms)
        a = np.array([t[1] for t in self.terms]); b = np.array([t[2] for t in self.terms])
        return float(p @ a / p.sum()), float(p @ b / p.sum()), float(p @ a ** 2 / (p @ a)) if (p @ a) > 0 else 0.0


class Curve:
    """Closed rocket mass m_r(z) = m_inf + B exp(-z/L) and the suspended mass (1 + g) m_r."""

    def __init__(self, minf, B, L, g=GON):
        self.minf, self.B, self.L, self.g = minf, B, L, g

    def ms(self, z):
        return (1 + self.g) * (self.minf + self.B * math.exp(-z / self.L))

    def dms(self, z):
        return -(1 + self.g) * self.B / self.L * math.exp(-z / self.L)


# geometric altitudes (km) at which the lapse rate of the standard atmosphere changes; H_rho is discontinuous there
KINKS = [rm.R0_GEOPOT * h / (rm.R0_GEOPOT - h) / 1e3 for h in (11e3, 20e3, 32e3)]


def smooth(z):
    """False within 60 m of a change of lapse rate, where the optimum condition holds only between one-sided values."""
    return min(abs(z - k) for k in KINKS) > 0.06


def state(cv, law, z):
    ms = cv.ms(z); s = law.at(ms, z)
    if s is None:
        return None
    V, me = s
    al, be, a2 = law.exponents(V, ms)
    return dict(z=z, ms=ms, me=me, V=V, eps=me / ms, Lam=ms / abs(cv.dms(z)), H=Hrho(z), alpha=al, beta=be, alpha2=a2)


def ratio(cv, law, z, mu=1.0):
    """Saving side of the optimum condition over its envelope side; 1 at the optimum, above 1 below it."""
    s = state(cv, law, z)
    return (1 + mu * s['beta'] * s['eps'] + (mu - 1) * s['alpha'] * s['eps'] / (1 + s['eps'])) / s['Lam'] / (mu * s['eps'] * s['alpha'] / s['H'])


def exact(cv, law, mu=1.0, lo=5.0, hi=46.0):
    """Direct minimisation of m_s + mu m_e on the fitted curve."""
    def J(z):
        s = law.at(cv.ms(z), z)
        return np.inf if s is None else cv.ms(z) + mu * s[1]
    zz = np.linspace(lo, hi, 165); jj = [J(z) for z in zz]; i = int(np.argmin(jj))
    r = minimize_scalar(J, bounds=(zz[max(i - 1, 0)], zz[min(i + 1, len(zz) - 1)]), method='bounded', options=dict(xatol=1e-8))
    return float(r.x)


def one_step(cv, law, zt, mu=1.0):
    s = state(cv, law, zt)
    Leff = 1.0 / (s['alpha2'] / s['H'] + 1.0 / cv.L)
    return zt + Leff * math.log(ratio(cv, law, zt, mu)), Leff


def grid_optimum(zs, mr, law, g=GON):
    """Optimum of the full model: launched mass at the computed altitudes, parabola through the lowest three."""
    tot = []
    for z, m in zip(zs, mr):
        s = law.at((1 + g) * m, z)
        tot.append(np.inf if s is None else (1 + GAS * rm.M_H2 / (rm.M_AIR - rm.M_H2)) * ((1 + g) * m + s[1]))
    keep = np.isfinite(tot)
    return optimum(np.asarray(zs)[keep], np.asarray(tot)[keep])


def latex_fit():
    m = np.array([b[0] for b in LATEX]) / 1e3; d = np.array([b[1] for b in LATEX]) * FT
    V = math.pi / 6 * d ** 3
    a, lnc = np.polyfit(np.log(V), np.log(m), 1)
    r = np.exp(lnc) * V ** a / m - 1
    return dict(c=float(np.exp(lnc)), alpha=float(a), rms=float(np.sqrt(np.mean(r ** 2))), worst=float(np.max(np.abs(r))), m=m.tolist(),
                V=V.tolist(), areal=(m / (math.pi * d ** 2)).tolist())


def slope(xs, zs):
    return float(np.polyfit(np.log10(xs), zs, 1)[0])


def main(put, RES, TEX):
    S = json.load(open(os.path.join(ROOT, 'results', 'similarity.json')))
    rev = RES['revision']
    a_fit, b_fit = rev['envelope_fit']['law2']['coef']
    lf = latex_fit()
    i3 = [b[0] for b in LATEX].index(CLUSTER_G); mu1 = lf['m'][i3] / lf['V'][i3]
    laws = [Law('area', 'film', [(0.030 * K36, 2 / 3, 0.0)]),
            Law('load', 'film + load', [(0.030 * K36, 2 / 3, 0.0), (K_LOAD, 0.0, 1.0)]),
            Law('fitted', 'fitted', [(a_fit * K36, 2 / 3, 0.0), (b_fit, 1 / 3, 1.0)]),
            Law('cluster', 'cluster', [(mu1, 1.0, 0.0)], lift=LIFT_LATEX),
            Law('latex', 'latex', [(lf['c'], lf['alpha'], 0.0)], lift=LIFT_LATEX)]
    grid = collections.defaultdict(list)
    for r in S['system']:
        if r['target'] == 100e3 and r['sigma'] == 0.030 and r['gf'] == 0.10 and r['m0'] is not None and r['z'] >= 5:
            grid[(r['mat'], r['mfix'])].append((r['z'], r['m0']))
    fams = {(r['mat'], r['mfix']): r for r in rev['analytic']}
    out = dict(latex=lf, mu1=mu1, laws={})
    worst_res, step1, step2, n_kink = 0.0, 0.0, 0.0, 0
    for law in laws:
        rows = []
        for mat in MATS:
            for mfix in MFIX:
                fr = fams[(mat, mfix)]; cv = Curve(fr['minf'], fr['B'], fr['L'])
                zs, mr = zip(*sorted(grid[(mat, mfix)]))
                zx = exact(cv, law); st = state(cv, law, zx)
                zg = grid_optimum(zs, mr, law)[0]
                res = ratio(cv, law, zx)
                z1, Leff = one_step(cv, law, 20.0); z2 = one_step(cv, law, z1)[0]
                if smooth(zx):
                    worst_res = max(worst_res, abs(res - 1))
                else:
                    n_kink += 1
                if law.tag in TABLE:
                    step1 = max(step1, abs(z1 - zx)); step2 = max(step2, abs(z2 - zx))
                Lx = 1.0 / (st['alpha2'] / st['H'] + 1.0 / cv.L)
                rows.append(dict(mat=mat, mfix=mfix, ground=fr['ground'], z_exact=zx, z_grid=zg, z_one=z1, eps=st['eps'], Lam=st['Lam'],
                                 alpha=st['alpha'], beta=st['beta'], alpha2=st['alpha2'], Leff=Lx, V=st['V'], ms=st['ms']))
            # curves of the two sides for the figure (baseline family)
        cv = Curve(*[fams[('Al6061T6', 80.0)][k] for k in ('minf', 'B', 'L')])
        zz = np.arange(8.0, 36.01, 0.5); rr = []
        for z in zz:
            rr.append(ratio(cv, law, z) if law.at(cv.ms(z), z) is not None else float('nan'))
        al = [r for r in rows if r['mat'] == 'Al6061T6']; ss = [r for r in rows if r['mat'] == 'SS321']
        out['laws'][law.tag] = dict(rows=rows, line=dict(z=zz.tolist(), ratio=rr),
                                    slope_al=slope([r['ground'] for r in al], [r['z_grid'] for r in al]),
                                    slope_ss=slope([r['ground'] for r in ss], [r['z_grid'] for r in ss]),
                                    slope_al_rule=slope([r['ground'] for r in al], [r['z_exact'] for r in al]))
    L_ = out['laws']
    # ---------------------------------------------------------- numbers for the text
    put('grStepOne', step1, '{:.1f}'); put('grStepTwo', step2, '{:.2f}')
    for tag, T in (('area', 'Area'), ('load', 'Load'), ('fitted', 'Fit'), ('latex', 'Ltx'), ('cluster', 'Clu')):
        rows = L_[tag]['rows']
        for key, nm, fmt, sc in (('z_grid', 'Z', '{:.0f}', 1), ('eps', 'Eps', '{:.0f}', 100), ('Leff', 'Leff', '{:.1f}', 1), ('Lam', 'Lam', '{:.0f}', 1),
                                 ('alpha', 'Al', '{:.2f}', 1), ('beta', 'Be', '{:.2f}', 1)):
            put('gr' + T + nm + 'Lo', sc * min(r[key] for r in rows), fmt); put('gr' + T + nm + 'Hi', sc * max(r[key] for r in rows), fmt)
        put('gr' + T + 'Diff', max(abs(r['z_exact'] - r['z_grid']) for r in rows), '{:.1f}')
        put('gr' + T + 'Slope', L_[tag]['slope_al'], '{:+.1f}'); put('gr' + T + 'SlopeSs', L_[tag]['slope_ss'], '{:+.1f}')
        put('gr' + T + 'SumLo', min(r['alpha'] + r['beta'] for r in rows), '{:.2f}'); put('gr' + T + 'SumHi', max(r['alpha'] + r['beta'] for r in rows), '{:.2f}')
    put('grFitSlopeAbs', abs(L_['fitted']['slope_al']), '{:.1f}'); put('grCluSlopeAbs', abs(L_['cluster']['slope_al']), '{:.1f}')
    # the fitted law family by family: where the optimum is highest, and the shift against the area law
    fit_al = [r for r in L_['fitted']['rows'] if r['mat'] == 'Al6061T6']; area_al = [r for r in L_['area']['rows'] if r['mat'] == 'Al6061T6']
    put('grFitAlSmall', fit_al[0]['z_grid'], '{:.1f}'); put('grFitAlLarge', fit_al[-1]['z_grid'], '{:.1f}')
    put('grAreaAlSmall', area_al[0]['z_grid'], '{:.1f}'); put('grAreaAlLarge', area_al[-1]['z_grid'], '{:.1f}')
    sh = [f['z_grid'] - a['z_grid'] for f, a in zip(L_['fitted']['rows'], L_['area']['rows'])]
    put('grFitShiftLo', min(sh), '{:+.1f}'); put('grFitShiftHi', max(sh), '{:+.1f}')
    # the same comparison in the full model with each vehicle designed for its own launch (scenarios of revision_analysis)
    sn, sf = rev['scenario_none'], rev['scenario_fittedNone']
    for mat, T in (('Al6061T6', 'Al'), ('SS321', 'Ss')):
        g_ = [rev['comparators'][i]['one'] for i, c in enumerate(rev['comparators']) if c['mat'] == mat]
        zf = [sf['%s|%g' % (mat, m)]['zstar'] for m in MFIX]; zn = [sn['%s|%g' % (mat, m)]['zstar'] for m in MFIX]
        put('grOwnFit' + T + 'Slope', slope(g_, zf), '{:+.1f}'); put('grOwnArea' + T + 'Slope', slope(g_, zn), '{:+.1f}')
        put('grOwnFit' + T + 'Small', zf[0], '{:.1f}'); put('grOwnFit' + T + 'Large', zf[-1], '{:.1f}')
        put('grOwnFit' + T + 'SlopeAbs', abs(slope(g_, zf)), '{:.1f}')
    shf = [sf[k]['zstar'] - sn[k]['zstar'] for k in sn]
    put('grOwnShiftLo', min(shf), '{:+.1f}'); put('grOwnShiftHi', max(shf), '{:+.1f}')
    shc = [rev['scenario_fitted'][k]['zstar'] - rev['scenario_central'][k]['zstar'] for k in sn]
    put('grCenShiftLo', min(shc), '{:+.1f}'); put('grCenShiftHi', max(shc), '{:+.1f}')
    # latex data
    put('ltxN', len(LATEX), '{:d}'); put('ltxC', 1000 * lf['c'], '{:.1f}'); put('ltxAlpha', lf['alpha'], '{:.2f}')
    put('ltxRms', 100 * lf['rms'], '{:.0f}'); put('ltxWorst', 100 * lf['worst'], '{:.0f}')
    put('ltxVLo', min(lf['V']), '{:.0f}'); TEX['ltxVHi'] = '{:,.0f}'.format(round(max(lf['V']), -1)).replace(',', '\\,')
    put('ltxArealLo', 1000 * min(lf['areal']), '{:.1f}'); put('ltxArealHi', 1000 * max(lf['areal']), '{:.1f}')
    put('cluMu', 1000 * mu1, '{:.1f}'); put('cluV', lf['V'][i3], '{:.0f}'); put('cluKg', CLUSTER_G / 1000, '{:.0f}')
    nb = [r['V'] / lf['V'][i3] for r in L_['cluster']['rows']]
    put('cluNLo', min(nb), '{:.1f}'); put('cluNHi', max(nb), '{:.0f}')
    put('cluCeil', brentq(lambda z: drho(z) - LIFT_LATEX * mu1, 20.0, 60.0), '{:.0f}')
    vl = [r['V'] for r in L_['latex']['rows']]
    TEX['ltxNeedLo'] = '{:,.0f}'.format(round(min(vl), -1)).replace(',', '\\,'); TEX['ltxNeedHi'] = '{:,.0f}'.format(round(max(vl), -3)).replace(',', '\\,')
    put('ltxFits', sum(v <= max(lf['V']) for v in vl), '{:d}')
    # ---------------------------------------------------------- sensitivities predicted by the rule (baseline family)
    b80 = fams[('Al6061T6', 80.0)]; cv = Curve(b80['minf'], b80['B'], b80['L'])
    for tag, T in (('area', 'Area'), ('fitted', 'Fit'), ('cluster', 'Clu')):
        law = [l for l in laws if l.tag == tag][0]
        dbl = Law(tag, '', [(2 * c, a, b) for c, a, b in law.terms], lift=law.lift)
        z1, z2 = exact(cv, law), exact(cv, dbl)
        hlf = Law(tag, '', [(math.sqrt(2) * c, a, b) for c, a, b in law.terms], lift=law.lift)
        st = state(cv, hlf, 0.5 * (z1 + z2)); Leff = 1.0 / (st['alpha2'] / st['H'] + 1.0 / cv.L)
        put('grDbl' + T + 'Pred', Leff * math.log(2), '{:.1f}'); put('grDbl' + T + 'Exact', z1 - z2, '{:.1f}')
        cg = Curve(b80['minf'], b80['B'], b80['L'], g=0.05); ch = Curve(b80['minf'], b80['B'], b80['L'], g=0.20)
        st0 = state(cv, law, z1); L0 = 1.0 / (st0['alpha2'] / st0['H'] + 1.0 / cv.L)
        put('grGon' + T + 'Pred', L0 * (1 - st0['alpha'] - st0['beta']) * math.log(1.20 / 1.05), '{:+.2f}')
        put('grGon' + T + 'Exact', exact(ch, law) - exact(cg, law), '{:+.2f}')
        put('grGon' + T + 'PredAbs', abs(RES['grGon' + T + 'Pred']), '{:.2f}')
    # ---------------------------------------------------------- cost weighting: the rule against Sec. 5.4
    from revision_analysis import Grid, CENTRAL
    from scipy.optimize import curve_fit
    R2 = json.load(open(os.path.join(ROOT, 'results', 'revision.json'))); G = Grid(R2)
    gam = rm.M_H2 / (rm.M_AIR - rm.M_H2); gc = CENTRAL['g']
    mu_of = lambda r: r * (1 + GAS * gam) * (1 + gc) / (1 + r * gc + r * GAS * gam * (1 + gc))
    central = Law('central', '', [(CENTRAL['env'][1] * K36, 2 / 3, 0.0), (CENTRAL['env'][2], 0.0, 1.0)])
    f = lambda z, a, b, L: a + b * np.exp(-z / L)
    cvs = {}
    for mat in MATS:
        for mfix in MFIX:
            m = np.array([G(mat, mfix, z, areal=CENTRAL['areal']) for z in ZB])
            p, _ = curve_fit(f, np.array(ZB), m, p0=[m[-1], m[0], 10.0], maxfev=20000)
            cvs[(mat, mfix)] = Curve(p[0], p[1], p[2], g=gc)
    diffs = []; cost = {}
    for r_c, tag in ((0.3, 'A'), (1.0, 'B'), (3.0, 'C'), (10.0, 'D')):
        zr = {k: exact(cv_, central, mu=mu_of(r_c)) for k, cv_ in cvs.items()}
        zm = {tuple([k.split('|')[0], float(k.split('|')[1])]): v['zstar'] for k, v in rev['cost_' + tag].items()}
        diffs += [abs(zr[k] - zm[k]) for k in zr]; cost[tag] = dict(mu=mu_of(r_c), rule=list(zr.values()))
        worst_res = max([worst_res] + [abs(ratio(cvs[k], central, zr[k], mu=mu_of(r_c)) - 1) for k in zr if smooth(zr[k])])
        n_kink += sum(not smooth(z_) for z_ in zr.values())
        put('grMu' + tag, mu_of(r_c), '{:.2f}')
    put('grCostDiff', max(diffs), '{:.1f}'); put('grCostMed', float(np.median(diffs)), '{:.1f}')
    mu_max = (1 + GAS * gam) * (1 + gc) / (gc + GAS * gam * (1 + gc))
    put('grMuMax', mu_max, '{:.1f}')
    zinf = [exact(cv_, central, mu=mu_max) for cv_ in cvs.values()]; zone = [exact(cv_, central, mu=1.0) for cv_ in cvs.values()]
    worst_res = max([worst_res] + [abs(ratio(cv_, central, z_, mu=mu_max) - 1) for cv_, z_ in zip(cvs.values(), zinf) if smooth(z_)])
    n_kink += sum(not smooth(z_) for z_ in zinf)
    put('grInfZLo', min(zinf), '{:.0f}'); put('grInfZHi', max(zinf), '{:.0f}')
    fall = [a - b for a, b in zip(zone, zinf)]
    put('grInfFallLo', min(fall), '{:.0f}'); put('grInfFallHi', max(fall), '{:.0f}')
    # a balloon system that is cheap per kilogram: latex cluster, envelope weighted by mu < 1
    clu = [l for l in laws if l.tag == 'cluster'][0]
    cheap = {}
    for mu, tag in ((0.1, 'Tenth'), (0.03, 'Cheap')):
        cc = [Curve(fr['minf'], fr['B'], fr['L']) for fr in fams.values()]
        zc = [exact(c_, clu, mu=mu, hi=47.0) for c_ in cc]
        worst_res = max([worst_res] + [abs(ratio(c_, clu, z_, mu=mu) - 1) for c_, z_ in zip(cc, zc) if smooth(z_)])
        n_kink += sum(not smooth(z_) for z_ in zc)
        cheap[tag] = zc; put('grClu' + tag + 'Lo', min(zc), '{:.0f}'); put('grClu' + tag + 'Hi', max(zc), '{:.0f}')
    TEX['grResid'] = '%d' % math.floor(-math.log10(max(worst_res, 1e-12))); put('grKink', n_kink, '{:d}')
    V_fs, me_fs, ms_fs = FARSIDE[0] * 0.0283168466, FARSIDE[1] * 0.45359237, FARSIDE[2] * 0.45359237
    TEX['fsVol'] = '{:,.0f}'.format(round(V_fs, -2)).replace(',', '\\,'); RES['fsVol'] = V_fs
    put('fsEnv', me_fs, '{:.0f}'); TEX['fsLoad'] = '{:,.0f}'.format(ms_fs).replace(',', '\\,'); RES['fsLoad'] = ms_fs
    put('fsPhi', 100 * me_fs / ms_fs, '{:.0f}'); put('fsSigma', 1000 * me_fs / (K36 * V_fs ** (2 / 3)), '{:.0f}')
    out.update(cost=cost, mu_max=mu_max, cheap=cheap, worst_res=worst_res)
    # ---------------------------------------------------------- table of envelope laws
    lab = {'area': 'Film law, $\\sigma S$', 'load': 'Film with a load term, $\\sigma S+%.1f\\,m_s$' % K_LOAD,
           'fitted': 'Fitted law, Eq.~(\\ref{eq:envfit})', 'cluster': 'Cluster of %d\\,kg latex balloons' % (CLUSTER_G / 1000)}
    lines = []
    for tag, T in (('area', 'Area'), ('load', 'Load'), ('fitted', 'Fit'), ('cluster', 'Clu')):
        g = lambda k: TEX['gr' + T + k]
        rng = lambda k: g(k + 'Lo') if g(k + 'Lo') == g(k + 'Hi') else g(k + 'Lo') + '--' + g(k + 'Hi')
        lines.append('%s & %s & %s & %s & %s & %s & %s & $%s$ \\\\' % (lab[tag], rng('Al'), rng('Be'), rng('Z'), g('Diff'), rng('Eps'), rng('Leff'),
                                                                 g('Slope').replace('-', '{-}')))
    with open(os.path.join(ROOT, 'paper', 'tab_general.tex'), 'w') as f_:
        f_.write('\n'.join(lines) + '\n')
    RES['general'] = out
    return out
