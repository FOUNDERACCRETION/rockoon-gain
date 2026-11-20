"""Run every calculation in the paper. Writes results/results.json and paper/results.tex."""
import json, math, os, sys
from multiprocessing import Pool
import numpy as np
from scipy.optimize import brentq, minimize_scalar
import rockoon_model as rm

ROOT = os.path.dirname(rm.HERE)
RES = {}
TEX = {}


def put(key, val, fmt='{:.1f}'):
    RES[key] = float(val) if isinstance(val, (int, float, np.floating)) else val
    TEX[key] = fmt.format(val) if not isinstance(val, str) else val


def thousands(x):
    return '{:,.0f}'.format(x).replace(',', '\\,')


M0, MP = 301.27, 151.27
CDREF = 0.40
LB, FT, MI = 0.45359237, 0.3048, 1609.344


def mc_one(args):
    cd, eta, tw, fdry, zf = args
    m0 = MP + (M0 - MP) * fdry
    b = rm.vehicle_run(m0, MP, zf, d=0.40, tw=tw, eta=eta, cd_scale=cd)['apogee']
    g = rm.vehicle_run(m0, MP, 0.0, d=0.40, tw=tw, eta=eta, cd_scale=cd)['apogee']
    return b, g


def mc_shap(args):
    cd, eta, tw, fdry, zf = args
    m0 = MP + (M0 - MP) * fdry
    s = rm.shapley(m0, MP, zf, d=0.40, tw=tw, eta=eta, cd_scale=cd)
    tot = s['balloon'] - s['ground']
    return s['head'] / tot, s['drag'] / tot, s['nozzle'] / tot


def main():
    t = rm.table()
    nzt = t['nozzle']
    # ---- atmosphere
    T, p, rho = rm.atmosphere(30e3)
    put('atmT', T); put('atmP', p, '{:.0f}'); put('atmRho', rho, '{:.4f}')
    Tp, pp, rp = rm.lapse_extrapolation(30e3)
    put('lapT', Tp, '{:.2f}'); put('lapP', pp, '{:.0f}'); put('lapRho', rp, '{:.4f}')
    put('lapZero', 288.15 / 6.49, '{:.0f}')
    put('propRhoOwn', 2654.98 * rm.M_AIR / (rm.RSTAR * Tp), '{:.3f}')
    put('rhoH2float', p * rm.M_H2 / (rm.RSTAR * T), '{:.5f}')
    # ---- balloon
    put('gasH2', rm.M_H2 / (rm.M_AIR - rm.M_H2), '{:.4f}')
    put('gasHe', rm.M_HE / (rm.M_AIR - rm.M_HE), '{:.3f}')
    Vbare = M0 * rm.RSTAR * T / (p * (rm.M_AIR - rm.M_H2))
    put('balBareV', Vbare, '{:.0f}'); RES['balBareV_tex'] = None
    TEX['balBareV'] = thousands(round(Vbare, -2))
    put('balBareH2', M0 * RES['gasH2'], '{:.1f}')
    V, me, mg = rm.balloon(M0 + 30.0, 30e3)
    TEX['balV'] = thousands(round(V, -2)); RES['balV'] = V
    put('balEnv', me, '{:.0f}'); put('balH2', mg, '{:.1f}')
    put('balArea', (36 * math.pi) ** (1 / 3) * V ** (2 / 3), '{:.0f}')
    for s, tag in ((0.020, 'Lo'), (0.040, 'Hi')):
        v_, _, _ = rm.balloon(M0 + 30.0, 30e3, sigma=s)
        TEX['balV' + tag] = thousands(round(v_, -2)); RES['balV' + tag] = v_
    Vh, meh, mgh = rm.balloon(M0 + 30.0, 30e3, M_gas=rm.M_HE)
    TEX['balVHe'] = thousands(round(Vh, -2)); RES['balVHe'] = Vh
    put('balHe', mgh, '{:.1f}'); put('balHeVolPct', 100 * (Vh / V - 1), '{:.0f}'); put('balHeMassRatio', mgh / mg, '{:.2f}')
    put('balLaunchH2', 1.1 * mg, '{:.0f}')
    put('propVolRatio', V / 2870.0, '{:.1f}'); put('propH2Ratio', 258.42 / mg, '{:.1f}')
    put('filmAreal', 20e-6 * 920 * 1000, '{:.1f}')
    # ---- thermochemistry
    put('Tc', nzt['Tc'], '{:.0f}'); put('cstar', nzt['cstar'], '{:.0f}')
    sw = t['of_sweep']
    b = max(sw, key=lambda s: s['shift']); put('ofShiftPeak', b['of'], '{:.2f}'); put('ispShiftPeak', b['shift'])
    b = max(sw, key=lambda s: s['frozen']); put('ofFrozenPeak', b['of'], '{:.2f}'); put('ispFrozenPeak', b['frozen'])
    b = [s for s in sw if abs(s['of'] - 2.4) < 1e-9][0]
    put('ispShift', b['shift']); put('ispFrozen', b['frozen']); put('ispDel', 0.92 * b['shift'])
    put('ispDelLo', 0.85 * b['shift'], '{:.0f}'); put('ispDelHi', 0.95 * b['shift'], '{:.0f}')
    fr = t['fuel_rich']
    put('frTc', fr['Tc'], '{:.0f}'); put('frShift', fr['shift']); put('frFrozen', fr['frozen']); put('frDel', 0.92 * fr['shift'], '{:.0f}')
    put('frDeficitPct', 100 * (1 - fr['shift'] / 316.7), '{:.0f}')
    bm = t['benchmark_68atm']
    put('bmFrozen', bm['2.1']['frozen']); put('bmShift', bm['2.6']['shift']); put('bmTc', bm['2.6']['Tc'], '{:.0f}'); put('bmCstar', bm['2.6']['cstar'], '{:.0f}')
    # ---- validation: Deacon ground rounds (NACA RM L8H26, Table I)
    rounds = [(2, 232.65, 3307, 20008), (3, 215.30, 3571, 19800), (4, 213.40, 3747, 20168), (5, 197.45, 4094, 19912),
              (6, 190.50, 4180, 19512), (7, 196.51, 4122, 19904), (8, 196.08, 3817, 19160), (9, 197.49, 3667, 18976)]
    mp_d, d_d, tb_d = 103.0 * LB, 6.75 * 0.0254, 2.9
    put('deaMp', mp_d, '{:.1f}'); put('deaD', d_d, '{:.3f}')
    put('deaIsp', 19680 * 4.4482216 / (mp_d * rm.G0), '{:.0f}'); put('deaImp', 19680 * 4.4482216 / 1e3, '{:.1f}')

    def deacon(wlb, imp, cd=1.0, z0=0.0, elev=70.0, coast=False):
        m0 = wlb * LB
        I = imp * 4.4482216
        return rm.fly(m0, mp_d, d_d, mp_d / tb_d, I / (mp_d * rm.G0), 0.0, z_kin=z0, cd_scale=cd, elev_deg=elev, coast=coast)

    def rms(cd, sel):
        return math.sqrt(np.mean([(deacon(w, I, cd)['v_max'] / (v * FT) - 1) ** 2 for (_, w, v, I) in sel]))
    warm = rounds[:6]
    val = []
    for (n, w, v, I) in rounds:
        r = deacon(w, I)
        val.append(dict(round=n, m0=w * LB, v_meas=v * FT, v_model=r['v_max'], err=100 * (r['v_max'] / (v * FT) - 1),
                        imp=I * 4.4482216 / 1e3, v_nodrag=deacon(w, I, 0.0)['v_max']))
    RES['deacon_rounds'] = val
    put('deaErrMean', np.mean([x['err'] for x in val[:6]]), '{:+.1f}')
    put('deaErrMaxAbs', max(abs(x['err']) for x in val[:6]), '{:.1f}')
    put('deaErrCold', np.mean([x['err'] for x in val[6:]]), '{:+.1f}')
    put('deaRmsUnit', 100 * rms(1.0, warm), '{:.1f}')
    fit = minimize_scalar(lambda c: rms(c, warm), bounds=(0.3, 3.0), method='bounded')
    put('deaCdFit', fit.x, '{:.2f}'); put('deaRmsFit', 100 * fit.fun, '{:.1f}')
    put('deaDragShare', 100 * np.mean([(x['v_nodrag'] - x['v_model']) / x['v_nodrag'] for x in val[:6]]), '{:.0f}')
    # rockoon SUI flight 3: fired 38,000 ft, summit about 200,000 ft (Ludwig 2004, p. 8)
    z_fire, z_top = 38000 * FT, 200000 * FT
    put('suiFire', z_fire / 1e3); put('suiTop', z_top / 1e3, '{:.0f}')
    sui = {}
    for cd in (0.8, 1.0, 1.2, RES['deaCdFit']):
        f = lambda w: deacon(w, 19680, cd, z0=z_fire, elev=90.0, coast=True)['apogee'] - z_top
        try:
            w = brentq(f, 120.0, 400.0)
        except ValueError:
            w = float('nan')
        sui[cd] = w * LB
    put('suiMassLo', sui[0.8], '{:.0f}'); put('suiMassMid', sui[1.0], '{:.0f}'); put('suiMassHi', sui[1.2], '{:.0f}')
    put('suiMassFit', sui[RES['deaCdFit']], '{:.0f}')
    aps = [deacon(w, 19680, 1.0, z0=z_fire, elev=90.0, coast=True)['apogee'] / 1e3 for w in (190.50, 232.65)]
    put('suiApLight', aps[0], '{:.0f}'); put('suiApHeavy', aps[1], '{:.0f}')
    put('nacaMassLo', 190.50 * LB, '{:.0f}'); put('nacaMassHi', 232.65 * LB, '{:.0f}')
    # Deacon rockoon firings of 1952-53 (NRL Report 4276, data sheets): launch altitude and zenith, statute miles
    SUI = [('SUI-3', 7.2, 37.9), ('SUI-4', 6.6, 36.9), ('SUI-7', 8.1, 39.8), ('SUI-13', 10.8, 54.8), ('SUI-15', 8.5, 40.7), ('SUI-23', 12.2, 61.6)]
    NRLF = [('NRL-1', 15.7, 49.2), ('NRL-4', 14.2, 50.2), ('NRL-5', 14.0, 41.7)]
    def rk_ap(m, zmi, cd=1.0):
        return deacon(m / LB, 19680, cd, z0=zmi * MI, elev=90.0, coast=True)['apogee'] / 1e3
    def rk_rms(m, D, cd=1.0):
        return math.sqrt(np.mean([(rk_ap(m, a_, cd) / (b_ * MI / 1e3) - 1) ** 2 for _, a_, b_ in D]))
    rk = {}
    for tag, D, cd in (('Sui', SUI, 1.0), ('Nrl', NRLF, 1.0), ('Early', SUI[:3], 1.0), ('SuiCd', SUI, 2.0)):
        r_ = minimize_scalar(lambda m: rk_rms(m, D, cd), bounds=(70.0, 160.0), method='bounded')
        rk[tag] = r_.x
        put('rk' + tag + 'Mass', r_.x, '{:.0f}'); put('rk' + tag + 'Rms', 100 * r_.fun, '{:.1f}')
    rows = []
    for grp, D, m_ in (('SUI', SUI, rk['Sui']), ('NRL', NRLF, rk['Nrl'])):
        for nm, a_, b_ in D:
            mod = rk_ap(m_, a_)
            rows.append(dict(name=nm, grp=grp, fire=a_ * MI / 1e3, obs=b_ * MI / 1e3, model=mod, err=100 * (mod / (b_ * MI / 1e3) - 1)))
    RES['rockoons'] = rows; RES['rockoon_mass'] = rk
    put('rkN', len(rows), '{:d}'); put('rkFireLo', min(r['fire'] for r in rows), '{:.0f}'); put('rkFireHi', max(r['fire'] for r in rows), '{:.0f}')
    put('rkObsLo', min(r['obs'] for r in rows), '{:.0f}'); put('rkObsHi', max(r['obs'] for r in rows), '{:.0f}')
    put('rkWorst', min(r['err'] for r in rows), '{:+.0f}'); put('rkBest', max(r['err'] for r in rows), '{:+.0f}')
    put('rkAllRms', 100 * math.sqrt(np.mean([(r['err'] / 100) ** 2 for r in rows])), '{:.0f}')
    with open(os.path.join(ROOT, 'paper', 'tab_rockoons.tex'), 'w') as f_:
        f_.write('\n'.join('%s & %.1f & %.1f & %.1f & $%+.0f$ \\\\' % (r['name'], r['fire'], r['obs'], r['model'], r['err']) for r in rows) + '\n')
    # the same flights with an independently reported gross mass (no fit): Encyclopedia Astronautica lists the
    # Deacon Rockoon at 205 lb gross with a 37 lb payload
    m_ax = 205.0 * LB
    ax_err = {grp: [100 * (rk_ap(m_ax, a_) / (b_ * MI / 1e3) - 1) for _, a_, b_ in D] for grp, D in (('Sui', SUI), ('Nrl', NRLF))}
    put('rkAxMass', m_ax, '{:.0f}')
    for grp in ('Sui', 'Nrl'):
        put('rkAx' + grp + 'Mean', np.mean(ax_err[grp]), '{:.0f}'); put('rkAx' + grp + 'Rms', math.sqrt(np.mean(np.square(ax_err[grp]))), '{:.0f}')
    allax = ax_err['Sui'] + ax_err['Nrl']
    put('rkAxRms', math.sqrt(np.mean(np.square(allax))), '{:.0f}'); put('rkAxLo', min(allax), '{:+.0f}'); put('rkAxHi', max(allax), '{:+.0f}')
    RES['rockoon_ax_err'] = ax_err
    fx = minimize_scalar(lambda c: rk_rms(m_ax, SUI, c), bounds=(0.5, 6.0), method='bounded')
    put('rkAxCdFit', fx.x, '{:.1f}'); put('rkAxCdRms', 100 * fx.fun, '{:.0f}')
    gnd = [deacon(w, 19680, 1.0, z0=0.0, elev=90.0, coast=True)['apogee'] / 1e3 for w in (190.50, 232.65)]
    put('deaGroundLight', gnd[0], '{:.0f}'); put('deaGroundHeavy', gnd[1], '{:.0f}')
    # ---- baseline
    cases = {'bal': rm.vehicle_run(M0, MP, 30e3), 'gnd': rm.vehicle_run(M0, MP, 0.0),
             'fall': rm.vehicle_run(M0, MP, 29e3, v0=-math.sqrt(2 * 9.72 * 1000.0))}
    put('fallSpeed', math.sqrt(2 * 9.72 * 1000.0), '{:.0f}')
    RES['cases'] = {k: {a: float(b) if not isinstance(b, str) else b for a, b in v.items()} for k, v in cases.items()}
    for k, r in cases.items():
        put(k + 'Ap', r['apogee'] / 1e3, '{:.1f}' if k == 'gnd' else '{:.0f}'); put(k + 'Tb', r['tb']); put(k + 'Hb', (r['apogee'] * 0 + r['h_b'] + {'bal': 30e3, 'gnd': 0.0, 'fall': 29e3}[k]) / 1e3)
        put(k + 'Vb', r['v_b'], '{:.0f}'); put(k + 'Dv', r['dv_thrust'], '{:.0f}'); put(k + 'Lg', r['loss_g'], '{:.0f}')
        put(k + 'Ld', r['loss_d'], '{:.0f}'); put(k + 'Q', r['q_max'] / 1e3, '{:.1f}' if r['q_max'] < 1e4 else '{:.0f}')
        put(k + 'Eps', r['eps'], '{:.1f}' if r['eps'] > 10 else '{:.2f}'); put(k + 'Isp', r['isp_ign']); put(k + 'IspVac', r['isp_vac'])
        put(k + 'F', r['F0'] / 1e3); put(k + 'Mdot', r['mdot'], '{:.2f}'); put(k + 'Pe', r['pe'] / 1e3, '{:.1f}'); put(k + 'Tap', r['t_apogee'], '{:.0f}')
    put('gainRatio', cases['bal']['apogee'] / cases['gnd']['apogee'], '{:.1f}')
    put('gainKm', (cases['bal']['apogee'] - cases['gnd']['apogee']) / 1e3, '{:.1f}')
    put('fallCost', (cases['bal']['apogee'] - cases['fall']['apogee']) / 1e3, '{:.0f}')
    put('vmin', math.sqrt(2 * rm.GM * (1 / (rm.R_EARTH + 30e3) - 1 / (rm.R_EARTH + 100e3))), '{:.0f}')
    put('peakAccel', cases['bal']['F0'] / (M0 - MP) / rm.G0 * (1 + rm.atmosphere(30e3)[1] * cases['bal']['Ae'] / cases['bal']['F0']), '{:.1f}')
    # one-at-a-time sensitivity of the baseline
    def ap(zf=30e3, **kw):
        a = dict(d=0.40); a.update(kw)
        return rm.vehicle_run(M0, MP, zf, **a)['apogee'] / 1e3
    put('sTwLo', ap(tw=2.0), '{:.0f}'); put('sTwHi', ap(tw=8.0), '{:.0f}')
    put('sCdLo', ap(cd_scale=0.5), '{:.0f}'); put('sCdHi', ap(cd_scale=2.0), '{:.0f}')
    put('sEtaLo', ap(eta=0.85), '{:.0f}'); put('sEtaHi', ap(eta=0.95), '{:.0f}')
    put('sGndCdLo', ap(0.0, cd_scale=0.8), '{:.0f}'); put('sGndCdHi', ap(0.0, cd_scale=1.2), '{:.0f}')
    put('sGndD30', ap(0.0, d=0.30), '{:.0f}'); put('sGndD35', ap(0.0, d=0.35), '{:.0f}')
    # ---- decomposition
    def dec(m0, mp, zf, **kw):
        s = rm.shapley(m0, mp, zf, **kw)
        tot = s['balloon'] - s['ground']
        return dict(m0=m0, zf=zf, ground=s['ground'] / 1e3, balloon=s['balloon'] / 1e3, gain=tot / 1e3,
                    head=s['head'] / 1e3, drag=s['drag'] / 1e3, nozzle=s['nozzle'] / 1e3,
                    f_head=s['head'] / tot, f_drag=s['drag'] / tot, f_nozzle=s['nozzle'] / tot,
                    runs={''.join(map(str, k)): v / 1e3 for k, v in s['runs'].items()}, **{k: v for k, v in kw.items()})
    base = dec(M0, MP, 30e3)
    RES['dec_base'] = base
    for k in ('head', 'drag', 'nozzle'):
        put('dec' + k.title(), base[k], '{:.1f}'); put('decF' + k.title(), 100 * base['f_' + k], '{:.0f}')
    put('decGain', base['gain'], '{:.1f}')
    for k, v in base['runs'].items():
        put('run' + k.replace('0', 'G').replace('1', 'B'), v, '{:.0f}')
    # spread of each contribution over the six switching orders
    import itertools
    rr = {tuple(int(c) for c in k): v for k, v in base['runs'].items()}
    for i_, nm in enumerate(('Head', 'Drag', 'Nozzle')):
        vals = []
        for order in itertools.permutations(range(3)):
            st = [0, 0, 0]
            for j_ in order:
                if j_ == i_:
                    t_ = list(st); t_[i_] = 1
                    vals.append(100 * (rr[tuple(t_)] - rr[tuple(st)]) / base['gain'])
                    break
                st[j_] = 1
        put('ord' + nm + 'Lo', min(vals), '{:.0f}'); put('ord' + nm + 'Hi', max(vals), '{:.0f}')
    # attribution measured in burnout velocity and in burnout energy, from the same eight runs
    runs8 = {(k, dr, n): rm.vehicle_run(M0, MP, 30e3, switches=(k, dr, n)) for k in (0, 1) for dr in (0, 1) for n in (0, 1)}
    def shap(metric):
        w = {0: 1 / 3, 1: 1 / 6, 2: 1 / 3}; phi = [0.0, 0.0, 0.0]
        for i in range(3):
            for s_ in runs8:
                if s_[i] == 0:
                    t_ = list(s_); t_[i] = 1
                    phi[i] += w[sum(s_)] * (metric(runs8[tuple(t_)]) - metric(runs8[s_]))
        return phi
    pv = shap(lambda r: r['v_b']); tv = sum(pv)
    put('dvbGain', tv, '{:.0f}')
    for i, nm in enumerate(('Head', 'Drag', 'Nozzle')):
        put('dvb' + nm, pv[i], '{:.0f}'); put('dvbF' + nm, 100 * pv[i] / tv, '{:.0f}')
    rb, rg = runs8[(1, 1, 1)], runs8[(0, 0, 0)]
    put('budDvThr', rb['dv_thrust'] - rg['dv_thrust'], '{:.0f}'); put('budLd', rg['loss_d'] - rb['loss_d'], '{:.0f}'); put('budLg', rb['loss_g'] - rg['loss_g'], '{:.0f}')
    # coast approximation: apogee ~ z0 + h_b + v_b^2/(2 g); the derivative v_b/g converts velocity to height
    put('coastFactor', rb['v_b'] / rm.G0, '{:.0f}'); put('coastFactorGnd', rg['v_b'] / rm.G0, '{:.0f}')
    put('coastApprox', (30e3 + rb['h_b'] + rb['v_b'] ** 2 / (2 * rm.G0 * (1 + (30e3 + rb['h_b']) / 6.371e6) ** -2)) / 1e3, '{:.0f}')
    fr_ = MP / M0
    put('propFrac', fr_, '{:.3f}')
    RES['dec_alt'] = [dec(M0, MP, z * 1e3) for z in (5, 10, 15, 20, 25, 30, 35, 40)]
    RES['dec_mass'] = [dec(m, fr_ * m, 30e3) for m in (25, 50, 100, 200, 301.27, 500, 1000, 2000)]
    RES['dec_diam'] = [dec(M0, MP, 30e3, d=d) for d in (0.25, 0.30, 0.35, 0.40, 0.50, 0.60)]
    for row, tag in ((RES['dec_alt'][1], 'AltTen'), (RES['dec_alt'][3], 'AltTwenty'), (RES['dec_mass'][1], 'MassFifty'), (RES['dec_mass'][-1], 'MassBig')):
        for k in ('head', 'drag', 'nozzle'):
            put('dec' + tag + k.title(), 100 * row['f_' + k], '{:.0f}')
        put('dec' + tag + 'Gnd', row['ground'], '{:.0f}'); put('dec' + tag + 'Bal', row['balloon'], '{:.0f}')
        put('dec' + tag + 'Ratio', row['balloon'] / row['ground'], '{:.1f}')
    for row in RES['dec_mass'] + RES['dec_diam']:
        dd = row.get('d', rm.diameter(row['m0']))
        row['beta'] = row['m0'] / (math.pi * dd * dd / 4)
    put('betaBase', M0 / (math.pi * 0.04), '{:.0f}')
    put('betaLo', RES['dec_mass'][0]['beta'], '{:.0f}'); put('betaHi', RES['dec_mass'][-1]['beta'], '{:.0f}')
    # ---- robustness: each vehicle designed for its own launch condition, and chamber pressure
    tw_max = 15.0 * (1 - fr_)
    put('twMax', tw_max, '{:.1f}')
    best = {}
    for tag, z in (('Gnd', 0.0), ('Bal', 30e3)):
        cand = []
        for fin in (10.0, 12.5, 15.0, 17.5, 20.0):
            for tw in np.arange(1.5, tw_max + 1e-9, 0.25):
                try:
                    cand.append((rm.vehicle_run(M0, MP, z, d=rm.diameter(M0, fineness=fin), tw=float(tw))['apogee'] / 1e3, fin, float(tw)))
                except ValueError:
                    pass
        best[tag] = max(cand)
        put('opt' + tag + 'Ap', best[tag][0], '{:.0f}'); put('opt' + tag + 'Fin', best[tag][1], '{:.1f}'); put('opt' + tag + 'Tw', best[tag][2], '{:.2f}')
        put('opt' + tag + 'D', rm.diameter(M0, fineness=best[tag][1]), '{:.2f}')
    put('optRatio', best['Bal'][0] / best['Gnd'][0], '{:.1f}'); put('optGain', best['Bal'][0] - best['Gnd'][0], '{:.0f}')
    slim = dec(M0, MP, 30e3, d=rm.diameter(M0, fineness=20.0))
    RES['dec_slim'] = slim
    put('slimD', rm.diameter(M0, fineness=20.0), '{:.2f}')
    for k in ('head', 'drag', 'nozzle'):
        put('slimF' + k.title(), 100 * slim['f_' + k], '{:.0f}')
    put('slimGnd', slim['ground'], '{:.0f}'); put('slimBal', slim['balloon'], '{:.0f}'); put('slimRatio', slim['balloon'] / slim['ground'], '{:.1f}')
    RES['dec_pc'] = []
    for key, tag, pcv in (('nozzle_1MPa', 'One', 1.0), ('nozzle', 'Two', 2.0), ('nozzle_4MPa', 'Four', 4.0)):
        row = dec(M0, MP, 30e3, key=key); row['pc'] = pcv
        g = rm.vehicle_run(M0, MP, 0.0, key=key); b = rm.vehicle_run(M0, MP, 30e3, key=key)
        row.update(eps_g=g['eps'], eps_b=b['eps'], isp_g=g['isp_ign'], isp_b=b['isp_ign'])
        RES['dec_pc'].append(row)
        for k in ('head', 'drag', 'nozzle'):
            put('pc' + tag + 'F' + k.title(), 100 * row['f_' + k], '{:.0f}')
        put('pc' + tag + 'Gnd', row['ground'], '{:.0f}'); put('pc' + tag + 'Bal', row['balloon'], '{:.0f}')
        put('pc' + tag + 'EpsG', g['eps'], '{:.1f}'); put('pc' + tag + 'EpsB', b['eps'], '{:.0f}')
        put('pc' + tag + 'IspG', g['isp_ign'], '{:.0f}'); put('pc' + tag + 'IspB', b['isp_ign'], '{:.0f}')
        put('pc' + tag + 'Ratio', row['balloon'] / row['ground'], '{:.1f}')
    # ---- scaling map and balloon cost of altitude
    zs = np.arange(0.0, 42.5, 2.5); ms = [50, 100, 200, 301.27, 500]
    grid = [[rm.vehicle_run(m, fr_ * m, z * 1e3)['apogee'] / 1e3 for z in zs] for m in ms]
    RES['map'] = dict(z=zs.tolist(), m=ms, apogee=grid)
    alt = []
    for z in np.arange(15.0, 42.5, 2.5):
        V_, me_, mg_ = rm.balloon(M0 + 30.0, z * 1e3)
        alt.append(dict(z=float(z), V=V_, env=me_, gas=mg_, apogee=rm.vehicle_run(M0, MP, z * 1e3)['apogee'] / 1e3))
    RES['alt_cost'] = alt
    a20 = [x for x in alt if x['z'] == 20.0][0]; a30 = [x for x in alt if x['z'] == 30.0][0]; a40 = [x for x in alt if x['z'] == 40.0][0]
    put('altApTwenty', a20['apogee'], '{:.0f}'); put('altApForty', a40['apogee'], '{:.0f}')
    TEX['altVTwenty'] = thousands(round(a20['V'], -2)); TEX['altVForty'] = thousands(round(a40['V'], -3))
    RES['altVTwenty'] = a20['V']; RES['altVForty'] = a40['V']
    put('altGainLow', (a30['apogee'] - a20['apogee']) / 10, '{:.1f}'); put('altGainHigh', (a40['apogee'] - a30['apogee']) / 10, '{:.1f}')
    put('altVolRatioLow', a30['V'] / a20['V'], '{:.1f}'); put('altVolRatioHigh', a40['V'] / a30['V'], '{:.1f}')
    put('gndFirstTen', (grid[3][4] - grid[3][0]) / 10, '{:.1f}')
    put('mapApTen', grid[3][4], '{:.0f}')
    # ---- tanks and closure
    ft_hp = rm.tank_fraction(20.7e6, 2 / 3, 8000.0, 205e6)
    put('ftHP', ft_hp, '{:.2f}'); put('tankHP', ft_hp * MP, '{:.0f}'); put('tankHPsf', rm.tank_fraction(20.7e6, 2 / 3, 8000.0, 205e6, sf=1.1) * MP, '{:.0f}')
    put('propVolL', rm.prop_volume(2 / 3) * MP * 1e3, '{:.0f}')
    put('mrCap', (1 + ft_hp) / ft_hp, '{:.2f}')
    put('dvCap', 0.92 * fr['shift'] * rm.G0 * math.log((1 + ft_hp) / ft_hp), '{:.0f}')
    put('dvCapOpt', RES['ispDel'] * rm.G0 * math.log((1 + ft_hp) / ft_hp), '{:.0f}')
    for pe_, tag in ((22e6, 'Lo'), (25e6, 'Mid'), (28e6, 'Hi')):
        mhe, mb, Vb = rm.pressurant(20.7e6, 2 / 3, 34.5e6, pe_, 8000.0, 205e6)
        put('heBottle' + tag, mb * MP / 1e3, '{:.1f}'); put('heVol' + tag, Vb * MP * 1e3, '{:.0f}'); put('heMass' + tag, mhe * MP, '{:.0f}')
    put('ftSS', rm.tank_fraction(3e6, 2.4, 8000.0, 205e6), '{:.3f}'); put('ftAl', rm.tank_fraction(3e6, 2.4, 2700.0, 276e6), '{:.3f}')
    mhe, mb, _ = rm.pressurant(3e6, 2.4, 30e6, 3.6e6, 2700.0, 276e6)
    put('heKg', mhe, '{:.4f}'); put('heBottleKg', mb, '{:.3f}')
    fAl, fSS = rm.feed_fraction('Al6061T6'), rm.feed_fraction('SS321')
    put('feedAl', fAl, '{:.3f}'); put('feedSS', fSS, '{:.3f}')
    clos = {}
    for name, f in (('Al', fAl), ('SS', fSS)):
        for mfix in (40.0, 80.0, 120.0):
            rows = []
            for z in (0, 5, 10, 15, 20, 25, 30, 35, 40):
                r = rm.close_vehicle(z * 1e3, mfix, f)
                if r is None:
                    rows.append(dict(z=z, m0=None)); continue
                row = dict(z=z, m0=r['m0'], mp=r['mp'], dry=r['m_dry'], F=r['F0'], tb=r['tb'], d=r['d'], eps=r['eps'])
                if z >= 15:
                    V_, me_, mg_ = rm.balloon(r['m0'] + 30.0, z * 1e3)
                    row.update(V=V_, env=me_, gas=mg_)
                rows.append(row)
            clos[f'{name}_{int(mfix)}'] = rows
    RES['closure'] = clos
    c = {r['z']: r for r in clos['Al_80']}
    put('cGndM', c[0]['m0'], '{:.0f}'); put('cGndMp', c[0]['mp'], '{:.0f}')
    put('cBalM', c[30]['m0'], '{:.0f}'); put('cBalMp', c[30]['mp'], '{:.1f}'); put('cBalDry', c[30]['dry'], '{:.1f}')
    put('cBalF', c[30]['F'] / 1e3); put('cBalTb', c[30]['tb']); put('cBalD', c[30]['d'], '{:.2f}')
    TEX['cBalV'] = thousands(round(c[30]['V'], -2)); RES['cBalV'] = c[30]['V']
    put('cBalEnv', c[30]['env'], '{:.0f}'); put('cBalGas', c[30]['gas'], '{:.1f}')
    put('cMassRatio', c[0]['m0'] / c[30]['m0'], '{:.1f}'); put('cPropRatio', c[0]['mp'] / c[30]['mp'], '{:.1f}')
    put('cSysBal', c[30]['m0'] + 30 + c[30]['env'] + 1.1 * c[30]['gas'], '{:.0f}')
    c20 = c[20]; put('cTwentyM', c20['m0'], '{:.0f}'); TEX['cTwentyV'] = thousands(round(c20['V'], -2)); RES['cTwentyV'] = c20['V']
    s = {r['z']: r for r in clos['SS_80']}
    put('cBalMSS', s[30]['m0'], '{:.0f}'); put('cGndMSS', s[0]['m0'], '{:.0f}')
    r83 = brentq(lambda mp: rm.vehicle_run(150 + mp, mp, 30e3, d=0.40)['apogee'] - 100e3, 30, 150)
    put('mpFor100', r83, '{:.0f}')
    d274 = brentq(lambda md: rm.vehicle_run(md + MP, MP, 30e3, d=0.40)['apogee'] - 100e3, 150, 500)
    put('dryFor100', d274, '{:.0f}')
    # ---- Monte Carlo
    rng = np.random.default_rng(20261118)
    N = 10000
    smp = np.column_stack([rng.uniform(0.8, 1.2, N), rng.uniform(0.85, 0.95, N), rng.uniform(3.0, 5.0, N),
                           rng.normal(1.0, 0.05, N).clip(0.85, 1.15), rng.uniform(28e3, 32e3, N)])
    with Pool() as pool:
        out = np.array(pool.map(mc_one, smp.tolist(), chunksize=100))
        sh = np.array(pool.map(mc_shap, smp[:1000].tolist(), chunksize=25))
    bal, gnd_ = out[:, 0] / 1e3, out[:, 1] / 1e3
    RES['mc'] = dict(N=N, bal=bal.tolist(), gnd=gnd_.tolist(), shares=sh.tolist())
    TEX['mcN'] = thousands(N); TEX['mcNshap'] = thousands(1000)
    for arr, tag, f in ((bal, 'Bal', '{:.0f}'), (gnd_, 'Gnd', '{:.0f}'), (bal / gnd_, 'Ratio', '{:.1f}')):
        for q, qt in ((5, 'Lo'), (50, 'Med'), (95, 'Hi')):
            put('mc' + tag + qt, np.percentile(arr, q), f)
    bsr = np.random.default_rng(5)
    med = [np.median(bsr.choice(bal, len(bal))) for _ in range(500)]; medg = [np.median(bsr.choice(gnd_, len(gnd_))) for _ in range(500)]
    put('mcBalMedSe', np.std(med), '{:.1f}'); put('mcGndMedSe', np.std(medg), '{:.1f}')
    put('mcBalCv', 100 * np.std(bal) / np.mean(bal), '{:.0f}'); put('mcGndCv', 100 * np.std(gnd_) / np.mean(gnd_), '{:.0f}')
    put('mcPbelow', 100 * np.mean(bal < 100.0), '{:.1f}')
    for i, tag in enumerate(('Head', 'Drag', 'Nozzle')):
        for q, qt in ((5, 'Lo'), (50, 'Med'), (95, 'Hi')):
            put('mcF' + tag + qt, 100 * np.percentile(sh[:, i], q), '{:.0f}')
    put('mcDragLargest', 100 * np.mean((sh[:, 1] > sh[:, 0]) & (sh[:, 1] > sh[:, 2])), '{:.0f}')
    # rank correlations of the balloon and ground apogee with the inputs
    from scipy.stats import spearmanr
    names = ['Cd', 'Eta', 'Tw', 'Dry', 'Zf']
    for i, n in enumerate(names):
        put('rhoBal' + n, spearmanr(smp[:, i], bal)[0], '{:+.2f}'); put('rhoGnd' + n, spearmanr(smp[:, i], gnd_)[0], '{:+.2f}')
    # ---- audit arithmetic (Appendix)
    put('audDvBack', 317.6 * 9.81 * math.log(301.27 / 150.0), '{:.0f}')
    put('audDvPct', 100 * (1 - 317.6 * 9.81 * math.log(301.27 / 150.0) / 2854.7), '{:.0f}')
    put('audLnA', math.log(2.5), '{:.3f}'); put('audLnB', math.log(3.5), '{:.3f}'); put('audDiv', 226.215 / 3.7419, '{:.2f}')
    put('audPratio', 2654.98 / pp, '{:.1f}'); put('audRhoRatio', 1.03 / RES['propRhoOwn'], '{:.1f}')
    put('audRhoLow', rho / rp, '{:.1f}'); put('audRhoTot', 1.03 / rho, '{:.0f}')
    put('audH2ratio', 0.0899 / RES['rhoH2float'], '{:.0f}'); put('audMolar', 301.27 * 9.81 / 1.03 / 22.414 * 7.4, '{:.0f}')
    put('audVolRatio', Vbare / 2870.0, '{:.1f}'); put('audGasRatio', 258.42 / RES['balBareH2'], '{:.1f}')
    v = math.sqrt(rm.GM / (rm.R_EARTH + 400e3))
    put('audEratio', (v * v / 2 + rm.GM * (1 / rm.R_EARTH - 1 / (rm.R_EARTH + 400e3))) / (rm.GM * (1 / rm.R_EARTH - 1 / (rm.R_EARTH + 100e3))), '{:.0f}')
    put('h2oRatio', 103.0 / 17.5, '{:.0f}')

    # ---- validation: Super Loki Dart (NASA CR-61238)
    F_sl, FTm = 4450 * 4.4482216, 0.3048
    def loki(z0, el, cds=1.0):
        return rm.fly2d(28.32, 16.87, 85.47e-4, F_sl, 2.0, rm.cd_superloki, el, z0=z0, m2=6.12, area2=9.84e-4, cd_fn2=rm.cd_dart, cd_scale=cds)
    a = loki(0.0, 80.0); b = loki(1219.0, 80.0)
    put('lokiVb', a['v_b'], '{:.0f}'); put('lokiVbRep', 5893 * FTm, '{:.0f}'); put('lokiVbErr', 100 * (a['v_b'] / (5893 * FTm) - 1), '{:+.1f}')
    put('lokiAp', a['apogee'] / 1e3, '{:.1f}'); put('lokiApErr', 100 * (a['apogee'] / 113.4e3 - 1), '{:+.1f}'); put('lokiTap', a['t_apogee'], '{:.0f}')
    put('lokiRange', a['range'] / 1e3, '{:.0f}')
    put('lokiWs', b['apogee'] / 1e3, '{:.0f}'); put('lokiWsErr', 100 * (b['apogee'] / 125e3 - 1), '{:+.0f}')
    put('lokiImp', 8900 * 4.4482216 / 1e3, '{:.1f}'); put('lokiF', F_sl / 1e3, '{:.1f}')
    A_l = 85.47e-4
    put('lokiPiZero', 0.5 * A_l * 101325.0 / (28.32 * rm.G0), '{:.1f}')
    put('lokiPiNom', 0.5 * math.pi * 0.1016 ** 2 / 4 * 101325.0 / (28.32 * rm.G0), '{:.1f}')
    # ---- robustness to the shape of the drag curve: NASA CR-61238 booster-plus-dart table in place of the generic curve
    alt = dec(M0, MP, 30e3, cd_fn=rm.cd_superloki)
    alt.pop('cd_fn', None)
    RES['dec_altcurve'] = alt
    put('acGnd', alt['ground'], '{:.0f}'); put('acBal', alt['balloon'], '{:.0f}'); put('acRatio', alt['balloon'] / alt['ground'], '{:.1f}')
    for k in ('head', 'drag', 'nozzle'):
        put('acF' + k.title(), 100 * alt['f_' + k], '{:.0f}')
    put('acBalChange', 100 * (alt['balloon'] / base['balloon'] - 1), '{:+.1f}'); put('acGndChange', 100 * (alt['ground'] / base['ground'] - 1), '{:+.0f}')
    # ---- non-vertical ignition (planar model)
    tilt = []
    nb = rm.Nozzle(30e3, 4 * M0 * rm.G0, 0.40); ng = rm.Nozzle(0.0, 4 * M0 * rm.G0, 0.40)
    for a_ in (0, 5, 10, 15, 20, 30):
        bf = rm.fly_tilted(M0, MP, 0.40, nb.mdot, nb.isp_vac, nb.Ae, 30e3, a_, 'fixed')
        bt = rm.fly_tilted(M0, MP, 0.40, nb.mdot, nb.isp_vac, nb.Ae, 30e3, a_, 'turn')
        gt = rm.fly_tilted(M0, MP, 0.40, ng.mdot, ng.isp_vac, ng.Ae, 0.0, a_, 'turn')
        tilt.append(dict(tilt=a_, bal_fixed=bf['apogee'] / 1e3, bal_fixed_range=bf['range'] / 1e3, bal_turn=bt['apogee'] / 1e3,
                         bal_turn_range=bt['range'] / 1e3, gnd_turn=gt['apogee'] / 1e3, gnd_range=gt['range'] / 1e3))
    RES['tilt'] = tilt
    put('tiltCheck', tilt[0]['bal_fixed'], '{:.0f}')
    for row, tag in ((tilt[1], 'Five'), (tilt[2], 'Ten'), (tilt[4], 'Twenty')):
        put('tiltFixed' + tag, row['bal_fixed'], '{:.0f}'); put('tiltTurn' + tag, row['bal_turn'], '{:.0f}'); put('tiltGnd' + tag, row['gnd_turn'], '{:.0f}')
        put('tiltFixedLoss' + tag, 100 * (1 - row['bal_fixed'] / tilt[0]['bal_fixed']), '{:.0f}'); put('tiltTurnLoss' + tag, 100 * (1 - row['bal_turn'] / tilt[0]['bal_turn']), '{:.0f}')
        put('tiltGndLoss' + tag, 100 * (1 - row['gnd_turn'] / tilt[0]['gnd_turn']), '{:.0f}')
        put('tiltRatio' + tag, row['bal_fixed'] / row['gnd_turn'], '{:.1f}'); put('tiltRange' + tag, row['bal_fixed_range'], '{:.0f}')
    with open(os.path.join(ROOT, 'paper', 'tab_tilt.tex'), 'w') as f:
        f.write('\n'.join('%d & %.0f & %.0f & %.0f & %.0f & %.0f & %.1f \\\\' % (r['tilt'], r['bal_fixed'], r['bal_fixed_range'], r['bal_turn'], r['bal_turn_range'], r['gnd_turn'], r['bal_fixed'] / r['gnd_turn']) for r in tilt) + '\n')
    # ---- drag curve computed by component build-up for a defined geometry
    bu = rm.DragBuildUp()
    Ms = [0.3, 0.8, 1.05, 2.0, 3.0, 5.0]
    T0_, p0_, r0_ = rm.atmosphere(0.0); a0_ = math.sqrt(1.4 * 287.05 * T0_)
    ratios = [bu(M_, r0_, T0_, M_ * a0_) / rm.cd_generic(M_) for M_ in np.linspace(0.1, 5.0, 99)]
    put('buRatioLo', min(ratios), '{:.2f}'); put('buRatioHi', max(ratios), '{:.2f}')
    sup = [bu(M_, r0_, T0_, M_ * a0_) / rm.cd_generic(M_) for M_ in np.linspace(1.5, 5.0, 36)]
    put('buSupLo', min(sup), '{:.2f}'); put('buSupHi', max(sup), '{:.2f}')
    put('buPeak', max(bu(M_, r0_, T0_, M_ * a0_) for M_ in np.linspace(0.9, 1.3, 81)), '{:.2f}')
    put('buSub', bu(0.3, r0_, T0_, 0.3 * a0_), '{:.2f}'); put('buFive', bu(5.0, r0_, T0_, 5.0 * a0_), '{:.2f}')
    bd = dec(M0, MP, 30e3, cd_fn=bu); bd.pop('cd_fn', None)
    RES['dec_buildup'] = bd
    put('buGnd', bd['ground'], '{:.0f}'); put('buBal', bd['balloon'], '{:.0f}'); put('buRatio', bd['balloon'] / bd['ground'], '{:.1f}')
    for k in ('head', 'drag', 'nozzle'):
        put('buF' + k.title(), 100 * bd['f_' + k], '{:.0f}')
    put('buGndChange', 100 * (bd['ground'] / base['ground'] - 1), '{:+.0f}'); put('buBalChange', 100 * (bd['balloon'] / base['balloon'] - 1), '{:+.2f}')
    put('buWetBody', bu.Awb, '{:.2f}'); put('buWetFins', bu.Awf, '{:.2f}')
    # ---- real-gas helium (NIST WebBook, 290 K): density 43.571 kg/m3 at 30 MPa
    Zhe = 30e6 / (43.571 * rm.R_HE * 290.0)
    put('heZ', Zhe, '{:.3f}')
    mhe_, mb_, _ = rm.pressurant(3e6, 2.4, 30e6, 3.6e6, 2700.0, 276e6)
    fAl_rg = 1.5 * rm.tank_fraction(3e6, 2.4, 2700.0, 276e6) + mhe_ + Zhe * mb_
    put('feedAlRg', fAl_rg, '{:.3f}')
    crg = rm.close_vehicle(30e3, 80.0, fAl_rg); crg0 = rm.close_vehicle(0.0, 80.0, fAl_rg)
    put('cBalMRg', crg['m0'], '{:.0f}'); put('cGndMRg', crg0['m0'], '{:.0f}')
    # ---- NASA zero-pressure balloons: effective areal density on the equal-volume sphere
    for tag, mcf, lb in (('A', 29.47, 3600.0), ('B', 59.84, 2751.0)):
        Vn = mcf * 1e6 * 0.0283168466; An = (36 * math.pi) ** (1 / 3) * Vn ** (2 / 3)
        put('nasaSig' + tag, 1000 * lb * LB / An, '{:.0f}'); put('nasaM' + tag, lb * LB, '{:.0f}')
        TEX['nasaV' + tag] = thousands(round(Vn, -3)); RES['nasaV' + tag] = Vn
    # ---- decomposition at the drag scale factors preferred by the two validation sets
    for cds_, tag in ((RES.get('deaCdFit', 1.16), 'Dea'), (2.0, 'Two')):
        dd = dec(M0, MP, 30e3, cd_scale=float(cds_)); RES['dec_cd' + tag] = dd
        put('cd' + tag + 'Gnd', dd['ground'], '{:.1f}'); put('cd' + tag + 'Bal', dd['balloon'], '{:.0f}')
        put('cd' + tag + 'Ratio', dd['balloon'] / dd['ground'], '{:.1f}')
        for k in ('head', 'drag', 'nozzle'):
            put('cd' + tag + 'F' + k.title(), 100 * dd['f_' + k], '{:.0f}')
    # ---- ground nozzle allowed to over-expand (exit pressure down to a fraction `sep` of ambient)
    SEPS = (1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4)
    sw = [(sp_, rm.vehicle_run(M0, MP, 0.0, sep=sp_)) for sp_ in SEPS]
    RES['sep_sweep'] = [dict(sep=a_, eps=b_['eps'], apogee=b_['apogee'] / 1e3) for a_, b_ in sw]
    sbest, gbest = max(sw, key=lambda x: x[1]['apogee'])
    sp = dec(M0, MP, 30e3, sep=sbest); RES['dec_sep_best'] = sp
    put('sepBest', sbest, '{:.1f}'); put('sepEps', gbest['eps'], '{:.1f}'); put('sepGnd', sp['ground'], '{:.1f}'); put('gndApOne', base['ground'], '{:.1f}')
    put('sepRatio', sp['balloon'] / sp['ground'], '{:.1f}')
    for k in ('head', 'drag', 'nozzle'):
        put('sepF' + k.title(), 100 * sp['f_' + k], '{:.0f}')
    lim = dec(M0, MP, 30e3, sep=0.4); RES['dec_sep_limit'] = lim
    put('sepLimEps', sw[-1][1]['eps'], '{:.1f}'); put('sepLimGnd', lim['ground'], '{:.1f}'); put('sepLimRatio', lim['balloon'] / lim['ground'], '{:.1f}')
    put('sepLimFNozzle', 100 * lim['f_nozzle'], '{:.0f}')
    # both design freedoms together: any over-expansion down to 0.4 and each vehicle optimised
    cand = []
    for sp_ in SEPS:
        for fin in (10.0, 12.5, 15.0, 17.5, 20.0):
            for tw in np.arange(1.5, tw_max + 1e-9, 0.25):
                try:
                    cand.append(rm.vehicle_run(M0, MP, 0.0, d=rm.diameter(M0, fineness=fin), tw=float(tw), sep=sp_)['apogee'] / 1e3)
                except ValueError:
                    pass
    put('sepOptGnd', max(cand), '{:.0f}'); put('sepOptRatio', best['Bal'][0] / max(cand), '{:.1f}')
    # ---- similarity study (results/similarity.json, written by similarity.py)
    S = json.load(open(os.path.join(ROOT, 'results', 'similarity.json')))
    E = S['ensemble']
    g = lambda k: np.array([e[k] for e in E])
    TEX['ensN'] = str(len(E)); RES['ensN'] = len(E)
    pi0 = g('pi0'); L = np.log
    put('ensPiLo', pi0.min(), '{:.2f}'); put('ensPiHi', pi0.max(), '{:.1f}')
    ldf = g('ld_g') / g('dv_g'); Rr = g('balloon') / g('ground')
    def fit(X, y):
        X = np.column_stack([np.ones(len(y))] + X); b_, *_ = np.linalg.lstsq(X, y, rcond=None); r = y - X @ b_
        return 1 - r.var() / y.var(), r.std(), b_, X @ b_
    one = lambda: [L(pi0), L(pi0) ** 2]
    two = lambda: one() + [L(g('psi')), g('xi')]
    full = lambda: two() + [L(g('pc') / 101325.0), g('z0') / 7000.0, L(g('pic'))]
    rows = []
    fits = {}
    for tag, lab, y, scale, unit in (('Ld', 'Ground drag loss / thrust $\\Delta v$ (log)', L(ldf), 100, '\\%'), ('Sd', 'Drag share $S_D$', g('f_drag'), 100, 'pp'),
                                     ('Sn', 'Nozzle share $S_N$', g('f_noz'), 100, 'pp'), ('Sh', 'Head-start share $S_H$', g('f_head'), 100, 'pp'),
                                     ('R', 'Apogee ratio $R$ (log)', L(Rr), 100, '\\%')):
        r1 = fit(one(), y); r2 = fit(two(), y); r3 = fit(full(), y)
        fits[tag] = (r1, r3)
        for nm, r in (('One', r1), ('Two', r2), ('Full', r3)):
            put('fit' + tag + nm, r[0], '{:.2f}'); put('fit' + tag + nm + 'Sd', scale * r[1], '{:.1f}')
        rows.append('%s & %.2f & %.1f\\,%s & %.2f & %.1f\\,%s & %.2f & %.1f\\,%s \\\\' % (lab, r1[0], scale * r1[1], unit, r2[0], scale * r2[1], unit, r3[0], scale * r3[1], unit))
    with open(os.path.join(ROOT, 'paper', 'tab_collapse.tex'), 'w') as f:
        f.write('\n'.join(rows) + '\n')
    # out-of-sample check of the fits: five-fold cross-validation
    def cv(Xf, y, k=5):
        X = np.column_stack([np.ones(len(y))] + Xf())
        idx = np.random.default_rng(11).permutation(len(y)); pred = np.empty(len(y))
        for f_ in range(k):
            te = idx[f_::k]; tr = np.setdiff1d(idx, te)
            b_, *_ = np.linalg.lstsq(X[tr], y[tr], rcond=None); pred[te] = X[te] @ b_
        return 1 - np.var(y - pred) / np.var(y), np.std(y - pred)
    c1 = cv(one, L(ldf)); c2 = cv(full, g('f_drag')); c3 = cv(full, L(Rr)); c4 = cv(one, g('f_drag'))
    put('cvLd', c1[0], '{:.2f}'); put('cvSd', c2[0], '{:.2f}'); put('cvSdSd', 100 * c2[1], '{:.1f}'); put('cvR', c3[0], '{:.2f}'); put('cvSdOne', c4[0], '{:.2f}')
    # numerical convergence of the integrator
    nzb = rm.Nozzle(30e3, 4 * M0 * rm.G0, 0.40); nzg = rm.Nozzle(0.0, 4 * M0 * rm.G0, 0.40)
    conv = []
    for nz_, z_ in ((nzb, 30e3), (nzg, 0.0)):
        a8 = rm.fly(M0, MP, 0.40, nz_.mdot, nz_.isp_vac, nz_.Ae, z_kin=z_)['apogee']
        a11 = rm.fly(M0, MP, 0.40, nz_.mdot, nz_.isp_vac, nz_.Ae, z_kin=z_, rtol=1e-11)['apogee']
        conv.append(abs(a8 - a11))
    put('convMax', max(conv), '{:.1f}')
    # coefficients of the six-group fits, for reuse
    names = ['1', '$\\ln\\Pi_D$', '$(\\ln\\Pi_D)^2$', '$\\ln\\psi$', '$\\xi$', '$\\ln(p_c/p_0)$', '$z_0/H$', '$\\ln\\Pi_c$']
    cols = [('$S_D$', g('f_drag')), ('$S_N$', g('f_noz')), ('$S_H$', g('f_head')), ('$\\ln R$', L(Rr))]
    coefs = [fit(full(), y_)[2] for _, y_ in cols]
    with open(os.path.join(ROOT, 'paper', 'tab_coef.tex'), 'w') as f_:
        f_.write('\n'.join('%s & %s \\\\' % (names[i], ' & '.join('$%+.4f$' % c_[i] for c_ in coefs)) for i in range(8)) + '\n')
    RES['fit_coefficients'] = {n_.strip('$'): [float(c_[i]) for c_ in coefs] for i, n_ in enumerate(names)}
    c1_ = fit(one(), L(ldf))[2]
    put('ldA', c1_[0], '{:+.3f}'); put('ldB', c1_[1], '{:+.3f}'); put('ldC', c1_[2], '{:+.3f}')
    put('rngPsiLo', g('psi').min(), '{:.1f}'); put('rngPsiHi', g('psi').max(), '{:.1f}'); put('rngPicLo', g('pic').min(), '{:.0f}'); put('rngPicHi', g('pic').max(), '{:.0f}')
    sl = np.polyfit(L(pi0), L(ldf), 1)
    put('ldSlope', sl[0], '{:.2f}'); put('ldCoef', math.exp(sl[1]), '{:.3f}')
    from scipy.stats import spearmanr as sp_
    put('ldRho', sp_(pi0, ldf)[0], '{:.2f}'); put('sdRho', sp_(pi0, g('f_drag'))[0], '{:.2f}'); put('rRho', sp_(pi0, Rr)[0], '{:.2f}')
    put('ensLdLo', 100 * np.percentile(ldf, 5), '{:.0f}'); put('ensLdHi', 100 * np.percentile(ldf, 95), '{:.0f}')
    put('ensLdBalMed', 100 * np.median(g('ld_b') / g('dv_b')), '{:.1f}')
    put('ensDragLargest', 100 * np.mean((g('f_drag') > g('f_noz')) & (g('f_drag') > g('f_head'))), '{:.0f}')
    put('ensSdLo', 100 * np.percentile(g('f_drag'), 5), '{:.0f}'); put('ensSdHi', 100 * np.percentile(g('f_drag'), 95), '{:.0f}')
    put('ensRLo', np.percentile(Rr, 5), '{:.1f}'); put('ensRHi', np.percentile(Rr, 95), '{:.0f}')
    sol = g('eng') == 3
    put('ensSolidOff', 100 * np.mean((g('f_drag') - fits['Sd'][1][3])[sol]), '{:+.1f}')
    RES['ens_fit_sd'] = fits['Sd'][1][3].tolist()
    # optimum float altitude
    import collections
    fam = collections.defaultdict(list)
    for r in S['system']:
        fam[(r['mat'], r['mfix'], r['target'], r['sigma'], r['gf'])].append(r)
    opt = []
    for (mat, mfix, tgt, sg_, gf_), v in sorted(fam.items()):
        v = sorted(v, key=lambda r: r['z'])
        if any(r['m0'] is None for r in v):
            continue
        z = np.array([r['z'] for r in v]); tot = np.array([r['total'] for r in v])
        i = int(np.argmin(tot[1:])) + 1
        if 1 < i < len(z) - 1:
            c = np.polyfit(z[i - 1:i + 2], tot[i - 1:i + 2], 2); zst = -c[1] / (2 * c[0]); tst = float(np.polyval(c, zst))
        else:
            zst, tst = z[i], tot[i]
        above = [zz for zz, tt in zip(z[1:], tot[1:]) if tt > tot[0]]
        within = z[1:][tot[1:] <= 1.05 * tot[i]]
        A = math.pi * v[0]['d'] ** 2 / 4
        opt.append(dict(mat=mat, mfix=mfix, target=tgt / 1e3, sigma=sg_, gf=gf_, zgrid=float(z[i]), total_grid=float(tot[i]), ground=float(tot[0]), zstar=float(zst), total=tst, rocket=v[i]['m0'], V=v[i]['V'],
                        pi0=CDREF * A * 101325.0 / (v[0]['m0'] * rm.G0), lo=float(within.min()), hi=float(within.max()),
                        breakeven=float(min(above)) if above else None, curve=[z.tolist(), (tot / tot[0]).tolist()]))
    var = [o for o in opt if (o['sigma'], o['gf']) != (0.030, 0.10)]
    opt = [o for o in opt if (o['sigma'], o['gf']) == (0.030, 0.10)]
    RES['optimum'] = opt
    RES['optimum_var'] = var
    put('zsVarLo', min(o['zstar'] for o in var), '{:.0f}'); put('zsVarHi', max(o['zstar'] for o in var), '{:.0f}')
    for sg_, gf_, tag in ((0.020, 0.10, 'SigLo'), (0.040, 0.10, 'SigHi'), (0.030, 0.05, 'GonLo'), (0.030, 0.20, 'GonHi')):
        o = [x for x in var if x['sigma'] == sg_ and x['gf'] == gf_ and x['mfix'] == 80.0][0]
        put('zs' + tag, o['zstar'], '{:.1f}'); put('zs' + tag + 'Save', 100 * (1 - o['total'] / o['ground']), '{:.0f}')
    o100 = [o for o in opt if o['target'] == 100.0]
    put('zsLo', min(o['zstar'] for o in o100), '{:.0f}'); put('zsHi', max(o['zstar'] for o in o100), '{:.0f}')
    put('zsAllLo', min(o['zstar'] for o in opt), '{:.0f}'); put('zsAllHi', max(o['zstar'] for o in opt), '{:.0f}')
    put('zsMassLo', min(o['ground'] for o in o100), '{:.0f}'); TEX['zsMassHi'] = thousands(max(o['ground'] for o in o100)); RES['zsMassHi'] = max(o['ground'] for o in o100)
    put('zsSaveLo', 100 * (1 - max(o['total'] / o['ground'] for o in o100)), '{:.0f}'); put('zsSaveHi', 100 * (1 - min(o['total'] / o['ground'] for o in o100)), '{:.0f}')
    put('zsBandLo', max(o['lo'] for o in o100), '{:.1f}'); put('zsBandHi', min(o['hi'] for o in o100), '{:.1f}')
    zfit = np.polyfit([math.log(o['ground']) for o in o100 if o['mat'] == 'Al6061T6'], [o['zstar'] for o in o100 if o['mat'] == 'Al6061T6'], 1)
    put('zsSlope', zfit[0] * math.log(10), '{:.1f}')
    b80 = [o for o in o100 if o['mat'] == 'Al6061T6' and o['mfix'] == 80.0][0]
    put('zsBase', b80['zstar'], '{:.0f}'); put('zsBaseTot', b80['total'], '{:.0f}'); put('zsBaseGnd', b80['ground'], '{:.0f}'); put('zsBaseRocket', b80['rocket'], '{:.0f}')
    TEX['zsBaseV'] = thousands(round(b80['V'], -2)); RES['zsBaseV'] = b80['V']
    put('zsBaseGrid', b80['zgrid'], '{:.1f}'); put('zsBaseTotGrid', b80['total_grid'], '{:.0f}')
    be = [o['breakeven'] for o in o100 if o['breakeven']]
    put('beLo', min(be), '{:.1f}') if be else None; put('beHi', max(be), '{:.1f}') if be else None
    TEX['beCount'] = str(len(be)); TEX['optCount'] = str(len(o100))
    rows = []
    for o in o100:
        rows.append('%s & %.0f & %s & %.1f & %.1f & %s & %.0f & %s \\\\' % ('Al' if o['mat'].startswith('Al') else 'steel', o['mfix'], thousands(o['ground']), o['zstar'], o['zgrid'], thousands(o['total_grid']), o['rocket'], thousands(round(o['V'], -1))))
    with open(os.path.join(ROOT, 'paper', 'tab_optimum.tex'), 'w') as f:
        f.write('\n'.join(rows) + '\n')
    # optimum under a load-scaled envelope and with a thermal jacket (recomputed from the closed vehicles)
    def fam_opt(v, tot):
        z = np.array([r['z'] for r in v]); tot = np.asarray(tot)
        i = int(np.argmin(tot[1:])) + 1
        if 1 < i < len(z) - 1:
            c = np.polyfit(z[i - 1:i + 2], tot[i - 1:i + 2], 2); return -c[1] / (2 * c[0]), float(np.polyval(c, -c[1] / (2 * c[0]))), tot[0]
        return z[i], tot[i], tot[0]
    base_fams = {k: sorted(v, key=lambda r: r['z']) for k, v in fam.items() if k[3:] == (0.030, 0.10) and k[2] == 100e3 and all(r['m0'] is not None for r in v)}
    def variant(k_load=0.0, thermal=False):
        out = []
        for (mat, mfix, *_), v in sorted(base_fams.items()):
            tot = []
            for r in v:
                if r['z'] == 0:
                    tot.append(r['m0']); continue
                ms = 1.10 * r['m0'] + (rm.insulation_mass(r['mp'], r['z'] * 1e3)[0] if thermal else 0.0)
                V_, me_, mg_ = rm.balloon(ms, r['z'] * 1e3, k_load=k_load)
                tot.append(ms + me_ + 1.1 * mg_)
            zs_, ts_, g_ = fam_opt(v, tot)
            out.append(dict(mat=mat, mfix=mfix, zstar=float(zs_), save=float(1 - ts_ / g_), ground=float(g_)))
        return out
    chk = variant()
    assert max(abs(a_['zstar'] - b_['zstar']) for a_, b_ in zip(chk, sorted(o100, key=lambda o: (o['mat'], o['mfix'])))) < 1e-6
    RES['optimum_kload'] = {}
    for kl, tag in ((0.05, 'KlA'), (0.10, 'KlB'), (0.20, 'KlC')):
        vv = variant(k_load=kl); RES['optimum_kload'][str(kl)] = vv
        put('zs' + tag + 'Lo', min(o['zstar'] for o in vv), '{:.0f}'); put('zs' + tag + 'Hi', max(o['zstar'] for o in vv), '{:.0f}')
        put('zs' + tag + 'SaveLo', 100 * min(o['save'] for o in vv), '{:.0f}'); put('zs' + tag + 'SaveHi', 100 * max(o['save'] for o in vv), '{:.0f}')
        al = [o for o in vv if o['mat'] == 'Al6061T6']
        put('zs' + tag + 'Slope', np.polyfit([math.log(o['ground']) for o in al], [o['zstar'] for o in al], 1)[0] * math.log(10), '{:.1f}')
    vv = variant(thermal=True); RES['optimum_thermal'] = vv
    put('zsThLo', min(o['zstar'] for o in vv), '{:.0f}'); put('zsThHi', max(o['zstar'] for o in vv), '{:.0f}')
    put('zsThShift', max(abs(a_['zstar'] - b_['zstar']) for a_, b_ in zip(vv, chk)), '{:.2f}')
    put('zsThSaveLo', 100 * min(o['save'] for o in vv), '{:.0f}'); put('zsThSaveHi', 100 * max(o['save'] for o in vv), '{:.0f}')
    jk = [(rm.insulation_mass(r['mp'], 20e3), r) for k_, v in base_fams.items() for r in v if r['z'] == 20.0]
    put('insFracHi', 100 * max(m_[0] / r['m0'] for m_, r in jk), '{:.1f}'); put('insFracLo', 100 * min(m_[0] / r['m0'] for m_, r in jk), '{:.2f}')
    put('insThickReq', 1e3 * max(rm.insulation_mass(r['mp'], 20e3, t_min=0.0)[1] for _, r in jk), '{:.0f}')
    # fair ground comparator for the mass optimum: each ground vehicle designed for ground launch
    fair_gnd = {}
    for (mat, mfix, *_), v in sorted(base_fams.items()):
        best_ = None
        for fin in (10.0, 15.0, 20.0):
            for tw in (2.0, 3.0, 4.0, 5.0, 6.0):
                for sp_ in (1.0, 0.7):
                    try:
                        r = rm.close_vehicle(0.0, mfix, rm.feed_fraction(mat), target=100e3, tw=tw, fineness=fin, sep=sp_)
                    except ValueError:
                        r = None
                    if r and (best_ is None or r['m0'] < best_[0]):
                        best_ = (r['m0'], fin, tw, sp_)
        fair_gnd[(mat, mfix)] = best_
    RES['fair_ground'] = {'%s/%g' % k: list(v) for k, v in fair_gnd.items()}
    # the balloon-launched vehicles redesigned for their own launch as well (body and thrust)
    fair_bal = {}
    for (mat, mfix, *_), v in sorted(base_fams.items()):
        for r in v:
            if r['z'] == 0:
                continue
            best_ = (r['m0'], 10.0, 4.0)
            for fin in (10.0, 15.0):
                for tw in (4.0, 5.0, 6.0):
                    if (fin, tw) == (10.0, 4.0):
                        continue
                    try:
                        c = rm.close_vehicle(r['z'] * 1e3, mfix, rm.feed_fraction(mat), target=100e3, tw=tw, fineness=fin)
                    except ValueError:
                        c = None
                    if c and c['m0'] < best_[0]:
                        best_ = (c['m0'], fin, tw)
            fair_bal[(mat, mfix, r['z'])] = best_
    RES['fair_balloon'] = {'%s/%g/%g' % k: list(v) for k, v in fair_bal.items()}
    fair = [dict(mat=o['mat'], mfix=o['mfix'], save=1 - (1 - o['save']) * o['ground'] / fair_gnd[(o['mat'], o['mfix'])][0]) for o in chk]
    fb = []
    for (mat, mfix, *_), v in sorted(base_fams.items()):
        tot = [fair_gnd[(mat, mfix)][0]]
        for r in v[1:]:
            ms = 1.10 * fair_bal[(mat, mfix, r['z'])][0]
            V_, me_, mg_ = rm.balloon(ms, r['z'] * 1e3); tot.append(ms + me_ + 1.1 * mg_)
        zs_, ts_, g_ = fam_opt(v, tot)
        fb.append(dict(mat=mat, mfix=mfix, zstar=float(zs_), save=float(1 - ts_ / g_)))
    RES['optimum_fair_both'] = fb
    put('fbLo', min(o['zstar'] for o in fb), '{:.0f}'); put('fbHi', max(o['zstar'] for o in fb), '{:.0f}')
    put('fbSaveLo', 100 * min(o['save'] for o in fb), '{:.0f}'); put('fbSaveHi', 100 * max(o['save'] for o in fb), '{:.0f}')
    put('fairSaveLo', 100 * min(o['save'] for o in fair), '{:.0f}'); put('fairSaveHi', 100 * max(o['save'] for o in fair), '{:.0f}')
    put('fairGndRedLo', 100 * min(1 - fair_gnd[(o['mat'], o['mfix'])][0] / o['ground'] for o in chk), '{:.0f}')
    put('fairGndRedHi', 100 * max(1 - fair_gnd[(o['mat'], o['mfix'])][0] / o['ground'] for o in chk), '{:.0f}')
    # hardware only the rockoon needs: divergent-section mass (both vehicles), and a heavier gondola
    # carrying the launcher and attitude-control system
    def hw_variant(areal=10.0, gf=0.10, k_load=0.0, fair_ground=False):
        out = []
        for (mat, mfix, *_), v in sorted(base_fams.items()):
            f = rm.feed_fraction(mat); tot = []
            for r in v:
                kw = {}
                if r['z'] == 0 and fair_ground:
                    m0g, fin, tw, sp_ = fair_gnd[(mat, mfix)]; kw = dict(tw=tw, fineness=fin, sep=sp_)
                elif fair_ground:
                    _, fin, tw = fair_bal[(mat, mfix, r['z'])]; kw = dict(tw=tw, fineness=fin)
                c = rm.close_vehicle_skirt(r['z'] * 1e3, mfix, f, target=100e3, areal=areal, **kw) if areal > 0 else rm.close_vehicle(r['z'] * 1e3, mfix, f, target=100e3, **kw)
                if r['z'] == 0:
                    tot.append(c['m0']); continue
                ms = (1.0 + gf) * c['m0']
                V_, me_, mg_ = rm.balloon(ms, r['z'] * 1e3, k_load=k_load)
                tot.append(ms + me_ + 1.1 * mg_)
            zs_, ts_, g_ = fam_opt(v, tot)
            out.append(dict(mat=mat, mfix=mfix, zstar=float(zs_), save=float(1 - ts_ / g_)))
        return out
    for tag, kw in (('Sk', dict(areal=10.0)), ('SkHi', dict(areal=20.0)), ('Gf', dict(areal=0.0, gf=0.30)),
                    ('All', dict(areal=20.0, gf=0.30, k_load=0.20, fair_ground=True))):
        vv = hw_variant(**kw); RES['optimum_hw_' + tag] = vv
        put('hw' + tag + 'Lo', min(o['zstar'] for o in vv), '{:.0f}'); put('hw' + tag + 'Hi', max(o['zstar'] for o in vv), '{:.0f}')
        put('hw' + tag + 'SaveLo', 100 * min(o['save'] for o in vv), '{:.0f}'); put('hw' + tag + 'SaveHi', 100 * max(o['save'] for o in vv), '{:.0f}')
    al_ = [o for o in RES['optimum_hw_All'] if o['mat'] == 'Al6061T6']; ss_ = [o for o in RES['optimum_hw_All'] if o['mat'] != 'Al6061T6']
    put('hwAllAlLo', 100 * min(o['save'] for o in al_), '{:.0f}'); put('hwAllAlHi', 100 * max(o['save'] for o in al_), '{:.0f}')
    put('hwAllSsLo', 100 * min(o['save'] for o in ss_), '{:.0f}'); put('hwAllSsHi', 100 * max(o['save'] for o in ss_), '{:.0f}')
    TEX['hwAllAlNeg'] = str(sum(o['save'] < 0 for o in al_))
    rr_ = rm.vehicle_run(M0, MP, 30e3); sk = rm.nozzle_skirt(rr_['Ae'], rr_['eps'])
    put('skBaseM', sk[0], '{:.1f}'); put('skBaseL', sk[1], '{:.2f}')
    gg_ = rm.vehicle_run(M0, MP, 0.0); skg = rm.nozzle_skirt(gg_['Ae'], gg_['eps'])
    put('skGndM', skg[0], '{:.1f}'); put('skGndL', skg[1], '{:.2f}')
    # a cost weighting: balloon-system mass (envelope, gas, gondola) priced at r times the rocket per kg
    for r_c, tag in ((0.3, 'A'), (1.0, 'B'), (3.0, 'C'), (10.0, 'D')):
        res_ = []
        for (mat, mfix, *_), v in sorted(base_fams.items()):
            cost = [r['m0'] if r['z'] == 0 else r['m0'] + r_c * (0.10 * r['m0'] + r['env'] + 1.1 * r['gas']) for r in v]
            zs_, ts_, g_ = fam_opt(v, cost)
            res_.append(dict(zstar=float(zs_), save=float(1 - ts_ / g_)))
        RES['optimum_cost_' + tag] = res_
        put('cost' + tag + 'Lo', min(o['zstar'] for o in res_), '{:.0f}'); put('cost' + tag + 'Hi', max(o['zstar'] for o in res_), '{:.0f}')
        put('cost' + tag + 'SaveLo', 100 * min(o['save'] for o in res_), '{:.0f}'); put('cost' + tag + 'SaveHi', 100 * max(o['save'] for o in res_), '{:.0f}')
    b_ = [r for r in base_fams[('Al6061T6', 80.0, 100e3, 0.030, 0.10)] if r['z'] == 17.5][0]
    put('insBase', rm.insulation_mass(b_['mp'], 17.5e3)[0], '{:.1f}')

    # further tests (results/extra.json, written by extra_tests.py)
    X = json.load(open(os.path.join(ROOT, 'results', 'extra.json')))
    xs = X['sobol_vehicles']; RES['sobol_vehicles'] = xs
    vnames = []
    for x in xs:
        if x['vehicle'] not in vnames:
            vnames.append(x['vehicle'])
    cross_z, cross_pi = [], []
    for n in vnames:
        v = [x for x in xs if x['vehicle'] == n]
        for p_, q_ in zip(v, v[1:]):
            if p_['ST_cd'] >= 0.01 > q_['ST_cd']:
                f_ = (math.log(p_['ST_cd']) - math.log(0.01)) / (math.log(p_['ST_cd']) - math.log(max(q_['ST_cd'], 1e-6)))
                cross_z.append(p_['z'] + f_ * (q_['z'] - p_['z']))
                cross_pi.append(math.exp(math.log(p_['piz']) + f_ * (math.log(q_['piz']) - math.log(p_['piz']))))
                break
    put('h3ZLo', min(cross_z), '{:.0f}'); put('h3ZHi', max(cross_z), '{:.0f}')
    put('h3PiLo', min(cross_pi), '{:.2f}'); put('h3PiHi', max(cross_pi), '{:.2f}')
    put('h3PiZeroLo', min(x['pi0'] for x in xs), '{:.1f}'); put('h3PiZeroHi', max(x['pi0'] for x in xs), '{:.1f}')
    yy = np.array([max(x['ST_cd'], 0) for x in xs])
    def r2q(xv):
        c = np.polyfit(xv, yy, 3); r = yy - np.polyval(c, xv); return 1 - r.var() / yy.var()
    put('h3RPi', r2q(np.log([x['piz'] for x in xs])), '{:.2f}'); put('h3RZ', r2q(np.array([x['z'] for x in xs])), '{:.2f}')
    put('h3NS', 512, '{:d}')
    for zf, tag in (('30', 'Thirty'), ('15', 'Fifteen')):
        E = X['fixed_float'][zf]
        lp = np.log([e['pi0'] for e in E])
        for key, nm in (('f_drag', 'Sd'), ('f_noz', 'Sn'), ('f_head', 'Sh')):
            yv = np.array([e[key] for e in E])
            Xm = np.column_stack([np.ones(len(yv)), lp, lp ** 2]); b = np.linalg.lstsq(Xm, yv, rcond=None)[0]; r = yv - Xm @ b
            put('fx' + tag + nm, 1 - r.var() / yv.var(), '{:.2f}'); put('fx' + tag + nm + 'Sd', 100 * r.std(), '{:.1f}')
        put('fx' + tag + 'N', len(E), '{:d}')

    # ---- cross-checks against NASA CEA run directly (cea_check.py) and against RocketPy (rocketpy_check.py)
    C = json.load(open(os.path.join(ROOT, 'results', 'cea_check.json')))
    cn = C['cantera']; ds, dfz = C['design_shift'], C['design_frozen']
    put('ceaTc', ds['T'][0], '{:.0f}'); put('ceaCstar', ds['cstar'], '{:.0f}')
    put('ceaIspG', ds['isp_vac'][2], '{:.1f}'); put('ceaIspB', ds['isp_vac'][3], '{:.1f}'); put('ceaIspFr', dfz['isp_vac'][3], '{:.1f}')
    put('cnIspG', cn['isp_vac_g'], '{:.1f}'); put('cnIspB', cn['isp_vac_b'], '{:.1f}')
    dif = [abs(ds['isp_vac'][2] / cn['isp_vac_g'] - 1), abs(ds['isp_vac'][3] / cn['isp_vac_b'] - 1), abs(dfz['isp_vac'][3] / cn['frozen_b'] - 1),
           abs(ds['cstar'] / cn['cstar'] - 1)]
    dif += [abs(s_['cea_shift'] / s_['cantera_shift'] - 1) for s_ in C['sweep'] if s_['cantera_shift']]
    dif += [abs(s_['cea_frozen'] / s_['cantera_frozen'] - 1) for s_ in C['sweep'] if s_['cantera_frozen']]
    put('ceaMaxDiff', 100 * max(dif), '{:.2f}')
    put('ceaBench', C['bench_shift']['isp'][-1], '{:.1f}'); put('ceaBenchFr', C['bench_frozen']['isp'][-1], '{:.1f}')
    put('ceaTcDiff', abs(ds['T'][0] - cn['Tc']), '{:.1f}')
    P_ = json.load(open(os.path.join(ROOT, 'results', 'rocketpy_check.json')))
    for tag, k in (('Gnd', 'ground'), ('Bal', 'balloon')):
        a3, a6 = P_[k + '_3DOF'], P_[k + '_6DOF']
        put('rp' + tag + 'Three', a3['apogee'] / 1e3, '{:.2f}' if tag == 'Gnd' else '{:.1f}'); put('rp' + tag + 'Six', a6['apogee'] / 1e3, '{:.2f}' if tag == 'Gnd' else '{:.1f}')
        put('rp' + tag + 'Ours', a3['ours_apogee'] / 1e3, '{:.2f}' if tag == 'Gnd' else '{:.1f}')
        put('rp' + tag + 'Diff', 100 * max(abs(a3['apogee'] / a3['ours_apogee'] - 1), abs(a6['apogee'] / a3['ours_apogee'] - 1)), '{:.2f}')
        put('rp' + tag + 'Vdiff', abs(a3['v_max'] - a3['ours_vb']), '{:.1f}')
    put('rpMaxDiff', max(RES['rpGndDiff'], RES['rpBalDiff']), '{:.2f}')
    # ---- tests added in the revision
    import revision_analysis
    revision_analysis.main(put, RES, TEX)
    # ---- the envelope rule for any envelope law; published dry-mass data, fineness-dependent drag and the Aerobee check
    import general_rule
    general_rule.main(put, RES, TEX)
    import realdata_analysis
    realdata_analysis.main(put, RES, TEX)

    # Sobol
    RES['sobol'] = S['sobol']
    sb = {x['z']: x for x in S['sobol']}
    for z, tag in ((0, 'Zero'), (5, 'Five'), (10, 'Ten'), (15, 'Fifteen'), (20, 'Twenty'), (30, 'Thirty')):
        put('sobCd' + tag, max(sb[z]['ST'][0], 0.0), '{:.2f}'); put('sobEta' + tag, sb[z]['ST'][1], '{:.2f}'); put('sobDry' + tag, sb[z]['ST'][3], '{:.2f}')
        put('sobPi' + tag, sb[z]['piz'], '{:.2f}' if sb[z]['piz'] > 0.1 else '{:.3f}')
    TEX['sobN'] = thousands(1024)
    put('sobCdZeroLo', sb[0]['ST_ci'][0][0], '{:.2f}'); put('sobCdZeroHi', sb[0]['ST_ci'][0][1], '{:.2f}')
    put('sobCdTenLo', max(sb[10]['ST_ci'][0][0], 0), '{:.2f}'); put('sobCdTenHi', sb[10]['ST_ci'][0][1], '{:.2f}')
    put('sobCdTwentyHi', sb[20]['ST_ci'][0][1], '{:.3f}')

    os.makedirs(os.path.join(ROOT, 'results'), exist_ok=True)
    with open(os.path.join(ROOT, 'results', 'results.json'), 'w') as f:
        json.dump(RES, f, indent=1, default=float)
    with open(os.path.join(ROOT, 'paper', 'results.tex'), 'w') as f:
        for k in sorted(TEX):
            f.write('\\expandafter\\def\\csname res@%s\\endcsname{%s}\n' % (k, TEX[k]))
    print('wrote', len(TEX), 'values')


if __name__ == '__main__':
    main()
