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
    gnd = [deacon(w, 19680, 1.0, z0=0.0, elev=90.0, coast=True)['apogee'] / 1e3 for w in (190.50, 232.65)]
    put('deaGroundLight', gnd[0], '{:.0f}'); put('deaGroundHeavy', gnd[1], '{:.0f}')
    # ---- baseline
    cases = {'bal': rm.vehicle_run(M0, MP, 30e3), 'gnd': rm.vehicle_run(M0, MP, 0.0),
             'fall': rm.vehicle_run(M0, MP, 29e3, v0=-math.sqrt(2 * 9.72 * 1000.0))}
    put('fallSpeed', math.sqrt(2 * 9.72 * 1000.0), '{:.0f}')
    RES['cases'] = {k: {a: float(b) if not isinstance(b, str) else b for a, b in v.items()} for k, v in cases.items()}
    for k, r in cases.items():
        put(k + 'Ap', r['apogee'] / 1e3, '{:.0f}'); put(k + 'Tb', r['tb']); put(k + 'Hb', (r['apogee'] * 0 + r['h_b'] + {'bal': 30e3, 'gnd': 0.0, 'fall': 29e3}[k]) / 1e3)
        put(k + 'Vb', r['v_b'], '{:.0f}'); put(k + 'Dv', r['dv_thrust'], '{:.0f}'); put(k + 'Lg', r['loss_g'], '{:.0f}')
        put(k + 'Ld', r['loss_d'], '{:.0f}'); put(k + 'Q', r['q_max'] / 1e3, '{:.1f}' if r['q_max'] < 1e4 else '{:.0f}')
        put(k + 'Eps', r['eps'], '{:.1f}' if r['eps'] > 10 else '{:.2f}'); put(k + 'Isp', r['isp_ign']); put(k + 'IspVac', r['isp_vac'])
        put(k + 'F', r['F0'] / 1e3); put(k + 'Mdot', r['mdot'], '{:.2f}'); put(k + 'Pe', r['pe'] / 1e3, '{:.1f}'); put(k + 'Tap', r['t_apogee'], '{:.0f}')
    put('gainRatio', cases['bal']['apogee'] / cases['gnd']['apogee'], '{:.1f}')
    put('gainKm', (cases['bal']['apogee'] - cases['gnd']['apogee']) / 1e3, '{:.0f}')
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
        put('dec' + k.title(), base[k], '{:.0f}'); put('decF' + k.title(), 100 * base['f_' + k], '{:.0f}')
    for k, v in base['runs'].items():
        put('run' + k.replace('0', 'G').replace('1', 'B'), v, '{:.0f}')
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

    os.makedirs(os.path.join(ROOT, 'results'), exist_ok=True)
    with open(os.path.join(ROOT, 'results', 'results.json'), 'w') as f:
        json.dump(RES, f, indent=1, default=float)
    with open(os.path.join(ROOT, 'paper', 'results.tex'), 'w') as f:
        for k in sorted(TEX):
            f.write('\\expandafter\\def\\csname res@%s\\endcsname{%s}\n' % (k, TEX[k]))
    print('wrote', len(TEX), 'values')


if __name__ == '__main__':
    main()
