"""Check against a flown pressure-fed liquid sounding rocket, the Aerobee 150 (run with Cantera installed).

1. Thermochemistry: ideal shifting-equilibrium performance of IRFNA with aniline/furfuryl alcohol at the published
   chamber pressure, mixture ratio and area ratio, by the method of rockoon_model.expansion, against the published
   sea-level specific impulse. The ratio is an estimate of the efficiency eta for an engine of this class.
2. Trajectory: the published masses, thrust and flow rate flown with the generic drag curve, against the
   performance stated in the NASA compendia and against flights of 1964 with stated net payload weights.

Sources: NASA TR R-226 (Busse and Leffler, 1966) for the propulsion data, propellant weights and nominal
performance; Busse et al. (1967), compendium for 1964, for the flights; the item-by-item masses are those of the
Open Aerospace compilation, which agrees with TR R-226 on the oxidizer weight and on the empty weight.

Writes results/aerobee_check.json."""
import json
import math
import os

import numpy as np

import rockoon_model as rm

ROOT = os.path.dirname(rm.HERE)
LB, LBF, FT, MI, PSI = 0.45359237, 4.4482216, 0.3048, 1609.344, 6894.757
# ---- propulsion (TR R-226, nominal liquid-propulsion performance)
PC, OF, EPS = 324.0 * PSI, 2.56, 4.6
F_SL, MDOT, T_BURN = 4100.0 * LBF, 20.71 * LB, 50.9
# ---- propellants: IRFNA (type IIIA: 83.4% HNO3, 14% NO2, 2% H2O, 0.6% HF counted as water), ANFA 65/35
#      composition (atoms C, H, N, O), molar mass (kg/mol), enthalpy of formation of the liquid (J/mol)
SPECIES = {'HNO3': ((0, 1, 1, 3), 63.013e-3, -174.10e3), 'N2O4': ((0, 0, 2, 4), 92.011e-3, -19.56e3), 'H2O': ((0, 2, 0, 1), 18.015e-3, -285.83e3),
           'aniline': ((6, 7, 1, 0), 93.128e-3, 31.3e3), 'furfuryl': ((5, 6, 0, 2), 98.100e-3, -276.2e3)}
OX = {'HNO3': 0.834, 'N2O4': 0.140, 'H2O': 0.026}
FUEL = {'aniline': 0.65, 'furfuryl': 0.35}
# ---- masses (kg): Aerobee 150A sustainer without payload, and booster
SUSTAINER_ITEMS = dict(nose_cone=7.711, nose_extension=6.35, tank_shroud=3.0, head_end_plumbing=1.814, tank_assembly=74.156, aft_shroud=7.031,
                       thrust_chamber=14.56, misc=5.625, reserve_propellant=6.35, fins=12.701)
PROP = 134.354 + 343.914
HELIUM = 5.15 * LB
BOOSTER, BOOSTER_PROP, F_BOOST, T_BOOST = 153.314 + 118.841, 118.841, 18600.0 * LBF, 2.5
D = 15.0 * 0.0254
WSMR = 1219.0
# ---- published performance: (label, net payload lb, launch elevation m, peak statute miles, burnout ft/s, burnout ft, source)
NOMINAL = [('Nominal, 150 lb', 150.0, 0.0, 158.0, 6580.0, 132000.0), ('Nominal, 250 lb', 250.0, 0.0, 128.0, 5870.0, 128000.0)]
FLIGHTS = [('Flight 4.109 GG', 231.6, WSMR, 130.7, 5900.0, None), ('Flight 4.110 GG', 236.5, WSMR, 129.0, 5600.0, 130000.0),
           ('Flight 4.118 NA', 314.8, WSMR, 97.5, None, None), ('Flight 4.120 CG', 356.9, WSMR, 89.0, 4650.0, 108000.0)]


def thermo():
    gas = rm._gas()
    atoms = np.zeros(4); h = 0.0
    for mix, w in ((OX, OF / (1 + OF)), (FUEL, 1 / (1 + OF))):
        for name, x in mix.items():
            comp, M, hf = SPECIES[name]
            n = w * x / M
            atoms += n * np.array(comp); h += n * hf
    X = {'CO': atoms[0], 'H2': atoms[1] / 2, 'N2': atoms[2] / 2, 'O2': (atoms[3] - atoms[0]) / 2}
    gas.TPX = 2500.0, PC, X
    gas.HP = h, PC
    gas.equilibrate('HP', solver='gibbs', max_steps=5000)
    Tc, s = gas.T, gas.entropy_mass
    ps = PC * np.exp(np.linspace(math.log(0.75), math.log(2e-3), 200))
    u, rho = [], []
    for p in ps:
        gas.SPX = s, p, gas.X
        gas.equilibrate('SP', solver='gibbs', max_steps=5000)
        u.append(math.sqrt(max(2.0 * (h - gas.enthalpy_mass), 0.0))); rho.append(gas.density)
    u, rho = np.array(u), np.array(rho); flux = rho * u
    it = int(np.argmax(flux)); a = np.polyfit(np.log(ps[it - 1:it + 2]), flux[it - 1:it + 2], 2)
    Gt = float(np.polyval(a, -a[1] / (2 * a[0])))
    sl = slice(it + 1, None)
    eps = Gt / flux[sl]; isp_vac = (u[sl] + ps[sl] / flux[sl]) / rm.G0
    cstar = PC / Gt
    iv = float(np.interp(EPS, eps, isp_vac)); pe = float(np.interp(EPS, eps, ps[sl]))
    isp_sl = iv - 101325.0 * EPS * cstar / (PC * rm.G0)
    return dict(Tc=Tc, cstar=cstar, isp_vac_ideal=iv, isp_sl_ideal=isp_sl, pe=pe)


def fly(payload_lb, z0, Ae, cd_scale=1.0):
    isp_vac = (F_SL + 101325.0 * Ae) / (MDOT * rm.G0)
    m_sus = sum(SUSTAINER_ITEMS.values()) + PROP + HELIUM + payload_lb * LB
    m0 = m_sus + BOOSTER
    mb = BOOSTER_PROP / T_BOOST
    # booster and sustainer burn together for 2.5 s: one engine with the summed flow and vacuum thrust
    md1 = MDOT + mb
    isp1 = (F_SL + 101325.0 * Ae + F_BOOST) / (md1 * rm.G0)
    r1 = rm.fly(m0, md1 * T_BOOST, D, md1, isp1, Ae, z_kin=z0, cd_scale=cd_scale, coast=False)
    m1 = m0 - md1 * T_BOOST - (BOOSTER - BOOSTER_PROP)
    r2 = rm.fly(m1, PROP - MDOT * T_BOOST, D, MDOT, isp_vac, Ae, z_kin=z0 + r1['h_b'], v0=r1['v_b'], cd_scale=cd_scale)
    return dict(apogee=r2['apogee'], v_b=r2['v_b'], z_b=z0 + r1['h_b'] + r2['h_b'], m0=m0, isp_vac=isp_vac, loss_d=r1['loss_d'] + r2['loss_d'])


def main():
    out = dict(thermo=thermo())
    th = out['thermo']
    isp_sl_pub = F_SL / (MDOT * rm.G0)
    At = MDOT * th['cstar'] / PC; Ae = EPS * At
    isp_vac_del = (F_SL + 101325.0 * Ae) / (MDOT * rm.G0)
    out.update(isp_sl_published=isp_sl_pub, Ae=Ae, isp_vac_delivered=isp_vac_del, eta_vac=isp_vac_del / th['isp_vac_ideal'],
               eta_sl=isp_sl_pub / th['isp_sl_ideal'], feed_fraction=(SUSTAINER_ITEMS['tank_assembly'] + HELIUM) / ((758.2 + 303.3) * LB),
               tank_pressure_ratio=500.0 / 324.0)
    rows = []
    for label, pay, z0, peak, vb, zb in NOMINAL + FLIGHTS:
        r = fly(pay, z0, Ae); r16 = fly(pay, z0, Ae, cd_scale=1.16)
        rows.append(dict(label=label, payload=pay * LB, z0=z0, peak=peak * MI, vb=None if vb is None else vb * FT, zb=None if zb is None else zb * FT,
                         model_peak=r['apogee'], model_vb=r['v_b'], model_zb=r['z_b'], model_peak_cd=r16['apogee'], m0=r['m0'], loss_d=r['loss_d'],
                         nominal=(label, pay, z0, peak, vb, zb) in NOMINAL))
    out['rows'] = rows
    # sensitivity of the 250 lb case to the exit area (throat from an ideal thrust coefficient instead of c*)
    out['ae_sens'] = [fly(250.0, 0.0, k * Ae)['apogee'] for k in (0.9, 1.0, 1.1)]
    with open(os.path.join(ROOT, 'results', 'aerobee_check.json'), 'w') as f:
        json.dump(out, f, indent=1, default=float)
    print('thermo: Tc %.0f K, c* %.0f m/s, ideal Isp vac %.1f s (eps %.1f), ideal sea level %.1f s; published sea level %.1f s' % (
        th['Tc'], th['cstar'], th['isp_vac_ideal'], EPS, th['isp_sl_ideal'], isp_sl_pub))
    print('eta (vacuum) %.3f, eta (sea level) %.3f, delivered vacuum Isp %.1f s, feed fraction %.3f' % (out['eta_vac'], out['eta_sl'], isp_vac_del, out['feed_fraction']))
    for r in rows:
        print('%-26s payload %5.1f kg, z0 %4.0f m: peak %5.1f km (published %5.1f, %+5.1f%%; drag x1.16 %+5.1f%%)%s' % (
            r['label'], r['payload'], r['z0'], r['model_peak'] / 1e3, r['peak'] / 1e3, 100 * (r['model_peak'] / r['peak'] - 1),
            100 * (r['model_peak_cd'] / r['peak'] - 1),
            '' if r['vb'] is None else '; burnout %.0f m/s (published %.0f, %+.1f%%)' % (r['model_vb'], r['vb'], 100 * (r['model_vb'] / r['vb'] - 1))))


if __name__ == '__main__':
    main()
