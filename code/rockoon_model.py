"""First-order model of a balloon-launched pressure-fed sounding rocket.

Atmosphere (U.S. Standard Atmosphere 1976), zero-pressure balloon sizing,
UDMH/N2O4 equilibrium thermochemistry (Cantera), nozzle selection,
one-dimensional trajectory with switchable environments, and tank masses.
"""
import json, math, os
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

G0 = 9.80665
RSTAR = 8.31432
M_AIR = 28.9644e-3
R0_GEOPOT = 6356766.0
R_EARTH = 6371000.0
GM = 3.986004418e14
M_H2 = 2.01588e-3
M_HE = 4.002602e-3
R_HE = 2077.1
HERE = os.path.dirname(os.path.abspath(__file__))

# ---------------------------------------------------------------- atmosphere
_HB = [0.0, 11000.0, 20000.0, 32000.0, 47000.0, 51000.0, 71000.0, 84852.0]
_LB = [-6.5e-3, 0.0, 1.0e-3, 2.8e-3, 0.0, -2.8e-3, -2.0e-3, 0.0]
_TB = [288.15]
_PB = [101325.0]
for _i in range(1, 8):
    _T = _TB[-1] + _LB[_i - 1] * (_HB[_i] - _HB[_i - 1])
    if _LB[_i - 1] != 0.0:
        _p = _PB[-1] * (_TB[-1] / _T) ** (G0 * M_AIR / (RSTAR * _LB[_i - 1]))
    else:
        _p = _PB[-1] * math.exp(-G0 * M_AIR * (_HB[_i] - _HB[_i - 1]) / (RSTAR * _TB[-1]))
    _TB.append(_T)
    _PB.append(_p)


def atmosphere(z):
    """Temperature (K), pressure (Pa), density (kg/m3) at geometric altitude z (m)."""
    if z < 0.0:
        z = 0.0
    H = R0_GEOPOT * z / (R0_GEOPOT + z)
    i = 7
    for j in range(1, 8):
        if H < _HB[j]:
            i = j - 1
            break
    if i < 7 and _LB[i] != 0.0:
        T = _TB[i] + _LB[i] * (H - _HB[i])
        p = _PB[i] * (_TB[i] / T) ** (G0 * M_AIR / (RSTAR * _LB[i]))
    else:
        T = _TB[i]
        p = _PB[i] * math.exp(-G0 * M_AIR * (H - _HB[i]) / (RSTAR * T))
    return T, p, p * M_AIR / (RSTAR * T)


def lapse_extrapolation(z, lapse=6.49e-3):
    """Tropospheric lapse rate extrapolated to any altitude (model of the 2022 proposal)."""
    T = 288.15 - lapse * z
    p = 101325.0 * (T / 288.15) ** (G0 * M_AIR / (RSTAR * lapse))
    return T, p, p * M_AIR / (RSTAR * T)


# ------------------------------------------------------------------- balloon
def balloon(m_susp, z, M_gas=M_H2, sigma=0.030, k_load=0.0):
    """Zero-pressure balloon at float: volume, envelope mass, gas mass.

    Envelope mass = sigma * (area of the equal-volume sphere) + k_load * m_susp; the second term
    stands for load tapes and film that are sized by the suspended load rather than by the area."""
    T, p, _ = atmosphere(z)
    dens = p * (M_AIR - M_gas) / (RSTAR * T)
    k = sigma * (36.0 * math.pi) ** (1.0 / 3.0)
    V = brentq(lambda V: V * dens - (1.0 + k_load) * m_susp - k * V ** (2.0 / 3.0), 1e-3, 1e9)
    m_env = k * V ** (2.0 / 3.0) + k_load * m_susp
    return V, m_env, (m_susp + m_env) * M_gas / (M_AIR - M_gas)


def insulation_mass(mp, z_float, of=2.4, rate=5.0, hold=1800.0, T0=293.15, dT=25.0, c=1500.0,
                    k_ins=0.030, rho_ins=35.0, t_min=0.025):
    """Mass of a foam jacket that keeps the propellant above freezing during ascent and float.

    N2O4 freezes at about -11 C, so the propellant may cool by dT from T0. The jacket covers two
    spherical tanks, conducts heat to the local air temperature (no outer film resistance, no
    radiation, propellant held at T0: all conservative), and is released with the gondola at ignition.
    Returns (mass, thickness, heat-leak integral in K s)."""
    vf = (1.0 / (1.0 + of)) / RHO_UDMH * mp
    vo = (of / (1.0 + of)) / RHO_NTO * mp
    A = sum((36.0 * math.pi) ** (1.0 / 3.0) * v ** (2.0 / 3.0) for v in (vf, vo))
    zs = np.linspace(0.0, z_float, 400)
    Ta = np.array([atmosphere(z)[0] for z in zs])
    I = np.trapezoid(T0 - Ta, zs) / rate + (T0 - Ta[-1]) * hold
    t = max(k_ins * A * I / (mp * c * dT), t_min)
    return rho_ins * A * t, t, I


# ------------------------------------------------------------ thermochemistry
M_UDMH = 60.0983e-3
M_NTO = 92.011e-3
HF_UDMH = 48.3e3     # J/mol, liquid
HF_NTO = -19.56e3    # J/mol, liquid
RHO_UDMH = 791.0
RHO_NTO = 1450.0
_GAS = None


def _gas():
    global _GAS
    if _GAS is None:
        import cantera as ct
        sp = [s for s in ct.Species.list_from_file('nasa_gas.yaml')
              if set(s.composition) <= {'C', 'H', 'N', 'O'} and s.composition.get('C', 0) <= 3]
        _GAS = ct.Solution(thermo='ideal-gas', species=sp)
    return _GAS


def chamber(of, pc):
    gas = _gas()
    nU = 1.0 / (1.0 + of) / M_UDMH
    nN = of / (1.0 + of) / M_NTO
    h = nU * HF_UDMH + nN * HF_NTO
    gas.TPX = 2500.0, pc, {'CH4': 2 * nU, 'N2': nU + nN, 'O2': 2 * nN}
    gas.HP = h, pc
    gas.equilibrate('HP', solver='gibbs', max_steps=5000)
    return h


def expansion(of, pc, shifting=True, n=160, pr_min=2.0e-5):
    """Isentropic nozzle expansion. Returns dict of arrays versus exit pressure."""
    gas = _gas()
    hc = chamber(of, pc)
    Tc, s, Xc = gas.T, gas.entropy_mass, gas.X.copy()
    ps = pc * np.exp(np.linspace(math.log(0.75), math.log(pr_min), n))
    u, rho = [], []
    for p in ps:
        if shifting:
            gas.SPX = s, p, gas.X
            gas.equilibrate('SP', solver='gibbs', max_steps=5000)
        else:
            gas.SPX = s, p, Xc
        u.append(math.sqrt(max(2.0 * (hc - gas.enthalpy_mass), 0.0)))
        rho.append(gas.density)
    u, rho = np.array(u), np.array(rho)
    flux = rho * u
    it = int(np.argmax(flux))
    # refine the throat with a parabola through the three points around the maximum
    x = np.log(ps[it - 1:it + 2]); y = flux[it - 1:it + 2]
    a = np.polyfit(x, y, 2)
    Gt = float(np.polyval(a, -a[1] / (2 * a[0])))
    sl = slice(it + 1, None)
    eps = Gt / flux[sl]
    isp_vac = (u[sl] + ps[sl] / flux[sl]) / G0
    isp_opt = u[sl] / G0
    return dict(of=of, pc=pc, Tc=Tc, cstar=pc / Gt, pe=ps[sl].tolist(), eps=eps.tolist(),
                isp_vac=isp_vac.tolist(), isp_opt=isp_opt.tolist())


def build_tables(path):
    out = {'nozzle': expansion(2.4, 2.0e6, True), 'nozzle_1MPa': expansion(2.4, 1.0e6, True),
           'nozzle_4MPa': expansion(2.4, 4.0e6, True)}
    ofs = [0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.05, 2.2, 2.4, 2.6, 2.75, 3.0, 3.3, 3.6, 4.0]
    sweep = []
    for of in ofs:
        a = expansion(of, 2.0e6, True); b = expansion(of, 2.0e6, False)
        sweep.append(dict(of=of, Tc=a['Tc'], cstar=a['cstar'],
                          shift=float(np.interp(36.4, a['eps'], a['isp_vac'])),
                          frozen=float(np.interp(36.4, b['eps'], b['isp_vac']))))
    out['of_sweep'] = sweep
    a = expansion(2 / 3, 2.0e6, True); b = expansion(2 / 3, 2.0e6, False)
    out['fuel_rich'] = dict(Tc=a['Tc'], shift=float(np.interp(36.4, a['eps'], a['isp_vac'])),
                            frozen=float(np.interp(36.4, b['eps'], b['isp_vac'])))
    pc = 68 * 101325.0
    bench = {}
    for of in (2.1, 2.6):
        a = expansion(of, pc, True); b = expansion(of, pc, False)
        bench[str(of)] = dict(shift=float(np.interp(-101325.0, -np.array(a['pe']), a['isp_opt'])),
                              frozen=float(np.interp(-101325.0, -np.array(b['pe']), b['isp_opt'])),
                              Tc=a['Tc'], cstar=a['cstar'])
    out['benchmark_68atm'] = bench
    with open(path, 'w') as f:
        json.dump(out, f, indent=1)
    return out


def ideal_gas_nozzle(gamma, cstar, pc, n=220, pr_min=2.0e-5):
    """Nozzle table for a constant-gamma gas (used for the generic solid-motor family)."""
    pr = np.exp(np.linspace(math.log(0.5), math.log(pr_min), n))
    g = gamma
    Me = np.sqrt(2.0 / (g - 1.0) * (pr ** (-(g - 1.0) / g) - 1.0))
    eps = (1.0 / Me) * ((2.0 / (g + 1.0)) * (1.0 + 0.5 * (g - 1.0) * Me ** 2)) ** ((g + 1.0) / (2.0 * (g - 1.0)))
    cf0 = np.sqrt(2.0 * g * g / (g - 1.0) * (2.0 / (g + 1.0)) ** ((g + 1.0) / (g - 1.0)) * (1.0 - pr ** ((g - 1.0) / g)))
    keep = Me > 1.0
    return dict(pc=pc, cstar=cstar, Tc=float('nan'), pe=(pr * pc)[keep], eps=eps[keep],
                isp_vac=(cstar * (cf0 + pr * eps) / G0)[keep], isp_opt=(cstar * cf0 / G0)[keep])


SOLID = dict(gamma=1.16, cstar=1500.0, pc=6.0e6)
_TABLE = None


def table():
    global _TABLE
    if _TABLE is None:
        path = os.path.join(HERE, 'thermo_table.json')
        if not os.path.exists(path):
            build_tables(path)
        with open(path) as f:
            _TABLE = json.load(f)
        for key in ('nozzle', 'nozzle_1MPa', 'nozzle_4MPa'):
            nz = _TABLE[key]
            for k in ('pe', 'eps', 'isp_vac', 'isp_opt'):
                nz[k] = np.array(nz[k])
        _TABLE['solid'] = ideal_gas_nozzle(**SOLID)
    return _TABLE


class Nozzle:
    """Nozzle for the baseline propellant state (O/F 2.4, 2.0 MPa), chosen for an ignition altitude.

    The expansion ratio is the largest one for which the exit pressure is not below `sep` times the
    ambient pressure at ignition (sep = 1: no over-expansion; sep = 0.4: over-expanded to a typical
    separation limit) and the exit area does not exceed the body cross-section.
    """

    def __init__(self, z_ign, thrust, d_body, eta=0.92, key='nozzle', sep=1.0):
        nz = table()[key]
        pa = atmosphere(z_ign)[1]
        A_body = math.pi * d_body ** 2 / 4.0
        cstar, pc = nz['cstar'], nz['pc']
        best = None
        for eps, pe, isp in zip(nz['eps'], nz['pe'], nz['isp_vac']):
            if pe < sep * pa:
                break
            mdot = thrust / (G0 * eta * isp - pa * eps * cstar / pc)
            Ae = eps * mdot * cstar / pc
            if Ae > A_body:
                break
            best = (eps, pe, isp, mdot, Ae)
        if best is None:
            raise ValueError('no admissible nozzle')
        self.eps, self.pe, isp, self.mdot, self.Ae = best
        self.isp_vac = eta * isp
        self.limit = 'pressure' if self.pe < 1.05 * pa else 'area'


# ---------------------------------------------------------------- trajectory
CD_SUB, CD_PEAK, CD_INF, CD_EXP = 0.30, 0.55, 0.25, 1.2


def cd_generic(M):
    if M <= 0.8:
        return CD_SUB
    if M <= 1.05:
        return CD_SUB + (CD_PEAK - CD_SUB) * (M - 0.8) / 0.25
    return CD_INF + (CD_PEAK - CD_INF) * (1.05 / M) ** CD_EXP


_SL_M = [0, 0.50, 0.80, 0.81, 1.40, 1.50, 1.75, 2.00, 2.25, 2.75, 3.00, 3.25, 3.50, 3.75, 4.00, 4.25, 4.50, 4.75, 5.00,
         5.25, 5.50, 5.75, 6.00, 7.00]
_SL_C = [0.400, 0.420, 0.430, 0.930, 0.930, 0.880, 0.770, 0.682, 0.620, 0.533, 0.500, 0.470, 0.440, 0.420, 0.397, 0.380,
         0.360, 0.345, 0.331, 0.320, 0.310, 0.300, 0.290, 0.260]
_DART_M = [0, 0.90, 1.00, 1.50, 1.75, 2.00, 2.25, 2.50, 2.75, 3.00, 3.50, 4.00, 4.50, 5.00, 6.00, 7.00, 8.00]
_DART_C = [0.350, 0.350, 0.576, 0.576, 0.489, 0.432, 0.393, 0.365, 0.345, 0.329, 0.308, 0.294, 0.284, 0.277, 0.277,
           0.262, 0.262]


def cd_superloki(M):
    """Booster plus dart, NASA CR-61238 Table 4.1 (reference area 85.47 cm2; the Mach 2.50 entry is not used)."""
    return float(np.interp(M, _SL_M, _SL_C))


def cd_dart(M):
    """Dart alone, NASA CR-61238 Table 4.1 (reference area 9.84 cm2)."""
    return float(np.interp(M, _DART_M, _DART_C))


def fly2d(m0, mp, area, thrust, tb, cd_fn, elev_deg, z0=0.0, m2=None, area2=None, cd_fn2=None, cd_scale=1.0):
    """Planar point-mass flight with thrust and drag along the velocity (gravity turn), constant
    thrust for tb seconds, and an optional change of mass, area and drag curve at burnout."""
    th = math.radians(elev_deg)
    mdot = mp / tb

    def rhs(t, y, burning, A, cdf):
        x, z, vx, vz, m = y
        v = math.hypot(vx, vz)
        T, _, rho = atmosphere(z)
        g = G0 * (R_EARTH / (R_EARTH + z)) ** 2
        D = 0.5 * rho * v * v * cd_scale * cdf(v / math.sqrt(1.4 * RSTAR / M_AIR * T)) * A
        a = ((thrust if burning else 0.0) - D) / m
        return [vx, vz, a * vx / v, a * vz / v - g, -mdot if burning else 0.0]

    y0 = [0.0, z0, math.cos(th), math.sin(th), m0]
    s1 = solve_ivp(lambda t, y: rhs(t, y, True, area, cd_fn), [0.0, tb], y0, rtol=1e-9, atol=1e-9, dense_output=True)
    yb = s1.y[:, -1].copy()
    out = dict(v_b=math.hypot(yb[2], yb[3]), z_b=yb[1], sol=s1.sol)
    if m2 is not None:
        yb[4] = m2
    ev = lambda t, y: y[3]
    ev.terminal = True
    ev.direction = -1
    s2 = solve_ivp(lambda t, y: rhs(t, y, False, area2 or area, cd_fn2 or cd_fn), [tb, tb + 3000.0], yb, rtol=1e-9, atol=1e-9, events=ev)
    out.update(apogee=s2.y[1, -1], t_apogee=s2.t[-1], range=s2.y[0, -1])
    return out


# ---------------------------------------------------- component drag build-up
class DragBuildUp:
    """Zero-lift drag by the component build-up of the OpenRocket technical documentation
    (Niskanen, version 13.05): skin friction Eqs. 3.78-3.85, tangent-ogive nose pressure drag
    from the conical profile (Appendix B), fin leading- and trailing-edge drag Eqs. 3.89-3.93,
    and base drag Eq. 3.94. The method is stated for Mach numbers up to about 5.
    All coefficients are referred to the body cross-section."""
    full = True

    def __init__(self, d=0.40, length=4.0, nose=1.2, n_fins=4, span=0.40, c_root=0.60, c_tip=0.30,
                 thick=0.010, sweep_deg=30.0, rough=20e-6):
        self.d, self.L, self.rough = d, length, rough
        self.Aref = math.pi * d * d / 4.0
        R = d / 2.0
        rho_o = (R * R + nose * nose) / (2.0 * R)          # ogive radius
        xs = np.linspace(0.0, nose, 2001)
        y = np.sqrt(rho_o ** 2 - (nose - xs) ** 2) + R - rho_o
        dy = np.gradient(y, xs)
        a_nose = float(np.trapezoid(2.0 * math.pi * y * np.sqrt(1.0 + dy ** 2), xs))
        self.Awb = a_nose + math.pi * d * (length - nose)
        self.Awf = 2.0 * n_fins * 0.5 * (c_root + c_tip) * span
        lam = c_tip / c_root
        self.mac = (2.0 / 3.0) * c_root * (1.0 + lam + lam * lam) / (1.0 + lam)
        self.fB = length / d
        self.t = thick
        self.fin_front = n_fins * thick * span
        self.cos2 = math.cos(math.radians(sweep_deg)) ** 2
        self.sin_eps = math.sin(math.atan(1.0 / (2.0 * nose / d)))

    def friction(self, M, rho, T, v):
        mu = 1.458e-6 * T ** 1.5 / (T + 110.4)             # Sutherland
        Re = max(rho * v * self.L / mu, 1.0)
        rel = self.rough / self.L
        Rcrit = 51.0 * rel ** (-1.039)
        if Re < 1e4:
            cf, rough = 1.48e-2, False
        elif Re < Rcrit:
            cf, rough = 1.0 / (1.50 * math.log(Re) - 5.6) ** 2, False
        else:
            cf, rough = 0.032 * rel ** 0.2, True
        if M < 1.0:
            cfc = cf * (1.0 - 0.1 * M * M)
        else:
            turb = (1.0 / (1.50 * math.log(max(Re, 1e4)) - 5.6) ** 2) / (1.0 + 0.15 * M * M) ** 0.58
            cfc = max(cf / (1.0 + 0.18 * M * M), turb) if rough else cf / (1.0 + 0.15 * M * M) ** 0.58
        return cfc * ((1.0 + 1.0 / (2.0 * self.fB)) * self.Awb + (1.0 + 2.0 * self.t / self.mac) * self.Awf) / self.Aref

    def nose(self, M):
        se = self.sin_eps
        c1 = se
        s1 = (4.0 / 2.4) * (1.0 - 0.5 * c1)
        f = lambda m: 2.1 * se * se + 0.5 * se / math.sqrt(m * m - 1.0)
        if M >= 1.3:
            return f(M)
        if M >= 1.0:                                        # cubic through values and slopes at 1 and 1.3
            c2 = f(1.3)
            s2 = -0.5 * se * 1.3 / (1.3 * 1.3 - 1.0) ** 1.5
            h = 0.3
            u = (M - 1.0) / h
            return ((2 * u ** 3 - 3 * u ** 2 + 1) * c1 + (u ** 3 - 2 * u ** 2 + u) * h * s1
                    + (-2 * u ** 3 + 3 * u ** 2) * c2 + (u ** 3 - u ** 2) * h * s2)
        return c1 * M ** (s1 / c1)                          # tangent ogive: zero at M = 0

    @staticmethod
    def base(M):
        return 0.12 + 0.13 * M * M if M < 1.0 else 0.25 / M

    def fins(self, M):
        if M < 0.9:
            le = (1.0 - M * M) ** (-0.417) - 1.0
        elif M < 1.0:
            le = 1.0 - 1.785 * (M - 0.9)
        else:
            le = 1.214 - 0.502 / M ** 2 + 0.1095 / M ** 4
        return (le * self.cos2 + self.base(M)) * self.fin_front / self.Aref

    def parts(self, M, rho, T, v, base_frac=1.0):
        return dict(friction=self.friction(M, rho, T, v), nose=self.nose(M), fins=self.fins(M), base=self.base(M) * base_frac)

    def __call__(self, M, rho, T, v, base_frac=1.0):
        return sum(self.parts(M, rho, T, v, base_frac).values())


def fly_tilted(m0, mp, d, mdot, isp_vac, Ae, z0, tilt_deg, mode='fixed', rail=10.0, cd_scale=1.0):
    """Planar flight of the liquid vehicle ignited at an angle tilt_deg from the vertical.

    mode 'fixed': the thrust keeps its initial direction (a spin-stabilised or attitude-held vehicle,
    the relevant case where there is no air to turn the body). mode 'turn': the thrust follows the
    velocity after the first `rail` metres of travel (gravity turn, aerodynamically stable vehicle).
    """
    A = math.pi * d * d / 4.0
    tb = mp / mdot
    th = math.radians(90.0 - tilt_deg)
    ux, uz = math.cos(th), math.sin(th)

    def rhs(t, y, burning):
        x, z, vx, vz, m, s = y
        v = math.hypot(vx, vz)
        T, p, rho = atmosphere(z)
        g = G0 * (R_EARTH / (R_EARTH + z)) ** 2
        D = 0.5 * rho * v * v * cd_scale * cd_generic(v / math.sqrt(1.4 * RSTAR / M_AIR * T)) * A if v > 0 else 0.0
        F = (mdot * G0 * isp_vac - p * Ae) if burning else 0.0
        if mode == 'fixed' or s < rail or v < 1e-6:
            tx, tz = ux, uz
        else:
            tx, tz = vx / v, vz / v
        dx, dz = (vx / v, vz / v) if v > 1e-9 else (0.0, 0.0)
        ax = (F * tx - D * dx) / m
        az = (F * tz - D * dz) / m - g
        if s < rail and mode == 'turn':            # on the rail: motion constrained to the rail direction
            a_r = ax * ux + az * uz
            ax, az = a_r * ux, a_r * uz
        return [vx, vz, ax, az, -mdot if burning else 0.0, v]

    s1 = solve_ivp(lambda t, y: rhs(t, y, True), [0.0, tb], [0.0, z0, 0.0, 0.0, m0, 0.0], rtol=1e-8, atol=1e-8, max_step=0.5)
    ev = lambda t, y: y[3]
    ev.terminal = True
    ev.direction = -1
    s2 = solve_ivp(lambda t, y: rhs(t, y, False), [tb, tb + 5000.0], s1.y[:, -1], rtol=1e-8, atol=1e-8, events=ev)
    return dict(apogee=s2.y[1, -1], range=s2.y[0, -1], v_b=math.hypot(s1.y[2, -1], s1.y[3, -1]))


def fly(m0, mp, d, mdot, isp_vac, Ae, z_kin=0.0, z_drag=None, z_noz=None, v0=0.0,
        cd_scale=1.0, elev_deg=90.0, drag=True, coast=True, trace=False, cd_fn=None, rtol=1e-8):
    """Point-mass flight along a straight line of elevation elev_deg.

    Altitude is z_kin + h. Drag is evaluated in the atmosphere at z_drag + h and the nozzle
    back-pressure at z_noz + h, so that the three effects of a high start can be switched
    independently. With the defaults all three offsets coincide.
    """
    z_drag = z_kin if z_drag is None else z_drag
    z_noz = z_kin if z_noz is None else z_noz
    A = math.pi * d * d / 4.0
    tb = mp / mdot
    se = math.sin(math.radians(elev_deg))
    Fv = mdot * G0 * isp_vac

    def rhs(t, y, burning):
        h, v, m = y[0], y[1], y[2]
        g = G0 * (R_EARTH / (R_EARTH + z_kin + h)) ** 2 * se
        T, _, rho = atmosphere(z_drag + h)
        D = 0.0
        if drag and v != 0.0:
            Mach = abs(v) / math.sqrt(1.4 * RSTAR / M_AIR * T)
            if getattr(cd_fn, 'full', False):
                cdv = cd_fn(Mach, rho, T, abs(v), max(1.0 - Ae / A, 0.0) if burning else 1.0)
            else:
                cdv = (cd_fn or cd_generic)(Mach)
            D = 0.5 * rho * v * abs(v) * cd_scale * cdv * A
        if burning:
            F = Fv - atmosphere(z_noz + h)[1] * Ae
            return [v * se, (F - D) / m - g, -mdot, g, D / m, F / m, 0.5 * rho * v * v, 0.0]
        return [v * se, -D / m - g, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    y0 = [0.0, v0, m0, 0.0, 0.0, 0.0, 0.0, 0.0]
    s1 = solve_ivp(lambda t, y: rhs(t, y, True), [0.0, tb], y0, rtol=rtol, atol=rtol,
                   dense_output=True)
    yb = s1.y[:, -1]
    ts = np.linspace(0.0, tb, 400)
    ys = s1.sol(ts)
    Ts_rho = [atmosphere(z_drag + hh)[2] for hh in ys[0]]
    q = 0.5 * np.array(Ts_rho) * ys[1] ** 2
    out = dict(tb=tb, h_b=yb[0], v_b=yb[1], v_max=float(np.max(ys[1])), t_vmax=float(ts[int(np.argmax(ys[1]))]),
               dv_thrust=yb[5], loss_g=yb[3], loss_d=yb[4], q_max=float(np.max(q)),
               F0=Fv - atmosphere(z_noz)[1] * Ae, isp_ign=(Fv - atmosphere(z_noz)[1] * Ae) / (mdot * G0))
    if coast:
        ev = lambda t, y: y[1]
        ev.terminal = True
        ev.direction = -1
        s2 = solve_ivp(lambda t, y: rhs(t, y, False), [tb, tb + 5000.0], yb, rtol=rtol, atol=rtol,
                       events=ev, dense_output=True)
        out['apogee'] = z_kin + s2.y[0, -1]
        out['t_apogee'] = s2.t[-1]
        if trace:
            tc = np.linspace(tb, s2.t[-1], 300)
            yc = s2.sol(tc)
            out['trace'] = (np.concatenate([ts, tc]), z_kin + np.concatenate([ys[0], yc[0]]),
                            np.concatenate([ys[1], yc[1]]))
    return out


def diameter(m0, d_ref=0.40, m_ref=301.27, fineness=10.0):
    """Body diameter at constant bulk density. The baseline (0.40 m) is taken to have a
    length-to-diameter ratio of 10; other values rescale the diameter at constant volume."""
    return d_ref * (m0 / m_ref) ** (1.0 / 3.0) * (10.0 / fineness) ** (1.0 / 3.0)


def vehicle_run(m0, mp, z_ign, d=None, tw=4.0, eta=0.92, cd_scale=1.0, switches=None, v0=0.0, key='nozzle', cd_fn=None, sep=1.0, fineness=10.0):
    """Design the nozzle for the ignition altitude and fly. switches = (kin, drag, nozzle) booleans
    selecting the balloon (True) or ground (False) environment for each effect."""
    d = diameter(m0, fineness=fineness) if d is None else d
    F = tw * m0 * G0
    if switches is None:
        nz = Nozzle(z_ign, F, d, eta, key, sep)
        r = fly(m0, mp, d, nz.mdot, nz.isp_vac, nz.Ae, z_kin=z_ign, cd_scale=cd_scale, v0=v0, cd_fn=cd_fn)
    else:
        k, dr, n = switches
        nz = Nozzle(z_ign if n else 0.0, F, d, eta, key, sep)
        r = fly(m0, mp, d, nz.mdot, nz.isp_vac, nz.Ae, z_kin=z_ign if k else 0.0,
                z_drag=z_ign if dr else 0.0, z_noz=z_ign if n else 0.0, cd_scale=cd_scale, cd_fn=cd_fn)
    r.update(eps=nz.eps, pe=nz.pe, Ae=nz.Ae, mdot=nz.mdot, isp_vac=nz.isp_vac, limit=nz.limit, d=d)
    return r


def shapley(m0, mp, z_f, **kw):
    """Order-averaged attribution of the apogee gain to head start, drag and nozzle."""
    ap = {}
    for k in (0, 1):
        for dr in (0, 1):
            for n in (0, 1):
                ap[(k, dr, n)] = vehicle_run(m0, mp, z_f, switches=(k, dr, n), **kw)['apogee']
    phi = [0.0, 0.0, 0.0]
    w = {0: 1 / 3, 1: 1 / 6, 2: 1 / 3}
    for i in range(3):
        for s in ap:
            if s[i] == 0:
                t = list(s); t[i] = 1
                phi[i] += w[sum(s)] * (ap[tuple(t)] - ap[s])
    return dict(ground=ap[(0, 0, 0)], balloon=ap[(1, 1, 1)], head=phi[0], drag=phi[1], nozzle=phi[2], runs=ap)


# --------------------------------------------------------------------- tanks
def prop_volume(of):
    return (1.0 / (1.0 + of)) / RHO_UDMH + (of / (1.0 + of)) / RHO_NTO


def tank_fraction(p, of, rho_m, sigma_y, sf=1.5, ullage=0.05):
    """Mass of ideal thin-walled spherical tanks per kg of propellant."""
    return 1.5 * p * (1.0 + ullage) * prop_volume(of) * rho_m / (sigma_y / sf)


def pressurant(p_tank, of, p_bottle, p_end, rho_m, sigma_y, sf=1.5, ullage=0.05, T0=290.0):
    """Helium and bottle mass per kg of propellant; adiabatic blow-down of an ideal gas."""
    Vt = (1.0 + ullage) * prop_volume(of)
    Vb = (5.0 / 3.0) * p_tank * Vt / (p_bottle - p_end)
    m_he = p_bottle * Vb / (R_HE * T0)
    m_b = 1.5 * p_bottle * Vb * rho_m / (sigma_y / sf)
    return m_he, m_b, Vb


MATERIALS = {'SS321': (8000.0, 205e6), 'Al6061T6': (2700.0, 276e6)}


def feed_fraction(material, p_tank=3.0e6, of=2.4, geom=1.5, p_bottle=30e6):
    rho_m, sy = MATERIALS[material]
    ft = tank_fraction(p_tank, of, rho_m, sy)
    m_he, m_b, _ = pressurant(p_tank, of, p_bottle, 1.2 * p_tank, 2700.0, 276e6)
    return geom * ft + m_he + m_b


def nozzle_skirt(Ae, eps, half_angle=15.0, areal=10.0):
    """Mass (kg) and length (m) of a conical divergent section from throat to exit.

    areal is the wall mass per unit area in kg/m^2 (10 kg/m^2 is about 1.3 mm of steel); it is an
    assumption, varied in the paper."""
    re_ = math.sqrt(Ae / math.pi)
    rt = re_ / math.sqrt(eps)
    a = math.radians(half_angle)
    L = (re_ - rt) / math.tan(a)
    return areal * math.pi * (re_ + rt) * L / math.cos(a), L


def close_vehicle_skirt(z_f, m_fix, f_feed, target=100e3, areal=10.0, **kw):
    """As close_vehicle, with the mass of the divergent nozzle section added to the fixed mass."""
    extra = 0.0
    r = None
    for _ in range(6):
        r = close_vehicle(z_f, m_fix + extra, f_feed, target=target, **kw)
        if r is None:
            return None
        new = nozzle_skirt(r['Ae'], r['eps'], areal=areal)[0]
        if abs(new - extra) < 0.01:
            break
        extra = new
    r['skirt'] = extra
    return r


def close_vehicle(z_f, m_fix, f_feed, target=100e3, **kw):
    """Propellant mass for which the vehicle reaches the target apogee from z_f."""
    def miss(mp):
        m0 = m_fix + (1.0 + f_feed) * mp
        return vehicle_run(m0, mp, z_f, **kw)['apogee'] - target
    lo, hi = 2.0, 20.0
    if miss(lo) > 0:
        return None
    while miss(hi) < 0:
        hi *= 1.6
        if hi > 5e4:
            return None
    mp = brentq(miss, lo, hi, xtol=1e-3)
    m0 = m_fix + (1.0 + f_feed) * mp
    r = vehicle_run(m0, mp, z_f, **kw)
    r.update(mp=mp, m0=m0, m_dry=m0 - mp)
    return r


if __name__ == '__main__':
    build_tables(os.path.join(HERE, 'thermo_table.json'))
