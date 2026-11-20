"""Independent checks of the model against values that can be worked out by hand."""
import math
import rockoon_model as rm


def test_atmosphere_sea_level_and_tropopause():
    T, p, rho = rm.atmosphere(0.0)
    assert abs(T - 288.15) < 1e-9 and abs(p - 101325.0) < 1e-6 and abs(rho - 1.2250) < 1e-3
    T, p, _ = rm.atmosphere(11019.0)          # 11 km geopotential
    assert abs(T - 216.65) < 0.02 and abs(p - 22632.0) < 15.0


def test_column_mass_equals_pressure_over_g():
    z0, dz, col = 20e3, 20.0, 0.0
    z = z0
    while z < 85e3:
        col += rm.atmosphere(z + dz / 2)[2] * dz * (rm.R_EARTH / (rm.R_EARTH + z)) ** 0
        z += dz
    assert abs(col / (rm.atmosphere(z0)[1] / rm.G0) - 1.0) < 0.02


def test_balloon_buoyancy_closes():
    V, me, mg = rm.balloon(331.27, 30e3)
    T, p, rho = rm.atmosphere(30e3)
    assert abs(rho * V - (331.27 + me + mg)) < 1e-6 * rho * V
    assert abs(mg / (331.27 + me) - 0.0748) < 2e-4


def test_vacuum_flight_matches_rocket_equation():
    m0, mp, isp = 300.0, 150.0, 300.0
    mdot = 5.0
    r = rm.fly(m0, mp, 0.4, mdot, isp, 0.0, z_kin=200e3, drag=False, coast=False)
    assert abs(r['dv_thrust'] - isp * rm.G0 * math.log(2.0)) < 0.5
    assert abs(r['v_b'] - (r['dv_thrust'] - r['loss_g'])) < 1e-3


def test_coast_obeys_energy_conservation():
    r = rm.fly(300.0, 150.0, 0.4, 50.0, 300.0, 0.0, z_kin=200e3, drag=False)
    rb = rm.R_EARTH + 200e3 + r['h_b']
    ra = 1.0 / (1.0 / rb - r['v_b'] ** 2 / (2.0 * rm.GM))
    assert abs((ra - rm.R_EARTH) - r['apogee']) / r['apogee'] < 3e-3


def test_shapley_shares_sum_to_gain():
    s = rm.shapley(301.27, 151.27, 30e3)
    assert abs(s['head'] + s['drag'] + s['nozzle'] - (s['balloon'] - s['ground'])) < 1e-6


def test_ideal_gas_nozzle_limits():
    nz = rm.ideal_gas_nozzle(1.2, 1500.0, 5e6)
    assert nz['eps'][0] < 1.2 and all(b > a for a, b in zip(nz['isp_vac'], nz['isp_vac'][1:]))


def test_tank_mass_scales_with_pressure():
    a = rm.tank_fraction(3e6, 2.4, 2700.0, 276e6)
    b = rm.tank_fraction(6e6, 2.4, 2700.0, 276e6)
    assert abs(b / a - 2.0) < 1e-12


def test_fineness_drag_vanishes_at_the_reference_fineness():
    T, p, rho = rm.atmosphere(5e3)
    a = math.sqrt(1.4 * rm.RSTAR / rm.M_AIR * T)
    same = rm.FinenessDrag(0.30, 10.0); slim = rm.FinenessDrag(0.30, 25.0)
    for M in (0.5, 2.0, 4.0):
        assert abs(same(M, rho, T, M * a) - rm.cd_generic(M)) < 1e-12
        assert slim(M, rho, T, M * a) > rm.cd_generic(M)


def test_envelope_rule_reduces_to_the_film_rule():
    import general_rule as gr
    cv = gr.Curve(130.0, 150.0, 9.5)
    film = gr.Law('area', '', [(0.030 * gr.K36, 2 / 3, 0.0)])
    z = gr.exact(cv, film)
    st = gr.state(cv, film, z)
    assert abs(st['me'] / (1.5 * st['H'] * abs(cv.dms(z))) - 1.0) < 1e-5
    assert abs(st['eps'] - 1.0 / (st['alpha'] * st['Lam'] / st['H'] - st['beta'])) < 1e-6


def test_envelope_rule_holds_for_a_two_term_law_and_a_price():
    import general_rule as gr
    cv = gr.Curve(130.0, 150.0, 9.5)
    law = gr.Law('fitted', '', [(0.0071 * gr.K36, 2 / 3, 0.0), (6.77e-3, 1 / 3, 1.0)])
    for mu in (0.4, 1.0, 3.0):
        z = gr.exact(cv, law, mu=mu)
        assert abs(gr.ratio(cv, law, z, mu=mu) - 1.0) < 1e-5
