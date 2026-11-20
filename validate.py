#!/usr/bin/env python3
"""One-command validity check for the rockoon paper.

    python validate.py            quick checks (about 20 seconds)
    python validate.py --full     also regenerates every result and compares it with the committed files

Exit status is 0 only if every check passes.
"""
import hashlib, json, math, os, shutil, subprocess, sys, tempfile

ROOT = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.join(ROOT, 'code')
sys.path.insert(0, CODE)
import rockoon_model as rm  # noqa: E402

ROWS = []


def check(group, name, got, want, tol, unit='', rel=True):
    err = abs(got / want - 1.0) if rel else abs(got - want)
    ok = err <= tol
    ROWS.append((group, name, '%.4g' % got, '%.4g' % want, ('%.2f%%' % (100 * err)) if rel else ('%.3g' % err),
                 ('%.1f%%' % (100 * tol)) if rel else ('%.3g' % tol), unit, ok))
    return ok


def flag(group, name, ok, detail=''):
    ROWS.append((group, name, detail, '', '', '', '', bool(ok)))
    return ok


def physics():
    g = 'Physics identities'
    T, p, rho = rm.atmosphere(0.0)
    check(g, 'sea-level density', rho, 1.2250, 1e-3, 'kg/m3')
    T, p, rho = rm.atmosphere(30e3)
    check(g, 'pressure at 30 km (USSA-76)', p, 1197.0, 2e-3, 'Pa')
    check(g, 'temperature at 30 km (USSA-76)', T, 226.51, 1e-3, 'K')
    z, col, dz = 20e3, 0.0, 20.0
    while z < 85e3:
        col += rm.atmosphere(z + dz / 2)[2] * dz
        z += dz
    check(g, 'column mass above 20 km = p/g', col, rm.atmosphere(20e3)[1] / rm.G0, 0.02, 'kg/m2')
    V, me, mg = rm.balloon(331.27, 30e3)
    check(g, 'buoyancy balance closes', rm.atmosphere(30e3)[2] * V, 331.27 + me + mg, 1e-6, 'kg')
    r = rm.fly(300.0, 150.0, 0.4, 5.0, 300.0, 0.0, z_kin=200e3, drag=False, coast=False)
    check(g, 'vacuum burn = rocket equation', r['dv_thrust'], 300.0 * rm.G0 * math.log(2.0), 1e-4, 'm/s')
    r = rm.fly(300.0, 150.0, 0.4, 50.0, 300.0, 0.0, z_kin=200e3, drag=False)
    rb = rm.R_EARTH + 200e3 + r['h_b']
    check(g, 'coast conserves energy', r['apogee'], 1.0 / (1.0 / rb - r['v_b'] ** 2 / (2.0 * rm.GM)) - rm.R_EARTH, 3e-3, 'm')
    s = rm.shapley(301.27, 151.27, 30e3)
    check(g, 'attribution sums to the gain', s['head'] + s['drag'] + s['nozzle'], s['balloon'] - s['ground'], 1e-9, 'm')


def external():
    g = 'External references'
    t = rm.table()
    b = t['benchmark_68atm']['2.6']
    check(g, 'UDMH/N2O4 Isp, 68 atm -> 1 atm, shifting (compiled value 285 s)', b['shift'], 285.0, 0.01, 's')
    check(g, 'UDMH/N2O4 c* (compiled value 1720 m/s)', b['cstar'], 1720.0, 0.01, 'm/s')
    LB, FT, LBF = 0.45359237, 0.3048, 4.4482216
    mp_d, d_d = 103.0 * LB, 6.75 * 0.0254
    rounds = [(2, 232.65, 3307, 20008), (3, 215.30, 3571, 19800), (4, 213.40, 3747, 20168), (5, 197.45, 4094, 19912),
              (6, 190.50, 4180, 19512), (7, 196.51, 4122, 19904)]
    for n, w, v, I in rounds:
        r = rm.fly(w * LB, mp_d, d_d, mp_d / 2.9, I * LBF / (mp_d * rm.G0), 0.0, elev_deg=70.0, coast=False)
        check(g, 'Deacon round %d maximum velocity (NACA RM L8H26)' % n, r['v_max'], v * FT, 0.03, 'm/s')
    for nm, zmi, obs in (('SUI-3', 7.2, 37.9), ('SUI-4', 6.6, 36.9), ('SUI-7', 8.1, 39.8)):
        rr = rm.fly(105.0, mp_d, d_d, mp_d / 2.9, 19680 * LBF / (mp_d * rm.G0), 0.0, z_kin=zmi * 1609.344)
        check(g, 'Deacon rockoon %s zenith, 105 kg (NRL Report 4276)' % nm, rr['apogee'], obs * 1609.344, 0.05, 'm')
    r = rm.fly2d(28.32, 16.87, 85.47e-4, 4450 * LBF, 2.0, rm.cd_superloki, 80.0, m2=6.12, area2=9.84e-4, cd_fn2=rm.cd_dart)
    check(g, 'Super Loki burnout velocity (NASA CR-61238)', r['v_b'], 5893 * FT, 0.05, 'm/s')
    check(g, 'Super Loki apogee (NASA CR-61238 trajectory)', r['apogee'], 113.4e3, 0.08, 'm')
    C = json.load(open(os.path.join(ROOT, 'results', 'cea_check.json')))
    check(g, 'vacuum Isp at the 30 km area ratio vs NASA CEA run directly', C['cantera']['isp_vac_b'], C['design_shift']['isp_vac'][3], 0.005, 's')
    check(g, 'vacuum Isp at the ground area ratio vs NASA CEA run directly', C['cantera']['isp_vac_g'], C['design_shift']['isp_vac'][2], 0.005, 's')
    check(g, 'c* vs NASA CEA run directly', C['cantera']['cstar'], C['design_shift']['cstar'], 0.005, 'm/s')
    P = json.load(open(os.path.join(ROOT, 'results', 'rocketpy_check.json')))
    for k in ('ground_3DOF', 'ground_6DOF', 'balloon_3DOF', 'balloon_6DOF'):
        check(g, 'baseline apogee vs RocketPy (%s)' % k.replace('_', ', '), P[k]['ours_apogee'], P[k]['apogee'], 0.002, 'm')
    A = json.load(open(os.path.join(ROOT, 'results', 'aerobee_check.json')))
    for r in A['rows']:
        if r['nominal']:
            check(g, 'Aerobee 150 nominal peak altitude (%s)' % r['label'], r['model_peak'], r['peak'], 0.05, 'm')
            check(g, 'Aerobee 150 nominal burnout velocity (%s)' % r['label'], r['model_vb'], r['vb'], 0.03, 'm/s')
    flag(g, 'Aerobee 150 engine: delivered/ideal vacuum specific impulse inside 0.85-0.95', 0.85 <= A['eta_vac'] <= 0.95, '%.3f' % A['eta_vac'])


def rule():
    g = 'Envelope rule (general_rule.py)'
    import general_rule as gr
    R = json.load(open(os.path.join(ROOT, 'results', 'results.json')))
    fr = [r for r in R['revision']['analytic'] if r['mat'] == 'Al6061T6' and r['mfix'] == 80.0][0]
    cv = gr.Curve(fr['minf'], fr['B'], fr['L'])
    film = gr.Law('area', '', [(0.030 * gr.K36, 2 / 3, 0.0)])
    zx = gr.exact(cv, film)
    check(g, 'film law: exact condition equals direct minimisation (two sides of Eq. general)', gr.ratio(cv, film, zx), 1.0, 1e-5)
    check(g, 'film law: optimum on the fitted curve equals the tabulated rule optimum', zx, fr['z_rule'], 2e-3, 'km')
    st = gr.state(cv, film, zx)
    check(g, 'film law: envelope = (3/2) H |dm_s/dz| at the optimum', st['me'], 1.5 * st['H'] * abs(cv.dms(zx)), 1e-5, 'kg')
    for mu in (0.3, 3.0):
        zm = gr.exact(cv, film, mu=mu)
        check(g, 'weighted objective, mu = %g: two sides of the condition' % mu, gr.ratio(cv, film, zm, mu=mu), 1.0, 1e-5)
    clu = gr.Law('cluster', '', [(R['general']['mu1'], 1.0, 0.0)], lift=gr.LIFT_LATEX)
    z5 = gr.exact(gr.Curve(fr['minf'], fr['B'], fr['L'], g=0.05), clu); z20 = gr.exact(gr.Curve(fr['minf'], fr['B'], fr['L'], g=0.20), clu)
    check(g, 'cluster law (alpha + beta = 1): optimum independent of the gondola fraction', z20 - z5 + 1.0, 1.0, 1e-4, 'km')
    z1, _ = gr.one_step(cv, film, 20.0)
    check(g, 'one-step estimate from 20 km against the exact optimum', z1, zx, 0.02, 'km')


def headline():
    g = 'Headline results vs committed results.json'
    R = json.load(open(os.path.join(ROOT, 'results', 'results.json')))
    b = rm.vehicle_run(301.27, 151.27, 30e3)
    gr = rm.vehicle_run(301.27, 151.27, 0.0)
    s = rm.shapley(301.27, 151.27, 30e3)
    tot = s['balloon'] - s['ground']
    check(g, 'balloon-launched apogee', b['apogee'] / 1e3, R['balAp'], 1e-6, 'km')
    check(g, 'ground-launched apogee', gr['apogee'] / 1e3, R['gndAp'], 1e-6, 'km')
    check(g, 'drag share', 100 * s['drag'] / tot, R['decFDrag'], 1e-6, '%')
    check(g, 'nozzle share', 100 * s['nozzle'] / tot, R['decFNozzle'], 1e-6, '%')
    c = rm.close_vehicle(30e3, 80.0, rm.feed_fraction('Al6061T6'))
    check(g, '100 km vehicle gross mass', c['m0'], R['cBalM'], 1e-5, 'kg')
    nb = rm.Nozzle(30e3, 4 * 301.27 * rm.G0, 0.40)
    p2 = rm.fly_tilted(301.27, 151.27, 0.40, nb.mdot, nb.isp_vac, nb.Ae, 30e3, 0.0, 'fixed')
    check(g, 'planar model at zero tilt = vertical model', p2['apogee'] / 1e3, b['apogee'] / 1e3, 1e-4, 'km')
    a11 = rm.fly(301.27, 151.27, 0.40, nb.mdot, nb.isp_vac, nb.Ae, z_kin=30e3, rtol=1e-11)['apogee']
    check(g, 'integrator converged (tolerance 1e-8 vs 1e-11)', b['apogee'], a11, 1e-6, 'm')
    bu = rm.DragBuildUp(); T0, p0, r0 = rm.atmosphere(0.0); a0 = math.sqrt(1.4 * 287.05 * T0)
    worst = max(abs(bu(m, r0, T0, m * a0) / rm.cd_generic(m) - 1.0) for m in (2.0, 3.0, 4.0, 5.0))
    flag(g, 'assumed drag curve within 20% of the computed build-up for Mach 2-5', worst <= 0.20, 'largest difference %.0f%%' % (100 * worst))
    sob = R['sobol']
    flag(g, 'Sobol first-order indices sum to <= 1 at every altitude', all(sum(x['S']) <= 1.02 for x in sob),
         'max sum %.2f' % max(sum(x['S']) for x in sob))
    flag(g, 'ensemble size', R['ensN'] >= 700, '%d vehicles' % R['ensN'])


def paper():
    g = 'Paper consistency'
    for name, cmd in (('unit tests (pytest)', [sys.executable, '-m', 'pytest', '-q', 'test_model.py']),
                      ('check_paper.py', [sys.executable, 'check_paper.py'])):
        p = subprocess.run(cmd, cwd=CODE, capture_output=True, text=True)
        flag(g, name, p.returncode == 0, (p.stdout.strip().splitlines() or [''])[-1][:60])


def sha(path):
    return hashlib.sha256(open(path, 'rb').read()).hexdigest()


def full():
    g = 'Full regeneration'
    tmp = tempfile.mkdtemp()
    dst = os.path.join(tmp, 'repo')
    shutil.copytree(ROOT, dst, ignore=shutil.ignore_patterns('.git', '__pycache__', '.pytest_cache', '*.pdf', '*.gz'))
    for sub in ('results', os.path.join('code', 'thermo_table.json'), os.path.join('paper', 'results.tex')):
        pth = os.path.join(dst, sub)
        shutil.rmtree(pth) if os.path.isdir(pth) else os.remove(pth)
    for script in ('rockoon_model.py', 'similarity.py', 'extra_tests.py', 'cea_check.py', 'rocketpy_check.py', 'revision_tests.py',
                   'realdata_tests.py', 'aerobee_check.py', 'run_analysis.py', 'make_tables.py'):
        p = subprocess.run([sys.executable, script], cwd=os.path.join(dst, 'code'), capture_output=True, text=True)
        if p.returncode != 0:
            flag(g, script, False, p.stderr.strip().splitlines()[-1][:60])
            return
    for f in ['results.tex'] + sorted(x for x in os.listdir(os.path.join(ROOT, 'paper')) if x.startswith('tab_')):
        flag(g, 'paper/%s identical after regeneration' % f, sha(os.path.join(ROOT, 'paper', f)) == sha(os.path.join(dst, 'paper', f)))
    shutil.rmtree(tmp)


def main():
    physics(); external(); headline(); rule(); paper()
    if '--full' in sys.argv:
        print('Regenerating every result in a temporary copy (about an hour on 2 cores)...')
        full()
    last = None
    w = max(len(r[1]) for r in ROWS)
    for grp, name, got, want, err, tol, unit, ok in ROWS:
        if grp != last:
            print('\n' + grp); last = grp
        tail = ('  got %s, expected %s %s, off by %s (limit %s)' % (got, want, unit, err, tol)) if want else ('  ' + got)
        print('  [%s] %s%s' % ('PASS' if ok else 'FAIL', name.ljust(w), tail))
    bad = sum(1 for r in ROWS if not r[-1])
    print('\n%d checks, %d failed' % (len(ROWS), bad))
    sys.exit(1 if bad else 0)


if __name__ == '__main__':
    main()
