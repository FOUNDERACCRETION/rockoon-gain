"""Cross-check of the Cantera thermochemistry against NASA CEA run directly.

Uses NASA's CEA (Python package `cea`, NASA Glenn, Apache-2.0) with its own thermodynamic library and its own
reactant data for C2H8N2(L),UDMH and N2O4(L). Writes results/cea_check.json."""
import json, os
import numpy as np
import cea
import rockoon_model as rm

ROOT = os.path.dirname(rm.HERE)
G0 = rm.G0


def run(of, pc_bar, supar=None, pi_p=None, frozen=False):
    reac = cea.Mixture(['C2H8N2(L),UDMH', 'N2O4(L)'])
    prod = cea.Mixture(['C2H8N2(L),UDMH', 'N2O4(L)'], products_from_reactants=True)
    solver = cea.RocketSolver(prod, reactants=reac)
    sol = cea.RocketSolution(solver)
    w = reac.of_ratio_to_weights(np.array([0.0, 1.0]), np.array([1.0, 0.0]), of)
    hc = reac.calc_property(cea.ENTHALPY, w, 298.15) / cea.R
    kw = dict(iac=True, hc=hc)
    if supar is not None:
        kw['supar'] = supar
    if pi_p is not None:
        kw['pi_p'] = pi_p
    if frozen:
        kw['n_frz'] = 1
    solver.solve(sol, w, pc_bar, **kw)
    return dict(T=np.array(sol.T).tolist(), P=np.array(sol.P).tolist(), eps=np.array(sol.ae_at).tolist(),
                cstar=float(np.array(sol.c_star)[0]), isp=(np.array(sol.Isp) / G0).tolist(), isp_vac=(np.array(sol.Isp_vacuum) / G0).tolist(),
                converged=bool(sol.converged))


def main():
    t = rm.table()
    nz = t['nozzle']
    out = {}
    # design point: 2.0 MPa, O/F 2.4, the two selected area ratios
    eps_g, eps_b = 3.6809, 36.4
    for tag, fr in (('shift', False), ('frozen', True)):
        r = run(2.4, 20.0, supar=[eps_g, eps_b], frozen=fr)
        out['design_' + tag] = r
    # Cantera values at the same points (shifting), from the stored nozzle table
    out['cantera'] = dict(Tc=nz['Tc'], cstar=nz['cstar'],
                          isp_vac_g=float(np.interp(eps_g, nz['eps'], nz['isp_vac'])),
                          isp_vac_b=float(np.interp(eps_b, nz['eps'], nz['isp_vac'])))
    sw = {s['of']: s for s in t['of_sweep']}
    out['cantera']['frozen_b'] = sw[2.4]['frozen']
    # mixture-ratio sweep at 2.0 MPa, area ratio 36.4
    ofs = [1.6, 2.0, 2.4, 2.75, 3.0, 3.6]
    out['sweep'] = []
    for of in ofs:
        a = run(of, 20.0, supar=[36.4]); b = run(of, 20.0, supar=[36.4], frozen=True)
        c = sw.get(of)
        out['sweep'].append(dict(of=of, cea_shift=a['isp_vac'][-1], cea_frozen=b['isp_vac'][-1], cea_Tc=a['T'][0],
                                 cantera_shift=c['shift'] if c else None, cantera_frozen=c['frozen'] if c else None,
                                 cantera_Tc=c['Tc'] if c else None))
    # benchmark: 68 atm to 1 atm
    pc = 68 * 1.01325
    r = run(2.6, pc, pi_p=[68.0])
    out['bench_shift'] = r
    out['bench_frozen'] = run(2.1, pc, pi_p=[68.0], frozen=True)
    with open(os.path.join(ROOT, 'results', 'cea_check.json'), 'w') as f:
        json.dump(out, f, indent=1)
    print(json.dumps({k: v for k, v in out.items() if k != 'sweep'}, indent=1)[:3000])
    for s in out['sweep']:
        print(s)


if __name__ == '__main__':
    main()
