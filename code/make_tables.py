"""Write the LaTeX table bodies from results/results.json."""
import json, os
import rockoon_model as rm
ROOT = os.path.dirname(rm.HERE)
R = json.load(open(os.path.join(ROOT, 'results', 'results.json')))
out = os.path.join(ROOT, 'paper')


def w(name, lines):
    with open(os.path.join(out, name), 'w') as f:
        f.write('\n'.join(lines) + '\n')


w('tab_deacon.tex', ['%d & %.1f & %.1f & %.0f & %.0f & %+.1f \\\\' % (x['round'], x['m0'], x['imp'], x['v_meas'], x['v_model'], x['err'])
                     for x in R['deacon_rounds']])
lab = {'000': 'ground & ground & ground', '100': '30 km & ground & ground', '010': 'ground & 30 km & ground', '001': 'ground & ground & 30 km',
       '110': '30 km & 30 km & ground', '101': '30 km & ground & 30 km', '011': 'ground & 30 km & 30 km', '111': '30 km & 30 km & 30 km'}
w('tab_runs.tex', ['%s & %.0f \\\\' % (lab[k], R['dec_base']['runs'][k]) for k in ('000', '100', '010', '001', '110', '101', '011', '111')])
rows = []
for r in R['dec_alt']:
    rows.append('%.0f & %.0f & %.0f & %.1f & %.0f & %.0f & %.0f \\\\' % (r['zf'] / 1e3, r['ground'], r['balloon'], r['balloon'] / r['ground'], 100 * r['f_drag'], 100 * r['f_nozzle'], 100 * r['f_head']))
w('tab_dec_alt.tex', rows)
rows = []
for r in R['dec_mass']:
    rows.append('%.0f & %.2f & %.0f & %.0f & %.0f & %.1f & %.0f & %.0f & %.0f \\\\' % (r['m0'], rm.diameter(r['m0']), r['beta'], r['ground'], r['balloon'], r['balloon'] / r['ground'], 100 * r['f_drag'], 100 * r['f_nozzle'], 100 * r['f_head']))
w('tab_dec_mass.tex', rows)
rows = []
al = {r['z']: r for r in R['closure']['Al_80']}; ss = {r['z']: r for r in R['closure']['SS_80']}
for z in (0, 10, 15, 20, 25, 30, 35, 40):
    a, s = al[z], ss[z]
    bal = ('{:,.0f}'.format(round(a['V'], -2)).replace(',', '\\,') + ' & %.0f & %.1f' % (a['env'], a['gas'])) if 'V' in a else '-- & -- & --'
    rows.append('%d & %.0f & %.0f & %.0f & %s \\\\' % (z, a['mp'], a['m0'], s['m0'], bal))
w('tab_closure.tex', rows)
rows = []
for x in R['alt_cost']:
    if x['z'] in (15.0, 20.0, 25.0, 30.0, 35.0, 40.0):
        rows.append('%.0f & %s & %.0f & %.1f & %.0f \\\\' % (x['z'], '{:,.0f}'.format(round(x['V'], -2)).replace(',', '\\,'), x['env'], x['gas'], x['apogee']))
w('tab_altcost.tex', rows)
print('tables written')
