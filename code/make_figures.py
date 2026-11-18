"""Draw every figure of the paper from results/results.json and the model."""
import json, math, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import rockoon_model as rm

ROOT = os.path.dirname(rm.HERE)
OUT = os.path.join(ROOT, 'paper', 'figures')
os.makedirs(OUT, exist_ok=True)
R = json.load(open(os.path.join(ROOT, 'results', 'results.json')))
plt.rcParams.update({'font.size': 8, 'axes.labelsize': 8, 'xtick.labelsize': 7.5, 'ytick.labelsize': 7.5,
                     'legend.fontsize': 7.5, 'pdf.fonttype': 42, 'axes.spines.top': False, 'axes.spines.right': False,
                     'axes.grid': True, 'grid.color': '0.9', 'grid.linewidth': 0.6, 'lines.linewidth': 1.4,
                     'font.family': 'DejaVu Sans', 'legend.frameon': False, 'figure.dpi': 300})
BLUE, ORANGE, GREEN, GREY, PURPLE = '#2166ac', '#d6604d', '#1b9e77', '#666666', '#7b3294'
W2 = 6.6


def panel(ax, s):
    ax.text(-0.14, 1.04, s, transform=ax.transAxes, fontweight='bold', fontsize=9, va='bottom')


def save(fig, name):
    fig.savefig(os.path.join(OUT, name + '.pdf'), bbox_inches='tight', pad_inches=0.02)
    plt.close(fig)


def fig_inputs():
    fig, ax = plt.subplots(1, 2, figsize=(W2, 2.5)); fig.subplots_adjust(wspace=0.33)
    M = np.linspace(0, 6, 400); cd = np.array([rm.cd_generic(m) for m in M])
    ax[0].fill_between(M, 0.8 * cd, 1.2 * cd, color=BLUE, alpha=0.18, lw=0, label='scale factor 0.8-1.2')
    ax[0].plot(M, cd, color=BLUE, label='assumed curve')
    ax[0].set_xlabel('Mach number'); ax[0].set_ylabel('Drag coefficient $C_D$'); ax[0].set_ylim(0, 0.7); ax[0].legend(loc='upper right')
    nz = rm.table()['nozzle']
    ax[1].semilogx(nz['eps'], 0.92 * nz['isp_vac'], color=BLUE)
    for key, col, lab in (('gnd', GREEN, 'ground'), ('bal', ORANGE, '30 km')):
        c = R['cases'][key]
        ax[1].plot(c['eps'], c['isp_vac'], 'o', color=col, ms=5)
        ax[1].annotate('%s: $\\varepsilon$ = %.1f' % (lab, c['eps']), (c['eps'], c['isp_vac']), xytext=(8, -10), textcoords='offset points', color=col)
    ax[1].set_xlim(1.5, 300); ax[1].set_ylim(220, 340)
    ax[1].set_xlabel('Nozzle area ratio $\\varepsilon$'); ax[1].set_ylabel('Delivered vacuum $I_{sp}$ (s)')
    panel(ax[0], '(a)'); panel(ax[1], '(b)'); save(fig, 'fig_inputs')


def fig_validation():
    fig, ax = plt.subplots(1, 2, figsize=(W2, 2.6)); fig.subplots_adjust(wspace=0.33)
    d = R['deacon_rounds']
    vm = np.array([x['v_meas'] for x in d]); vp = np.array([x['v_model'] for x in d]); vn = np.array([x['v_nodrag'] for x in d])
    lim = [950, 1450]
    ax[0].plot(lim, lim, color=GREY, lw=0.8)
    ax[0].plot(vm[:6], vp[:6], 'o', color=BLUE, ms=4.5, label='model, rounds 2-7')
    ax[0].plot(vm[6:], vp[6:], 's', color=ORANGE, ms=4.5, label='model, rounds 8-9')
    ax[0].plot(vm, vn, 'x', color=GREY, ms=4, label='model without drag')
    ax[0].set_xlim(lim); ax[0].set_ylim(lim); ax[0].set_aspect('equal')
    ax[0].set_xlabel('Measured maximum velocity (m s$^{-1}$)'); ax[0].set_ylabel('Computed maximum velocity (m s$^{-1}$)'); ax[0].legend(loc='lower right', fontsize=6.5)
    LB = 0.45359237
    mp_d, d_d = 103.0 * LB, 6.75 * 0.0254
    ms = np.linspace(80, 115, 36)
    def ap(m, cd):
        return rm.fly(m, mp_d, d_d, mp_d / 2.9, 19680 * 4.4482216 / (mp_d * rm.G0), 0.0, z_kin=38000 * 0.3048, cd_scale=cd)['apogee'] / 1e3
    lo = [ap(m, 1.2) for m in ms]; hi = [ap(m, 0.8) for m in ms]; mid = [ap(m, 1.0) for m in ms]
    ax[1].axvspan(R['nacaMassLo'], R['nacaMassHi'], color='0.92', lw=0)
    ax[1].fill_between(ms, lo, hi, color=BLUE, alpha=0.18, lw=0)
    ax[1].plot(ms, mid, color=BLUE)
    ax[1].axhline(R['suiTop'], color=ORANGE, lw=1.0)
    ax[1].text(81, R['suiTop'] + 1.5, 'reported summit, about 61 km', color=ORANGE, fontsize=7)
    ax[1].text(0.5 * (R['nacaMassLo'] + R['nacaMassHi']), 108, 'gross masses of the\nNACA ground rounds', ha='center', va='top', fontsize=7, color=GREY)
    ax[1].set_xlim(80, 115); ax[1].set_ylim(40, 105)
    ax[1].set_xlabel('Gross mass at ignition (kg)'); ax[1].set_ylabel('Computed apogee (km)')
    panel(ax[0], '(a)'); panel(ax[1], '(b)'); save(fig, 'fig_validation')


def fig_trajectory():
    fig, ax = plt.subplots(1, 2, figsize=(W2, 2.5)); fig.subplots_adjust(wspace=0.33)
    M0, MP = 301.27, 151.27
    specs = [('Ignited at rest at 30 km', BLUE, 30e3, 0.0), ('Ignited at 29 km after a 1 km fall', ORANGE, 29e3, -R['fallSpeed']), ('Ground launch', GREEN, 0.0, 0.0)]
    for lab, col, z, v0 in specs:
        nz = rm.Nozzle(z, 4 * M0 * rm.G0, 0.40)
        r = rm.fly(M0, MP, 0.40, nz.mdot, nz.isp_vac, nz.Ae, z_kin=z, v0=v0, trace=True)
        t, zz, vv = r['trace']
        ax[0].plot(t, zz / 1e3, color=col, label=lab); ax[1].plot(t, vv, color=col)
        ax[0].plot(r['tb'], (z + r['h_b']) / 1e3, 'o', color=col, ms=3.5); ax[1].plot(r['tb'], r['v_b'], 'o', color=col, ms=3.5)
        ax[0].annotate('%.0f km' % (r['apogee'] / 1e3), (t[-1], zz[-1] / 1e3), xytext=(-4, 4), textcoords='offset points', ha='right', color=col, fontsize=7)
    ax[0].axhline(100, color=GREY, lw=0.8); ax[0].text(2, 103, '100 km', color=GREY, fontsize=7)
    ax[0].set_xlabel('Time from ignition (s)'); ax[0].set_ylabel('Altitude (km)'); ax[0].set_ylim(0, 240)
    ax[1].set_xlabel('Time from ignition (s)'); ax[1].set_ylabel('Vertical velocity (m s$^{-1}$)')
    fig.legend(*ax[0].get_legend_handles_labels(), loc='lower center', ncol=3, bbox_to_anchor=(0.5, -0.2))
    panel(ax[0], '(a)'); panel(ax[1], '(b)'); save(fig, 'fig_trajectory')


def stack(ax, x, rows, width, xl):
    h = np.array([r['head'] for r in rows]); d = np.array([r['drag'] for r in rows]); n = np.array([r['nozzle'] for r in rows]); g = np.array([r['ground'] for r in rows])
    ax.bar(x, g, width, color='0.75', label='ground-launch apogee')
    ax.bar(x, d, width, bottom=g, color=BLUE, label='drag removed')
    ax.bar(x, n, width, bottom=g + d, color=ORANGE, label='nozzle')
    ax.bar(x, h, width, bottom=g + d + n, color=GREEN, label='head start')
    ax.set_xlabel(xl)


def fig_decomposition():
    fig, ax = plt.subplots(1, 3, figsize=(W2, 2.5), gridspec_kw={'width_ratios': [1.15, 1, 1]}); fig.subplots_adjust(wspace=0.42)
    rows = R['dec_alt']
    stack(ax[0], [r['zf'] / 1e3 for r in rows], rows, 3.6, 'Float altitude (km)')
    ax[0].set_ylabel('Apogee (km)'); ax[0].set_ylim(0, 330)
    ax[0].legend(loc='upper left', fontsize=6.5, handlelength=1.0, labelspacing=0.25, borderaxespad=0.1)
    for a, key, xl in ((ax[1], 'dec_alt', 'Float altitude (km)'), (ax[2], None, 'Ballistic parameter $m_0/A$ (kg m$^{-2}$)')):
        if key:
            rr = R[key]; x = [r['zf'] / 1e3 for r in rr]
            for k, c, lab in (('f_drag', BLUE, 'drag'), ('f_nozzle', ORANGE, 'nozzle'), ('f_head', GREEN, 'head start')):
                a.plot(x, [100 * r[k] for r in rr], 'o-', color=c, ms=3, label=lab)
        else:
            for kk, ls, mk in (('dec_mass', '-', 'o'), ('dec_diam', '--', 's')):
                rr = sorted(R[kk], key=lambda r: r['beta'])
                for k, c in (('f_drag', BLUE), ('f_nozzle', ORANGE), ('f_head', GREEN)):
                    a.plot([r['beta'] for r in rr], [100 * r[k] for r in rr], ls=ls, marker=mk, color=c, ms=3, mfc='w' if kk == 'dec_diam' else c)
            a.set_xscale('log'); a.set_xticks([1000, 2000, 4000]); a.set_xticklabels(['1000', '2000', '4000']); a.minorticks_off()
        a.set_ylim(0, 75); a.set_xlabel(xl); a.set_ylabel('Share of apogee gain (%)')
    ax[1].text(22, 57, 'drag', color=BLUE, fontsize=7); ax[1].text(22, 34, 'nozzle', color=ORANGE, fontsize=7); ax[1].text(22, 7, 'head start', color=GREEN, fontsize=7)
    for a, s in zip(ax, ('(a)', '(b)', '(c)')):
        a.text(-0.2, 1.04, s, transform=a.transAxes, fontweight='bold', fontsize=9, va='bottom')
    save(fig, 'fig_decomposition')


def fig_scaling():
    fig, ax = plt.subplots(1, 2, figsize=(W2, 2.6)); fig.subplots_adjust(wspace=0.36)
    m = R['map']; cols = plt.cm.viridis(np.linspace(0.05, 0.85, len(m['m'])))
    for mm, row, c in zip(m['m'], m['apogee'], cols):
        ax[0].plot(m['z'], row, color=c, label='%.0f kg' % mm)
    ax[0].axhline(100, color=GREY, lw=0.8)
    ax[0].set_xlabel('Ignition altitude (km)'); ax[0].set_ylabel('Apogee (km)'); ax[0].set_ylim(0, 240); ax[0].legend(title='Gross mass', title_fontsize=7.5, loc='lower right')
    a = R['alt_cost']; z = [x['z'] for x in a]
    ax[1].semilogy(z, [x['V'] for x in a], color=BLUE)
    ax[1].set_xlabel('Float altitude (km)'); ax[1].set_ylabel('Balloon volume (m$^3$)', color=BLUE); ax[1].tick_params(axis='y', colors=BLUE)
    b = ax[1].twinx(); b.spines['right'].set_visible(True); b.grid(False)
    b.plot(z, [x['apogee'] for x in a], color=ORANGE); b.set_ylabel('Apogee (km)', color=ORANGE); b.tick_params(axis='y', colors=ORANGE); b.set_ylim(150, 240)
    panel(ax[0], '(a)'); panel(ax[1], '(b)'); save(fig, 'fig_scaling')


def fig_closure():
    fig, ax = plt.subplots(1, 2, figsize=(W2, 2.6)); fig.subplots_adjust(wspace=0.33)
    for key, col, ls, lab in (('Al_40', GREEN, '-', 'Al, 40 kg'), ('Al_80', BLUE, '-', 'Al, 80 kg'), ('Al_120', PURPLE, '-', 'Al, 120 kg'),
                              ('SS_40', GREEN, '--', 'steel, 40 kg'), ('SS_80', BLUE, '--', 'steel, 80 kg'), ('SS_120', PURPLE, '--', 'steel, 120 kg')):
        rows = [r for r in R['closure'][key] if r['m0']]
        ax[0].semilogy([r['z'] for r in rows], [r['m0'] for r in rows], color=col, ls=ls, marker='o', ms=2.5, label=lab)
        rb = [r for r in rows if 'V' in r]
        ax[1].semilogy([r['z'] for r in rb], [r['V'] for r in rb], color=col, ls=ls, marker='o', ms=2.5)
    ax[0].set_xlabel('Ignition altitude (km)'); ax[0].set_ylabel('Gross mass for a 100 km apogee (kg)')
    fig.legend(*ax[0].get_legend_handles_labels(), title='Tank material and non-tank dry mass', title_fontsize=7.5, ncol=6, loc='lower center', bbox_to_anchor=(0.5, -0.2), columnspacing=1.0, handlelength=1.6, fontsize=7)
    ax[0].set_yticks([50, 100, 200, 500, 1000]); ax[0].set_yticklabels(['50', '100', '200', '500', '1000']); ax[0].minorticks_off()
    ax[1].set_xlabel('Float altitude (km)'); ax[1].set_ylabel('Balloon volume (m$^3$)'); ax[1].set_xlim(0, 41)
    panel(ax[0], '(a)'); panel(ax[1], '(b)'); save(fig, 'fig_closure')


def fig_montecarlo():
    fig, ax = plt.subplots(1, 2, figsize=(W2, 2.4)); fig.subplots_adjust(wspace=0.33)
    mc = R['mc']
    ax[0].hist(mc['gnd'], bins=40, color=GREEN, alpha=0.8, label='ground launch')
    ax[0].hist(mc['bal'], bins=60, color=BLUE, alpha=0.8, label='ignition at 28-32 km')
    ax[0].axvline(100, color=GREY, lw=0.8); ax[0].set_xlabel('Apogee (km)'); ax[0].set_ylabel('Samples'); ax[0].legend(loc='upper center')
    sh = 100 * np.array(mc['shares'])
    for i, (c, lab) in enumerate(((GREEN, 'head start'), (BLUE, 'drag'), (ORANGE, 'nozzle'))):
        ax[1].hist(sh[:, i], bins=30, color=c, alpha=0.8, label=lab)
    ax[1].set_xlabel('Share of apogee gain (%)'); ax[1].set_ylabel('Samples'); ax[1].legend(loc='upper center', ncol=3, columnspacing=0.8, handlelength=1.0)
    panel(ax[0], '(a)'); panel(ax[1], '(b)'); save(fig, 'fig_montecarlo')


def fig_appendix():
    fig, ax = plt.subplots(1, 3, figsize=(W2, 2.3)); fig.subplots_adjust(wspace=0.45)
    z = np.linspace(0, 44, 300)
    ax[0].plot([rm.atmosphere(x * 1e3)[0] for x in z], z, color=BLUE, label='USSA 1976')
    ax[0].plot([rm.lapse_extrapolation(x * 1e3)[0] for x in z], z, color=ORANGE, label='lapse rate\nextrapolated')
    ax[0].axhline(30, color=GREY, lw=0.8); ax[0].set_xlabel('Temperature (K)'); ax[0].set_ylabel('Altitude (km)'); ax[0].set_xlim(0, 300); ax[0].legend(loc='upper left', fontsize=6.5)
    sw = R_T['of_sweep']; of = [s['of'] for s in sw]
    ax[1].plot(of, [s['shift'] for s in sw], color=BLUE, label='shifting'); ax[1].plot(of, [s['frozen'] for s in sw], color=ORANGE, label='frozen')
    ax[1].axhline(316.7, color=GREY, lw=0.8); ax[1].text(2.4, 318, '316.7 s (proposal)', fontsize=6.5, color=GREY)
    ax[1].set_xlabel('Mixture ratio O/F'); ax[1].set_ylabel('Ideal vacuum $I_{sp}$ (s), $\\varepsilon$ = 36.4'); ax[1].legend(loc='lower right')
    ax[2].plot(of, [s['Tc'] for s in sw], color=BLUE); ax[2].set_xlabel('Mixture ratio O/F'); ax[2].set_ylabel('Chamber temperature (K)')
    for a, s in zip(ax, ('(a)', '(b)', '(c)')):
        a.text(-0.24, 1.04, s, transform=a.transAxes, fontweight='bold', fontsize=9, va='bottom')
    save(fig, 'fig_appendix')


def fig_collapse():
    S = json.load(open(os.path.join(ROOT, 'results', 'similarity.json')))['ensemble']
    g = lambda k: np.array([e[k] for e in S])
    fig, ax = plt.subplots(1, 3, figsize=(W2, 2.3)); fig.subplots_adjust(wspace=0.55)
    eng = g('eng'); names = ['liquid, 1 MPa', 'liquid, 2 MPa', 'liquid, 4 MPa', 'solid, 6 MPa']; cols = [GREEN, BLUE, PURPLE, ORANGE]
    for i in range(4):
        m = eng == i
        ax[0].loglog(g('pi0')[m], (g('ld_g') / g('dv_g'))[m], '.', ms=2.5, color=cols[i], label=names[i])
    x = np.array([0.3, 5.0]); ax[0].loglog(x, R['ldCoef'] * x ** R['ldSlope'], color='k', lw=0.8)
    ax[0].set_xticks([0.3, 1, 3]); ax[0].set_xticklabels(['0.3', '1', '3']); ax[0].set_yticks([0.1, 0.2, 0.5]); ax[0].set_yticklabels(['0.1', '0.2', '0.5']); ax[0].minorticks_off()
    ax[0].set_xlabel('$\\Pi_D(0)$'); ax[0].set_ylabel('Ground drag loss / thrust $\\Delta v$'); ax[0].legend(loc='lower right', fontsize=6, handletextpad=0.1, markerscale=2)
    sc = ax[1].scatter(g('pi0'), 100 * g('f_drag'), c=g('z0') / 1e3, s=3, cmap='viridis', lw=0)
    ax[1].set_xscale('log'); ax[1].set_xticks([0.3, 1, 3]); ax[1].set_xticklabels(['0.3', '1', '3']); ax[1].minorticks_off()
    ax[1].set_xlabel('$\\Pi_D(0)$'); ax[1].set_ylabel('Drag share of gain (%)')
    cax = ax[1].inset_axes([0.05, 0.9, 0.45, 0.04]); cb = fig.colorbar(sc, cax=cax, orientation='horizontal'); cb.set_ticks([10, 25, 40]); cb.ax.tick_params(labelsize=6); cax.set_title('$z_0$ (km)', fontsize=6.5, pad=1)
    pred = 100 * np.array(R['ens_fit_sd'])
    ax[2].plot([15, 80], [15, 80], color=GREY, lw=0.8)
    for i in range(4):
        m = eng == i
        ax[2].plot(pred[m], 100 * g('f_drag')[m], '.', ms=2.5, color=cols[i])
    ax[2].set_xlabel('Drag share from six-group fit (%)'); ax[2].set_ylabel('Drag share of gain (%)'); ax[2].set_xlim(15, 80); ax[2].set_ylim(15, 80)
    for a, t in zip(ax, ('(a)', '(b)', '(c)')):
        a.text(-0.26, 1.04, t, transform=a.transAxes, fontweight='bold', fontsize=9, va='bottom')
    save(fig, 'fig_collapse')


def fig_optimum():
    fig, ax = plt.subplots(1, 2, figsize=(W2, 2.6)); fig.subplots_adjust(wspace=0.33)
    opt = [o for o in R['optimum'] if o['target'] == 100.0]
    al = [o for o in opt if o['mat'] == 'Al6061T6']; ss = [o for o in opt if o['mat'] == 'SS321']
    cols = plt.cm.viridis(np.linspace(0.05, 0.9, len(al)))
    for o, c in zip(al, cols):
        ax[0].plot(o['curve'][0], o['curve'][1], color=c, label='%.0f kg' % o['mfix'])
    ax[0].axhline(1.0, color=GREY, lw=0.8)
    ax[0].set_xlabel('Float altitude (km)'); ax[0].set_ylabel('Launched mass / ground-launched mass'); ax[0].set_ylim(0.3, 1.4)
    ax[0].legend(title='$m_{fix}$ (Al tanks)', title_fontsize=7, fontsize=6.5, ncol=2, loc='upper center', columnspacing=0.8)
    for fam, col, lab, mk in ((al, BLUE, 'aluminium tanks', 'o'), (ss, ORANGE, 'steel tanks', 's')):
        x = [o['ground'] for o in fam]; z = [o['zstar'] for o in fam]
        ax[1].errorbar(x, z, yerr=[[zz - o['lo'] for zz, o in zip(z, fam)], [o['hi'] - zz for zz, o in zip(z, fam)]], fmt=mk, color=col, ms=3.5, lw=0.7, capsize=1.5, label=lab)
    oth = [o for o in R['optimum'] if o['target'] != 100.0]
    ax[1].plot([o['ground'] for o in oth], [o['zstar'] for o in oth], '^', color=GREEN, ms=3.5, label='60 and 150 km targets')
    ax[1].set_xscale('log'); ax[1].set_xlabel('Gross mass of the ground-launched vehicle (kg)'); ax[1].set_ylabel('Mass-optimal float altitude (km)'); ax[1].set_ylim(5, 35)
    ax[1].legend(loc='upper left', fontsize=6.5)
    panel(ax[0], '(a)'); panel(ax[1], '(b)'); save(fig, 'fig_optimum')


def fig_sobol():
    fig, ax = plt.subplots(1, 2, figsize=(W2, 2.4)); fig.subplots_adjust(wspace=0.33)
    sb = R['sobol']; z = [x['z'] for x in sb]
    for i, c, lab in ((0, BLUE, 'drag scale factor'), (1, ORANGE, 'efficiency $\\eta$'), (3, GREEN, 'dry mass'), (2, PURPLE, 'thrust-to-weight')):
        ax[0].plot(z, [max(x['ST'][i], 0) for x in sb], 'o-', color=c, ms=3, label=lab)
    ax[0].set_xlabel('Ignition altitude (km)'); ax[0].set_ylabel('Total-effect Sobol index'); ax[0].set_ylim(0, 0.85); ax[0].legend(loc='upper right', fontsize=6.5, ncol=2, columnspacing=0.8)
    ax[1].semilogx([x['piz'] for x in sb], [max(x['ST'][0], 0) for x in sb], 'o-', color=BLUE, ms=3)
    ax[1].set_xlabel('$\\Pi_D(z_0)$'); ax[1].set_ylabel('Sobol index of the drag scale factor'); ax[1].set_ylim(0, 0.55)
    panel(ax[0], '(a)'); panel(ax[1], '(b)'); save(fig, 'fig_sobol')


if __name__ == '__main__':
    R_T = json.load(open(os.path.join(rm.HERE, 'thermo_table.json')))
    fig_inputs(); fig_validation(); fig_trajectory(); fig_decomposition(); fig_scaling(); fig_closure(); fig_montecarlo(); fig_appendix(); fig_collapse(); fig_optimum(); fig_sobol()
    print('figures written to', OUT)
