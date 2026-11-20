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
    bu = rm.DragBuildUp()
    for z_, ls, lab in ((0.0, '-', 'build-up, sea level'), (30e3, '--', 'build-up, 30 km')):
        T_, p_, r_ = rm.atmosphere(z_); a_ = math.sqrt(1.4 * 287.05 * T_)
        ax[0].plot(M[5:], [bu(m, r_, T_, m * a_) for m in M[5:]], color=ORANGE, ls=ls, lw=1.1, label=lab)
    ax[0].plot(M, [rm.cd_superloki(m) for m in M], color=GREY, lw=0.9, ls=':', label='Super Loki table')
    ax[0].set_xlabel('Mach number'); ax[0].set_ylabel('Drag coefficient $C_D$'); ax[0].set_ylim(0, 1.0); ax[0].legend(loc='upper right', fontsize=6.3)
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
    def ap(m, zkm, cd=1.0):
        return rm.fly(m, mp_d, d_d, mp_d / 2.9, 19680 * 4.4482216 / (mp_d * rm.G0), 0.0, z_kin=zkm * 1e3, cd_scale=cd)['apogee'] / 1e3
    zz = np.linspace(8, 27, 20)
    for grp, col, mk, key in (('SUI', BLUE, 'o', 'Sui'), ('NRL', ORANGE, 's', 'Nrl')):
        m_ = R['rockoon_mass'][key]
        ax[1].fill_between(zz, [ap(m_, z, 1.2) for z in zz], [ap(m_, z, 0.8) for z in zz], color=col, alpha=0.15, lw=0)
        ax[1].plot(zz, [ap(m_, z) for z in zz], color=col, lw=1.1, label='model, %.0f kg' % m_)
        pts = [r for r in R['rockoons'] if r['grp'] == grp]
        ax[1].plot([r['fire'] for r in pts], [r['obs'] for r in pts], mk, color=col, ms=4.5, mfc='w', mew=1.2, label='%s flights' % ('Iowa' if grp == 'SUI' else 'NRL'))
    ax[1].set_xlabel('Firing altitude (km)'); ax[1].set_ylabel('Apogee (km)'); ax[1].set_ylim(40, 110); ax[1].legend(loc='upper left', fontsize=6.3, ncol=2, columnspacing=0.8)
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
    fig, ax = plt.subplots(1, 2, figsize=(W2, 2.4)); fig.subplots_adjust(wspace=0.35)
    sw = R_T['of_sweep']; of = [s['of'] for s in sw]
    ax[0].plot(of, [s['shift'] for s in sw], color=BLUE, label='shifting (Cantera)'); ax[0].plot(of, [s['frozen'] for s in sw], color=ORANGE, label='frozen (Cantera)')
    C = json.load(open(os.path.join(ROOT, 'results', 'cea_check.json')))['sweep']
    ax[0].plot([c['of'] for c in C], [c['cea_shift'] for c in C], 'o', mfc='none', color=BLUE, ms=4.5, label='shifting (NASA CEA)')
    ax[0].plot([c['of'] for c in C], [c['cea_frozen'] for c in C], 's', mfc='none', color=ORANGE, ms=4.5, label='frozen (NASA CEA)')
    ax[0].set_xlabel('Mixture ratio O/F'); ax[0].set_ylabel('Ideal vacuum $I_{sp}$ (s), $\\varepsilon$ = 36.4'); ax[0].legend(loc='lower right', fontsize=6.3)
    ax[1].plot(of, [s['Tc'] for s in sw], color=BLUE, label='Cantera')
    ax[1].plot([c['of'] for c in C], [c['cea_Tc'] for c in C], 'o', mfc='none', color=BLUE, ms=4.5, label='NASA CEA')
    ax[1].set_xlabel('Mixture ratio O/F'); ax[1].set_ylabel('Chamber temperature (K)'); ax[1].legend(loc='lower right', fontsize=6.5)
    panel(ax[0], '(a)'); panel(ax[1], '(b)'); save(fig, 'fig_appendix')


def fig_analytic():
    """The envelope rule: at the optimum the envelope weighs 3H/2 times the suspended-mass saving per kilometre."""
    S = json.load(open(os.path.join(ROOT, 'results', 'similarity.json')))['system']
    fig, ax = plt.subplots(1, 2, figsize=(W2, 2.6)); fig.subplots_adjust(wspace=0.33)
    cols = {20.0: GREEN, 80.0: BLUE, 320.0: PURPLE}
    for mfix, c in cols.items():
        v = sorted([r for r in S if r['mat'] == 'Al6061T6' and r['mfix'] == mfix and r['target'] == 100e3 and r['sigma'] == 0.030
                    and r['gf'] == 0.10 and r['z'] > 0], key=lambda r: r['z'])
        z = np.array([r['z'] for r in v]); ms = 1.1 * np.array([r['m0'] for r in v]); me = np.array([r['env'] for r in v])
        dms = -np.gradient(ms, z)
        H = np.array([-0.1 / math.log(rm.atmosphere(zz * 1e3 + 50)[2] / rm.atmosphere(zz * 1e3 - 50)[2]) for zz in z])
        sel = (z >= 7.5) & (z <= 35)
        ax[0].semilogy(z[sel], me[sel], color=c, ls='--', lw=1.1)
        ax[0].semilogy(z[sel], 1.5 * H[sel] * dms[sel], color=c, lw=1.4, label='$m_{fix}$ = %.0f kg' % mfix)
        o = [x for x in R['optimum'] if x['mat'] == 'Al6061T6' and x['mfix'] == mfix and x['target'] == 100.0][0]
        ax[0].axvline(o['zstar'], color=c, lw=0.6, ls=':')
    ax[0].set_xlabel('Float altitude (km)'); ax[0].set_ylabel('Mass (kg)')
    ax[0].plot([], [], color=GREY, lw=1.4, label='$\\frac{3}{2}H\\,|dm_s/dz|$'); ax[0].plot([], [], color=GREY, ls='--', lw=1.1, label='envelope $m_e$')
    ax[0].legend(loc='lower center', bbox_to_anchor=(0.5, 1.0), fontsize=6, ncol=3, columnspacing=0.8, handlelength=1.6)
    an = R['revision']['analytic']
    for mat, col, mk, lab in (('Al6061T6', BLUE, 'o', 'aluminium tanks'), ('SS321', ORANGE, 's', 'steel tanks')):
        a = [r for r in an if r['mat'] == mat]
        ax[1].plot([r['z_model'] for r in a], [r['z_rule'] for r in a], mk, color=col, ms=4, label=lab)
    ax[1].plot([15, 25], [15, 25], color=GREY, lw=0.8)
    ax[1].set_xlim(15.5, 24.5); ax[1].set_ylim(15.5, 24.5)
    ax[1].set_xlabel('Optimum of the full model, $z^*$ (km)'); ax[1].set_ylabel('Optimum from the envelope rule (km)'); ax[1].legend(loc='upper left', fontsize=6.5)
    panel(ax[0], '(a)'); panel(ax[1], '(b)'); save(fig, 'fig_analytic')


def fig_saving():
    """Monte Carlo over the rockoon-only penalties: the location of the optimum against the size of the saving."""
    mc = R['revision']['mc_saving']
    fig, ax = plt.subplots(1, 2, figsize=(W2, 2.5)); fig.subplots_adjust(wspace=0.3)
    labels, pos = [], []
    for j, (mat, col) in enumerate((('Al6061T6', BLUE), ('SS321', ORANGE))):
        for i, mfix in enumerate((10, 20, 40, 80, 160, 320, 640)):
            k = '%s|%g' % (mat, mfix); x = i + 8 * j
            for a, key, sc in ((ax[0], 'save', 100.0), (ax[1], 'z', 1.0)):
                b = a.boxplot([sc * np.array(mc[k][key])], positions=[x], widths=0.6, whis=(5, 95), showfliers=False, patch_artist=True)
                for p in b['boxes']:
                    p.set_facecolor(col); p.set_alpha(0.55); p.set_edgecolor(col)
                for part in ('whiskers', 'caps', 'medians'):
                    for l in b[part]:
                        l.set_color(col if part != 'medians' else 'k'); l.set_linewidth(0.9)
            labels.append('%d' % mfix); pos.append(x)
    for a in ax:
        a.set_xticks(pos); a.set_xticklabels(labels, rotation=90, fontsize=6.3); a.set_xlim(-1, 15)
        a.set_xlabel('$m_{fix}$ (kg)')
        a.text(3, 0.97, 'aluminium tanks', transform=a.get_xaxis_transform(), ha='center', va='top', fontsize=7, color=BLUE)
        a.text(11, 0.97, 'steel tanks', transform=a.get_xaxis_transform(), ha='center', va='top', fontsize=7, color=ORANGE)
    ax[0].axhline(0, color=GREY, lw=0.8); ax[0].set_ylabel('Saving in launched mass (%)')
    ax[0].set_ylim(-30, 85)
    ax[1].set_ylabel('Mass-optimal float altitude (km)'); ax[1].set_ylim(10, 33)
    panel(ax[0], '(a)'); panel(ax[1], '(b)'); save(fig, 'fig_saving')


def fig_lapse():
    """Supplement: the proposal's extrapolated lapse rate against the standard atmosphere."""
    fig, ax = plt.subplots(figsize=(3.3, 2.6))
    z = np.linspace(0, 44, 300)
    ax.plot([rm.atmosphere(x * 1e3)[0] for x in z], z, color=BLUE, label='USSA 1976')
    ax.plot([rm.lapse_extrapolation(x * 1e3)[0] for x in z], z, color=ORANGE, label='lapse rate\nextrapolated')
    ax.axhline(30, color=GREY, lw=0.8); ax.set_xlabel('Temperature (K)'); ax.set_ylabel('Altitude (km)'); ax.set_xlim(0, 300); ax.legend(loc='upper left', fontsize=6.5)
    save(fig, 'fig_lapse')


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
    fig, ax = plt.subplots(1, 3, figsize=(W2, 2.3)); fig.subplots_adjust(wspace=0.45)
    sb = R['sobol']; z = [x['z'] for x in sb]
    for i, c, lab in ((0, BLUE, 'drag scale factor'), (1, ORANGE, 'efficiency $\\eta$'), (3, GREEN, 'dry mass'), (2, PURPLE, 'thrust-to-weight')):
        ax[0].plot(z, [max(x['ST'][i], 0) for x in sb], 'o-', color=c, ms=3, label=lab)
    ax[0].set_xlabel('Ignition altitude (km)'); ax[0].set_ylabel('Total-effect Sobol index'); ax[0].set_ylim(0, 0.85); ax[0].legend(loc='upper right', fontsize=6, ncol=1)
    X = json.load(open(os.path.join(ROOT, 'results', 'extra.json')))['sobol_vehicles']
    names = []
    for x in X:
        if x['vehicle'] not in names:
            names.append(x['vehicle'])
    cols = [GREEN, BLUE, PURPLE, ORANGE]
    for n, c in zip(names, cols):
        v = [x for x in X if x['vehicle'] == n]
        ax[1].plot([x['z'] for x in v], [max(x['ST_cd'], 0) for x in v], 'o-', color=c, ms=3, label=n)
        ax[2].semilogx([x['piz'] for x in v], [max(x['ST_cd'], 0) for x in v], 'o-', color=c, ms=3, label=n)
    ax[1].set_xlabel('Ignition altitude (km)'); ax[1].set_ylabel('Index of the drag scale factor'); ax[1].legend(loc='upper right', fontsize=6)
    ax[2].set_xlabel('$\\Pi_D(z_0)$'); ax[2].set_ylabel('Index of the drag scale factor')
    for a in ax[1:]:
        a.set_ylim(0, 0.75)
    panel(ax[0], '(a)'); panel(ax[1], '(b)'); panel(ax[2], '(c)'); save(fig, 'fig_sobol')


def fig_general():
    """The envelope rule for several envelope laws: the ratio of the two sides of the optimum condition against altitude,
    the optimum against vehicle size, and the one-step estimate against the exact optimum."""
    Gn = R['general']['laws']
    sty = {'area': (BLUE, '-', 'o', 'film'), 'load': (GREY, '--', 'v', 'film + load'),
           'fitted': (ORANGE, '-', 's', 'fitted'), 'cluster': (GREEN, '-.', '^', 'cluster')}
    fig, ax = plt.subplots(1, 3, figsize=(W2, 2.45)); fig.subplots_adjust(wspace=0.42)
    for tag, (c, ls, mk, lab) in sty.items():
        ln = Gn[tag]['line']; z = np.array(ln['z']); r = np.array(ln['ratio'], dtype=float)
        ax[0].semilogy(z, r, color=c, ls=ls, lw=1.3, label=lab)
        zx = [x for x in Gn[tag]['rows'] if x['mat'] == 'Al6061T6' and x['mfix'] == 80.0][0]['z_exact']
        ax[0].plot([zx], [1.0], mk, color=c, ms=4, mec='white', mew=0.6, zorder=5)
    ax[0].axhline(1.0, color='0.4', lw=0.6)
    ax[0].set_xlim(8, 36); ax[0].set_ylim(0.03, 30)
    ax[0].set_xlabel('Float altitude (km)'); ax[0].set_ylabel('Ratio of the two sides, $\\mathcal{R}$')
    ax[0].legend(loc='lower left', fontsize=6, handlelength=2.4, borderaxespad=0.2)
    for tag in ('area', 'fitted', 'cluster'):
        c, ls, mk, lab = sty[tag]
        for mat, fill in (('Al6061T6', c), ('SS321', 'white')):
            v = [x for x in Gn[tag]['rows'] if x['mat'] == mat]
            ax[1].semilogx([x['ground'] for x in v], [x['z_exact'] for x in v], color=c, ls=ls, lw=0.9)
            ax[1].semilogx([x['ground'] for x in v], [x['z_grid'] for x in v], mk, color=c, mfc=fill, ms=3.6, mew=0.9,
                           label=lab if mat == 'Al6061T6' else None)
    ax[1].set_xlabel('Ground-launched mass (kg)'); ax[1].set_ylabel('Optimum float altitude $z^*$ (km)')
    ax[1].set_ylim(15.5, 27.0); ax[1].legend(loc='lower center', fontsize=6, ncol=3, columnspacing=0.8, handletextpad=0.1, bbox_to_anchor=(0.58, 0.99))
    lo, hi = 15.5, 27.0
    for tag, (c, ls, mk, lab) in sty.items():
        v = Gn[tag]['rows']
        ax[2].plot([x['z_exact'] for x in v], [x['z_one'] for x in v], mk, color=c, ms=3.4, mfc='none' if tag == 'load' else c, mew=0.8, label=lab)
    ax[2].plot([lo, hi], [lo, hi], color='0.4', lw=0.7)
    ax[2].set_xlim(lo, hi); ax[2].set_ylim(lo, hi)
    ax[2].set_xlabel('Exact optimum (km)'); ax[2].set_ylabel('One step from 20 km (km)')
    panel(ax[0], '(a)'); panel(ax[1], '(b)'); panel(ax[2], '(c)'); save(fig, 'fig_general')


if __name__ == '__main__':
    R_T = json.load(open(os.path.join(rm.HERE, 'thermo_table.json')))
    fig_inputs(); fig_validation(); fig_trajectory(); fig_decomposition(); fig_scaling(); fig_closure(); fig_montecarlo(); fig_appendix(); fig_lapse(); fig_collapse(); fig_optimum(); fig_sobol(); fig_analytic(); fig_saving(); fig_general()
    print('figures written to', OUT)
