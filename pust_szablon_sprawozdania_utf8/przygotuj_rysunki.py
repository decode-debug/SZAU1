"""Wektorowe rysunki i tabele z wyników run_part1.m (bez ponownej symulacji)."""
from pathlib import Path
import numpy as np
from scipy.io import loadmat
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / 'wyniki_czesc1' / 'dane_sprawozdanie.mat'
OUT = HERE / 'rysunki' / 'czesc1'
TABLES = HERE / 'dane'
OUT.mkdir(parents=True, exist_ok=True)
TABLES.mkdir(exist_ok=True)
d = loadmat(DATA, simplify_cells=True)
plt.rcParams.update({'font.size': 8, 'axes.labelsize': 9, 'axes.titlesize': 9,
                     'legend.fontsize': 7, 'lines.linewidth': 1.0,
                     'pdf.fonttype': 42, 'savefig.bbox': 'tight'})
colors = ['#1171BE', '#DD5400', '#EDB120']
fmt = FuncFormatter(lambda x, pos: f'{x:g}'.replace('.', ','))

def axes_style(ax):
    ax.grid(True, alpha=.25)
    ax.xaxis.set_major_formatter(fmt)
    ax.yaxis.set_major_formatter(fmt)
    ax.set_xlabel(r'$t$ [s]')
    ax.set_xlim(0, 8000)

def save(fig, name):
    fig.savefig(OUT / (name + '.pdf'))
    png_names = {'skoki_F1in': 'skoki_F1in', 'skoki_FD': 'skoki_FD',
                 'DMC_porownanie_y_u': 'regulacja_DMC_porownanie',
                 'DMC_przypadek_4cm': 'regulacja_DMC'}
    if name.startswith('DMC_scenariusz_'):
        png_names[name] = name
    if name in png_names:
        fig.savefig(DATA.parent / (png_names[name] + '.png'), dpi=160)
    plt.close(fig)

for kind, key, factor in [('F1in', 'stepTracesU', 1), ('FD', 'stepTracesD', .2)]:
    fig, axs = plt.subplots(1, 2, figsize=(6.4, 3.1), layout='constrained')
    for sign, ax in enumerate(axs):
        handles = []
        for c, j in zip(colors, range(sign*3, sign*3+3)):
            ax.step(d['t'], d[key][:, 2*j], where='post', color=c)
            ax.step(d['t'], d[key][:, 2*j+1], where='post', linestyle='--', color=c)
            handles.append(Line2D([], [], color=c, label=f'{d["steps"][j]*factor:+g} cm³/s'))
        handles.extend([Line2D([], [], color='.25', label='model nieliniowy'),
                        Line2D([], [], color='.25', linestyle='--', label='model liniowy')])
        ax.legend(handles=handles, loc='best')
        axes_style(ax)
        ax.set_xlim(0, 6000)
        ax.set_ylabel(r'$y$ [cm]')
        ax.set_title('Skoki ujemne' if sign == 0 else 'Skoki dodatnie')
    save(fig, 'skoki_' + kind)

traces = d['dmcTraces']
for indices, name in [([0], 'DMC_scenariusz_1'),
                      ([1], 'DMC_scenariusz_2'),
                      ([2], 'DMC_scenariusz_3'),
                      ([3], 'DMC_scenariusz_4'),
                      ([0, 1, 2, 3], 'DMC_porownanie_y_u')]:
    rows, cols = (2, 1) if len(indices) == 1 else (4, 2)
    fig, axs = plt.subplots(rows, cols, figsize=(6.4, 4.2 if rows == 2 else 8.5),
                            layout='constrained', squeeze=False)
    for offset, j in enumerate(indices):
        tr, sc = traces[j], d['scenarios'][j]
        row, col = 2*(offset//cols), offset % cols
        ay, au = axs[row, col], axs[row+1, col]
        au.sharex(ay)
        ay.step(d['tc'], tr['yNL'], where='post', color=colors[0], label='obiekt nieliniowy')
        ay.step(d['tc'], tr['yLin'], where='post', linestyle='--', color=colors[1], label='obiekt liniowy')
        ay.step(d['tc'], tr['ref'], where='post', color='.2', linestyle=':',
                label=r'$y_{\mathrm{zad}}$')
        # Ostatnie sterowanie obowiązuje do końca ostatniego przedziału.
        au.step(d['tc'], np.r_[tr['uNL'], tr['uNL'][-1]], where='post',
                color=colors[0], label='obiekt nieliniowy')
        au.step(d['tc'], np.r_[tr['uLin'], tr['uLin'][-1]], '--', where='post',
                color=colors[1], label='obiekt liniowy')
        for ax in (ay, au):
            axes_style(ax)
            ax.axvline(400, color='.6', linewidth=.6, linestyle=':')
            ax.axvline(3800, color='.6', linewidth=.6, linestyle=':')
            ax.set_xticks([0, 2000, 4000, 6000, 8000])
        ay.set_xlabel('')
        ay.tick_params(labelbottom=False)
        ay.set_ylabel(r'$y$ [cm]')
        au.set_ylabel(r'$u$ [cm³/s]')
        ay.set_title(rf'({chr(97+j)}) $\Delta y_{{\mathrm{{zad}}}}={sc[0]:+g}$ cm, $\Delta F_D={sc[1]:+g}$ cm³/s', fontsize=8)
        if offset == 0:
            ay.legend(loc='best', ncol=3 if cols == 1 else 1)
            au.legend(loc='best', ncol=2 if cols == 1 else 1)
    save(fig, name)

tr = traces[1]
fig, axs = plt.subplots(3, 1, figsize=(6.4, 4.9), sharex=True, layout='constrained')
axs[0].step(d['tc'], tr['yNL'], where='post', color=colors[0], label='obiekt nieliniowy')
axs[0].step(d['tc'], tr['yLin'], where='post', linestyle='--', color=colors[1], label='obiekt liniowy')
axs[0].step(d['tc'], tr['ref'], where='post', color='.2', linestyle=':', label=r'$y_{\mathrm{zad}}$')
axs[0].set_ylabel(r'$y$ [cm]')
axs[0].legend(ncol=3, loc='upper right')
axs[1].step(d['tc'], np.r_[tr['uNL'], tr['uNL'][-1]], where='post', color=colors[0])
axs[1].step(d['tc'], np.r_[tr['uLin'], tr['uLin'][-1]], '--', where='post', color=colors[1])
axs[1].set_ylabel(r'$u$ [cm³/s]')
axs[2].step(d['tc'], tr['dist'], where='post', color=colors[2])
axs[2].set_ylabel(r'$F_D$ [cm³/s]')
for ax in axs:
    axes_style(ax)
    ax.axvline(3800, color='.6', linewidth=.6, linestyle=':')
    if ax is not axs[-1]:
        ax.set_xlabel('')
save(fig, 'DMC_przypadek_4cm')

def row(values, decimals):
    return ' & '.join(r'\num{' + f'{float(v):.{n}f}' + '}' for v, n in zip(values, decimals)) + r' \\'

for key, name in [('stepResults', 'skoki_F1in'), ('distResults', 'skoki_FD')]:
    (TABLES / (name + '.tex')).write_text('\n'.join(row(r, [0, 4, 4, 4, 4]) for r in d[key]) + '\n', encoding='utf-8')

(TABLES / 'DMC_wskazniki.tex').write_text('\n'.join(row(r[:4], [0, 0, 3, 3]) for r in d['dmcResults']) + '\n', encoding='utf-8')
control_rows = []
extra = []
for sc, tr in zip(d['scenarios'], traces):
    residual = tr['yNL'][-1] - tr['ref'][-1]
    control_rows.append(row([sc[0], sc[1], np.min(tr['uNL']), np.max(tr['uNL']), residual], [0, 0, 2, 2, 4]))
    mask_r = (d['tc'] >= 400) & (d['tc'] < 3800)
    mask_d = d['tc'] >= 3800
    extra.append([sc[0], sc[1], np.mean((tr['yLin'][mask_r]-tr['ref'][mask_r])**2),
                  np.mean((tr['yLin'][mask_d]-tr['ref'][mask_d])**2)])
(TABLES / 'DMC_sterowanie.tex').write_text('\n'.join(control_rows) + '\n', encoding='utf-8')
print('Utworzono rysunki PDF z parami y(k), u(k), obrazy PNG oraz 4 tabele LaTeX.')
print('MSE modelu liniowego [skok y_zad, skok FD, MSE zadanie, MSE FD]:')
print(np.array(extra))
