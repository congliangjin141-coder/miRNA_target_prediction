# ==============================================================================
# Fig4.py - Figure 4: power
# ==============================================================================
# Author : Congliang Jin
# Date   : 2026-10-23
#
# Description:
#   Power = 1 - Type II error = TP / (TP + FN); dashed line = 0.80.
#   Fig4A - all predictions of each tool, truth signal c = 1, 2, 3 (weak to strong)
#   Fig4B - truth c = 1, model signal = top 10%-100% of each tool's predictions
#
# Input:
#   analysis/results/refscore_metrics.csv
#
# Output (figure/Fig4/results/):
#   fig4.{png,pdf}
#
# Prerequisite:
#   analysis/03_compute_confusion_metrics.py (its output is tracked)
#
# Run from the repository root:
#   python figure/Fig4/Fig4.py
# ==============================================================================
import os

import pandas as pd
import matplotlib.pyplot as plt

plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans'],
                     'pdf.fonttype': 42, 'ps.fonttype': 42, 'axes.linewidth': 1.4})


METRICS = 'analysis/results/refscore_metrics.csv'
OUT_DIR = 'figure/Fig4/results'
TOOLS   = ['DIANA-microT', 'miRDB', 'mirDIP', 'miRmap', 'miRWalk', 'RNA22', 'TargetScan']
COLORS  = ['#E07B39', '#E8C030', '#7B5EA7', '#C94040', '#2AADAD', '#6AAF3D', '#A0345A']
MARKERS = ['o', 's', '^', 'D', 'v', '*', 'P']
# legend order: tools by number of predicted pairs, fewest to most (as in Fig. 2C)
LEGEND_ORDER = ['TargetScan', 'miRDB', 'miRmap', 'DIANA-microT', 'RNA22', 'miRWalk', 'mirDIP']
PCTS    = list(range(10, 101, 10))
TRUTH_CUTOFF   = 1
WEAK_TO_STRONG = [1, 2, 3]
POWER_BASELINE = 0.80


FIG_HEIGHT = 5.4
PANEL_WIDTH_TRUTH, PANEL_WIDTH_MODEL = 4.6, 7.0
LABEL_FS, TICK_FS, NOTE_FS, LEGEND_FS, LETTER_FS = 19, 16.5, 16, 16, 26
LINE_W, MARKER_S, STAR_S, MARKER_EDGE = 2.8, 10, 14, 2.0
RED = '#C03030'

X_TRUTH = 'Truth signal (weak to strong)'
X_MODEL = 'Model signal (strong to weak)'


m = pd.read_csv(METRICS)
m['power'] = 1 - m['type2']
BY_TRUTH = m[m['pct'] == 100]
BY_MODEL = m[m['truth_cutoff'] == TRUTH_CUTOFF]


def draw_lines(ax, data, x_col, y_col, x_pos, transform=lambda y: y):
    for tool, color, marker in zip(TOOLS, COLORS, MARKERS):
        y = transform(data[data['tool'] == tool].sort_values(x_col)[y_col].to_numpy())
        ax.plot(x_pos, y, color=color, marker=marker, linewidth=LINE_W,
                markersize=STAR_S if marker == '*' else MARKER_S,
                markerfacecolor='white', markeredgewidth=MARKER_EDGE, label=tool, zorder=3)


def style(ax, xlabel, ylabel):
    ax.set_xlabel(xlabel, fontsize=LABEL_FS, labelpad=6)
    ax.set_ylabel(ylabel, fontsize=LABEL_FS, labelpad=6)
    ax.tick_params(labelsize=TICK_FS, width=1.4, length=6)
    ax.spines[['top', 'right']].set_visible(False)
    ax.grid(True, which='major', linestyle='--', linewidth=0.6, color='#d4d4d4', zorder=0)
    ax.set_axisbelow(True)


def truth_positions():
    return [WEAK_TO_STRONG.index(c) for c in sorted(WEAK_TO_STRONG)]


def truth_axis(ax):
    ax.set_xticks(range(len(WEAK_TO_STRONG))); ax.set_xticklabels([str(c) for c in WEAK_TO_STRONG])
    ax.set_xlim(-0.25, len(WEAK_TO_STRONG) - 0.75)


def model_axis(ax, right=104):
    ax.set_xticks(PCTS); ax.set_xticklabels([str(p) for p in PCTS]); ax.set_xlim(6, right)


def baseline(ax, y, text, x, ha='left', above=True, offset=0.012):
    ax.axhline(y, color=RED, linestyle='--', linewidth=2.0, zorder=2)
    ax.text(x, y + offset if above else y - offset, text, color=RED, fontsize=NOTE_FS, ha=ha,
            va='bottom' if above else 'top')


def panel_power_truth(ax):
    draw_lines(ax, BY_TRUTH, 'truth_cutoff', 'power', truth_positions())
    baseline(ax, POWER_BASELINE, f'{POWER_BASELINE:.2f}', -0.22)
    truth_axis(ax); ax.set_ylim(-0.03, 1.04)
    style(ax, X_TRUTH, 'Power')


def panel_power_model(ax):
    draw_lines(ax, BY_MODEL, 'pct', 'power', PCTS)
    baseline(ax, POWER_BASELINE, f'{POWER_BASELINE:.2f}', 7, above=False)
    model_axis(ax); ax.set_ylim(-0.03, 1.04)
    style(ax, X_MODEL, 'Power')


def make_figure(panels, widths, out, letters='ABCDEFG'):
    fig, axes = plt.subplots(1, len(panels), figsize=(sum(widths), FIG_HEIGHT), facecolor='white',
                             gridspec_kw={'width_ratios': widths})
    for ax, panel, letter in zip(axes, panels, letters):
        panel(ax)
        ax.text(-0.2, 1.03, letter, transform=ax.transAxes, fontsize=LETTER_FS, fontweight='bold',
                va='bottom', ha='left')
    handles, labels = axes[0].get_legend_handles_labels()
    handles = [handles[labels.index(t)] for t in LEGEND_ORDER]; labels = LEGEND_ORDER
    fig.legend(handles, labels, loc='lower center', ncol=4, fontsize=LEGEND_FS, frameon=False,
               handlelength=2.0, handletextpad=0.5, columnspacing=1.4, bbox_to_anchor=(0.5, 0.0))
    fig.tight_layout(rect=(0, 0.15, 1, 1), w_pad=1.5)
    os.makedirs(OUT_DIR, exist_ok=True)
    fig.savefig(os.path.join(OUT_DIR, f'{out}.png'), dpi=300)
    fig.savefig(os.path.join(OUT_DIR, f'{out}.pdf'))
    plt.close(fig)
    print('saved', os.path.join(OUT_DIR, out) + '.png / .pdf')


make_figure([panel_power_truth, panel_power_model], [PANEL_WIDTH_TRUTH, PANEL_WIDTH_MODEL], 'fig4')
