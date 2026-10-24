# ==============================================================================
# Fig2.py - Figure 2: reference datasets and tool coverage
# ==============================================================================
# Author : Congliang Jin
# Date   : 2026-10-23
#
# Description:
#   Fig2A - pairs, miRNAs and genes of miRTarBase, TarBase, the iPSC dataset and their union
#   Fig2B - Venn diagram of the three references
#   Fig2C - pairs, miRNAs and genes predicted by each tool
#   Fig2D - genes and miRNAs covered by the top 10%-100% of each tool, with the
#           reference score 3 / 2 / 1 pairs
#   Panels A-C and the D reference points use counts written below (from the reference
#   table and the tool files); the D tool curves are read from fig2_curve_data.json, or
#   recomputed from data/predictions/ when the json is absent.
#
# Input:
#   figure/Fig2/results/fig2_curve_data.json  (or data/predictions/)
#
# Output (figure/Fig2/results/):
#   fig2.{png,pdf}
#
# Memory:
#   large only when the curves are recomputed
#
# Run from the repository root:
#   python figure/Fig2/Fig2.py
# ==============================================================================
import json, os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Circle
from matplotlib.gridspec import GridSpecFromSubplotSpec
from matplotlib.ticker import FixedLocator, FixedFormatter, NullFormatter, NullLocator

OUT_DIR = 'figure/Fig2/results'
os.makedirs(OUT_DIR, exist_ok=True)

plt.rcParams.update({
    'font.family'    : 'sans-serif',
    'font.sans-serif': ['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans'],
    'pdf.fonttype'   : 42,
    'ps.fonttype'    : 42,
    'axes.linewidth' : 0.8,
})
FS_AXLABEL = 14
FS_TICK    = 13
FS_LETTER  = 19
FS_VTITLE  = 16
FS_LEG     = 12

def draw_bars_panel(fig, subspec, *, labels, leg_labels, colors, edges,
                    pairs, mirnas, genes,
                    pairs_ylim, pairs_yticks,
                    mir_ylim, mir_yticks, mir_yticklab,
                    gene_ylim, gene_yticks, gene_yticklab,
                    bar_w, pairs_mult, pairs_fs, pairs_rot,
                    mir_off, mir_fs, gene_off, gene_fs, leg_ncol,
                    xtick_fs=None, xtick_rot=0):
    N = len(labels)
    x = np.arange(N)

    pgs = GridSpecFromSubplotSpec(2, 1, subplot_spec=subspec,
                                  height_ratios=[10, 1.4], hspace=0.42)
    cgs = GridSpecFromSubplotSpec(2, 2, subplot_spec=pgs[0, 0],
                                  width_ratios=[1, 1.0],
                                  height_ratios=[2.2, 1.0],
                                  hspace=0.07, wspace=0.34)
    ax1    = fig.add_subplot(cgs[:, 0])
    ax_top = fig.add_subplot(cgs[0, 1])
    ax_bot = fig.add_subplot(cgs[1, 1])

    b = ax1.bar(x, pairs, width=bar_w, color=colors, edgecolor=edges,
                linewidth=0.8, zorder=3)
    for bar, v in zip(b, pairs):
        _ha = 'left' if pairs_rot > 30 else 'center'
        ax1.text(bar.get_x()+bar.get_width()/2, v*pairs_mult, f'{v:,}',
                 ha=_ha, va='bottom', fontsize=pairs_fs,
                 color='#222', rotation=pairs_rot, fontweight='bold')
    ax1.set_yscale('log'); ax1.set_ylim(*pairs_ylim)
    ax1.set_yticks(pairs_yticks)
    ax1.set_yticklabels([f'$10^{{{int(np.log10(t))}}}$' for t in pairs_yticks],
                        fontsize=FS_TICK)
    _xfs = xtick_fs if xtick_fs is not None else FS_TICK
    ax1.set_xticks(x); ax1.set_xticklabels(labels, fontsize=_xfs,
                                            rotation=xtick_rot,
                                            ha='right' if xtick_rot > 0 else 'center')
    ax1.set_ylabel('# of pairs (log scale)', fontsize=FS_AXLABEL)
    ax1.spines[['top', 'right']].set_visible(False)
    ax1.yaxis.grid(True, which='both', linestyle='--', linewidth=0.5,
                   color='#ddd', zorder=0); ax1.set_axisbelow(True)

    bt = ax_top.bar(x, genes, width=bar_w, color=colors, edgecolor=edges,
                    linewidth=0.8, zorder=3)
    for bar, v in zip(bt, genes):
        ax_top.text(bar.get_x()+bar.get_width()/2, v+gene_off, f'{v:,}',
                    ha='center', va='bottom', fontsize=gene_fs, color='#222', fontweight='bold')
    ax_top.set_xlim(-0.5, N-0.5); ax_top.set_ylim(*gene_ylim)
    ax_top.set_yticks(gene_yticks); ax_top.set_yticklabels(gene_yticklab, fontsize=FS_TICK)
    ax_top.set_ylabel('# of genes', fontsize=FS_AXLABEL, labelpad=6)
    ax_top.set_xticks(x); ax_top.set_xticklabels([])
    ax_top.spines[['top', 'right', 'bottom']].set_visible(False)
    ax_top.tick_params(bottom=False)
    ax_top.yaxis.grid(True, linestyle='--', linewidth=0.5, color='#ddd', zorder=0)
    ax_top.set_axisbelow(True); ax_top.axhline(0, color='#555', linewidth=1.0)

    bb = ax_bot.bar(x, mirnas, width=bar_w, color=colors, edgecolor=edges,
                    linewidth=0.8, zorder=3)
    for bar, v, ec in zip(bb, mirnas, edges):
        ax_bot.text(bar.get_x()+bar.get_width()/2, v+mir_off, f'{v:,}',
                    ha='center', va='top', fontsize=mir_fs, color='#222',
                    fontweight='bold', zorder=5)
    ax_bot.set_xlim(-0.5, N-0.5); ax_bot.set_ylim(*mir_ylim); ax_bot.invert_yaxis()
    _mir_labs = ['' if t == 0 else l for t, l in zip(mir_yticks, mir_yticklab)]
    ax_bot.set_yticks(mir_yticks); ax_bot.set_yticklabels(_mir_labs, fontsize=FS_TICK)
    ax_bot.set_ylabel('# of miRNA', fontsize=FS_AXLABEL, labelpad=6)
    ax_bot.set_xticks(x); ax_bot.set_xticklabels(labels, fontsize=_xfs,
                                                  rotation=xtick_rot,
                                                  ha='right' if xtick_rot > 0 else 'center')
    ax_bot.spines[['top', 'right']].set_visible(False)
    ax_bot.yaxis.grid(True, linestyle='--', linewidth=0.5, color='#ddd', zorder=0)
    ax_bot.set_axisbelow(True)

    legax = fig.add_subplot(pgs[1, 0]); legax.axis('off')
    patches = [mpatches.Patch(facecolor=c, edgecolor=e, linewidth=0.8, label=l)
               for c, e, l in zip(colors, edges, leg_labels)]
    legax.legend(handles=patches, loc='center', ncol=leg_ncol, fontsize=FS_LEG,
                 frameon=True, framealpha=0.9, edgecolor='#ccc',
                 handlelength=1.4, columnspacing=1.0, borderpad=0.6)


def draw_venn(fig, subspec, real_center=(0.06, -0.80), max_fs=13):
    C_MTB, E_MTB = '#85B7EB', '#185FA5'
    C_TB,  E_TB  = '#5DCAA5', '#0F6E56'
    C_RD,  E_RD  = '#F4A078', '#C05020'
    PINK = '#F4A7BB'
    n_mtb, n_tb, n_rd = 951_010, 1_271_467, 256_876
    base = np.sqrt(n_tb)
    r_mtb = 1.55*np.sqrt(n_mtb)/base
    r_tb  = 1.55
    r_rd  = 1.55*np.sqrt(n_rd)/base
    cx1, cy1 = -0.40, 0.50
    cx2, cy2 =  0.90, 0.50
    cx3, cy3 = real_center
    ax = fig.add_subplot(subspec)
    ax.set_facecolor('white'); ax.set_aspect('equal'); ax.axis('off')
    ax.set_xlim(-1.9, 2.6); ax.set_ylim(-2.0, 2.55)
    for cx, cy, r in [(cx1, cy1, r_mtb), (cx2, cy2, r_tb), (cx3, cy3, r_rd)]:
        ax.add_patch(Circle((cx, cy), r, fc='white', ec='none', zorder=1))
    for cx, cy, r, fc in [(cx1, cy1, r_mtb, C_MTB), (cx2, cy2, r_tb, C_TB),
                          (cx3, cy3, r_rd, C_RD)]:
        ax.add_patch(Circle((cx, cy), r, fc=fc, ec='none', alpha=0.40, zorder=2))
    def pink_clip(ca, ra, cb, rb, z=3, color=None):
        fc = color if color else PINK
        p = Circle(ca, ra, fc=fc, ec='none', alpha=0.85, zorder=z)
        ax.add_patch(p); p.set_clip_path(Circle(cb, rb, transform=ax.transData))
    pink_clip((cx1, cy1), r_mtb, (cx2, cy2), r_tb, color='#A8C4DC')
    pink_clip((cx1, cy1), r_mtb, (cx3, cy3), r_rd)
    pink_clip((cx2, cy2), r_tb,  (cx3, cy3), r_rd)

    for cx, cy, r, ec in [(cx1, cy1, r_mtb, E_MTB), (cx2, cy2, r_tb, E_TB),
                          (cx3, cy3, r_rd, E_RD)]:
        ax.add_patch(Circle((cx, cy), r, fc='none', ec=ec, lw=1.8, zorder=5))

    region_labels = [((True,  False, False), '862,478'),
                     ((False, True,  False), '1,164,907'),
                     ((True,  True,  False), '80,836'),
                     ((True,  True,  True),  '1,683'),
                     ((True,  False, True),  '6,013'),
                     ((False, True,  True),  '24,041'),
                     ((False, False, True),  '225,139')]
    step, clearance = 0.004, 0.022
    gx, gy = np.meshgrid(np.arange(-2.4, 3.1, step), np.arange(-2.2, 2.6, step))
    circles = [(cx1, cy1, r_mtb), (cx2, cy2, r_tb), (cx3, cy3, r_rd)]
    signed = [r - np.hypot(gx - cx, gy - cy) for cx, cy, r in circles]

    def place(flags, w, h):
        mask = np.ones(gx.shape, dtype=bool)
        depth = np.full(gx.shape, np.inf)
        for d, inside in zip(signed, flags):
            mask &= (d > clearance) if inside else (d < -clearance)
            depth = np.minimum(depth, np.abs(d))
        bw, bh = int(np.ceil(w / step)), int(np.ceil(h / step))
        I = np.pad(mask.astype(np.int32).cumsum(0).cumsum(1), ((1, 0), (1, 0)))
        ok = (I[bh:, bw:] - I[:-bh, bw:] - I[bh:, :-bw] + I[:-bh, :-bw]) == bw * bh
        if not ok.any():
            return None
        centre_depth = depth[bh // 2: bh // 2 + ok.shape[0], bw // 2: bw // 2 + ok.shape[1]]
        i, j = np.unravel_index(np.argmax(np.where(ok, centre_depth, -1.0)), ok.shape)
        return gx[0, j] + w / 2, gy[i, 0] + h / 2

    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    to_data = ax.transData.inverted()
    for fs in np.arange(max_fs, 7.5, -0.5):
        spots = []
        for flags, t in region_labels:
            probe = ax.text(0, 0, t, fontsize=fs, fontweight='bold')
            box = probe.get_window_extent(renderer).transformed(to_data)
            probe.remove()
            spot = place(flags, box.width, box.height)
            if spot is None:
                break
            spots.append(spot)
        else:
            break
    else:
        raise RuntimeError('A region is too small to hold its label; move the iPSC-dataset circle.')
    for (x, y), (_, t) in zip(spots, region_labels):
        ax.text(x, y, t, ha='center', va='center', fontsize=fs,
                color='#111', fontweight='bold', zorder=7)
    print(f'Venn region labels drawn at {fs:g} pt')

    ax.text(cx1-0.25, cy1+r_mtb+0.15, 'miRTarBase', fontsize=FS_VTITLE,
            ha='center', va='bottom', color='#111', fontweight='bold')
    ax.text(cx2+0.40, cy2+r_tb+0.15, 'TarBase', fontsize=FS_VTITLE,
            ha='center', va='bottom', color='#111', fontweight='bold')
    ax.text(cx3, cy3-r_rd-0.18, 'iPSC dataset', fontsize=FS_VTITLE,
            ha='center', va='top', color='#111', fontweight='bold')


CURVE_COLORS  = ['#E07B39', '#E8C030', '#7B5EA7', '#C94040', '#2AADAD', '#6AAF3D', '#A0345A']
CURVE_MARKERS = ['o', 's', '^', 'D', 'v', '*', 'P']
DISPLAY_NAMES = ['DIANA-microT', 'miRDB', 'mirDIP',
                 'miRmap', 'miRWalk', 'RNA22', 'TargetScan']


CURVE_CACHE = f'{OUT_DIR}/fig2_curve_data.json'


REF_SCORE_POINTS = {
    3: (1_122,  193,   '*', 16),
    2: (12_941, 1_025, 'D', 9),
    1: (19_811, 2_011, 's', 9),
}
CURVE_X_TICKS = [5_000, 10_000, 15_000, 20_000, 25_000]
CURVE_Y_TICKS = [200, 500, 1_000, 1_500, 2_000, 2_500]
CURVE_X_LIM = (850, 33_000)
CURVE_Y_LIM = (140, 3_400)


def load_curve_data():
    import pandas as pd
    base = 'data/predictions'
    configs = {
        'DIANA-microT': {'path':f'{base}/diana_microT.csv',   'mi':'mirna_name','gene':'gene_name','score':'score','higher_better':True},
        'miRDB':        {'path':f'{base}/miRDB v6.0.csv',     'mi':'mirna_name','gene':'gene_name','score':'score','higher_better':True},
        'mirDIP':       {'path':f'{base}/mirDIP v5.2.csv',    'mi':'mirna_name','gene':'gene_name','score':'score','higher_better':True},
        'miRmap':       {'path':f'{base}/mirmap.csv',          'mi':'mirna_name','gene':'gene_name','score':'score','higher_better':False},
        'miRWalk':      {'path':f'{base}/mirwalk.csv',         'mi':'mirna_name','gene':'gene_name','score':'score','higher_better':True},
        'RNA22':        {'path':f'{base}/rna22.csv',           'mi':'mirna_name','gene':'gene_name','score':'score','higher_better':False,'chunksize':500_000},
        'TargetScan':   {'path':f'{base}/TargetScan v8.0.csv','mi':'mirna_name','gene':'gene_name','score':'score','higher_better':False},
    }
    pcts = list(range(10, 101, 10))
    cum = {}
    print("Reading data (RNA22 is slow)...")
    for name, cfg in configs.items():
        print(f"  {name} ...")
        mi, ge, sc, hb = cfg['mi'], cfg['gene'], cfg['score'], cfg['higher_better']
        if cfg.get('chunksize'):
            rows = []
            for chunk in pd.read_csv(cfg['path'], usecols=[mi,ge,sc],
                                     chunksize=cfg['chunksize'], low_memory=False):
                rows.append(chunk.dropna(subset=[mi,ge,sc]))
            df = pd.concat(rows, ignore_index=True)
        else:
            df = pd.read_csv(cfg['path'], usecols=[mi,ge,sc], low_memory=False)
            df = df.dropna(subset=[mi,ge,sc])
        df[sc] = pd.to_numeric(df[sc], errors='coerce')
        df = df.dropna(subset=[sc])

        # several rows of one pair: keep the strongest score, as in the reference table
        df = df.groupby([mi, ge], as_index=False, sort=False)[sc].agg('max' if hb else 'min')
        g_list, m_list = [], []
        for p in pcts:
            if hb:
                thr = df[sc].quantile(1-p/100); sub = df[df[sc]>=thr]
            else:
                thr = df[sc].quantile(p/100);   sub = df[df[sc]<=thr]
            g_list.append(sub[ge].nunique()); m_list.append(sub[mi].nunique())
        cum[name] = {'genes': g_list, 'mirna': m_list}
    return cum

def draw_curve(fig, subspec, cum_data):
    inner = GridSpecFromSubplotSpec(2, 1, subplot_spec=subspec, height_ratios=[10, 1.0], hspace=0)
    ax = fig.add_subplot(inner[0, 0])
    # tools in the same order as panel C (fewest to most predicted pairs); colour, marker
    # and drawing layer stay with each tool
    for name in c_tools:
        i = DISPLAY_NAMES.index(name)
        cd = cum_data[name]
        ax.plot(cd['genes'], cd['mirna'], color=CURVE_COLORS[i],
                marker=CURVE_MARKERS[i], linewidth=1.8, markersize=5,
                label=name, markerfacecolor='white', markeredgewidth=1.4, zorder=3+i)

    for score, (n_genes, n_mirna, marker, size) in REF_SCORE_POINTS.items():
        ax.plot(n_genes, n_mirna, linestyle='none', marker=marker, markersize=size,
                markerfacecolor='#111', markeredgecolor='#111', zorder=20,
                label=f'Ref score = {score}')
        ax.annotate(f'score = {score}', (n_genes, n_mirna), xytext=(0, -14),
                    textcoords='offset points', ha='center', va='top',
                    fontsize=FS_LEG, color='#111')

    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlim(*CURVE_X_LIM); ax.set_ylim(*CURVE_Y_LIM)
    ax.xaxis.set_major_locator(FixedLocator(CURVE_X_TICKS))
    ax.xaxis.set_major_formatter(FixedFormatter([f'{t:,}' for t in CURVE_X_TICKS]))
    ax.yaxis.set_major_locator(FixedLocator(CURVE_Y_TICKS))
    ax.yaxis.set_major_formatter(FixedFormatter([f'{t:,}' for t in CURVE_Y_TICKS]))
    for axis in (ax.xaxis, ax.yaxis):
        axis.set_minor_locator(NullLocator()); axis.set_minor_formatter(NullFormatter())
    ax.tick_params(labelsize=FS_TICK)
    plt.setp(ax.get_xticklabels(), rotation=60, ha='right', rotation_mode='anchor')
    ax.set_xlabel('# of genes (log scale)', fontsize=FS_AXLABEL)
    ax.set_ylabel('# of miRNA (log scale)', fontsize=FS_AXLABEL)
    ax.spines[['top', 'right']].set_visible(False)
    ax.grid(True, linestyle='--', linewidth=0.35, color='#e0e0e0', zorder=0)
    ax.set_axisbelow(True)
    ax.legend(loc='upper left', fontsize=FS_LEG, frameon=True)

def add_letter(fig, subspec, letter):
    pos = subspec.get_position(fig)
    fig.text(pos.x0 - 0.012, pos.y1 + 0.012, letter,
             fontsize=FS_LETTER, fontweight='bold', va='top', ha='left')

fig = plt.figure(figsize=(20, 10.4), facecolor='white')
outer = fig.add_gridspec(2, 2, width_ratios=[1.75, 1], height_ratios=[1, 1],
                         hspace=0.08, wspace=0.12,
                         left=0.045, right=0.98, top=0.95, bottom=0.055)

draw_bars_panel(
    fig, outer[0, 0],
    labels=['miRTarBase', 'TarBase', 'iPSC dataset', 'Union'],
    leg_labels=['miRTarBase', 'TarBase', 'iPSC dataset', 'Union'],
    colors=['#85B7EB', '#5DCAA5', '#F4A078', '#B9A3E3'],
    edges=['#185FA5', '#0F6E56', '#C05020', '#5B3E99'],
    pairs=[951_010, 1_271_467, 256_876, 2_365_097],
    mirnas=[1964, 1132, 1265, 2011],
    genes=[16250, 17844, 13844, 19877],
    pairs_ylim=(5e4, 2.5e7), pairs_yticks=[1e5, 1e6, 1e7],
    mir_ylim=(0, 3000), mir_yticks=[0, 1000, 2000], mir_yticklab=['0', '1,000', '2,000'],
    gene_ylim=(0, 24000), gene_yticks=[0, 5000, 10000, 15000, 20000],
    gene_yticklab=['0', '5,000', '10,000', '15,000', '20,000'],
    bar_w=0.55, pairs_mult=1.06, pairs_fs=11, pairs_rot=0,
    mir_off=40, mir_fs=11, gene_off=600, gene_fs=11, leg_ncol=4,
    xtick_fs=11.5)

draw_venn(fig, outer[0, 1])

# Panel C: tools ordered by the number of predicted pairs, fewest to most
PANEL_C = [  # tool, colour, edge colour, pairs, miRNAs, genes
    ('DIANA-microT', '#E07B39', '#A04A10', 24_787_164, 2654, 19404),
    ('miRDB',        '#E8C030', '#B08000', 1_607_873,  2656, 17454),
    ('mirDIP',       '#7B5EA7', '#4A3070', 45_783_772, 2656, 27936),
    ('miRmap',       '#C94040', '#8A1010', 16_783_062, 2656, 18909),
    ('miRWalk',      '#2AADAD', '#107070', 31_807_243, 2656, 19342),
    ('RNA22',        '#6AAF3D', '#3A7010', 30_435_856, 2588, 18669),
    ('TargetScan',   '#A0345A', '#601030', 202_584,    321,  13035),
]
PANEL_C.sort(key=lambda row: row[3])
c_tools, c_colors, c_edges, c_pairs, c_mirnas, c_genes = (list(col) for col in zip(*PANEL_C))

draw_bars_panel(
    fig, outer[1, 0],
    labels=c_tools,
    leg_labels=c_tools,
    colors=c_colors,
    edges=c_edges,
    pairs=c_pairs,
    mirnas=c_mirnas,
    genes=c_genes,
    pairs_ylim=(1e5, 2e8), pairs_yticks=[1e5, 1e6, 1e7, 1e8],
    mir_ylim=(0, 3600), mir_yticks=[0, 1000, 2000, 3000],
    mir_yticklab=['0', '1,000', '2,000', '3,000'],
    gene_ylim=(0, 36000), gene_yticks=[0, 10000, 20000, 30000],
    gene_yticklab=['0', '10,000', '20,000', '30,000'],
    bar_w=0.6, pairs_mult=1.06, pairs_fs=9.5, pairs_rot=65,
    mir_off=30, mir_fs=9, gene_off=700, gene_fs=10, leg_ncol=4,
    xtick_fs=10.5, xtick_rot=25)

if os.path.exists(CURVE_CACHE):
    with open(CURVE_CACHE) as fh:
        cum_data = json.load(fh)
    print(f'Tool curves read from the cache {CURVE_CACHE}')
else:
    cum_data = load_curve_data()
    with open(CURVE_CACHE, 'w') as fh:
        json.dump({k: {a: [int(x) for x in b] for a, b in v.items()} for k, v in cum_data.items()}, fh, indent=1)
    print(f'Tool curves saved to {CURVE_CACHE}')
draw_curve(fig, outer[1, 1], cum_data)

for cell, letter in [(outer[0, 0], 'A'), (outer[0, 1], 'B'),
                     (outer[1, 0], 'C'), (outer[1, 1], 'D')]:
    add_letter(fig, cell, letter)

out = f'{OUT_DIR}/fig2'
fig.savefig(out+'.pdf', dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(out+'.png', dpi=300, bbox_inches='tight', facecolor='white')
print(f'Saved {out}.pdf and {out}.png')
