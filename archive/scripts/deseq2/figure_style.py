import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

WIDTHS = {
    'single': 3.50,   # 89 mm
    '1.5col': 5.04,   # 128 mm
    'double': 7.09,   # 180 mm
}
MAX_HEIGHT = 6.69  # 170 mm

PALETTE_QUALITATIVE = [
    '#332288', '#88CCEE', '#44AA99', '#117733', '#999933',
    '#DDCC77', '#CC6677', '#882255', '#AA4499',
]

PALETTE_TWO = {
    'blue_orange': ['#0077BB', '#EE7733'],
    'blue_red': ['#4477AA', '#CC3311'],
    'teal_wine': ['#44AA99', '#882255'],
    'blue_gold': ['#2E86AB', '#D4A03C'],
}

PALETTE_GRAY = ['#333333', '#666666', '#999999', '#CCCCCC']



def set_style(context='paper'):
    scale = 1.0 if context == 'paper' else 1.4
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'font.size': 7 * scale,
        'axes.titlesize': 8 * scale,
        'axes.labelsize': 7 * scale,
        'xtick.labelsize': 6 * scale,
        'ytick.labelsize': 6 * scale,
        'legend.fontsize': 5 * scale,
        'legend.title_fontsize': 6 * scale,
        'figure.titlesize': 8 * scale,
        'axes.linewidth': 0.5,
        'xtick.major.width': 0.5,
        'ytick.major.width': 0.5,
        'xtick.minor.width': 0.3,
        'ytick.minor.width': 0.3,
        'lines.linewidth': 0.75,
        'patch.linewidth': 0.5,
        'hatch.linewidth': 0.5,
        'grid.linewidth': 0.3,
        'xtick.major.size': 3,
        'ytick.major.size': 3,
        'xtick.minor.size': 1.5,
        'ytick.minor.size': 1.5,
        'xtick.direction': 'out',
        'ytick.direction': 'out',
        'xtick.major.pad': 2,
        'ytick.major.pad': 2,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'axes.grid': False,
        'axes.axisbelow': True,
        'legend.frameon': False,
        'legend.borderpad': 0.3,
        'legend.handlelength': 1.2,
        'legend.handletextpad': 0.4,
        'legend.labelspacing': 0.3,
        'legend.columnspacing': 1.0,
        'figure.dpi': 150,
        'savefig.dpi': 450,
        'savefig.format': 'pdf',
        'savefig.bbox': 'tight',
        'savefig.pad_inches': 0.02,
        'savefig.transparent': False,
        'figure.facecolor': 'white',
        'axes.facecolor': 'white',
        'figure.constrained_layout.use': True,
        'figure.constrained_layout.h_pad': 0.02,
        'figure.constrained_layout.w_pad': 0.02,
        'figure.autolayout': False,
        'scatter.edgecolors': 'white',
        'errorbar.capsize': 2,
        'mathtext.fontset': 'dejavusans',
        'image.interpolation': 'nearest',
        'image.cmap': 'viridis',
    })
    plt.rcParams['axes.prop_cycle'] = plt.cycler(color=PALETTE_QUALITATIVE)


def create_figure(width='single', height=None, nrows=1, ncols=1,
                  aspect_ratio=None, **kwargs):
    if isinstance(width, str):
        w = WIDTHS.get(width, WIDTHS['single'])
    else:
        w = width
    if height is None:
        if aspect_ratio is None:
            aspect_ratio = 0.618
        height = w * aspect_ratio * (nrows / ncols) ** 0.5
        height = min(height, MAX_HEIGHT)
    return plt.subplots(nrows, ncols, figsize=(w, height), **kwargs)


def add_panel_label(ax, label, x=-0.12, y=1.08, **kwargs):
    defaults = {'fontsize': 8, 'fontweight': 'bold', 'va': 'top',
                'ha': 'right', 'transform': ax.transAxes}
    defaults.update(kwargs)
    ax.text(x, y, label, **defaults)


def add_panel_labels(axes, labels=None, **kwargs):
    import string
    if hasattr(axes, 'flatten'):
        axes = axes.flatten()
    if labels is None:
        labels = list(string.ascii_lowercase[:len(axes)])
    for ax, label in zip(axes, labels):
        add_panel_label(ax, label, **kwargs)


def add_zero_line(ax, orientation='horizontal', **kwargs):
    defaults = {'color': 'gray', 'linestyle': '--', 'linewidth': 0.5,
                'alpha': 0.5, 'zorder': 0}
    defaults.update(kwargs)
    if orientation == 'horizontal':
        ax.axhline(y=0, **defaults)
    else:
        ax.axvline(x=0, **defaults)


def add_identity_line(ax, **kwargs):
    defaults = {'color': 'gray', 'linestyle': '--', 'linewidth': 0.5,
                'alpha': 0.7, 'zorder': 0}
    defaults.update(kwargs)
    lims = [min(ax.get_xlim()[0], ax.get_ylim()[0]),
            max(ax.get_xlim()[1], ax.get_ylim()[1])]
    ax.plot(lims, lims, **defaults)
    ax.set_xlim(lims)
    ax.set_ylim(lims)


def save_figure(fig, filename, formats=None, dpi=450):
    if formats is None:
        formats = ['pdf', 'png']
    for fmt in formats:
        save_dpi = dpi if fmt in ['png', 'jpg', 'jpeg', 'tiff'] else None
        filepath = f'{filename}.{fmt}'
        fig.savefig(filepath, dpi=save_dpi, bbox_inches='tight',
                    pad_inches=0.02, format=fmt, facecolor='white')


def get_palette(n_colors=None, palette_type='qualitative'):
    palettes = {'qualitative': PALETTE_QUALITATIVE, 'gray': PALETTE_GRAY, **PALETTE_TWO}
    colors = palettes.get(palette_type, PALETTE_QUALITATIVE)
    if n_colors is not None:
        colors = colors[:n_colors]
    return colors

get_palette_compat = get_palette
