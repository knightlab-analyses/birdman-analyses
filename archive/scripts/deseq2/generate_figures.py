#!/usr/bin/env python
"""generate figure 5: birdman vs deseq2 concordance (scatter + forest plot)."""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from figure_style import (
    set_style, PALETTE_TWO, PALETTE_QUALITATIVE,
    add_panel_label, add_zero_line, save_figure,
)

set_style()
os.makedirs('figures', exist_ok=True)

# load data
results = pd.read_csv('deseq2_table/comparison_results.tsv', sep='\t', index_col=0)
results.index.name = 'Feature'
results['deseq2_ln'] = results['deseq2_lfc']

spearman_r, _ = stats.spearmanr(results['birdman_lfc'], results['deseq2_lfc'])
pearson_r, _ = stats.pearsonr(results['birdman_lfc'], results['deseq2_lfc'])
print(f"spearman r = {spearman_r:.3f}, pearson r = {pearson_r:.3f}")

# ci error bars
results['lower_err'] = results['birdman_lfc'] - results['birdman_lfc_q025']
results['upper_err'] = results['birdman_lfc_q975'] - results['birdman_lfc']
results['credible'] = np.where(
    (results['birdman_lfc_q025'] > 0) | (results['birdman_lfc_q975'] < 0), 'yes', 'no')
print(f"credible: {(results['credible'] == 'yes').sum()} / {len(results)}")

# colors
scatter_color = PALETTE_QUALITATIVE[0]
birdman_color = PALETTE_TWO['blue_red'][0]
deseq2_color = PALETTE_TWO['blue_red'][1]

# select top/bottom 5 features by deseq2 lfc for forest plot
sorted_by_lfc = results.sort_values('deseq2_lfc')
all_selected = sorted_by_lfc.head(5).index.tolist() + sorted_by_lfc.tail(5).index.tolist()
selected_data = results.loc[all_selected].sort_values('deseq2_lfc')

transition_idx = None
for i, (_, row) in enumerate(selected_data.iterrows()):
    if row['deseq2_lfc'] > 0:
        transition_idx = i
        break

panel_size = 2.6
fig, axes = plt.subplots(1, 2, figsize=(panel_size * 2 + 0.4, panel_size))

# panel a: scatter
ax = axes[0]
ax.scatter(results['deseq2_ln'], results['birdman_lfc'],
           s=18, alpha=0.75, c=scatter_color, edgecolors='white', linewidth=0.3, zorder=5)
lims = [min(results['deseq2_ln'].min(), results['birdman_lfc'].min()) - 0.5,
        max(results['deseq2_ln'].max(), results['birdman_lfc'].max()) + 0.5]
ax.plot(lims, lims, '--', color='#888888', linewidth=0.5, alpha=0.6, zorder=0)
ax.set_xlim(lims)
ax.set_ylim(lims)
ax.text(0.05, 0.95, f'Spearman \u03c1 = {spearman_r:.2f}',
        transform=ax.transAxes, fontsize=6, fontweight='bold', va='top')
ax.set_xlabel('DESeq2 log$_2$ fold change')
ax.set_ylabel('BIRDMAn log$_2$ fold change')
ax.set_aspect('equal')
add_panel_label(ax, 'a')

# panel b: forest plot
ax = axes[1]
y_positions = np.arange(len(selected_data))
offset = 0.15

for i, (feature, row) in enumerate(selected_data.iterrows()):
    y = y_positions[i]
    ax.scatter(row['deseq2_ln'], y + offset, marker='s', s=22,
               color=deseq2_color, edgecolors='white', linewidth=0.3, zorder=5,
               label='DESeq2' if i == 0 else '')
    ax.errorbar(row['birdman_lfc'], y - offset,
                xerr=[[row['lower_err']], [row['upper_err']]],
                fmt='o', markersize=4, color=birdman_color,
                markeredgecolor='white', markeredgewidth=0.3,
                capsize=1.5, capthick=0.5, elinewidth=0.5, zorder=5,
                label='BIRDMAn (95% CI)' if i == 0 else '')

add_zero_line(ax, orientation='vertical')

if transition_idx:
    ax.axhline(y=transition_idx - 0.5, color='#cccccc', linestyle='-', linewidth=0.4, zorder=0)

ax.set_yticks(y_positions)
display_names = [f.replace('_', ' ').replace('Lactobacillus ', 'L. ') for f in selected_data.index]
ax.set_yticklabels(display_names)
ax.set_xlabel('Log$_2$ fold change (BV vs healthy)')
ax.set_ylabel('')
ax.legend(loc='lower right', fontsize=5, handlelength=1.2,
          handletextpad=0.3, borderpad=0.2, labelspacing=0.3)

# directional arrows
arrow_y = -0.18
ax.annotate('', xy=(0.08, arrow_y), xytext=(0.35, arrow_y),
            xycoords='axes fraction', textcoords='axes fraction',
            arrowprops=dict(arrowstyle='->', color='#666666', lw=0.75),
            annotation_clip=False)
ax.text(0.21, arrow_y - 0.06, 'Health-\nassociated',
        ha='center', va='top', fontsize=5, color='#666666',
        linespacing=0.85, transform=ax.transAxes)
ax.annotate('', xy=(0.92, arrow_y), xytext=(0.65, arrow_y),
            xycoords='axes fraction', textcoords='axes fraction',
            arrowprops=dict(arrowstyle='->', color='#666666', lw=0.75),
            annotation_clip=False)
ax.text(0.79, arrow_y - 0.06, 'BV-\nassociated',
        ha='center', va='top', fontsize=5, color='#666666',
        linespacing=0.85, transform=ax.transAxes)

add_panel_label(ax, 'b')

plt.tight_layout()
plt.subplots_adjust(bottom=0.18)
save_figure(fig, 'figures/figure5', formats=['pdf', 'png'], dpi=450)
plt.close()

