import io
import os
import re
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
from psynlig import plot_correlation_heatmap
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = 'Helvetica'
mpl.rcParams['pdf.use14corefonts'] = True
mpl.rcParams['pdf.fonttype'] = 42
%config InlineBackend.figure_format = 'retina'
from matplotlib.colors import LinearSegmentedColormap

## Figure 5c <Dendograms>

ct = pd.read_csv('../01_data/TNBC_Times_MTX.csv')
ct = ct.drop(columns=['Mutation'])

g = sns.clustermap(ct.corr(),
                   vmin=-1, vmax=1,
                   cmap=custom_cmap,
                   dendrogram_ratio=(.1, .1),
                   cbar_pos=(.02, .32, .03, .2),
                   linewidths=.5, figsize=(10, 8))
plt.subplots_adjust(left = 0.5)
for label in g.ax_heatmap.get_xticklabels():
    label.set_rotation(90)
do = '../03_outs/Dendos_Time.pdf'
plt.savefig(do, bbox_inches='tight', pad_inches=1)
plt.show()