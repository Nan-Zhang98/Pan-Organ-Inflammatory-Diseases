#### import packages -------------------------------
import numpy as np
import pandas as pd
import scanpy as sc
import omicverse as ov
from omicverse.external import VIA

import matplotlib.pyplot as plt

#### global ----------------------------------------
momac_color = {
  "IL1B+NLRP3+_Macro":"#FA8072",
  "MMP3+CXCL8+_Macro":"#FFC5BF",
  "CCL13+_Complement-associated_Macro" : "#B6D0E2",
  "SPP1+CCL2+_Macro":"#ff9d5c",
  "CD1C+_DC-like_Macro":"#53738c"
}
#### load data ------------------------------------
Mac_10 = sc.read_h5ad("/storage/data/KAI/Pan_Inflammation/Homo/Mac2025/analysis/SCENIC/data/Mac_10.h5ad")


#### VIA analysis ---------------------------------
## preporcess

## parameters setting
adata = Mac_10
ncomps=30
knn=15
v0_random_seed=4
root_user = ['IL1B+NLRP3+_Macro'] #the index of a cell belonging to the HSC cell type
memory = 10
dataset = 'group'
use_rep = 'X_scVI'
clusters = 'cell_type_level3'
basis='X_umap'

## staVIA running

v0 = VIA.core.VIA(data=adata.obsm[use_rep][:, 0:ncomps],
             true_label=adata.obs[clusters],
             edgepruning_clustering_resolution=0.15, cluster_graph_pruning=0.15,
             knn=knn,  root_user=root_user, resolution_parameter=1.5,
             dataset=dataset, random_seed=v0_random_seed, memory=memory)#, do_compute_embedding=True, embedding_type='via-atlas')

v0.run_VIA()


fig, ax, ax1 = VIA.core.plot_piechart_viagraph_ov(adata,clusters='Organ',dpi=80,
                                             via_object=v0, ax_text=False,show_legend=False)
fig.set_size_inches(8,4)


adata.obs['pt_via']=v0.single_cell_pt_markov
ov.pl.embedding(adata,basis='X_umap',
                   color=['pt_via'],
                   frameon='small',cmap='Reds')


## Visual
adata.uns['cell_type_level3_colors'] = [
    momac_color[cell] for cell in adata.obs['cell_type_level3']
]



# streamplot colored by level3
fig, ax = VIA.core.via_streamplot_ov(adata,clusters='cell_type_level3',
                                     via_object=v0, embedding=adata.obsm['X_umap'], dpi=80,
                             density_grid=0.8, scatter_size=30,
                             scatter_alpha=0.3, linewidth=0.5,
                             labels=adata.obs['cell_type_level3'],
                             color_dict=momac_color,
                                     )
fig.set_size_inches(5,5)
fig.savefig("/storage/data/KAI/Pan_Inflammation/Homo/Mac2025/analysis/staVIA/result/streamplot_Mac2025_level3.svg")


# streamplot colored by pseudotime
v0.embedding = adata.obsm['X_umap']
fig, ax = VIA.core.via_streamplot_ov(adata,clusters='cell_type_level3',
                             via_object=v0,density_grid=0.8, scatter_size=30, color_scheme='time', linewidth=0.5,
                             min_mass = 1, cutoff_perc = 5, scatter_alpha=0.3, marker_edgewidth=0.1,
                             density_stream = 2, smooth_transition=1, smooth_grid=0.5,dpi=80,
                             cmap="plasma"
                                     )
fig.set_size_inches(5,5)
cbar = plt.colorbar(ax.collections[0], ax=ax)
cbar.set_label('Pseudotime', rotation=270, labelpad=15)  # 为颜色条加标签
cbar.set_ticks([0, 0.25, 0.5, 0.75, 1])  # 设置颜色条的刻度
cbar.set_ticklabels(['Start', '25%', '50%', '75%', 'End'])  # 自定义刻度标签

fig.savefig("/storage/data/KAI/Pan_Inflammation/Homo/Mac2025/analysis/staVIA/result/streamplot_Mac2025_pseudotime.svg")
