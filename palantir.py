### Module ----------------------------
import palantir
import scanpy as sc
import pandas as pd
import os

import matplotlib
import matplotlib.pyplot as plt
import scipy.sparse as sp

### golbal ----------------------------
momac_color = {
  "IL1B+NLRP3+_Macro":"#A77980",
  "MMP3+CXCL8+_Macro":"#BF94B2",
  "CCL13+_Complement-associated_Macro" : "#4d75be",
  "SPP1+CCL2+_Macro":"#D9D19B",
  "CD1C+_DC-like_Macro":"#83978C"
}





### load data -------------------------
## (0.1)data prepare
# all moMac 10% downsampled
Mac_10 = sc.read_h5ad("/storage/data/KAI/Pan_Inflammation/Mac2025/analysis/cytotrace/data/Mac_10.h5ad")
# common cells filtered
dm_eigen = pd.read_csv("/storage/data/KAI/Pan_Inflammation/Mac2025/analysis/palantir/data/dm_eigen.txt", sep='\t')
# pd.set_option('display.max_columns', 10)
dm_eigen.index.intersection(Mac_10.obs.index) # 匹配细胞数27886，同R中结果
Mac_palafiltered = Mac_10[Mac_10.obs.index.isin(dm_eigen.index.intersection(Mac_10.obs.index))].copy()

# embedding DCs & dpt
common_cells = dm_eigen.index.intersection(Mac_10.obs.index)
Mac_palafiltered.obs['DPT'] = dm_eigen.loc[common_cells, 'pseudotime']
Mac_palafiltered.obsm['X_diffusion_map'] = dm_eigen.loc[common_cells, ['DC1','DC2']].values
adata = Mac_palafiltered.copy()

DM_EigenValues = pd.read_csv('/storage/data/KAI/Pan_Inflammation/Mac2025/analysis/palantir/data/DM_EigenValues.txt',sep='\t')
DM_EigenVectors = pd.read_csv('/storage/data/KAI/Pan_Inflammation/Mac2025/analysis/palantir/data/DM_EigenVectors.txt',sep='\t')

adata.obsm['DM_EigenVectors'] = DM_EigenVectors.reindex(adata.obs.index).values
adata.uns['DM_EigenValues'] = DM_EigenValues.values

## (1)Run diffusion maps
# dm_res = palantir.utils.run_diffusion_maps(adata, n_components=5)
ms_data = palantir.utils.determine_multiscale_space(adata)



imputed_X = palantir.utils.run_magic_imputation(adata)
# palantir.plot.plot_diffusion_components(adata)
# plt.show()



# 不指定分化终点，进行轨迹推断
# GSM7883944_ACCCACTGTCCTGCTT-1-9-4-1 作为分化起点: dpt = 0
pr_res = palantir.core.run_palantir(adata,
                                    early_cell='GSM7883944_ACCCACTGTCCTGCTT-1-9-4-1',
                                    num_waypoints=500,
                                    # terminal_states=terminal_states
                                    )

Pyroptosis_genes = ['BAK1','BAX','CASP1','CASP3','CASP4','CASP5','CHMP2A','CHMP2B','CHMP3','CHMP4A','CHMP4B','CHMP4C','CHMP6','CHMP7','CYCS','ELANE','GSDMD','GSDME','GZMB','HMGB1','IL18','IL1A','IL1B','IRF1','IRF2','TP53','TP63']

sc.pl.embedding(
    adata,
    basis="X_diffusion_map",
    # layer="MAGIC_imputed_data",
    color=Pyroptosis_genes,
    frameon=False,
    use_raw=False
)

print(os.getcwd()) # 查找工作路径
# os.chdir('/storage/data/KAI/Pan_Inflammation/Infla_py') # 更改工作路径
fig, ax = plt.subplots(figsize=[12, 12])
palantir.plot.plot_palantir_results(adata, s=3, embedding_basis='X_diffusion_map')
plt.savefig('Mac2025_palantir_pseudotime.pdf', format="pdf", bbox_inches="tight")
plt.close()
# plt.show()

masks = palantir.presults.select_branch_cells(adata, q=.01, eps=.01)
fig, ax = plt.subplots(figsize=[12, 8])
palantir.plot.plot_branch_selection(adata,embedding_basis='X_diffusion_map')
plt.savefig('Mac2025_palantir_trajectory.pdf', format="pdf", bbox_inches="tight")
plt.close()
# plt.show()

fig, ax = plt.subplots(figsize=[5, 4])
palantir.plot.plot_trajectory(adata, "GSM7103304_CTAAGACTCAGGCAAG-1-2-5-0",embedding_basis='X_diffusion_map')
plt.savefig('Mac2025_palantir_SPP1_trajectory.pdf', format="pdf", bbox_inches="tight")
plt.close()
# plt.show()

# log
sc.pp.normalize_total(adata, target_sum=1e4)  # 标准化到 10,000
sc.pp.log1p(adata)  # 对数变换
# save anndata
adata.write('/storage/data/KAI/Pan_Inflammation/Mac2025/analysis/palantir/data/Mac2025_palantir.h5ad')


palantir.plot.plot_terminal_state_probs(adata, )
plt.show()
palantir.plot.highlight_cells_on_umap(adata, ['GSM7103304_CTAAGACTCAGGCAAG-1-2-5-0','GSM7457969_AATAACGCGACG-1-1'], embedding_basis='X_diffusion_map')
plt.show()

palantir.plot.plot_trajectory(adata, "GSM7103304_CTAAGACTCAGGCAAG-1-2-5-0", embedding_basis='X_diffusion_map')
plt.show()

palantir.plot.plot_trajectory(
    adata, # your anndata
    "GSM7103304_CTAAGACTCAGGCAAG-1-2-5-0", # the branch to plot
    cell_color="palantir_entropy", # the ad.obs colum to color the cells by
    n_arrows=10, # the number of arrow heads along the path
    color="red", # the color of the path and arrow heads
    scanpy_kwargs=dict(cmap="viridis"), # arguments passed to scanpy.pl.embedding
    arrowprops=dict(arrowstyle="->,head_length=.5,head_width=.5", lw=3), # appearance of the arrow heads
    lw=3, # thickness of the path
    pseudotime_interval=(0, .9), # interval of the pseudotime to cover with the path
    embedding_basis='X_diffusion_map',
)
plt.show()

adata = sc.read_h5ad('/storage/data/KAI/Pan_Inflammation/Mac2025/analysis/palantir/data/Mac2025_palantir.h5ad')

DM_Similarity = pd.read_csv("/storage/data/KAI/Pan_Inflammation/Mac2025/analysis/palantir/data/DM_Similarity.txt", sep='\t')
adata.obsp['DM_Similarity'] = sp.csr_matrix(DM_Similarity)

imputed_X = palantir.utils.run_magic_imputation(adata)

# save anndata
adata.write('/storage/data/KAI/Pan_Inflammation/Mac2025/analysis/palantir/data/Mac2025_palantir.h5ad')


gene_trends = palantir.presults.compute_gene_trends(
    adata,
    expression_key="MAGIC_imputed_data",
)

genes = ["IL1B", "NLRP3", "VEGFA", "STAT3",
         "CASP1","CASP3","CASP4","CASP5",
         "CASP6","CASP8","CASP9","ELANE",
         "GSDMA","GSDMB","GSDMC","GSDMD",
         ]
palantir.plot.plot_gene_trends(adata, genes)
plt.show()

palantir.plot.plot_gene_trend_heatmaps(adata, genes)
plt.show()


