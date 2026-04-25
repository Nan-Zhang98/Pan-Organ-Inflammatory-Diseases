# Module
import matplotlib.pyplot as plt
import scanpy as sc
import pandas as pd
import numpy as np
import subprocess
import os
import re
import scvi
import matplotlib as mlt
import matplotlib.pyplot as plt
from scib_metrics.benchmark import Benchmarker
import bbknn


#import omicverse as ov

### Mouth Inflammation scRNA-seq

adata_list = []

# Chronic_apical_periodontitis

preprecess_path = "/storage/data/KAI/Pan_Inflammation/Mouth/Chronic_apical_periodontitis/1.GSE181688/GSE181688_py/GSE181688_preprocess.h5ad"
adata = sc.read_h5ad(preprecess_path)
print(f'GSE181688 dim: {adata.shape}')
adata_list.append(adata)


preprecess_path = "/storage/data/KAI/Pan_Inflammation/Mouth/Chronic_apical_periodontitis/2.GSE197680/GSE197680_py/GSE197680_preprocess.h5ad"
adata = sc.read_h5ad(preprecess_path)
print(f'GSE197680 dim: {adata.shape}')
adata_list.append(adata)


# Oral_lichen_planus
# preprecess_path = "/storage/data/KAI/Pan_Inflammation/Mouth/Oral_lichen_planus/2.GSE211630/GSE211630_py/GSE211630_preprocess.h5ad"
# adata = sc.read_h5ad(preprecess_path)
# print(f'GSE211630 dim: {adata.shape}')
# adata_list.append(adata)


# Periodontitis
preprecess_path = "/storage/data/KAI/Pan_Inflammation/Mouth/Periodontitis/1.GSE164241/GSE164241_py/GSE164241_preprocess.h5ad"
adata = sc.read_h5ad(preprecess_path)
print(f'GSE164241 dim: {adata.shape}')
adata_list.append(adata)


preprecess_path = "/storage/data/KAI/Pan_Inflammation/Mouth/Periodontitis/2.GSE152042/GSE152042_py/GSE152042_preprocess.h5ad"
adata = sc.read_h5ad(preprecess_path)
print(f'GSE152042 dim: {adata.shape}')
adata_list.append(adata)


preprecess_path = "/storage/data/KAI/Pan_Inflammation/Mouth/Periodontitis/5.GSE171213/GSE171213_py/GSE171213_preprocess.h5ad"
adata = sc.read_h5ad(preprecess_path)
print(f'GSE171213 dim: {adata.shape}')   #该数据集平台为BD，整合时缺少很多基因的表达，如EPCAM
adata_list.append(adata)


# 检查交集基因
len(adata_list)

var_list = []
for i in range(0,len(adata_list)):
    var_gse = list(adata_list[i].var_names)
    var_list.append(var_gse)


inte_adata = []
for i in range(1,len(var_list)+1):
    intesect = set.intersection(*[set(v) for v in var_list[:i]])
    inte_adata.append(intesect)

for i, inter in enumerate(inte_adata):
    print(f"交集（前 {i+1} 个集合）: {len(inter)}")



# Integration
combined_adata = adata_list[0].concatenate(*adata_list[1:], join='outer', batch_key=None, batch_categories=None) # 整合后 (291114, 8161)

temp = combined_adata

# bool转0-1
#temp.obs = temp.obs.astype({'predicted_doublet': 'int'})

# 转换object类型
temp.obs = temp.obs.astype({'Age': 'str'}) #数据中不能含有缺失值(NaN)

# bool转0-1
#temp.var = temp.var.astype({'mt': 'int'})
#temp.var = temp.var.astype({'ribo': 'int'})
#temp.var = temp.var.astype({'hb': 'int'})
#temp.var = temp.var.astype({'highly_variable': 'int'})

# 转换object类型
temp.var = temp.var.astype({'Gene_ids-0': 'str'})
temp.var = temp.var.astype({'gene_ids-1': 'str'})
temp.var = temp.var.astype({'gene_ids-3': 'str'})



# 删除var中多余的属性
temp.var.drop(columns=['gene_ids-1', 'gene_ids-3'], inplace=True)

rubb_var = ['n_cells_by_counts-0', 'mean_counts-0', 'log1p_mean_counts-0', 'pct_dropout_by_counts-0', 'total_counts-0', 'log1p_total_counts-0', 'n_cells-0', 'n_cells_by_counts-1', 'mean_counts-1', 'log1p_mean_counts-1', 'pct_dropout_by_counts-1', 'total_counts-1', 'log1p_total_counts-1', 'n_cells-1', 'n_cells_by_counts-2', 'mean_counts-2', 'log1p_mean_counts-2', 'pct_dropout_by_counts-2', 'total_counts-2', 'log1p_total_counts-2', 'n_cells-2', 'n_cells_by_counts-3', 'mean_counts-3', 'log1p_mean_counts-3', 'pct_dropout_by_counts-3', 'total_counts-3', 'log1p_total_counts-3', 'n_cells-3', 'n_cells_by_counts-4', 'mean_counts-4', 'log1p_mean_counts-4', 'pct_dropout_by_counts-4', 'total_counts-4', 'log1p_total_counts-4', 'n_cells-4']
combined_adata.var.drop(columns=rubb_var, inplace=True)

temp.var.rename(columns={'Gene_ids-0':'Gene_ids'}, inplace=True)


# 写出整合后文件(未进行scvi)
#temp.write("/storage/data/KAI/Pan_Inflammation/Mouth/Mouth_py/scvi/" + "Mouth_Inte.h5ad")


# scvi-----------------------
combined_adata = temp

#combined_adata = sc.read_h5ad("/storage/data/KAI/Pan_Inflammation/Mouth/Mouth_py/scvi/" + "Mouth_Inte.h5ad")

combined_adata.raw = combined_adata  # keep full dimension safe

#combined_adata.layers["counts"] = combined_adata.X
'''
sc.pp.highly_variable_genes(
    adata,
    flavor="seurat_v3",
    #n_top_genes=2000,
    min_mean=0.0125,
    max_mean=3,
    min_disp=0.5,
    layer="counts",
    #batch_key="GEO",
    subset=True
)
'''
sc.pp.highly_variable_genes(combined_adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

scvi.model.SCVI.setup_anndata(combined_adata, layer="counts", categorical_covariate_keys=["Sample_geo_accession", "Technology", "Sample_Site"])
vae = scvi.model.SCVI(combined_adata)
vae.train()
# 保存scvi模型
model_dir = os.path.join("/storage/data/KAI/Pan_Inflammation/Mouth/Mouth_py/scvi/model/mouth_allcells/", "mouth_allcells_scvi_model")
vae.save(model_dir, overwrite=True, prefix="mouth_allcells_scvi_")

vae = scvi.model.SCVI.load("/storage/data/KAI/Pan_Inflammation/Mouth/Mouth_py/scvi/model/mouth_allcells/mouth_allcells_scvi_model/", adata=combined_adata, prefix="mouth_allcells_scvi_")

combined_adata.obsm["X_scVI"] = vae.get_latent_representation()
combined_adata.layers["X_normalized_scVI"] = vae.get_normalized_expression(library_size=10e4)

# run PcA then generate UMAp plots
sc.tl.pca(combined_adata)
sc.pl.pca_variance_ratio(combined_adata,log=True)
# use scvI latent space for UMAp generation
sc.pp.neighbors(combined_adata,use_rep="X_scVI")
sc.tl.umap(combined_adata,min_dist=0.3)

#result_folder = "/storage/data/KAI/Pan_Inflammation/Mouth/Mouth_py/scvi/result/"

sc.pl.umap(combined_adata,color=["GEO","Technology","Group","Sample_Site"], ncols=2, frameon=False, save = "Mouth_scvi_batch.png")

#写出
combined_adata.write("/storage/data/KAI/Pan_Inflammation/Mouth/Mouth_py/scvi/" + "Mouth_scvi.h5ad")

#combined_adata = sc.read_h5ad("/storage/data/KAI/Pan_Inflammation/Mouth/Mouth_py/scvi/" + "Mouth_scvi.h5ad")


# BBKNN ------------------------------
combined_adata
sc.external.pp.bbknn(combined_adata, batch_key="Sample_geo_accession")

sc.tl.umap(combined_adata,min_dist=0.3)

sc.pl.umap(combined_adata,color=["GEO","Technology","Group","Sample_Site"],
           ncols=2,
           frameon=False,
           #save = "Mouth_scvi_batch.png"
           show=True
           )




combined_adata.raw.var.index = combined_adata.raw.var.index.astype(str)  # 确保索引为字符串
combined_adata.raw.var = combined_adata.raw.var.astype(str)  # 确保所有列的类型为字符串
rubb_var = ['n_cells_by_counts-0', 'mean_counts-0', 'log1p_mean_counts-0', 'pct_dropout_by_counts-0', 'total_counts-0', 'log1p_total_counts-0', 'n_cells-0', 'n_cells_by_counts-1', 'mean_counts-1', 'log1p_mean_counts-1', 'pct_dropout_by_counts-1', 'total_counts-1', 'log1p_total_counts-1', 'n_cells-1', 'n_cells_by_counts-2', 'mean_counts-2', 'log1p_mean_counts-2', 'pct_dropout_by_counts-2', 'total_counts-2', 'log1p_total_counts-2', 'n_cells-2', 'n_cells_by_counts-3', 'mean_counts-3', 'log1p_mean_counts-3', 'pct_dropout_by_counts-3', 'total_counts-3', 'log1p_total_counts-3', 'n_cells-3', 'n_cells_by_counts-4', 'mean_counts-4', 'log1p_mean_counts-4', 'pct_dropout_by_counts-4', 'total_counts-4', 'log1p_total_counts-4', 'n_cells-4']
rubb_var = ['mt-0', 'ribo-0', 'hb-0', 'mt-1', 'ribo-1', 'hb-1',
       'gene_ids-1', 'mt-2', 'ribo-2', 'hb-2', 'gene_ids-2', 'mt-3', 'ribo-3',
       'hb-3', 'gene_ids-3', 'feature_types-3', 'genome-3', 'mt-4', 'ribo-4',
       'hb-4']
combined_adata.var.drop(columns=rubb_var, inplace=True)

combined_adata.obs = combined_adata.obs.astype({'Age': 'str'})
combined_adata.write("/storage/data/KAI/Pan_Inflammation/new_inflammation/Mouth/data/Mouth_all_gene.h5ad")






# 细胞注释-----------------------
sc.tl.leiden(combined_adata, resolution=0.8)
fig, ax = plt.subplots(figsize=(10, 8), dpi=300)
sc.pl.umap(combined_adata, color=["leiden"],
           legend_loc='on data',
           ax=ax,
           save="_mouth_leiden.png"
           )


Epi_marker = ["EPCAM","KRT6B", "KRT5", "SFN", "KRT14", "SPRR1B"]

Fib_marker = ["DCN","LUM","FAP","COL1A2", "COL3A1", "COL6A1", "PDPN"]

Myo_marker = ["ACTA1", "ACTA2", "ACTN2", "MYL2", "MPZ", "PLP1"]

Mural_marker = ["PDGFRB", "RGS5", "KCNJ8"]

Osteocyte_marker = ["OMD", "COL1A1", "TNFRSF11B", "SPP1",  "PHEX",  "ALPL"]

Endo_marker = ["PECAM1","VWF","PLVAP","CDH5"]

Mac_marker = ["AIF1","CSF1R","LYZ", "CD163", "FCGR2A", "CSF1R"]

Mono_marker = ["S100A9","S100A8","CD14","CD68"]

Mast_marker = ["TPSB2","CPA3","TPSAB1","KIT", "CMA1", "MS4A2", "TPSB2"]

Dendritic_marker = ["CLEC9A", "IRF8", "CD40", "CD80", "CD83", "CCR7"]

Neutrophil_marker =["G0S2", "CXCL8", "SOD2", "NAMPT", "S100A8"]

TNK_marker =  ["TRAC", "CD3D", "CD3E", "CD3G", "CD2", "NKG7","GNLY","NCAM1","KLRD1", "KLRF1", "TRDC", "XCL2", "XCL1", "CCL4", "KLRC1"]

B_marker = ["MS4A1","CD79A","CD79B", "CD19", "BLNK", "FCRL5", "CD24", "IGHG1", "CD27", "BANK1"]

Plasma_marker = ["IGHA1", "IGHG1", "IGHG2", "IGKC", "IGLC2", "JCHAIN", "CD19", "SDC1", "CD79A", "CD24"]

Neural_marker = ["NRXN1","S100B","NGFR"]


sc.pl.violin(combined_adata, Epi_marker, groupby="leiden", stripplot=False)
sc.pl.violin(combined_adata, Fib_marker, groupby="leiden", stripplot=False)
sc.pl.violin(combined_adata, Myo_marker, groupby="leiden", stripplot=False)
sc.pl.violin(combined_adata, Mural_marker, groupby="leiden", stripplot=False)
sc.pl.violin(combined_adata, Osteocyte_marker, groupby="leiden", stripplot=False)
sc.pl.violin(combined_adata, Endo_marker, groupby="leiden", stripplot=False)
sc.pl.violin(combined_adata, Mac_marker, groupby="leiden", stripplot=False)
sc.pl.violin(combined_adata, Mono_marker, groupby="leiden", stripplot=False)
sc.pl.violin(combined_adata, Mast_marker, groupby="leiden", stripplot=False)
sc.pl.violin(combined_adata, Dendritic_marker, groupby="leiden", stripplot=False)
sc.pl.violin(combined_adata, Neutrophil_marker, groupby="leiden", stripplot=False)
sc.pl.violin(combined_adata, TNK_marker, groupby="leiden", stripplot=False)
sc.pl.violin(combined_adata, B_marker, groupby="leiden", stripplot=False)
sc.pl.violin(combined_adata, Plasma_marker, groupby="leiden", stripplot=False)
sc.pl.violin(combined_adata, Neural_marker, groupby="leiden", stripplot=False)


sc.tl.rank_genes_groups(combined_adata, "leiden", groups=["16"], method="wilcoxon")
result = combined_adata.uns['rank_genes_groups']
cluster_16_genes = pd.DataFrame({
    'gene': result['names']['16'],
    'logfoldchanges': result['logfoldchanges']['16'],
    'pvals': result['pvals']['16'],
    'pvals_adj': result['pvals_adj']['16']
})
cluster_16_genes['gene'].head(20).to_list()
sc.pl.rank_genes_groups(combined_adata, groups=["16"], n_genes=10, sharey=False)

sc.tl.rank_genes_groups(combined_adata, "leiden", groups=["19"], method="wilcoxon")
result = combined_adata.uns['rank_genes_groups']
cluster_19_genes = pd.DataFrame({
    'gene': result['names']['19'],
    'logfoldchanges': result['logfoldchanges']['19'],
    'pvals': result['pvals']['19'],
    'pvals_adj': result['pvals_adj']['19']
})
cluster_19_genes['gene'].head(20).to_list()

sc.tl.rank_genes_groups(combined_adata, "leiden", groups=["26"], method="wilcoxon")
result = combined_adata.uns['rank_genes_groups']
cluster_26_genes = pd.DataFrame({
    'gene': result['names']['26'],
    'logfoldchanges': result['logfoldchanges']['26'],
    'pvals': result['pvals']['26'],
    'pvals_adj': result['pvals_adj']['26']
})
cluster_26_genes['gene'].head(30).to_list()

sc.pl.umap(combined_adata, color=["leiden","TRAC"])

cell_annotation_level1 = {
                            "0":"Endothelial cell",
                            "1":"T/NK cell",
                            "2":"Fibroblast",
                            "3":"T/NK cell",
                            "4":"Fibroblast",
                            "5":"Macrophage",
                            "6":"Plasma cell",
                            "7":"Fibroblast",
                            "8":"Plasma cell",
                            "9":"Endothelial cell",
                            "10":"Epithelial cell",
                            "11":"Plasma cell",
                            "12":"B cell",
                            "13":"Plasma cell",
                            "14":"T/NK cell",
                            "15":"Mast cell",
                            "16":"Neutrophil",
                            "17":"Fibroblast",
                            "18":"Epithelial cell",
                            "19":"Endothelial cell",
                            "20":"T/NK cell",
                            "21":"Fibroblast",
                            "22":"Dendritic cell",
                            "23":"Dendritic cell",
                            "24":"Endothelial cell",
                            "25":"Endothelial cell",
                            "26":"T/NK cell",
                            "27":"Neural cell",
                            "28":"Endothelial cell"
}
leiden_annotation_level1 = pd.Series(cell_annotation_level1)
combined_adata.obs['cell_type_level1'] = combined_adata.obs['leiden'].map(leiden_annotation_level1)

cell_color = {
    "Epithelial cell":"#e43030",
    "Fibroblast":"#ff8831",
    "Endothelial cell":"#704ba3",
    "Macrophage":"#67a8cd",
    "Mast cell":"#ff9d9f",
    "Dendritic cell":"#cf9f88",
    "Neutrophil":"#ffc17f",
    "T/NK cell":"#50aa4b",
    "B cell":"#b3e19b",
    "Plasma cell":"#ab3181",
    "Neural cell":"#00afba"
}
#combined_adata.write("/storage/data/KAI/Pan_Inflammation/Mouth/Mouth_py/scvi/" + "Mouth_scvi_cell_level1.h5ad")

combined_adata = sc.read_h5ad("/storage/data/KAI/Pan_Inflammation/Mouth/Mouth_py/scvi/" + "Mouth_scvi_cell_level1.h5ad")


# ov.pl.embedding(combined_adata,
#     basis="X_umap",
#     color=['cell_type_level1'],title='Mouth_scvi_cell_annotation',
#     show=True,
#     frameon=False,
#     size=10,
# )

sc.pl.umap(combined_adata, color="cell_type_level1",
           legend_loc="on data",
           palette=cell_color,
           add_outline=False,
           legend_fontsize=8,
           legend_fontoutline=0.1,
           legend_fontweight="normal",
           title="Mouth_scvi_cell_annotation",
           frameon=False,
           save = "_Mouth_scvi_cell_level1_clear.pdf"
           )
sc.pl.umap(combined_adata, color="cell_type_level1",
           #legend_loc="on data",
           palette=cell_color,
           add_outline=False,
           legend_fontsize=8,
           legend_fontoutline=0.1,
           legend_fontweight="normal",
           title="Mouth_scvi_cell_annotation",
           frameon=False,
           save = "_Mouth_scvi_cell_level1_legend.pdf"
           )

# scanvi-----------------------

lvae = scvi.model.SCANVI.from_scvi_model(
    vae,
    adata=combined_adata,
    labels_key="cell_type_level1",
    unlabeled_category="Unknown",
)
lvae.train()
model_dir2 = os.path.join("/storage/data/KAI/Pan_Inflammation/Mouth/Mouth_py/scanvi/model/mouth_allcells/", "mouth_allcells_scanvi_model")
lvae.save(model_dir2, overwrite=True, prefix="mouth_allcells_scanvi_")

lvae = scvi.model.SCANVI.load("/storage/data/KAI/Pan_Inflammation/Mouth/Mouth_py/scanvi/model/mouth_allcells/mouth_allcells_scanvi_model/", adata=combined_adata, prefix="mouth_allcells_scanvi_")


combined_adata.obsm["X_scANVI"] = lvae.get_latent_representation()

# run PcA then generate UMAp plots
sc.tl.pca(combined_adata)
sc.pl.pca_variance_ratio(combined_adata,log=True)
# use scANVI latent space for UMAp generation
sc.pp.neighbors(combined_adata,use_rep="X_scANVI")
sc.tl.umap(combined_adata,min_dist=0.3)

# ov.pl.embedding(combined_adata,
#     basis="X_umap",
#     color=['cell_type_level1'],title='Mouth_scvi_cell_annotation',
#     show=True,
#     frameon=False,
#     size=10,
# )
fig, ax=plt.subplots(1,1,figsize=(5,4))
sc.pl.umap(combined_adata, color="cell_type_level1",
           legend_loc="right margin",
           palette=cell_color,
           add_outline=False,
           legend_fontsize=8,
           legend_fontoutline=0.1,
           legend_fontweight="normal",
           title="Mouth_scanvi_cell_annotation",
           frameon=False,
           show=False,
           ax=ax
           #save = "_Mouth_scanvi_cell_level1_legend.pdf"
           )
# fig = plt.gcf()
for collection in ax.collections:
    collection.set_rasterized(True)
# for ax in fig.get_axes():

plt.savefig("/storage/data/KAI/Pan_Inflammation/Infla_py/figures/umap_Mouth_scanvi_cell_level1_legend.pdf", format="pdf", dpi=600)
plt.close()

fig, ax=plt.subplots(1,1,figsize=(5,4))
sc.pl.umap(combined_adata, color="cell_type_level1",
           legend_loc="on data",
           palette=cell_color,
           add_outline=False,
           legend_fontsize=8,
           legend_fontoutline=0.1,
           legend_fontweight="normal",
           title="Mouth_scanvi_cell_annotation",
           frameon=False,
           show=False,
           ax=ax
           #save = "_Mouth_scanvi_cell_level1_clear.pdf"
           )
# fig = plt.gcf()
for collection in ax.collections:
    collection.set_rasterized(True)
# for ax in fig.get_axes():

plt.savefig("/storage/data/KAI/Pan_Inflammation/Infla_py/figures/umap_Mouth_scanvi_cell_level1_clear.pdf", format="pdf", dpi=600)
plt.close()





# Evaluate
bm = Benchmarker(
    combined_adata,
    batch_key="GEO",
    label_key="cell_type_level1",
    embedding_obsm_keys=["X_pca", "X_scVI", "X_scANVI"],
    n_jobs=-1,
)
bm.benchmark()
fig, ax=plt.subplots(1,1,figsize=(4,8))
bm.plot_results_table(min_max_scale=False, show=False)
plt.savefig("/storage/data/KAI/Pan_Inflammation/Infla_py/figures/mouth_scvi_scanvi_eva.pdf")

#combined_adata.write("/storage/data/KAI/Pan_Inflammation/Mouth/Mouth_py/scanvi/Mouth_scanvi.h5ad")


#combined_adata = sc.read_h5ad("/storage/data/KAI/Pan_Inflammation/Mouth/Mouth_py/scanvi/Mouth_scanvi.h5ad")


# marker violin

violin_markers = ["KRT5", "SFN", "KRT14",
                  "DCN","LUM","COL1A2",
                  "PECAM1","VWF","PLVAP",
                  "AIF1","LYZ","C1QA",
                  "IRF8", "CD83", "CLEC9A",
                  "G0S2", "CXCL8", "SOD2",
                  "CPA3","TPSAB1","TPSB2",
                  "TRAC", "CD3D", "CD2",
                  "MS4A1","BANK1", "CD79A",
                  "IGHA1", "IGHG1", "IGHG2",
                  "S100B","PLP1"
                  ]
sc.pl.StackedViolin(combined_adata,var_names=violin_markers,groupby="cell_type_level1").style(row_palette=None, linewidth=0,).show()

# 可视化效果不佳，将h5ad转换成rds在R中绘制小提琴图展示marker表达情况。

sc.pl.tracksplot(combined_adata,violin_markers,"cell_type_level1") #tracksplot展示marker表达效果不佳


sc.pl.umap(combined_adata,color=["GEO"], ncols=2, frameon=False, save = "Mouth_scanvi_batch.png")




# TO R
combined_adata = sc.read_h5ad("/storage/data/KAI/Pan_Inflammation/Mouth/Mouth_py/scanvi/Mouth_scanvi.h5ad")
combined_adata.X = combined_adata.layers['counts']
combined_adata.write("/storage/data/KAI/Pan_Inflammation/Mouth/Mouth_py/scanvi/Mouth_scanvi_toR.h5ad")





### Skin Inflammation scRNA-seq
adata_list = []

# Atopic_dermatitis

preprecess_path = "/storage/data/KAI/Pan_Inflammation/Skin/Atopic dermatitis/4.GSE222840/GSE222840_py/GSE222840_preprocess.h5ad"
adata = sc.read_h5ad(preprecess_path)
adata.obs['GEO']
print(f'GSE222840 dim: {adata.shape}')
adata_list.append(adata)

preprecess_path = "/storage/data/KAI/Pan_Inflammation/Skin/Atopic dermatitis/7.GSE230575/GSE230575_py/GSE230575_preprocess.h5ad"
adata = sc.read_h5ad(preprecess_path)
adata.obs['GEO']
print(f'GSE230575 dim: {adata.shape}')
adata_list.append(adata)

# Hidradenitis suppurativa

preprecess_path = "/storage/data/KAI/Pan_Inflammation/Skin/Hidradenitis suppurativa/1.GSE154775/GSE154775_py/GSE154775_preprocess.h5ad"
adata = sc.read_h5ad(preprecess_path)
adata.obs['GEO']
print(f'GSE154775 dim: {adata.shape}')
adata_list.append(adata)

preprecess_path = "/storage/data/KAI/Pan_Inflammation/Skin/Hidradenitis suppurativa/2.GSE175990/GSE175990_py/GSE175990_preprocess.h5ad"
adata = sc.read_h5ad(preprecess_path)
adata.obs['GEO']
print(f'GSE175990 dim: {adata.shape}')
adata_list.append(adata)

preprecess_path = "/storage/data/KAI/Pan_Inflammation/Skin/Hidradenitis suppurativa/3.GSE220116/GSE220116_py/GSE220116_preprocess.h5ad"
adata = sc.read_h5ad(preprecess_path)
adata.obs['GEO']
print(f'GSE220116 dim: {adata.shape}')
adata_list.append(adata)


# Prurigo nodularis
preprecess_path = "/storage/data/KAI/Pan_Inflammation/Skin/Prurigo nodularis/1.GSE233280/GSE233280_py/GSE233280_preprocess.h5ad"
adata = sc.read_h5ad(preprecess_path)
adata.obs['GEO']
print(f'GSE233280 dim: {adata.shape}')
adata_list.append(adata)

preprecess_path = "/storage/data/KAI/Pan_Inflammation/Skin/Prurigo nodularis/3.GSE273559/GSE273559_py/GSE273559_preprocess.h5ad"
adata = sc.read_h5ad(preprecess_path)
adata.obs['GEO']
print(f'GSE273559 dim: {adata.shape}')
adata_list.append(adata)



# Sarcoidosis
preprecess_path = "/storage/data/KAI/Pan_Inflammation/Skin/Sarcoidosis/1.GSE234901/GSE234901_py/GSE234901_preprocess.h5ad"
adata = sc.read_h5ad(preprecess_path)
adata.obs['GEO']
print(f'GSE234901 dim: {adata.shape}')
adata_list.append(adata)

preprecess_path = "/storage/data/KAI/Pan_Inflammation/Skin/Sarcoidosis/2.GSE192456/GSE192456_py/GSE192456_preprocess.h5ad"
adata = sc.read_h5ad(preprecess_path)
adata.obs['GEO'] = "GSE192456"
print(f'GSE192456 dim: {adata.shape}')
adata_list.append(adata)

# Sample_Info = Read_Sample_Info("/storage/data/KAI/Pan_Inflammation/Skin/Sarcoidosis/")
# original_obs_names = adata.obs_names.copy()
# adata.obs = adata.obs.merge(Sample_Info, on='Sample_geo_accession', how='left')
# adata.obs_names = original_obs_names
# adata.write("/storage/data/KAI/Pan_Inflammation/Skin/Sarcoidosis/2.GSE192456/GSE192456_py/GSE192456_scanpy.h5ad")



# Psoriasis
preprecess_path = "/storage/data/KAI/Pan_Inflammation/Skin/Psoriasis/1.GSE173706/GSE173706_py/GSE173706_preprocess_sp.h5ad"
adata = sc.read_h5ad(preprecess_path)
adata.obs['GEO']
#adata.obs['GEO'] = 'GSE173706'
#adata.write("/storage/data/KAI/Pan_Inflammation/Skin/Psoriasis/1.GSE173706/GSE173706_py/GSE173706_preprocess_sp.h5ad")
print(f'GSE173706 dim: {adata.shape}')
adata_list.append(adata)


preprecess_path = "/storage/data/KAI/Pan_Inflammation/Skin/Psoriasis/2.GSE198805/GSE198805_py/GSE198805_preprocess.h5ad"
adata = sc.read_h5ad(preprecess_path)
adata.obs['GEO']
print(f'GSE198805 dim: {adata.shape}')
adata_list.append(adata)


preprecess_path = "/storage/data/KAI/Pan_Inflammation/Skin/Psoriasis/5.GSE221648/GSE221648_py/GSE221648_preprocess.h5ad"
adata = sc.read_h5ad(preprecess_path)
adata.obs['GEO']
print(f'GSE221648 dim: {adata.shape}')
adata_list.append(adata)


preprecess_path = "/storage/data/KAI/Pan_Inflammation/Skin/Psoriasis/7.GSE162183/GSE162183_py/GSE162183_preprocess.h5ad"
adata = sc.read_h5ad(preprecess_path)
adata.obs['GEO']
print(f'GSE162183 dim: {adata.shape}')
adata_list.append(adata)


preprecess_path = "/storage/data/KAI/Pan_Inflammation/Skin/Psoriasis/8.GSE151177/GSE151177_py/GSE151177_preprocess.h5ad"
adata = sc.read_h5ad(preprecess_path)
adata.obs['GEO']
print(f'GSE151177 dim: {adata.shape}')
adata_list.append(adata)


preprecess_path = "/storage/data/KAI/Pan_Inflammation/Skin/Psoriasis/9.GSE230842/GSE230842_py/GSE230842_preprocess.h5ad"
adata = sc.read_h5ad(preprecess_path)
adata.obs['GEO']
print(f'GSE230842 dim: {adata.shape}')
adata_list.append(adata)

# 检查交集基因
len(adata_list)

var_list = []
for i in range(0,len(adata_list)):
    var_gse = list(adata_list[i].var_names)
    var_list.append(var_gse)


inte_adata = []
for i in range(1,len(var_list)+1):
    intesect = set.intersection(*[set(v) for v in var_list[:i]])
    inte_adata.append(intesect)

for i, inter in enumerate(inte_adata):
    print(f"交集（前 {i+1} 个集合）: {len(inter)}")





# Integration
combined_adata = adata_list[0].concatenate(*adata_list[1:], join='outer', batch_key=None, batch_categories=None) # 整合后 (1084964, 13591)

temp = combined_adata

# bool转0-1
#temp.obs = temp.obs.astype({'predicted_doublet': 'int'})

# 转换object类型
temp.obs = temp.obs.astype({'Age': 'str'}) #数据中不能含有缺失值(NaN)

# bool转0-1
#temp.var = temp.var.astype({'mt': 'int'})
#temp.var = temp.var.astype({'ribo': 'int'})
#temp.var = temp.var.astype({'hb': 'int'})
#temp.var = temp.var.astype({'highly_variable': 'int'})

# 转换object类型
temp.var = temp.var.astype({'Gene_ids-0': 'str'})
temp.var = temp.var.astype({'gene_ids-1': 'str'})
temp.var = temp.var.astype({'gene_ids-3': 'str'})


# 删除var中多余的属性
temp.var.drop(columns=['gene_ids-1', 'gene_ids-3'], inplace=True)

rubb_var = ['n_cells_by_counts-0', 'mean_counts-0', 'log1p_mean_counts-0', 'pct_dropout_by_counts-0', 'total_counts-0', 'log1p_total_counts-0', 'n_cells-0', 'gene_ids-1', 'n_cells_by_counts-1', 'mean_counts-1', 'log1p_mean_counts-1', 'pct_dropout_by_counts-1', 'total_counts-1', 'log1p_total_counts-1', 'n_cells-1', 'n_cells_by_counts-10', 'mean_counts-10', 'log1p_mean_counts-10', 'pct_dropout_by_counts-10', 'total_counts-10', 'log1p_total_counts-10', 'n_cells-10', 'n_cells_by_counts-11', 'mean_counts-11', 'log1p_mean_counts-11', 'pct_dropout_by_counts-11', 'total_counts-11', 'log1p_total_counts-11', 'n_cells-11', 'Gene_ids-11', 'n_cells_by_counts-12', 'mean_counts-12', 'log1p_mean_counts-12', 'pct_dropout_by_counts-12', 'total_counts-12', 'log1p_total_counts-12', 'n_cells-12', 'gene_ids-13', 'n_cells_by_counts-13', 'mean_counts-13', 'log1p_mean_counts-13', 'pct_dropout_by_counts-13', 'total_counts-13', 'log1p_total_counts-13', 'n_cells-13', 'gene_ids-14', 'n_cells_by_counts-14', 'mean_counts-14', 'log1p_mean_counts-14', 'pct_dropout_by_counts-14', 'total_counts-14', 'log1p_total_counts-14', 'n_cells-14', 'gene_ids-2', 'n_cells_by_counts-2', 'mean_counts-2', 'log1p_mean_counts-2', 'pct_dropout_by_counts-2', 'total_counts-2', 'log1p_total_counts-2', 'n_cells-2', 'gene_ids-3', 'n_cells_by_counts-3', 'mean_counts-3', 'log1p_mean_counts-3', 'pct_dropout_by_counts-3', 'total_counts-3', 'log1p_total_counts-3', 'n_cells-3', 'gene_ids-4', 'n_cells_by_counts-4', 'mean_counts-4', 'log1p_mean_counts-4', 'pct_dropout_by_counts-4', 'total_counts-4', 'log1p_total_counts-4', 'n_cells-4', 'gene_ids-5', 'n_cells_by_counts-5', 'mean_counts-5', 'log1p_mean_counts-5', 'pct_dropout_by_counts-5', 'total_counts-5', 'log1p_total_counts-5', 'n_cells-5', 'gene_ids-6', 'n_cells_by_counts-6', 'mean_counts-6', 'log1p_mean_counts-6', 'pct_dropout_by_counts-6', 'total_counts-6', 'log1p_total_counts-6', 'n_cells-6', 'feature_types-6', 'genome-6', 'gene_ids-7', 'n_cells_by_counts-7', 'mean_counts-7', 'log1p_mean_counts-7', 'pct_dropout_by_counts-7', 'total_counts-7', 'log1p_total_counts-7', 'n_cells-7', 'n_cells_by_counts-8', 'mean_counts-8', 'log1p_mean_counts-8', 'pct_dropout_by_counts-8', 'total_counts-8', 'log1p_total_counts-8', 'n_cells-8', 'n_cells_by_counts-9', 'mean_counts-9', 'log1p_mean_counts-9', 'pct_dropout_by_counts-9', 'total_counts-9', 'log1p_total_counts-9', 'n_cells-9']
temp.var.drop(columns=rubb_var, inplace=True)

temp.var.rename(columns={'Gene_ids-0':'Gene_ids'}, inplace=True)


# 写出整合后文件(未进行scvi)
#temp.write("/storage/data/KAI/Pan_Inflammation/Skin/Skin_py/scvi/" + "Skin_Inte.h5ad")

# scvi-----------------------
combined_adata = temp

#combined_adata = sc.read_h5ad("/storage/data/KAI/Pan_Inflammation/Skin/Skin_py/scvi/" + "Skin_Inte.h5ad")

combined_adata.raw = combined_adata  # keep full dimension safe

combined_adata.layers["counts"] = combined_adata.X
'''
sc.pp.highly_variable_genes(
    adata,
    flavor="seurat_v3",
    #n_top_genes=2000,
    min_mean=0.0125,
    max_mean=3,
    min_disp=0.5,
    layer="counts",
    #batch_key="GEO",
    subset=True
)
'''
sc.pp.highly_variable_genes(combined_adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

scvi.model.SCVI.setup_anndata(combined_adata, layer="counts", categorical_covariate_keys=["Sample_geo_accession", "Sample_Type", "Technology", "Sample_Site", "Group"])
vae = scvi.model.SCVI(combined_adata)
vae.train()
# 保存scvi模型
model_dir = os.path.join("/storage/data/KAI/Pan_Inflammation/Skin/Skin_py/scvi/model/skin_allcells/", "skin_allcells_scvi_model")
vae.save(model_dir, overwrite=True, prefix="skin_allcells_scvi_")

vae = scvi.model.SCVI.load("/storage/data/KAI/Pan_Inflammation/Skin/Skin_py/scvi/model/skin_allcells/skin_allcells_scvi_model/", adata=combined_adata, prefix="skin_allcells_scvi_")

combined_adata.obsm["X_scVI"] = vae.get_latent_representation()
combined_adata.layers["X_normalized_scVI"] = vae.get_normalized_expression(library_size=10e4)

# run PcA then generate UMAp plots
sc.tl.pca(combined_adata)
sc.pl.pca_variance_ratio(combined_adata,log=True)
# use scvI latent space for UMAp generation
sc.pp.neighbors(combined_adata,use_rep="X_scVI")
sc.tl.umap(combined_adata,min_dist=0.3)

#result_folder = "/storage/data/KAI/Pan_Inflammation/Skin/Skin_py/scvi/result/"

sc.pl.umap(combined_adata,color=["GEO","Technology","Group","Sample_Site"], ncols=2, frameon=False, save = "Skin_scvi_batch.png")

#写出
combined_adata.write("/storage/data/KAI/Pan_Inflammation/Skin/Skin_py/scvi/" + "Skin_scvi.h5ad")

#combined_adata = sc.read_h5ad("/storage/data/KAI/Pan_Inflammation/Skin/Skin_py/scvi/" + "Skin_scvi.h5ad")


# full gene ---------------------------
combined_adata.raw.var.index = combined_adata.raw.var.index.astype(str)  # 确保索引为字符串
combined_adata.raw.var = combined_adata.raw.var.astype(str)  # 确保所有列的类型为字符串
rubb_var = ['mt-0', 'ribo-0', 'hb-0', 'n_cells_by_counts-0',
       'mean_counts-0', 'log1p_mean_counts-0', 'pct_dropout_by_counts-0',
       'total_counts-0', 'log1p_total_counts-0','mt-9', 'ribo-9', 'hb-9', 'n_cells_by_counts-9', 'mean_counts-9',
       'log1p_mean_counts-9', 'pct_dropout_by_counts-9', 'total_counts-9',
       'log1p_total_counts-9', 'n_cells-9']
combined_adata.var.drop(columns=rubb_var, inplace=True)
rubb_var = [
'total_counts-1', 'log1p_total_counts-1', 'n_cells-1', 'mt-10', 'ribo-10', 'hb-10', 'n_cells_by_counts-10', 'mean_counts-10', 'log1p_mean_counts-10', 'pct_dropout_by_counts-10', 'total_counts-10', 'log1p_total_counts-10', 'n_cells-10', 'mt-11', 'ribo-11', 'hb-11', 'n_cells_by_counts-11', 'mean_counts-11', 'log1p_mean_counts-11', 'pct_dropout_by_counts-11', 'total_counts-11', 'log1p_total_counts-11', 'n_cells-11', 'Gene_ids-11', 'mt-12', 'ribo-12', 'hb-12', 'n_cells_by_counts-12', 'mean_counts-12', 'log1p_mean_counts-12', 'pct_dropout_by_counts-12', 'total_counts-12', 'log1p_total_counts-12', 'n_cells-12', 'gene_ids-13', 'mt-13', 'ribo-13', 'hb-13', 'n_cells_by_counts-13', 'mean_counts-13', 'log1p_mean_counts-13', 'pct_dropout_by_counts-13', 'total_counts-13', 'log1p_total_counts-13', 'n_cells-13', 'gene_ids-14', 'mt-14', 'ribo-14', 'hb-14', 'n_cells_by_counts-14', 'mean_counts-14', 'log1p_mean_counts-14', 'pct_dropout_by_counts-14', 'total_counts-14', 'log1p_total_counts-14', 'n_cells-14', 'gene_ids-2', 'mt-2', 'ribo-2', 'hb-2', 'n_cells_by_counts-2', 'mean_counts-2', 'log1p_mean_counts-2', 'pct_dropout_by_counts-2', 'total_counts-2', 'log1p_total_counts-2', 'n_cells-2', 'gene_ids-3', 'mt-3', 'ribo-3', 'hb-3', 'n_cells_by_counts-3', 'mean_counts-3', 'log1p_mean_counts-3', 'pct_dropout_by_counts-3', 'total_counts-3', 'log1p_total_counts-3', 'n_cells-3', 'gene_ids-4', 'mt-4', 'ribo-4', 'hb-4', 'n_cells_by_counts-4', 'mean_counts-4', 'log1p_mean_counts-4', 'pct_dropout_by_counts-4', 'total_counts-4', 'log1p_total_counts-4', 'n_cells-4', 'gene_ids-5', 'mt-5', 'ribo-5', 'hb-5', 'n_cells_by_counts-5', 'mean_counts-5', 'log1p_mean_counts-5', 'pct_dropout_by_counts-5', 'total_counts-5', 'log1p_total_counts-5', 'n_cells-5', 'gene_ids-6', 'mt-6', 'ribo-6', 'hb-6', 'n_cells_by_counts-6', 'mean_counts-6', 'log1p_mean_counts-6', 'pct_dropout_by_counts-6', 'total_counts-6', 'log1p_total_counts-6', 'n_cells-6', 'feature_types-6', 'genome-6', 'gene_ids-7', 'mt-7', 'ribo-7', 'hb-7', 'n_cells_by_counts-7', 'mean_counts-7', 'log1p_mean_counts-7', 'pct_dropout_by_counts-7', 'total_counts-7', 'log1p_total_counts-7', 'n_cells-7'
]
combined_adata.var.drop(columns=rubb_var, inplace=True)

combined_adata.obs = combined_adata.obs.astype({'Age': 'str'})
combined_adata.write("/storage/data/KAI/Pan_Inflammation/new_inflammation/Skin/data/Skin_all_gene.h5ad")



# 细胞注释-----------------------
sc.tl.leiden(combined_adata, resolution=0.8)
fig, ax = plt.subplots(figsize=(10, 8), dpi=300)
sc.pl.umap(combined_adata, color=["leiden"],
           legend_loc='on data',
           ax=ax,
           save="_skin_leiden.png"
           )


Keratinocyte_marker = ["KRT1", "KRT14", "DMKN", "KRT5", "KRT10", "CDSN"]

Melanocyte_marker = ["DCT","TYRP1","PMEL","MLANA","QPCT","MITF"]

Eccrine_gland_marker = ["MUCL1"]

Fib_marker = ["DCN","CFD","FAP","COL1A2", "COL3A1", "COL6A1", "APOD"]

Langerhans_marker = ["CD1C", "CD4"]

Mural_marker = ["RGS5", "ACTA2", "TAGLN", "TPM2", "MYL9", "NOTCH3"]

Endo_marker = ["VWF","PLVAP","CDH5", "CLDN5"]

Mac_marker = ["AIF1","CSF1R","LYZ", "CD163", "FCGR2A", "CSF1R", "CD163"]

Mono_marker = ["S100A9","S100A8","CD14","CD68"]

Mast_marker = ["CPA3","TPSAB1","KIT", "CMA1", "MS4A2", "CTSG"]

Dendritic_marker = ["ITGAX", "CLEC9A", "IRF8", "CD40", "CD80", "CD83", "CCR7"]

Neutrophil_marker =["G0S2", "NAMPT", "S100A8"]

TNK_marker =  ["CD3D", "CD3E", "CD3G","IL32", "CD52", "CXCR4" ,"CD2", "NKG7","GNLY","NCAM1","KLRD1", "KLRF1", "XCL2", "XCL1", "CCL4", "KLRC1"]

B_marker = ["MS4A1","CD79A","CD79B", "CD19", "BLNK", "CD27", "BANK1"]

Plasma_marker = ["CD19", "SDC1", "CD79A"]

Erythrocyte_marker = ["HBA1","HBA2","HBB"]

Neural_marker = ["MPZ", "PLP1","S100B", "MBP", "PMP22", "CRYAB"]

sc.pl.violin(combined_adata, Keratinocyte_marker, groupby="leiden", stripplot=False, show=True, use_raw=False, rotation=90)
sc.pl.violin(combined_adata, Melanocyte_marker, groupby="leiden", stripplot=False,show=True, use_raw=False, rotation=90, multi_panel= True)
#sc.pl.violin(combined_adata, Eccrine_gland_marker, groupby="leiden", stripplot=False, show=True, use_raw=False, rotation=90)
sc.pl.violin(combined_adata, Fib_marker, groupby="leiden", stripplot=False, show=True, use_raw=False, rotation=90)
sc.pl.violin(combined_adata, Mural_marker, groupby="leiden", stripplot=False, show=True, use_raw=False, rotation=90)
sc.pl.violin(combined_adata, Endo_marker, groupby="leiden", stripplot=False, show=True, use_raw=False, rotation=90)
sc.pl.violin(combined_adata, Langerhans_marker, groupby="leiden", stripplot=False, show=True, use_raw=False, rotation=90)
sc.pl.violin(combined_adata, Mac_marker, groupby="leiden", stripplot=False, show=True, use_raw=False, rotation=90)
#sc.pl.violin(combined_adata, Mono_marker, groupby="leiden", stripplot=False, show=True, use_raw=False, rotation=90)
sc.pl.violin(combined_adata, Mast_marker, groupby="leiden", stripplot=False, show=True, use_raw=False, rotation=90)
sc.pl.violin(combined_adata, Dendritic_marker, groupby="leiden", stripplot=False, show=True, use_raw=False, rotation=90)
sc.pl.violin(combined_adata, Neutrophil_marker, groupby="leiden", stripplot=False, show=True, use_raw=False, rotation=90)
sc.pl.violin(combined_adata, TNK_marker, groupby="leiden", stripplot=False, show=True, use_raw=False, rotation=90)
sc.pl.violin(combined_adata, B_marker, groupby="leiden", stripplot=False, show=True, use_raw=False, rotation=90)
sc.pl.violin(combined_adata, Plasma_marker, groupby="leiden", stripplot=False, show=True, use_raw=False, rotation=90)
sc.pl.violin(combined_adata, Erythrocyte_marker, groupby="leiden", stripplot=False, show=True, use_raw=False, rotation=90)
sc.pl.violin(combined_adata, Neural_marker, groupby="leiden", stripplot=False, show=True, use_raw=False, rotation=90)


sc.tl.rank_genes_groups(combined_adata, "leiden", groups=["20"], method="wilcoxon")
result = combined_adata.uns['rank_genes_groups']
cluster_20_genes = pd.DataFrame({
    'gene': result['names']['20'],
    'logfoldchanges': result['logfoldchanges']['20'],
    'pvals': result['pvals']['20'],
    'pvals_adj': result['pvals_adj']['20']
})
cluster_20_genes['gene'].head(20).to_list()

sc.tl.rank_genes_groups(combined_adata, "leiden", groups=["27"], method="wilcoxon")
result = combined_adata.uns['rank_genes_groups']
cluster_27_genes = pd.DataFrame({
    'gene': result['names']['27'],
    'logfoldchanges': result['logfoldchanges']['27'],
    'pvals': result['pvals']['27'],
    'pvals_adj': result['pvals_adj']['27']
})
cluster_27_genes['gene'].head(20).to_list()

sc.tl.rank_genes_groups(combined_adata, "leiden", groups=["37"], method="wilcoxon")
result = combined_adata.uns['rank_genes_groups']
cluster_37_genes = pd.DataFrame({
    'gene': result['names']['37'],
    'logfoldchanges': result['logfoldchanges']['37'],
    'pvals': result['pvals']['37'],
    'pvals_adj': result['pvals_adj']['37']
})
cluster_37_genes['gene'].head(20).to_list()



cell_annotation_level1 = {
                        "0":"Endothelial cell",
                        "1":"Fibroblast",
                        "2":"Macrophage",
                        "3":"Keratinocyte",
                        "4":"Keratinocyte",
                        "5":"T/NK cell",
                        "6":"Fibroblast",
                        "7":"Keratinocyte",
                        "8":"Keratinocyte",
                        "9":"Keratinocyte",
                        "10":"Mural cell",
                        "11":"T/NK cell",
                        "12":"Keratinocyte",
                        "13":"Keratinocyte",
                        "14":"Keratinocyte",
                        "15":"Keratinocyte",
                        "16":"Fibroblast",
                        "17":"Mural cell",
                        "18":"T/NK cell",
                        "19":"T/NK cell",
                        "20":"Endothelial cell",
                        "21":"T/NK cell",
                        "22":"Melanocyte",
                        "23":"Keratinocyte",
                        "24":"Dendritic cell",
                        "25":"Endothelial cell",
                        "26":"Endothelial cell",
                        "27":"Keratinocyte",
                        "28":"Mast cell",
                        "29":"B cell",
                        "30":"B cell",
                        "31":"Keratinocyte",
                        "32":"Neural cell",
                        "33":"Neutrophil",
                        "34":"Macrophage",
                        "35":"Keratinocyte",
                        "36":"Keratinocyte",
                        "37":"Erythrocyte",
                        "38":"Fibroblast",
                        "39":"T/NK cell",
                        "40":"Keratinocyte",
                        "41":"Keratinocyte",
                        "42":"Keratinocyte",
                        "43":"Keratinocyte",
                        "44":"Keratinocyte",
                        "45":"Keratinocyte",
                        "46":"Macrophage",
                        "47":"Keratinocyte",
                        "48":"Keratinocyte",
                        "49":"Keratinocyte",

}
leiden_annotation_level1 = pd.Series(cell_annotation_level1)
combined_adata.obs['cell_type_level1'] = combined_adata.obs['leiden'].map(leiden_annotation_level1)

cell_color = {
    "Keratinocyte":"#4d79a6",
    "Melanocyte":"#37d5ab",
    "Mural cell":"#20452e",
    "Erythrocyte":"#a1acbd",
    "Epithelial cell":"#e43030",
    "Fibroblast":"#ff8831",
    "Endothelial cell":"#704ba3",
    "Macrophage":"#67a8cd",
    "Mast cell":"#ff9d9f",
    "Dendritic cell":"#cf9f88",
    "Neutrophil":"#ffc17f",
    "T/NK cell":"#50aa4b",
    "B cell":"#b3e19b",
    "Plasma cell":"#ab3181",
    "Neural cell":"#00afba"
}
#combined_adata.write("/storage/data/KAI/Pan_Inflammation/Skin/Skin_py/scvi/" + "Skin_scvi_cell_level1.h5ad")

combined_adata = sc.read_h5ad("/storage/data/KAI/Pan_Inflammation/Skin/Skin_py/scvi/" + "Skin_scvi_cell_level1.h5ad")


sc.pl.umap(combined_adata, color="cell_type_level1",
           legend_loc="on data",
           palette=cell_color,
           add_outline=False,
           legend_fontsize=8,
           legend_fontoutline=0.1,
           legend_fontweight="normal",
           title="Skin_scvi_cell_annotation",
           frameon=False,
           save = "_Skin_scvi_cell_level1_clear.pdf"
           )
sc.pl.umap(combined_adata, color="cell_type_level1",
           #legend_loc="on data",
           palette=cell_color,
           add_outline=False,
           legend_fontsize=8,
           legend_fontoutline=0.1,
           legend_fontweight="normal",
           title="Skin_scvi_cell_annotation",
           frameon=False,
           save = "_Skin_scvi_cell_level1_legend.pdf"
           )


# scanvi-----------------------

lvae = scvi.model.SCANVI.from_scvi_model(
    vae,
    adata=combined_adata,
    labels_key="cell_type_level1",
    unlabeled_category="Unknown",
)
lvae.train()
model_dir2 = os.path.join("/storage/data/KAI/Pan_Inflammation/Skin/Skin_py/scanvi/model/skin_allcells/", "skin_allcells_scanvi_model")
lvae.save(model_dir2, overwrite=True, prefix="skin_allcells_scanvi_")

lvae = scvi.model.SCANVI.load("/storage/data/KAI/Pan_Inflammation/Skin/Skin_py/scanvi/model/skin_allcells/skin_allcells_scanvi_model/", adata=combined_adata, prefix="skin_allcells_scanvi_")


combined_adata.obsm["X_scANVI"] = lvae.get_latent_representation()


# use scANVI latent space for UMAp generation
sc.pp.neighbors(combined_adata,use_rep="X_scANVI")
sc.tl.umap(combined_adata,min_dist=0.3)


fig, ax=plt.subplots(1,1,figsize=(5,4))
sc.pl.umap(combined_adata, color="cell_type_level1",
           legend_loc="right margin",
           palette=cell_color,
           add_outline=False,
           legend_fontsize=8,
           legend_fontoutline=0.1,
           legend_fontweight="normal",
           title="Skin_scanvi_cell_annotation",
           frameon=False,
           show=False,
           ax=ax
           #save = "_Skin_scanvi_cell_level1_legend.pdf"
           )
# fig = plt.gcf()
for collection in ax.collections:
    collection.set_rasterized(True)
# for ax in fig.get_axes():

plt.savefig("/storage/data/KAI/Pan_Inflammation/Infla_py/figures/umap_Skin_scanvi_cell_level1_legend.pdf", format="pdf", dpi=600)
plt.close()

fig, ax=plt.subplots(1,1,figsize=(5,4))
sc.pl.umap(combined_adata, color="cell_type_level1",
           legend_loc="on data",
           palette=cell_color,
           add_outline=False,
           legend_fontsize=8,
           legend_fontoutline=0.1,
           legend_fontweight="normal",
           title="Skin_scanvi_cell_annotation",
           frameon=False,
           show=False,
           ax=ax
           #save = "_Skin_scanvi_cell_level1_clear.pdf"
           )
# fig = plt.gcf()
for collection in ax.collections:
    collection.set_rasterized(True)
# for ax in fig.get_axes():

plt.savefig("/storage/data/KAI/Pan_Inflammation/Infla_py/figures/umap_Skin_scanvi_cell_level1_clear.pdf", format="pdf", dpi=600)
plt.close()



# Evaluate
bm = Benchmarker(
    combined_adata,
    batch_key="GEO",
    label_key="cell_type_level1",
    embedding_obsm_keys=["X_pca", "X_scVI", "X_scANVI"],
    n_jobs=-1,
)
bm.benchmark()
fig, ax=plt.subplots(1,1,figsize=(4,8))
bm.plot_results_table(min_max_scale=False, show=False)
plt.savefig("/storage/data/KAI/Pan_Inflammation/Infla_py/figures/Skin_scvi_scanvi_eva.pdf")

#combined_adata.write("/storage/data/KAI/Pan_Inflammation/Skin/Skin_py/scanvi/Skin_scanvi.h5ad")


#combined_adata = sc.read_h5ad("/storage/data/KAI/Pan_Inflammation/Skin/Skin_py/scanvi/Skin_scanvi.h5ad")

violin_markers = [
  "KRT1", "KRT14", "DMKN",
  "DCT","PMEL","MLANA",
  "DCN","COL1A2", "COL6A1",
  "ACTA2", "TAGLN", "MYL9",
  "VWF","PLVAP", "CLDN5",
  "AIF1", "CSF1R", "CD68",
  "CD83", "CCR7",
  "G0S2", "NAMPT", "S100A8",
  "CPA3","TPSAB1", "CTSG",
  "CD3D", "CD3E", "CD3G",
  "MS4A1","CD79A","CD79B",
  "MPZ", "PLP1","S100B",
  "HBA1","HBA2","HBB"
                  ]

sc.pl.umap(combined_adata,color=["GEO"], ncols=2, frameon=False, save = "Skin_scanvi_batch.png")



# TO R
combined_adata = sc.read_h5ad("/storage/data/KAI/Pan_Inflammation/Skin/Skin_py/scanvi/Skin_scanvi.h5ad")
combined_adata.X = combined_adata.layers['counts']
combined_adata.write("/storage/data/KAI/Pan_Inflammation/Skin/Skin_py/scanvi/Skin_scanvi_toR.h5ad")





## Bladder Inflammation scRNA-seq

adata_list = []

# Interstitial_cystitis

preprecess_path = "/storage/data/KAI/Pan_Inflammation/Bladder/Cystitis_glandularis/1.GSE225190/GSE225190_py/GSE225190_preprocess.h5ad"
adata = sc.read_h5ad(preprecess_path)
print(f'GSE225190 dim: {adata.shape}')
adata_list.append(adata)

# Interstitial_cystitis

preprecess_path = "/storage/data/KAI/Pan_Inflammation/Bladder/Interstitial_cystitis/2.GSE175526/GSE175526_py/GSE175526_preprocess.h5ad"
adata = sc.read_h5ad(preprecess_path)
print(f'GSE175526 dim: {adata.shape}')
adata_list.append(adata)


# 检查交集基因
len(adata_list)

var_list = []
for i in range(0,len(adata_list)):
    var_gse = list(adata_list[i].var_names)
    var_list.append(var_gse)


inte_adata = []
for i in range(1,len(var_list)+1):
    intesect = set.intersection(*[set(v) for v in var_list[:i]])
    inte_adata.append(intesect)

for i, inter in enumerate(inte_adata):
    print(f"交集（前 {i+1} 个集合）: {len(inter)}")


# Integration
combined_adata = adata_list[0].concatenate(*adata_list[1:], join='outer', batch_key=None, batch_categories=None) # 整合后 (159514, 23755)

temp = combined_adata

# bool转0-1
#temp.obs = temp.obs.astype({'predicted_doublet': 'int'})

# 转换object类型
temp.obs = temp.obs.astype({'Age': 'str'}) #数据中不能含有缺失值(NaN)

# bool转0-1
#temp.var = temp.var.astype({'mt': 'int'})
#temp.var = temp.var.astype({'ribo': 'int'})
#temp.var = temp.var.astype({'hb': 'int'})
#temp.var = temp.var.astype({'highly_variable': 'int'})

# 转换object类型
temp.var = temp.var.astype({'Gene_ids-0': 'str'})



# 删除var中多余的属性

rubb_var = ['n_cells_by_counts-0', 'mean_counts-0', 'log1p_mean_counts-0', 'pct_dropout_by_counts-0', 'total_counts-0', 'log1p_total_counts-0', 'n_cells-0', 'gene_ids-1', 'n_cells_by_counts-1', 'mean_counts-1', 'log1p_mean_counts-1', 'pct_dropout_by_counts-1', 'total_counts-1', 'log1p_total_counts-1', 'n_cells-1']
temp.var.drop(columns=rubb_var, inplace=True)

temp.var.rename(columns={'Gene_ids-0':'Gene_ids'}, inplace=True)


# 写出整合后文件(未进行scvi)
#temp.write("/storage/data/KAI/Pan_Inflammation/Bladder/Bladder_py/scvi/" + "Bladder_Inte.h5ad")


# scvi-----------------------
combined_adata = temp

#combined_adata = sc.read_h5ad("/storage/data/KAI/Pan_Inflammation/Bladder/Bladder_py/scvi/" + "Bladder_Inte.h5ad")

combined_adata.raw = combined_adata  # keep full dimension safe

combined_adata.layers["counts"] = combined_adata.X
'''
sc.pp.highly_variable_genes(
    adata,
    flavor="seurat_v3",
    #n_top_genes=2000,
    min_mean=0.0125,
    max_mean=3,
    min_disp=0.5,
    layer="counts",
    #batch_key="GEO",
    subset=True
)
'''
sc.pp.highly_variable_genes(combined_adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

scvi.model.SCVI.setup_anndata(combined_adata, layer="counts", categorical_covariate_keys=["Sample_geo_accession"])
vae = scvi.model.SCVI(combined_adata)
vae.train()
# 保存scvi模型
model_dir = os.path.join("/storage/data/KAI/Pan_Inflammation/Bladder/Bladder_py/scvi/model/bladder_allcells/", "bladder_allcells_scvi_model")
vae.save(model_dir, overwrite=True, prefix="bladder_allcells_scvi_")

vae = scvi.model.SCVI.load("/storage/data/KAI/Pan_Inflammation/Bladder/Bladder_py/scvi/model/bladder_allcells/bladder_allcells_scvi_model/", adata=combined_adata, prefix="mouth_allcells_scvi_")

combined_adata.obsm["X_scVI"] = vae.get_latent_representation()
combined_adata.layers["X_normalized_scVI"] = vae.get_normalized_expression(library_size=10e4)

# run PcA then generate UMAp plots
sc.tl.pca(combined_adata)
sc.pl.pca_variance_ratio(combined_adata,log=True)
# use scvI latent space for UMAp generation
sc.pp.neighbors(combined_adata,use_rep="X_scVI")
sc.tl.umap(combined_adata,min_dist=0.3)


#result_folder = "/storage/data/KAI/Pan_Inflammation/Bladder/Bladder_py/scvi/result/"

sc.pl.umap(combined_adata,color=["GEO","Technology","Group","Sample_Site"], ncols=2, frameon=False, save = "Bladder_scvi_batch.png")

#写出
combined_adata.write("/storage/data/KAI/Pan_Inflammation/Bladder/Bladder_py/scvi/" + "Bladder_scvi.h5ad")

#combined_adata = sc.read_h5ad("/storage/data/KAI/Pan_Inflammation/Bladder/Bladder_py/scvi/" + "Bladder_scvi.h5ad")

# full gene ---------------------------
combined_adata.raw.var.index = combined_adata.raw.var.index.astype(str)  # 确保索引为字符串
combined_adata.raw.var = combined_adata.raw.var.astype(str)  # 确保所有列的类型为字符串

rubb_var = [
'mt-0', 'ribo-0', 'hb-0', 'n_cells_by_counts-0', 'mean_counts-0', 'log1p_mean_counts-0', 'pct_dropout_by_counts-0', 'total_counts-0', 'log1p_total_counts-0', 'n_cells-0', 'gene_ids-1', 'mt-1', 'ribo-1', 'hb-1', 'n_cells_by_counts-1', 'mean_counts-1', 'log1p_mean_counts-1', 'pct_dropout_by_counts-1', 'total_counts-1', 'log1p_total_counts-1', 'n_cells-1'
]
combined_adata.var.drop(columns=rubb_var, inplace=True)

combined_adata.obs = combined_adata.obs.astype({'Age': 'str'})
combined_adata.write("/storage/data/KAI/Pan_Inflammation/new_inflammation/Bladder/data/Bladder_all_gene.h5ad")





# 细胞注释-----------------------
sc.tl.leiden(combined_adata, resolution=0.8)
fig, ax = plt.subplots(figsize=(10, 8), dpi=300)
sc.pl.umap(combined_adata, color=["leiden"],
           legend_loc='on data',
           ax=ax,
           save="_bladder_leiden.png"
           )



Epi_marker = ["EPCAM","KRT8","KRT18", "KRT19", "KRT13", "KRT17", "UPK1A"]

Fib_marker = ["DCN","LUM","COL1A2", "COL3A1", "COL6A1", "PDPN"]

Mural_marker = ["PDGFRB", "RGS5", "KCNJ8"]

Endo_marker = ["PECAM1","VWF","PLVAP","CLDN5"]

Mac_marker = ["AIF1","CSF1R","LYZ", "CD163", "FCGR2A", "CSF1R", "C1QA", "C1QB"]

Mono_marker = ["S100A9","S100A8"]

Mast_marker = ["TPSB2","CPA3","TPSAB1"]

Dendritic_marker = ["CLEC9A", "IRF8", "CD40", "CD80", "CD83", "CCR7"]

Neutrophil_marker =["G0S2", "CXCL8", "SOD2", "NAMPT", "S100A8"]

TNK_marker =  ["TRAC", "CD3D", "CD3E", "CD3G", "CD2", "NKG7", "TRDC", "XCL2", "XCL1"]

B_marker = ["MS4A1","CD79A","CD79B", "CD19", "BLNK", "FCRL5", "CD24", "IGHG1", "CD27", "BANK1"]

Plasma_marker = ["IGHA1", "IGHG1","IGKC", "IGLC2", "JCHAIN", "IGHM"]

Erythrocyte_marker = ["HBA1","HBA2","HBB"]

Neural_marker = ["NRXN1","S100B","NGFR"]


sc.pl.violin(combined_adata, Epi_marker, groupby="leiden", stripplot=False)
sc.pl.violin(combined_adata, Fib_marker, groupby="leiden", stripplot=False)
sc.pl.violin(combined_adata, Mural_marker, groupby="leiden", stripplot=False)
sc.pl.violin(combined_adata, Endo_marker, groupby="leiden", stripplot=False)
sc.pl.violin(combined_adata, Mac_marker, groupby="leiden", stripplot=False)
sc.pl.violin(combined_adata, Mono_marker, groupby="leiden", stripplot=False)
sc.pl.violin(combined_adata, Mast_marker, groupby="leiden", stripplot=False)
#sc.pl.violin(combined_adata, Dendritic_marker, groupby="leiden", stripplot=False)
#sc.pl.violin(combined_adata, Neutrophil_marker, groupby="leiden", stripplot=False)
sc.pl.violin(combined_adata, TNK_marker, groupby="leiden", stripplot=False)
sc.pl.violin(combined_adata, B_marker, groupby="leiden", stripplot=False)
sc.pl.violin(combined_adata, Plasma_marker, groupby="leiden", stripplot=False)
sc.pl.violin(combined_adata, Erythrocyte_marker, groupby="leiden", stripplot=False)
sc.pl.violin(combined_adata, Neural_marker, groupby="leiden", stripplot=False)


sc.tl.rank_genes_groups(combined_adata, "leiden", groups=["0"], method="wilcoxon")
result = combined_adata.uns['rank_genes_groups']
cluster_0_genes = pd.DataFrame({
    'gene': result['names']['0'],
    'logfoldchanges': result['logfoldchanges']['0'],
    'pvals': result['pvals']['0'],
    'pvals_adj': result['pvals_adj']['0']
})
cluster_0_genes['gene'].head(30).to_list()
sc.pl.umap(combined_adata, color=["DCN"],
           legend_loc='on data'
           )

sc.tl.rank_genes_groups(combined_adata, "leiden", groups=["4"], method="wilcoxon")
result = combined_adata.uns['rank_genes_groups']
cluster_4_genes = pd.DataFrame({
    'gene': result['names']['4'],
    'logfoldchanges': result['logfoldchanges']['4'],
    'pvals': result['pvals']['4'],
    'pvals_adj': result['pvals_adj']['4']
})
cluster_4_genes['gene'].head(20).to_list()

sc.tl.rank_genes_groups(combined_adata, "leiden", groups=["7"], method="wilcoxon")
result = combined_adata.uns['rank_genes_groups']
cluster_7_genes = pd.DataFrame({
    'gene': result['names']['7'],
    'logfoldchanges': result['logfoldchanges']['7'],
    'pvals': result['pvals']['7'],
    'pvals_adj': result['pvals_adj']['7']
})
cluster_7_genes['gene'].head(20).to_list()


sc.tl.rank_genes_groups(combined_adata, "leiden", groups=["12"], method="wilcoxon")
result = combined_adata.uns['rank_genes_groups']
cluster_12_genes = pd.DataFrame({
    'gene': result['names']['12'],
    'logfoldchanges': result['logfoldchanges']['12'],
    'pvals': result['pvals']['12'],
    'pvals_adj': result['pvals_adj']['12']
})
cluster_12_genes['gene'].head(20).to_list()


sc.tl.rank_genes_groups(combined_adata, "leiden", groups=["15"], method="wilcoxon")
result = combined_adata.uns['rank_genes_groups']
cluster_15_genes = pd.DataFrame({
    'gene': result['names']['15'],
    'logfoldchanges': result['logfoldchanges']['15'],
    'pvals': result['pvals']['15'],
    'pvals_adj': result['pvals_adj']['15']
})
cluster_15_genes['gene'].head(20).to_list()


sc.tl.rank_genes_groups(combined_adata, "leiden", groups=["26"], method="wilcoxon")
result = combined_adata.uns['rank_genes_groups']
cluster_26_genes = pd.DataFrame({
    'gene': result['names']['26'],
    'logfoldchanges': result['logfoldchanges']['26'],
    'pvals': result['pvals']['26'],
    'pvals_adj': result['pvals_adj']['26']
})
cluster_26_genes['gene'].head(50).to_list()

sc.tl.rank_genes_groups(combined_adata, "leiden", groups=["31"], method="wilcoxon")
result = combined_adata.uns['rank_genes_groups']
cluster_31_genes = pd.DataFrame({
    'gene': result['names']['31'],
    'logfoldchanges': result['logfoldchanges']['31'],
    'pvals': result['pvals']['31'],
    'pvals_adj': result['pvals_adj']['31']
})
cluster_31_genes['gene'].head(20).to_list()

cell_annotation_level1 = {
"0"	:"Fibroblast",
"1" :"Epithelial cell",
"2"	:"Epithelial cell",
"3"	:"T/NK cell",
"4"	:"Macrophage",
"5"	:"Fibroblast",
"6"	:"B cell",
"7"	:"Plasma cell",
"8"	:"Epithelial cell",
"9"	:"T/NK cell",
"10":"T/NK cell",
"11":"Plasma cell",
"12":"Endothelial cell",
"13":"Epithelial cell",
"14":"Endothelial cell",
"15":"Plasma cell",
"16":"Mural cell",
"17":"Macrophage",
"18":"B cell",
"19":"T/NK cell",
"20":"Fibroblast",
"21":"Epithelial cell",
"22":"Fibroblast",
"23":"Monocyte",
"24":"Mast cell",
"25":"Epithelial cell",
"26":"Macrophage",
"27":"Neural cell",
"28":"T/NK cell",
"29":"Macrophage",
"30":"Epithelial cell",
"31":"Erythrocyte",
"32":"Fibroblast"

}
leiden_annotation_level1 = pd.Series(cell_annotation_level1)
combined_adata.obs['cell_type_level1'] = combined_adata.obs['leiden'].map(leiden_annotation_level1)

cell_color = {
    "Epithelial cell":"#e43030",
    "Fibroblast":"#ff8831",
    "Mural cell":"#20452e",
    "Endothelial cell":"#704ba3",
    "Macrophage":"#67a8cd",
    "Monocyte":"#8184e2",
    "Mast cell":"#ff9d9f",
    "Dendritic cell":"#cf9f88",
    "Neutrophil":"#ffc17f",
    "T/NK cell":"#50aa4b",
    "B cell":"#b3e19b",
    "Plasma cell":"#ab3181",
    "Erythrocyte":"#a1acbd",
    "Neural cell":"#00afba"
}
#combined_adata.write("/storage/data/KAI/Pan_Inflammation/Bladder/Bladder_py/scvi/" + "Bladder_scvi_cell_level1.h5ad")

combined_adata = sc.read_h5ad("/storage/data/KAI/Pan_Inflammation/Bladder/Bladder_py/scvi/" + "Bladder_scvi_cell_level1.h5ad")


sc.pl.umap(combined_adata, color="cell_type_level1",
           legend_loc="on data",
           palette=cell_color,
           add_outline=False,
           legend_fontsize=8,
           legend_fontoutline=0.1,
           legend_fontweight="normal",
           title="Bladder_scvi_cell_annotation",
           frameon=False,
           save = "_Bladder_scvi_cell_level1_clear.pdf"
           )
sc.pl.umap(combined_adata, color="cell_type_level1",
           #legend_loc="on data",
           palette=cell_color,
           add_outline=False,
           legend_fontsize=8,
           legend_fontoutline=0.1,
           legend_fontweight="normal",
           title="Bladder_scvi_cell_annotation",
           frameon=False,
           save = "_Bladder_scvi_cell_level1_legend.pdf"
           )

# scanvi-----------------------

lvae = scvi.model.SCANVI.from_scvi_model(
    vae,
    adata=combined_adata,
    labels_key="cell_type_level1",
    unlabeled_category="Unknown",
)
lvae.train()
model_dir2 = os.path.join("/storage/data/KAI/Pan_Inflammation/Bladder/Bladder_py/scanvi/model/bladder_allcells/", "bladder_allcells_scanvi_model")
lvae.save(model_dir2, overwrite=True, prefix="bladder_allcells_scanvi_")

lvae = scvi.model.SCANVI.load("/storage/data/KAI/Pan_Inflammation/Bladder/Bladder_py/scanvi/model/bladder_allcells/bladder_allcells_scanvi_model/", adata=combined_adata, prefix="bladder_allcells_scanvi_")


combined_adata.obsm["X_scANVI"] = lvae.get_latent_representation()

# run PcA then generate UMAp plots
sc.tl.pca(combined_adata)
sc.pl.pca_variance_ratio(combined_adata,log=True)
# use scANVI latent space for UMAp generation
sc.pp.neighbors(combined_adata,use_rep="X_scANVI")
sc.tl.umap(combined_adata,min_dist=0.3)


fig, ax=plt.subplots(1,1,figsize=(5,4))
sc.pl.umap(combined_adata, color="cell_type_level1",
           legend_loc="right margin",
           palette=cell_color,
           add_outline=False,
           legend_fontsize=8,
           legend_fontoutline=0.1,
           legend_fontweight="normal",
           title="Bladder_scanvi_cell_annotation",
           frameon=False,
           show=False,
           ax=ax
           #save = "_Bladder_scanvi_cell_level1_legend.pdf"
           )
# fig = plt.gcf()
for collection in ax.collections:
    collection.set_rasterized(True)
# for ax in fig.get_axes():

plt.savefig("/storage/data/KAI/Pan_Inflammation/Infla_py/figures/umap_Bladder_scanvi_cell_level1_legend.pdf", format="pdf", dpi=600)
plt.close()



fig, ax=plt.subplots(1,1,figsize=(5,4))
sc.pl.umap(combined_adata, color="cell_type_level1",
           legend_loc="on data",
           palette=cell_color,
           add_outline=False,
           legend_fontsize=8,
           legend_fontoutline=0.1,
           legend_fontweight="normal",
           title="Bladder_scanvi_cell_annotation",
           frameon=False,
           show=False,
           ax=ax
           #save = "_Bladder_scanvi_cell_level1_clear.pdf"
           )
# fig = plt.gcf()
for collection in ax.collections:
    collection.set_rasterized(True)
# for ax in fig.get_axes():

plt.savefig("/storage/data/KAI/Pan_Inflammation/Infla_py/figures/umap_Bladder_scanvi_cell_level1_clear.pdf", format="pdf", dpi=600)
plt.close()


#combined_adata.write("/storage/data/KAI/Pan_Inflammation/Bladder/Bladder_py/scanvi/Bladder_scanvi.h5ad")


#combined_adata = sc.read_h5ad("/storage/data/KAI/Pan_Inflammation/Bladder/Bladder_py/scanvi/Bladder_scanvi.h5ad")


sc.pl.umap(combined_adata,color=["GEO"], ncols=2, frameon=False, save = "Bladder_scanvi_batch.png")


# TO R
combined_adata.X = combined_adata.layers['counts']
combined_adata.write("/storage/data/KAI/Pan_Inflammation/Bladder/Bladder_py/scanvi/Bladder_scanvi_toR.h5ad")


### Uterus Inflammation scRNA-seq

adata_list = []

# Endometriosis

preprecess_path = "/storage/data/KAI/Pan_Inflammation/Uterus/Endometriosis/1.GSE213216/GSE213216_py/GSE213216_preprocess.h5ad"
adata = sc.read_h5ad(preprecess_path)
print(f'GSE213216 dim: {adata.shape}')
adata_list.append(adata)


preprecess_path = "/storage/data/KAI/Pan_Inflammation/Uterus/Endometriosis/2.GSE179640/GSE179640_py/GSE179640_preprocess.h5ad"
adata = sc.read_h5ad(preprecess_path)
print(f'GSE179640 dim: {adata.shape}')
adata_list.append(adata)


preprecess_path = "/storage/data/KAI/Pan_Inflammation/Uterus/Endometriosis/4.GSE214411/GSE214411_py/GSE214411_preprocess.h5ad"
adata = sc.read_h5ad(preprecess_path)
print(f'GSE214411 dim: {adata.shape}')
adata_list.append(adata)


# 检查交集基因
len(adata_list)

var_list = []
for i in range(0,len(adata_list)):
    var_gse = list(adata_list[i].var_names)
    var_list.append(var_gse)


inte_adata = []
for i in range(1,len(var_list)+1):
    intesect = set.intersection(*[set(v) for v in var_list[:i]])
    inte_adata.append(intesect)

for i, inter in enumerate(inte_adata):
    print(f"交集（前 {i+1} 个集合）: {len(inter)}")


# Integration
combined_adata = adata_list[0].concatenate(*adata_list[1:], join='outer', batch_key=None, batch_categories=None) # 整合后 (291114, 8161)

temp = combined_adata

# bool转0-1
#temp.obs = temp.obs.astype({'predicted_doublet': 'int'})

# 转换object类型
temp.obs = temp.obs.astype({'Age': 'str'}) #数据中不能含有缺失值(NaN)
temp.obs = temp.obs.astype({'Donor': 'str'})
# bool转0-1
#temp.var = temp.var.astype({'mt': 'int'})
#temp.var = temp.var.astype({'ribo': 'int'})
#temp.var = temp.var.astype({'hb': 'int'})
#temp.var = temp.var.astype({'highly_variable': 'int'})

# 转换object类型
temp.var = temp.var.astype({'Gene_ids-0': 'str'})
temp.var = temp.var.astype({'gene_ids-1': 'str'})
temp.var = temp.var.astype({'gene_ids-3': 'str'})


# 删除var中多余的属性
rubb_var = ['gene_ids-1-1',
       'gene_ids-10-1', 'gene_ids-11-1', 'gene_ids-12-1', 'gene_ids-13-1',
       'gene_ids-14-1', 'gene_ids-15-1', 'gene_ids-16-1', 'gene_ids-17-1',
       'gene_ids-18-1', 'gene_ids-19-1', 'gene_ids-2-1', 'gene_ids-20-1',
       'gene_ids-21-1', 'gene_ids-22-1', 'gene_ids-23-1', 'gene_ids-24-1',
       'gene_ids-25-1', 'gene_ids-26-1', 'gene_ids-27-1', 'gene_ids-28-1',
       'gene_ids-29-1', 'gene_ids-3-1', 'gene_ids-30-1', 'gene_ids-4-1',
       'gene_ids-5-1', 'gene_ids-6-1', 'gene_ids-7-1', 'gene_ids-8-1',
       'gene_ids-9-1', 'n_cells-2', 'n_cells_by_counts-2', 'mean_counts-2',
       'log1p_mean_counts-2', 'pct_dropout_by_counts-2', 'total_counts-2',
       'log1p_total_counts-2', 'gene_ids-2']
rubb_var = ['n_cells-0', 'n_cells_by_counts-0',
       'mean_counts-0', 'log1p_mean_counts-0', 'pct_dropout_by_counts-0',
       'total_counts-0', 'log1p_total_counts-0', 'n_cells-1',
       'n_cells_by_counts-1', 'mean_counts-1', 'log1p_mean_counts-1',
       'pct_dropout_by_counts-1', 'total_counts-1', 'log1p_total_counts-1',
       'feature_types-1', 'genome-1', 'gene_ids-0-1']
temp.var.drop(columns=rubb_var, inplace=True)

temp.var.rename(columns={'Gene_ids-0':'Gene_ids'}, inplace=True)




combined_adata = temp
combined_adata = combined_adata[combined_adata.obs['Group'] != 'Endometrioma'].copy()
combined_adata.obs['Sample_Site'] = combined_adata.obs['Sample_Site'].replace('Endometirum', 'Endometrium')


rubb_var = ['gene_ids-1-1', 'gene_ids-10-1', 'gene_ids-11-1', 'gene_ids-12-1', 'gene_ids-13-1', 'gene_ids-14-1', 'gene_ids-15-1', 'gene_ids-16-1', 'gene_ids-17-1', 'gene_ids-18-1', 'gene_ids-19-1', 'gene_ids-2-1', 'gene_ids-20-1', 'gene_ids-21-1', 'gene_ids-22-1', 'gene_ids-23-1', 'gene_ids-24-1', 'gene_ids-25-1', 'gene_ids-26-1', 'gene_ids-27-1', 'gene_ids-28-1', 'gene_ids-29-1', 'gene_ids-3-1', 'gene_ids-30-1', 'gene_ids-4-1', 'gene_ids-5-1', 'gene_ids-6-1', 'gene_ids-7-1', 'gene_ids-8-1', 'gene_ids-9-1', 'n_cells-2', 'n_cells_by_counts-2', 'mean_counts-2', 'log1p_mean_counts-2', 'pct_dropout_by_counts-2', 'total_counts-2', 'log1p_total_counts-2', 'gene_ids-2']
combined_adata.var.drop(columns=rubb_var, inplace=True)
# 写出整合后文件(未进行scvi)
#combined_adata.write("/storage/data/KAI/Pan_Inflammation/Uterus/Uterus_py/scvi/" + "Uterus_Inte.h5ad")



# scvi-----------------------
#combined_adata = sc.read_h5ad("/storage/data/KAI/Pan_Inflammation/Uterus/Uterus_py/scvi/" + "Uterus_Inte.h5ad")

combined_adata.raw = combined_adata  # keep full dimension safe

#combined_adata.layers["counts"] = combined_adata.X

sc.pp.highly_variable_genes(combined_adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

scvi.model.SCVI.setup_anndata(combined_adata, layer="counts", categorical_covariate_keys=["Sample_geo_accession", "Sample_Site","GEO" ])
vae = scvi.model.SCVI(combined_adata, n_latent=20, n_layers=2, dropout_rate=0.1)
vae.train(max_epochs=30, early_stopping=True, batch_size=256)
# 保存scvi模型
model_dir = os.path.join("/storage/data/KAI/Pan_Inflammation/Uterus/Uterus_py/scvi/model/uterus_allcells/", "uterus_allcells_scvi_model")
vae.save(model_dir, overwrite=True, prefix="uterus_allcells_scvi_")

vae = scvi.model.SCVI.load("/storage/data/KAI/Pan_Inflammation/Uterus/Uterus_py/scvi/model/uterus_allcells/uterus_allcells_scvi_model/", adata=combined_adata, prefix="uterus_allcells_scvi_")


combined_adata.obsm["X_scVI"] = vae.get_latent_representation()
combined_adata.layers["X_normalized_scVI"] = vae.get_normalized_expression(library_size=10e4)

combined_adata = combined_adata[combined_adata.obs['Group'] != 'Endometrioma'].copy()
combined_adata.obs['Sample_Site'] = combined_adata.obs['Sample_Site'].replace('Endometirum', 'Endometrium')
# run PcA then generate UMAp plots
sc.tl.pca(combined_adata)
sc.pl.pca_variance_ratio(combined_adata,log=True)
# use scvI latent space for UMAp generation
sc.pp.neighbors(combined_adata,use_rep="X_scVI",n_neighbors=30)
sc.tl.umap(combined_adata,min_dist=0.3)

#result_folder = "/storage/data/KAI/Pan_Inflammation/Uterus/Uterus_py/scvi/result/"

sc.pl.umap(combined_adata,color=["GEO","Technology","Group","Sample_Site"], ncols=2, frameon=False, save = "Uterus_scvi_batch.png")

#写出
combined_adata.write("/storage/data/KAI/Pan_Inflammation/Uterus/Uterus_py/scvi/" + "Uterus_scvi.h5ad")

#combined_adata = sc.read_h5ad("/storage/data/KAI/Pan_Inflammation/Uterus/Uterus_py/scvi/" + "Uterus_scvi.h5ad")

# full gene ---------------------------
combined_adata.raw.var.index = combined_adata.raw.var.index.astype(str)  # 确保索引为字符串
combined_adata.raw.var = combined_adata.raw.var.astype(str)  # 确保所有列的类型为字符串

rubb_var = [
'mt-0', 'ribo-0', 'hb-0', 'n_cells-0', 'n_cells_by_counts-0', 'mean_counts-0', 'log1p_mean_counts-0', 'pct_dropout_by_counts-0', 'total_counts-0', 'log1p_total_counts-0', 'mt-1', 'ribo-1', 'hb-1', 'n_cells-1', 'n_cells_by_counts-1', 'mean_counts-1', 'log1p_mean_counts-1', 'pct_dropout_by_counts-1', 'total_counts-1', 'log1p_total_counts-1', 'feature_types-1', 'genome-1', 'gene_ids-0-1', 'gene_ids-1-1', 'gene_ids-10-1', 'gene_ids-11-1', 'gene_ids-12-1', 'gene_ids-13-1', 'gene_ids-14-1', 'gene_ids-15-1', 'gene_ids-16-1', 'gene_ids-17-1', 'gene_ids-18-1', 'gene_ids-19-1', 'gene_ids-2-1', 'gene_ids-20-1', 'gene_ids-21-1', 'gene_ids-22-1', 'gene_ids-23-1', 'gene_ids-24-1', 'gene_ids-25-1', 'gene_ids-26-1', 'gene_ids-27-1', 'gene_ids-28-1', 'gene_ids-29-1', 'gene_ids-3-1', 'gene_ids-30-1', 'gene_ids-4-1', 'gene_ids-5-1', 'gene_ids-6-1', 'gene_ids-7-1', 'gene_ids-8-1', 'gene_ids-9-1', 'mt-2', 'ribo-2', 'hb-2', 'n_cells-2', 'n_cells_by_counts-2', 'mean_counts-2', 'log1p_mean_counts-2', 'pct_dropout_by_counts-2', 'total_counts-2', 'log1p_total_counts-2', 'gene_ids-2'
]
combined_adata.var.drop(columns=rubb_var, inplace=True)

combined_adata.obs = combined_adata.obs.astype({'Age': 'str'})
combined_adata.obs = combined_adata.obs.astype({'Donor': 'str'})
combined_adata.write("/storage/data/KAI/Pan_Inflammation/new_inflammation/Uterus/data/Uterus_all_gene.h5ad")






# 细胞注释-----------------------
sc.tl.leiden(combined_adata, resolution=0.8)
fig, ax = plt.subplots(figsize=(10, 8), dpi=300)
sc.pl.umap(combined_adata, color=["leiden"],
           legend_loc='on data',
           ax=ax,
           save="_uterus_leiden.png"
           )


Epi_marker = ["EPCAM","KRT8","KRT18", "KRT19"]

Fib_marker = ["DCN","LUM","COL1A2", "COL3A1", "COL6A1", "PDPN"]

Mural_marker = ["PDGFRB", "RGS5", "KCNJ8"]

Endo_marker = ["PECAM1","VWF","PLVAP","CLDN5"]

Mac_marker = ["AIF1","CSF1R","LYZ", "CD163", "FCGR2A", "CSF1R", "C1QA", "C1QB"]

Mono_marker = ["S100A9","S100A8"]

Mast_marker = ["TPSB2","CPA3","TPSAB1"]

Dendritic_marker = ["ITGAX","CLEC9A", "IRF8", "CD40", "CD80", "CD83", "CCR7"]

Neutrophil_marker =["G0S2", "CXCL8", "SOD2", "NAMPT", "S100A8"]

TNK_marker =  ["TRAC", "CD3D", "CD3E", "CD3G", "CD2", "NKG7", "TRDC", "XCL2", "XCL1"]

B_marker = ["MS4A1","CD79A","CD79B", "CD19", "BLNK", "FCRL5", "CD24", "IGHG1", "CD27", "BANK1"]

Plasma_marker = ["IGHA1", "IGHG1","IGKC", "IGLC2", "JCHAIN", "IGHM"]

Erythrocyte_marker = ["HBA1","HBA2","HBB"]

Neural_marker = ["NRXN1","S100B","NGFR"]


sc.pl.violin(combined_adata, Epi_marker, groupby="leiden", stripplot=False)
sc.pl.violin(combined_adata, Fib_marker, groupby="leiden", stripplot=False)
sc.pl.violin(combined_adata, Mural_marker, groupby="leiden", stripplot=False)
sc.pl.violin(combined_adata, Endo_marker, groupby="leiden", stripplot=False)
sc.pl.violin(combined_adata, Mac_marker, groupby="leiden", stripplot=False)
sc.pl.violin(combined_adata, Mono_marker, groupby="leiden", stripplot=False)
sc.pl.violin(combined_adata, Mast_marker, groupby="leiden", stripplot=False)
#sc.pl.violin(combined_adata, Dendritic_marker, groupby="leiden", stripplot=False)
#sc.pl.violin(combined_adata, Neutrophil_marker, groupby="leiden", stripplot=False)
sc.pl.violin(combined_adata, TNK_marker, groupby="leiden", stripplot=False)
sc.pl.violin(combined_adata, B_marker, groupby="leiden", stripplot=False)
sc.pl.violin(combined_adata, Plasma_marker, groupby="leiden", stripplot=False)
sc.pl.violin(combined_adata, Erythrocyte_marker, groupby="leiden", stripplot=False)
sc.pl.violin(combined_adata, Neural_marker, groupby="leiden", stripplot=False)


sc.tl.rank_genes_groups(combined_adata, "leiden", groups=["15"], method="wilcoxon")
result = combined_adata.uns['rank_genes_groups']
cluster_15_genes = pd.DataFrame({
    'gene': result['names']['15'],
    'logfoldchanges': result['logfoldchanges']['15'],
    'pvals': result['pvals']['15'],
    'pvals_adj': result['pvals_adj']['15']
})
cluster_15_genes['gene'].head(30).to_list()

cell_annotation_level1 = {
"0":	"Fibroblast",
"1":	"Fibroblast",
"2":	"Fibroblast",
"3":	"T/NK cell",
"4":	"Fibroblast",
"5":	"Fibroblast",
"6":	"Fibroblast",
"7":	"Epithelial cell",
"8":	"T/NK cell",
"9":	"Endothelial cell",
"10":	"Fibroblast",
"11":	"T/NK cell",
"12":	"Mural cell",
"13":	"Macrophage",
"14":	"Endothelial cell",
"15":   "Macrophage",
"16":	"Fibroblast",
"17":	"Epithelial cell",
"18":	"B cell",
"19":	"Fibroblast",
"20":	"Epithelial cell",
"21":	"Fibroblast",
"22":	"Monocyte",
"23":	"T/NK cell",
"24":	"Epithelial cell",
"25":	"Fibroblast",
"26":	"T/NK cell",
"27":	"Endothelial cell",
"28":	"T/NK cell",
"29":	"Mast cell",
"30":	"Fibroblast",
"31":	"T/NK cell",
"32":	"Plasma cell",
"33":	"Epithelial cell",
"34":	"Fibroblast",
"35":	"Epithelial cell",
"36":	"Neural cell"

}
leiden_annotation_level1 = pd.Series(cell_annotation_level1)
combined_adata.obs['cell_type_level1'] = combined_adata.obs['leiden'].map(leiden_annotation_level1)

cell_color = {
    "Epithelial cell":"#e43030",
    "Fibroblast":"#ff8831",
    "Mural cell":"#20452e",
    "Endothelial cell":"#704ba3",
    "Macrophage":"#67a8cd",
    "Monocyte":"#8184e2",
    "Mast cell":"#ff9d9f",
    "Dendritic cell":"#cf9f88",
    "Neutrophil":"#ffc17f",
    "T/NK cell":"#50aa4b",
    "B cell":"#b3e19b",
    "Plasma cell":"#ab3181",
    "Erythrocyte":"#a1acbd",
    "Neural cell":"#00afba"
}
#combined_adata.write("/storage/data/KAI/Pan_Inflammation/Uterus/Uterus_py/scvi/" + "Uterus_scvi_cell_level1.h5ad")

combined_adata = sc.read_h5ad("/storage/data/KAI/Pan_Inflammation/Uterus/Uterus_py/scvi/" + "Uterus_scvi_cell_level1.h5ad")


sc.pl.umap(combined_adata, color="cell_type_level1",
           legend_loc="on data",
           palette=cell_color,
           add_outline=False,
           legend_fontsize=8,
           legend_fontoutline=0.1,
           legend_fontweight="normal",
           title="Uterus_scvi_cell_annotation",
           frameon=False,
           save = "_Uterus_scvi_cell_level1_clear.pdf"
           )
sc.pl.umap(combined_adata, color="cell_type_level1",
           #legend_loc="on data",
           palette=cell_color,
           add_outline=False,
           legend_fontsize=8,
           legend_fontoutline=0.1,
           legend_fontweight="normal",
           title="Uterus_scvi_cell_annotation",
           frameon=False,
           save = "_Uterus_scvi_cell_level1_legend.pdf"
           )

# scanvi-----------------------

lvae = scvi.model.SCANVI.from_scvi_model(
    vae,
    adata=combined_adata,
    labels_key="cell_type_level1",
    unlabeled_category="Unknown",
)
lvae.train()
model_dir2 = os.path.join("/storage/data/KAI/Pan_Inflammation/Uterus/uterus_py/scanvi/model/uterus_allcells/", "uterus_allcells_scanvi_model")
lvae.save(model_dir2, overwrite=True, prefix="bladder_allcells_scanvi_")

lvae = scvi.model.SCANVI.load("/storage/data/KAI/Pan_Inflammation/Uterus/Uterus_py/scanvi/model/uterus_allcells/uterus_allcells_scanvi_model/", adata=combined_adata, prefix="uterus_allcells_scanvi_")


combined_adata.obsm["X_scANVI"] = lvae.get_latent_representation()

# run PcA then generate UMAp plots
sc.tl.pca(combined_adata)
sc.pl.pca_variance_ratio(combined_adata,log=True)
# use scANVI latent space for UMAp generation
sc.pp.neighbors(combined_adata,use_rep="X_scANVI", n_neighbors=30)
sc.tl.umap(combined_adata,min_dist=0.3)


fig, ax=plt.subplots(1,1,figsize=(5,4))
sc.pl.umap(combined_adata, color="cell_type_level1",
           legend_loc="right margin",
           palette=cell_color,
           add_outline=False,
           legend_fontsize=8,
           legend_fontoutline=0.1,
           legend_fontweight="normal",
           title="Uterus_scanvi_cell_annotation",
           frameon=False,
           show=False,
           ax=ax
           #save = "_Uterus_scanvi_cell_level1_legend.pdf"
           )
# fig = plt.gcf()
for collection in ax.collections:
    collection.set_rasterized(True)
# for ax in fig.get_axes():

plt.savefig("/storage/data/KAI/Pan_Inflammation/Infla_py/figures/umap_Uterus_scanvi_cell_level1_legend.pdf", format="pdf", dpi=600)
plt.close()



fig, ax=plt.subplots(1,1,figsize=(5,4))
sc.pl.umap(combined_adata, color="cell_type_level1",
           legend_loc="on data",
           palette=cell_color,
           add_outline=False,
           legend_fontsize=8,
           legend_fontoutline=0.1,
           legend_fontweight="normal",
           title="Uterus_scanvi_cell_annotation",
           frameon=False,
           show=False,
           ax=ax
           #save = "_Uterus_scanvi_cell_level1_clear.pdf"
           )
# fig = plt.gcf()
for collection in ax.collections:
    collection.set_rasterized(True)
# for ax in fig.get_axes():

plt.savefig("/storage/data/KAI/Pan_Inflammation/Infla_py/figures/umap_Uterus_scanvi_cell_level1_clear.pdf", format="pdf", dpi=600)
plt.close()


#combined_adata.write("/storage/data/KAI/Pan_Inflammation/Uterus/Uterus_py/scanvi/Uterus_scanvi.h5ad")


#combined_adata = sc.read_h5ad("/storage/data/KAI/Pan_Inflammation/Uterus/Uterus_py/scanvi/Uterus_scanvi.h5ad")


#violin_marker
violin_markers = [
  "EPCAM", "KRT8","KRT18",
  "DCN", "LUM","COL3A1",
  "PDGFRB", "RGS5",
  "PECAM1","VWF","CLDN5",
  "C1QA","AIF1","LYZ",
  "S100A8","S100A9",
  "CPA3","TPSAB1", "CTSG",
  "TRAC","CD3D", "CD3E",
  "MS4A1","CD79A","BANK1",
  "IGKC","IGHG1", "JCHAIN",
  "NRXN1", "PLP1","S100B"
                  ]


#batch effector

sc.pl.umap(combined_adata,color=["GEO"], ncols=2, frameon=False, save = "Uterus_scanvi_batch.png")


# TO R
combined_adata = sc.read_h5ad("/storage/data/KAI/Pan_Inflammation/Uterus/Uterus_py/scanvi/Uterus_scanvi.h5ad")
combined_adata.X = combined_adata.layers['counts']
combined_adata.write("/storage/data/KAI/Pan_Inflammation/Uterus/Uterus_py/scanvi/Uterus_scanvi_toR.h5ad")








###Colon Inflammation scRNA-seq

adata_list = []

# Crohn’s disease

preprecess_path = "/storage/data/KAI/Pan_Inflammation/Colorectum/Inflammatory bowel disease/Crohn’s disease/1.GSE214695/GSE214695_py/GSE214695_preprocess.h5ad"
adata = sc.read_h5ad(preprecess_path)
print(f'GSE214695 dim: {adata.shape}')
adata_list.append(adata)


preprecess_path = "/storage/data/KAI/Pan_Inflammation/Colorectum/Inflammatory bowel disease/Crohn’s disease/8.GSE225199/GSE225199_py/GSE225199_preprocess.h5ad"
adata = sc.read_h5ad(preprecess_path)
print(f'GSE225199 dim: {adata.shape}')
adata_list.append(adata)


preprecess_path = "/storage/data/KAI/Pan_Inflammation/Colorectum/Inflammatory bowel disease/Crohn’s disease/9.GSE215001/GSE215001_py/GSE215001_preprocess.h5ad"
adata = sc.read_h5ad(preprecess_path)
print(f'GSE215001 dim: {adata.shape}')
adata_list.append(adata)

preprecess_path = "/storage/data/KAI/Pan_Inflammation/Colorectum/Inflammatory bowel disease/Crohn’s disease/10.GSE202052/GSE202052_py/GSE202052_preprocess.h5ad"
adata = sc.read_h5ad(preprecess_path)
print(f'GSE202052 dim: {adata.shape}')
adata_list.append(adata)


# preprecess_path = "/storage/data/KAI/Pan_Inflammation/Colorectum/Inflammatory bowel disease/Crohn’s disease/13.GSE156776/GSE156776_py/GSE156776_preprocess.h5ad"
# adata = sc.read_h5ad(preprecess_path)
# print(f'GSE156776 dim: {adata.shape}')
# adata_list.append(adata)


preprecess_path = "/storage/data/KAI/Pan_Inflammation/Colorectum/Inflammatory bowel disease/Crohn’s disease/16.GSE246987/GSE246987_py/GSE246987_preprocess.h5ad"
adata = sc.read_h5ad(preprecess_path)
print(f'GSE246987 dim: {adata.shape}')
adata_list.append(adata)


preprecess_path = "/storage/data/KAI/Pan_Inflammation/Colorectum/Inflammatory bowel disease/Crohn’s disease/19.GSE134809/GSE134809_py/GSE134809_preprocess.h5ad"
adata = sc.read_h5ad(preprecess_path)
print(f'GSE134809 dim: {adata.shape}')
adata_list.append(adata)


# Ulcerative colitis
preprecess_path = "/storage/data/KAI/Pan_Inflammation/Colorectum/Inflammatory bowel disease/Ulcerative colitis/4.GSE231993/GSE231993_py/GSE231993_preprocess.h5ad"
adata = sc.read_h5ad(preprecess_path)
print(f'GSE231993 dim: {adata.shape}')
adata_list.append(adata)

preprecess_path = "/storage/data/KAI/Pan_Inflammation/Colorectum/Inflammatory bowel disease/Ulcerative colitis/8.GSE182270/GSE182270_py/GSE182270_preprocess.h5ad"
adata = sc.read_h5ad(preprecess_path)
print(f'GSE182270 dim: {adata.shape}')
adata_list.append(adata)

preprecess_path = "/storage/data/KAI/Pan_Inflammation/Colorectum/Inflammatory bowel disease/Ulcerative colitis/15.GSE242086/GSE242086_py/GSE242086_preprocess.h5ad"
adata = sc.read_h5ad(preprecess_path)
print(f'GSE242086 dim: {adata.shape}')
adata_list.append(adata)

# preprecess_path = "/storage/data/KAI/Pan_Inflammation/Colorectum/Inflammatory bowel disease/Ulcerative colitis/17.GSE150115/GSE150115_py/GSE150115_preprocess.h5ad"
# adata = sc.read_h5ad(preprecess_path)
# print(f'GSE150115 dim: {adata.shape}')
# adata_list.append(adata)


# Immune-checkpoint inhibitor colitis
preprecess_path = "/storage/data/KAI/Pan_Inflammation/Colorectum/Inflammatory bowel disease/Immune-checkpoint inhibitor colitis/1.GSE253720/GSE253720_py/GSE253720_preprocess.h5ad"
adata = sc.read_h5ad(preprecess_path)
print(f'GSE253720 dim: {adata.shape}')
adata_list.append(adata)


# 检查交集基因
len(adata_list)

var_list = []
for i in range(0,len(adata_list)):
    var_gse = list(adata_list[i].var_names)
    var_list.append(var_gse)


inte_adata = []
for i in range(1,len(var_list)+1):
    intesect = set.intersection(*[set(v) for v in var_list[:i]])
    inte_adata.append(intesect)

for i, inter in enumerate(inte_adata):
    print(f"交集（前 {i+1} 个集合）: {len(inter)}")


# Integration
combined_adata = adata_list[0].concatenate(*adata_list[1:], join='outer', batch_key=None, batch_categories=None) # 整合后 (736281, 16331)

temp = combined_adata

# bool转0-1
#temp.obs = temp.obs.astype({'predicted_doublet': 'int'})

# 转换object类型
temp.obs = temp.obs.astype({'Age': 'str'}) #数据中不能含有缺失值(NaN)
temp.obs = temp.obs.astype({'Sample_title': 'str'})
temp.obs = temp.obs.astype({'Donor': 'str'})
# bool转0-1
#temp.var = temp.var.astype({'mt': 'int'})
#temp.var = temp.var.astype({'ribo': 'int'})
#temp.var = temp.var.astype({'hb': 'int'})
#temp.var = temp.var.astype({'highly_variable': 'int'})

# 转换object类型
temp.var = temp.var.astype({'Gene_ids-0': 'str'})



# 删除var中多余的属性
rubb_var = ['n_cells_by_counts-0', 'mean_counts-0', 'log1p_mean_counts-0', 'pct_dropout_by_counts-0', 'total_counts-0', 'log1p_total_counts-0', 'n_cells-0', 'gene_ids-1', 'n_cells_by_counts-1', 'mean_counts-1', 'log1p_mean_counts-1', 'pct_dropout_by_counts-1', 'total_counts-1', 'log1p_total_counts-1', 'n_cells-1', 'gene_ids-2', 'n_cells_by_counts-2', 'mean_counts-2', 'log1p_mean_counts-2', 'pct_dropout_by_counts-2', 'total_counts-2', 'log1p_total_counts-2', 'n_cells-2', 'n_cells_by_counts-3', 'mean_counts-3', 'log1p_mean_counts-3', 'pct_dropout_by_counts-3', 'total_counts-3', 'log1p_total_counts-3', 'n_cells-3', 'Gene_ids-3', 'gene_ids-4', 'n_cells_by_counts-4', 'mean_counts-4', 'log1p_mean_counts-4', 'pct_dropout_by_counts-4', 'total_counts-4', 'log1p_total_counts-4', 'n_cells-4', 'gene_ids-5', 'n_cells_by_counts-5', 'mean_counts-5', 'log1p_mean_counts-5', 'pct_dropout_by_counts-5', 'total_counts-5', 'log1p_total_counts-5', 'n_cells-5', 'gene_ids-6', 'n_cells_by_counts-6', 'mean_counts-6', 'log1p_mean_counts-6', 'pct_dropout_by_counts-6', 'total_counts-6', 'log1p_total_counts-6', 'n_cells-6', 'gene_ids-7', 'n_cells_by_counts-7', 'mean_counts-7', 'log1p_mean_counts-7', 'pct_dropout_by_counts-7', 'total_counts-7', 'log1p_total_counts-7', 'n_cells-7', 'gene_ids-8', 'n_cells_by_counts-8', 'mean_counts-8', 'log1p_mean_counts-8', 'pct_dropout_by_counts-8', 'total_counts-8', 'log1p_total_counts-8', 'n_cells-8', 'n_cells_by_counts-9', 'mean_counts-9', 'log1p_mean_counts-9', 'pct_dropout_by_counts-9', 'total_counts-9', 'log1p_total_counts-9', 'n_cells-9', 'assay-9']
temp.var.drop(columns=rubb_var, inplace=True)

temp.var.rename(columns={'Gene_ids-0':'Gene_ids'}, inplace=True)





# scvi-----------------------
combined_adata = temp



combined_adata.obs.loc[combined_adata.obs['Sample_geo_accession'] == 'GSM7748909', 'Group'] = 'Normal'
combined_adata.obs.loc[combined_adata.obs['Sample_geo_accession'] == 'GSM7748910', 'Group'] = 'Normal'
combined_adata.obs.loc[combined_adata.obs['Sample_geo_accession'] == 'GSM7748911', 'Group'] = 'Normal'
combined_adata.obs.loc[combined_adata.obs['Sample_geo_accession'] == 'GSM7748912', 'Group'] = 'Normal'
combined_adata.obs.loc[combined_adata.obs['Sample_geo_accession'] == 'GSM7748913', 'Group'] = 'Normal'
combined_adata.obs.loc[combined_adata.obs['Sample_geo_accession'] == 'GSM7748914', 'Group'] = 'Normal'


combined_adata.obs.loc[combined_adata.obs['Sample_geo_accession'] == 'GSM7748927', 'Group'] = 'Normal'
combined_adata.obs.loc[combined_adata.obs['Sample_geo_accession'] == 'GSM7748928', 'Group'] = 'Normal'
combined_adata.obs.loc[combined_adata.obs['Sample_geo_accession'] == 'GSM7748929', 'Group'] = 'Normal'
combined_adata.obs.loc[combined_adata.obs['Sample_geo_accession'] == 'GSM7748930', 'Group'] = 'Normal'
combined_adata.obs.loc[combined_adata.obs['Sample_geo_accession'] == 'GSM7748931', 'Group'] = 'Normal'
combined_adata.obs.loc[combined_adata.obs['Sample_geo_accession'] == 'GSM7748932', 'Group'] = 'Normal'


combined_adata = combined_adata[combined_adata.obs['Group'] != 'UC-IPAA'].copy()
# 写出整合后文件(未进行scvi)
#combined_adata.write("/storage/data/KAI/Pan_Inflammation/Colorectum/Colorectum_py/scvi/" + "Colorectum_Inte.h5ad")

#combined_adata = sc.read_h5ad("/storage/data/KAI/Pan_Inflammation/Colorectum/Colorectum_py/scvi/" + "Colorectum_Inte.h5ad")

combined_adata.raw = combined_adata  # keep full dimension safe

#combined_adata.layers["counts"] = combined_adata.X

sc.pp.highly_variable_genes(combined_adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

scvi.model.SCVI.setup_anndata(combined_adata, layer="counts", categorical_covariate_keys=["Sample_geo_accession", "Sample_Type", "Technology", "Sample_Site", "GEO"])

vae = scvi.model.SCVI(combined_adata, n_latent=20, n_layers=2, dropout_rate=0.1)
vae.train(max_epochs=30, early_stopping=True, batch_size=256)
# 保存scvi模型
model_dir = os.path.join("/storage/data/KAI/Pan_Inflammation/Colorectum/Colorectum_py/scvi/model/colorectum_allcells/", "colorectum_allcells_scvi_model")
vae.save(model_dir, overwrite=True, prefix="colorectum_allcells_scvi_")

vae = scvi.model.SCVI.load("/storage/data/KAI/Pan_Inflammation/Colorectum/Colorectum_py/scvi/model/colorectum_allcells/colorectum_allcells_scvi_model/", adata=combined_adata, prefix="colorectum_allcells_scvi_")

combined_adata.obsm["X_scVI"] = vae.get_latent_representation()
combined_adata.layers["X_normalized_scVI"] = vae.get_normalized_expression(library_size=10e4)


# run PcA then generate UMAp plots
sc.tl.pca(combined_adata)
sc.pl.pca_variance_ratio(combined_adata,log=True)
# use scvI latent space for UMAp generation
sc.pp.neighbors(combined_adata,use_rep="X_scVI")
sc.tl.umap(combined_adata,min_dist=0.3)


#result_folder = "/storage/data/KAI/Pan_Inflammation/Colorectum/Colorectum_py/scvi/result/"

sc.pl.umap(combined_adata,color=["GEO","Technology","Group","Sample_Site"], ncols=2, frameon=False, save = "Colorectum_scvi_batch.png")

#写出
combined_adata.write("/storage/data/KAI/Pan_Inflammation/Colorectum/Colorectum_py/scvi/" + "Colorectum_scvi.h5ad")

#combined_adata = sc.read_h5ad("/storage/data/KAI/Pan_Inflammation/Colorectum/Colorectum_py/scvi/" + "Colorectum_scvi.h5ad")

# full gene ---------------------------
combined_adata.raw.var.index = combined_adata.raw.var.index.astype(str)  # 确保索引为字符串
combined_adata.raw.var = combined_adata.raw.var.astype(str)  # 确保所有列的类型为字符串

rubb_var = [
'mt-0', 'ribo-0', 'hb-0', 'n_cells_by_counts-0', 'mean_counts-0', 'log1p_mean_counts-0', 'pct_dropout_by_counts-0', 'total_counts-0', 'log1p_total_counts-0', 'n_cells-0', 'gene_ids-1', 'mt-1', 'ribo-1', 'hb-1', 'n_cells_by_counts-1', 'mean_counts-1', 'log1p_mean_counts-1', 'pct_dropout_by_counts-1', 'total_counts-1', 'log1p_total_counts-1', 'n_cells-1', 'gene_ids-2', 'mt-2', 'ribo-2', 'hb-2', 'n_cells_by_counts-2', 'mean_counts-2', 'log1p_mean_counts-2', 'pct_dropout_by_counts-2', 'total_counts-2', 'log1p_total_counts-2', 'n_cells-2', 'mt-3', 'ribo-3', 'hb-3', 'n_cells_by_counts-3', 'mean_counts-3', 'log1p_mean_counts-3', 'pct_dropout_by_counts-3', 'total_counts-3', 'log1p_total_counts-3', 'n_cells-3', 'Gene_ids-3', 'gene_ids-4', 'mt-4', 'ribo-4', 'hb-4', 'n_cells_by_counts-4', 'mean_counts-4', 'log1p_mean_counts-4', 'pct_dropout_by_counts-4', 'total_counts-4', 'log1p_total_counts-4', 'n_cells-4', 'gene_ids-5', 'mt-5', 'ribo-5', 'hb-5', 'n_cells_by_counts-5', 'mean_counts-5', 'log1p_mean_counts-5', 'pct_dropout_by_counts-5', 'total_counts-5', 'log1p_total_counts-5', 'n_cells-5', 'gene_ids-6', 'mt-6', 'ribo-6', 'hb-6', 'n_cells_by_counts-6', 'mean_counts-6', 'log1p_mean_counts-6', 'pct_dropout_by_counts-6', 'total_counts-6', 'log1p_total_counts-6', 'n_cells-6', 'gene_ids-7', 'mt-7', 'ribo-7', 'hb-7', 'n_cells_by_counts-7', 'mean_counts-7', 'log1p_mean_counts-7', 'pct_dropout_by_counts-7', 'total_counts-7', 'log1p_total_counts-7', 'n_cells-7', 'gene_ids-8', 'mt-8', 'ribo-8', 'hb-8', 'n_cells_by_counts-8', 'mean_counts-8', 'log1p_mean_counts-8', 'pct_dropout_by_counts-8', 'total_counts-8', 'log1p_total_counts-8', 'n_cells-8', 'mt-9', 'ribo-9', 'hb-9', 'n_cells_by_counts-9', 'mean_counts-9', 'log1p_mean_counts-9', 'pct_dropout_by_counts-9', 'total_counts-9', 'log1p_total_counts-9', 'n_cells-9', 'assay-9'
]
combined_adata.var.drop(columns=rubb_var, inplace=True)

combined_adata.obs = combined_adata.obs.astype({'Age': 'str'})
combined_adata.obs = combined_adata.obs.astype({'Donor': 'str'})
combined_adata.obs = combined_adata.obs.astype({'Sample_title': 'str'})
combined_adata.write("/storage/data/KAI/Pan_Inflammation/new_inflammation/Colorectum/data/Colorectum_all_gene.h5ad")

# 细胞注释-----------------------
sc.tl.leiden(combined_adata, resolution=0.8)
fig, ax = plt.subplots(figsize=(10, 8), dpi=300)
sc.pl.umap(combined_adata, color=["leiden"],
           legend_loc='on data',
           ax=ax,
           save="_colorectum_leiden.png"
           )



Epi_marker = ["EPCAM","KRT8","KRT18", "KRT19"]

Fib_marker = ["DCN","LUM","COL1A2", "COL3A1", "COL6A1", "PDPN"]

Mural_marker = ["PDGFRB", "RGS5", "KCNJ8"]

Endo_marker = ["PECAM1","VWF","PLVAP","CLDN5"]

Mac_marker = ["AIF1","CSF1R","LYZ", "CD163", "FCGR2A", "CSF1R", "C1QA", "C1QB"]

Mono_marker = ["S100A9","S100A8"]

Mast_marker = ["TPSB2","CPA3","TPSAB1"]

Dendritic_marker = ["ITGAX","CLEC9A", "IRF8", "CD40", "CD80", "CD83", "CCR7"]

Neutrophil_marker =["G0S2", "CXCL8", "SOD2", "NAMPT", "S100A8"]

TNK_marker =  ["TRAC", "CD3D", "CD3E", "CD3G", "CD2", "NKG7", "TRDC", "XCL2", "XCL1"]

B_marker = ["MS4A1","CD79A","CD79B", "CD19", "BLNK", "FCRL5", "CD24", "IGHG1", "CD27", "BANK1"]

Plasma_marker = ["IGHA1","IGKC", "JCHAIN"]

#Erythrocyte_marker = ["HBA1","HBA2","HBB"]

Neural_marker = ["NRXN1","S100B","NGFR"]


sc.pl.violin(combined_adata, Epi_marker, groupby="leiden", stripplot=False, show=True, use_raw=False, rotation=90)
sc.pl.violin(combined_adata, Fib_marker, groupby="leiden", stripplot=False, show=True, use_raw=False, rotation=90)
sc.pl.violin(combined_adata, Mural_marker, groupby="leiden", stripplot=False, show=True, use_raw=False, rotation=90)
sc.pl.violin(combined_adata, Endo_marker, groupby="leiden", stripplot=False, show=True, use_raw=False, rotation=90)
sc.pl.violin(combined_adata, Mac_marker, groupby="leiden", stripplot=False, show=True, use_raw=False, rotation=90)
sc.pl.violin(combined_adata, Mono_marker, groupby="leiden", stripplot=False, show=True, use_raw=False, rotation=90)
sc.pl.violin(combined_adata, Mast_marker, groupby="leiden", stripplot=False, show=True, use_raw=False, rotation=90)
sc.pl.violin(combined_adata, Dendritic_marker, groupby="leiden", stripplot=False, show=True, use_raw=False, rotation=90)
sc.pl.violin(combined_adata, Neutrophil_marker, groupby="leiden", stripplot=False, show=True, use_raw=False, rotation=90)
sc.pl.violin(combined_adata, TNK_marker, groupby="leiden", stripplot=False, show=True, use_raw=False, rotation=90)
sc.pl.violin(combined_adata, B_marker, groupby="leiden", stripplot=False, show=True, use_raw=False, rotation=90)
sc.pl.violin(combined_adata, Plasma_marker, groupby="leiden", stripplot=False, show=True, use_raw=False, rotation=90)
#sc.pl.violin(combined_adata, Erythrocyte_marker, groupby="leiden", stripplot=False, show=True, use_raw=False, rotation=90)
sc.pl.violin(combined_adata, Neural_marker, groupby="leiden", stripplot=False, show=True, use_raw=False, rotation=90)



sc.tl.rank_genes_groups(combined_adata, "leiden", groups=["38"], method="wilcoxon")
result = combined_adata.uns['rank_genes_groups']
cluster_38_genes = pd.DataFrame({
    'gene': result['names']['38'],
    'logfoldchanges': result['logfoldchanges']['38'],
    'pvals': result['pvals']['38'],
    'pvals_adj': result['pvals_adj']['38']
})
cluster_38_genes['gene'].head(30).to_list()


cell_annotation_level1 = {
"0":	"T/NK cell",
"1":	"Plasma cell",
"2":	"B cell",
"3":	"T/NK cell",
"4":	"B cell",
"5":	"T/NK cell",
"6":	"Epithelial cell",
"7":	"Fibroblast",
"8":	"Macrophage",
"9":	"T/NK cell",
"10":	"Plasma cell",
"11":	"Epithelial cell",
"12":	"Epithelial cell",
"13":	"Epithelial cell",
"14":	"Epithelial cell",
"15":	"Plasma cell",
"16":	"Endothelial cell",
"17":	"B cell",
"18":	"Plasma cell",
"19":	"B cell",
"20":	"B cell",
"21":	"Neutrophil",
"22":	"Mast cell",
"23":	"Mural cell",
"24":	"B cell",
"25":	"T/NK cell",
"26":	"Epithelial cell",
"27":	"T/NK cell",
"28":	"Plasma cell",
"29":	"Epithelial cell",
"30":	"Epithelial cell",
"31":	"Neural cell",
"32":	"Epithelial cell",
"33":	"Endothelial cell",
"34":	"Epithelial cell",
"35":	"Dendritic cell",
"36":	"B cell",
"37":	"T/NK cell",
"38":	"B cell"

}
leiden_annotation_level1 = pd.Series(cell_annotation_level1)
combined_adata.obs['cell_type_level1'] = combined_adata.obs['leiden'].map(leiden_annotation_level1)

cell_color = {
    "Epithelial cell":"#e43030",
    "Fibroblast":"#ff8831",
    "Mural cell":"#20452e",
    "Endothelial cell":"#704ba3",
    "Macrophage":"#67a8cd",
    "Monocyte":"#8184e2",
    "Mast cell":"#ff9d9f",
    "Dendritic cell":"#cf9f88",
    "Neutrophil":"#ffc17f",
    "T/NK cell":"#50aa4b",
    "B cell":"#b3e19b",
    "Plasma cell":"#ab3181",
    "Erythrocyte":"#a1acbd",
    "Neural cell":"#00afba"
}

combined_adata = sc.read_h5ad("/storage/data/KAI/Pan_Inflammation/Colorectum/Colorectum_py/scvi/" + "Colorectum_scvi_cell_level1.h5ad")

#由于10x测序平台合理性，将"Neutrophil"替换为"Monocyte"
print(combined_adata.obs['cell_type_level1'].unique())
combined_adata.obs['cell_type_level1'] = combined_adata.obs['cell_type_level1'].replace("Neutrophil", "Monocyte")
print(combined_adata.obs['cell_type_level1'].unique())

#combined_adata.write("/storage/data/KAI/Pan_Inflammation/Colorectum/Colorectum_py/scvi/" + "Colorectum_scvi_cell_level1.h5ad")






sc.pl.umap(combined_adata, color="cell_type_level1",
           legend_loc="on data",
           palette=cell_color,
           add_outline=False,
           legend_fontsize=8,
           legend_fontoutline=0.1,
           legend_fontweight="normal",
           title="Colorectum_scvi_cell_annotation",
           frameon=False,
           save = "_Colorectum_scvi_cell_level1_clear.pdf"
           )
sc.pl.umap(combined_adata, color="cell_type_level1",
           #legend_loc="on data",
           palette=cell_color,
           add_outline=False,
           legend_fontsize=8,
           legend_fontoutline=0.1,
           legend_fontweight="normal",
           title="Colorectum_scvi_cell_annotation",
           frameon=False,
           save = "_Colorectum_scvi_cell_level1_legend.pdf"
           )

# scanvi-----------------------

lvae = scvi.model.SCANVI.from_scvi_model(
    vae,
    adata=combined_adata,
    labels_key="cell_type_level1",
    unlabeled_category="Unknown",
)
lvae.train()
model_dir2 = os.path.join("/storage/data/KAI/Pan_Inflammation/Colorectum/Colorectum_py/scanvi/model/colorectum_allcells/", "colorectum_allcells_scanvi_model")
lvae.save(model_dir2, overwrite=True, prefix="colorectum_allcells_scanvi_")

lvae = scvi.model.SCANVI.load("/storage/data/KAI/Pan_Inflammation/Colorectum/Colorectum_py/scanvi/model/colorectum_allcells/colorectum_allcells_scanvi_model/", adata=combined_adata, prefix="colorectum_allcells_scanvi_")


combined_adata.obsm["X_scANVI"] = lvae.get_latent_representation()

# run PcA then generate UMAp plots
sc.tl.pca(combined_adata)
sc.pl.pca_variance_ratio(combined_adata,log=True)
# use scANVI latent space for UMAp generation
sc.pp.neighbors(combined_adata,use_rep="X_scANVI", n_neighbors=30)
sc.tl.umap(combined_adata,min_dist=0.3)



#combined_adata = sc.read_h5ad("/storage/data/KAI/Pan_Inflammation/Colorectum/Colorectum_py/scanvi/"+"Colorectum_scanvi.h5ad")

#由于10x测序平台合理性，将"Neutrophil"替换为"Monocyte"
print(combined_adata.obs['cell_type_level1'].unique().tolist())
combined_adata.obs['cell_type_level1'] = combined_adata.obs['cell_type_level1'].replace("Neutrophil", "Monocyte")
print(combined_adata.obs['cell_type_level1'].unique().tolist())

combined_adata.write("/storage/data/KAI/Pan_Inflammation/Colorectum/Colorectum_py/scanvi/"+"Colorectum_scanvi.h5ad")






fig, ax=plt.subplots(1,1,figsize=(5,4))
sc.pl.umap(combined_adata, color="cell_type_level1",
           legend_loc="right margin",
           palette=cell_color,
           add_outline=False,
           legend_fontsize=8,
           legend_fontoutline=0.1,
           legend_fontweight="normal",
           title="Colorectum_scanvi_cell_annotation",
           frameon=False,
           show=False,
           ax=ax
           #save = "_Colorectum_scanvi_cell_level1_legend.pdf"
           )
# fig = plt.gcf()
for collection in ax.collections:
    collection.set_rasterized(True)
# for ax in fig.get_axes():

plt.savefig("/storage/data/KAI/Pan_Inflammation/Infla_py/figures/umap_Colorectum_scanvi_cell_level1_legend.pdf", format="pdf", dpi=600)
plt.close()



fig, ax=plt.subplots(1,1,figsize=(5,4))
sc.pl.umap(combined_adata, color="cell_type_level1",
           legend_loc="on data",
           palette=cell_color,
           add_outline=False,
           legend_fontsize=8,
           legend_fontoutline=0.1,
           legend_fontweight="normal",
           title="Colorectum_scanvi_cell_annotation",
           frameon=False,
           show=False,
           ax=ax
           #save = "_Colorectum_scanvi_cell_level1_clear.pdf"
           )
# fig = plt.gcf()
for collection in ax.collections:
    collection.set_rasterized(True)
# for ax in fig.get_axes():

plt.savefig("/storage/data/KAI/Pan_Inflammation/Infla_py/figures/umap_Colorectum_scanvi_cell_level1_clear.pdf", format="pdf", dpi=600)
plt.close()



combined_adata.write("/storage/data/KAI/Pan_Inflammation/Colorectum/Colorectum_py/scanvi/"+"Colorectum_scanvi.h5ad")


#batch effector

sc.pl.umap(combined_adata,color=["GEO"], ncols=2, frameon=False, save = "Colorectum_scanvi_batch.png")

#violin_marker
violin_markers = [
  "EPCAM", "KRT8","KRT18",
  "DCN", "LUM","COL3A1",
  "PDGFRB", "RGS5",
  "PECAM1","VWF","CLDN5",
  "C1QA","AIF1","LYZ",
  "S100A8","S100A9",
  "IRF8", "PLD4", "LILRA4",
  "CPA3","TPSAB1", "TPSB2",
  "TRAC","CD3D", "CD3E",
  "MS4A1","CD79A","BANK1",
  "IGHA1","IGKC", "JCHAIN",
  "NRXN1", "PLP1","S100B"
                  ]

combined_adata = sc.read_h5ad("/storage/data/KAI/Pan_Inflammation/Colorectum/Colorectum_py/scanvi/"+"Colorectum_scanvi.h5ad")
combined_adata.X = combined_adata.layers['counts']
combined_adata.write("/storage/data/KAI/Pan_Inflammation/Colorectum/Colorectum_py/scanvi/"+"Colorectum_scanvi_toR.h5ad")



# Module
import matplotlib.pyplot as plt
import scanpy as sc
import pandas as pd
import numpy as np
import subprocess
import os
import re
import scvi
import matplotlib.pyplot as plt
from scib_metrics.benchmark import Benchmarker

### Joint Inflammation scRNA-seq
adata_list = []

# Rheumatoid_arthritis
preprecess_path = "/storage/data/KAI/Pan_Inflammation/Joint/Rheumatoid_arthritis/4.GSE200815/GSE200815_py/GSE200815_preprocess.h5ad"
adata = sc.read_h5ad(preprecess_path)
print(f'GSE200815 dim: {adata.shape}')
adata_list.append(adata)

# IA中已涵盖该数据集
# preprecess_path = "/storage/data/KAI/Pan_Inflammation/Joint/Rheumatoid_arthritis/5.GSE181082/GSE181082_py/GSE181082_preprocess.h5ad"
# adata = sc.read_h5ad(preprecess_path)
# print(f'GSE181082 dim: {adata.shape}')
# adata_list.append(adata)


# Osteoarthritis
preprecess_path = "/storage/data/KAI/Pan_Inflammation/Joint/Osteoarthritis/3.GSE152805/GSE152805_py/GSE152805_preprocess.h5ad"
adata = sc.read_h5ad(preprecess_path)
print(f'GSE152805 dim: {adata.shape}')
adata_list.append(adata)

preprecess_path = "/storage/data/KAI/Pan_Inflammation/Joint/Osteoarthritis/5.GSE216651/GSE216651_py/GSE216651_preprocess.h5ad"
adata = sc.read_h5ad(preprecess_path)
unique_genes = pd.Series(adata.var_names).drop_duplicates(keep="first")
adata = adata[adata.obs.index, unique_genes.index]
adata.var_names = unique_genes.values
print(f'GSE216651 dim: {adata.shape}')
adata_list.append(adata)


# Inflammatory arthritis
preprecess_path = "/storage/data/KAI/Pan_Inflammation/Joint/Inflammatory arthritis/1.E-MTAB-11791/E-MTAB-11791_py/E-MTAB-11791_preprocess.h5ad"
adata = sc.read_h5ad(preprecess_path)
unique_genes = pd.Series(adata.var_names).drop_duplicates(keep="first")
adata = adata[adata.obs.index, unique_genes.index]
adata.var_names = unique_genes.values
print(f'E-MTAB-11791 dim: {adata.shape}')
adata_list.append(adata)

# Normal
preprecess_path = "/storage/data/KAI/Pan_Inflammation/Joint/Normal/1.E-MTAB-14339/E-MTAB-14339_py/E-MTAB-14339_preprocess.h5ad"
adata = sc.read_h5ad(preprecess_path)
print(f'E-MTAB-14339 dim: {adata.shape}')
adata_list.append(adata)




len(adata_list)

var_list = []
for i in range(0,len(adata_list)):
    var_gse = list(adata_list[i].var_names)
    var_list.append(var_gse)


inte_adata = []
for i in range(1,len(var_list)+1):
    intesect = set.intersection(*[set(v) for v in var_list[:i]])
    inte_adata.append(intesect)

for i, inter in enumerate(inte_adata):
    print(f"交集（前 {i+1} 个集合）: {len(inter)}")






# Integration
combined_adata = adata_list[0].concatenate(adata_list[1:], join='outer', batch_key=None, batch_categories=None) # 整合后 (386769, 14452)


temp = combined_adata



# 转换object类型
temp.obs = temp.obs.astype({'Age': 'str'}) #数据中不能含有缺失值(NaN)

temp.obs = temp.obs.astype({'Donor': 'str'})


# 删除var中多余的属性

rubb_var = ['n_cells_by_counts-0', 'mean_counts-0', 'log1p_mean_counts-0', 'pct_dropout_by_counts-0', 'total_counts-0', 'log1p_total_counts-0', 'n_cells-0', 'gene_ids-1', 'n_cells_by_counts-1', 'mean_counts-1', 'log1p_mean_counts-1', 'pct_dropout_by_counts-1', 'total_counts-1', 'log1p_total_counts-1', 'n_cells-1', 'n_cells_by_counts-2', 'mean_counts-2', 'log1p_mean_counts-2', 'pct_dropout_by_counts-2', 'total_counts-2', 'log1p_total_counts-2', 'n_cells-2', 'Gene_ids-2', 'n_cells_by_counts-3', 'mean_counts-3', 'log1p_mean_counts-3', 'pct_dropout_by_counts-3', 'total_counts-3', 'log1p_total_counts-3', 'n_cells-3', 'gene_ids-4', 'n_cells_by_counts-4', 'mean_counts-4', 'log1p_mean_counts-4', 'pct_dropout_by_counts-4', 'total_counts-4', 'log1p_total_counts-4', 'n_cells-4']
temp.var.drop(columns=rubb_var, inplace=True)

temp.var.rename(columns={'gene_ids-0':'Gene_ids'}, inplace=True)


# 写出整合后文件(未进行scvi)
#temp.write("/storage/data/KAI/Pan_Inflammation/Joint/Joint_py/scvi/" + "Joint_Inte.h5ad")


# scvi-----------------------
combined_adata = temp

#combined_adata = sc.read_h5ad("/storage/data/KAI/Pan_Inflammation/Joint/Joint_py/scvi/" + "Joint_Inte.h5ad")

combined_adata.raw = combined_adata  # keep full dimension safe


sc.pp.highly_variable_genes(combined_adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

scvi.model.SCVI.setup_anndata(combined_adata, layer="counts", categorical_covariate_keys=["Sample_geo_accession", "Sample_Type", "Sample_Site","GEO"])
vae = scvi.model.SCVI(combined_adata, n_latent=20, n_layers=2, dropout_rate=0.1)
vae.train(max_epochs=50, early_stopping=True, batch_size=256)
# 保存scvi模型
model_dir = os.path.join("/storage/data/KAI/Pan_Inflammation/Joint/Joint_py/scvi/model/joint_allcells/", "joint_allcells_scvi_model")
vae.save(model_dir, overwrite=True, prefix="joint_allcells_scvi_")

vae = scvi.model.SCVI.load("/storage/data/KAI/Pan_Inflammation/Joint/Joint_py/scvi/model/joint_allcells/joint_allcells_scvi_model/", adata=combined_adata, prefix="joint_allcells_scvi_")

combined_adata.obsm["X_scVI"] = vae.get_latent_representation()
combined_adata.layers["X_normalized_scVI"] = vae.get_normalized_expression(library_size=10e4)

# run PcA then generate UMAp plots
sc.tl.pca(combined_adata)
sc.pl.pca_variance_ratio(combined_adata,log=True)
# use scvI latent space for UMAp generation
sc.pp.neighbors(combined_adata,use_rep="X_scVI")
sc.tl.umap(combined_adata,min_dist=0.3)

#result_folder = "/storage/data/KAI/Pan_Inflammation/Joint/Joint_py/scvi/result/"

sc.pl.umap(combined_adata,color=["GEO","Technology","Group","Sample_Site"], ncols=2, frameon=False, save = "Joint_scvi_batch.png")

#写出
combined_adata.write("/storage/data/KAI/Pan_Inflammation/Joint/Joint_py/scvi/" + "Joint_scvi.h5ad")

combined_adata = sc.read_h5ad("/storage/data/KAI/Pan_Inflammation/Joint/Joint_py/scvi/" + "Joint_scvi.h5ad")


# full gene ---------------------------
combined_adata.raw.var.index = combined_adata.raw.var.index.astype(str)  # 确保索引为字符串
combined_adata.raw.var = combined_adata.raw.var.astype(str)  # 确保所有列的类型为字符串

rubb_var = [
'mt-0', 'ribo-0', 'hb-0', 'n_cells_by_counts-0', 'mean_counts-0', 'log1p_mean_counts-0', 'pct_dropout_by_counts-0', 'total_counts-0', 'log1p_total_counts-0', 'n_cells-0', 'gene_ids-1', 'mt-1', 'ribo-1', 'hb-1', 'n_cells_by_counts-1', 'mean_counts-1', 'log1p_mean_counts-1', 'pct_dropout_by_counts-1', 'total_counts-1', 'log1p_total_counts-1', 'n_cells-1', 'mt-2', 'ribo-2', 'hb-2', 'n_cells_by_counts-2', 'mean_counts-2', 'log1p_mean_counts-2', 'pct_dropout_by_counts-2', 'total_counts-2', 'log1p_total_counts-2', 'n_cells-2', 'Gene_ids-2', 'mt-3', 'ribo-3', 'hb-3', 'n_cells_by_counts-3', 'mean_counts-3', 'log1p_mean_counts-3', 'pct_dropout_by_counts-3', 'total_counts-3', 'log1p_total_counts-3', 'n_cells-3', 'gene_ids-4', 'mt-4', 'ribo-4', 'hb-4', 'n_cells_by_counts-4', 'mean_counts-4', 'log1p_mean_counts-4', 'pct_dropout_by_counts-4', 'total_counts-4', 'log1p_total_counts-4', 'n_cells-4'
]
combined_adata.var.drop(columns=rubb_var, inplace=True)

combined_adata.obs = combined_adata.obs.astype({'Age': 'str'})
combined_adata.obs = combined_adata.obs.astype({'Donor': 'str'})
combined_adata.obs = combined_adata.obs.astype({'Sample_title': 'str'})
combined_adata.write("/storage/data/KAI/Pan_Inflammation/new_inflammation/Joint/data/Joint_all_gene.h5ad")










# cell annotation ------------------------------

sc.tl.leiden(combined_adata, resolution=0.8)
fig, ax = plt.subplots(figsize=(10, 8), dpi=300)
sc.pl.umap(combined_adata, color=["leiden"],
           legend_loc='on data',
           ax=ax,
           save="_joint_leiden.png"
           )
#Epi_marker = ["EPCAM","KRT8","KRT18", "KRT19", "KRT17"]

Fib_marker = ["DCN","LUM","COL1A2", "COL3A1", "COL6A1", "PDPN"]

Mural_marker = ["PDGFRB", "RGS5", "KCNJ8"]

Endo_marker = ["VWF","PLVAP","CLDN5"]

Mac_marker = ["AIF1","CSF1R","LYZ", "CD163", "FCGR2A", "CSF1R", "C1QA", "C1QB","CD68"]

Mono_marker = ["S100A9","S100A8"]

Mast_marker = ["TPSAB1", "KIT","GATA2"]

Dendritic_marker = ["CLEC9A", "IRF8", "CD40", "CD80", "CD83", "CCR7"]

#Neutrophil_marker =["G0S2", "CXCL8", "SOD2", "NAMPT", "S100A8"]

TNK_marker =  [ "CD3D", "CD3E", "CD3G", "CD2", "NKG7", "XCL2", "XCL1"]

B_marker = ["MS4A1","CD79A","CD79B", "CD19", "BLNK", "FCRL5", "CD27", "BANK1"]

Plasma_marker = ["SDC1", "MZB1", "XBP1"]

#Erythrocyte_marker = ["HBA1","HBA2","HBB"]

#Neural_marker = ["NRXN1","S100B","NGFR"]


#sc.pl.violin(combined_adata, Epi_marker, groupby="leiden", stripplot=False)
sc.pl.violin(combined_adata, Fib_marker, groupby="leiden", stripplot=False)
sc.pl.violin(combined_adata, Mural_marker, groupby="leiden", stripplot=False)
sc.pl.violin(combined_adata, Endo_marker, groupby="leiden", stripplot=False)
sc.pl.violin(combined_adata, Mac_marker, groupby="leiden", stripplot=False)
sc.pl.violin(combined_adata, Mono_marker, groupby="leiden", stripplot=False)
sc.pl.violin(combined_adata, Mast_marker, groupby="leiden", stripplot=False)
sc.pl.violin(combined_adata, Dendritic_marker, groupby="leiden", stripplot=False)
#sc.pl.violin(combined_adata, Neutrophil_marker, groupby="leiden", stripplot=False)
sc.pl.violin(combined_adata, TNK_marker, groupby="leiden", stripplot=False)
sc.pl.violin(combined_adata, B_marker, groupby="leiden", stripplot=False)
sc.pl.violin(combined_adata, Plasma_marker, groupby="leiden", stripplot=False)
#sc.pl.violin(combined_adata, Erythrocyte_marker, groupby="leiden", stripplot=False)
#sc.pl.violin(combined_adata, Neural_marker, groupby="leiden", stripplot=False)


sc.tl.rank_genes_groups(combined_adata, "leiden", groups=["0"], method="wilcoxon")
result = combined_adata.uns['rank_genes_groups']
cluster_0_genes = pd.DataFrame({
    'gene': result['names']['0'],
    'logfoldchanges': result['logfoldchanges']['0'],
    'pvals': result['pvals']['0'],
    'pvals_adj': result['pvals_adj']['0']
})
cluster_0_genes['gene'].head(30).to_list()


sc.tl.rank_genes_groups(combined_adata, "leiden", groups=["9"], method="wilcoxon")
result = combined_adata.uns['rank_genes_groups']
cluster_9_genes = pd.DataFrame({
    'gene': result['names']['9'],
    'logfoldchanges': result['logfoldchanges']['9'],
    'pvals': result['pvals']['9'],
    'pvals_adj': result['pvals_adj']['9']
})
cluster_9_genes['gene'].head(30).to_list()

sc.tl.rank_genes_groups(combined_adata, "leiden", groups=["13"], method="wilcoxon")
result = combined_adata.uns['rank_genes_groups']
cluster_13_genes = pd.DataFrame({
    'gene': result['names']['13'],
    'logfoldchanges': result['logfoldchanges']['13'],
    'pvals': result['pvals']['13'],
    'pvals_adj': result['pvals_adj']['13']
})
cluster_13_genes['gene'].head(30).to_list()

sc.tl.rank_genes_groups(combined_adata, "leiden", groups=["18"], method="wilcoxon")
result = combined_adata.uns['rank_genes_groups']
cluster_18_genes = pd.DataFrame({
    'gene': result['names']['18'],
    'logfoldchanges': result['logfoldchanges']['18'],
    'pvals': result['pvals']['18'],
    'pvals_adj': result['pvals_adj']['18']
})
cluster_18_genes['gene'].head(30).to_list()



cell_annotation_level1 = {
    "0": "Fibroblast",
    "1": "Fibroblast",
    "2": "Endothelial cell",
    "3": "T/NK cell",
    "4": "Macrophage",
    "5": "Fibroblast",
    "6": "Fibroblast",
    "7": "Fibroblast",
    "8": "Fibroblast",
    "9": "Macrophage",
    "10": "Fibroblast",
    "11": "Mural cell",
    "12": "Fibroblast",
    "13": "Plasma cell",
    "14": "Monocyte",
    "15": "Fibroblast",
    "16": "B cell",
    "17": "Mast cell",
    "18": "Monocyte",
    "19": "Dendritic cell",
    "20": "Endothelial cell",
    "21": "Macrophage"

}
leiden_annotation_level1 = pd.Series(cell_annotation_level1)
combined_adata.obs['cell_type_level1'] = combined_adata.obs['leiden'].map(leiden_annotation_level1)

cell_color = {
    "Epithelial cell": "#e43030",
    "Cholangiocyte": "#c67e82",
    "Hepatocyte": "#ddcccc",
    "Fibroblast": "#ff8831",
    "Mural cell": "#20452e",
    "Endothelial cell": "#704ba3",
    "Macrophage": "#67a8cd",
    "Monocyte": "#8184e2",
    "Mast cell": "#ff9d9f",
    "Dendritic cell": "#cf9f88",
    "Neutrophil": "#ffc17f",
    "T/NK cell": "#50aa4b",
    "B cell": "#b3e19b",
    "Plasma cell": "#ab3181",
    "Erythrocyte": "#a1acbd",
    "Neural cell": "#00afba"
}
# combined_adata.write("/storage/data/KAI/Pan_Inflammation/Joint/Joint_py/scvi/" + "Joint_scvi_cell_level1.h5ad")

combined_adata = sc.read_h5ad("/storage/data/KAI/Pan_Inflammation/Joint/Joint_py/scvi/" + "Joint_scvi_cell_level1.h5ad")

sc.pl.umap(combined_adata, color="cell_type_level1",
           legend_loc="on data",
           palette=cell_color,
           add_outline=False,
           legend_fontsize=8,
           legend_fontoutline=0.1,
           legend_fontweight="normal",
           title="Joint_scvi_cell_annotation",
           frameon=False,
           save="_Joint_scvi_cell_level1_clear.pdf"
           )
sc.pl.umap(combined_adata, color="cell_type_level1",
           # legend_loc="on data",
           palette=cell_color,
           add_outline=False,
           legend_fontsize=8,
           legend_fontoutline=0.1,
           legend_fontweight="normal",
           title="Joint_scvi_cell_annotation",
           frameon=False,
           save="_Joint_scvi_cell_level1_legend.pdf"
           )

# scanvi-----------------------

lvae = scvi.model.SCANVI.from_scvi_model(
    vae,
    adata=combined_adata,
    labels_key="cell_type_level1",
    unlabeled_category="Unknown",
)
lvae.train()
model_dir2 = os.path.join("/storage/data/KAI/Pan_Inflammation/Joint/Joint_py/scanvi/model/Joint_allcells/",
                          "Joint_allcells_scanvi_model")
lvae.save(model_dir2, overwrite=True, prefix="Joint_allcells_scanvi_")

lvae = scvi.model.SCANVI.load(
    "/storage/data/KAI/Pan_Inflammation/Joint/Joint_py/scanvi/model/Joint_allcells/Joint_allcells_scanvi_model/",
    adata=combined_adata, prefix="Joint_allcells_scanvi_")

combined_adata.obsm["X_scANVI"] = lvae.get_latent_representation()

# run PcA then generate UMAp plots
sc.tl.pca(combined_adata)
sc.pl.pca_variance_ratio(combined_adata, log=True)
# use scANVI latent space for UMAp generation
sc.pp.neighbors(combined_adata, use_rep="X_scANVI")
sc.tl.umap(combined_adata, min_dist=0.3)

fig, ax = plt.subplots(1, 1, figsize=(5, 4))
sc.pl.umap(combined_adata, color="cell_type_level1",
           legend_loc="right margin",
           palette=cell_color,
           add_outline=False,
           legend_fontsize=8,
           legend_fontoutline=0.1,
           legend_fontweight="normal",
           title="Joint_scanvi_cell_annotation",
           frameon=False,
           show=False,
           ax=ax
           # save = "_Joint_scanvi_cell_level1_legend.pdf"
           )
# fig = plt.gcf()
for collection in ax.collections:
    collection.set_rasterized(True)
# for ax in fig.get_axes():

plt.savefig("/storage/data/KAI/Pan_Inflammation/Infla_py/figures/umap_Joint_scanvi_cell_level1_legend.pdf",
            format="pdf", dpi=600)
plt.close()

fig, ax = plt.subplots(1, 1, figsize=(5, 4))
sc.pl.umap(combined_adata, color="cell_type_level1",
           legend_loc="on data",
           palette=cell_color,
           add_outline=False,
           legend_fontsize=8,
           legend_fontoutline=0.1,
           legend_fontweight="normal",
           title="Joint_scanvi_cell_annotation",
           frameon=False,
           show=False,
           ax=ax
           # save = "_joint_scanvi_cell_level1_clear.pdf"
           )
# fig = plt.gcf()
for collection in ax.collections:
    collection.set_rasterized(True)
# for ax in fig.get_axes():

plt.savefig("/storage/data/KAI/Pan_Inflammation/Infla_py/figures/umap_Joint_scanvi_cell_level1_clear.pdf", format="pdf",
            dpi=600)
plt.close()

combined_adata.write("/storage/data/KAI/Pan_Inflammation/Joint/Joint_py/scanvi/Joint_scanvi.h5ad")

sc.pl.umap(combined_adata, color=["GEO"], ncols=2, frameon=False, save="Joint_scanvi_batch.png")


violin_markers = [
  "DCN", "LUM","COL3A1",
  "PDGFRB", "RGS5","ACTA2",
  "VWF","PLVAP","CD34",
  "C1QA","AIF1", "CD68",
  "S100A8","S100A9",
  "IRF8", "CLEC9A","IDO1",
  "KIT","GATA2","HPGD",
  "CD3D", "CD3E", "NKG7",
  "MS4A1","BANK1","CD79A",
  "MZB1", "XBP1"
                  ]

# TO R
combined_adata = sc.read_h5ad("/storage/data/KAI/Pan_Inflammation/Joint/Joint_py/scanvi/Joint_scanvi.h5ad")
combined_adata.X = combined_adata.layers['counts']
combined_adata.write("/storage/data/KAI/Pan_Inflammation/Joint/Joint_py/scanvi/Joint_scanvi_toR.h5ad")







### Liver Inflammation scRNA-seq
adata_list = []

# Chronic hepatitis B
preprecess_path = "/storage/data/KAI/Pan_Inflammation/Liver/Viral hepatitis/Chronic hepatitis B/1.GSE148881/GSE148881_py/GSE148881_preprocess.h5ad"
adata = sc.read_h5ad(preprecess_path)
unique_genes = pd.Series(adata.var_names).drop_duplicates(keep="first")
adata = adata[adata.obs.index, unique_genes.index]
adata.var_names = unique_genes.values
print(f'GSE148881 dim: {adata.shape}')
adata_list.append(adata)

preprecess_path = "/storage/data/KAI/Pan_Inflammation/Liver/Viral hepatitis/Chronic hepatitis B/2.GSE234241/GSE234241_py/GSE234241_preprocess.h5ad"
adata = sc.read_h5ad(preprecess_path)
unique_genes = pd.Series(adata.var_names).drop_duplicates(keep="first")
adata = adata[adata.obs.index, unique_genes.index]
adata.var_names = unique_genes.values
print(f'GSE234241 dim: {adata.shape}')
adata_list.append(adata)

preprecess_path = "/storage/data/KAI/Pan_Inflammation/Liver/Viral hepatitis/Chronic hepatitis B/5.GSE186343/GSE186343_py/GSE186343_preprocess.h5ad"
adata = sc.read_h5ad(preprecess_path)
unique_genes = pd.Series(adata.var_names).drop_duplicates(keep="first")
adata = adata[adata.obs.index, unique_genes.index]
adata.var_names = unique_genes.values
print(f'GSE186343 dim: {adata.shape}')
adata_list.append(adata)



# Alcohol-associated hepatitis
preprecess_path = "/storage/data/KAI/Pan_Inflammation/Liver/Non-viral hepatitis/Alcohol-associated hepatitis/2.GSE255772/GSE255772_py/GSE255772_preprocess.h5ad"
adata = sc.read_h5ad(preprecess_path)
unique_genes = pd.Series(adata.var_names).drop_duplicates(keep="first")
adata = adata[adata.obs.index, unique_genes.index]
adata.var_names = unique_genes.values
print(f'GSE255772 dim: {adata.shape}')
adata_list.append(adata)



# Primary biliary cholangitis
preprecess_path = "/storage/data/KAI/Pan_Inflammation/Liver/Non-viral hepatitis/Primary biliary cholangitis/1.GSE247128/GSE247128_py/GSE247128_preprocess.h5ad"
adata = sc.read_h5ad(preprecess_path)
unique_genes = pd.Series(adata.var_names).drop_duplicates(keep="first")
adata = adata[adata.obs.index, unique_genes.index]
adata.var_names = unique_genes.values
print(f'GSE247128 dim: {adata.shape}')
adata_list.append(adata)


# 检查交集基因
len(adata_list)

var_list = []
for i in range(0,len(adata_list)):
    var_gse = list(adata_list[i].var_names)
    var_list.append(var_gse)


inte_adata = []
for i in range(1,len(var_list)+1):
    intesect = set.intersection(*[set(v) for v in var_list[:i]])
    inte_adata.append(intesect)

for i, inter in enumerate(inte_adata):
    print(f"交集（前 {i+1} 个集合）: {len(inter)}")




for i, ad in enumerate(adata_list):
    print(f"AnnData {i} - Duplicate obs columns: {ad.obs.columns.duplicated().sum()}")
    print(f"AnnData {i} - Duplicate var columns: {ad.var.columns.duplicated().sum()}")


for ad in adata_list:
    ad.obs_names_make_unique()


all_obs_names=[]
for ad in adata_list:
    all_obs_names.extend(ad.obs.index.to_list())




# combined_adata = adata_list[0].concatenate(adata_list[1:], batch_key=None, batch_categories=["1","2","3","4","5"], join='inner')
#
# combined1 = adata.concatenate(adata_list[1], join='inner', batch_key=None, batch_categories=None)
# combined2 = adata.concatenate(adata_list[2], join='inner', batch_key=None, batch_categories=None)
# combined3 = adata.concatenate(adata_list[3], join='inner', batch_key=None, batch_categories=None)
# combined4 = adata.concatenate(adata_list[4], join='inner', batch_key=None, batch_categories=None)


# Integration
combined_adata = adata_list[0].concatenate(adata_list[1:], join='outer', batch_key=None, batch_categories=None) # 整合后 (165820, 15730)

temp = combined_adata

# bool转0-1
#temp.obs = temp.obs.astype({'predicted_doublet': 'int'})

# 转换object类型
temp.obs = temp.obs.astype({'Age': 'str'}) #数据中不能含有缺失值(NaN)

# bool转0-1
#temp.var = temp.var.astype({'mt': 'int'})
#temp.var = temp.var.astype({'ribo': 'int'})
#temp.var = temp.var.astype({'hb': 'int'})
#temp.var = temp.var.astype({'highly_variable': 'int'})

# 转换object类型
temp.var = temp.var.astype({'Gene_ids-0': 'str'})
temp.var = temp.var.astype({'gene_ids-1': 'str'})
temp.var = temp.var.astype({'gene_ids-3': 'str'})


# 删除var中多余的属性

rubb_var = ['feature_types-0', 'genome-0', 'n_cells_by_counts-0', 'mean_counts-0', 'log1p_mean_counts-0', 'pct_dropout_by_counts-0', 'total_counts-0', 'log1p_total_counts-0', 'n_cells-0', 'gene_ids-1', 'n_cells_by_counts-1', 'mean_counts-1', 'log1p_mean_counts-1', 'pct_dropout_by_counts-1', 'total_counts-1', 'log1p_total_counts-1', 'n_cells-1', 'gene_ids-2', 'feature_types-2', 'genome-2', 'n_cells_by_counts-2', 'mean_counts-2', 'log1p_mean_counts-2', 'pct_dropout_by_counts-2', 'total_counts-2', 'log1p_total_counts-2', 'n_cells-2', 'gene_ids-3', 'feature_types-3', 'genome-3', 'n_cells_by_counts-3', 'mean_counts-3', 'log1p_mean_counts-3', 'pct_dropout_by_counts-3', 'total_counts-3', 'log1p_total_counts-3', 'n_cells-3', 'n_cells_by_counts-4', 'mean_counts-4', 'log1p_mean_counts-4', 'pct_dropout_by_counts-4', 'total_counts-4', 'log1p_total_counts-4', 'n_cells-4', 'Gene_ids-4']
temp.var.drop(columns=rubb_var, inplace=True)

temp.var.rename(columns={'gene_ids-0':'Gene_ids'}, inplace=True)


# 写出整合后文件(未进行scvi)
#temp.write("/storage/data/KAI/Pan_Inflammation/Liver/Liver_py/scvi/" + "Liver_Inte.h5ad")


# scvi-----------------------
combined_adata = temp

#combined_adata = sc.read_h5ad("/storage/data/KAI/Pan_Inflammation/Liver/Liver_py/scvi/" + "Liver_Inte.h5ad")

combined_adata.raw = combined_adata  # keep full dimension safe


sc.pp.highly_variable_genes(combined_adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

scvi.model.SCVI.setup_anndata(combined_adata, layer="counts", categorical_covariate_keys=["Sample_geo_accession", "Sample_Type", "Technology", "Sample_Site"])
vae = scvi.model.SCVI(combined_adata, n_latent=20, n_layers=2, dropout_rate=0.1)
vae.train(early_stopping=True, batch_size=256)
# 保存scvi模型
model_dir = os.path.join("/storage/data/KAI/Pan_Inflammation/Liver/Liver_py/scvi/model/liver_allcells/", "liver_allcells_scvi_model")
vae.save(model_dir, overwrite=True, prefix="liver_allcells_scvi_")

vae = scvi.model.SCVI.load("/storage/data/KAI/Pan_Inflammation/Liver/Liver_py/scvi/model/liver_allcells/liver_allcells_scvi_model/", adata=combined_adata, prefix="liver_allcells_scvi_")

combined_adata.obsm["X_scVI"] = vae.get_latent_representation()
combined_adata.layers["X_normalized_scVI"] = vae.get_normalized_expression(library_size=10e4)

# run PcA then generate UMAp plots
sc.tl.pca(combined_adata)
sc.pl.pca_variance_ratio(combined_adata,log=True)
# use scvI latent space for UMAp generation
sc.pp.neighbors(combined_adata,use_rep="X_scVI")
sc.tl.umap(combined_adata,min_dist=0.3)

#result_folder = "/storage/data/KAI/Pan_Inflammation/Liver/Liver_py/scvi/result/"

sc.pl.umap(combined_adata,color=["GEO","Technology","Group","Sample_Site"], ncols=2, frameon=False, save = "Liver_scvi_batch.png")

#写出
combined_adata.write("/storage/data/KAI/Pan_Inflammation/Liver/Liver_py/scvi/" + "Liver_scvi.h5ad")

#combined_adata = sc.read_h5ad("/storage/data/KAI/Pan_Inflammation/Liver/Liver_py/scvi/" + "Liver_scvi.h5ad")



# full gene ---------------------------
combined_adata.raw.var.index = combined_adata.raw.var.index.astype(str)  # 确保索引为字符串
combined_adata.raw.var = combined_adata.raw.var.astype(str)  # 确保所有列的类型为字符串

rubb_var = [
'feature_types-0', 'genome-0', 'mt-0', 'ribo-0', 'hb-0', 'n_cells_by_counts-0', 'mean_counts-0', 'log1p_mean_counts-0', 'pct_dropout_by_counts-0', 'total_counts-0', 'log1p_total_counts-0', 'n_cells-0', 'gene_ids-1', 'mt-1', 'ribo-1', 'hb-1', 'n_cells_by_counts-1', 'mean_counts-1', 'log1p_mean_counts-1', 'pct_dropout_by_counts-1', 'total_counts-1', 'log1p_total_counts-1', 'n_cells-1', 'gene_ids-2', 'feature_types-2', 'genome-2', 'mt-2', 'ribo-2', 'hb-2', 'n_cells_by_counts-2', 'mean_counts-2', 'log1p_mean_counts-2', 'pct_dropout_by_counts-2', 'total_counts-2', 'log1p_total_counts-2', 'n_cells-2', 'gene_ids-3', 'feature_types-3', 'genome-3', 'mt-3', 'ribo-3', 'hb-3', 'n_cells_by_counts-3', 'mean_counts-3', 'log1p_mean_counts-3', 'pct_dropout_by_counts-3', 'total_counts-3', 'log1p_total_counts-3', 'n_cells-3', 'mt-4', 'ribo-4', 'hb-4', 'n_cells_by_counts-4', 'mean_counts-4', 'log1p_mean_counts-4', 'pct_dropout_by_counts-4', 'total_counts-4', 'log1p_total_counts-4', 'n_cells-4', 'Gene_ids-4'
]
combined_adata.var.drop(columns=rubb_var, inplace=True)

combined_adata.obs = combined_adata.obs.astype({'Age': 'str'})
combined_adata.obs = combined_adata.obs.astype({'Donor': 'str'})
combined_adata.obs = combined_adata.obs.astype({'Sample_title': 'str'})
combined_adata.write("/storage/data/KAI/Pan_Inflammation/new_inflammation/Liver/data/Liver_all_gene.h5ad")








# 细胞注释 ---------------------------------

sc.tl.leiden(combined_adata, resolution=0.8)
fig, ax = plt.subplots(figsize=(10, 8), dpi=300)
sc.pl.umap(combined_adata, color=["leiden"],
           legend_loc='on data',
           ax=ax,
           save="_liver_leiden.png"
           )



Cholangiocyte_marker = ["EPCAM","KRT19"]

Hepatocyte_marker = ['ALB','TF','TTR','HNF4A','CYP2A6']

Fib_marker = ["DCN","LUM","COL1A2", "COL3A1", "COL6A1", "PDPN"]

Mural_marker = ["PDGFRB", "RGS5", "KCNJ8"]

Endo_marker = ["VWF","PLVAP","CLDN5"]

Mac_marker = ["AIF1","CSF1R","LYZ", "CD163", "FCGR2A", "CSF1R", "C1QA", "C1QB"]

Mono_marker = ["S100A9","S100A8"]

Mast_marker = ["KIT","TPSAB1","GATA2","HPGD"]

Dendritic_marker = ["CLEC9A", "IRF8", "CD40", "CD80", "CD83", "CCR7"]

Neutrophil_marker =["G0S2", "SOD2", "NAMPT", "S100A8"]

TNK_marker =  ["TRAC", "CD3D", "CD3E", "CD3G", "CD2", "NKG7", "TRDC", "XCL2", "XCL1"]

B_marker = ["MS4A1","CD79A","CD79B", "CD19", "BLNK", "FCRL5", "IGHG1", "CD27", "BANK1"]

Plasma_marker = ["IGHA1", "IGHG1","IGKC", "IGLC2", "IGHM"]

Erythrocyte_marker = ["HBA1","HBA2","HBB"]

Neural_marker = ["S100B","NGFR", "RBFOX3","MAPT","MAP2","ENO2","UCHL1"]


sc.pl.violin(combined_adata, Cholangiocyte_marker, groupby="leiden", stripplot=False)
sc.pl.violin(combined_adata, Hepatocyte_marker, groupby="leiden", stripplot=False)
sc.pl.violin(combined_adata, Fib_marker, groupby="leiden", stripplot=False)
sc.pl.violin(combined_adata, Mural_marker, groupby="leiden", stripplot=False)
sc.pl.violin(combined_adata, Endo_marker, groupby="leiden", stripplot=False)
sc.pl.violin(combined_adata, Mac_marker, groupby="leiden", stripplot=False)
sc.pl.violin(combined_adata, Mono_marker, groupby="leiden", stripplot=False)
sc.pl.violin(combined_adata, Mast_marker, groupby="leiden", stripplot=False)
sc.pl.violin(combined_adata, Dendritic_marker, groupby="leiden", stripplot=False)
sc.pl.violin(combined_adata, Neutrophil_marker, groupby="leiden", stripplot=False)
sc.pl.violin(combined_adata, TNK_marker, groupby="leiden", stripplot=False)
sc.pl.violin(combined_adata, B_marker, groupby="leiden", stripplot=False)
sc.pl.violin(combined_adata, Plasma_marker, groupby="leiden", stripplot=False)
sc.pl.violin(combined_adata, Erythrocyte_marker, groupby="leiden", stripplot=False)
#sc.pl.violin(combined_adata, Neural_marker, groupby="leiden", stripplot=False)


sc.tl.rank_genes_groups(combined_adata, "leiden", groups=["18"], method="wilcoxon")
result = combined_adata.uns['rank_genes_groups']
cluster_18_genes = pd.DataFrame({
    'gene': result['names']['18'],
    'logfoldchanges': result['logfoldchanges']['18'],
    'pvals': result['pvals']['18'],
    'pvals_adj': result['pvals_adj']['18']
})
cluster_18_genes['gene'].head(30).to_list()


sc.tl.rank_genes_groups(combined_adata, "leiden", groups=["21"], method="wilcoxon")
result = combined_adata.uns['rank_genes_groups']
cluster_21_genes = pd.DataFrame({
    'gene': result['names']['21'],
    'logfoldchanges': result['logfoldchanges']['21'],
    'pvals': result['pvals']['21'],
    'pvals_adj': result['pvals_adj']['21']
})
cluster_21_genes['gene'].head(50).to_list()


sc.tl.rank_genes_groups(combined_adata, "leiden", groups=["11"], method="wilcoxon")
result = combined_adata.uns['rank_genes_groups']
cluster_11_genes = pd.DataFrame({
    'gene': result['names']['11'],
    'logfoldchanges': result['logfoldchanges']['11'],
    'pvals': result['pvals']['11'],
    'pvals_adj': result['pvals_adj']['11']
})
cluster_11_genes['gene'].head(20).to_list()


sc.tl.rank_genes_groups(combined_adata, "leiden", groups=["28"], method="wilcoxon")
result = combined_adata.uns['rank_genes_groups']
cluster_28_genes = pd.DataFrame({
    'gene': result['names']['28'],
    'logfoldchanges': result['logfoldchanges']['28'],
    'pvals': result['pvals']['28'],
    'pvals_adj': result['pvals_adj']['28']
})
cluster_28_genes['gene'].head(50).to_list()




sc.pl.umap(combined_adata, color=["DCN"],
           legend_loc='on data'
           )


cell_annotation_level1 = {
"0"	:"T/NK cell",
"1" :"T/NK cell",
"2"	:"Macrophage",
"3"	:"Macrophage",
"4"	:"T/NK cell",
"5"	:"T/NK cell",
"6"	:"T/NK cell",
"7"	:"T/NK cell",
"8"	:"Macrophage",
"9"	:"B cell",
"10":"Endothelial cell",
"11":"Neutrophil",
"12":"T/NK cell",
"13":"Macrophage",
"14":"Macrophage",
"15":"T/NK cell",
"16":"T/NK cell",
"17":"Hepatocyte",
"18":"Neutrophil",
"19":"Cholangiocyte",
"20":"Fibroblast",
"21":"Monocyte",
"22":"T/NK cell",
"23":"Plasma cell",
"24":"Hepatocyte",
"25":"Erythrocyte",
"26":"Mast cell",
"27":"Dendritic cell",
"28":"Dendritic cell",
"29":"T/NK cell",

}
leiden_annotation_level1 = pd.Series(cell_annotation_level1)
combined_adata.obs['cell_type_level1'] = combined_adata.obs['leiden'].map(leiden_annotation_level1)

cell_color = {
    "Epithelial cell":"#e43030",
    "Cholangiocyte":"#c67e82",
    "Hepatocyte":"#ddcccc",
    "Fibroblast":"#ff8831",
    "Mural cell":"#20452e",
    "Endothelial cell":"#704ba3",
    "Macrophage":"#67a8cd",
    "Monocyte":"#8184e2",
    "Mast cell":"#ff9d9f",
    "Dendritic cell":"#cf9f88",
    "Neutrophil":"#ffc17f",
    "T/NK cell":"#50aa4b",
    "B cell":"#b3e19b",
    "Plasma cell":"#ab3181",
    "Erythrocyte":"#a1acbd",
    "Neural cell":"#00afba"
}
#combined_adata.write("/storage/data/KAI/Pan_Inflammation/Liver/Liver_py/scvi/" + "Liver_scvi_cell_level1.h5ad")

combined_adata = sc.read_h5ad("/storage/data/KAI/Pan_Inflammation/Liver/Liver_py/scvi/" + "Liver_scvi_cell_level1.h5ad")

sc.pl.umap(combined_adata, color="cell_type_level1",
           legend_loc="on data",
           palette=cell_color,
           add_outline=False,
           legend_fontsize=8,
           legend_fontoutline=0.1,
           legend_fontweight="normal",
           title="Liver_scvi_cell_annotation",
           frameon=False,
           save = "_Liver_scvi_cell_level1_clear.pdf"
           )
sc.pl.umap(combined_adata, color="cell_type_level1",
           #legend_loc="on data",
           palette=cell_color,
           add_outline=False,
           legend_fontsize=8,
           legend_fontoutline=0.1,
           legend_fontweight="normal",
           title="Liver_scvi_cell_annotation",
           frameon=False,
           save = "_Liver_scvi_cell_level1_legend.pdf"
           )


# scanvi-----------------------

lvae = scvi.model.SCANVI.from_scvi_model(
    vae,
    adata=combined_adata,
    labels_key="cell_type_level1",
    unlabeled_category="Unknown",
)
lvae.train()
model_dir2 = os.path.join("/storage/data/KAI/Pan_Inflammation/Liver/Liver_py/scanvi/model/liver_allcells/", "liver_allcells_scanvi_model")
lvae.save(model_dir2, overwrite=True, prefix="liver_allcells_scanvi_")

lvae = scvi.model.SCANVI.load("/storage/data/KAI/Pan_Inflammation/Liver/Liver_py/scanvi/model/liver_allcells/liver_allcells_scanvi_model/", adata=combined_adata, prefix="liver_allcells_scanvi_")


combined_adata.obsm["X_scANVI"] = lvae.get_latent_representation()

# run PcA then generate UMAp plots
sc.tl.pca(combined_adata)
sc.pl.pca_variance_ratio(combined_adata,log=True)
# use scANVI latent space for UMAp generation
sc.pp.neighbors(combined_adata,use_rep="X_scANVI")
sc.tl.umap(combined_adata,min_dist=0.3)


fig, ax=plt.subplots(1,1,figsize=(5,4))
sc.pl.umap(combined_adata, color="cell_type_level1",
           legend_loc="right margin",
           palette=cell_color,
           add_outline=False,
           legend_fontsize=8,
           legend_fontoutline=0.1,
           legend_fontweight="normal",
           title="Liver_scanvi_cell_annotation",
           frameon=False,
           show=False,
           ax=ax
           #save = "_Liver_scanvi_cell_level1_legend.pdf"
           )
# fig = plt.gcf()
for collection in ax.collections:
    collection.set_rasterized(True)
# for ax in fig.get_axes():

plt.savefig("/storage/data/KAI/Pan_Inflammation/Infla_py/figures/umap_Liver_scanvi_cell_level1_legend.pdf", format="pdf", dpi=600)
plt.close()

fig, ax=plt.subplots(1,1,figsize=(5,4))
sc.pl.umap(combined_adata, color="cell_type_level1",
           legend_loc="on data",
           palette=cell_color,
           add_outline=False,
           legend_fontsize=8,
           legend_fontoutline=0.1,
           legend_fontweight="normal",
           title="Liver_scanvi_cell_annotation",
           frameon=False,
           show=False,
           ax=ax
           #save = "_Liver_scanvi_cell_level1_clear.pdf"
           )
# fig = plt.gcf()
for collection in ax.collections:
    collection.set_rasterized(True)
# for ax in fig.get_axes():

plt.savefig("/storage/data/KAI/Pan_Inflammation/Infla_py/figures/umap_Liver_scanvi_cell_level1_clear.pdf", format="pdf", dpi=600)
plt.close()


combined_adata.write("/storage/data/KAI/Pan_Inflammation/Liver/Liver_py/scanvi/Liver_scanvi.h5ad")
sc.pl.umap(combined_adata,color=["GEO"], ncols=2, frameon=False, save = "liver_scanvi_batch.png")



combined_adata = sc.read_h5ad("/storage/data/KAI/Pan_Inflammation/Liver/Liver_py/scanvi/Liver_scanvi.h5ad")

# combined_adata.obs['cell_type_level1'] = combined_adata.obs['cell_type_level1'].astype('category')
# sc.tl.rank_genes_groups(combined_adata, 'cell_type_level1', groups=['Monocyte'], reference='rest', method='wilcoxon')
#
# result = combined_adata.uns['rank_genes_groups']
# cluster_mono_genes = pd.DataFrame({
#     'gene': result['names']['Monocyte'],
#     'logfoldchanges': result['logfoldchanges']['Monocyte'],
#     'pvals': result['pvals']['Monocyte'],
#     'pvals_adj': result['pvals_adj']['Monocyte']
# })
# cluster_mono_genes['gene'].head(50).to_list()


violin_markers = [
  "EPCAM", "KRT7", "KRT18",
  'ALB','APOC3','CYP3A4',
  "DCN", "LUM","COL3A1",
  "VWF","PLVAP","CLDN5",
  "C1QA","CD68","CD163",
  "S100A9",
  "IRF8", "PLD4", "LILRA4",
  "KIT","GATA2","HPGD",
  "ALPL","CXCR1", "CXCR2",
  "TRAC","CD3D", "CD3E",
  "MS4A1","BANK1","CD79A",
  "IGHA1", "IGHG1","IGKC",
  "HBA1","HBA2","ALAS2"
                  ]

# TO R
combined_adata = sc.read_h5ad("/storage/data/KAI/Pan_Inflammation/Liver/Liver_py/scanvi/Liver_scanvi.h5ad")
combined_adata.X = combined_adata.layers['counts']
combined_adata.write("/storage/data/KAI/Pan_Inflammation/Liver/Liver_py/scanvi/Liver_scanvi_toR.h5ad")





### Lung Inflammation scRNA-seq
adata_list = []

# Chronic obstructive pulmonary disease
preprecess_path = "/storage/data/KAI/Pan_Inflammation/Lung/Chronic obstructive pulmonary disease/1.GSE173896/GSE173896_py/GSE173896_preprocess.h5ad"
adata = sc.read_h5ad(preprecess_path)
unique_genes = pd.Series(adata.var_names).drop_duplicates(keep="first")
adata = adata[adata.obs.index, unique_genes.index]
adata.var_names = unique_genes.values
print(f'GSE173896 dim: {adata.shape}')
adata_list.append(adata)

preprecess_path = "/storage/data/KAI/Pan_Inflammation/Lung/Chronic obstructive pulmonary disease/2.GSE136831/GSE136831_py/GSE136831_preprocess.h5ad"
adata = sc.read_h5ad(preprecess_path)
unique_genes = pd.Series(adata.var_names).drop_duplicates(keep="first")
adata = adata[adata.obs.index, unique_genes.index]
adata.var_names = unique_genes.values
print(f'GSE136831 dim: {adata.shape}')
adata_list.append(adata)


preprecess_path = "/storage/data/KAI/Pan_Inflammation/Lung/Chronic obstructive pulmonary disease/6.GSE270667/GSE270667_py/GSE270667_preprocess.h5ad"
adata = sc.read_h5ad(preprecess_path)
unique_genes = pd.Series(adata.var_names).drop_duplicates(keep="first")
adata = adata[adata.obs.index, unique_genes.index]
adata.var_names = unique_genes.values
print(f'GSE270667 dim: {adata.shape}')
adata_list.append(adata)


preprecess_path = "/storage/data/KAI/Pan_Inflammation/Lung/Chronic obstructive pulmonary disease/7.GSE171541/GSE171541_py/GSE171541_preprocess.h5ad"
adata = sc.read_h5ad(preprecess_path)
unique_genes = pd.Series(adata.var_names).drop_duplicates(keep="first")
adata = adata[adata.obs.index, unique_genes.index]
adata.var_names = unique_genes.values
print(f'GSE171541 dim: {adata.shape}')
adata_list.append(adata)



# Hypersensitivity pneumonitis
preprecess_path = "/storage/data/KAI/Pan_Inflammation/Lung/Hypersensitivity pneumonitis/1.GSE122960/GSE122960_py/GSE122960_preprocess.h5ad"
adata = sc.read_h5ad(preprecess_path)
unique_genes = pd.Series(adata.var_names).drop_duplicates(keep="first")
adata = adata[adata.obs.index, unique_genes.index]
adata.var_names = unique_genes.values
print(f'GSE122960 dim: {adata.shape}')
adata_list.append(adata)


preprecess_path = "/storage/data/KAI/Pan_Inflammation/Lung/Hypersensitivity pneumonitis/2.GSE227136/GSE227136_py/GSE227136_preprocess.h5ad"
adata = sc.read_h5ad(preprecess_path)
unique_genes = pd.Series(adata.var_names).drop_duplicates(keep="first")
adata = adata[adata.obs.index, unique_genes.index]
adata.var_names = unique_genes.values
print(f'GSE227136 dim: {adata.shape}')
adata_list.append(adata)



# COVID-19
preprecess_path = "/storage/data/KAI/Pan_Inflammation/Lung/COVID-19/1.GSE171668/GSE171668_py/GSE171668_preprocess.h5ad"
adata = sc.read_h5ad(preprecess_path)
unique_genes = pd.Series(adata.var_names).drop_duplicates(keep="first")
adata = adata[adata.obs.index, unique_genes.index]
adata.var_names = unique_genes.values
print(f'GSE171668 dim: {adata.shape}')
adata_list.append(adata)


preprecess_path = "/storage/data/KAI/Pan_Inflammation/Lung/COVID-19/2.GSE163919/GSE163919_py/GSE163919_preprocess.h5ad"
adata = sc.read_h5ad(preprecess_path)
unique_genes = pd.Series(adata.var_names).drop_duplicates(keep="first")
adata = adata[adata.obs.index, unique_genes.index]
adata.var_names = unique_genes.values
print(f'GSE163919 dim: {adata.shape}')
adata_list.append(adata)


preprecess_path = "/storage/data/KAI/Pan_Inflammation/Lung/COVID-19/3.GSE149878/GSE149878_py/GSE149878_preprocess.h5ad"
adata = sc.read_h5ad(preprecess_path)
unique_genes = pd.Series(adata.var_names).drop_duplicates(keep="first")
adata = adata[adata.obs.index, unique_genes.index]
adata.var_names = unique_genes.values
print(f'GSE149878 dim: {adata.shape}')
adata_list.append(adata)


# 检查交集基因
len(adata_list)

var_list = []
for i in range(0,len(adata_list)):
    var_gse = list(adata_list[i].var_names)
    var_list.append(var_gse)


inte_adata = []
for i in range(1,len(var_list)+1):
    intesect = set.intersection(*[set(v) for v in var_list[:i]])
    inte_adata.append(intesect)

for i, inter in enumerate(inte_adata):
    print(f"交集（前 {i+1} 个集合）: {len(inter)}")


# Integration
combined_adata = adata_list[0].concatenate(adata_list[1:], join='outer', batch_key=None, batch_categories=None) # 整合后 (432357,14753)

temp = combined_adata

# bool转0-1
#temp.obs = temp.obs.astype({'predicted_doublet': 'int'})

# 转换object类型
temp.obs = temp.obs.astype({'Age': 'str'}) #数据中不能含有缺失值(NaN)

# bool转0-1
#temp.var = temp.var.astype({'mt': 'int'})
#temp.var = temp.var.astype({'ribo': 'int'})
#temp.var = temp.var.astype({'hb': 'int'})
#temp.var = temp.var.astype({'highly_variable': 'int'})

# 转换object类型
temp.var = temp.var.astype({'gene_ids-1': 'str'})
temp.var = temp.var.astype({'gene_ids-3': 'str'})


# 删除var中多余的属性

rubb_var = ['n_cells_by_counts-0', 'mean_counts-0', 'log1p_mean_counts-0', 'pct_dropout_by_counts-0', 'total_counts-0', 'log1p_total_counts-0', 'n_cells-0', 'n_cells_by_counts-1', 'mean_counts-1', 'log1p_mean_counts-1', 'pct_dropout_by_counts-1', 'total_counts-1', 'log1p_total_counts-1', 'n_cells-1', 'n_cells_by_counts-2', 'mean_counts-2', 'log1p_mean_counts-2', 'pct_dropout_by_counts-2', 'total_counts-2', 'log1p_total_counts-2', 'n_cells-2', 'Gene_ids-2', 'n_cells_by_counts-3', 'mean_counts-3', 'log1p_mean_counts-3', 'pct_dropout_by_counts-3', 'total_counts-3', 'log1p_total_counts-3', 'n_cells-3', 'gene_ids-4', 'n_cells_by_counts-4', 'mean_counts-4', 'log1p_mean_counts-4', 'pct_dropout_by_counts-4', 'total_counts-4', 'log1p_total_counts-4', 'n_cells-4', 'gene_ids-5', 'n_cells_by_counts-5', 'mean_counts-5', 'log1p_mean_counts-5', 'pct_dropout_by_counts-5', 'total_counts-5', 'log1p_total_counts-5', 'n_cells-5', 'gene_ids-6', 'n_cells_by_counts-6', 'mean_counts-6', 'log1p_mean_counts-6', 'pct_dropout_by_counts-6', 'total_counts-6', 'log1p_total_counts-6', 'n_cells-6', 'feature_types-6', 'genome-6', 'gene_ids-7', 'n_cells_by_counts-7', 'mean_counts-7', 'log1p_mean_counts-7', 'pct_dropout_by_counts-7', 'total_counts-7', 'log1p_total_counts-7', 'n_cells-7', 'feature_types-7', 'genome-7', 'gene_ids-8', 'n_cells_by_counts-8', 'mean_counts-8', 'log1p_mean_counts-8', 'pct_dropout_by_counts-8', 'total_counts-8', 'log1p_total_counts-8', 'n_cells-8', 'feature_types-8', 'genome-8']
temp.var.drop(columns=rubb_var, inplace=True)

temp.var.rename(columns={'gene_ids-0':'Gene_ids'}, inplace=True)


# 写出整合后文件(未进行scvi)
#temp.write("/storage/data/KAI/Pan_Inflammation/Lung/Lung_py/scvi/" + "Lung_Inte.h5ad")


# scvi-----------------------
combined_adata = temp

#combined_adata = sc.read_h5ad("/storage/data/KAI/Pan_Inflammation/Lung/Lung_py/scvi/" + "Lung_Inte.h5ad")

combined_adata.raw = combined_adata  # keep full dimension safe


sc.pp.highly_variable_genes(combined_adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

scvi.model.SCVI.setup_anndata(combined_adata, layer="counts", categorical_covariate_keys=["Sample_geo_accession", "Sample_Type", "Technology", "Sample_Site","GEO"])
vae = scvi.model.SCVI(combined_adata, n_latent=20, n_layers=2, dropout_rate=0.1)
vae.train(max_epochs=50, early_stopping=True, batch_size=256)
# 保存scvi模型
model_dir = os.path.join("/storage/data/KAI/Pan_Inflammation/Lung/Lung_py/scvi/model/Lung_allcells/", "Lung_allcells_scvi_model")
vae.save(model_dir, overwrite=True, prefix="Lung_allcells_scvi_")

vae = scvi.model.SCVI.load("/storage/data/KAI/Pan_Inflammation/Lung/Lung_py/scvi/model/Lung_allcells/Lung_allcells_scvi_model/", adata=combined_adata, prefix="Lung_allcells_scvi_")

combined_adata.obsm["X_scVI"] = vae.get_latent_representation()
combined_adata.layers["X_normalized_scVI"] = vae.get_normalized_expression(library_size=10e4)

# run PcA then generate UMAp plots
sc.tl.pca(combined_adata)
sc.pl.pca_variance_ratio(combined_adata,log=True)
# use scvI latent space for UMAp generation
sc.pp.neighbors(combined_adata,use_rep="X_scVI")
sc.tl.umap(combined_adata,min_dist=0.3)


#sc.pl.umap(combined_adata,color=["GEO","Technology","Group","Sample_Site"], ncols=2, frameon=False, save = "Lung_scvi_batch.png")

fig, ax=plt.subplots(2,2,figsize=(10,8))
axes = ax.flatten()
for i, color in enumerate(["GEO", "Technology", "Group", "Sample_Site"]):
    sc.pl.umap(
        combined_adata,
        color=color,
        legend_loc="right margin",
        add_outline=False,
        legend_fontsize=8,
        legend_fontoutline=0.1,
        legend_fontweight="normal",
        title=f"Lung_scvi_batch - {color}",
        frameon=False,
        show=False,
        ax=axes[i],  # 指定子图
    )
    # fig = plt.gcf()
for ax_i in axes:
    for collection in ax_i.collections:
        collection.set_rasterized(True)

plt.savefig("/storage/data/KAI/Pan_Inflammation/Infla_py/figures/Lung_scvi_batch.pdf", format="pdf", dpi=600)
plt.close()



#写出
combined_adata.write("/storage/data/KAI/Pan_Inflammation/Lung/Lung_py/scvi/" + "Lung_scvi.h5ad")

#combined_adata = sc.read_h5ad("/storage/data/KAI/Pan_Inflammation/Lung/Lung_py/scvi/" + "Lung_scvi.h5ad")





# full gene ---------------------------
combined_adata.raw.var.index = combined_adata.raw.var.index.astype(str)  # 确保索引为字符串
combined_adata.raw.var = combined_adata.raw.var.astype(str)  # 确保所有列的类型为字符串

rubb_var = [
'mt-0', 'ribo-0', 'hb-0', 'n_cells_by_counts-0', 'mean_counts-0', 'log1p_mean_counts-0', 'pct_dropout_by_counts-0', 'total_counts-0', 'log1p_total_counts-0', 'n_cells-0', 'mt-1', 'ribo-1', 'hb-1', 'n_cells_by_counts-1', 'mean_counts-1', 'log1p_mean_counts-1', 'pct_dropout_by_counts-1', 'total_counts-1', 'log1p_total_counts-1', 'n_cells-1', 'mt-2', 'ribo-2', 'hb-2', 'n_cells_by_counts-2', 'mean_counts-2', 'log1p_mean_counts-2', 'pct_dropout_by_counts-2', 'total_counts-2', 'log1p_total_counts-2', 'n_cells-2', 'Gene_ids-2', 'mt-3', 'ribo-3', 'hb-3', 'n_cells_by_counts-3', 'mean_counts-3', 'log1p_mean_counts-3', 'pct_dropout_by_counts-3', 'total_counts-3', 'log1p_total_counts-3', 'n_cells-3', 'gene_ids-4', 'mt-4', 'ribo-4', 'hb-4', 'n_cells_by_counts-4', 'mean_counts-4', 'log1p_mean_counts-4', 'pct_dropout_by_counts-4', 'total_counts-4', 'log1p_total_counts-4', 'n_cells-4', 'gene_ids-5', 'mt-5', 'ribo-5', 'hb-5', 'n_cells_by_counts-5', 'mean_counts-5', 'log1p_mean_counts-5', 'pct_dropout_by_counts-5', 'total_counts-5', 'log1p_total_counts-5', 'n_cells-5', 'gene_ids-6', 'mt-6', 'ribo-6', 'hb-6', 'n_cells_by_counts-6', 'mean_counts-6', 'log1p_mean_counts-6', 'pct_dropout_by_counts-6', 'total_counts-6', 'log1p_total_counts-6', 'n_cells-6', 'feature_types-6', 'genome-6', 'gene_ids-7', 'mt-7', 'ribo-7', 'hb-7', 'n_cells_by_counts-7', 'mean_counts-7', 'log1p_mean_counts-7', 'pct_dropout_by_counts-7', 'total_counts-7', 'log1p_total_counts-7', 'n_cells-7', 'feature_types-7', 'genome-7', 'gene_ids-8', 'mt-8', 'ribo-8', 'hb-8', 'n_cells_by_counts-8', 'mean_counts-8', 'log1p_mean_counts-8', 'pct_dropout_by_counts-8', 'total_counts-8', 'log1p_total_counts-8', 'n_cells-8', 'feature_types-8', 'genome-8'
]
combined_adata.var.drop(columns=rubb_var, inplace=True)

combined_adata.obs = combined_adata.obs.astype({'Age': 'str'})
combined_adata.obs = combined_adata.obs.astype({'Donor': 'str'})
combined_adata.obs = combined_adata.obs.astype({'Sample_title': 'str'})
combined_adata.write("/storage/data/KAI/Pan_Inflammation/new_inflammation/Lung/data/Lung_all_gene.h5ad")





# 细胞注释-------------------------------------
sc.tl.leiden(combined_adata, resolution=0.8)
fig, ax = plt.subplots(figsize=(10, 8), dpi=300)
sc.pl.umap(combined_adata, color=["leiden"],
           legend_loc='on data',
           ax=ax,
           save="_lung_leiden.png"
           )

Epi_marker = ["EPCAM","KRT18","KRT19", "CDH1","CAPS"]

Fib_marker = ["DCN","LUM","COL1A2", "COL3A1", "COL1A1", "THY1"]

#Mural_marker = ["PDGFRB", "RGS5", "CSPG4", "ACTA2","TRPC6"]

Endo_marker = ["PECAM1","VWF","PLVAP","CLDN5","CDH5"]

Mac_marker = ["AIF1","CSF1R","LYZ", "CD163", "FCGR2A", "CSF1R", "C1QA", "C1QB","MARCO"]

Mono_marker = ["S100A9","S100A8"]

Mast_marker = ["TPSAB1", "KIT","CPA3","GATA2"]

Dendritic_marker = ["CLEC9A", "IRF8", "CD1C", "CCL17", "CD83", "ITGAX"]

#Neutrophil_marker =["G0S2", "CXCL8", "SOD2", "NAMPT"]

TNK_marker =  ["CD3D", "CD3E", "CD3G", "TRAC", "NKG7", "XCL2", "XCL1"]

B_marker = ["MS4A1","CD79A","CD79B", "CD19", "BLNK", "FCRL5", "CD27", "BANK1"]

Plasma_marker = ["SDC1", "MZB1", "XBP1","JCHAIN","IGHG1"]

Erythrocyte_marker = ["HBA1","HBA2","HBB"]

#Neural_marker = ["S100B","NGFR","TUBB2B"]


sc.pl.violin(combined_adata, Epi_marker, groupby="leiden", stripplot=False)
sc.pl.violin(combined_adata, Fib_marker, groupby="leiden", stripplot=False)
#sc.pl.violin(combined_adata, Mural_marker, groupby="leiden", stripplot=False)
sc.pl.violin(combined_adata, Endo_marker, groupby="leiden", stripplot=False)
sc.pl.violin(combined_adata, Mac_marker, groupby="leiden", stripplot=False)
sc.pl.violin(combined_adata, Mono_marker, groupby="leiden", stripplot=False)
sc.pl.violin(combined_adata, Mast_marker, groupby="leiden", stripplot=False)
sc.pl.violin(combined_adata, Dendritic_marker, groupby="leiden", stripplot=False)
#sc.pl.violin(combined_adata, Neutrophil_marker, groupby="leiden", stripplot=False)
sc.pl.violin(combined_adata, TNK_marker, groupby="leiden", stripplot=False)
sc.pl.violin(combined_adata, B_marker, groupby="leiden", stripplot=False)
sc.pl.violin(combined_adata, Plasma_marker, groupby="leiden", stripplot=False)
sc.pl.violin(combined_adata, Erythrocyte_marker, groupby="leiden", stripplot=False)
#sc.pl.violin(combined_adata, Neural_marker, groupby="leiden", stripplot=False)


sc.tl.rank_genes_groups(combined_adata, "leiden", groups=["27"], method="wilcoxon")
result = combined_adata.uns['rank_genes_groups']
cluster_27_genes = pd.DataFrame({
    'gene': result['names']['27'],
    'logfoldchanges': result['logfoldchanges']['27'],
    'pvals': result['pvals']['27'],
    'pvals_adj': result['pvals_adj']['27']
})
cluster_27_genes['gene'].head(30).to_list()


cell_annotation_level1 = {
'0':	'Macrophage',
'1':	'T/NK cell',
'2':	'Macrophage',
'3':	'T/NK cell',
'4':	'Macrophage',
'5':	'Endothelial cell',
'6':	'Macrophage',
'7':	'Epithelial cell',
'8':	'Macrophage',
'9':	'Epithelial cell',
'10':	'Epithelial cell',
'11':	'T/NK cell',
'12':	'Fibroblast',
'13':	'Monocyte',
'14':	'Macrophage',
'15':	'Macrophage',
'16':	'Epithelial cell',
'17':	'Mast cell',
'18':	'B cell',
'19':	'Epithelial cell',
'20':	'T/NK cell',
'21':	'Plasma cell',
'22':	'Endothelial cell',
'23':	'Epithelial cell',
'24':	'Macrophage',
'25':	'Macrophage',
'26':	'Macrophage',
'27':	'Epithelial cell',
'28':	'Dendritic cell',
'29':	'Erythrocyte',
'30':	'Monocyte',
'31':	'Mast cell',
'32':	'T/NK cell',

}
leiden_annotation_level1 = pd.Series(cell_annotation_level1)
combined_adata.obs['cell_type_level1'] = combined_adata.obs['leiden'].map(leiden_annotation_level1)

cell_color = {
    "Epithelial cell":"#e43030",
    "Cholangiocyte":"#c67e82",
    "Hepatocyte":"#ddcccc",
    "Fibroblast":"#ff8831",
    "Mural cell":"#20452e",
    "Endothelial cell":"#704ba3",
    "Macrophage":"#67a8cd",
    "Monocyte":"#8184e2",
    "Mast cell":"#ff9d9f",
    "Dendritic cell":"#cf9f88",
    "Neutrophil":"#ffc17f",
    "T/NK cell":"#50aa4b",
    "B cell":"#b3e19b",
    "Plasma cell":"#ab3181",
    "Erythrocyte":"#a1acbd",
    "Neural cell":"#00afba"
}
#combined_adata.write("/storage/data/KAI/Pan_Inflammation/Lung/Lung_py/scvi/" + "Lung_scvi_cell_level1.h5ad")

combined_adata = sc.read_h5ad("/storage/data/KAI/Pan_Inflammation/Lung/Lung_py/scvi/" + "Lung_scvi_cell_level1.h5ad")

fig, ax=plt.subplots(1,1,figsize=(5,4))
sc.pl.umap(combined_adata, color="cell_type_level1",
           legend_loc="right margin",
           palette=cell_color,
           add_outline=False,
           legend_fontsize=8,
           legend_fontoutline=0.1,
           legend_fontweight="normal",
           title="Lung_scvi_cell_annotation",
           frameon=False,
           show=False,
           ax=ax
           #save = "_Lung_scvi_cell_level1_legend.pdf"
           )
# fig = plt.gcf()
for collection in ax.collections:
    collection.set_rasterized(True)
# for ax in fig.get_axes():

plt.savefig("/storage/data/KAI/Pan_Inflammation/Infla_py/figures/umap_Lung_scvi_cell_level1_legend.pdf", format="pdf", dpi=600)
plt.close()

fig, ax=plt.subplots(1,1,figsize=(5,4))
sc.pl.umap(combined_adata, color="cell_type_level1",
           legend_loc="on data",
           palette=cell_color,
           add_outline=False,
           legend_fontsize=8,
           legend_fontoutline=0.1,
           legend_fontweight="normal",
           title="Lung_scvi_cell_annotation",
           frameon=False,
           show=False,
           ax=ax
           #save = "_Lung_scvi_cell_level1_clear.pdf"
           )
# fig = plt.gcf()
for collection in ax.collections:
    collection.set_rasterized(True)
# for ax in fig.get_axes():

plt.savefig("/storage/data/KAI/Pan_Inflammation/Infla_py/figures/umap_Lung_scvi_cell_level1_clear.pdf", format="pdf", dpi=600)
plt.close()



# scanvi-----------------------

lvae = scvi.model.SCANVI.from_scvi_model(
    vae,
    adata=combined_adata,
    labels_key="cell_type_level1",
    unlabeled_category="Unknown",
)
lvae.train()
model_dir2 = os.path.join("/storage/data/KAI/Pan_Inflammation/Lung/Lung_py/scanvi/model/lung_allcells/", "lung_allcells_scanvi_model")
lvae.save(model_dir2, overwrite=True, prefix="lung_allcells_scanvi_")

lvae = scvi.model.SCANVI.load("/storage/data/KAI/Pan_Inflammation/Lung/Lung_py/scanvi/model/lung_allcells/lung_allcells_scanvi_model/", adata=combined_adata, prefix="lung_allcells_scanvi_")


combined_adata.obsm["X_scANVI"] = lvae.get_latent_representation()

# run PcA then generate UMAp plots
sc.tl.pca(combined_adata)
sc.pl.pca_variance_ratio(combined_adata,log=True)
# use scANVI latent space for UMAp generation
sc.pp.neighbors(combined_adata,use_rep="X_scANVI")
sc.tl.umap(combined_adata,min_dist=0.3)


fig, ax=plt.subplots(1,1,figsize=(5,4))
sc.pl.umap(combined_adata, color="cell_type_level1",
           legend_loc="right margin",
           palette=cell_color,
           add_outline=False,
           legend_fontsize=8,
           legend_fontoutline=0.1,
           legend_fontweight="normal",
           title="Lung_scanvi_cell_annotation",
           frameon=False,
           show=False,
           ax=ax
           #save = "_Lung_scanvi_cell_level1_legend.pdf"
           )
# fig = plt.gcf()
for collection in ax.collections:
    collection.set_rasterized(True)
# for ax in fig.get_axes():

plt.savefig("/storage/data/KAI/Pan_Inflammation/Infla_py/figures/umap_Lung_scanvi_cell_level1_legend.pdf", format="pdf", dpi=600)
plt.close()

fig, ax=plt.subplots(1,1,figsize=(5,4))
sc.pl.umap(combined_adata, color="cell_type_level1",
           legend_loc="on data",
           palette=cell_color,
           add_outline=False,
           legend_fontsize=8,
           legend_fontoutline=0.1,
           legend_fontweight="normal",
           title="Lung_scanvi_cell_annotation",
           frameon=False,
           show=False,
           ax=ax
           #save = "_Lung_scanvi_cell_level1_clear.pdf"
           )
# fig = plt.gcf()
for collection in ax.collections:
    collection.set_rasterized(True)
# for ax in fig.get_axes():

plt.savefig("/storage/data/KAI/Pan_Inflammation/Infla_py/figures/umap_Lung_scanvi_cell_level1_clear.pdf", format="pdf", dpi=600)
plt.close()


combined_adata.write("/storage/data/KAI/Pan_Inflammation/Lung/Lung_py/scanvi/Lung_scanvi.h5ad")


fig, ax=plt.subplots(1,1,figsize=(5,4))
sc.pl.umap(combined_adata,color=["GEO"], ncols=2, frameon=False)
# fig = plt.gcf()
for collection in ax.collections:
    collection.set_rasterized(True)
# for ax in fig.get_axes():

plt.savefig("/storage/data/KAI/Pan_Inflammation/Infla_py/figures/lung_scanvi_batch.pdf", format="pdf", dpi=600)
plt.close()





violin_markers = [
  "EPCAM","KRT18","KRT19",
  "DCN", "LUM","COL1A1",
  "VWF","CDH5","PECAM1",
  "C1QA","AIF1", "MARCO",
  "S100A8","S100A9",
  "IRF8","CD83","LILRA4",
  "TPSAB1", "KIT","CPA3",
  "CD3D", "CD3E", "NKG7",
  "MS4A1","BANK1","CD79A",
  "MZB1","JCHAIN","IGHG1",
  "HBA1","HBA2","HBB"
                  ]

# TO R
combined_adata = sc.read_h5ad("/storage/data/KAI/Pan_Inflammation/Lung/Lung_py/scanvi/Lung_scanvi.h5ad")

combined_adata.X = combined_adata.layers['counts']
combined_adata.write("/storage/data/KAI/Pan_Inflammation/Lung/Lung_py/scanvi/Lung_scanvi_toR.h5ad")

