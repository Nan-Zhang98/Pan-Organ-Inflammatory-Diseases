#### Moudle --------------------------
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
import skimage

from typing import Sequence,List,Union
import gseapy
from gseapy import Biomart

import scarches as sca


#### load data -----------------------
file_path = '/storage/data/KAI/Pan_Inflammation/Mouse/Pneumonia/LC/counts/'
control = list(filter(lambda file: file.startswith(("con")), os.listdir(file_path)))
pl1 = list(filter(lambda file: file.startswith(("PL1")), os.listdir(file_path)))
pl2 = list(filter(lambda file: file.startswith(("PL2")), os.listdir(file_path)))

def gunzip_all_gz_files(file_path):
    # 获取当前路径下所有文件
    files = os.listdir(file_path)

    # 遍历文件列表，查找 .gz 文件
    for file in files:
        if file.endswith('.gz'):
            # 调用 gunzip 解压文件
            subprocess.run(['gunzip', os.path.join(file_path, file)])

def rename_features_to_genes(file_path):
    for filename in os.listdir(file_path):
        if filename.startswith("features"):
            # 分离文件名和扩展名
            name, ext = os.path.splitext(filename)
            new_filename = f"genes{ext}"
            # 重命名文件
            os.rename(os.path.join(file_path, filename), os.path.join(file_path, new_filename))
            print(f'Renamed: {filename} to {new_filename}')

adata_list = []

### con
for gsm_folder in control:
    gsm_path = os.path.join(file_path, gsm_folder)

    if os.path.isdir(gsm_path):
        if any(file.endswith('.gz') for file in os.listdir(gsm_path)):
            # 识别以.gz结尾的matrix、features和barcodes文件，并调用linux命令进行解压
            gunzip_all_gz_files(gsm_path)

        if any(file.startswith('features') for file in os.listdir(gsm_path)):
            # 将features.tsv更名为genes.tsv，防止报错
            rename_features_to_genes(gsm_path)

    adata = sc.read_10x_mtx(gsm_path, var_names="gene_symbols")
    adata.obs_names = gsm_folder + '_' + adata.obs_names

    # 删去重复的基因
    unique_genes = pd.Series(adata.var_names).drop_duplicates(keep="first")
    adata = adata[adata.obs.index, unique_genes.index]
    adata.var_names = unique_genes.values
    adata.obs['Inflammation'] = 'Normal'
    adata.obs['Sample'] = gsm_folder
    adata_list.append(adata)

### pl1
for gsm_folder in pl1:
    gsm_path = os.path.join(file_path, gsm_folder)

    if os.path.isdir(gsm_path):
        if any(file.endswith('.gz') for file in os.listdir(gsm_path)):
            # 识别以.gz结尾的matrix、features和barcodes文件，并调用linux命令进行解压
            gunzip_all_gz_files(gsm_path)

        if any(file.startswith('features') for file in os.listdir(gsm_path)):
            # 将features.tsv更名为genes.tsv，防止报错
            rename_features_to_genes(gsm_path)

    adata = sc.read_10x_mtx(gsm_path, var_names="gene_symbols")
    adata.obs_names = gsm_folder + '_' + adata.obs_names

    # 删去重复的基因
    unique_genes = pd.Series(adata.var_names).drop_duplicates(keep="first")
    adata = adata[adata.obs.index, unique_genes.index]
    adata.var_names = unique_genes.values
    adata.obs['Inflammation'] = 'PL_1W'
    adata.obs['Sample'] = gsm_folder
    adata_list.append(adata)

### pl2
for gsm_folder in pl2:
    gsm_path = os.path.join(file_path, gsm_folder)

    if os.path.isdir(gsm_path):
        if any(file.endswith('.gz') for file in os.listdir(gsm_path)):
            # 识别以.gz结尾的matrix、features和barcodes文件，并调用linux命令进行解压
            gunzip_all_gz_files(gsm_path)

        if any(file.startswith('features') for file in os.listdir(gsm_path)):
            # 将features.tsv更名为genes.tsv，防止报错
            rename_features_to_genes(gsm_path)

    adata = sc.read_10x_mtx(gsm_path, var_names="gene_symbols")
    adata.obs_names = gsm_folder + '_' + adata.obs_names

    # 删去重复的基因
    unique_genes = pd.Series(adata.var_names).drop_duplicates(keep="first")
    adata = adata[adata.obs.index, unique_genes.index]
    adata.var_names = unique_genes.values
    adata.obs['Inflammation'] = 'PL_2W'
    adata.obs['Sample'] = gsm_folder
    adata_list.append(adata)

### Integrate
combined_adata = adata_list[0].concatenate(*adata_list[1:], join='outer', batch_key=None)
combined_adata.write("/storage/data/KAI/Pan_Inflammation/Mouse/Pneumonia/LC/pl_py/plmouse_scanpy.h5ad")

#### preprocess ----------------------
adata = sc.read_h5ad("/storage/data/KAI/Pan_Inflammation/Mouse/Pneumonia/LC/pl_py/plmouse_scanpy.h5ad")
# mitochondrial genes, "MT-" for human, "Mt-" for mouse
adata.var["mt"] = adata.var_names.str.startswith("MT-")
# ribosomal genes
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
# hemoglobin genes
adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")
sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True
)
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata = adata[adata.obs['pct_counts_mt'] < 20, :]

sc.pp.scrublet(adata, batch_key="Sample")
adata = adata[~adata.obs['predicted_doublet'],:]
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
adata.obs['Organ'] = 'Lung'

adata.write("/storage/data/KAI/Pan_Inflammation/Mouse/Pneumonia/LC/pl_py/plmouse_preprocess.h5ad")

# run PcA then generate UMAp plots
sc.tl.pca(adata)
sc.pl.pca_variance_ratio(adata,log=True)
# use scvI latent space for UMAp generation
sc.pp.neighbors(adata,use_rep="X_pca")
sc.tl.umap(adata,min_dist=0.3)
sc.pl.umap(adata,color=["Sample","Inflammation"], ncols=2, frameon=False, save = "PLmouse_scvi_batch.png")



#### scvi ---------------
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

scvi.model.SCVI.setup_anndata(adata, layer="counts", categorical_covariate_keys=["Sample"])
vae = scvi.model.SCVI(adata, n_latent=20, n_layers=1, dropout_rate=0.1)
vae.train(max_epochs=50, early_stopping=True, batch_size=128,early_stopping_patience=5)

# 保存scvi模型
model_dir = os.path.join()
vae.save(model_dir, overwrite=True, prefix="plmouse_scvi_")

vae = scvi.model.SCVI.load("/storage/data/KAI/Pan_Inflammation/Mouse/Pneumonia/LC/pl_py/scvi/model/", adata=adata, prefix="plmouse_scvi_")

adata.obsm["X_scVI"] = vae.get_latent_representation()
adata.layers["X_normalized_scVI"] = vae.get_normalized_expression(library_size=10e4)


# use scvI latent space for UMAp generation
sc.pp.neighbors(adata,use_rep="X_scVI")
sc.tl.umap(adata,min_dist=0.3)

# PLmouse_scvi_script.py
# nohup /storage/data/KAI/miniconda3/envs/XKscvi/bin/python /storage/data/KAI/Pan_Inflammation/PyProject/pythonProject/Python_inflammation/PLmouse_scvi_script.py > /storage/data/KAI/Pan_Inflammation/Mouse/Pneumonia/LC/pl_py/PLmouse_scvi_log.txt 2>&1 &

# PID =  2150223
# Time 5.4 9:40
# done


# load scvi anndata
adata = sc.read_h5ad("/storage/data/KAI/Pan_Inflammation/Mouse/Pneumonia/LC/pl_py/scvi/plmouse_scvi.h5ad")

sc.pl.umap(adata,color=["Sample","Inflammation"], ncols=2, frameon=False, save = "PLmouse_afterscvi_batch.pdf"
            )

print(os.getcwd()) # 查找工作路径
# os.chdir('/storage/data/KAI/Pan_Inflammation/PyProject/') # 更改工作路径


#### annotation cell_type_level1
sc.tl.leiden(adata, resolution=0.8)
sc.pl.umap(adata, color=["leiden"],
           legend_loc='on data',
           # save="PLmouse_leiden.png"
           )

Macrophage_marker = ["Cd68","Adgre1","Csf1r","Lyz2","Itgam","Fcgr1"]
Neutrophil_marker = ["S100a8","S100a9","Mmp8","Cxcr2","Camp","Ngp"]
Monocyte_marker = ["Plac8","Ly6c2","Ccr2"] # done
Mast_marker = ["Cpa3","Fcer1a","Mcpt8"] # done
Dendritic_marker = ["Ccl17","Vdr","Ccl22","Fcrls"]

# ILC_marker = ['Il7r','Cxcr6']

Endothelial_marker = ["Pecam1","Cldn5","Kdr",'Cdh5','Vwf','Eng'] # done
Epithelial_marker = ["Epcam","Cdh1",'Scgb1a1','Foxj1','Sftpc'] # done
Fibroblast_marker = ["Mgp","Gsn","Inmt","Col3a1","Dcn","Col1a1","Col1a2"]
# Mural_marker = ['Rgs5','Pdgfrb','Kcnj8','Acta2','Vtn']

TNK_marker = ["Cd3g","Cd3e","Cd3d","Ccl5","Gzma","Gzmb","Nkg7"]

B_marker = ["Cd79a","Cd79b","Ms4a1","Cd19"]
Plasma_marker = ["Mzb1","Jchain","Cd27","Xbp1"]

# Neural_marker = ['Tubb3','Snap25','Rbfox3','Meg3','Syt1','Map2','Stmn2']
# Erythrocytes_marker = []


sc.pl.violin(adata, Macrophage_marker, groupby="leiden", stripplot=False, show=True, use_raw=False, rotation=90)
sc.pl.violin(adata, Monocyte_marker, groupby="leiden", stripplot=False, show=True, use_raw=False, rotation=90)
sc.pl.violin(adata, Neutrophil_marker, groupby="leiden", stripplot=False, show=True, use_raw=False, rotation=90)
sc.pl.violin(adata, Dendritic_marker, groupby="leiden", stripplot=False, show=True, use_raw=False, rotation=90)
sc.pl.violin(adata, Mast_marker, groupby="leiden", stripplot=False, show=True, use_raw=False, rotation=90)

# sc.pl.violin(adata, ILC_marker, groupby="leiden", stripplot=False, show=True, use_raw=False, rotation=90)

sc.pl.violin(adata, Endothelial_marker, groupby="leiden", stripplot=False, show=True, use_raw=False, rotation=90)
sc.pl.violin(adata, Epithelial_marker, groupby="leiden", stripplot=False, show=True, use_raw=False, rotation=90)
sc.pl.violin(adata, Fibroblast_marker, groupby="leiden", stripplot=False, show=True, use_raw=False, rotation=90)
# sc.pl.violin(adata, Mural_marker, groupby="leiden", stripplot=False, show=True, use_raw=False, rotation=90)


sc.pl.violin(adata, TNK_marker, groupby="leiden", stripplot=False, show=True, use_raw=False, rotation=90)
sc.pl.violin(adata, B_marker, groupby="leiden", stripplot=False, show=True, use_raw=False, rotation=90)
sc.pl.violin(adata, Plasma_marker, groupby="leiden", stripplot=False, show=True, use_raw=False, rotation=90)

# sc.pl.violin(adata, Neural_marker, groupby="leiden", stripplot=False, show=True, use_raw=False, rotation=90)


sc.pl.umap(adata,
           color="Cd68",
           use_raw = False,
           #palette=macro_color,
           show=True,
           size=5
           )
sc.pl.umap(adata,
           color="Cd14",
           use_raw = False,
           #palette=macro_color,
           show=True,
           size=5
           )


sc.tl.rank_genes_groups(adata, "leiden", groups=["11"], method="wilcoxon")
result = adata.uns['rank_genes_groups']
cluster_11_genes = pd.DataFrame({
    'gene': result['names']['11'],
    'logfoldchanges': result['logfoldchanges']['11'],
    'pvals': result['pvals']['11'],
    'pvals_adj': result['pvals_adj']['11']
})
cluster_11_genes['gene'].head(50).to_list()


adata.write("/storage/data/KAI/Pan_Inflammation/Mouse/Pneumonia/LC/pl_py/scvi/plmouse_scvi.h5ad")

# Annotation
cell_annotation_level1 = {
'0':	'Neutrophil',
'1':	'Macrophage',
'2':	'T cell',
'3':	'B cell',
'4':	'Monocyte',
'5':	'Endothelial cell',
'6':	'Macrophage',
'7':	'Fibroblast',
'8':	'T cell',
'9':	'Dendritic cell',
'10':	'Endothelial cell',
'11':	'Macrophage',
'12':	'NK cell',
'13':	'T cell',
'14':	'Macrophage',
'15':	'Epithelial cell',
'16':	'Macrophage',
'17':	'T cell',
'18':	'Epithelial cell',
'19':	'Endothelial cell',
'20':	'Fibroblast',
'21':	'Endothelial cell',
'22':	'Neutrophil',
'23':	'T cell',
'24':	'Plasma cell',
'25':	'Neutrophil',
'26':	'Fibroblast',
'27':	'Endothelial cell',
'28':	'Mast cell',
'29':	'Endothelial cell',
'30':	'Fibroblast',
}
leiden_annotation_level1 = pd.Series(cell_annotation_level1)
adata.obs['cell_type_level1'] = adata.obs['leiden'].map(leiden_annotation_level1)

#### scANVI -----------------------
lvae = scvi.model.SCANVI.from_scvi_model(
    vae,
    adata=adata,
    labels_key="cell_type_level1",
    unlabeled_category="Unknown",
)
lvae.train()

lvae.save('/storage/data/KAI/Pan_Inflammation/Mouse/Pneumonia/LC/pl_py/scanvi/', overwrite=True, prefix="plmouse_scanvi_")

# umap
adata.obsm["X_scANVI"] = lvae.get_latent_representation()

# use scvI latent space for UMAp generation
sc.pp.neighbors(adata,use_rep="X_scANVI")
sc.tl.umap(adata,min_dist=0.3)

adata.write("/storage/data/KAI/Pan_Inflammation/Mouse/Pneumonia/LC/pl_py/scanvi/plmouse_scanvi.h5ad")


# PLmouse_scvi_script.py
# nohup /storage/data/KAI/miniconda3/envs/XKscvi/bin/python /storage/data/KAI/Pan_Inflammation/PyProject/pythonProject/Python_inflammation/PLmouse_scvi_script.py > /storage/data/KAI/Pan_Inflammation/Mouse/Pneumonia/LC/pl_py/PLmouse_scanvi_log.txt 2>&1 &

# PID = 3084803
# Time 5.7 13:55
# done

# scanvi visual
adata = sc.read_h5ad("/storage/data/KAI/Pan_Inflammation/Mouse/Pneumonia/LC/pl_py/scanvi/plmouse_scanvi.h5ad")
print(os.getcwd()) # 查找工作路径
# os.chdir('/storage/data/KAI/Pan_Inflammation/PyProject/') # 更改工作路径

sc.pl.umap(adata,color=["Sample","Inflammation"], ncols=2, frameon=False,  save = "PLmouse_afterscanvi_batch.pdf"
            )

cell_color = {
    "Epithelial cell":"#e43030",
    "Fibroblast":"#ff8831",
    "Endothelial cell":"#704ba3",
    "Macrophage":"#67a8cd",
    "Monocyte":"#8184e2",
    "Mast cell":"#ff9d9f",
    "Dendritic cell":"#cf9f88",
    "Neutrophil":"#ffc17f",
    "T cell":"#50aa4b",
    "B cell":"#b3e19b",
    "Plasma cell":"#ab3181",
    "NK cell":"#00afba"
}


fig, ax=plt.subplots(1,1,figsize=(5,4))
sc.pl.umap(adata, color="cell_type_level1",
           legend_loc="on data",
           palette=cell_color,
           add_outline=False,
           legend_fontsize=8,
           legend_fontoutline=0.1,
           legend_fontweight="normal",
           title="Mouse_Lung_scanvi_cell_annotation",
           frameon=False,
           show=False,
           ax=ax
           #save = "MouseLung_scanvi_cell_level1.pdf"
           )
# fig = plt.gcf()
for collection in ax.collections:
    collection.set_rasterized(True)
# for ax in fig.get_axes():

plt.savefig("/storage/data/KAI/Pan_Inflammation/Mouse/Pneumonia/LC/pl_py/scanvi/umap_Mouse_Lung_scanvi_cell_level1.pdf", format="pdf", dpi=600)
plt.close()

#### marker 展示 ----------------------------
level1_marker = [
"Cd68","Adgre1","Mrc1", # Macrophage
"Plac8","Ly6c2","Ccr2", # Monocyte
"S100a8","S100a9","Cxcr2", # Neutrophil
"Ccl17","Vdr","Ccl22", # Dendritic
"Cpa3","Fcer1a","Mcpt8", # Mast

"Pecam1",'Cdh5','Vwf', # Endothelial
"Epcam",'Scgb1a1','Foxj1','Sftpc', # Epithelial
"Mgp","Gsn","Col3a1", # Fibroblast
"Cd3g","Cd3e","Cd3d", # T/NK
"Gzma","Gzmb","Nkg7", # T/NK
"Cd79a","Cd79b","Cd19", # B
"Mzb1","Jchain","Xbp1", # Plasma
]

# 为dotplot重新排序
desired_order = [
"Macrophage",
"Monocyte",
"Neutrophil",
"Dendritic cell",
"Mast cell",
"Endothelial cell",
"Epithelial cell",
"Fibroblast",
"T cell",
"NK cell",
"B cell",
"Plasma cell",
]
adata.obs['cell_type_level1'] = pd.Categorical(
    adata.obs['cell_type_level1'],
    categories=desired_order,
    ordered=True
)
print(os.getcwd()) # 查找工作路径
# os.chdir('/storage/data/KAI/Pan_Inflammation/PyProject/') # 更改工作路径
sc.pl.dotplot(
    adata,
    var_names=level1_marker,       # marker 基因列表
    groupby='cell_type_level1',          # 根据细胞类型分组
    cmap='pink_r',               # 色彩映射，可选 'viridis', 'magma' 等
    dendrogram=False,              # 是否显示分组的树状图
    swap_axes=True,               # 是否交换行列以适应展示
    standard_scale="var",
    show=True,
    save='dotplot_Mouse_Lung_level1_markers.pdf'
)

adata.X = adata.layers['counts'].copy()

adata.write("/storage/data/KAI/Pan_Inflammation/Mouse/Pneumonia/LC/pl_py/scanvi/plmouse_scanvi_count.h5ad")



#### 提取巨噬细胞和单核细胞 --------------------------
adata = sc.read_h5ad("/storage/data/KAI/Pan_Inflammation/Mouse/Pneumonia/LC/pl_py/scanvi/plmouse_scanvi_count.h5ad")

monoMac = adata[adata.obs['cell_type_level1'].isin(['Monocyte', 'Macrophage'])].copy()
monoMac.write("/storage/data/KAI/Pan_Inflammation/Mouse/Pneumonia/LC/analysis/scArches/data/plmouse_monomac.h5ad")

### 整合重聚类
monoMac = sc.read_h5ad("/storage/data/KAI/Pan_Inflammation/Mouse/Pneumonia/LC/analysis/scArches/data/plmouse_monomac.h5ad")
sc.pp.normalize_total(monoMac, target_sum=1e4)  # 标准化到 10,000
sc.pp.log1p(monoMac)  # 对数变换

scvi.model.SCVI.setup_anndata(monoMac, layer="counts", categorical_covariate_keys=["Sample","Inflammation"])
vae = scvi.model.SCVI(monoMac, n_latent=30, n_layers=2, dropout_rate=0.1)
vae.train(max_epochs=80, early_stopping=True, batch_size=128,early_stopping_patience=5)
vae.save('/storage/data/KAI/Pan_Inflammation/Mouse/Pneumonia/LC/analysis/scArches/data/scvi/model', overwrite=True, prefix="plmouse_monoMac_scvi_")


monoMac.obsm["X_scVI"] = vae.get_latent_representation()
monoMac.layers["X_normalized_scVI"] = vae.get_normalized_expression(library_size=10e4)
sc.pp.neighbors(monoMac,use_rep="X_scVI")
sc.tl.umap(monoMac,min_dist=0.3)

sc.pl.umap(monoMac, color="cell_type_level1",
           legend_loc="on data",
           palette=cell_color,
           add_outline=False,
           legend_fontsize=8,
           legend_fontoutline=0.1,
           legend_fontweight="normal",
           title="Mouse_Lung_scanvi_cell_annotation",
           frameon=False,
           show=True,
           #save = "MouseLung_monoMac_scvi_cell_level1.pdf"
           )

sc.pl.umap(monoMac,color=["Sample","Inflammation"], ncols=2, frameon=False,  # save = "PLmouse_monoMac_afterscvi_batch.pdf"
          )



# 提高resolution，删除不高表达成纤维细胞marker和高表达其他细胞类型marker的簇
sc.tl.leiden(monoMac, resolution=6.0, key_added='leiden_res_6.0')
sc.tl.leiden(monoMac, resolution=10.0, key_added='leiden_res_10.0')

monoMac.write("/storage/data/KAI/Pan_Inflammation/Mouse/Pneumonia/LC/analysis/scArches/data/scvi/plmouse_monomac_all.h5ad")

# PLmouse_scvi_script.py
# nohup /storage/data/KAI/miniconda3/envs/XKscvi/bin/python /storage/data/KAI/Pan_Inflammation/PyProject/pythonProject/Python_inflammation/PLmouse_scvi_script.py > /storage/data/KAI/Pan_Inflammation/Mouse/Pneumonia/LC/pl_py/PLmouse_monoMac_scvi_log.txt 2>&1 &
# PID = 3835309
# Time 5.8 20:49
# done

monoMac = sc.read_h5ad("/storage/data/KAI/Pan_Inflammation/Mouse/Pneumonia/LC/analysis/scArches/data/scvi/plmouse_monomac_all.h5ad")


# leiden res = 10.0
sc.pl.umap(monoMac, color=["leiden_res_10.0"],
           size = 3,
           legend_loc='on data',
           )

# 过滤低表达monocyte和macrophage marker和表达其他细胞marker的cluster
cell_marker = [
"Cd68","Adgre1","Mrc1", # Macrophage
"Plac8","Ly6c2","Ccr2", # Monocyte
"S100a8","S100a9","Cxcr2", # Neutrophil
"Ccl17","Vdr","Ccl22", # Dendritic
"Cpa3","Fcer1a","Mcpt8", # Mast

"Pecam1",'Cdh5','Vwf', # Endothelial
"Epcam",'Scgb1a1','Foxj1','Sftpc', # Epithelial
"Mgp","Gsn","Col3a1", # Fibroblast
"Cd3g","Cd3e","Cd3d", # T/NK
"Gzma","Gzmb","Nkg7", # T/NK
"Cd79a","Cd79b","Cd19", # B
"Mzb1","Jchain","Xbp1", # Plasma
"Hba-a1","Hba-a2","Hbb-bs","Hbb-bt" # Erythrocyte
]

sc.set_figure_params(figsize=(25, 5))
dp = sc.pl.dotplot(monoMac, cell_marker, groupby="leiden_res_10.0", return_fig=True, swap_axes=False, standard_scale="var",title="")
# 添加总计和样式设置
dp.add_totals().style(dot_edge_color='black', dot_edge_lw=0.5).show()

# 待删除cluster
del_leiden_100_1 = [
    '21','26','42','77','86','97','126','139','143','160','162' # low Mac & mono
    '15','17','24','85','90','134', # high T doublet
    '128', # high B doublet
    '34', # high fibroblast
    '159', # high Erythrocyte
]

monoMac = monoMac[~monoMac.obs['leiden_res_10.0'].isin(del_leiden_100_1)]
monoMac = monoMac.copy()
monoMac.write("/storage/data/KAI/Pan_Inflammation/Mouse/Pneumonia/LC/analysis/scArches/data/scvi/plmouse_monomac_filtered.h5ad")

# 重新进行scVI整合
monoMac = sc.read_h5ad("/storage/data/KAI/Pan_Inflammation/Mouse/Pneumonia/LC/analysis/scArches/data/scvi/plmouse_monomac_filtered.h5ad")

sc.pp.normalize_total(monoMac, target_sum=1e4)  # 标准化到 10,000
sc.pp.log1p(monoMac)  # 对数变换

scvi.model.SCVI.setup_anndata(monoMac, layer="counts", categorical_covariate_keys=["Sample","Inflammation"])
vae = scvi.model.SCVI(monoMac, n_latent=30, n_layers=2, dropout_rate=0.1)
vae.train(max_epochs=80, early_stopping=True, batch_size=128,early_stopping_patience=5)
vae.save('/storage/data/KAI/Pan_Inflammation/Mouse/Pneumonia/LC/analysis/scArches/data/filtered_scvi/model/', overwrite=True, prefix="plmouse_monoMac_filtered_scvi_")

monoMac.obsm["X_scVI"] = vae.get_latent_representation()
monoMac.layers["X_normalized_scVI"] = vae.get_normalized_expression(library_size=10e4)
sc.pp.neighbors(monoMac,use_rep="X_scVI")
sc.tl.umap(monoMac,min_dist=0.1)

monoMac.X = monoMac.layers['counts'].copy()

monoMac.write("/storage/data/KAI/Pan_Inflammation/Mouse/Pneumonia/LC/analysis/scArches/data/filtered_scvi/plmouse_monomac_filtered_scvi.h5ad")

# scANVI调整
# lvae = scvi.model.SCANVI.from_scvi_model(
#     vae,
#     adata=monoMac,
#     labels_key="cell_type_level1",
#     unlabeled_category="Unknown",
# )
# lvae.train()
# lvae.save('/storage/data/KAI/Pan_Inflammation/Mouse/Pneumonia/LC/analysis/scArches/data/filtered_scanvi/model/', overwrite=True, prefix="plmouse_monoMac_filtered_scanvi_")
#
# # umap
# monoMac.obsm["X_scANVI"] = lvae.get_latent_representation()
#
# # use scvI latent space for UMAp generation
# sc.pp.neighbors(monoMac,use_rep="X_scANVI")
# sc.tl.umap(monoMac,min_dist=0.1)
#
# monoMac.write("/storage/data/KAI/Pan_Inflammation/Mouse/Pneumonia/LC/analysis/scArches/data/filtered_scanvi/plmouse_monomac_filtered_scanvi.h5ad")


# PLmouse_scvi_script.py
# nohup /storage/data/KAI/miniconda3/envs/XKscvi/bin/python /storage/data/KAI/Pan_Inflammation/PyProject/pythonProject/Python_inflammation/PLmouse_scvi_script.py > /storage/data/KAI/Pan_Inflammation/Mouse/Pneumonia/LC/pl_py/PLmouse_monoMac_scvi_log.txt 2>&1 &
# PID = 117079
# Time 5.10 10:49
# done

monoMac = sc.read_h5ad("/storage/data/KAI/Pan_Inflammation/Mouse/Pneumonia/LC/analysis/scArches/data/filtered_scvi/plmouse_monomac_filtered_scvi.h5ad")


fig, ax=plt.subplots(1,1,figsize=(5,4))
sc.pl.umap(monoMac, color="cell_type_level1",
           legend_loc="on data",
           palette=cell_color,
           add_outline=False,
           legend_fontsize=8,
           legend_fontoutline=0.1,
           legend_fontweight="normal",
           title="Mouse_Lung_scvi_monocyte&macrophage",
           frameon=False,
           show=False,
           ax=ax
           #save = "MouseLung_scanvi_cell_level1.pdf"
           )
# fig = plt.gcf()
for collection in ax.collections:
    collection.set_rasterized(True)
# for ax in fig.get_axes():

plt.savefig("/storage/data/KAI/Pan_Inflammation/Mouse/Pneumonia/LC/analysis/scArches/data/filtered_scvi/umap_Mouse_Lung_monoMac_scvi_level1.pdf", format="pdf", dpi=600)
plt.close()


sc.pl.umap(monoMac,
           color=["Sample","Inflammation"],
           ncols=2,
           frameon=False,
           # save = "PLmouse_monoMac_afterscvi_batch.pdf"
          )

# Alveolar macrophage marker
sc.pl.umap(monoMac,
           color="Itgax",
           use_raw = False,
           #palette=macro_color,
           show=True,
           size=5,
           )
sc.pl.umap(monoMac,
           color="Siglecf",
           use_raw = False,
           #palette=macro_color,
           show=True,
           size=5,
           )
sc.pl.umap(monoMac,
           color="Marco",
           use_raw = False,
           #palette=macro_color,
           show=True,
           size=5,
           )

Mac = monoMac[monoMac.obs['cell_type_level1'].isin(['Macrophage'])].copy()
Mono = monoMac[monoMac.obs['cell_type_level1'].isin(['Monocyte'])].copy()
Mono.write("/storage/data/KAI/Pan_Inflammation/Mouse/Pneumonia/LC/analysis/scArches/data/filtered_scvi/plmouse_mono_filtered_scvi.h5ad")

sc.pl.umap(Mac, color="cell_type_level1",
           legend_loc="on data",
           palette=cell_color,
           add_outline=False,
           legend_fontsize=8,
           legend_fontoutline=0.1,
           legend_fontweight="normal",
           title="Mouse_Lung_scvi_macrophage",
           frameon=False,
           show=True,
           #save = "MouseLung_scanvi_cell_level1.pdf"
           )
Mac.write("/storage/data/KAI/Pan_Inflammation/Mouse/Pneumonia/LC/analysis/scArches/data/filtered_scvi/plmouse_mac_filtered_scvi.h5ad")
Mac = sc.read_h5ad("/storage/data/KAI/Pan_Inflammation/Mouse/Pneumonia/LC/analysis/scArches/data/filtered_scvi/plmouse_mac_filtered_scvi.h5ad")

sc.tl.leiden(Mac, resolution=0.8)
sc.pl.umap(Mac, color=["leiden"],
           legend_loc='on data',
           # save="PLmouse_leiden.png"
           )
pd.crosstab(Mac.obs['leiden'], Mac.obs['Inflammation'])

AM_marker = ['Chil3','Krt79','Plet1','Lpl', # LC's marker
             'Itgax','Siglecf','Marco','Ear1','Ear2', # ACT
             'Mertk',
             ]

sc.pl.violin(Mac, AM_marker, groupby="leiden", stripplot=False, show=True, use_raw=False, rotation=90)

# Macro annotation level2
cell_annotation_level2 = {
'0':	'Mono-derived Macro',
'1':	'Mono-derived Macro',
'2':	'Mono-derived Macro',
'3':	'Alveolar Macro',
'4':	'Mono-derived Macro',
'5':	'Mono-derived Macro',
'6':	'Mono-derived Macro',
'7':	'Mono-derived Macro',
'8':	'Mono-derived Macro',
'9':	'Mono-derived Macro',
'10':	'Alveolar Macro',
'11':	'Mono-derived Macro',
'12':	'Mono-derived Macro',
}

leiden_annotation_level2 = pd.Series(cell_annotation_level2)
Mac.obs['cell_type_level2'] = Mac.obs['leiden'].map(leiden_annotation_level2)

level2_color = {
    "Monocyte":"#8184e2",
    "Alveolar Macro":"#79B472",
    "Mono-derived Macro":"#544E70"
}

fig, ax=plt.subplots(1,1,figsize=(5,4))
sc.pl.umap(Mac, color="cell_type_level2",
           legend_loc="on data",
           palette=level2_color,
           add_outline=False,
           legend_fontsize=8,
           legend_fontoutline=0.1,
           legend_fontweight="normal",
           title="Mouse_Lung_scvi_Macrophage",
           frameon=False,
           show=False,
           ax=ax
           #save = "MouseLung_scanvi_cell_level1.pdf"
           )
# fig = plt.gcf()
for collection in ax.collections:
    collection.set_rasterized(True)
# for ax in fig.get_axes():

plt.savefig("/storage/data/KAI/Pan_Inflammation/Mouse/Pneumonia/LC/analysis/scArches/data/filtered_scvi/umap_Mouse_Lung_scvi_macrophage_level2.pdf", format="pdf", dpi=600)
plt.close()

# AM marker 展示验证
AM_marker = ['Chil3','Krt79','Lpl', # LC's marker
             'Itgax','Siglecf','Marco', # ACT
             'Mertk',]

print(os.getcwd()) # 查找工作路径
# os.chdir('/storage/data/KAI/Pan_Inflammation/PyProject/') # 更改工作路径
sc.pl.dotplot(
    Mac,
    var_names=AM_marker,       # marker 基因列表
    groupby='cell_type_level2',          # 根据细胞类型分组
    cmap='pink_r',               # 色彩映射，可选 'viridis', 'magma' 等
    dendrogram=False,              # 是否显示分组的树状图
    swap_axes=True,               # 是否交换行列以适应展示
    standard_scale="var",
    show=True,
    save='dotplot_Mouse_Macrophage_level2_markers.pdf'
)
Mac.write("/storage/data/KAI/Pan_Inflammation/Mouse/Pneumonia/LC/analysis/scArches/data/filtered_scvi/plmouse_mac_filtered_scvi.h5ad")

level2_meta = Mac.obs['cell_type_level2']
level2_meta.to_csv('/storage/data/KAI/Pan_Inflammation/Mouse/Pneumonia/LC/analysis/annotation/momac_level3/data/level2_meta.txt',sep='\t')




# 保存非AM巨噬细胞（Mono-derived Macro）
moMac = Mac[Mac.obs['cell_type_level2'].isin(['Mono-derived Macro'])].copy()
moMac.X = moMac.layers['counts'].copy()
moMac.write("/storage/data/KAI/Pan_Inflammation/Mouse/Pneumonia/LC/analysis/scArches/data/filtered_scvi/plmouse_moMac_filtered.h5ad")










#### scArches 进行样本迁移映射 -------------------------
### 人鼠同源基因映射 (Biomart based)
def hsa2mmu(genelist:Sequence[str],drop:bool=False) -> Union[pd.DataFrame,List[str]]:
    q = Biomart().query(dataset='hsapiens_gene_ensembl',
        attributes=['ensembl_gene_id','external_gene_name',
                    'mmusculus_homolog_ensembl_gene',
                    'mmusculus_homolog_associated_gene_name'])
    if drop:
        return q.loc[q.external_gene_name.isin(genelist),'mmusculus_homolog_associated_gene_name'].dropna().tolist()
    else:
        return q.loc[q.external_gene_name.isin(genelist),:]

def mmu2hsa(genelist:Sequence[str],drop:bool=False) -> Union[pd.DataFrame,List[str]]:
    from gseapy import Biomart
    q = Biomart().query(dataset='hsapiens_gene_ensembl',
        attributes=['ensembl_gene_id','external_gene_name',
                    'hsapiens_homolog_ensembl_gene',
                    'hsapiens_homolog_associated_gene_name'])
    if drop:
        return q.loc[q.external_gene_name.isin(genelist),'hsapiens_homolog_associated_gene_name'].dropna().tolist()
    else:
        return q.loc[q.external_gene_name.isin(genelist),:]

# pd.set_option('display.max_columns', 10)
hsa2mmu(['IL1B','NLRP3','SPP1','CCL2'],drop=False)


### reference data
hs_moMac = sc.read_h5ad("/storage/data/KAI/Pan_Inflammation/Homo/Mac2025/analysis/CellPhoneDB/data/moMac.h5ad") # counts
hs_moMac.obs['species'] = 'human'
hs_moMac.write("/storage/data/KAI/Pan_Inflammation/Homo/Mac2025/analysis/CellPhoneDB/data/moMac.h5ad")

### query data
mm_moMac = sc.read_h5ad("/storage/data/KAI/Pan_Inflammation/Mouse/Pneumonia/LC/analysis/scArches/data/filtered_scvi/plmouse_moMac_filtered.h5ad") # counts
mm_moMac.obs['species'] = 'mouse'
mm_moMac.write("/storage/data/KAI/Pan_Inflammation/Mouse/Pneumonia/LC/analysis/scArches/data/filtered_scvi/plmouse_moMac_filtered.h5ad")

### 人鼠同源基因映射 (OrthoIntegrate based)
# 基于OrthoIntegrate包生成人-鼠symbol转换dataframe
Orthologue_DF = pd.read_csv("/storage/data/KAI/Pan_Inflammation/Mouse/Pneumonia/LC/analysis/OrthoIntegrate/data/hm2mm_Orthologue_DF.txt",sep='\t')

mouse_to_human = dict(zip(Orthologue_DF['mouse'], Orthologue_DF['human']))
mapped_mouse_genes = [gene for gene in mm_moMac.var_names if gene in mouse_to_human] # 16532
mm_moMac_Orth = mm_moMac[:, mapped_mouse_genes].copy()
mm_moMac_Orth.var['mouse_gene'] = mm_moMac_Orth.var_names
new_var_names = [mouse_to_human[gene] for gene in mm_moMac_Orth.var_names]
mm_moMac_Orth.var_names = new_var_names
mm_moMac_Orth.write("/storage/data/KAI/Pan_Inflammation/Mouse/Pneumonia/LC/analysis/scArches/data/scArches_scvi/input_data/mm_moMac_Orth.h5ad") # 共同同源基因数 16532

hs_moMac_Orth_common = hs_moMac[:, hs_moMac.var_names.isin(new_var_names)].copy()
hs_moMac_Orth_common.write("/storage/data/KAI/Pan_Inflammation/Mouse/Pneumonia/LC/analysis/scArches/data/scArches_scvi/input_data/hs_moMac_Orth.h5ad") # 共同同源基因数 16332
hs_moMac_Orth = sc.read_h5ad("/storage/data/KAI/Pan_Inflammation/Mouse/Pneumonia/LC/analysis/scArches/data/scArches_scvi/input_data/hs_moMac_Orth.h5ad")

# 人类-小鼠(16532-16332)同源基因数不同，过滤小鼠基因
mm_moMac_Orth_common = mm_moMac_Orth[:,mm_moMac_Orth.var_names.isin(hs_moMac_Orth_common.var_names)].copy() # 共同同源基因数 16332
# 检查人-鼠同源基因数是否相同
set(mm_moMac_Orth_common.var_names) == set(hs_moMac_Orth_common.var_names) # True

(mm_moMac_Orth_common.var_names == hs_moMac_Orth_common.var_names).all() # False 顺序不相同
mm_moMac_Orth_common.write("/storage/data/KAI/Pan_Inflammation/Mouse/Pneumonia/LC/analysis/scArches/data/scArches_scvi/input_data/mm_moMac_Orth.h5ad") # 共同同源基因数 16332


### scArches
# 训练scVI模型
sca.models.SCVI.setup_anndata(hs_moMac_Orth, batch_key='species', labels_key='cell_type_level3')
vae = sca.models.SCVI(
    hs_moMac_Orth,
    n_layers=2,
    encode_covariates=True,
    deeply_inject_covariates=False,
    use_layer_norm="both",
    use_batch_norm="none",
)
vae.train(max_epochs=400, early_stopping=True)
vae.save('/storage/data/KAI/Pan_Inflammation/Mouse/Pneumonia/LC/analysis/scArches/data/scArches_scvi/model/', overwrite=True, prefix="hs_moMac_scvi_")


# PLmouse_scvi_script.py
# nohup /storage/data/KAI/miniconda3/envs/XKscarches/bin/python /storage/data/KAI/Pan_Inflammation/PyProject/pythonProject/Python_inflammation/PLmouse_scvi_script.py > /storage/data/KAI/Pan_Inflammation/Mouse/Pneumonia/LC/analysis/scArches/data/scArches_scvi/scArches_scvi_log.txt 2>&1 &
# PID = 3398282
# Time 5.13 13:41



vae = sca.models.SCVI.load('/storage/data/KAI/Pan_Inflammation/Mouse/Pneumonia/LC/analysis/scArches/data/scArches_scvi/model/', adata=hs_moMac, prefix="hs_moMac_scvi_")
# 创建scANVI模型
scanvae = sca.models.SCANVI.from_scvi_model(vae, unlabeled_category = "Unknown")



#### macrophage level3预先注释 ----------------------------
mm_moMac = sc.read_h5ad("/storage/data/KAI/Pan_Inflammation/Mouse/Pneumonia/LC/analysis/scArches/data/filtered_scvi/plmouse_moMac_filtered.h5ad")
sc.pp.normalize_total(mm_moMac)
sc.pp.log1p(mm_moMac)
sc.pl.umap(mm_moMac, color=["leiden"],
           legend_loc='on data',
           # save="PLmouse_leiden.png"
           )

sc.pl.umap(mm_moMac,
           color=["Nlrp3"],
           # save="PLmouse_leiden.png",
           use_raw = False,
           )
sc.pl.umap(mm_moMac,
           color=["Il1b"],
           # save="PLmouse_leiden.png",
           use_raw = False,
           )
sc.pl.umap(mm_moMac,
           color=["Inhba"],
           # save="PLmouse_leiden.png",
           use_raw = False,
           )
Orthologue_DF[Orthologue_DF['human']=='MMP3']

# SPP1+CCL2+_Macro
sc.pl.umap(mm_moMac,
           color=["Spp1"],
           # save="PLmouse_leiden.png",
           use_raw = False,
           )
sc.pl.umap(mm_moMac,
           color=["Ccl2"],
           # save="PLmouse_leiden.png",
           use_raw = False,
           )
sc.pl.umap(mm_moMac,
           color=["Hif1a"],
           # save="PLmouse_leiden.png",
           use_raw = False,
           )

# CCL13+_Complement-associated_Macro
sc.pl.umap(mm_moMac,
           color=["C1qa"],
           # save="PLmouse_leiden.png",
           use_raw = False,
           )
sc.pl.umap(mm_moMac,
           color=["C1qb"],
           # save="PLmouse_leiden.png",
           use_raw = False,
           )
sc.pl.umap(mm_moMac,
           color=["C1qc"],
           # save="PLmouse_leiden.png",
           use_raw = False,
           )
sc.pl.umap(mm_moMac,
           color=["Ccl12"],
           # save="PLmouse_leiden.png",
           use_raw = False,
           )

# CD1C+_DC-like_Macro
sc.pl.umap(mm_moMac,
           color=["Csf1r"],
           # save="PLmouse_leiden.png",
           use_raw = False,
           )


# T cell attraction
sc.pl.umap(mm_moMac,
           color=["Cxcl9","Cxcl10","Ccl3","Ccl5","Ccl4","Ccl19","Ccl22","Ccl17",],
           # save="PLmouse_leiden.png",
           use_raw = False,
           )

sc.pl.umap(mm_moMac,
           color=['Plau','S100a10','Mafb'],
           # save="PLmouse_leiden.png",
           use_raw = False,
           )

sc.pl.umap(mm_moMac,
           color=['Cd163'],
           # save="PLmouse_leiden.png",
           use_raw = False,
           )








# 去除只有1个细胞的12簇（leiden 0.8）
mm_moMac = mm_moMac[mm_moMac.obs['leiden']!='12'].copy()
mm_moMac.write("/storage/data/KAI/Pan_Inflammation/Mouse/Pneumonia/LC/analysis/scArches/data/filtered_scvi/plmouse_moMac_filtered.h5ad")
sc.tl.rank_genes_groups(mm_moMac, groupby='leiden', method='wilcoxon')

result = mm_moMac.uns['rank_genes_groups']

clusters = result['names'].dtype.names
top_genes = {}
for cluster in clusters:
    # 提取每个聚类的 top 100 基因名
    top_genes[cluster] = result['names'][cluster][:100]
top_genes_df = pd.DataFrame(top_genes)
top_genes_df.to_csv('/storage/data/KAI/Pan_Inflammation/Mouse/Pneumonia/LC/analysis/annotation/momac_level3/moMac_clustertop_100.txt', sep='\t', index=True)

cluster_1_genes = pd.DataFrame({
    'gene': result['names']['1'],
    'logfoldchanges': result['logfoldchanges']['1'],
    'pvals': result['pvals']['1'],
    'pvals_adj': result['pvals_adj']['1']
})
cluster_1_genes['gene'].head(100).to_list()










