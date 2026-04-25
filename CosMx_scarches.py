#### package -------------------------------------
from pathlib import Path
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import scarches as sca
from scarches.dataset.trvae.data_handling import remove_sparsity
import matplotlib.pyplot as plt

#### load data -------------------------------------
source_adata = sc.read_h5ad("/storage/data/KAI/Pan_Inflammation/STRNA/adata_ref/Colorectum/scarches/data/source_adata.h5ad")
source_adata.obs['cell_type_TL'] = source_adata.obs['cell_type_TL'].cat.add_categories(['Plasma cell'])
source_adata.obs.loc[
    source_adata.obs['cell_type_level3'].isin(['Plasma', 'Plasmablast']),
    'cell_type_TL'
] = 'Plasma cell'

source_adata.write("/storage/data/KAI/Pan_Inflammation/STRNA/adata_ref/Colorectum/scarches/data/source_adata.h5ad")

source_adata = sc.read_h5ad("/storage/data/KAI/Pan_Inflammation/STRNA/adata_ref/Colorectum/scarches/data/source_adata.h5ad")
target_adata = sc.read_h5ad("/storage/data/KAI/Pan_Inflammation/STRNA/adata_ref/Colorectum/scarches/data/target_adata.h5ad")

### scArches (all_cell_type)

sca.models.SCVI.setup_anndata(source_adata, batch_key='Omics', labels_key='cell_type_TL')
vae = sca.models.SCVI(
    source_adata,
    n_layers=2,
    encode_covariates=True,
    deeply_inject_covariates=False,
    use_layer_norm="both",
    use_batch_norm="none",
)
vae.train(max_epochs=400, early_stopping=True)

vae.save('/storage/data/KAI/Pan_Inflammation/STRNA/adata_ref/Colorectum/scarches/model/scvi/', overwrite=True, prefix="colonref_scvi_",)
vae = sca.models.SCVI.load("/storage/data/KAI/Pan_Inflammation/STRNA/adata_ref/Colorectum/scarches/model/scvi/",
                           adata=source_adata,
                           prefix="colonref_scvi_",
                           )

# model train (scANVI)
scanvae = sca.models.SCANVI.from_scvi_model(vae, unlabeled_category = "Unknown")
scanvae.train(max_epochs=20)
scanvae.save('/storage/data/KAI/Pan_Inflammation/STRNA/adata_ref/Colorectum/scarches/model/', overwrite=True,prefix="colonref_scanvi_" )
scanvae = sca.models.SCANVI.load("/storage/data/KAI/Pan_Inflammation/STRNA/adata_ref/Colorectum/scarches/model/",
                                 adata=source_adata,
                                 prefix="colonref_scanvi_",
                                 )

# Create anndata file of latent representation and compute UMAP
reference_latent = sc.AnnData(scanvae.get_latent_representation())
reference_latent.obs["cell_type"] = source_adata.obs['cell_type_TL'].tolist()
reference_latent.obs["batch"] = source_adata.obs['Omics'].tolist()
reference_latent.obs['predictions'] = scanvae.predict()
print("Acc: {}".format(np.mean(reference_latent.obs.predictions == reference_latent.obs.cell_type)))

# Surgery model
model = sca.models.SCANVI.load_query_data(
    target_adata,
    reference_model="/storage/data/KAI/Pan_Inflammation/STRNA/adata_ref/Colorectum/scarches/model/scanvi/",
    freeze_dropout = False,           # 允许 Dropout 层参与训练
    freeze_batchnorm_encoder = False, # 放开编码器 BatchNorm 训练
    freeze_classifier = False         # 解冻分类器层
)

model.train(
    max_epochs=400,
    plan_kwargs=dict(weight_decay=0.0),
    check_val_every_n_epoch=10,
)

model.save('/storage/data/KAI/Pan_Inflammation/STRNA/adata_ref/Colorectum/scarches/model/surgery/epoch400/',overwrite=True)

query_latent = sc.AnnData(model.get_latent_representation())
query_latent.obs['cell_type'] = target_adata.obs['cell_type_TL'].tolist()
query_latent.obs['batch'] = target_adata.obs['Omics'].tolist()
query_latent.obs['predictions'] = model.predict()


## 参考原有meta信息中的细胞大类注释信息，计算模型不同迭代次数的迁移标签准确率
# Epoch = 400
target_adata.obs['cell_type_TL400'] = query_latent.obs['predictions'].values

sc.pl.umap(target_adata,
           color='cell_type_TL400',
           use_raw = False,
           legend_loc='right margin',
           legend_fontsize='small',
           show=False,
           size=3
           )
plt.gcf().set_size_inches(12, 6)  # 调整画布尺寸，确保 legend 有空间
plt.tight_layout()  # 自动调整布局避免裁剪
plt.show()

# marker visual


sc.pl.umap(target_adata,
           color=['EPCAM','RGS5','DCN','PECAM1','CD68','CD3E','CD79A','NRXN1'],
           use_raw = False,
           # legend_loc='right margin',
           legend_fontsize='small',
           show=False,
           size=3
           )
plt.gcf().set_size_inches(12, 6)  # 调整画布尺寸，确保 legend 有空间
plt.tight_layout()  # 自动调整布局避免裁剪
plt.show()

sc.pl.umap(target_adata,
           color=['JCHAIN','FKBP11'],
           use_raw = False,
           # legend_loc='right margin',
           legend_fontsize='small',
           show=True,
           size=3
           )


target_adata.write('/storage/data/KAI/Pan_Inflammation/STRNA/adata_ref/Colorectum/scarches/data/target_adata.h5ad')

### scArches (only for Macrophage)

## 导入原文注释信息观察
cosmx_meta = pd.read_csv("/storage/data/KAI/Pan_Inflammation/STRNA/adata_ref/Colorectum/scarches/data/meta/GSE234713_CosMx_annotation.csv")
adata = sc.read_h5ad('/storage/data/KAI/Pan_Inflammation/STRNA/adata_ref/Colorectum/scarches/data/target_adata.h5ad')

# 修改adata.obs index实现匹配
import re
def insert_underscore(index_str):
    return re.sub(r'^([A-Z]{2})([a-z])', r'\1_\2', index_str)
adata.obs["matching_id"] = adata.obs.index.to_series().apply(insert_underscore)

cosmx_meta = cosmx_meta.set_index("id")
adata.obs = adata.obs.join(cosmx_meta, on="matching_id")

sc.pl.umap(adata,
           color='subset',
           # use_raw = False,
           legend_loc='right margin',
           legend_fontsize='small',
           show=False,
           size=3
           )
plt.gcf().set_size_inches(8, 6)  # 调整画布尺寸，确保 legend 有空间
plt.tight_layout()  # 自动调整布局避免裁剪
plt.show()
sc.pl.umap(adata,
           color='SingleR2',
           # use_raw = False,
           legend_loc='right margin',
           legend_fontsize='small',
           show=False,
           size=3
           )
plt.gcf().set_size_inches(12, 6)  # 调整画布尺寸，确保 legend 有空间
plt.tight_layout()  # 自动调整布局避免裁剪
plt.show()

unmatched_ids = set(adata.obs['matching_id']) - set(cosmx_meta.index) # 不是所有细胞都能在cosmx_meta中匹配上
adata.obs['matching_id'].isin(unmatched_ids).sum()
adata.obs['subset'].isna().sum() # cosmx_meta中本身也存在注释信息为NA的细胞

# 删除原文中细胞注释类型为空的细胞
adata = adata[~adata.obs['subset'].isna()].copy() # 458983,966
adata.write("/storage/data/KAI/Pan_Inflammation/STRNA/adata_ref/Colorectum/scarches/data/annotation/GSE234713_anno.h5ad")

## 细胞大类注释
sc.pl.umap(adata,
           color='subset',
           # use_raw = False,
           legend_loc='right margin',
           legend_fontsize='small',
           show=False,
           size=3
           )
plt.gcf().set_size_inches(8, 6)  # 调整画布尺寸，确保 legend 有空间
plt.tight_layout()  # 自动调整布局避免裁剪
plt.show()

## 提取髓系细胞进行再分型
adata_mye = adata[adata.obs['subset']=='myeloids'].copy()
import scvi
import torch

scvi.model.SCVI.setup_anndata(adata_mye, layer="counts", categorical_covariate_keys=["batch_fov"])
vae = scvi.model.SCVI(adata_mye,n_layers=2, n_latent=30,gene_likelihood="nb")
vae.train()
model_dir = os.path.join("/storage/data/KAI/Pan_Inflammation/STRNA/adata_ref/Colorectum/scarches/data/annotation/mye/scvi_model/", "colonmye_scvi_model")
vae.save(model_dir, overwrite=True, prefix="colonmye_scvi_")

adata_mye.obsm["X_scVI"] = vae.get_latent_representation()
adata_mye.layers["X_normalized_scVI"] = vae.get_normalized_expression(library_size=10e4)

sc.pp.neighbors(adata_mye,use_rep="X_scVI")
sc.tl.umap(adata_mye,min_dist=0.3)

sc.pl.umap(adata_mye,
           color=["Sample_geo_accession","Inflammation","batch_fov"],
           ncols=3,
           frameon=False
           )

adata_mye.write("/storage/data/KAI/Pan_Inflammation/STRNA/adata_ref/Colorectum/scarches/data/annotation/mye/colonmye_scvi.h5ad")

## 准备colon myeloid参考单细胞数据
source_adata = sc.read_h5ad("/storage/data/KAI/Pan_Inflammation/STRNA/adata_ref/Colorectum/scarches/data/source_adata.h5ad")
source_mye = source_adata[source_adata.obs['cell_type_TL'].isin(['Monocyte','CD1C+_DC-like_Macro',
                                                                 'CCL13+_Complement-associated_Macro',
                                                                 'IL1B+NLRP3+_Macro','SPP1+CCL2+_Macro',
                                                                 'MMP3+CXCL8+_Macro','Intestinal resident Macro',
                                                                 'Dendritic cell','Mast cell'
                                                                 ])].copy()


sca.models.SCVI.setup_anndata(source_mye, batch_key='Omics', labels_key='cell_type_TL')
vae = sca.models.SCVI(
    source_mye,
    n_layers=2,
    encode_covariates=True,
    deeply_inject_covariates=False,
    use_layer_norm="both",
    use_batch_norm="none",
)
vae.train(max_epochs=400, early_stopping=True)

vae.save('/storage/data/KAI/Pan_Inflammation/STRNA/adata_ref/Colorectum/scarches/data/annotation/mye/ref_model/scvi/', overwrite=True, )
vae = sca.models.SCVI.load("/storage/data/KAI/Pan_Inflammation/STRNA/adata_ref/Colorectum/scarches/data/annotation/mye/ref_model/scvi/",
                           adata=source_mye,
                           )
# model train (scANVI)
scanvae = sca.models.SCANVI.from_scvi_model(vae, unlabeled_category = "Unknown")
scanvae.train(max_epochs=20)
scanvae.save('/storage/data/KAI/Pan_Inflammation/STRNA/adata_ref/Colorectum/scarches/data/annotation/mye/ref_model/scanvi/', overwrite=True, )
scanvae = sca.models.SCANVI.load("/storage/data/KAI/Pan_Inflammation/STRNA/adata_ref/Colorectum/scarches/data/annotation/mye/ref_model/scanvi/",
                                 adata=source_mye,
                                 )

reference_latent = sc.AnnData(scanvae.get_latent_representation())
reference_latent.obs["cell_type"] = source_mye.obs['cell_type_TL'].tolist()
reference_latent.obs["batch"] = source_mye.obs['Omics'].tolist()

reference_latent.obs['predictions'] = scanvae.predict()
print("Acc: {}".format(np.mean(reference_latent.obs.predictions == reference_latent.obs.cell_type)))

# epoch = 400
model = sca.models.SCANVI.load_query_data(
    adata_mye,
    "/storage/data/KAI/Pan_Inflammation/STRNA/adata_ref/Colorectum/scarches/data/annotation/mye/ref_model/scanvi/",
    freeze_dropout = False,           # 允许 Dropout 层参与训练
    freeze_batchnorm_encoder = False, # 放开编码器 BatchNorm 训练
    freeze_classifier = False         # 解冻分类器层
)

model.train(
    max_epochs=400,
    plan_kwargs=dict(weight_decay=0.0),
    check_val_every_n_epoch=10,
)

model.save('/storage/data/KAI/Pan_Inflammation/STRNA/adata_ref/Colorectum/scarches/data/annotation/mye/sur_model/epoch400/',overwrite=True)

query_latent = sc.AnnData(model.get_latent_representation())
query_latent.obs['cell_type'] = adata_mye.obs['cell_type_TL'].tolist()
query_latent.obs['batch'] = adata_mye.obs['Omics'].tolist()
query_latent.obs['predictions'] = model.predict()

query_latent.obs['predictions'].value_counts()

adata_mye.obs['cell_type_TL400'] = query_latent.obs['predictions'].values

adata_mye.write("/storage/data/KAI/Pan_Inflammation/STRNA/adata_ref/Colorectum/scarches/data/annotation/mye/colonmye_scvi.h5ad")

sc.pl.umap(adata_mye,
           color='cell_type_TL400',
           # use_raw = False,
           legend_loc='right margin',
           legend_fontsize='small',
           show=False,
           size=3
           )
plt.gcf().set_size_inches(8, 6)  # 调整画布尺寸，确保 legend 有空间
plt.tight_layout()  # 自动调整布局避免裁剪
plt.show()

# epoch = 100
model = sca.models.SCANVI.load_query_data(
    adata_mye,
    "/storage/data/KAI/Pan_Inflammation/STRNA/adata_ref/Colorectum/scarches/data/annotation/mye/ref_model/scanvi/",
    freeze_dropout = False,           # 允许 Dropout 层参与训练
    freeze_batchnorm_encoder = False, # 放开编码器 BatchNorm 训练
    freeze_classifier = False         # 解冻分类器层
)

model.train(
    max_epochs=100,
    plan_kwargs=dict(weight_decay=0.0),
    check_val_every_n_epoch=10,
)

model.save('/storage/data/KAI/Pan_Inflammation/STRNA/adata_ref/Colorectum/scarches/data/annotation/mye/sur_model/epoch100/',overwrite=True)
query_latent = sc.AnnData(model.get_latent_representation())
query_latent.obs['cell_type'] = adata_mye.obs['cell_type_TL'].tolist()
query_latent.obs['batch'] = adata_mye.obs['Omics'].tolist()
query_latent.obs['predictions'] = model.predict()

query_latent.obs['predictions'].value_counts()

adata_mye.obs['cell_type_TL100'] = query_latent.obs['predictions'].values

adata_mye.write("/storage/data/KAI/Pan_Inflammation/STRNA/adata_ref/Colorectum/scarches/data/annotation/mye/colonmye_scvi.h5ad")


# epoch = 1000
model = sca.models.SCANVI.load_query_data(
    adata_mye,
    "/storage/data/KAI/Pan_Inflammation/STRNA/adata_ref/Colorectum/scarches/data/annotation/mye/ref_model/scanvi/",
    freeze_dropout = False,           # 允许 Dropout 层参与训练
    freeze_batchnorm_encoder = False, # 放开编码器 BatchNorm 训练
    freeze_classifier = False         # 解冻分类器层
)

model.train(
    max_epochs=1000,
    plan_kwargs=dict(weight_decay=0.0),
    check_val_every_n_epoch=10,
)

model.save('/storage/data/KAI/Pan_Inflammation/STRNA/adata_ref/Colorectum/scarches/data/annotation/mye/sur_model/epoch1000/',overwrite=True)

query_latent = sc.AnnData(model.get_latent_representation())
query_latent.obs['cell_type'] = adata_mye.obs['cell_type_TL'].tolist()
query_latent.obs['batch'] = adata_mye.obs['Omics'].tolist()
query_latent.obs['predictions'] = model.predict()

query_latent.obs['predictions'].value_counts()

adata_mye.obs['cell_type_TL1000'] = query_latent.obs['predictions'].values

sc.pl.umap(adata_mye,
           color='cell_type_TL1000',
           # use_raw = False,
           legend_loc='right margin',
           legend_fontsize='small',
           show=False,
           size=3
           )
plt.gcf().set_size_inches(8, 6)  # 调整画布尺寸，确保 legend 有空间
plt.tight_layout()  # 自动调整布局避免裁剪
plt.show()



adata_mye.write("/storage/data/KAI/Pan_Inflammation/STRNA/adata_ref/Colorectum/scarches/data/annotation/mye/colonmye_scvi.h5ad")


# epoch = 2000
model = sca.models.SCANVI.load_query_data(
    adata_mye,
    "/storage/data/KAI/Pan_Inflammation/STRNA/adata_ref/Colorectum/scarches/data/annotation/mye/ref_model/scanvi/",
    freeze_dropout = False,           # 允许 Dropout 层参与训练
    freeze_batchnorm_encoder = False, # 放开编码器 BatchNorm 训练
    freeze_classifier = False         # 解冻分类器层
)

model.train(
    max_epochs=2000,
    plan_kwargs=dict(weight_decay=0.0),
    check_val_every_n_epoch=10,
)

model.save('/storage/data/KAI/Pan_Inflammation/STRNA/adata_ref/Colorectum/scarches/data/annotation/mye/sur_model/epoch2000/',overwrite=True)


query_latent = sc.AnnData(model.get_latent_representation())
query_latent.obs['cell_type'] = adata_mye.obs['cell_type_TL'].tolist()
query_latent.obs['batch'] = adata_mye.obs['Omics'].tolist()
query_latent.obs['predictions'] = model.predict()

query_latent.obs['predictions'].value_counts()
adata_mye.obs['cell_type_TL2000'] = query_latent.obs['predictions'].values

sc.pl.umap(adata_mye,
           color='cell_type_TL2000',
           # use_raw = False,
           legend_loc='right margin',
           legend_fontsize='small',
           show=False,
           size=5
           )
plt.gcf().set_size_inches(8, 6)  # 调整画布尺寸，确保 legend 有空间
plt.tight_layout()  # 自动调整布局避免裁剪
plt.show()

sc.pl.umap(adata_mye,
           color='SingleR2',
           # use_raw = False,
           legend_loc='right margin',
           legend_fontsize='small',
           show=False,
           size=5
           )
plt.gcf().set_size_inches(8, 6)  # 调整画布尺寸，确保 legend 有空间
plt.tight_layout()  # 自动调整布局避免裁剪
plt.show()

adata_mye.X = adata_mye.layers['counts'].copy()

adata_mye.write("/storage/data/KAI/Pan_Inflammation/STRNA/adata_ref/Colorectum/scarches/data/annotation/mye/colonmye_scvi.h5ad")


adata_mye.X = adata_mye.layers['counts'].copy()
sc.pp.normalize_total(adata_mye, inplace=True, target_sum=1e4)
sc.pp.log1p(adata_mye)

sc.pl.umap(adata_mye,
           color=['IL1B','NLRP3','SPP1','CCL2'],
           use_raw = False,
           legend_loc='right margin',
           legend_fontsize='small',
           show=True,
           size=5
           )
plt.gcf().set_size_inches(8, 6)  # 调整画布尺寸，确保 legend 有空间
plt.tight_layout()  # 自动调整布局避免裁剪
plt.show()



### scArches (only for Macrophage in myeloids)
## construct source_adata
source_adata = sc.read_h5ad("/storage/data/KAI/Pan_Inflammation/STRNA/adata_ref/Colorectum/scarches/data/source_adata.h5ad")
source_macro = source_adata[source_adata.obs['cell_type_level3'].isin(['CD1C+_DC-like_Macro',
                                                                       'CCL13+_Complement-associated_Macro',
                                                                       'IL1B+NLRP3+_Macro',
                                                                       'SPP1+CCL2+_Macro',
                                                                       'MMP3+CXCL8+_Macro',
                                                                       'Intestinal resident Macro'
                                                                       ])].copy()
source_macro.obs['cell_type_TL'] = source_macro.obs['cell_type_level3'].copy()
source_macro.write("/storage/data/KAI/Pan_Inflammation/STRNA/adata_ref/Colorectum/scarches/data/annotation/macro/data/source_macro.h5ad")
## construct target_adata
target_adata = sc.read_h5ad("/storage/data/KAI/Pan_Inflammation/STRNA/adata_ref/Colorectum/scarches/data/annotation/mye/colonmye_scvi.h5ad")
target_macro = target_adata[target_adata.obs['SingleR2'].isin(['Macrophage NRG1','M0','M1','M2'])].copy()

target_macro.write("/storage/data/KAI/Pan_Inflammation/STRNA/adata_ref/Colorectum/scarches/data/annotation/macro/data/target_macro.h5ad")


## scVI model (ref_model)
sca.models.SCVI.setup_anndata(source_macro, batch_key='Omics', labels_key='cell_type_TL')
vae = sca.models.SCVI(
    source_macro,
    n_layers=2,
    encode_covariates=True,
    deeply_inject_covariates=False,
    use_layer_norm="both",
    use_batch_norm="none",
)
vae.train(max_epochs=400, early_stopping=True)
vae.save('/storage/data/KAI/Pan_Inflammation/STRNA/adata_ref/Colorectum/scarches/data/annotation/macro/ref_model/scvi/', overwrite=True, )

vae = sca.models.SCVI.load("/storage/data/KAI/Pan_Inflammation/STRNA/adata_ref/Colorectum/scarches/data/annotation/macro/ref_model/scvi/",
                           adata=source_macro,
                           )
# scANVI model (ref_model)
scanvae = sca.models.SCANVI.from_scvi_model(vae, unlabeled_category = "Unknown")
scanvae.train(max_epochs=20)
scanvae.save('/storage/data/KAI/Pan_Inflammation/STRNA/adata_ref/Colorectum/scarches/data/annotation/macro/ref_model/scanvi/', overwrite=True, )
scanvae = sca.models.SCANVI.load("/storage/data/KAI/Pan_Inflammation/STRNA/adata_ref/Colorectum/scarches/data/annotation/macro/ref_model/scanvi/",
                                 adata=source_macro,
                                 )

reference_latent = sc.AnnData(scanvae.get_latent_representation())
reference_latent.obs["cell_type"] = source_macro.obs['cell_type_TL'].tolist()
reference_latent.obs["batch"] = source_macro.obs['Omics'].tolist()

reference_latent.obs['predictions'] = scanvae.predict()
print("Acc: {}".format(np.mean(reference_latent.obs.predictions == reference_latent.obs.cell_type)))

# epoch = 100
model = sca.models.SCANVI.load_query_data(
    target_macro,
    "/storage/data/KAI/Pan_Inflammation/STRNA/adata_ref/Colorectum/scarches/data/annotation/macro/ref_model/scanvi/",
    freeze_dropout = False,           # 允许 Dropout 层参与训练
    freeze_batchnorm_encoder = False, # 放开编码器 BatchNorm 训练
    freeze_classifier = False         # 解冻分类器层
)

model.train(
    max_epochs=100,
    plan_kwargs=dict(weight_decay=0.0),
    check_val_every_n_epoch=10,
)

model.save('/storage/data/KAI/Pan_Inflammation/STRNA/adata_ref/Colorectum/scarches/data/annotation/macro/surg_model/epoch100/',overwrite=True)
query_latent = sc.AnnData(model.get_latent_representation())
query_latent.obs['cell_type'] = target_macro.obs['cell_type_TL'].tolist()
query_latent.obs['batch'] = target_macro.obs['Omics'].tolist()
query_latent.obs['predictions'] = model.predict()

query_latent.obs['predictions'].value_counts()

target_macro.obs['cell_type_TL100'] = query_latent.obs['predictions'].values

target_macro.write("/storage/data/KAI/Pan_Inflammation/STRNA/adata_ref/Colorectum/scarches/data/annotation/macro/data/colonmacro_scvi.h5ad")

# epoch = 2000
model = sca.models.SCANVI.load_query_data(
    target_macro,
    "/storage/data/KAI/Pan_Inflammation/STRNA/adata_ref/Colorectum/scarches/data/annotation/macro/ref_model/scanvi/",
    freeze_dropout = False,           # 允许 Dropout 层参与训练
    freeze_batchnorm_encoder = False, # 放开编码器 BatchNorm 训练
    freeze_classifier = False         # 解冻分类器层
)

model.train(
    max_epochs=2000,
    plan_kwargs=dict(weight_decay=0.0),
    check_val_every_n_epoch=10,
)

model.save('/storage/data/KAI/Pan_Inflammation/STRNA/adata_ref/Colorectum/scarches/data/annotation/macro/surg_model/epoch2000/',overwrite=True)
query_latent = sc.AnnData(model.get_latent_representation())
query_latent.obs['cell_type'] = target_macro.obs['cell_type_TL'].tolist()
query_latent.obs['batch'] = target_macro.obs['Omics'].tolist()
query_latent.obs['predictions'] = model.predict()

query_latent.obs['predictions'].value_counts()



### scArches (only for Macrophage, IL1B/SPP1/Other)
source_macro = sc.read_h5ad("/storage/data/KAI/Pan_Inflammation/STRNA/adata_ref/Colorectum/scarches/data/annotation/macro/data/source_macro.h5ad")
target_macro = sc.read_h5ad("/storage/data/KAI/Pan_Inflammation/STRNA/adata_ref/Colorectum/scarches/data/annotation/macro/data/target_macro.h5ad")

## 调整source数据巨噬细胞标签
source_macro.obs['cell_type_level3'] = source_macro.obs['cell_type_level3'].cat.add_categories('Other Macro')
source_macro.obs['cell_scarches'] = source_macro.obs['cell_type_level3'].where(
    source_macro.obs['cell_type_level3'].isin(['IL1B+NLRP3+_Macro', 'SPP1+CCL2+_Macro']),
    'Other Macro'
)
# 将未使用的 category 从分类列表中移除
source_macro.obs['cell_type_level3'] = source_macro.obs['cell_type_level3'].cat.remove_unused_categories()
source_macro.obs['cell_scarches'] = source_macro.obs['cell_scarches'].cat.remove_unused_categories()

source_macro.obs['cell_scarches'].value_counts()
# Other Macro          15960
# IL1B+NLRP3+_Macro     5951
# SPP1+CCL2+_Macro      2452
# Name: cell_scarches, dtype: int64

## scVI model (ref_model)
sca.models.SCVI.setup_anndata(source_macro, batch_key='Omics', labels_key='cell_scarches')
vae = sca.models.SCVI(
    source_macro,
    n_layers=2,
    encode_covariates=True,
    deeply_inject_covariates=False,
    use_layer_norm="both",
    use_batch_norm="none",
)
vae.train(max_epochs=400, early_stopping=True)
vae.save('/storage/data/KAI/Pan_Inflammation/STRNA/adata_ref/Colorectum/scarches/data/annotation/macro/new_ref_model/scvi/', overwrite=True, )

vae = sca.models.SCVI.load("/storage/data/KAI/Pan_Inflammation/STRNA/adata_ref/Colorectum/scarches/data/annotation/macro/new_ref_model/scvi/",
                           adata=source_macro,
                           )
# scANVI model (ref_model)
scanvae = sca.models.SCANVI.from_scvi_model(vae, unlabeled_category = "Unknown")
scanvae.train(max_epochs=20)
scanvae.save('/storage/data/KAI/Pan_Inflammation/STRNA/adata_ref/Colorectum/scarches/data/annotation/macro/new_ref_model/scanvi/', overwrite=True, )
scanvae = sca.models.SCANVI.load("/storage/data/KAI/Pan_Inflammation/STRNA/adata_ref/Colorectum/scarches/data/annotation/macro/new_ref_model/scanvi/",
                                 adata=source_macro,
                                 )
reference_latent = sc.AnnData(scanvae.get_latent_representation())
reference_latent.obs["cell_type"] = source_macro.obs['cell_scarches'].tolist()
reference_latent.obs["batch"] = source_macro.obs['Omics'].tolist()

reference_latent.obs['predictions'] = scanvae.predict()
print("Acc: {}".format(np.mean(reference_latent.obs.predictions == reference_latent.obs.cell_type)))
# Acc: 0.9848130361613923

model = sca.models.SCANVI.load_query_data(
    target_macro,
    "/storage/data/KAI/Pan_Inflammation/STRNA/adata_ref/Colorectum/scarches/data/annotation/macro/new_ref_model/scanvi/",
    freeze_dropout = False,           # 允许 Dropout 层参与训练
    freeze_batchnorm_encoder = False, # 放开编码器 BatchNorm 训练
    freeze_classifier = False         # 解冻分类器层
)
# epoch = 400
weights = {
    "IL1B+NLRP3+_Macro": 5.0,
    "SPP1+CCL2+_Macro": 1.0,
    "Other Macro": 1.0
}

model.train(
    max_epochs=400,
    plan_kwargs=dict(weight_decay=0.0),
    check_val_every_n_epoch=10,
)
model.save('/storage/data/KAI/Pan_Inflammation/STRNA/adata_ref/Colorectum/scarches/data/annotation/macro/new_surg_model/epoch400/',overwrite=True)
query_latent = sc.AnnData(model.get_latent_representation())
query_latent.obs['cell_type'] = target_macro.obs['cell_type_TL'].tolist()
query_latent.obs['batch'] = target_macro.obs['Omics'].tolist()
query_latent.obs['predictions'] = model.predict()


query_latent.obs['predictions'].value_counts()

target_macro.obs['cell_arches'] = query_latent.obs['predictions'].tolist()
sc.pl.violin(target_macro, keys=["SPP1", "CCL2"], groupby="cell_arches", stripplot=False)



# epoch = 100
model.train(
    max_epochs=100,
    plan_kwargs=dict(weight_decay=0.0),
    check_val_every_n_epoch=10,
    early_stopping=True,
)
model.save('/storage/data/KAI/Pan_Inflammation/STRNA/adata_ref/Colorectum/scarches/data/annotation/macro/new_surg_model/epoch100/',overwrite=True)
query_latent = sc.AnnData(model.get_latent_representation())
query_latent.obs['cell_type'] = target_macro.obs['cell_type_TL'].tolist()
query_latent.obs['batch'] = target_macro.obs['Omics'].tolist()
query_latent.obs['predictions'] = model.predict()

query_latent.obs['predictions'].value_counts()

target_macro.obs['cell_arches'] = query_latent.obs['predictions'].tolist()

target_macro.write("/storage/data/KAI/Pan_Inflammation/STRNA/adata_ref/Colorectum/scarches/data/annotation/macro/data/target_macro.h5ad")

sc.pp.normalize_total(target_macro, inplace=True, target_sum=1e4)
sc.pp.log1p(target_macro)


sc.pl.violin(target_macro, keys=["SPP1", "CCL2"], groupby="cell_arches", stripplot=False,use_raw=False)
sc.pl.violin(target_macro, keys=["IL1B", "NLRP3"], groupby="cell_arches", stripplot=False,use_raw=False)


