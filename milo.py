#### packages ---------------------------
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import scipy

import milopy
import milopy.core as milo
import milopy.plot as milopl

import pertpy

#### load anndata ------------------------
Mac2025 = sc.read_h5ad("/storage/data/KAI/Pan_Inflammation/Homo/Mac2025/scVI/new/pan_Mac_counts.h5ad")

#### build kNN graph ---------------------
### run PCA
sc.tl.pca(Mac2025, n_comps=50)
sc.pp.neighbors(Mac2025, n_neighbors=30, n_pcs=30)
### kNN graph
milo.make_nhoods(Mac2025, prop=0.1)


#### count cells -------------------
# Mac2025 = sc.read_h5ad("/storage/data/KAI/Pan_Inflammation/Homo/Mac2025/analysis/milopy/data/allMac2025_milopy.h5ad")

Mac2025.obsm["nhoods"]

Mac2025[Mac2025.obs['nhood_ixs_refined'] != 0].obs[['nhood_ixs_refined', 'nhood_kth_distance']]
nhood_size = np.array(Mac2025.obsm["nhoods"].sum(0)).ravel()
plt.hist(nhood_size, bins=100);

milo.count_nhoods(Mac2025, sample_col="Sample_geo_accession")
Mac2025.uns["nhood_adata"]

#### Differential abundance ---------------------------------
Mac2025.obs["Inflammation_continuous"] = Mac2025.obs["Inflammation"].cat.codes
milo.DA_nhoods(Mac2025, design="~Inflammation_continuous")

Mac2025.uns["nhood_adata"].obs

### spatialFDR plot
old_figsize = plt.rcParams["figure.figsize"]
plt.rcParams["figure.figsize"] = [10,5]
plt.subplot(1,2,1)
plt.hist(Mac2025.uns["nhood_adata"].obs.PValue, bins=50);
plt.xlabel("P-Vals");
plt.subplot(1,2,2)
plt.plot(Mac2025.uns["nhood_adata"].obs.logFC, -np.log10(Mac2025.uns["nhood_adata"].obs.SpatialFDR), '.');
plt.xlabel("log-Fold Change");
plt.ylabel("- log10(Spatial FDR)");
plt.tight_layout()
plt.rcParams["figure.figsize"] = old_figsize


####
import milopy.utils
milopy.utils.build_nhood_graph(Mac2025)
plt.rcParams["figure.figsize"] = [10,10]
milopl.plot_nhood_graph(Mac2025,
                        alpha=0.01, ## SpatialFDR level (1%)
                        min_size=2 ## Size of smallest dot
                       )

milopy.plot.plot_nhood_graph(Mac2025, alpha=0.2, min_size=5)

milopy.utils.annotate_nhoods(Mac2025, anno_col='cell_type_level3')

### Mixed identify
# Mac2025.uns['nhood_adata'].obs.loc[Mac2025.uns['nhood_adata'].obs["nhood_annotation_frac"] < 0.6, "nhood_annotation"] = "Mixed"

nhood_adata = Mac2025.uns["nhood_adata"]
obs = nhood_adata.obs

# 确保该列存在且是分类类型；若不是，则转为分类
if "nhood_annotation" not in obs.columns:
    raise KeyError("obs 中没有列 'nhood_annotation'")

if not pd.api.types.is_categorical_dtype(obs["nhood_annotation"]):
    obs["nhood_annotation"] = obs["nhood_annotation"].astype("category")

# 将新类别加入分类再赋值
obs["nhood_annotation"] = obs["nhood_annotation"].cat.add_categories(["Mixed"])
mask = obs["nhood_annotation_frac"] < 0.6
obs.loc[mask, "nhood_annotation"] = "Mixed"

# 移除未使用类别（可选）
obs["nhood_annotation"] = obs["nhood_annotation"].cat.remove_unused_categories()

# 回写到 Mac2025（稳妥）
nhood_adata.obs = obs
Mac2025.uns["nhood_adata"] = nhood_adata






sc.pl.violin(Mac2025.uns['nhood_adata'], "logFC", groupby="nhood_annotation", rotation=90, show=False);
plt.axhline(y=0, color='black', linestyle='--');
plt.show()




# Mac2025.write("/storage/data/KAI/Pan_Inflammation/Homo/Mac2025/analysis/milopy/data/allMac2025_milopy.h5ad")


##### save data to R --------------------------------------

nh = Mac2025.uns["nhood_adata"]     # 邻域级 AnnData
Nhood = np.array(nh.obs_names)      # 邻域ID（主键）
Cells = np.array(Mac2025.obs_names) # 细胞条码（与 R 侧 SCE 匹配）

# 1) 邻域DA结果（给 plotDAbeeswarm 用，也会给 plotNhoodGraphDA 着色）
da_cols = ["logFC","logCPM","F","PValue","FDR","SpatialFDR",
           "Nhood_size","nhood_annotation","nhood_annotation_frac"]
da_df = nh.obs[da_cols].copy()
da_df.insert(0, "Nhood", Nhood)
da_df.to_csv("/storage/data/KAI/Pan_Inflammation/Homo/Mac2025/analysis/milopy/data/nhood_da.csv", index=False)


# 2) 邻域二维坐标（milopy的图坐标，一般是力导FR布局；列名用 x,y）
coords = pd.DataFrame(nh.obsm["X_milo_graph"], columns=["x","y"])
coords.insert(0, "Nhood", Nhood)
coords.to_csv("/storage/data/KAI/Pan_Inflammation/Homo/Mac2025/analysis/milopy/data/nhood_coords.csv", index=False)

# 3) 邻域图（边）——从稀疏邻接转边表
from scipy import sparse
from scipy.io import mmwrite
# 确保 A 是 CSR 或 CSC 类型
A = nh.obsp["nhood_connectivities"]
# 取上三角，避免重复边
A = sparse.triu(A, k=1).tocsr()   # ✅ 关键改动：.tocsr()
# 获取非零边（source-target）
src, dst = A.nonzero()
# 提取权重（连接强度）
w = np.asarray(A[src, dst]).ravel()
# 转成 DataFrame
Nhood = np.array(nh.obs_names)
edges = pd.DataFrame({
    "source": Nhood[src],
    "target": Nhood[dst],
    "weight": w
})
edges.to_csv("/storage/data/KAI/Pan_Inflammation/Homo/Mac2025/analysis/milopy/data/nhood_edges.csv", index=False)

# 4) 邻域成员矩阵（cells × nhoods，0/1），供 miloR::nhoods(milo) 使用
#   milopy 一般把成员矩阵放在 adata.obsm["nhoods"]；如果你的对象没有，请先在 Python 侧生成/取到它
M = Mac2025.obsm["nhoods"]  # 稀疏 csr_matrix，形状 = n_cells × n_nhoods
assert M.shape == (Mac2025.n_obs, len(Nhood))

# 保存为Matrix Market + 行/列名（R端很好读）
mmwrite("/storage/data/KAI/Pan_Inflammation/Homo/Mac2025/analysis/milopy/data/nhood_membership.mtx", M.tocoo())
pd.Series(Cells).to_csv("/storage/data/KAI/Pan_Inflammation/Homo/Mac2025/analysis/milopy/data/nhood_cells.txt", index=False, header=False)
pd.Series(Nhood).to_csv("/storage/data/KAI/Pan_Inflammation/Homo/Mac2025/analysis/milopy/data/nhood_ids.txt", index=False, header=False)



