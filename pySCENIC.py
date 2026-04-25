import os,sys
import loompy as lp
import numpy as np
import pandas as pd
import scanpy as sc
from dask.distributed import Client
from scipy.sparse import csr_matrix

#--------------------- full size macrophage  ------------------------
pan_Macro_filtered = sc.read_h5ad("/storage/data/KAI/Pan_Inflammation/Homo/Infla_R/Macrophage/pan_Macro_nodoublets_scanvi.h5ad")

row_attrs = {"Gene": np.array(pan_Macro_filtered.var_names), };
col_attrs = {"CellID": np.array(pan_Macro_filtered.obs_names)};


lp.create("/storage/data/KAI/Pan_Inflammation/Homo/Infla_R/Macrophage/pyscenic/input/sce.loom", pan_Macro_filtered.X.transpose(), row_attrs, col_attrs)

# 抽样检测pyscenic是否正常运行
np.random.seed(77)
n_cells = pan_Macro_filtered.n_obs
sample_size = int(n_cells * 0.01)

selected_indices = np.random.choice(n_cells, sample_size, replace=False)
pan_Macro_subset = pan_Macro_filtered[selected_indices, :]

row_attrs = {"Gene": np.array(pan_Macro_subset.var_names), };
col_attrs = {"CellID": np.array(pan_Macro_subset.obs_names)};

lp.create("/storage/data/KAI/Pan_Inflammation/Homo/Infla_R/Macrophage/pyscenic/input/sce_downsample.loom", pan_Macro_subset.X.transpose(), row_attrs, col_attrs)

# 检查写出的loom文件
with lp.connect("/storage/data/KAI/Pan_Inflammation/Homo/Infla_R/Macrophage/pyscenic/input/sce_downsample.loom") as ds:
    print("Shape:", ds.shape)  # 打印数据维度
    print("Genes:", ds.ra.Gene[:5])  # 打印前5个基因
    print("Cells:", ds.ca.CellID[:5])  # 打印前5个细胞ID


client = Client()  # 创建 Dask 客户端
print(client)


#--------------------- downsampled macrophage  ------------------------

np.random.seed(10)  # 设置随机种子以确保可重复性
sample_indices = pan_Macro_filtered.obs.groupby('cell_type_level2').apply(
    lambda x: x.sample(frac=0.1, random_state=10).index
).explode()

Macro_downsampled = pan_Macro_filtered[sample_indices]
Macro_downsampled.write("/storage/data/KAI/Pan_Inflammation/Homo/Infla_R/Macrophage/Macrophage_downsampled_10/Macro_downsampled_10.h5ad")
X_f32 = csr_matrix(Macro_downsampled.X.toarray().astype('float32'))
Macro_downsampled.X = X_f32

Macro_downsampled.write("/storage/data/KAI/Pan_Inflammation/Homo/Infla_R/Macrophage/Macrophage_downsampled_10/Macro_downsampled_10_f32.h5ad")


row_attrs = {"Gene": np.array(Macro_downsampled.var_names), };
col_attrs = {"CellID": np.array(Macro_downsampled.obs_names)};

lp.create("/storage/data/KAI/Pan_Inflammation/Homo/Infla_R/Macrophage/pyscenic/input/sce_downsampled.loom", Macro_downsampled.X.transpose(), row_attrs, col_attrs)



############################################# Mac 2025 #############################################
pan_Mac = sc.read_h5ad("/storage/data/KAI/Pan_Inflammation/Homo/Mac2025/scVI/new/pan_Mac_counts.h5ad")
Mac10_barcodes = pd.read_csv("/storage/data/KAI/Pan_Inflammation/Homo/Mac2025/analysis/SCENIC/data/Mac10_barcodes.txt",sep='\t')
Mac10_barcodes = Mac10_barcodes['x'].to_list()

Mac_10 = pan_Mac[pan_Mac.obs_names.isin(Mac10_barcodes)].copy()

Mac_10.write("/storage/data/KAI/Pan_Inflammation/Homo/Mac2025/analysis/SCENIC/data/Mac_10.h5ad")


### write loom
row_attrs = {"Gene": np.array(Mac_10.var_names), };
col_attrs = {"CellID": np.array(Mac_10.obs_names)};

lp.create("/storage/data/KAI/Pan_Inflammation/Homo/Mac2025/analysis/SCENIC/data/sce.loom", Mac_10.X.transpose(), row_attrs, col_attrs)


############################################# Fibro 2025 #############################################
Fibro = sc.read_h5ad("/storage/data/KAI/Pan_Inflammation/Homo/Fibro/scVI/new/Fibro_2025/pan_Fibro.h5ad")
Fibro10_barcodes = pd.read_csv("/storage/data/KAI/Pan_Inflammation/Homo/Fibro/analysis/SCENIC/data/Fibro10_barcodes.txt",sep='\t')
Fibro10_barcodes = Fibro10_barcodes['x'].to_list()
Fibro_10 = Fibro[Fibro.obs_names.isin(Fibro10_barcodes)].copy()
Fibro_10.write("/storage/data/KAI/Pan_Inflammation/Homo/Fibro/analysis/SCENIC/data/Fibro_10.h5ad")

### write loom
row_attrs = {"Gene": np.array(Fibro_10.var_names), };
col_attrs = {"CellID": np.array(Fibro_10.obs_names)};

lp.create("/storage/data/KAI/Pan_Inflammation/Homo/Fibro/analysis/SCENIC/data/sce.loom", Fibro_10.X.transpose(), row_attrs, col_attrs)




# bash---------------------------------
nohup pyscenic grn   --output /storage/data/KAI/Pan_Inflammation/Homo/Fibro/analysis/SCENIC/data/sce_adj_loomR.csv   --method grnboost2   /storage/data/KAI/Pan_Inflammation/Homo/Fibro/analysis/SCENIC/data/sce.loom   /storage/data/KAI/Pan_Inflammation/Homo/Infla_R/Macrophage/pyscenic/input_re-annotation/allTFs_hg38.txt > /storage/data/KAI/Pan_Inflammation/Homo/Fibro/analysis/SCENIC/data/pyscenic_macro.log 2>&1 &
# PID 604852
# done


nohup pyscenic ctx --num_workers 20 \
--output /storage/data/KAI/Pan_Inflammation/Homo/Fibro/analysis/SCENIC/data/sce.regulons.csv \
--expression_mtx_fname /storage/data/KAI/Pan_Inflammation/Homo/Fibro/analysis/SCENIC/data/sce.loom \
--all_modules \
--mask_dropouts \
--mode "dask_multiprocessing" \
--min_genes 10 \
--annotations_fname /storage/data/KAI/Pan_Inflammation/Homo/Infla_R/Macrophage/pyscenic/input_re-annotation/motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
/storage/data/KAI/Pan_Inflammation/Homo/Fibro/analysis/SCENIC/data/sce_adj_loomR.csv \
/storage/data/KAI/Pan_Inflammation/Homo/Infla_R/Macrophage/pyscenic/input_re-annotation/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather \
> /storage/data/KAI/Pan_Inflammation/Homo/Fibro/analysis/SCENIC/data/pyscenic_macro_ctx.log 2>&1 &
# PID
#进行TF-motif富集分析，识别直接靶标
#得到转录因子(TF)与其对应的直接作用的靶点,称为regulon(每一个regulon是1个TF和其调控的靶基因)
# done


nohup pyscenic aucell --num_workers 3 --output /storage/data/KAI/Pan_Inflammation/Homo/Mac2025/analysis/SCENIC/data/sce_SCENIC.loom \
/storage/data/KAI/Pan_Inflammation/Homo/Mac2025/analysis/SCENIC/data/sce.loom \
/storage/data/KAI/Pan_Inflammation/Homo/Mac2025/analysis/SCENIC/data/sce.regulons.csv \
> /storage/data/KAI/Pan_Inflammation/Homo/Mac2025/analysis/SCENIC/data/pyscenic_macro_aucell.log 2>&1 &
# PID 780948
#使用AUCell对每个细胞的每个regulon活性进行评分。




# --------------------------------------





