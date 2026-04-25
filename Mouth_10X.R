#### Package ------------------------------
library(dplyr)
library(ggplot2)
library(readxl)
library(patchwork)
library(tidyr)

library(Startrac)
library(tictoc)
library(ggpubr)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(tidyverse)
library(sscVis)
library(Seurat)
library(SeuratData)
library(tidyverse)
library(readr)
library(qs)
library(BiocParallel)
library(coin)
library(Deducer)
library(boot)

library(ggunchull)
library(tidydr)

### specific
library(spacexr)
library(Matrix)


#### GSM6258251 -------------------------------------------------

### GSE206621 数据中不提供表达谱的h5文件
### 自行编写读取函数

LoadSpatial <- function(data.dir, assay = 'Spatial', slice = data.dir){
  data <- Read10X(data.dir = file.path(data.dir, 'filtered_feature_bc_matrix'))
  tissue <- read.csv(file.path(data.dir, 'spatial', "tissue_positions_list.csv"), header = F)
  tissuetype <- tissue[,2]
  names(tissuetype) <- tissue[,1]
  object <- CreateSeuratObject(counts = data, assay = assay, project = data.dir)
  object <- AddMetaData(object, metadata = tissuetype, col.name = "tissue")
  image <- Read10X_Image( image.dir = file.path(data.dir, 'spatial'), filter.matrix = TRUE )
  image <- image[Cells(x = object)]
  DefaultAssay(object = image) <- assay
  object[[slice]] <- image
  object <- subset(object, tissue == 1)
  return(object)
}


GSM6258251 <- LoadSpatial(
  data.dir = "/storage/data/KAI/Pan_Inflammation/Homo/STRNA/Mouth/Periodontitis/1.GSE206621/GSE206621_RAW/GSM6258251/"
  )
Idents(GSM6258251) <- 'GSM6258251'

### QC

VlnPlot(GSM6258251, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()


#### scRNA-seq reference -------------------------------------
Mouth_ref <- readRDS("/storage/data/KAI/Pan_Inflammation/Homo/STRNA/Mouth/RCTD_ref/MouthRCTD_ref.RDS")

counts <- as.matrix(GetAssayData(Mouth_ref, assay = "RNA", layer = "counts"))
cluster <- as.factor(Mouth_ref$cell_type_RCTD)
levels(cluster) <- gsub("/", "_", levels(cluster))
nUMI <- Matrix::colSums(counts)
set.seed(123)
reference <- Reference(counts, cluster, nUMI, n_max_cells = 200)


saveRDS(reference,'/storage/data/KAI/Pan_Inflammation/Homo/STRNA/Mouth/RCTD_ref/reference.RDS')


#### 批量读取GSE206621 ---------------------------------------
setwd('/storage/data/KAI/Pan_Inflammation/Homo/STRNA/Mouth/Periodontitis/1.GSE206621/GSE206621_RAW/')
# LoadSpatial <- function(data.dir, assay = 'Spatial', slice = data.dir){
#   data <- Read10X(data.dir = file.path(data.dir, 'filtered_feature_bc_matrix'))
#   tissue <- read.csv(file.path(data.dir, 'spatial', "tissue_positions_list.csv"), header = F)
#   tissuetype <- tissue[,2]
#   names(tissuetype) <- tissue[,1]
#   object <- CreateSeuratObject(counts = data, assay = assay, project = data.dir)
#   object <- AddMetaData(object, metadata = tissuetype, col.name = "tissue")
#   image <- Read10X_Image( image.dir = file.path(data.dir, 'spatial'), filter.matrix = TRUE )
#   image <- image[Cells(x = object)]
#   DefaultAssay(object = image) <- assay
#   object[[slice]] <- image
#   object <- subset(object, tissue == 1)
#   return(object)
# }

LoadSpatial <- function(data.dir, assay = "Spatial") {
  sample_id <- basename(data.dir)                 # e.g., 'GSM6258251'
  stopifnot(nzchar(sample_id))
  
  counts <- Read10X(data.dir = file.path(data.dir, "filtered_feature_bc_matrix"))
  obj <- CreateSeuratObject(counts = counts,
                            assay = assay,
                            project = sample_id)  # 项目名=短ID，而不是长路径
  
  # 组织掩膜（可选）
  tissue <- read.csv(file.path(data.dir, "spatial", "tissue_positions_list.csv"),
                     header = FALSE)
  tissuetype <- tissue[, 2]; names(tissuetype) <- tissue[, 1]
  obj <- AddMetaData(obj, metadata = tissuetype, col.name = "tissue")
  
  # 读取图像；切片名用 sample_id，避免后面 slot 名重复且难看
  img <- Read10X_Image(image.dir = file.path(data.dir, "spatial"),
                       filter.matrix = TRUE)
  img <- img[Cells(obj)]
  DefaultAssay(img) <- assay
  obj[[sample_id]] <- img
  
  # 只保留 tissue==1 的 spot（按需）
  obj <- subset(obj, tissue == 1)
  
  # ——v5：保证 Spatial assay 的默认层为 "counts"——
  # 通常会有一个 "counts" 层；稳妥起见检查一下
  lyr <- Layers(obj[[assay]])
  if (!"counts" %in% lyr) {
    # 找到第一个 counts-like 层，复制为 "counts" 并删除旧层（安全重命名）
    counts_like <- lyr[grepl("^counts", lyr)][1]
    if (is.na(counts_like)) counts_like <- lyr[1]
    m <- GetAssayData(obj, assay = assay, layer = counts_like)
    obj[[assay]] <- AddLayer(obj[[assay]], layer = m, name = "counts")
    obj[[assay]] <- DropLayers(obj[[assay]], layers = counts_like)
  }
  try(DefaultLayer(obj[[assay]]) <- "counts", silent = TRUE)
  
  # 给细胞条形码加前缀，避免后续合并重名
  obj <- RenameCells(obj, new.names = paste0(sample_id, "_", colnames(obj)))
  
  # 记录样本 ID
  obj$orig.ident <- sample_id
  
  return(obj)
}


GSE206621_list <- as.list(list.files())
GSE206621 <- lapply(GSE206621_list, LoadSpatial)

scttransform <- function(object){
  return(SCTransform(object, assay = "Spatial", verbose = FALSE))
}

GSE206621 <- lapply(GSE206621, scttransform)
GSE206621_vars <- lapply(GSE206621, VariableFeatures)


GSE206621_all <- merge(x = GSE206621[[1]], y = GSE206621[-1])


SpatialFeaturePlot(GSE206621_all, features = "nCount_Spatial", pt.size.factor = 3.5) + theme(legend.position = "right")


#### GSE206621每个slice进行RCTD -----------------------------

run_rctd_for_sample <- function(obj_merged, sample_id, reference,
                                assay = "Spatial",
                                max_cores = 4,
                                save_prefix = "rctd_result"
                                ) {
  message(">>> Running RCTD for sample: ", sample_id)
  
  # 1) 取该样本的 counts（在合并对象的 layer 里）
  layer_name <- paste0("counts.", sample_id)
  stopifnot(layer_name %in% Layers(obj_merged[[assay]]))
  
  counts <- GetAssayData(obj_merged, assay = assay, layer = layer_name)
  # 稀疏格式（一般已经是 dgCMatrix）
  counts <- as(counts, "dgCMatrix")
  
  # 2) 取该样本的空间坐标（Seurat v5：image 名就是 sample_id）
  coords_df <- GetTissueCoordinates(obj_merged, image = sample_id)
  # 统一列名为 x/y
  coords_df <- coords_df[, c("x", "y"), drop = FALSE]
  # colnames(coords_df) <- c("x", "y")
  
  # 3) 对齐（barcodes）
  common <- intersect(rownames(coords_df), colnames(counts))
  coords <- as.matrix(coords_df[common, , drop = FALSE])
  counts <- counts[, common, drop = FALSE]
  stopifnot(identical(rownames(coords), colnames(counts)))
  
  # 4) 构建 SpatialRNA
  spatialRNA <- SpatialRNA(coords %>% as.data.frame(), counts)
  
  # 5) 运行 RCTD（先 bulk 后 pixels，更稳）
  set.seed(123)
  rctd <- create.RCTD(spatialRNA, reference,
                      max_cores = max_cores
                      )
  
  rctd <- run.RCTD(rctd, doublet_mode = "full") 
  
  saveRDS(rctd,paste0('/storage/data/KAI/Pan_Inflammation/Homo/STRNA/Mouth/Periodontitis/1.GSE206621/GSE206621_RCTD/', sprintf("%s_%s.RDS", sample_id, save_prefix )))
  message("<<< Done: ", sample_id)
  return(list(obj = obj_merged, rctd = rctd))
}

sample_ids <- sub("^counts\\.", "", Layers(GSE206621_all[["Spatial"]]))


for (sid in sample_ids) {
  out <- run_rctd_for_sample(GSE206621_all, sample_id = sid, reference = reference,
                             assay = "Spatial",
                             max_cores = 4
                             )
  GSE206621_all <- out$obj
}



#### 提取RCTD结果 -----------------------------
result_dir <- "/storage/data/KAI/Pan_Inflammation/Homo/STRNA/Mouth/Periodontitis/1.GSE206621/GSE206621_RCTD/"

# 找到所有 rctd_result.RDS 文件
files <- list.files(result_dir, pattern = "_rctd_result\\.RDS$", full.names = TRUE)

# 批量读取并放进一个 list
rctd_results <- lapply(files, readRDS)

# 给 list 命名（提取 GSM 编号）
names(rctd_results) <- sub("_rctd_result\\.RDS$", "", basename(files))

# 检查
names(rctd_results)







### GSM6258251
GSM6258251 <- subset(GSE206621_all, subset = orig.ident == 'GSM6258251')

GSM6258251@meta.data[setdiff(colnames(GSM6258251),rownames(weights)), c("nCount_Spatial","nFeature_Spatial","tissue")]

# GSM6258251.RCTD <- readRDS("/storage/data/KAI/Pan_Inflammation/Homo/STRNA/Mouth/Periodontitis/1.GSE206621/GSE206621_RCTD/GSM6258251_rctd_result.RDS")
GSM6258251.RCTD<- rctd_results[["GSM6258251"]]
GSM6258251.RCTD.result <- GSM6258251.RCTD@results
weights <- GSM6258251.RCTD.result$weights %>% as.data.frame()
norm_weights = normalize_weights(GSM6258251.RCTD.result$weights) 



cell_type_names <- GSM6258251.RCTD@cell_type_info$info[[2]] #list of cell type names
spatialRNA <- GSM6258251.RCTD@spatialRNA

coords <- GetTissueCoordinates(GSE206621_all, image = "GSM6258251")
plot_cell_type_weights <- function(weights, coords, celltype){
  spatial_df <- data.frame(
    x = coords$x,
    y = coords$y,
    weight = weights[,celltype]
  )
  ggplot(spatial_df, aes(x=x, y=y, color = weight)) + 
    geom_point(size = 3) + 
    scale_color_gradientn(colors = c('#5C569E','#3D87AF','#68BEA2', '#A5D5A2',
                                     '#ece399','#ffdb92','#F6C37E','#F18D5E','#DC5C52'),
                          limits = c(0,1)
                          ) + 
    labs(title = paste(celltype, " Proportion")) +
    theme_bw() +
    coord_fixed()
}

plot_cell_type_weights(weights,coords[intersect(colnames(GSM6258251),rownames(weights)),],"SPP1+CCL2+_Macro")



### GSM6258256
GSM6258256 <- subset(GSE206621_all, subset = orig.ident == 'GSM6258256')


# GSM6258251.RCTD <- readRDS("/storage/data/KAI/Pan_Inflammation/Homo/STRNA/Mouth/Periodontitis/1.GSE206621/GSE206621_RCTD/GSM6258251_rctd_result.RDS")
GSM6258256.RCTD<- rctd_results[["GSM6258256"]]
GSM6258256.RCTD.result <- GSM6258256.RCTD@results
weights <- GSM6258256.RCTD.result$weights %>% as.data.frame()
norm_weights = normalize_weights(GSM6258256.RCTD.result$weights) %>% as.data.frame() 

GSM6258256@meta.data[setdiff(colnames(GSM6258256),rownames(weights)), c("nCount_Spatial","nFeature_Spatial","tissue")]

cell_type_names <- GSM6258256.RCTD@cell_type_info$info[[2]] #list of cell type names
spatialRNA <- GSM6258256.RCTD@spatialRNA

coords <- GetTissueCoordinates(GSE206621_all, image = "GSM6258256")

plot_cell_type_weights(weights,coords[intersect(colnames(GSM6258256),rownames(weights)),],"T_NK cell")
 
pearson.corr <- cor.test(norm_weights$`T_NK cell`,norm_weights$`SPP1+CCL2+_Macro`,method="pearson")


### GSM6258257
GSM6258257 <- subset(GSE206621_all, subset = orig.ident == 'GSM6258257')


# GSM6258251.RCTD <- readRDS("/storage/data/KAI/Pan_Inflammation/Homo/STRNA/Mouth/Periodontitis/1.GSE206621/GSE206621_RCTD/GSM6258251_rctd_result.RDS")
GSM6258257.RCTD<- rctd_results[["GSM6258257"]]
GSM6258257.RCTD.result <- GSM6258257.RCTD@results
weights <- GSM6258257.RCTD.result$weights %>% as.data.frame()
norm_weights = normalize_weights(GSM6258257.RCTD.result$weights) %>% as.data.frame() 

GSM6258257@meta.data[setdiff(colnames(GSM6258257),rownames(weights)), c("nCount_Spatial","nFeature_Spatial","tissue")]

cell_type_names <- GSM6258257.RCTD@cell_type_info$info[[2]] #list of cell type names
spatialRNA <- GSM6258257.RCTD@spatialRNA

coords <- GetTissueCoordinates(GSE206621_all, image = "GSM6258257")

plot_cell_type_weights(weights,coords[intersect(colnames(GSM6258257),rownames(weights)),],"SPP1+CCL2+_Macro")

cor(norm_weights, method="pearson")
pearson.corr <- cor.test(norm_weights$`T_NK cell`,norm_weights$`SPP1+CCL2+_Macro`,method="pearson")


### GSM6258258
GSM6258258 <- subset(GSE206621_all, subset = orig.ident == 'GSM6258258')


# GSM6258251.RCTD <- readRDS("/storage/data/KAI/Pan_Inflammation/Homo/STRNA/Mouth/Periodontitis/1.GSE206621/GSE206621_RCTD/GSM6258251_rctd_result.RDS")
GSM6258258.RCTD<- rctd_results[["GSM6258258"]]
GSM6258258.RCTD.result <- GSM6258258.RCTD@results
weights <- GSM6258258.RCTD.result$weights %>% as.data.frame()
norm_weights = normalize_weights(GSM6258258.RCTD.result$weights) %>% as.data.frame() 

GSM6258258@meta.data[setdiff(colnames(GSM6258258),rownames(weights)), c("nCount_Spatial","nFeature_Spatial","tissue")]

cell_type_names <- GSM6258258.RCTD@cell_type_info$info[[2]] #list of cell type names
spatialRNA <- GSM6258258.RCTD@spatialRNA

coords <- GetTissueCoordinates(GSE206621_all, image = "GSM6258258")

plot_cell_type_weights(weights,coords[intersect(colnames(GSM6258258),rownames(weights)),],"T_NK cell")

cor(norm_weights, method="pearson")
pearson.corr <- cor.test(norm_weights$`T_NK cell`,norm_weights$`SPP1+CCL2+_Macro`,method="pearson")


#### 将normalize_weights导入meta.data ---------------------------------------------

samples <- paste0("GSM625825", 1:8)

for (s in samples) {
  # 1. 提取当前样本 Seurat 子集
  seurat_sub <- subset(GSE206621_all, subset = orig.ident == s)
  
  # 2. 提取 RCTD normalize_weights (spot x celltypes)
  norm_w <- normalize_weights(rctd_results[[s]]@results$weights) %>% as.data.frame()   # 请确认名字对，可能是 norm_weights
  
  # 3. 确认行名一致（spot ID），确保和 Seurat 对应
  # 如果 RCTD 的行名包含 sample 前缀，可以统一一下
  common_spots <- intersect(colnames(seurat_sub), rownames(norm_w))
  norm_w <- norm_w[common_spots, , drop = FALSE]
  
  # 4. 将 weights 加入到 metadata
  seurat_sub <- AddMetaData(seurat_sub, metadata = norm_w)
  
  # 5. 更新回大对象（替换子集）
  GSE206621_all@meta.data[common_spots, colnames(norm_w)] <- seurat_sub@meta.data[common_spots, colnames(norm_w)]
}


SpatialFeaturePlot(
  subset(GSE206621_all, subset = orig.ident %in% c("GSM6258256","GSM6258257","GSM6258258")),
  features = c("SPP1+CCL2+_Macro", "T_NK cell"),
  ncol = 3,
  pt.size.factor = 3.5
) # size 10*8







#### SPP1+CCL2+_Macro - T/NK cell pearson -----------------------------------

### 提取所有样本的RCTD weights
samples <- paste0("GSM625825", 1:8)


res_list <- list()

for (s in samples) {
  # 提取当前样本的 RCTD 权重矩阵（spot x celltypes）
  w <- normalize_weights(rctd_results[[s]]@results$weights) %>% as.data.frame()  # 确保这里是权重矩阵
  # 确认有SPP1+CCL2+_Macro
  if(!"SPP1+CCL2+_Macro" %in% colnames(w)) next
  
  macro <- w[, "SPP1+CCL2+_Macro"]
  
  # 遍历其他细胞类型
  for (ct in setdiff(colnames(w), "SPP1+CCL2+_Macro")) {
    test <- cor.test(macro, w[, ct], method = "pearson")
    
    res_list[[length(res_list) + 1]] <- data.frame(
      sample = s,
      celltype = ct,
      cor = test$estimate,
      pval = test$p.value
    )
  }
}

res_df <- do.call(rbind, res_list)
res_df <- res_df %>%
  mutate(
    size_cat = case_when(
      pval < 0.01 ~ "p < 0.01",
      pval < 0.05 ~ "p < 0.05",
      TRUE ~ "p > 0.05"
    )
  )
# plot_df <- res_df %>% filter(size_cat != "ns")


res_df$celltype <- factor(res_df$celltype, levels = rev(c("Epithelial cell","Endothelial cell","Fibroblast",
                                                      "IL1B+NLRP3+_Macro","Other Macro","Dendritic cell",
                                                      "Mast cell","Neutrophil","T_NK cell","B cell","Plasma cell",
                                                      "Neural cell"
                                                      )) )
saveRDS(res_df,"/storage/data/KAI/Pan_Inflammation/Homo/STRNA/Mouth/Periodontitis/1.GSE206621/GSE206621_R/analysis/RCTD/PearsonRes_SPP1_other.RDS")
ggplot(res_df, aes(x = sample, y = celltype)) +
  # 用 tile 画网格，调整宽高贴近点
  geom_tile(color = "black", fill = "white", width = 0.9, height = 0.9) +
  # 叠加 bubble
  geom_point(aes(color = cor, size = size_cat)) +
  scale_size_manual(values = c("p > 0.05" = 2,"p < 0.05" = 4, "p < 0.01" = 6)) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  coord_fixed(ratio = 1) +  # 保持格子正方形
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),        # 去掉背景灰色网格
    axis.text.x = element_text(angle = 90, hjust = 0.5)
  ) +
  labs(x = "Sample", y = "Cell type", color = "Pearson r", size = "P value scale")



ggplot(res_df, aes(x = sample, y = celltype)) +
  # 用 tile 画网格
  geom_tile(color = "black", fill = "white", width = 0.9, height = 0.9) +
  # 叠加 bubble
  geom_point(aes(color = cor, size = size_cat)) +
  scale_size_manual(values = c("p > 0.05" = 2,"p < 0.05" = 4, "p < 0.01" = 6)) +
  # 使用自定义调色板
  scale_color_gradientn(
    colours = c("#1f294e","#5390b5", "#eaebea","#d56e5e","#57121d"),
    limits = c(-1, 1)   # Pearson 相关范围
  ) +
  coord_fixed(ratio = 1) +  
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 0.5)
  ) +
  labs(x = "Sample", y = "Cell type", color = "Pearson r", size = "P value scale") # size 7*5



