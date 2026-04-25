#### Package ------------------------------
### global
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
library(future)
library(ggalluvial)
library(AUCell)
library(upstartr)
library(multcompView)
library(reticulate)
library(NMF)

#### Global -----------------------------
level3_color <- c(
  'Epithelial cell'='#a78982',
  'Endothelium'='#cc7f73',
  'Fibroblasts'='#c4612f',
  'Myofibroblasts'='#837d52',
  'Inflammatory fibroblasts'='#d19f9f',
  'S1'='#00c2b2',
  'S2a'='#6c408e',
  'S2b'='#9a70a8',
  'S3'='#814e3a',
  'Pericytes'='#efd2c9',
  'FRCs'='#6aa57e',
  'Glia'='#b7deea',
  'Inflammatory monocytes'='#cea5c7',
  'IL1B+NLRP3+_Macro'='#FA8072',
  'SPP1+CCL2+_Macro'='#ff9d5c',
  'Other Macro'='#FFC5BF',
  'Mast'='#e29eaf',
  'DCs'='#927c9a',
  'Eosinophils'='#83ab8e',
  'N1'='#d25774',
  'N2'='#e6e2a3',
  'N3'='#c6adb0',
  'Cycling myeloid'='#a5a9b0',
  'T cell'='#2d3462',
  'B/Plasma cell'='#4b6aa8',
  "Macrophage NRG1"='#25b1bf',
  "M0"='#737b5e',        
  "M1"='#98585c',
  "M2"='#463e4b'
)



#### load data ---------------------------
cosmx_mye <- readRDS("/storage/data/KAI/Pan_Inflammation/Homo/STRNA/Colorectum/1.GSE234713(CosMx)/squidpy/mye/scvi/colonmye_scvi.RDS")
cosmx_all <- readRDS("/storage/data/KAI/Pan_Inflammation/Homo/STRNA/Colorectum/1.GSE234713(CosMx)/squidpy/scarches/data/GSE234713_anno.RDS")
## 等待修改cell_type_level3和cell_arches之间的映射关系

#### marker_correlation (Pearson, no good) -----------------
### only macrophage subtypes in colorectum
cosmx_macro <- readRDS("/storage/data/KAI/Pan_Inflammation/Homo/STRNA/Colorectum/1.GSE234713(CosMx)/squidpy/macro/data/target_macro.RDS")
scrna_macro <- readRDS("/storage/data/KAI/Pan_Inflammation/Homo/STRNA/Colorectum/1.GSE234713(CosMx)/squidpy/macro/data/source_macro.RDS")

### normalized
cosmx_macro <- NormalizeData(cosmx_macro, normalization.method = "LogNormalize", scale.factor = 10000)
scrna_macro <- NormalizeData(scrna_macro, normalization.method = "LogNormalize", scale.factor = 10000)

Idents(cosmx_macro) <- cosmx_macro$cell_arches
Idents(scrna_macro) <- scrna_macro$cell_type_level3

### common markers
cosmx.markers <- FindAllMarkers(cosmx_macro, only.pos = TRUE)
scrna.markers <- FindAllMarkers(scrna_macro, only.pos = TRUE)
saveRDS(cosmx.markers,'/storage/data/KAI/Pan_Inflammation/Homo/Result_Fig/Fig6/material/colon_cosmx/marker_correlation/data/cosmx_topmarkers.RDS')
saveRDS(scrna.markers,'/storage/data/KAI/Pan_Inflammation/Homo/Result_Fig/Fig6/material/colon_cosmx/marker_correlation/data/scrna_topmarkers.RDS')

common_markers <- intersect(cosmx.markers$gene, scrna.markers$gene) # 642

### pseudobulk 
scrna.avg <- AggregateExpression(scrna_macro, return.seurat = TRUE)
cosmx.avg <- AggregateExpression(cosmx_macro, return.seurat = TRUE)

scrna.expr <- GetAssayData(scrna.avg, slot = "data")
cosmx.expr <- GetAssayData(cosmx.avg, slot = "data")

colnames(scrna.expr) <- c("CCL13+_C1Q+_Macro","CD1C+_DC-like_Macro","IL1B+NLRP3+_Macro",
                          "Intestinal resident Macro","MMP3+CXCL8+_Macro","SPP1+CCL2+_Macro")
colnames(cosmx.expr) <- c("IL1B+NLRP3+_Macro","Other Macro","SPP1+CCL2+_Macro")

# 统一前缀
colnames(scrna.expr) <- paste0("scRNA_", colnames(scrna.expr))
colnames(cosmx.expr) <- paste0("CosMx_", colnames(cosmx.expr))


### Pearson correlation (no common markers, all panel genes)
expr_mat <- cbind(scrna.expr[common_markers,], cosmx.expr[common_markers,])

sc_cols <- grep("^scRNA_", colnames(expr_mat), value = TRUE)
cosmx_cols <- grep("^CosMx_", colnames(expr_mat), value = TRUE)

# 初始化矩阵
cor_mat <- matrix(NA, nrow = length(cosmx_cols), ncol = length(sc_cols),
                  dimnames = list(cosmx_cols, sc_cols))
pval_mat <- matrix(NA, nrow = length(cosmx_cols), ncol = length(sc_cols),
                   dimnames = list(cosmx_cols, sc_cols))

# 循环计算
for (i in cosmx_cols) {
  for (j in sc_cols) {
    test <- cor.test(expr_mat[, i], expr_mat[, j], method = "spearman")
    cor_mat[i, j] <- test$estimate
    pval_mat[i, j] <- test$p.value
  }
}


Heatmap(cor_mat,
        name = "Pearson\nCorrelation",
        col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
        cell_fun = function(j, i, x, y, width, height, fill) {
          if (pval_mat[i, j] < 0.05) {
            grid.text(sprintf("%.2f*", cor_mat[i, j]), x, y, gp = gpar(fontsize = 8))
          } else {
            grid.text(sprintf("%.2f", cor_mat[i, j]), x, y, gp = gpar(fontsize = 8))
          }
        },
        row_names_side = "left",
        column_names_side = "top",
        column_names_rot = 45)

#### top marker AUCell score ----------------------------

pan_Macro <- readRDS("/storage/data/KAI/Pan_Inflammation/Homo/Mac2025/scVI/new/pan_Mac_counts.RDS") 
pan_Macro <- NormalizeData(pan_Macro, normalization.method = "LogNormalize", scale.factor = 10000)
Idents(pan_Macro) <- pan_Macro@meta.data$cell_type_level3

scrna_panmacro.markers <- FindAllMarkers(pan_Macro, only.pos = TRUE)
saveRDS(scrna_panmacro.markers,'/storage/data/KAI/Pan_Inflammation/Homo/Result_Fig/Fig6/material/colon_cosmx/marker_correlation/data/scrna_panmacro.markers.RDS')


markers_IL1B <- scrna_panmacro.markers %>%
  filter(cluster == "IL1B+NLRP3+_Macro") %>%
  filter(gene %in% rownames(cosmx_macro)) %>%
  arrange(desc(avg_log2FC)) %>%
  # slice_head(n = 50) %>%
  pull(gene)

markers_SPP1 <- scrna_panmacro.markers %>%
  filter(cluster == "SPP1+CCL2+_Macro") %>%
  filter(gene %in% rownames(cosmx_macro)) %>%
  arrange(desc(avg_log2FC)) %>%
  # slice_head(n = 50) %>%
  pull(gene)

geneSets <- list(
  IL1B_NLRP3_Score = markers_IL1B,
  SPP1_CCL2_Score = markers_SPP1
)


Mac_rankings <- AUCell_buildRankings(GetAssayData(cosmx_macro, layer = "counts"))

Mac_AUC <- AUCell_calcAUC(geneSets, Mac_rankings)
# 通过 AUC 评分直方图 选择一个合适的阈值，用于区分基因集活跃的细胞与不活跃的细胞。
Mac_assignment <- AUCell_exploreThresholds(Mac_AUC, plotHist = TRUE, assign=TRUE)

### 结果存入seurat对象 Mac中
all(colnames(cosmx_macro) == colnames(Mac_AUC)) # TRUE

Macro_scores <- as.data.frame(t(getAUC(Mac_AUC)))

cosmx_macro[["IL1B_NLRP3_Score"]] <- NULL
cosmx_macro[["SPP1_CCL2_Score"]] <- NULL

cosmx_macro@meta.data <- cbind(cosmx_macro@meta.data, Macro_scores)
saveRDS(cosmx_macro,"/storage/data/KAI/Pan_Inflammation/Homo/STRNA/Colorectum/1.GSE234713(CosMx)/squidpy/macro/data/target_macro.RDS")

### 小提琴图可视化结果
VlnPlot(
  cosmx_macro,
  features = "IL1B_NLRP3_Score",
  group.by = "cell_arches", # 或你的迁移标签字段
  pt.size = 0
) 


df <- FetchData(cosmx_macro, vars = c("IL1B_NLRP3_Score", "cell_arches"))
wilcox.test(
  IL1B_NLRP3_Score ~ cell_arches,
  data = subset(df, cell_arches %in% c("IL1B+NLRP3+_Macro", "SPP1+CCL2+_Macro"))
)

### 云雨图可视化结果
library(tidyverse)
library(gghalves)
library(Hmisc)
library(cowplot)
library(ggdist)

visual_df <- cosmx_macro@meta.data[, c("cell_arches", "IL1B_NLRP3_Score", "SPP1_CCL2_Score")]

## AUCell score

visual_df$cell_arches <- factor(visual_df$cell_arches, 
                                levels = c("IL1B+NLRP3+_Macro", "SPP1+CCL2+_Macro", "Other Macro"))
long_df <- visual_df %>%
  pivot_longer(cols = c(IL1B_NLRP3_Score, SPP1_CCL2_Score),
               names_to = "score_type",
               values_to = "score")

saveRDS(long_df,'/storage/data/KAI/Pan_Inflammation/Homo/Result_Fig/Fig6/material/colon_cosmx/marker_correlation/data/cosmx_il1bspp1_aucell_longdf.RDS')

long_df <- readRDS('/storage/data/KAI/Pan_Inflammation/Homo/Result_Fig/Fig6/material/colon_cosmx/marker_correlation/data/cosmx_il1bspp1_aucell_longdf.RDS')

## IL1B_NLRP3_Score

long_df_il1b <- long_df %>% filter(score_type == "IL1B_NLRP3_Score")
long_df_il1b$cell_arches <- factor(long_df_il1b$cell_arches, levels = c("IL1B+NLRP3+_Macro", "SPP1+CCL2+_Macro", "Other Macro"))

ggplot(long_df_il1b, aes(x = cell_arches, y = score, fill = cell_arches)) +  
  stat_halfeye(adjust = 0.6, width = 0.4, .width = 0, justification = -0.2, point_colour = NA) +   
  geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.5, color = "black") +  
  # stat_dots(side = "left", justification = 1.01, binwidth = 1.01, dotsize = 0.1, stroke = 0.1) +   
  stat_compare_means(
    method = "wilcox.test",
    comparisons = list(
      c("IL1B+NLRP3+_Macro", "SPP1+CCL2+_Macro"),
      c("IL1B+NLRP3+_Macro", "Other Macro")
    ),
    label = "p.format",
    vjust = 0.1,
    size = 6,
    tip.length = 0.015,
    bracket.size = 0.8,
    method.args = list(alternative = "greater"),
    y.position = c(1.05, 1.15)  # 两个显著性标签的高度
  ) +
  scale_fill_manual(values = level3_color) +  
  labs(x = "CosMx Macrophages", y = "AUCell (IL1B+NLRP3+_Macro)") +  
  theme_classic(base_size = 20) +  
  theme(axis.ticks.length = unit(0.3, "cm"),
        axis.text = element_text(size = 18, color = "black"),
        axis.line = element_line(linewidth = 0.8, color = "black"),  # 加粗轴线        
        axis.ticks = element_line(linewidth = 0.8, color = "black"),  # 加粗刻度线        
        legend.position = "none") # size 5*4


visual_df$cell_arches <- factor(visual_df$cell_arches, 
                                levels = c( "SPP1+CCL2+_Macro","IL1B+NLRP3+_Macro", "Other Macro"))
long_df <- visual_df %>%
  pivot_longer(cols = c(IL1B_NLRP3_Score, SPP1_CCL2_Score),
               names_to = "score_type",
               values_to = "score")

## IL1B_NLRP3_Score

long_df_spp1 <- long_df %>% filter(score_type == "SPP1_CCL2_Score")

long_df_spp1$cell_arches <- factor(long_df_spp1$cell_arches, levels = c( "SPP1+CCL2+_Macro","IL1B+NLRP3+_Macro", "Other Macro"))
ggplot(long_df_spp1, aes(x = cell_arches, y = score, fill = cell_arches)) +  
  stat_halfeye(adjust = 0.6, width = 0.4, .width = 0, justification = -0.2, point_colour = NA) +   
  geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.5, color = "black") +  
  # stat_dots(side = "left", justification = 1.01, binwidth = 1.01, dotsize = 0.1, stroke = 0.1) +   
  stat_compare_means(
    method = "wilcox.test",
    comparisons = list(
      c("SPP1+CCL2+_Macro","IL1B+NLRP3+_Macro"),
      c("SPP1+CCL2+_Macro", "Other Macro")
    ),
    label = "p.format",
    vjust = 0.1,
    size = 6,
    tip.length = 0.015,
    bracket.size = 0.8,
    method.args = list(alternative = "greater"),
    y.position = c(1.05, 1.15)  # 两个显著性标签的高度
  ) +
  scale_fill_manual(values = level3_color) +  
  labs(x = "CosMx Macrophages", y = "AUCell (SPP1+CCL2+_Macro)") +  
  theme_classic(base_size = 20) +  
  theme(axis.ticks.length = unit(0.3, "cm"),
        axis.text = element_text(size = 18, color = "black"),
        axis.line = element_line(linewidth = 0.8, color = "black"),  # 加粗轴线        
        axis.ticks = element_line(linewidth = 0.8, color = "black"),  # 加粗刻度线        
        legend.position = "none") 






#### proportion (all cells)--------------------------
cos_meta <- cosmx_all@meta.data[,c('fov','Sample_geo_accession','cell_type_level3','Donor','Group','Inflammation')]

cell_counts <- cos_meta %>%
  group_by(Donor, Group, cell_type_level3) %>%
  summarise(count = n(), .groups = 'drop')

cell_props <- cell_counts %>%
  group_by(Donor) %>%
  mutate(freq = count / sum(count))


cell_props$cell_type_level3 <- factor(cell_props$cell_type_level3,
                                      levels = c('Epithelial cell','Stromal cell','IL1B+NLRP3+_Macro','SPP1+CCL2+_Macro','Other Macro',
                                                 'Inflammatory monocytes',
                                                 'DCs','Mast','N1','N2','N3','Cycling myeloid','Eosinophils',
                                                 'T cell','B/Plasma cell'
                                                 )
                                      )
cell_props$Group <- factor(cell_props$Group, levels = c('HC','UC','CD'))
cell_props$Donor <- factor(cell_props$Donor, levels = c('HCa','HCb','HCc','UCa','UCb','UCc','CDa','CDb','CDc'))

ggplot(cell_props, aes(x = Donor, y = freq, fill = cell_type_level3)) +
  geom_bar(stat = "identity", width = 0.6) +
  facet_wrap(~ Group, scales = "free_x") +
  scale_fill_manual(values = level3_color) + 
  labs(
    x = "Donor",
    y = "Proportion",
    fill = "Cell Type",
    title = "Cell type composition per Donor"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5),
        panel.grid.major = element_blank(),     # 删除主网格线
        panel.grid.minor = element_blank(),     # 删除次网格线
        panel.background = element_blank()      # 可选：移除背景
        )

#### proportion (macros)--------------------------
cos_meta_macro <- cos_meta[grep('.*Macro$',cos_meta$cell_type_level3,value = F),]
cos_meta_macro <- droplevels(cos_meta_macro)
macro_counts <- cos_meta_macro %>%
  group_by(Donor, Group, cell_type_level3) %>%
  dplyr::summarise(count = n(), .groups = 'drop')

macro_props <- macro_counts %>%
  group_by(Donor, cell_type_level3) %>%
  dplyr::summarise(n = n(), .groups = "drop_last") %>%
  mutate(freq = n / sum(n))
macro_props <- macro_counts %>% 
  dplyr::group_by(Donor) %>%
  dplyr::mutate(freq = count / sum(count))



macro_props$cell_type_level3 <- factor(macro_props$cell_type_level3, 
                                       levels = c('IL1B+NLRP3+_Macro','SPP1+CCL2+_Macro','Other Macro'
                                                 )
                                       )

macro_props$Group <- factor(macro_props$Group, levels = c('HC','UC','CD'))
macro_props$Donor <- factor(macro_props$Donor, levels = c('HCa','HCb','HCc','UCa','UCb','UCc','CDa','CDb','CDc'))
ggplot(macro_props, aes(x = Donor, y = freq, fill = cell_type_level3,stratum = cell_type_level3, alluvium = cell_type_level3)) +
  geom_bar(stat = "identity", width = 0.6) +
  # facet_wrap(~ Group, scales = "free_x") +
  scale_fill_manual(values = level3_color) + 
  geom_flow(width=0.4,alpha=0.3,knot.pos=0) +
  labs(
    x = "Donor",
    y = "Proportion",
    fill = "Cell Type",
    title = "Cell type composition per Donor"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5),
        panel.grid.major = element_blank(),     # 删除主网格线
        panel.grid.minor = element_blank(),     # 删除次网格线
        panel.background = element_blank()      # 可选：移除背景
  )




#### umap --------------------------------

### 展示髓系细胞
cosmx_mye$SingleR2 <- factor(cosmx_mye$SingleR2, levels = c("Inflammatory monocytes","Macrophage NRG1","M0",        
                                                            "M1","M2","Mast","DCs","Eosinophils","N1",                      
                                                            "N2","N3","Cycling myeloid"))
DimPlot(cosmx_mye,
        cols = level3_color,
        group.by = 'SingleR2',
        raster = F,
        pt.size = 0.8
)+ theme_void() #size 5*4





cosmx_mye$cell_type_level3 <- factor(cosmx_mye$cell_type_level3, levels = c("Epithelial cell","Endothelium","Fibroblasts","Myofibroblasts","Inflammatory fibroblasts",
                                                                            "S1","S2a","S2b","S3","Pericytes",               
                                                                            "FRCs","Glia","Inflammatory monocytes","IL1B+NLRP3+_Macro","SPP1+CCL2+_Macro",        
                                                                            "Other Macro","Mast","DCs","Eosinophils","N1",                      
                                                                            "N2","N3","Cycling myeloid","T cell","B/Plasma cell"
                                                                            ))
DimPlot(cosmx_mye,
        cols = level3_color,
        group.by = 'cell_type_level3',
        raster = F,
        pt.size = 0.8
)+ theme_void()  # size 5*4



#### dotplot marker ----------------------

Idents(cosmx_macro) <- cosmx_macro$cell_arches
macro_marker_grouped <- list(
  Intestinal_Macro = c('CD68','CD163','MS4A4A','MRC1'),
  Complement_Macro = c('CCL13', 'C1QA', 'C1QB', 'C1QC'),
  Inflammasome_Macro = c('IL1B','NLRP3',),
  SPP1_Macro = c('SPP1','CCL2','HIF1A')
)
Cell_Type_Rank <- c( "Other Macro",  
                     "IL1B+NLRP3+_Macro",
                     "SPP1+CCL2+_Macro"
)
Cell_Type_Rank <- rev(Cell_Type_Rank)

cosmx_macro@active.ident <- factor(cosmx_macro@active.ident,
                                   levels = Cell_Type_Rank)
# cosmx_macro <- NormalizeData(cosmx_macro, normalization.method = "LogNormalize", scale.factor = 10000)


DotPlot(cosmx_macro, features = macro_marker_grouped) + 
  #RotatedAxis() + 
  theme(panel.border = element_rect(color = "black", size=0.8),
        panel.spacing = unit(2,"mm"),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 90,hjust = 1),
        panel.grid = element_line(color = "grey",linetype = 3,linewidth = 0.1),
        axis.title = element_blank()) +
  scale_color_gradientn(colours = c('#5C569E','#3D87AF','#68BEA2','#A5D5A2',
                                    '#EBF0B1','#F6C37E','#F18D5E','#DC5C52'),
                        limits = c(-2, 2),
                        oob = scales::squish
  ) 










#### neighbor vector normalized -------------------------

### load neighbor_matrix
## HCa
HCa_dir <- "/storage/data/KAI/Pan_Inflammation/Homo/STRNA/Colorectum/1.GSE234713(CosMx)/squidpy/GSM7473682/neighbor_matrix"

fov_dirs <- list.dirs(HCa_dir, recursive = FALSE)

# 初始化一个空的列表用于存储每个fov的矩阵
neighbor_list <- list()

# 遍历每个fov文件夹
for (fov_path in fov_dirs) {
  # 构造 neighbor matrix 文件路径
  neighbor_file <- list.files(fov_path, pattern = "*_neighbor_matrix.txt$", full.names = TRUE)
  
  if (length(neighbor_file) == 1) {
    neighbor_mat <- read.delim(neighbor_file, row.names = 1, check.names = FALSE)
    # neighbor_mat$fov <- basename(fov_path) # 添加一列记录fov来源
    
    neighbor_list[[basename(fov_path)]] <- neighbor_mat
  }
}
# 合并所有邻域矩阵（按行合并）
HCa.neighbor_matrix <- do.call(rbind, neighbor_list)

rownames(HCa.neighbor_matrix) <- sub("^fov\\d+\\.", "", rownames(HCa.neighbor_matrix))

write.table(HCa.neighbor_matrix, "/storage/data/KAI/Pan_Inflammation/Homo/STRNA/Colorectum/1.GSE234713(CosMx)/squidpy/GSM7473682/allmacro_neighbor/HCa.neighbor_matrix.txt", sep = '\t',quote = F,row.names = T,col.names = T)
HCa.neighbor_matrix <- read.table("/storage/data/KAI/Pan_Inflammation/Homo/STRNA/Colorectum/1.GSE234713(CosMx)/squidpy/GSM7473682/allmacro_neighbor/HCa.neighbor_matrix.txt", sep = '\t',header = T)

HCa.neighbor_matrix['HCa_15_6',]
## HCb
HCb_dir <- "/storage/data/KAI/Pan_Inflammation/Homo/STRNA/Colorectum/1.GSE234713(CosMx)/squidpy/GSM7473683/neighbor_matrix"

fov_dirs <- list.dirs(Hcb_dir, recursive = FALSE)

# 初始化一个空的列表用于存储每个fov的矩阵
neighbor_list <- list()

# 遍历每个fov文件夹
for (fov_path in fov_dirs) {
  # 构造 neighbor matrix 文件路径
  neighbor_file <- list.files(fov_path, pattern = "*_neighbor_matrix.txt$", full.names = TRUE)
  
  if (length(neighbor_file) == 1) {
    neighbor_mat <- read.delim(neighbor_file, row.names = 1, check.names = FALSE)
    # neighbor_mat$fov <- basename(fov_path) # 添加一列记录fov来源
    
    neighbor_list[[basename(fov_path)]] <- neighbor_mat
  }
}
# 合并所有邻域矩阵（按行合并）
HCb.neighbor_matrix <- do.call(rbind, neighbor_list)
rownames(HCb.neighbor_matrix) <- sub("^fov\\d+\\.", "", rownames(Hcb.neighbor_matrix))

write.table(HCb.neighbor_matrix, "/storage/data/KAI/Pan_Inflammation/Homo/STRNA/Colorectum/1.GSE234713(CosMx)/squidpy/GSM7473683/allmacro_neighbor/HCb.neighbor_matrix.txt", sep = '\t',quote = F,row.names = T,col.names = T)

HCb.neighbor_matrix <- read.table("/storage/data/KAI/Pan_Inflammation/Homo/STRNA/Colorectum/1.GSE234713(CosMx)/squidpy/GSM7473683/allmacro_neighbor/HCb.neighbor_matrix.txt", sep = '\t',header = T)



## HCc
HCc_dir <- "/storage/data/KAI/Pan_Inflammation/Homo/STRNA/Colorectum/1.GSE234713(CosMx)/squidpy/GSM7473684/neighbor_matrix"

fov_dirs <- list.dirs(HCc_dir, recursive = FALSE)

# 初始化一个空的列表用于存储每个fov的矩阵
neighbor_list <- list()

# 遍历每个fov文件夹
for (fov_path in fov_dirs) {
  # 构造 neighbor matrix 文件路径
  neighbor_file <- list.files(fov_path, pattern = "*_neighbor_matrix.txt$", full.names = TRUE)
  
  if (length(neighbor_file) == 1) {
    neighbor_mat <- read.delim(neighbor_file, row.names = 1, check.names = FALSE)
    # neighbor_mat$fov <- basename(fov_path) # 添加一列记录fov来源
    
    neighbor_list[[basename(fov_path)]] <- neighbor_mat
  }
}
# 合并所有邻域矩阵（按行合并）
HCc.neighbor_matrix <- do.call(rbind, neighbor_list)
rownames(HCc.neighbor_matrix) <- sub("^fov\\d+\\.", "", rownames(HCc.neighbor_matrix))

write.table(HCc.neighbor_matrix, "/storage/data/KAI/Pan_Inflammation/Homo/STRNA/Colorectum/1.GSE234713(CosMx)/squidpy/GSM7473684/allmacro_neighbor/HCc.neighbor_matrix.txt", sep = '\t',quote = F,row.names = T,col.names = T)



## UCa
UCa_dir <- "/storage/data/KAI/Pan_Inflammation/Homo/STRNA/Colorectum/1.GSE234713(CosMx)/squidpy/GSM7473685/neighbor_matrix"

fov_dirs <- list.dirs(UCa_dir, recursive = FALSE)

# 初始化一个空的列表用于存储每个fov的矩阵
neighbor_list <- list()

# 遍历每个fov文件夹
for (fov_path in fov_dirs) {
  # 构造 neighbor matrix 文件路径
  neighbor_file <- list.files(fov_path, pattern = "*_neighbor_matrix.txt$", full.names = TRUE)
  
  if (length(neighbor_file) == 1) {
    neighbor_mat <- read.delim(neighbor_file, row.names = 1, check.names = FALSE)
    # neighbor_mat$fov <- basename(fov_path) # 添加一列记录fov来源
    
    neighbor_list[[basename(fov_path)]] <- neighbor_mat
  }
}
# 合并所有邻域矩阵（按行合并）
UCa.neighbor_matrix <- do.call(rbind, neighbor_list)
rownames(UCa.neighbor_matrix) <- sub("^fov\\d+\\.", "", rownames(UCa.neighbor_matrix))

write.table(UCa.neighbor_matrix, "/storage/data/KAI/Pan_Inflammation/Homo/STRNA/Colorectum/1.GSE234713(CosMx)/squidpy/GSM7473685/allmacro_neighbor/UCa.neighbor_matrix.txt", sep = '\t',quote = F,row.names = T,col.names = T)



## UCb
UCb_dir <- "/storage/data/KAI/Pan_Inflammation/Homo/STRNA/Colorectum/1.GSE234713(CosMx)/squidpy/GSM7473686/neighbor_matrix"

fov_dirs <- list.dirs(UCb_dir, recursive = FALSE)

# 初始化一个空的列表用于存储每个fov的矩阵
neighbor_list <- list()

# 遍历每个fov文件夹
for (fov_path in fov_dirs) {
  # 构造 neighbor matrix 文件路径
  neighbor_file <- list.files(fov_path, pattern = "*_neighbor_matrix.txt$", full.names = TRUE)
  
  if (length(neighbor_file) == 1) {
    neighbor_mat <- read.delim(neighbor_file, row.names = 1, check.names = FALSE)
    # neighbor_mat$fov <- basename(fov_path) # 添加一列记录fov来源
    
    neighbor_list[[basename(fov_path)]] <- neighbor_mat
  }
}
# 合并所有邻域矩阵（按行合并）
UCb.neighbor_matrix <- do.call(rbind, neighbor_list)
rownames(UCb.neighbor_matrix) <- sub("^fov\\d+\\.", "", rownames(UCb.neighbor_matrix))
write.table(UCb.neighbor_matrix, "/storage/data/KAI/Pan_Inflammation/Homo/STRNA/Colorectum/1.GSE234713(CosMx)/squidpy/GSM7473686/allmacro_neighbor/UCb.neighbor_matrix.txt", sep = '\t',quote = F,row.names = T,col.names = T)





## UCc
UCc_dir <- "/storage/data/KAI/Pan_Inflammation/Homo/STRNA/Colorectum/1.GSE234713(CosMx)/squidpy/GSM7473687/neighbor_matrix"

fov_dirs <- list.dirs(UCc_dir, recursive = FALSE)

# 初始化一个空的列表用于存储每个fov的矩阵
neighbor_list <- list()
# 遍历每个fov文件夹
for (fov_path in fov_dirs) {
  # 构造 neighbor matrix 文件路径
  neighbor_file <- list.files(fov_path, pattern = "*_neighbor_matrix.txt$", full.names = TRUE)
  
  if (length(neighbor_file) == 1) {
    neighbor_mat <- read.delim(neighbor_file, row.names = 1, check.names = FALSE)
    # neighbor_mat$fov <- basename(fov_path) # 添加一列记录fov来源
    
    neighbor_list[[basename(fov_path)]] <- neighbor_mat
  }
}
# 合并所有邻域矩阵（按行合并）
UCc.neighbor_matrix <- do.call(rbind, neighbor_list)
rownames(UCc.neighbor_matrix) <- sub("^fov\\d+\\.", "", rownames(UCc.neighbor_matrix))
write.table(UCc.neighbor_matrix, "/storage/data/KAI/Pan_Inflammation/Homo/STRNA/Colorectum/1.GSE234713(CosMx)/squidpy/GSM7473687/allmacro_neighbor/UCc.neighbor_matrix.txt", sep = '\t',quote = F,row.names = T,col.names = T)



## CDa
CDa_dir <- "/storage/data/KAI/Pan_Inflammation/Homo/STRNA/Colorectum/1.GSE234713(CosMx)/squidpy/GSM7473688/neighbor_matrix"

fov_dirs <- list.dirs(CDa_dir, recursive = FALSE)

# 初始化一个空的列表用于存储每个fov的矩阵
neighbor_list <- list()
# 遍历每个fov文件夹
for (fov_path in fov_dirs) {
  # 构造 neighbor matrix 文件路径
  neighbor_file <- list.files(fov_path, pattern = "*_neighbor_matrix.txt$", full.names = TRUE)
  
  if (length(neighbor_file) == 1) {
    neighbor_mat <- read.delim(neighbor_file, row.names = 1, check.names = FALSE)
    # neighbor_mat$fov <- basename(fov_path) # 添加一列记录fov来源
    
    neighbor_list[[basename(fov_path)]] <- neighbor_mat
  }
}

# 合并所有邻域矩阵（按行合并）
CDa.neighbor_matrix <- do.call(rbind, neighbor_list)
rownames(CDa.neighbor_matrix) <- sub("^fov\\d+\\.", "", rownames(CDa.neighbor_matrix))
write.table(CDa.neighbor_matrix, "/storage/data/KAI/Pan_Inflammation/Homo/STRNA/Colorectum/1.GSE234713(CosMx)/squidpy/GSM7473688/allmacro_neighbor/CDa.neighbor_matrix.txt", sep = '\t',quote = F,row.names = T,col.names = T)


## CDb
CDb_dir <- "/storage/data/KAI/Pan_Inflammation/Homo/STRNA/Colorectum/1.GSE234713(CosMx)/squidpy/GSM7473689/neighbor_matrix"

fov_dirs <- list.dirs(CDb_dir, recursive = FALSE)

# 初始化一个空的列表用于存储每个fov的矩阵
neighbor_list <- list()
# 遍历每个fov文件夹
for (fov_path in fov_dirs) {
  # 构造 neighbor matrix 文件路径
  neighbor_file <- list.files(fov_path, pattern = "*_neighbor_matrix.txt$", full.names = TRUE)
  
  if (length(neighbor_file) == 1) {
    neighbor_mat <- read.delim(neighbor_file, row.names = 1, check.names = FALSE)
    # neighbor_mat$fov <- basename(fov_path) # 添加一列记录fov来源
    
    neighbor_list[[basename(fov_path)]] <- neighbor_mat
  }
}

# 合并所有邻域矩阵（按行合并）
CDb.neighbor_matrix <- do.call(rbind, neighbor_list)
rownames(CDb.neighbor_matrix) <- sub("^fov\\d+\\.", "", rownames(CDb.neighbor_matrix))
write.table(CDb.neighbor_matrix, "/storage/data/KAI/Pan_Inflammation/Homo/STRNA/Colorectum/1.GSE234713(CosMx)/squidpy/GSM7473689/allmacro_neighbor/CDb.neighbor_matrix.txt", sep = '\t',quote = F,row.names = T,col.names = T)



## CDc
CDc_dir <- "/storage/data/KAI/Pan_Inflammation/Homo/STRNA/Colorectum/1.GSE234713(CosMx)/squidpy/GSM7473690/neighbor_matrix"

fov_dirs <- list.dirs(CDc_dir, recursive = FALSE)

# 初始化一个空的列表用于存储每个fov的矩阵
neighbor_list <- list()
# 遍历每个fov文件夹
for (fov_path in fov_dirs) {
  # 构造 neighbor matrix 文件路径
  neighbor_file <- list.files(fov_path, pattern = "*_neighbor_matrix.txt$", full.names = TRUE)
  
  if (length(neighbor_file) == 1) {
    neighbor_mat <- read.delim(neighbor_file, row.names = 1, check.names = FALSE)
    # neighbor_mat$fov <- basename(fov_path) # 添加一列记录fov来源
    
    neighbor_list[[basename(fov_path)]] <- neighbor_mat
  }
}

# 合并所有邻域矩阵（按行合并）
CDc.neighbor_matrix <- do.call(rbind, neighbor_list)
rownames(CDc.neighbor_matrix) <- sub("^fov\\d+\\.", "", rownames(CDc.neighbor_matrix))
write.table(CDc.neighbor_matrix, "/storage/data/KAI/Pan_Inflammation/Homo/STRNA/Colorectum/1.GSE234713(CosMx)/squidpy/GSM7473690/allmacro_neighbor/CDc.neighbor_matrix.txt", sep = '\t',quote = F,row.names = T,col.names = T)


### 批量对每个样本的neighbor_matrix转换为neighbor_prop后进行标准化

donor_map <- c(
  "HCa" = "GSM7473682", "HCb" = "GSM7473683", "HCc" = "GSM7473684",
  "UCa" = "GSM7473685", "UCb" = "GSM7473686", "UCc" = "GSM7473687",
  "CDa" = "GSM7473688", "CDb" = "GSM7473689", "CDc" = "GSM7473690"
)

cell_types <- sort(unique(cosmx_all$cell_type_level3))
normalized_list <- list()

for (donor in names(donor_map)) {
  gsm_id <- donor_map[[donor]]
  file_path <- paste0("/storage/data/KAI/Pan_Inflammation/Homo/STRNA/Colorectum/1.GSE234713(CosMx)/squidpy/",
                      gsm_id, "/allmacro_neighbor/", donor, ".neighbor_matrix.txt")
  # 读取邻域计数矩阵
  neighbor_matrix <- read.table(file_path, sep = "\t", header = TRUE, row.names = 1,check.names = F)
  
  # 确保列顺序一致
  neighbor_matrix <- neighbor_matrix[, cell_types %>% as.character()]
  
  # Step 1: 转为邻域比例（行归一化）
  neighbor_prop <- neighbor_matrix / rowSums(neighbor_matrix)
  
  # Step 2: 计算 donor 的组织整体组成比例
  donor_subset <- cosmx_all@meta.data %>% filter(Donor == donor)
  donor_counts <- table(factor(donor_subset$cell_type_level3, levels = cell_types))
  donor_prop <- donor_counts / sum(donor_counts)
  
  # Step 3: 标准化（邻域比例 / 组织比例）
  donor_norm <- sweep(neighbor_prop, 2, donor_prop, FUN = "/")
  
  # 存入列表
  normalized_list[[donor]] <- donor_norm
}

all_normalized_neighbor <- do.call(rbind, normalized_list)
all_normalized_neighbor[is.na(all_normalized_neighbor)] <- 0

rownames(all_normalized_neighbor) <- sub("^\\w+\\.", "", rownames(all_normalized_neighbor))

write.table(all_normalized_neighbor,"/storage/data/KAI/Pan_Inflammation/Homo/STRNA/Colorectum/1.GSE234713(CosMx)/squidpy/macro_80umneighbor/IBDCosMX_normalized_80umneighbor.txt",sep = "\t",quote = F, row.names = T, col.names = T)

all_normalized_neighbor <- read.table("/storage/data/KAI/Pan_Inflammation/Homo/STRNA/Colorectum/1.GSE234713(CosMx)/squidpy/macro_80umneighbor/IBDCosMX_normalized_80umneighbor.txt",sep = "\t",header = T, check.names = F)


#### NMF identify spatial macrophage subtypes -------------------------------
### NMF
# k= 5, CancerCell Linghua Wang etal.
nmf_result <- nmf(all_normalized_neighbor, rank = 5, method = "brunet", nrun = 30, seed = 123)

# > which(rowSums(all_normalized_neighbor)==0)
# HCa_156_13   HCa_15_6 UCa_926_15 
# 414       1333       6398 

filtered_matrix <- all_normalized_neighbor[rowSums(all_normalized_neighbor != 0, na.rm = TRUE) != 0, ] # 去除3个空的邻域
estimate <- nmf(filtered_matrix, 2:10, nrun = 5)


estimate <- readRDS("/storage/data/KAI/Pan_Inflammation/Homo/STRNA/Colorectum/1.GSE234713(CosMx)/squidpy/macro_80umneighbor/NMF/nmf_k2-10.RDS")

plot(estimate)

nmf_result <- nmf(filtered_matrix, rank = 6, method = "brunet", nrun = 50, seed = 110)

saveRDS(nmf_result,"/storage/data/KAI/Pan_Inflammation/Homo/STRNA/Colorectum/1.GSE234713(CosMx)/squidpy/macro_80umneighbor/NMF/k=6/NMF_result.RDS")
nmf_result <- readRDS("/storage/data/KAI/Pan_Inflammation/Homo/STRNA/Colorectum/1.GSE234713(CosMx)/squidpy/macro_80umneighbor/NMF/k=6/NMF_result.RDS")


W <- basis(nmf_result)
H <- coef(nmf_result)

### Seurat Leiden
W_seurat <- CreateSeuratObject(counts = t(W))  # 注意转置，行为因子，列为细胞
W_seurat@assays$RNA$data <- W_seurat@assays$RNA$counts
# 归一化和降维
W_seurat <- ScaleData(W_seurat)
W_seurat <- RunPCA(W_seurat, features = rownames(W_seurat))
W_seurat <- FindNeighbors(W_seurat, dims = 1:5, k.param = 30, prune.SNN = 0)
use_condaenv("/storage/data/KAI/miniconda3/envs/sceasy",required = T)

W_seurat <- FindClusters(W_seurat, resolution = 0.1, algorithm = 4)  # Leiden resolution = 0.1, CancerCell Linghua Wang etal.

# 加入巨噬细胞注释信息
macro_label <- read.table('/storage/data/KAI/Pan_Inflammation/Homo/STRNA/Colorectum/1.GSE234713(CosMx)/squidpy/macro_80umneighbor/Niches/macro_level3.txt',sep='\t',header = T,row.names = 1,check.names = F)
lbl <- setNames(as.character(macro_label$cell_type_level3), rownames(macro_label))
W_seurat <- AddMetaData(W_seurat, metadata = lbl, col.name = "cell_type_level3")


# umap可视化
W_seurat <- RunUMAP(W_seurat, dims = 1:5, reduction = "pca", n.neighbors = 30, min.dist = 0.3)


DimPlot(W_seurat, reduction = "umap",
        group.by = "RNA_snn_res.0.1", label = TRUE, repel = TRUE,pt.size = 0.5
        ) +
  NoLegend()

DimPlot(W_seurat, reduction = "umap",
        group.by = "cell_type_level3", label = TRUE, repel = TRUE,pt.size = 0.5,cols = level3_color
) +
  NoLegend()


df <- FetchData(W_seurat, vars = c("RNA_snn_res.0.1", "cell_type_level3")) %>%
  tibble::as_tibble(rownames = "barcode") %>%
  transmute(
    barcode,
    cluster   = as.character(`RNA_snn_res.0.1`),
    cell_type = as.character(cell_type_level3)
  ) %>%
  filter(!is.na(cell_type))

prop_tbl <- df %>%
  dplyr::count(cluster, cell_type, name = "n") %>%   # 显式用 dplyr::count
  group_by(cluster) %>%
  mutate(total = sum(n), prop = n/total, percent = round(prop*100, 1)) %>%
  ungroup()

W_meta <- data.frame(
  row.names = Cells(W_seurat),
  leiden_cluster = W_seurat$RNA_snn_res.0.1,
  cell_type_level3 = W_seurat$cell_type_level3
)

write.table(W_meta, file = "/storage/data/KAI/Pan_Inflammation/Homo/STRNA/Colorectum/1.GSE234713(CosMx)/squidpy/macro_80umneighbor/NMF/k=6/leiden_meta.txt",
            sep = "\t", quote = FALSE, row.names = T)



