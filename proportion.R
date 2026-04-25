### Package ------------------------------
library(dplyr)
library(ggplot2)
library(readxl)
library(patchwork)
library(tidyr)

library(Startrac)
library(ggplot2)
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


library(lsa)
library(proxy)


### load data -------------------------------
Mac <- readRDS("/storage/data/KAI/Pan_Inflammation/Homo/Mac2025/scVI/new/pan_Mac_filtered.RDS") 


### global ----------------------------------
mac_color = c(
  "Synovial resident Macro"="#F7BC99",
  "Alveolar Macro"="#79B472",
  "Intestinal resident Macro"="#B79981",
  "Kupffer cell"="#553E2E",
  "Langerhans cell" = "#99A4BC",
  "Uterine resident Macro"="#E2CECE",
  "IL1B+NLRP3+_Macro"="#A77980",
  "MMP3+CXCL8+_Macro"="#BF94B2",
  "CCL13+_Complement-associated_Macro" = "#4d75be",
  "SPP1+CCL2+_Macro"="#D9D19B",
  "CD1C+_DC-like_Macro"="#83978C"
)


## dotplot 分组展示
DimPlot(object = Mac, 
        #split.by = 'Inflammation', 
        group.by = 'cell_type_level3',
        cols = mac_color,
        raster = F
)

Idents(Mac) <- Mac$cell_type_level3
DEG_Mac <- FindAllMarkers(Mac, log2FC.threshold = 1, test.use = "wilcox",
                                min.pct = 0.25, min.diff.pct = 0.1, only.pos = TRUE,
                                assay = "RNA")

DEG_Mac[DEG_Mac$cluster=='SPP1+_Macro',] %>% arrange(desc(avg_log2FC)) %>%        # 按log2FC降序排序
  slice_head(n = 50) %>%               # 每组取前50个
  ungroup()

VEGFA_marker <- FindMarkers(Mac, ident.1 = 'VEGFA+_Proangiogenic_Macro', logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
VEGFA_t50 <- VEGFA_marker %>% 
  arrange(desc(avg_diff)) %>%        # 按log2FC降序排序
  slice_head(n = 50) %>%               # 每组取前50个
  ungroup()


F13A1_marker <- FindMarkers(Mac, ident.1 = 'F13A1+_Macro', logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

F13A1_t50 <- F13A1_marker %>% 
  arrange(desc(avg_log2FC)) %>%        # 按log2FC降序排序
  slice_head(n = 50) %>%               # 每组取前50个
  ungroup()

top50_per_cluster <- DEG_Mac %>%
  group_by(cluster) %>%                # 按簇分组
  arrange(desc(avg_log2FC)) %>%        # 按log2FC降序排序
  slice_head(n = 50) %>%               # 每组取前50个
  ungroup()

################################################# Cosine Similarity ######################################################

Mac$cell_type_level3 %>% table()


IL1B_Macro <- subset(Mac, subset = cell_type_level3 == "IL1B+NLRP3+_Macro")
IL1B_Macro$leiden_res_0.8 <- droplevels(IL1B_Macro$leiden_res_0.8)

Alveolar_Macro <- subset(Mac, subset = cell_type_level3 == "Alveolar Macro")
Alveolar_Macro$leiden_res_0.8 <- droplevels(Alveolar_Macro$leiden_res_0.8)

Kupffer_cell <- subset(Mac, subset = cell_type_level3 == "Kupffer cell")
Kupffer_cell$leiden_res_0.8 <- droplevels(Kupffer_cell$leiden_res_0.8)

Langerhans_cell <- subset(Mac, subset = cell_type_level3 == "Langerhans cell")
Langerhans_cell$leiden_res_0.8 <- droplevels(Langerhans_cell$leiden_res_0.8)

Synovial_Macro <- subset(Mac, subset = cell_type_level3 == "Synovial resident Macro")
Synovial_Macro$leiden_res_0.8 <- droplevels(Synovial_Macro$leiden_res_0.8)

SPP1_Macro <- subset(Mac, subset = cell_type_level3 == "SPP1+CCL2+_Macro")
SPP1_Macro$leiden_res_0.8 <- droplevels(SPP1_Macro$leiden_res_0.8)

MMP3_Macro <- subset(Mac, subset = cell_type_level3 == "MMP3+CXCL8+_Macro")
MMP3_Macro$leiden_res_0.8 <- droplevels(MMP3_Macro$leiden_res_0.8)

CD1C_Macro <- subset(Mac, subset = cell_type_level3 == "CD1C+_DC-like_Macro")
CD1C_Macro$leiden_res_0.8 <- droplevels(CD1C_Macro$leiden_res_0.8)

CCL13_Macro <- subset(Mac, subset = cell_type_level3 == "CCL13+_Complement-associated_Macro")
CCL13_Macro$leiden_res_0.8 <- droplevels(CCL13_Macro$leiden_res_0.8)

Macro_list <- c('IL1B_Macro','Alveolar_Macro','Kupffer_cell',
                'Langerhans_cell','Synovial_Macro','SPP1_Macro',
                'MMP3_Macro','CD1C_Macro','CCL13_Macro'
                )

## Try1
expr_matrix <- GetAssayData(IL1B_Macro, layer = "counts")
group_info <- IL1B_Macro$leiden_res_0.8
avg_expr <- aggregate(t(as.matrix(expr_matrix)), by = list(group_info), FUN = mean)
avg_expr <- na.omit(avg_expr) # 行是分类，列是基因
rownames(avg_expr) <- avg_expr$Group.1
avg_expr <- avg_expr[, -1]
cos_sim <- proxy::simil(as.matrix(avg_expr), method = "cosine")
cossim_matrix <- as.matrix(cos_sim)
diag(cossim_matrix) <- 1




for (i in 1:length(Macro_list)) {
  temp_seu <- get(Macro_list[i])
  expr_matrix <- GetAssayData(temp_seu, layer = "counts")
  group_info <- temp_seu$leiden_res_0.8
  avg_expr <- aggregate(t(as.matrix(expr_matrix)), by = list(group_info), FUN = mean)
  avg_expr <- na.omit(avg_expr) # 行是分类，列是基因
  rownames(avg_expr) <- avg_expr$Group.1
  avg_expr <- avg_expr[, -1]
  cos_sim <- proxy::simil(as.matrix(avg_expr), method = "cosine")
  cossim_matrix <- as.matrix(cos_sim)
  diag(cossim_matrix) <- 1
  assign(paste0(Macro_list[i],'_cosmatrix'), cossim_matrix)
}


expr_matrix <- GetAssayData(Mac, layer = "counts")
group_info <- IL1B_Macro$leiden_res_0.8
avg_expr <- aggregate(t(as.matrix(expr_matrix)), by = list(group_info), FUN = mean)
avg_expr <- na.omit(avg_expr) # 行是分类，列是基因
rownames(avg_expr) <- avg_expr$Group.1
avg_expr <- avg_expr[, -1]
cos_sim <- proxy::simil(as.matrix(avg_expr), method = "cosine")
cossim_matrix <- as.matrix(cos_sim)
diag(cossim_matrix) <- 1
















### Pearson correlation ---------------------
hvgs <- VariableFeatures(Mac)
expression_matrix <- as.matrix(Mac@assays$RNA@data[hvgs, ])
clusters <- Mac$leiden_res_0.8
# 按簇计算平均表达值
cluster_means <- t(aggregate(t(expression_matrix), by = list(clusters), FUN = mean, na.rm = TRUE))

colnames(cluster_means) <- cluster_means[1, ]  # 第一行是簇名称
cluster_means <- cluster_means[-1, ]  # 移除第一行（簇名称）
cluster_means <- apply(cluster_means, 2, as.numeric)  # 转换为数值矩阵
rownames(cluster_means) <- rownames(expression_matrix)  # 设置基因名称

correlation_matrix <- cor(cluster_means, method = "pearson") # Pearson correlation


blue_gradient <- colorRampPalette(c("darkblue", "#79a7e8"))(50)  # 0 到 0.85
yellow_color <- rep("#ffc6a2", 10)  # 0.85 到 0.9，深黄色
red_gradient <- colorRampPalette(c("#e8798a", "darkred"))(40)  # 0.9 到 1
custom_colors <- c(blue_gradient, yellow_color, red_gradient)

pdf(file = "/storage/data/KAI/Pan_Inflammation/Homo/Mac2025/analysis/statistics/pearson/result/heatmap_mono_mac2025_pearson.pdf",
    width = 10, height = 8)  # 调整宽度和高度
pheatmap(
  correlation_matrix,
  color = custom_colors,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  scale = "none",  # 不需要额外标准化
  main = "Pearson Correlation Between Clusters (leiden_res_0.8)",
  fontsize = 10,
  fontsize_row = 8,
  fontsize_col = 8,
  border_color = NA,
  show_rownames = TRUE,
  show_colnames = TRUE
  #filename = "/storage/data/KAI/Pan_Inflammation/Homo/Mac2025/analysis/statistics/pearson/result/heatmap_mac2025_pearson.pdf"  # 保存为 PDF
)
dev.off()

saveRDS(correlation_matrix,"/storage/data/KAI/Pan_Inflammation/Homo/Mac2025/analysis/statistics/pearson/data/pearson_Mac_res0.8.RDS")


################################# Macrophage boxplot by Organ (sample) ###################################
Mac_meta_sample <- Mac@meta.data[, c("Organ", "Inflammation", "cell_type_level3","Sample_geo_accession")]
total_cells_sample <- Mac_meta_sample %>%
  group_by(Organ, Inflammation, Sample_geo_accession) %>%
  summarise(total = n(), .groups = "drop")

cell_type_counts_sample <- Mac_meta_sample %>%
  group_by(cell_type_level3, Organ, Inflammation, Sample_geo_accession) %>%
  summarise(count = n(), .groups = "drop")

Mac_proportion_sample <- cell_type_counts_sample %>%
  left_join(total_cells_sample, by = c("Organ", "Inflammation", "Sample_geo_accession")) %>%
  mutate(proportion = count / total) %>%
  select(cell_type_level3, Organ, Inflammation, Sample_geo_accession, proportion)

# 去除组织特异的组织驻留巨噬细胞
Mac_proportion_sample <- Mac_proportion_sample[which(Mac_proportion_sample$cell_type_level3 != 'Intestinal resident Macro'),]
Mac_proportion_sample <- Mac_proportion_sample[which(Mac_proportion_sample$cell_type_level3 != 'Uterine resident Macro'),]
Mac_proportion_sample <- Mac_proportion_sample[which(Mac_proportion_sample$cell_type_level3 != 'Alveolar Macro'),]
Mac_proportion_sample <- Mac_proportion_sample[which(Mac_proportion_sample$cell_type_level3 != 'Synovial resident Macro'),]
Mac_proportion_sample <- Mac_proportion_sample[which(Mac_proportion_sample$cell_type_level3 != 'Kupffer cell'),]
Mac_proportion_sample <- Mac_proportion_sample[which(Mac_proportion_sample$cell_type_level3 != 'Langerhans cell'),]




Mac_proportion_sample_boxplot  <- ggplot(Mac_proportion_sample, aes(x = Inflammation, y = proportion)) +
  geom_point(aes(color = Organ),position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.4), size = 0.5) +
  #geom_jitter(aes(color = Organ), width = 0.3, size = 0.5, alpha = 0.3) +
  geom_boxplot(aes(y = proportion), position = position_dodge(width = 0.3), alpha = 0, outlier.shape = NA) +  
  
  facet_wrap(~ cell_type_level3, nrow = 1) +
  theme_classic() + # theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 12),
        axis.text.y = element_text(color = "black", size = 12),
        plot.title = element_text(size = 14, face = "bold"),
        strip.text = element_text(size = 10, face = "bold")) +
  labs(title = "Macrphage cell by Organ (Sample)",
       #x = "Inflammation",
       y = "Proportion of Macrphage",
       color = "Organ") +
  scale_color_brewer(palette = "Set1")

# wilcoxon Test 
wc_results <- Mac_proportion_sample %>%
  group_by(cell_type_level3) %>%
  summarise(p_value = wilcox.test(proportion[Inflammation == "Inflamed"], 
                                  proportion[Inflammation == "Normal"],alternative = "greater",exact = T)$p.value, .groups = "drop") %>%
  mutate(significance = case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01 ~ "**",
    p_value < 0.05 ~ "*",
    TRUE ~ "ns"  # not significant
  ))
Mac_proportion_sample_boxplot + geom_text(data = wc_results, 
                                            aes(x = 1.5, y = max(Mac_proportion_sample$proportion) + 0.1, 
                                                label = significance),
                                            inherit.aes = FALSE, size = 4) +
  facet_wrap(~ cell_type_level3, nrow = 1)

################################# Macrophage proportion by Organ ###################################
# 计算构成比例
organ_composition <- cell_type_counts_sample %>%
  group_by(cell_type_level3, Organ) %>%
  summarise(total_count = sum(count), .groups = "drop") %>%
  group_by(cell_type_level3) %>%
  mutate(proportion = total_count / sum(total_count)) %>%
  ungroup()

# 绘制堆叠条形图，每一行是不同的 cell_type_level3
ggplot(organ_composition, aes(x = cell_type_level3, y = proportion, fill = Organ)) +
  geom_bar(stat = "identity", position = "fill") +  # 横向堆叠
  coord_flip() +  # 翻转坐标轴
  labs(x = "Cell Type", y = "Proportion", fill = "Organ",
       title = "Organ Composition in Different Cell Types") +
  scale_fill_brewer(palette = "Set1") + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1),  # 旋转 x 轴标签
        strip.text = element_text(face = "bold"),
        legend.position = "bottom")

################################ Ro/e ##################################
Mac_roe <- Mac@meta.data[, c("Organ", "Inflammation", "cell_type_level3","Sample_geo_accession")] %>% group_by(cell_type_level3, Organ, Inflammation) %>% summarise(total = n(), .groups = "drop")


Mac_repeated <- Mac_roe %>% uncount(total)# 为每行按 total 列的值重复行

# 去除组织特异的组织驻留巨噬细胞
Mac_repeated <- Mac_repeated[which(Mac_repeated$cell_type_level3 != 'Intestinal resident Macro'),]
Mac_repeated <- Mac_repeated[which(Mac_repeated$cell_type_level3 != 'Uterine resident Macro'),]
Mac_repeated <- Mac_repeated[which(Mac_repeated$cell_type_level3 != 'Alveolar Macro'),]
Mac_repeated <- Mac_repeated[which(Mac_repeated$cell_type_level3 != 'Synovial resident Macro'),]
Mac_repeated <- Mac_repeated[which(Mac_repeated$cell_type_level3 != 'Kupffer cell'),]
Mac_repeated <- Mac_repeated[which(Mac_repeated$cell_type_level3 != 'Langerhans cell'),]



Roe <- calTissueDist(Mac_repeated,
                     byPatient = F,
                     colname.cluster = "cell_type_level3", # 不同细胞亚群
                     colname.patient = "Organ", # 不同样本
                     colname.tissue = "Inflammation", # 不同组织
                     method = "chisq", # "chisq", "fisher", and "freq" 
                     min.rowSum = 0)

Roe <- Roe[c('CCL13+_Complement-associated_Macro','CD1C+_DC-like_Macro','IL1B+NLRP3+_Macro','MMP3+CXCL8+_Macro','SPP1+CCL2+_Macro'),]














