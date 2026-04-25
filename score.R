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


library(AUCell)
library(upstartr)
library(multcompView)

### load data -------------------------------
Mac <- readRDS("/storage/data/KAI/Pan_Inflammation/Homo/Mac2025/scVI/new/pan_Mac_counts.RDS") 

##为moMac 保存counts矩阵
moMac <- readRDS("/storage/data/KAI/Pan_Inflammation/Homo/Mac2025/scVI/new/mono_mac.RDS")

IL1BSPP1_Mac <- subset(moMac, subset = cell_type_level3 %in% c("IL1B+NLRP3+_Macro", "SPP1+CCL2+_Macro"))


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
momac_color = c(
  "IL1B+NLRP3+_Macro"="#A77980",
  "MMP3+CXCL8+_Macro"="#BF94B2",
  "CCL13+_Complement-associated_Macro" = "#4d75be",
  "SPP1+CCL2+_Macro"="#D9D19B",
  "CD1C+_DC-like_Macro"="#83978C"
)

Mac_color <- c(
  "IL1B+NLRP3+_Macro"="#FA8072",
  "SPP1+CCL2+_Macro"="#ff9d5c",
  "CCL13+_Complement-associated_Macro"="#B6D0E2",
  "MMP3+CXCL8+_Macro"="#FFC5BF",
  "CD1C+_DC-like_Macro"="#53738c",
  "Alveolar Macro"="#d4c2db",
  "Intestinal resident Macro"="#5599C8",
  "Uterine resident Macro"="#FFBF00",
  "Synovial resident Macro"="#bc9a7f",
  "Kupffer cell"="#CAA7DD",
  "Langerhans cell"="#97C1A9"
)


########################################## T-attraction/suppression AUCell #####################################################
### load genesets 
## 2023.01_Nature communications_A comprehensive single-cell map of T cell exhaustion-associated immune environments in human breast cancer
T_att <- c('CCL21','CCL17','CCL2','CXCL9','CXCL10','CXCL11','CXCL12','CCL3','CCL4','CCL5','CXCL16')

T_supp <- c('IDO1','CD274','PDCD1LG2','TNFSF10','HLA-E','HLA-G','VTCN1','IL10','TGFB1','TGFB2','PTGS2','PTGES', 'LGALS9','CCL22','CD80','CD86')

geneSets <- list(
  T_attraction = T_att,
  T_suppression = T_supp
)

### All Mac AUCell score
Mac_rankings <- AUCell_buildRankings(GetAssayData(Mac, layer = "counts"))

Mac_AUC <- AUCell_calcAUC(geneSets, Mac_rankings)
# 通过 AUC 评分直方图 选择一个合适的阈值，用于区分基因集活跃的细胞与不活跃的细胞。
Mac_assignment <- AUCell_exploreThresholds(Mac_AUC, plotHist = TRUE, assign=TRUE)
# thr阈值用于筛选出 AUC 评分高于阈值的细胞，表示它们在该基因集上是“活跃”的。
thr_T_att <-  Mac_assignment$T_attraction$aucThr$selected
thr_T_supp <-  Mac_assignment$T_suppression$aucThr$selected

### 结果存入seurat对象 Mac中
all(colnames(Mac) == colnames(Mac_AUC)) # TRUE

T_scores <- as.data.frame(t(getAUC(Mac_AUC)))
Mac@meta.data <- cbind(Mac@meta.data, T_scores)

### 按cell_type_level3中各巨噬细胞亚型计算AUCell均值/中值

Mac_level3_means <- Mac@meta.data %>%
  group_by(cell_type_level3) %>%
  summarise(
    T_attraction_mean = mean(T_attraction, na.rm = TRUE),
    T_suppression_mean = mean(T_suppression, na.rm = TRUE)
  )

Mac_level3_medians <- Mac@meta.data %>%
  group_by(cell_type_level3) %>%
  summarise(
    T_attraction_median = median(T_attraction, na.rm = TRUE),
    T_suppression_median = median(T_suppression, na.rm = TRUE)
  )


moMac_level3_means <- Mac_level3_means[Mac_level3_means$cell_type_level3 %in% c('CCL13+_Complement-associated_Macro','MMP3+CXCL8+_Macro','CD1C+_DC-like_Macro','IL1B+NLRP3+_Macro','SPP1+CCL2+_Macro'),]
moMac_level3_means_scale <- moMac_level3_means %>%
  mutate(
    T_attraction_scaled = scale(T_attraction_mean),
    T_suppression_scaled = scale(T_suppression_mean)
  )

### 可视化
ggplot(Mac_level3_means, aes(x = T_attraction_mean, y = T_suppression_mean, label = cell_type_level3)) +
  geom_point(aes(color = cell_type_level3), size = 4) +  # 颜色区分不同细胞亚型
  geom_text(vjust = -0.5, hjust = 0.5, size = 4) +  # 添加细胞亚型标签
  theme_minimal() +
  labs(
    x = "T-cell attraction score (mean)",
    y = "T-cell suppression score (mean)"
  ) +
  theme(legend.position = "none")  # 移除图例

## 标准化后可视化
ggplot(moMac_level3_means_scale, aes(x = T_attraction_scaled, y = T_suppression_scaled, color = cell_type_level3, label = cell_type_level3)) +
  geom_point(size = 8) +  # 调整点的大小
  geom_text(vjust = -1, hjust = 0.5, size = 5) +
  scale_color_manual(values = momac_color) +  # 自定义颜色
  labs(
    x = "T-cell attraction score (mean)",
    y = "T-cell suppression score (mean)",
    title = "T-cell Associated Signatures (Scale)"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", size = 1.5, fill = NA),
    legend.title = element_blank(),
    legend.position = 'none',
    text = element_text(size = 14)
  )+
  coord_fixed(ratio = 1)





##优先使用mean
ggplot(Mac_level3_medians, aes(x = T_attraction_median, y = T_suppression_median, label = cell_type_level3)) +
  geom_point(aes(color = cell_type_level3), size = 4) +  # 颜色区分不同细胞亚型
  geom_text(vjust = -0.5, hjust = 0.5, size = 4) +  # 添加细胞亚型标签
  theme_minimal() +
  labs(
    x = "T-cell attraction score (median)",
    y = "T-cell suppression score (median)"
  ) +
  theme(legend.position = "none")  # 移除图例






################################################# MSigDB genesets #####################################################

### load all genesets
score_dir <- "/storage/data/KAI/Pan_Inflammation/Homo/Mac2025/analysis/score"
score_list <- dir(score_dir)[-c(4,9,11,14)]  

for (i in 1:length(score_list)) {
  file <- dir(paste(score_dir,score_list[i],'data/',sep = '/'))
  path <- paste(score_dir,score_list[i],'data',file,sep = '/')
  temp <- clusterProfiler::read.gmt(path)
  assign(score_list[i],temp$gene)
}

Inflammasome <- clusterProfiler::read.gmt("/storage/data/KAI/Pan_Inflammation/Homo/Mac2025/analysis/score/Inflammasome/data/REACTOME_INFLAMMASOMES.v2024.1.Hs.gmt")
Inflammasome <- Inflammasome$gene

NLRP3_Inflammasome <- clusterProfiler::read.gmt("/storage/data/KAI/Pan_Inflammation/Homo/Mac2025/analysis/score/Inflammasome/data/REACTOME_THE_NLRP3_INFLAMMASOME.v2024.1.Hs.gmt")
NLRP3_Inflammasome <- NLRP3_Inflammasome$gene


geneSets <- list(
  cGAS = cGAS,
  Angiogenesis = Angiogenesis,
  Complement = Complement,
  Hypoxia = Hypoxia,
  IFNa = IFNa,
  IFNg = IFNg,
  IL6_STAT3 = IL6_STAT3,
  JAK_STAT = JAK_STAT,
  PI3K_AKT = PI3K_AKT,
  Inflammasome = Inflammasome,
  NLRP3_Inflammasome = NLRP3_Inflammasome,
  NFKB = NFKB,
  TGFb = TGFb,
  TNFa = TNFa
)

### All Mac AUCell score
Mac_rankings <- AUCell_buildRankings(GetAssayData(moMac, layer = "counts"))

Mac_AUC <- AUCell_calcAUC(geneSets, Mac_rankings)


### 结果存入seurat对象 Mac中
all(colnames(moMac) == colnames(Mac_AUC)) # TRUE

AUCell_scores <- as.data.frame(t(getAUC(Mac_AUC)))
moMac@meta.data <- cbind(moMac@meta.data, AUCell_scores)





################################################# AUCell MSigDB #####################################################

###可视化

score_sublist <- c("Angiogenesis",
                   "cGAS",
                   "Complement",
                   "Hypoxia", 
                   "IFNa",
                   "IFNg",
                   "IL6_STAT3",
                   "JAK_STAT",
                   "Inflammasome",
                   "NLRP3_Inflammasome",   
                   "NFKB", 
                   "PI3K_AKT",
                   "TGFb",
                   "TNFa") 

p <- VlnPlot(moMac,
        features = score_sublist, 
        pt.size = 0, 
        adjust = 2,
        group.by = "cell_type_level3",
        cols = momac_color,
        raster = F,
        flip = T,
        stack = T
       )
ggplot(p$data, aes(x = ident, y = expression, fill = ident)) +
  geom_violin(adjust = 2) +
  facet_wrap(~feature, scales = "free_y", strip.position = "left") +  # 堆叠特征
  scale_fill_manual(values = momac_color) +  # 使用group.by颜色
  coord_flip() +  # 翻转
  theme_minimal() +
  labs(x = "Cell Type", y = "AUCell") +
  theme(legend.position = "none")




### 炎症相关通路打分可视化

for (i in 1: length(score_sublist)){
  median_scores <- tapply(moMac@meta.data[,score_sublist[i]], moMac$cell_type_level3, median, na.rm = TRUE)
  # 按中位数降序排列
  sorted_levels <- names(sort(median_scores, decreasing = TRUE))
  # 重新设置 cell_type_level3 的因子顺序
  moMac$cell_type_level3 <- factor(moMac$cell_type_level3, levels = sorted_levels)
  p <- VlnPlot(moMac,
              features = score_sublist[i], 
              pt.size = 0, 
              adjust = 2,
              group.by = "cell_type_level3",
              cols = momac_color,
              raster = F
              )
  vln_data <- p$data
  scorename <- colnames(vln_data)[1]
  stats_data <- vln_data %>%
    group_by(ident) %>%
    summarise(mean_value = mean(.data[[scorename]], na.rm = TRUE),
              q25 = quantile(.data[[scorename]], 0.25, na.rm = TRUE),
              q75 = quantile(.data[[scorename]], 0.75, na.rm = TRUE))
  # ANOVA
  anova_result <- aov(as.formula(paste(scorename, "~ ident")), data = vln_data)
  anova_p <- summary(anova_result)[[1]]["Pr(>F)"][1]  # 提取 p 值
  
  anova_p_text <- ifelse(anova_p < 2.2e-16, "< 2.2e-16", signif(anova_p, 3))
  
  temp_vln  <- p + 
        geom_point(data = stats_data, aes(x = ident, y = mean_value), 
                   color = "black", size = 2, shape = 19) +  # 均值点
        geom_segment(data = stats_data, aes(x = ident, xend = ident, 
                                            y = q25, yend = q75), 
                     color = "black", size = 0.7)+  # 25%-75%分位数垂线
        annotate("text", 
                 x = 3,  # 确保 x 轴正确
                 y = max(vln_data[,1], na.rm = TRUE) , 
                 label = paste0("ANOVA p ", anova_p_text[1]), 
                 size = 5, hjust = 0.5)
  assign(paste0(score_sublist[i],"_Vlnplot"),temp_vln)
}



median_scores <- tapply(moMac$Inflammasome, moMac$cell_type_level3, median, na.rm = TRUE)
# 按中位数降序排列
sorted_levels <- names(sort(median_scores, decreasing = TRUE))
# 重新设置 cell_type_level3 的因子顺序
moMac$cell_type_level3 <- factor(moMac$cell_type_level3, levels = sorted_levels)

p<- VlnPlot(moMac,
        features = 'Inflammasome', 
        pt.size = 0, 
        adjust = 2,
        group.by = "cell_type_level3",
        cols = momac_color,
        raster = F

 )

vln_data <- p$data
stats_data <- vln_data %>%
  group_by(ident) %>%
  summarise(mean_value = mean(Inflammasome, na.rm = TRUE),
            q25 = quantile(Inflammasome, 0.25, na.rm = TRUE),
            q75 = quantile(Inflammasome, 0.75, na.rm = TRUE))

anova_result <- aov(Inflammasome ~ ident, data = vln_data)
anova_p <- summary(anova_result)[[1]]["Pr(>F)"][1]  # 提取 p 值

anova_p_text <- ifelse(anova_p < 2.2e-16, "< 2.2e-16", signif(anova_p, 3))

p + 
  geom_point(data = stats_data, aes(x = ident, y = mean_value), 
             color = "black", size = 2, shape = 19) +  # 均值点
  geom_segment(data = stats_data, aes(x = ident, xend = ident, 
                                      y = q25, yend = q75), 
               color = "black", size = 0.7)+  # 25%-75%分位数垂线
  annotate("text", 
           x = median(as.numeric(vln_data$ident)),  # 确保 x 轴正确
           y = max(vln_data[,1], na.rm = TRUE) , 
           label = paste0("ANOVA p = ", anova_p_text[1]), 
           size = 5, hjust = 0.5)
  
# 进行 Tukey HSD 事后检验
tukey_result <- TukeyHSD(anova_result)
tukey_df <- as.data.frame(tukey_result$ident) %>%
  rownames_to_column("comparison") %>%
  mutate(p_signif = case_when(
    `p adj` < 0.001 ~ "***",
    `p adj` < 0.01  ~ "**",
    `p adj` < 0.05  ~ "*",
    TRUE            ~ "ns"
  ))



# IL1B+NLRP3+_Macro
ggarrange(NFKB_Vlnplot, TNFa_Vlnplot, IL6_STAT3_Vlnplot,  
          Inflammasome_Vlnplot,IFNa_Vlnplot,IFNg_Vlnplot,
          JAK_STAT_Vlnplot,Hypoxia_Vlnplot,TGFb_Vlnplot,
          ncol=3,
          nrow=3
          )

# SPP1+CCL2+_Macro
Angiogenesis_Vlnplot + ggtitle("HALLMARK_ANGIOGENESIS")



############################################# Cytokines&chemokines score ###################################################

### load genesets
STTT <- read_xlsx("/storage/data/KAI/Pan_Inflammation/Homo/Mac2025/analysis/score/cytokine/data/STTT2025_TableS2.xlsx")

cytokine <- STTT$`Cytokine scores` %>% na.omit()
Infla_genes <- STTT$`Inflammatory genes` %>% na.omit()

chemokine_prod <- clusterProfiler::read.gmt("/storage/data/KAI/Pan_Inflammation/Homo/Mac2025/analysis/score/chemokine/data/GOBP_CHEMOKINE_PRODUCTION.v2024.1.Hs.gmt")
chemokine_prod <- chemokine_prod$gene
pos_chemokine_prod <- clusterProfiler::read.gmt("/storage/data/KAI/Pan_Inflammation/Homo/Mac2025/analysis/score/chemokine/data/GOBP_POSITIVE_REGULATION_OF_CHEMOKINE_PRODUCTION.v2024.1.Hs.gmt")
pos_chemokine_prod <- pos_chemokine_prod$gene


geneSets <- list(
  Infla_genes = Infla_genes,
  cytokine = cytokine,
  chemokine_prod = chemokine_prod,
  pos_chemokine_prod = pos_chemokine_prod
  
)

### All moMac AUCell score
Mac_rankings <- AUCell_buildRankings(GetAssayData(moMac, layer = "counts"))

Mac_AUC <- AUCell_calcAUC(geneSets, Mac_rankings)

### 结果存入seurat对象 Mac中
all(colnames(moMac) == colnames(Mac_AUC)) # TRUE

AUCell_scores <- as.data.frame(t(getAUC(Mac_AUC)))
moMac@meta.data <- cbind(moMac@meta.data, AUCell_scores)

score_sublist <- c("Infla_genes",
                   "cytokine",
                   "pos_chemokine_prod",
                   "chemokine_prod"
                   ) 

for (i in 1: length(score_sublist)){
  median_scores <- tapply(moMac@meta.data[,score_sublist[i]], moMac$cell_type_level3, median, na.rm = TRUE)
  # 按中位数降序排列
  sorted_levels <- names(sort(median_scores, decreasing = TRUE))
  # 重新设置 cell_type_level3 的因子顺序
  moMac$cell_type_level3 <- factor(moMac$cell_type_level3, levels = sorted_levels)
  p <- VlnPlot(moMac,
               features = score_sublist[i], 
               pt.size = 0, 
               adjust = 2,
               group.by = "cell_type_level3",
               cols = momac_color,
               raster = F
  )
  vln_data <- p$data
  scorename <- colnames(vln_data)[1]
  stats_data <- vln_data %>%
    group_by(ident) %>%
    summarise(mean_value = mean(.data[[scorename]], na.rm = TRUE),
              q25 = quantile(.data[[scorename]], 0.25, na.rm = TRUE),
              q75 = quantile(.data[[scorename]], 0.75, na.rm = TRUE))
  # ANOVA
  anova_result <- aov(as.formula(paste(scorename, "~ ident")), data = vln_data)
  anova_p <- summary(anova_result)[[1]]["Pr(>F)"][1]  # 提取 p 值
  
  anova_p_text <- ifelse(anova_p < 2.2e-16, "< 2.2e-16", signif(anova_p, 3))
  
  temp_vln  <- p + 
    geom_point(data = stats_data, aes(x = ident, y = mean_value), 
               color = "black", size = 2, shape = 19) +  # 均值点
    geom_segment(data = stats_data, aes(x = ident, xend = ident, 
                                        y = q25, yend = q75), 
                 color = "black", size = 0.7)+  # 25%-75%分位数垂线
    annotate("text", 
             x = 3,  # 确保 x 轴正确
             y = max(vln_data[,1], na.rm = TRUE) , 
             label = paste0("ANOVA p ", anova_p_text[1]), 
             size = 5, hjust = 0.5)
  assign(paste0(score_sublist[i],"_Vlnplot"),temp_vln)
}

Infla_genes_Vlnplot + ggtitle("Inflammatory Genes")
cytokine_Vlnplot
pos_chemokine_prod_Vlnplot
chemokine_prod_Vlnplot


############################################# M1/M2 score Addmodulescore #####################################################

# 标准化后再进行打分
IL1BSPP1_Mac <- NormalizeData(IL1BSPP1_Mac, normalization.method = "LogNormalize", scale.factor = 10000)

M1_gene <- c('IL23','TNF','CXCL9','CXCL10','CXCL11','CD86','IL1A','IL1B','IL6','CCL5','IRF5', 'IRF1', 'CD40', 'IDO1', 'KYNU', 'CCR7')

# M2 signature
M2_gene <- c('IL4R', 'CCL4', 'CCL13', 'CCL20', 'CCL17', 'CCL18', 'CCL22', 'CCL24', 'LYVE1', 'VEGFA', 'VEGFB', 'VEGFC', 'VEGFD', 'EGF','CTSA','CTSB' ,'CTSC','CTSD' ,'TGFB1','TGFB2','TGFB3','MMP14','MMP19','MMP9','CLEC7A','WNT7B','FASL','TNFSF12','TNFSF8','CD276','VTCN1','MSR1', 'FN1','IRF4' )

# AddModuleScore
IL1BSPP1_Mac <- AddModuleScore(IL1BSPP1_Mac, features = list(M1_gene), name = "M1_score")
IL1BSPP1_Mac <- AddModuleScore(IL1BSPP1_Mac, features = list(M2_gene), name = "M2_score")



plot_data <- data.frame(
  M1_score = IL1BSPP1_Mac$M1_score1,
  M2_score = IL1BSPP1_Mac$M2_score1,
  celltype = IL1BSPP1_Mac$cell_type_level3
)

p <- ggplot(plot_data, aes(x = M1_score, y = M2_score, color = celltype)) +
  geom_point(size = 1.5, alpha = 0.9) +  # 设置点的大小和透明度
  scale_color_manual(values = Mac_color) +  # 使用自定义颜色映射
  labs(
    x = "M1 Score",  # x 轴标签
    y = "M2 Score",  # y 轴标签
    title = "M1 vs M2 Scores by Cell Type",  # 图表标题
    color = "Macrophage subtypes"  # 图例标题
  ) +
  stat_ellipse(aes(group = celltype), linetype = "dotted", color = "black") +
  theme_minimal() +  # 使用简洁的主题
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # 标题居中加粗
    legend.position = "right",  # 图例放在右侧
    legend.text = element_text(size = 10),  # 图例文字大小
    panel.grid.major = element_blank(),  # 删除主要网格线
    panel.grid.minor = element_blank(),  # 删除次要网格线
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # 添加边框
    axis.line = element_line(color = "black", size = 0.5)  # 添加坐标轴线
  ) +
  guides(color = guide_legend(override.aes = list(size = 4))) +  # 调整图例中点的大小
  geom_hline(yintercept = 0, color = "black", linetype = "dashed", size = 0.5) +  # 添加 y=0 的水平线
  geom_vline(xintercept = 0, color = "black", linetype = "dashed", size = 0.5) +  # 添加 x=0 的垂直线
  coord_cartesian(xlim = c(min(plot_data$M1_score), max(plot_data$M1_score)),  # 设置 x 轴范围
                  ylim = c(min(plot_data$M2_score), max(plot_data$M2_score)))  # 设置 y 轴范围



plot_data$shape <- ifelse(plot_data$celltype == "IL1B+NLRP3+_Macro", 16, 3)  # 17 = triangle, 16 = solid circle

ggplot(plot_data, aes(x = M1_score, y = M2_score)) +
  # 添加密度线（类似等高线）
  geom_density_2d(aes(x = M1_score, y = M2_score), bins = 10, size = 0.3, colour = "black") + 
  # 添加点，按颜色和形状映射
  geom_point(aes(color = celltype, shape = celltype), size = 0.75, alpha = 1) +
  scale_color_manual(values = Mac_color) +
  scale_shape_manual(values = c("IL1B+NLRP3+_Macro" = 16, "SPP1+CCL2+_Macro" = 3)) +
  labs(
    x = "M1 Score",
    y = "M2 Score",
    title = "M1 vs M2 Scores by Macrophage Subtypes",
    color = "Cell Type",
    shape = "Cell Type"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "right",
    legend.text = element_text(size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.line = element_line(color = "black", size = 0.5)
  ) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed", size = 0.5) +
  geom_vline(xintercept = 0, color = "black", linetype = "dashed", size = 0.5) +
  coord_cartesian(xlim = range(plot_data$M1_score),
                  ylim = range(plot_data$M2_score))




