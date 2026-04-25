### Package ------------------------------
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

library(purrr)
library(limma)
library(org.Hs.eg.db)
library(BisqueRNA)



### load scRNA-seq reference ----------------------------
sc.lung <- readRDS("/storage/data/KAI/Pan_Inflammation/Homo/STRNA/adata_ref/Lung/Lung_ref.RDS")

sc.joint <- readRDS("/storage/data/KAI/Pan_Inflammation/Homo/STRNA/adata_ref/Joint/Joint_ref.RDS")


sc.bladder <- readRDS("/storage/data/KAI/Pan_Inflammation/Homo/STRNA/adata_ref/Bladder/Bladder_common.RDS")

### prepare sc.dat --------------------------------------

# Lung
Idents(sc.lung) <- sc.lung@meta.data$cell_type_level3
sc.lung <- NormalizeData(sc.lung, normalization.method = "LogNormalize", scale.factor = 10000)
sc.markers <- FindAllMarkers(sc.lung, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1)
saveRDS(sc.markers,"/storage/data/KAI/Pan_Inflammation/Homo/STRNA/adata_ref/Lung/Lung_allmarkers.RDS")
sc.markers <- readRDS("/storage/data/KAI/Pan_Inflammation/Homo/STRNA/adata_ref/Lung/Lung_allmarkers.RDS")
features <- sc.markers$gene
ref_matrix <- AverageExpression(sc.lung, group.by = "cell_type_level3", assays = "RNA", slot = "counts")$RNA



# Joint
Idents(sc.joint) <- sc.joint@meta.data$cell_type_level3
sc.joint <- NormalizeData(sc.joint, normalization.method = "LogNormalize", scale.factor = 10000)
sc.markers <- FindAllMarkers(sc.joint, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1)

ref_matrix <- AverageExpression(sc.joint, group.by = "cell_type_level3", assays = "RNA", slot = "counts")$RNA



# Bladder 
Idents(sc.bladder) <- sc.bladder@meta.data$cell_type_level3
sc.bladder <- NormalizeData(sc.bladder, normalization.method = "LogNormalize", scale.factor = 10000)
sc.markers <- FindAllMarkers(sc.bladder, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1)
saveRDS(sc.markers,"/storage/data/KAI/Pan_Inflammation/Homo/STRNA/adata_ref/Bladder/Bladder_allmarkers.RDS")
# sc.markers <- readRDS("/storage/data/KAI/Pan_Inflammation/Homo/STRNA/adata_ref/Lung/Lung_allmarkers.RDS")
features <- sc.markers$gene
ref_matrix <- AverageExpression(sc.bladder, group.by = "cell_type_level3", assays = "RNA", slot = "counts")$RNA




### BisqueRNA------------------------------------

## GSE239897 COPD --------------------------
GSE239897_tpm <- read.table("/storage/data/KAI/Pan_Inflammation/Homo_bulk/Lung/COPD/1.GSE239897/GSE239897_RAW/GSE239897_TPM.txt",sep = "\t")

bulk.eset <- Biobase::ExpressionSet(assayData = GSE239897_tpm %>% as.matrix())

scrna_eset <- SeuratToExpressionSet(sc.lung, delimiter = "_", position = 1, version = "v3")
scrna_eset@phenoData$celltype <- sc.lung$cell_type_level3

res <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset, scrna_eset, markers = NULL,use.overlap = FALSE)
ref.based.estimates <- res$bulk.props 
write.table(ref.based.estimates,"/storage/data/KAI/Pan_Inflammation/Homo_bulk/Lung/COPD/1.GSE239897/GSE239897_bisqueRNA_level3.txt", quote = F, sep = '\t')

GSE239897_Sample_Info <- read_xlsx("/storage/data/KAI/Pan_Inflammation/Homo_bulk/Lung/COPD/1.GSE239897/GSE239897_Sample_Info.xlsx")

GSE239897_SI_sub <- GSE239897_Sample_Info[GSE239897_Sample_Info$Sample_geo_accession %in% colnames(ref.based.estimates),]

GSE239897_decon <- data.frame(
  Group = GSE239897_SI_sub$Group[match(GSE239897_SI_sub$Sample_geo_accession,colnames(ref.based.estimates))] %>% sub(pattern = 'Control \\(donor\\)',replacement = 'Normal'),
  surface_density = GSE239897_SI_sub$`surface density`[match(GSE239897_SI_sub$Sample_geo_accession,colnames(ref.based.estimates))],
  IL1B_fraction = ref.based.estimates["IL1B+NLRP3+_Macro",]
)

GSE239897_decon$surface_density <- as.numeric(GSE239897_decon$surface_density)
GSE239897_decon$IL1B_fraction <- as.numeric(GSE239897_decon$IL1B_fraction)


ggscatter(GSE239897_decon, x = "IL1B_fraction", y = "surface_density", 
          color = "Group",fill = "lightgray",
          add = "reg.line", conf.int = TRUE, 
          add.params = list( color = "black",fill = "lightgray"),
          cor.coef = T,
          cor.method = "pearson"
) +
  ggtitle("GSE239897") +  # 设置标题
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),  # 添加黑色边框
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),     # 设置标题样式
    axis.text = element_text(color = "black")                             # 坐标轴文字颜色
  ) 

copd_df <- GSE239897_decon[GSE239897_decon$Group == "COPD", ]

median_sd <- median(copd_df$surface_density, na.rm = TRUE)
copd_df$SD_group <- ifelse(copd_df$surface_density >= median_sd, "High", "Low")
wilcox_res <- wilcox.test(IL1B_fraction ~ SD_group, data = copd_df)

ggboxplot(copd_df, x = "SD_group", y = "IL1B_fraction", fill = "SD_group",
          palette = c("skyblue", "tomato")) +
  stat_compare_means(method = "wilcox.test") + 
  labs(title = "IL1B+ Macro in COPD by surface_density",
       x = "Surface Density Group",
       y = "IL1B_fraction") +
  theme_minimal()



## GSE184316 HP --------------------------
GSE184316_tpm <- read_tsv("/storage/data/KAI/Pan_Inflammation/Homo_bulk/Lung/Hypersensitivity_pneumonitis/1.GSE184316/GSE184316_RAW/GSE184316_norm_counts_TPM_GRCh38.p13_NCBI.tsv")

gene_ids <- GSE184316_tpm$GeneID %>% as.character()
gene_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = gene_ids,
  columns = c("SYMBOL"),
  keytype = "ENTREZID"
)

gene_map <- gene_map %>% 
  filter(SYMBOL != "") %>% 
  distinct(ENTREZID, .keep_all = TRUE)

GSE184316_tpm <- GSE184316_tpm[GSE184316_tpm$GeneID %in% gene_map$ENTREZID, ]
GSE184316_tpm <- as.matrix(GSE184316_tpm)
rownames(GSE184316_tpm) <- gene_map$SYMBOL[match(GSE184316_tpm[,1], gene_map$ENTREZID)]
# GSE184316_tpm <- GSE184316_tpm[,-1]
GSE184316_tpm <- GSE184316_tpm[!duplicated(rownames(GSE184316_tpm)), ]

write.table(GSE184316_tpm,"/storage/data/KAI/Pan_Inflammation/Homo_bulk/Lung/Hypersensitivity_pneumonitis/1.GSE184316/GSE184316_RAW/GSE184316_TPM.txt",quote = F,sep = "\t")


bulk.eset <- Biobase::ExpressionSet(assayData = GSE184316_tpm %>% as.matrix())

res <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset, scrna_eset, markers = NULL,use.overlap = FALSE)
saveRDS(res,"/storage/data/KAI/Pan_Inflammation/Homo_bulk/Lung/Hypersensitivity_pneumonitis/1.GSE184316/bisqueRNA_level3.RDS")

ref.based.estimates <- res$bulk.props 
write.table(ref.based.estimates,"/storage/data/KAI/Pan_Inflammation/Homo_bulk/Lung/Hypersensitivity_pneumonitis/1.GSE184316/GSE184316_bisqueRNA_level3.txt", quote = F, sep = '\t')


GSE184316_Sample_Info <- read_xlsx("/storage/data/KAI/Pan_Inflammation/Homo_bulk/Lung/Hypersensitivity_pneumonitis/1.GSE184316/GSE184316_Sample_Info.xlsx")

GSE184316_SI_sub <- GSE184316_Sample_Info[GSE184316_Sample_Info$Sample_geo_accession %in% colnames(ref.based.estimates),]

GSE184316_decon <- data.frame(
  Group = GSE184316_SI_sub$Group[match(GSE184316_SI_sub$Sample_geo_accession,colnames(ref.based.estimates))] %>% sub(pattern = 'DONOR',replacement = 'Normal'),
  surface_density = GSE184316_SI_sub$`surface density`[match(GSE184316_SI_sub$Sample_geo_accession,colnames(ref.based.estimates))],
  IL1B_fraction = ref.based.estimates["IL1B+NLRP3+_Macro",]
)

GSE184316_decon$surface_density <- as.numeric(GSE184316_decon$surface_density)
GSE184316_decon$IL1B_fraction <- as.numeric(GSE184316_decon$IL1B_fraction)


ggscatter(GSE184316_decon, x = "IL1B_fraction", y = "surface_density", 
          color = "Group",fill = "lightgray",
          add = "reg.line", conf.int = TRUE, 
          add.params = list( color = "black",fill = "lightgray"),
          cor.coef = T,
          cor.method = "pearson"
) +
  ggtitle("GSE239897") +  # 设置标题
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),  # 添加黑色边框
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),     # 设置标题样式
    axis.text = element_text(color = "black")                             # 坐标轴文字颜色
  ) 


hp_df <- GSE184316_decon[GSE184316_decon$Group == "HP", ]

ggscatter(hp_df, x = "IL1B_fraction", y = "surface_density", 
          color = "black",fill = NULL,
          size = 2.5, shape = 16, stroke = 2,
          add = "reg.line", conf.int = TRUE, 
          add.params = list( color = "#FA8072",fill = "#ecced0"),
          cor.coeff.args = list(
            label.x = 0.05,        # 控制相关系数的位置 (x轴)
            label.y = 0.0155,       # 控制相关系数的位置 (y轴)
            size = 5,              # 字体大小
            color = "#88878d",        # 字体颜色
            fontface = "bold"      # 字体样式
          ),
          cor.coef = T,
          cor.method = "pearson"
) +
  ggtitle("GSE184316") +  # 设置标题
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # 添加黑色边框
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),     # 设置标题样式
    axis.text = element_text(color = "black")                             # 坐标轴文字颜色
  ) +
  labs(y = "Surface Density",
       x = "IL1B+NLRP3+_Macro fraction"
       )


median_sd <- median(copd_df$surface_density, na.rm = TRUE)
copd_df$SD_group <- ifelse(copd_df$surface_density >= median_sd, "High", "Low")
wilcox_res <- wilcox.test(IL1B_fraction ~ SD_group, data = copd_df)

ggboxplot(copd_df, x = "SD_group", y = "IL1B_fraction", fill = "SD_group",
          palette = c("skyblue", "tomato")) +
  stat_compare_means(method = "wilcox.test") + 
  labs(title = "IL1B+ Macro in COPD by surface_density",
       x = "Surface Density Group",
       y = "IL1B_fraction") +
  theme_minimal()


## GSE225731  RA ----------------------------
sc.joint_decon <- sc.joint_decon
scrna_eset <- SeuratToExpressionSet(sc.joint, delimiter = "_", position = 1, version = "v3")
scrna_eset@phenoData$celltype <- sc.joint$cell_type_level3



GSE225731_tpm <- read_tsv("/storage/data/KAI/Pan_Inflammation/Homo_bulk/Joint/Rheumatoid_arthritis/1.GSE225731/GSE225731_RAW/GSE225731_norm_counts_TPM_GRCh38.p13_NCBI.tsv")

gene_ids <- GSE225731_tpm$GeneID %>% as.character()
gene_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = gene_ids,
  columns = c("SYMBOL"),
  keytype = "ENTREZID"
)

gene_map <- gene_map %>% 
  filter(SYMBOL != "") %>% 
  distinct(ENTREZID, .keep_all = TRUE)

GSE225731_tpm <- GSE225731_tpm[GSE225731_tpm$GeneID %in% gene_map$ENTREZID, ]
GSE225731_tpm <- as.matrix(GSE225731_tpm)
rownames(GSE225731_tpm) <- gene_map$SYMBOL[match(GSE225731_tpm[,1], gene_map$ENTREZID)]
# GSE225731_tpm <- GSE225731_tpm[,-1]
GSE225731_tpm <- GSE225731_tpm[!duplicated(rownames(GSE225731_tpm)), ]

write.table(GSE225731_tpm,"/storage/data/KAI/Pan_Inflammation/Homo_bulk/Joint/Rheumatoid_arthritis/1.GSE225731/GSE225731_RAW/GSE225731_TPM.txt",quote = F,sep = "\t")

bulk.eset <- Biobase::ExpressionSet(assayData = GSE225731_tpm %>% as.matrix())

res <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset, scrna_eset, markers = NULL,use.overlap = FALSE)
saveRDS(res,"/storage/data/KAI/Pan_Inflammation/Homo_bulk/Joint/Rheumatoid_arthritis/1.GSE225731/bisqueRNA/bisqueRNA_level3.RDS")

ref.based.estimates <- res$bulk.props 
write.table(ref.based.estimates,"/storage/data/KAI/Pan_Inflammation/Homo_bulk/Joint/Rheumatoid_arthritis/1.GSE225731/bisqueRNA/GSE225731_bisqueRNA_level3.txt", quote = F, sep = '\t')






GSE225731_Sample_Info <- read_xlsx("/storage/data/KAI/Pan_Inflammation/Homo_bulk/Joint/Rheumatoid_arthritis/1.GSE225731/GSE225731_Sample_Info.xlsx")


GSE225731_SI_sub <- GSE225731_Sample_Info[GSE225731_Sample_Info$Sample_geo_accession %in% colnames(ref.based.estimates),]

GSE225731_decon <- data.frame(
  Group = GSE225731_SI_sub$`disease state`[match(GSE225731_SI_sub$Sample_geo_accession,colnames(ref.based.estimates))] %>% sub(pattern = 'rheumatoid arthritis',replacement = 'RA'),
  haq_bl = GSE225731_SI_sub$haq_bl[match(GSE225731_SI_sub$Sample_geo_accession,colnames(ref.based.estimates))],
  das28 = GSE225731_SI_sub$das28crp_bl[match(GSE225731_SI_sub$Sample_geo_accession,colnames(ref.based.estimates))],
  IL1B_fraction = ref.based.estimates["IL1B+NLRP3+_Macro",]
)

GSE225731_decon$haq_bl <- as.numeric(GSE225731_decon$haq_bl)
GSE225731_decon$IL1B_fraction <- as.numeric(GSE225731_decon$IL1B_fraction)


ggscatter(GSE225731_decon, x = "IL1B_fraction", y = "surface_density", 
          color = "Group",fill = "lightgray",
          add = "reg.line", conf.int = TRUE, 
          add.params = list( color = "black",fill = "lightgray"),
          cor.coef = T,
          cor.method = "pearson"
) +
  ggtitle("GSE239897") +  # 设置标题
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),  # 添加黑色边框
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),     # 设置标题样式
    axis.text = element_text(color = "black")                             # 坐标轴文字颜色
  ) 


ra_df <- GSE225731_decon[GSE225731_decon$Group == "RA", ]

ggscatter(ra_df, x = "IL1B_fraction", y = "haq_bl", 
          color = "black",fill = NULL,
          size = 2.5, shape = 16, stroke = 2,
          add = "reg.line", conf.int = TRUE, 
          add.params = list( color = "#FA8072",fill = "#ecced0"),
          cor.coeff.args = list(
            label.x = 0.05,        # 控制相关系数的位置 (x轴)
            label.y = 0.0155,       # 控制相关系数的位置 (y轴)
            size = 5,              # 字体大小
            color = "#88878d",        # 字体颜色
            fontface = "bold"      # 字体样式
          ),
          cor.coef = T,
          cor.method = "pearson"
) +
  ggtitle("GSE225731") +  # 设置标题
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # 添加黑色边框
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),     # 设置标题样式
    axis.text = element_text(color = "black")                             # 坐标轴文字颜色
  ) +
  labs(y = "Surface Density",
       x = "IL1B+NLRP3+_Macro fraction"
  )




## GSE238208 IC ----------------------------




# tpm preprocess

GSE238208_tpm <- read_tsv("/storage/data/KAI/Pan_Inflammation/Homo_bulk/Bladder/1.GSE238208/GSE238208_RAW/GSE238208_norm_counts_TPM_GRCh38.p13_NCBI.tsv")

gene_ids <- GSE238208_tpm$GeneID %>% as.character()
gene_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = gene_ids,
  columns = c("SYMBOL"),
  keytype = "ENTREZID"
)

gene_map <- gene_map %>% 
  filter(SYMBOL != "") %>% 
  distinct(ENTREZID, .keep_all = TRUE)

GSE238208_tpm <- GSE238208_tpm[GSE238208_tpm$GeneID %in% gene_map$ENTREZID, ]
GSE238208_tpm <- as.matrix(GSE238208_tpm)
rownames(GSE238208_tpm) <- gene_map$SYMBOL[match(GSE238208_tpm[,1], gene_map$ENTREZID)]
# GSE238208_tpm <- GSE238208_tpm[,-1]
GSE238208_tpm <- GSE238208_tpm[!duplicated(rownames(GSE238208_tpm)), ]

# write.table(GSE238208_tpm,"/storage/data/KAI/Pan_Inflammation/Homo_bulk/Bladder/1.GSE238208/GSE238208_RAW//GSE238208_TPM.txt",quote = F,sep = "\t")

bulk.eset <- Biobase::ExpressionSet(assayData = GSE238208_tpm %>% as.matrix())

# scrna ref data
scrna_eset <- SeuratToExpressionSet(sc.bladder, delimiter = "_", position = 1, version = "v3")
scrna_eset@phenoData$celltype <- sc.bladder$cell_type_level3
res <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset, scrna_eset, markers = NULL,use.overlap = FALSE)
saveRDS(res,"/storage/data/KAI/Pan_Inflammation/Homo_bulk/Bladder/1.GSE238208/bisqueRNA/bisqueRNA_level3.RDS")

ref.based.estimates <- res$bulk.props 
write.table(ref.based.estimates,"/storage/data/KAI/Pan_Inflammation/Homo_bulk/Bladder/1.GSE238208/bisqueRNA/GSE238208_bisqueRNA_level3.txt", quote = F, sep = '\t')

GSE238208_Sample_Info <- read_xlsx("/storage/data/KAI/Pan_Inflammation/Homo_bulk/Bladder/1.GSE238208/GSE238208_RAW/GSE238208_Sample_Info.xlsx")


GSE238208_SI_sub <- GSE238208_Sample_Info[GSE238208_Sample_Info$Sample_geo_accession %in% colnames(ref.based.estimates),]

GSE238208_decon <- data.frame(
  Group = GSE238208_SI_sub$Grade[match(GSE238208_SI_sub$Sample_geo_accession,colnames(ref.based.estimates))] ,
  IL1B_fraction = ref.based.estimates["IL1B+NLRP3+_Macro",]
)

df <- GSE238208_decon
df$Group <- factor(
  df$Group,
  levels = c(
    "Non_Hunner_Lesion",
    "reddened mucosal flat lesion",
    "Hunner_Lesion"
  )
)
ggplot(df, aes(x = Group, y = IL1B_fraction, fill = Group)) +
  geom_boxplot(
    width = 0.6,
    outlier.shape = NA,
    alpha = 0.7
  ) +
  geom_jitter(
    width = 0.15,
    size = 1.8,
    alpha = 0.8
  ) +
  labs(
    x = NULL,
    y = "IL1B fraction"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

df2 <- subset(
  GSE238208_decon,
  Group %in% c("Non_Hunner_Lesion", "Hunner_Lesion")
)

df2$Group <- factor(
  df2$Group,
  levels = c("Non_Hunner_Lesion", "Hunner_Lesion")
)

group_cols <- c(
  "Non_Hunner_Lesion" = "#7F8C8D",  # 灰蓝
  "Hunner_Lesion"     = "#C0392B"   # 深红
)


ggplot(df2, aes(x = Group, y = IL1B_fraction, fill = Group)) +
  geom_boxplot(
    width = 0.4,
    outlier.shape = NA,
    alpha = 0.5
  ) +
  geom_jitter(
    width = 0.2,
    size = 1.2,
    alpha = 0.75
  ) +
  stat_compare_means(
    method = "wilcox.test",
    alternative = "less",   
    label = "p.format"
  ) +
  scale_fill_manual(values = group_cols) +
  scale_color_manual(values = group_cols) +
  labs(
    x = NULL,
    y = "IL1B fraction"
  ) +
  theme_classic(base_size = 14) +
  ggtitle("GSE238208") +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 13)
  )



wilcox.test(
  IL1B_fraction ~ Group,
  data = subset(df, Group %in% c("Hunner_Lesion", "Non_Hunner_Lesion")),
  alternative = "less"
)
