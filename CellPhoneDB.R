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
library(data.table)


library(CellChat)
library(scde)
library(DEsingle)
library(MAST)
library(scater)

library(iTALK)
library(circlize)
library(igraph)

library(cpplot)

library(ktplots)
library(SingleCellExperiment)


### load data -------------------------------
Mac_10 <- readRDS("/storage/data/KAI/Pan_Inflammation/Homo/Mac2025/analysis/monocle2/data/Mac_10.RDS") # 同monocle2抽样输入数据
Idents(Mac_10) <- Mac_10@meta.data$cell_type_level3

T_10 <- readRDS("/storage/data/KAI/Pan_Inflammation/Homo/Mac2025/analysis/CellChat/data/T_10.RDS")
Endo_10 <- readRDS("/storage/data/KAI/Pan_Inflammation/Homo/Mac2025/analysis/CellChat/data/Endo_10.RDS")
B_10 <- readRDS("/storage/data/KAI/Pan_Inflammation/Homo/Mac2025/analysis/CellChat/data/B_10.RDS")
Fibro_10 <- readRDS("/storage/data/KAI/Pan_Inflammation/Homo/Mac2025/analysis/CellChat/data/Fibro_10.RDS")


### CellPhoneDB input -----------------------
merge_seurat <- merge(Mac_10, y = list(Endo_10, T_10, B_10, Fibro_10), 
                      add.cell.ids = c("Macro", "Endo", "T", "B", "Fibro")) 
write.table(as.matrix(merge_seurat@assays$RNA$counts),file = "/storage/data/KAI/Pan_Inflammation/Homo/Mac2025/analysis/CellPhoneDB/data/input/panCell.txt",sep = '\t', quote = F)



meta_data <- cbind(rownames(merge_seurat@meta.data),merge_seurat@meta.data[,'cell_type_level3', drop=F])
meta_data <- as.matrix(meta_data)
colnames(meta_data) <- c('Cell','cell_type_level3')
is.na(meta_data) %>% sum() # 没有NA值

write.table(meta_data, "/storage/data/KAI/Pan_Inflammation/Homo/Mac2025/analysis/CellPhoneDB/data/input/cellphonedb_meta.txt",sep = '\t', quote = F)


merge_infla <- subset(merge_seurat, Inflammation == "Inflamed")
merge_infla <- NormalizeData(merge_infla, normalization.method = "LogNormalize", scale.factor = 10000)
write.table(as.matrix(merge_infla@assays$RNA$data),file = "/storage/data/KAI/Pan_Inflammation/Homo/Mac2025/analysis/CellPhoneDB/data/input/panCell_Inflamed.txt",sep = '\t', quote = F)
# fwrite(as.matrix(merge_infla@assays$RNA$data), file = "/storage/data/KAI/Pan_Inflammation/Homo/Mac2025/analysis/CellPhoneDB/data/input/panCell_Inflamed.txt", sep = "\t", quote = FALSE)



meta_infla <- cbind(rownames(merge_infla@meta.data),merge_infla@meta.data[,'cell_type_level3', drop=F])
meta_infla <- as.matrix(meta_infla)
colnames(meta_infla) <- c('Cell','cell_type_level3')
is.na(meta_infla) %>% sum()
dim(meta_infla)
write.table(meta_infla, "/storage/data/KAI/Pan_Inflammation/Homo/Mac2025/analysis/CellPhoneDB/data/input/cellphonedb_meta_infla.txt",sep = '\t', quote = F)
# fwrite(meta_infla, file = "/storage/data/KAI/Pan_Inflammation/Homo/Mac2025/analysis/CellPhoneDB/data/input/cellphonedb_meta_infla.txt", sep = "\t", quote = FALSE)


merge_norm <- subset(merge_seurat, Inflammation == "Normal")
merge_norm <- NormalizeData(merge_norm, normalization.method = "LogNormalize", scale.factor = 10000)
# write.table(as.matrix(merge_norm@assays$RNA$data),file = "/storage/data/KAI/Pan_Inflammation/Homo/Mac2025/analysis/CellPhoneDB/data/input/panCell_Normal.txt",sep = '\t', quote = F)
fwrite(as.matrix(merge_norm@assays$RNA$data), file = "/storage/data/KAI/Pan_Inflammation/Homo/Mac2025/analysis/CellPhoneDB/data/input/panCell_Normal.txt", sep = "\t", quote = FALSE)


meta_norm <- cbind(rownames(merge_norm@meta.data),merge_norm@meta.data[,'cell_type_level3', drop=F])
meta_norm <- as.matrix(meta_norm)
colnames(meta_norm) <- c('Cell','cell_type_level3')
is.na(meta_norm) %>% sum()
dim(meta_norm)
# write.table(meta_norm, "/storage/data/KAI/Pan_Inflammation/Homo/Mac2025/analysis/CellPhoneDB/data/input/cellphonedb_meta_norm.txt",sep = '\t', quote = F)
fwrite(meta_norm, file = "/storage/data/KAI/Pan_Inflammation/Homo/Mac2025/analysis/CellPhoneDB/data/input/cellphonedb_meta_norm.txt", sep = "\t", quote = FALSE)





### CellPhoneDB visual (igraph based)------------------
setwd("/storage/data/KAI/Pan_Inflammation/Homo/Mac2025/analysis/CellPhoneDB/result/infla_output/")
all_files <- list.files('.', full.names = TRUE)

# 解析cellphonedb输出文件：means文件
network_dat <- all_files  %>% 
  str_subset('_analysis_means_') %>% 
  read_tsv() %>% 
  dplyr::select(interacting_pair, where(is.numeric)) %>%
  pivot_longer(
    -interacting_pair,
    names_to = 'cell_pair',
    values_to = 'p'
  ) %>% 
  dplyr::filter(p < 0.05) %>% 
  group_by(cell_pair) %>% 
  summarise(counts = n()) %>% 
  mutate(SOURCE = str_split(cell_pair, "\\|") %>% map_chr(1),
         TARGET = str_split(cell_pair, "\\|") %>% map_chr(2)) %>% 
  dplyr::select(SOURCE, TARGET, counts) %>% 
  dplyr::filter(counts != 0) 

# 整体网络
all_nodes <- unique(c(network_dat$SOURCE, network_dat$TARGET))

net <- graph_from_data_frame(network_dat)

E(net)$width <- E(net)$counts %>% {(. - min(.))/diff(range(.)) + 1}

# E(net)$color <- attr(E(net), 'vnames' ) %>% str_remove("\\|.*$")  %>% { color_match_table[.]}
# V(net)$color <- color_match_table[names(V(net))]

karate_groups <- cluster_optimal(net)
coords        <- layout_in_circle(net, order = order(membership(karate_groups)))  # 设置网络布局

plot(
  net,
  edge.curved=0.3,
  edge.arrow.size=1
  #layout = coords,

)


### CellPhoneDB visual (cpplot based)------------------

decon = read.delim(dir(pattern="deconvoluted_\\d{2}")[1], check.names = FALSE)


ccc_number_heatmap1(pfile = dir(pattern="analysis_pvalues")) #ggplot对象
ccc_number_heatmap2(pfile = dir(pattern="analysis_pvalues")) #ggplot对象
ccc_number_line(pfile = "test/pvalues.txt",vertex.size = 20) #不是ggplot对象，不能用ggsave保存
# Error in log10(pvalues.df1$pvalue) : 数学函数中用了非数值参数 

ccc_bubble(
  pfile=dir(pattern="analysis_pvalues"),
  mfile=dir(pattern="analysis_means"),
  cell.pair=c("IL1B+NLRP3+_Macro|Venous_EC",
              "IL1B+NLRP3+_Macro|Arterial_EC",
              "IL1B+NLRP3+_Macro|Lymphatic_EC",
              "IL1B+NLRP3+_Macro|Capillary_EC"
              )
  # 下面这些是默认参数，可以不变
  # neg_log10_th = -log10(0.05),
  # means_exp_log2_th = 1,
  # notused.cell = NULL,
  # used.cell = NULL,
  # neg_log10_th2 = 3,
  # means_exp_log2_th2 = c(-4, 6),
  # cell.pair = NULL,
  # gene.pair = NULL,
  # color_palette = c("#313695", "#4575B4", "#ABD9E9", "#FFFFB3", "#FDAE61", "#F46D43","#D73027", "#A50026"),
  # text_size = 12
)
# Error in log10(pvalues.df1$pvalue) : 数学函数中用了非数值参数 







#### cellphonedb level1 analysis input -----------------------
panCell_Inflamed10 <- readRDS("/storage/data/KAI/Pan_Inflammation/Homo/Mac2025/analysis/CellPhoneDB/level1_analysis/data/panCell_Inflamed10.RDS")




### CellPhoneDB visual (ktplots based)------------------
# load result data
setwd("/storage/data/KAI/Pan_Inflammation/Homo/Mac2025/analysis/CellPhoneDB/level1_analysis/result/infla_output/")
all_files <- list.files('.', full.names = TRUE)

pvals <- read.delim(dir(pattern="analysis_pvalues"), check.names = FALSE)
means <- read.delim(dir(pattern="analysis_means"), check.names = FALSE)
decon <- read.delim(dir(pattern = "deconvoluted")[1], check.names = FALSE)

# singlecellexperiment object
infla_10 <- panCell_Inflamed10 %>% as.SingleCellExperiment()

saveRDS(infla_10,"/storage/data/KAI/Pan_Inflammation/Homo/Mac2025/analysis/CellPhoneDB/data/infla_10_sce.RDS")


# Heatmap
plot_cpdb_heatmap(pvals = pvals, cellheight = 15, cellwidth = 15)

# Mac_T/NK
plot_cpdb_heatmap(pvals = pvals, cell_types = c("IL1B+NLRP3+_Macro","SPP1+CCL2+_Macro", "Mast cell", "T/NK cell",), cellheight = 20, cellwidth = 20)





# Mac_Endo
plot_cpdb(
  scdata = infla_10,
  cell_type1 = "IL1B+NLRP3+_Macro|SPP1+CCL2+_Macro",
  cell_type2 = "Endothelial cell", # "." this means all cell-types
  celltype_key = "cell_type_CC",
  means = means,
  pvals = pvals,
  genes = c('VEGFA','VEGFB','ICAM1','VCAM1','ANGPT1','PDGFA','PDGFB'),
  title = "interacting interactions",
  keep_id_cp_interaction = F
)

# Mac_Fibro
plot_cpdb(
  scdata = infla_10,
  cell_type1 = "IL1B+NLRP3+_Macro|SPP1+CCL2+_Macro",
  cell_type2 = "Fibroblast", # "." this means all cell-types
  celltype_key = "cell_type_CC",
  means = means,
  pvals = pvals,
  genes = c('TGFB1','TGFB2','TGFB3','TNF','LGALS3','AREG'),
  title = "interacting interactions",
  keep_id_cp_interaction = F
)


# Mac_T/NK
plot_cpdb(
  scdata = infla_10,
  cell_type1 = "IL1B+NLRP3+_Macro|SPP1+CCL2+_Macro",
  cell_type2 = "T/NK cell", # "." this means all cell-types
  celltype_key = "cell_type_CC",
  means = means,
  pvals = pvals,
  # gene_family = "chemokines",
  title = "interacting interactions",
  keep_id_cp_interaction = F
)

# Mac_B/Plasma
plot_cpdb(
  scdata = infla_10,
  cell_type1 = "IL1B+NLRP3+_Macro|SPP1+CCL2+_Macro",
  cell_type2 = "B cell", # "." this means all cell-types
  celltype_key = "cell_type_CC",
  means = means,
  pvals = pvals,
  gene_family = "chemokines",
  title = "interacting interactions",
  keep_id_cp_interaction = F
)


plot_cpdb3(
  scdata = infla_10,
  cell_type1 = "SPP1+CCL2+_Macro",
  cell_type2 = "Plasmablast",
  celltype_key = "cell_type_level3", # column name where the cell ids are located in the metadata
  means = means,
  pvals = pvals,
  deconvoluted = decon, # new options from here on specific to plot_cpdb3
  keep_significant_only = TRUE
)




#### Cross-organs cell-cell interactions ------------------------------------------
### sepcific packages
library(tidyverse)
library(fs)
library(pheatmap)
library(ktplots)
library(ComplexHeatmap)


library(dplyr)
library(tidyr)
library(readr)   # 或 vroom
library(stringr)
library(ggplot2)
library(forcats)


### load cellphonedb result

root <- "/storage/data/KAI/Pan_Inflammation/Homo/Mac2025/analysis/CellPhoneDB/result_new"
org_dirs <- fs::dir_ls(root, type = "directory")

cpdb_count_one_organ <- function(org_dir,
                                 prefer = c("significant_means","pvalues"),
                                 p_cut = 0.05,   # 回退到 pvalues 时使用
                                 collapse_symmetric = TRUE,  # A|B 与 B|A 合并
                                 verbose = TRUE) {
  prefer <- match.arg(prefer)
  organ  <- basename(org_dir)
  
  files <- list.files(org_dir, full.names = TRUE)
  # 两类文件（允许 .txt 或 .txt.gz）
  f_sig <- files[grepl("statistical_analysis_significant_means_.*\\.txt(\\.gz)?$", basename(files), ignore.case = TRUE)]
  f_p   <- files[grepl("statistical_analysis_pvalues_.*\\.txt(\\.gz)?$",         basename(files), ignore.case = TRUE)]
  
  # 选择使用哪个
  chosen <- character(0); mode <- NA_character_
  if (prefer == "significant_means" && length(f_sig) > 0) {
    chosen <- f_sig[1]; mode <- "significant_means"
  } else if (length(f_p) > 0) {
    chosen <- f_p[1];   mode <- "pvalues"
  } else if (length(f_sig) > 0) {  # 兜底
    chosen <- f_sig[1]; mode <- "significant_means"
  } else {
    if (verbose) message("No usable CPDB file in: ", org_dir)
    return(tibble())
  }
  
  # 读取
  df <- readr::read_tsv(chosen, show_col_types = FALSE)
  if (!nrow(df)) return(tibble())
  
  # 找出形如 "A|B" 的列（每列代表一个细胞对）
  pair_cols <- names(df)[grepl("\\|", names(df))]
  if (!length(pair_cols)) return(tibble())
  
  # 拉成长表 + 选择显著
  if (mode == "significant_means") {
    long <- df %>%
      tidyr::pivot_longer(dplyr::all_of(pair_cols), names_to = "pair", values_to = "val") %>%
      dplyr::filter(!is.na(val)) %>%
      dplyr::select(pair)
  } else {
    long <- df %>%
      tidyr::pivot_longer(dplyr::all_of(pair_cols), names_to = "pair", values_to = "p") %>%
      dplyr::filter(!is.na(p), p < p_cut) %>%
      dplyr::select(pair)
  }
  if (!nrow(long)) return(tibble())
  
  # 规范“无向对”：A|B 与 B|A 视作同一列；去多余空格
  long <- long %>%
    tidyr::separate(pair, c("A","B"), sep = "\\|", remove = FALSE) %>%
    dplyr::mutate(A = stringr::str_squish(A),
                  B = stringr::str_squish(B)) %>%
    dplyr::mutate(pair_canon = if (collapse_symmetric) {
      purrr::pmap_chr(list(A,B), ~ paste(sort(c(..1, ..2)), collapse = " \u2194 ")) # ↔
    } else {
      paste(A, B, sep = " \u2194 ")
    })
  
  # 计数：该器官里，每个“细胞对”的显著条目数
  out <- long %>%
    dplyr::count(organ = organ, pair = pair_canon, name = "n")
  
  if (verbose) message(sprintf("[%-20s] %s: %d significant entries, %d unique pairs",
                               organ, mode, sum(out$n), nrow(out)))
  out
}
res_all <- purrr::map_dfr(org_dirs, cpdb_count_one_organ)

CCI_all <- res_all %>%
  tidyr::pivot_wider(names_from = pair, values_from = n, values_fill = 0) %>%
  dplyr::arrange(organ) %>%
  tibble::column_to_rownames("organ") %>%
  as.matrix()

rownames(CCI_all) <- gsub("_inflamed$","", rownames(CCI_all))

### 细胞互作类型过多 拆分为巨噬细胞亚型与免疫细胞/基质细胞互作两个热图

## 基质互作
incl_macros <- c(
  "IL1B+NLRP3+_Macro",
  "SPP1+CCL2+_Macro",
  "CD1C+_DC-like_Macro",
  "CCL13+_Complement-associated_Macro",
  "MMP3+CXCL8+_Macro"
)

partners_stromal <- c("Epithelial cell", "Mural cell", "Fibroblast", "Endothelial cell", "Neural cell")

organ_order <- rownames(CCI_all)
arrow <- "\u2194"



# 拆分细胞互作对名称
pairs_df <- tibble(col = colnames(CCI_all)) %>%
  tidyr::separate(col, c("A","B"), sep = "\\s*↔\\s*", remove = FALSE) %>%
  mutate(A = stringr::str_squish(A), B = stringr::str_squish(B))

CCI_long <- as.data.frame(CCI_all) %>%
  tibble::rownames_to_column("organ") %>%
  tidyr::pivot_longer(-organ, names_to = "col", values_to = "n") %>%
  left_join(pairs_df, by = "col")

# 去除不含指定Macro的互作对
CCI_long <- CCI_long %>%
  mutate(
    Macro   = dplyr::case_when(A %in% incl_macros ~ A,
                               B %in% incl_macros ~ B,
                               TRUE ~ NA_character_),
    Partner = dplyr::case_when(!is.na(Macro) & Macro == A ~ B,
                               !is.na(Macro) & Macro == B ~ A,
                               TRUE ~ NA_character_)
  ) %>%
  filter(!is.na(Macro))   # 丢弃两侧都不是 5 个 Macro 的列


# 将Skin中的Keratinocyte转换成Epithelial cell进行处理
CCI_long <- CCI_long %>%
  mutate(
    Partner2 = if_else(grepl("^skin$", organ, ignore.case = TRUE) & Partner == "Keratinocyte",
                       "Epithelial cell", Partner)
  )

CCI_grp <- CCI_long %>%
  group_by(organ, Macro, Partner = Partner2) %>%
  summarise(n = sum(n), .groups = "drop")

# 提取与Macro互作为Stromal cell的互作类型
stromal_df <- CCI_grp %>%
  filter(Partner %in% partners_stromal)


stromal_cols_order <- tidyr::expand_grid(Macro = incl_macros, Partner = partners_stromal) %>%
  mutate(pair = paste(Macro, '↔', Partner)) %>%
  pull(pair)

mat_stromal <- stromal_df %>%
  dplyr::mutate(pair = paste(Macro, arrow, Partner)) %>%
  dplyr::select(organ, pair, n) %>%                                            
  tidyr::pivot_wider(names_from = pair, values_from = n, values_fill = 0) %>%
  dplyr::right_join(tibble::tibble(organ = organ_order), by = "organ") %>%     
  dplyr::arrange(factor(organ, levels = organ_order)) %>%
  tibble::column_to_rownames("organ") %>%
  as.matrix()

mat_stromal <- mat_stromal[, intersect(stromal_cols_order, colnames(mat_stromal)), drop = FALSE]


ann_stromal <- tibble::tibble(col = colnames(mat_stromal)) %>%
  tidyr::separate(col, c("Macro","Partner"), sep = paste0("\\s*", '↔', "\\s*"), remove = FALSE) %>%
  dplyr::select(-col) %>%                 # 显式 dplyr::
  as.data.frame()
rownames(ann_stromal) <- colnames(mat_stromal)

desired_cols <- tidyr::expand_grid(
  Partner = partners_stromal,   # 先遍历 Partner
  Macro   = incl_macros         # Partner 内的 Macro 顺序
) %>%
  dplyr::mutate(pair = paste(Macro, arrow, Partner)) %>%
  dplyr::pull(pair)

present_cols   <- intersect(desired_cols, colnames(mat_stromal))
mat_stromal_re <- mat_stromal[, present_cols, drop = FALSE]

ann_stromal_re <- ann_stromal[present_cols, , drop = FALSE]

partner_block_sizes <- table(ann_stromal_re$Partner)          # 以当前顺序统计每块列数
gaps_col <- cumsum(as.integer(partner_block_sizes))
if (length(gaps_col) > 0) gaps_col <- gaps_col[-length(gaps_col)]  # 最后一块后面不加

# 5) 重新作图
stromal_hp <- pheatmap(
  mat_stromal_re,
  color = colorRampPalette(c("#F7FBFF","#dcffd7","#ADC6A9","#7DB391","#47A183","#156C73","#1B5364"))(100),
  cluster_rows = FALSE, cluster_cols = FALSE,
  border_color = NA,
  # annotation_col = ann_stromal_re,
  gaps_col = gaps_col,                # ← 分块可视化（可选）
  fontsize_row = 10, fontsize_col = 8, angle_col = 90,
  cellwidth = 16,                 # ← 每格宽度（单位为 mm/像素样式，按经验调）
  cellheight = 20,               # ← 每格高度
  main = "Macrophage subtypes vs Stromal cells"
) # size 10*8


g <- stromal_hp$gtable
mat_id <- which(g$layout$name == "matrix")
g$grobs[[mat_id]] <- grobTree(
  g$grobs[[mat_id]],
  rectGrob(gp = gpar(col = "black", lwd = 0.8, fill = NA))
)

grid.newpage(); grid.draw(g)


### 免疫细胞互作

partners_immune  <- c(incl_macros, "Monocyte", "Neutrophil", "Dendritic cell", "Mast cell", "T/NK cell", "B cell", "Plasma cell")

# 提取与Macro互作为Stromal cell的互作类型
immune_df <- CCI_grp %>%
  filter(Partner %in% partners_immune)
# 剔除Macro亚型的内部互作 
immune_df <- immune_df[which(immune_df$Macro != immune_df$Partner),]


immune_cols_order <- tidyr::expand_grid(
  Macro   = incl_macros,
  Partner = partners_immune
) %>%
  dplyr::filter(Macro != Partner) %>%        # ← 去掉 “X ↔ X”
  dplyr::mutate(pair = paste(Macro, '↔', Partner)) %>%
  dplyr::pull(pair)



mat_immune <- immune_df %>%
  dplyr::mutate(pair = paste(Macro, '↔', Partner)) %>%
  dplyr::select(organ, pair, n) %>%                                            
  tidyr::pivot_wider(names_from = pair, values_from = n, values_fill = 0) %>%
  dplyr::right_join(tibble::tibble(organ = organ_order), by = "organ") %>%     
  dplyr::arrange(factor(organ, levels = organ_order)) %>%
  tibble::column_to_rownames("organ") %>%
  as.matrix()

mat_immune <- mat_immune[, intersect(immune_cols_order, colnames(mat_immune)), drop = FALSE]


ann_immune <- tibble::tibble(col = colnames(mat_immune)) %>%
  tidyr::separate(col, c("Macro","Partner"), sep = paste0("\\s*", '↔', "\\s*"), remove = FALSE) %>%
  dplyr::select(-col) %>%                 # 显式 dplyr::
  as.data.frame()
rownames(ann_immune) <- colnames(mat_immune)

desired_cols <- tidyr::expand_grid(
  Partner = partners_immune,   # 先遍历 Partner
  Macro   = incl_macros         # Partner 内的 Macro 顺序
) %>%
  dplyr::mutate(pair = paste(Macro, '↔', Partner)) %>%
  dplyr::pull(pair)

present_cols   <- intersect(desired_cols, colnames(mat_immune))

mat_immune_re <- mat_immune[, present_cols, drop = FALSE]

ann_immune_re <- ann_immune[present_cols, , drop = FALSE]

partner_block_sizes <- table(ann_immune_re$Partner)          # 以当前顺序统计每块列数
gaps_col <- c(10,15,20,25,30,35,40)

immune_hp <- pheatmap(
  mat_immune_re,
  color = colorRampPalette(c("#F7FBFF","#DBC6D9","#D4A7C4","#D48CB1","#9F6693","#765985","#4F3E6A"))(100),
  cluster_rows = FALSE, cluster_cols = FALSE,
  border_color = NA,
  # annotation_col = ann_stromal_re,
  gaps_col = gaps_col,                # ← 分块可视化（可选）
  fontsize_row = 10, fontsize_col = 8, angle_col = 90,
  cellwidth = 16,                 # ← 每格宽度（单位为 mm/像素样式，按经验调）
  cellheight = 20,               # ← 每格高度
  main = "Macrophage subtypes vs Immune cells"
) # size 10*8

g <- immune_hp$gtable
mat_id <- which(g$layout$name == "matrix")
g$grobs[[mat_id]] <- grobTree(
  g$grobs[[mat_id]],
  rectGrob(gp = gpar(col = "black", lwd = 0.8, fill = NA))
)

grid.newpage(); grid.draw(g)




#### Macro to T/NK cell L-R pairs分析 -----------------------------------------------------------------------------------

## load data
root <- "/storage/data/KAI/Pan_Inflammation/Homo/Mac2025/analysis/CellPhoneDB/result_new"
org_dirs <- fs::dir_ls(root, type = "directory")


incl_macros <- c(
  "IL1B+NLRP3+_Macro",
  "SPP1+CCL2+_Macro",
  "CD1C+_DC-like_Macro",
  "CCL13+_Complement-associated_Macro",
  "MMP3+CXCL8+_Macro"
)

target_immune <- "T/NK cell"

organ_names <- basename(org_dirs) |>
  str_replace("_inflamed$","")

read_one_organ <- function(dir_path, organ_label) {
  # 精确匹配文件名（只取 overall means，不取 significant_means）
  f_pval <- list.files(dir_path, "^statistical_analysis_pvalues_.*\\.txt$", full.names = TRUE)
  f_mean <- list.files(dir_path, "^statistical_analysis_means_.*\\.txt$",  full.names = TRUE)
  
  # 若同类出现多个文件，取最新一个
  pick_latest <- function(paths) if (length(paths) <= 1) paths else paths[which.max(file.info(paths)$mtime)]
  if (length(f_pval) != 1) f_pval <- pick_latest(f_pval)
  if (length(f_mean) != 1) f_mean <- pick_latest(f_mean)
  
  # 若缺文件，返回空 tibble
  if (length(f_pval) == 0 || length(f_mean) == 0) {
    message("Skip (missing files): ", organ_label)
    return(tibble(
      organ=character(), interacting_pair=character(),
      gene_a=character(), gene_b=character(),
      A=character(), B=character(),
      pval=double(), log10p=double(), expr_mean=double(),
      LR=character()
    ))
  }
  
  pval <- read_tsv(f_pval, show_col_types = FALSE)
  mean <- read_tsv(f_mean, show_col_types = FALSE)
  
  # 找到 "clusterA|clusterB" 的列
  pair_cols_p <- grep("\\|", names(pval), value = TRUE)
  pair_cols_m <- grep("\\|", names(mean), value = TRUE)
  if (length(pair_cols_p) == 0 || length(pair_cols_m) == 0) {
    message("No pair columns in: ", organ_label)
    return(tibble(
      organ=character(), interacting_pair=character(),
      gene_a=character(), gene_b=character(),
      A=character(), B=character(),
      pval=double(), log10p=double(), expr_mean=double(),
      LR=character()
    ))
  }
  
  # ---- pvalues: 拉长 + 让键唯一（对重复键取最小 p） ----
  p_long <- pval %>%
    pivot_longer(all_of(pair_cols_p), names_to = "pair", values_to = "pval") %>%
    separate(pair, c("A","B"), sep="\\|", fill="right") %>%
    mutate(A = str_squish(A), B = str_squish(B)) %>%
    # 保留 gene_a/gene_b（若不存在也没关系）
    dplyr::select(interacting_pair,
           gene_a = any_of("gene_a"),
           gene_b = any_of("gene_b"),
           A, B, pval) %>%
    group_by(interacting_pair, A, B) %>%
    summarise(
      pval   = suppressWarnings(min(as.numeric(pval), na.rm = TRUE)),
      gene_a = dplyr::first(gene_a),
      gene_b = dplyr::first(gene_b),
      .groups = "drop"
    )
  
  # ---- means: 拉长 + 让键唯一（对重复键取均值） ----
  m_long <- mean %>%
    pivot_longer(all_of(pair_cols_m), names_to = "pair", values_to = "expr_mean") %>%
    separate(pair, c("A","B"), sep="\\|", fill="right") %>%
    mutate(A = str_squish(A), B = str_squish(B)) %>%
    dplyr::select(interacting_pair,
           gene_a_m = any_of("gene_a"),
           gene_b_m = any_of("gene_b"),
           A, B, expr_mean) %>%
    group_by(interacting_pair, A, B) %>%
    summarise(
      expr_mean = suppressWarnings(mean(as.numeric(expr_mean), na.rm = TRUE)),
      gene_a_m  = dplyr::first(gene_a_m),
      gene_b_m  = dplyr::first(gene_b_m),
      .groups = "drop"
    )
  
  # ---- 合并（此时每个键唯一，不再 many-to-many）----
  out <- p_long %>%
    left_join(m_long, by = c("interacting_pair","A","B")) %>%
    mutate(
      organ   = organ_label,
      # 兜底基因名：优先用 pvals 的，若无，用 means 的；再无，用 interacting_pair
      gene_a2 = dplyr::coalesce(gene_a, gene_a_m),
      gene_b2 = dplyr::coalesce(gene_b, gene_b_m),
      LR      = dplyr::if_else(!is.na(gene_a2) & !is.na(gene_b2),
                               paste0(gene_a2, " – ", gene_b2),
                               interacting_pair),
      pval    = as.numeric(pval),
      log10p  = -log10(pmax(pval, 1e-300))
    ) %>%
    dplyr::select(organ, interacting_pair, gene_a = gene_a2, gene_b = gene_b2,
           A, B, pval, log10p, expr_mean, LR)
  
  out
}

df_all <- purrr::map2_dfr(org_dirs, organ_names, read_one_organ)

# 判定是否为 CellPhoneDB 里那类“复合体”命名
.is_complex <- function(x) {
  !is.na(x) & str_detect(x, "(?i)(^integrin_|_complex\\b)")
  # 说明：
  #  - ^integrin_    : 以 integrin_ 开头（CPDB 常见）
  #  - _complex\b    : 以 _complex 结尾（或单词边界），如 *_complex
}

# 过滤：去掉 gene_a/gene_b 缺失 & 去掉复合体命名
drop_complex_or_na <- function(df) {
  df %>%
    # 1) 基因名必须都有
    filter(!is.na(gene_a), !is.na(gene_b)) %>%
    # 2) 两侧都不能是“复合体”命名
    filter(!.is_complex(gene_a),
           !.is_complex(gene_b)) %>%
    # 3) 兜底：如果某些版本没有 gene_a/gene_b，或者列名不全，
    #    但 interacting_pair 里仍带有 complex，也顺手去掉
    filter(!str_detect(interacting_pair, "(?i)integrin_.*_complex|_complex\\b"))
}
df_all_clean <- df_all %>% drop_complex_or_na()

## 卡显著表达L-R pairs
df_mt <- df_all_clean %>%
  filter(A %in% incl_macros, B == "T/NK cell") %>%  # 你的方向性筛选
  filter(pval < 0.05) 


topN <- 20
order_by_macro <- df_mt %>%
  group_by(A, LR) %>%
  summarise(
    total_log10p = sum(log10p, na.rm = TRUE),
    mean_expr     = mean(expr_mean, na.rm = TRUE),
    n_organs      = n_distinct(organ),
    .groups = "drop"
  ) %>%
  arrange(A, desc(total_log10p), desc(mean_expr), desc(n_organs)) %>%
  group_by(A) %>% slice_head(n = topN) %>% ungroup() %>%
  mutate(LR_f = factor(LR, levels = unique(LR)))

df_mt <- df_mt %>%
  semi_join(order_by_macro, by = c("A","LR")) %>%
  left_join(order_by_macro, by = c("A","LR"))


df_mt$A <- factor(df_mt$A, levels = incl_macros)


ggplot(df_mt, aes(x = organ, y = LR_f)) +
  geom_point(aes(#size = log10p, 
                 colour = expr_mean), 
             alpha = 1, 
             na.rm = TRUE) +
  # scale_size_area(name = expression(-log[10](p)), max_size = 6,
  #                 breaks = c(1,2,3,5,10)) +
  scale_colour_gradientn(colors = c('#77C3BC','#4EA9C6','#2B84B2','#0E63A0',
                                    '#073E79','#0C203B'), 
                         name = "Expression (mean)") +
  facet_grid(A ~ ., scales = "free_y", space = "free_y",
             labeller = labeller(A = label_value)) +
  labs(x = NULL, y = NULL,
       title = "Macrophage subtypes (sender) → T/NK cell (receiver)",
       # subtitle = "Colour: mean expression; Size: -log10(p-value); p < 0.05"
       ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid = element_blank(), 
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text.y = element_text(face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.6),
    axis.text.x = element_text(angle = 90,hjust = 1),
    legend.position = "right"
  )


#### Macro to Monocyte L-R pairs分析 -----------------------------------------------------------------------------------

## load data
root <- "/storage/data/KAI/Pan_Inflammation/Homo/Mac2025/analysis/CellPhoneDB/result_new"
org_dirs <- fs::dir_ls(root, type = "directory")


incl_macros <- c(
  "IL1B+NLRP3+_Macro",
  "SPP1+CCL2+_Macro",
  "CD1C+_DC-like_Macro",
  "CCL13+_Complement-associated_Macro",
  "MMP3+CXCL8+_Macro"
)

target_immune <- "Monocyte"

organ_names <- basename(org_dirs) |>
  str_replace("_inflamed$","")

read_one_organ <- function(dir_path, organ_label) {
  # 精确匹配文件名（只取 overall means，不取 significant_means）
  f_pval <- list.files(dir_path, "^statistical_analysis_pvalues_.*\\.txt$", full.names = TRUE)
  f_mean <- list.files(dir_path, "^statistical_analysis_means_.*\\.txt$",  full.names = TRUE)
  
  # 若同类出现多个文件，取最新一个
  pick_latest <- function(paths) if (length(paths) <= 1) paths else paths[which.max(file.info(paths)$mtime)]
  if (length(f_pval) != 1) f_pval <- pick_latest(f_pval)
  if (length(f_mean) != 1) f_mean <- pick_latest(f_mean)
  
  # 若缺文件，返回空 tibble
  if (length(f_pval) == 0 || length(f_mean) == 0) {
    message("Skip (missing files): ", organ_label)
    return(tibble(
      organ=character(), interacting_pair=character(),
      gene_a=character(), gene_b=character(),
      A=character(), B=character(),
      pval=double(), log10p=double(), expr_mean=double(),
      LR=character()
    ))
  }
  
  pval <- read_tsv(f_pval, show_col_types = FALSE)
  mean <- read_tsv(f_mean, show_col_types = FALSE)
  
  # 找到 "clusterA|clusterB" 的列
  pair_cols_p <- grep("\\|", names(pval), value = TRUE)
  pair_cols_m <- grep("\\|", names(mean), value = TRUE)
  if (length(pair_cols_p) == 0 || length(pair_cols_m) == 0) {
    message("No pair columns in: ", organ_label)
    return(tibble(
      organ=character(), interacting_pair=character(),
      gene_a=character(), gene_b=character(),
      A=character(), B=character(),
      pval=double(), log10p=double(), expr_mean=double(),
      LR=character()
    ))
  }
  
  # ---- pvalues: 拉长 + 让键唯一（对重复键取最小 p） ----
  p_long <- pval %>%
    pivot_longer(all_of(pair_cols_p), names_to = "pair", values_to = "pval") %>%
    separate(pair, c("A","B"), sep="\\|", fill="right") %>%
    mutate(A = str_squish(A), B = str_squish(B)) %>%
    # 保留 gene_a/gene_b（若不存在也没关系）
    dplyr::select(interacting_pair,
                  gene_a = any_of("gene_a"),
                  gene_b = any_of("gene_b"),
                  A, B, pval) %>%
    group_by(interacting_pair, A, B) %>%
    summarise(
      pval   = suppressWarnings(min(as.numeric(pval), na.rm = TRUE)),
      gene_a = dplyr::first(gene_a),
      gene_b = dplyr::first(gene_b),
      .groups = "drop"
    )
  
  # ---- means: 拉长 + 让键唯一（对重复键取均值） ----
  m_long <- mean %>%
    pivot_longer(all_of(pair_cols_m), names_to = "pair", values_to = "expr_mean") %>%
    separate(pair, c("A","B"), sep="\\|", fill="right") %>%
    mutate(A = str_squish(A), B = str_squish(B)) %>%
    dplyr::select(interacting_pair,
                  gene_a_m = any_of("gene_a"),
                  gene_b_m = any_of("gene_b"),
                  A, B, expr_mean) %>%
    group_by(interacting_pair, A, B) %>%
    summarise(
      expr_mean = suppressWarnings(mean(as.numeric(expr_mean), na.rm = TRUE)),
      gene_a_m  = dplyr::first(gene_a_m),
      gene_b_m  = dplyr::first(gene_b_m),
      .groups = "drop"
    )
  
  # ---- 合并（此时每个键唯一，不再 many-to-many）----
  out <- p_long %>%
    left_join(m_long, by = c("interacting_pair","A","B")) %>%
    mutate(
      organ   = organ_label,
      # 兜底基因名：优先用 pvals 的，若无，用 means 的；再无，用 interacting_pair
      gene_a2 = dplyr::coalesce(gene_a, gene_a_m),
      gene_b2 = dplyr::coalesce(gene_b, gene_b_m),
      LR      = dplyr::if_else(!is.na(gene_a2) & !is.na(gene_b2),
                               paste0(gene_a2, " – ", gene_b2),
                               interacting_pair),
      pval    = as.numeric(pval),
      log10p  = -log10(pmax(pval, 1e-300))
    ) %>%
    dplyr::select(organ, interacting_pair, gene_a = gene_a2, gene_b = gene_b2,
                  A, B, pval, log10p, expr_mean, LR)
  
  out
}

df_all <- purrr::map2_dfr(org_dirs, organ_names, read_one_organ)

# 判定是否为 CellPhoneDB 里那类“复合体”命名
.is_complex <- function(x) {
  !is.na(x) & str_detect(x, "(?i)(^integrin_|_complex\\b)")
  # 说明：
  #  - ^integrin_    : 以 integrin_ 开头（CPDB 常见）
  #  - _complex\b    : 以 _complex 结尾（或单词边界），如 *_complex
}

# 过滤：去掉 gene_a/gene_b 缺失 & 去掉复合体命名
drop_complex_or_na <- function(df) {
  df %>%
    # 1) 基因名必须都有
    filter(!is.na(gene_a), !is.na(gene_b)) %>%
    # 2) 两侧都不能是“复合体”命名
    filter(!.is_complex(gene_a),
           !.is_complex(gene_b)) %>%
    # 3) 兜底：如果某些版本没有 gene_a/gene_b，或者列名不全，
    #    但 interacting_pair 里仍带有 complex，也顺手去掉
    filter(!str_detect(interacting_pair, "(?i)integrin_.*_complex|_complex\\b"))
}
df_all_clean <- df_all %>% drop_complex_or_na()

## 卡显著表达L-R pairs
df_mt <- df_all_clean %>%
  filter(A %in% incl_macros, B == target_immune) %>%  # 你的方向性筛选
  filter(pval < 0.05) 

TopN_per_organ <- 10
Final_cap <- 20

# 每个巨噬细胞A -> target_immune 取各器官top10
top_per_org <- df_mt %>%
  group_by(organ, A) %>%
  arrange(desc(expr_mean), .by_group = TRUE) %>%
  slice_head(n = TopN_per_organ) %>%
  ungroup()

# 各器官top10取并集
# 4) 对每个巨噬亚型 A：汇总“并集”并二次排序（频次→平均表达→统计强度）
rank_union_per_A <- top_per_org %>%
  group_by(A, organ, LR) %>%
  summarise(expr_mean = max(expr_mean, na.rm = TRUE),
            log10p = max(log10p, na.rm = TRUE),
            .groups = "drop") %>%
  group_by(A, LR) %>%
  summarise(
    freq_org   = n_distinct(organ),          # 跨器官出现频次
    mean_expr  = mean(expr_mean, na.rm = TRUE),
    sum_log10p = sum(log10p, na.rm = TRUE),  # 累积统计强度（可换成 mean 或 max）
    .groups = "drop"
  ) %>%
  arrange(A, desc(freq_org), desc(mean_expr), desc(sum_log10p)) %>%
  group_by(A) %>%
  mutate(rank = row_number()) %>%
  ungroup()

# 超过20选top20
keep_LR <- rank_union_per_A %>%
  group_by(A) %>%
  slice_head(n = Final_cap) %>%
  ungroup() %>%
  dplyr::select(A, LR)

# 制成绘图dataframe
plot_df <- df_mt %>%
  semi_join(keep_LR, by = c("A","LR")) %>%
  group_by(organ, A, LR) %>%
  summarise(
    expr_mean = mean(expr_mean, na.rm = TRUE),
    log10p    = max(log10p, na.rm = TRUE),   # 多条同 LR 取更强的统计
    .groups = "drop"
  ) %>%
  complete(organ = organ_names, nesting(A, LR), fill = list(expr_mean = 0, log10p = 0)) %>%
  mutate(
    organ = factor(organ, levels = organ_names),
    ## Row 排序：按 A 分面内，使用 4) 中的 rank
    LR = factor(LR, levels = rank_union_per_A %>%
                  semi_join(keep_LR, by = c("A","LR")) %>%
                  arrange(A, rank) %>%
                  pull(LR) %>% unique())
  )


eps <- 1e-8  # 判定“为 0”的容差

## 先找出“有任何非零表达”的器官
valid_orgs <- plot_df %>%
  dplyr::group_by(organ) %>%
  dplyr::summarise(any_pos = any(expr_mean > eps, na.rm = TRUE), .groups = "drop") %>%
  dplyr::filter(any_pos) %>%
  dplyr::pull(organ)

plot_df2 <- plot_df %>%
  dplyr::filter(organ %in% valid_orgs, expr_mean > eps) %>%
  dplyr::mutate(
    # 重新设定因子顺序，防止掉列后顺序错乱
    organ = factor(organ, levels = intersect(levels(plot_df$organ), valid_orgs))
  ) %>%
  droplevels()  # 清理掉未用到的水平

plot_df2$A <- factor(plot_df2$A, levels = incl_macros)

ggplot(plot_df2, aes(x = organ, y = LR)) +
  geom_point(aes(# size = log10p, 
                 colour = expr_mean),
             alpha = 1, na.rm = TRUE) +
  # scale_size_area(name = expression(-log[10](p)), max_size = 6,
  #                 breaks = c(1,2,3,5,10), limits = c(0, NA)) +
  scale_colour_gradientn(colors = c('white','#77C3BC','#4EA9C6','#2B84B2','#0E63A0',
                                    '#073E79','#0C203B'), 
                         name = "Expression (mean)") +
  facet_grid(A ~ ., scales = "free_y", space = "free_y") +
  labs(x = NULL, y = NULL,
       title    = "Macrophage subtypes (sender) → Monocyte (receiver)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid = element_blank(), 
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text.y  = element_text(face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    axis.text.x   = element_text(angle = 0, hjust = 0.5),
    legend.position = "right"
  ) # size 8*12


#### Macro to T/NK cell L-R pairs分析 ---------------------------------------------------
## load data
root <- "/storage/data/KAI/Pan_Inflammation/Homo/Mac2025/analysis/CellPhoneDB/result_new"
org_dirs <- fs::dir_ls(root, type = "directory")

incl_macros <- c(
  "IL1B+NLRP3+_Macro",
  "SPP1+CCL2+_Macro",
  "CD1C+_DC-like_Macro",
  "CCL13+_Complement-associated_Macro",
  "MMP3+CXCL8+_Macro"
)

target_immune <- "T/NK cell"

organ_names <- basename(org_dirs) |>
  str_replace("_inflamed$","")

read_one_organ <- function(dir_path, organ_label) {
  # 精确匹配文件名（只取 overall means，不取 significant_means）
  f_pval <- list.files(dir_path, "^statistical_analysis_pvalues_.*\\.txt$", full.names = TRUE)
  f_mean <- list.files(dir_path, "^statistical_analysis_means_.*\\.txt$",  full.names = TRUE)
  
  # 若同类出现多个文件，取最新一个
  pick_latest <- function(paths) if (length(paths) <= 1) paths else paths[which.max(file.info(paths)$mtime)]
  if (length(f_pval) != 1) f_pval <- pick_latest(f_pval)
  if (length(f_mean) != 1) f_mean <- pick_latest(f_mean)
  
  # 若缺文件，返回空 tibble
  if (length(f_pval) == 0 || length(f_mean) == 0) {
    message("Skip (missing files): ", organ_label)
    return(tibble(
      organ=character(), interacting_pair=character(),
      gene_a=character(), gene_b=character(),
      A=character(), B=character(),
      pval=double(), log10p=double(), expr_mean=double(),
      LR=character()
    ))
  }
  
  pval <- read_tsv(f_pval, show_col_types = FALSE)
  mean <- read_tsv(f_mean, show_col_types = FALSE)
  
  # 找到 "clusterA|clusterB" 的列
  pair_cols_p <- grep("\\|", names(pval), value = TRUE)
  pair_cols_m <- grep("\\|", names(mean), value = TRUE)
  if (length(pair_cols_p) == 0 || length(pair_cols_m) == 0) {
    message("No pair columns in: ", organ_label)
    return(tibble(
      organ=character(), interacting_pair=character(),
      gene_a=character(), gene_b=character(),
      A=character(), B=character(),
      pval=double(), log10p=double(), expr_mean=double(),
      LR=character()
    ))
  }
  
  # ---- pvalues: 拉长 + 让键唯一（对重复键取最小 p） ----
  p_long <- pval %>%
    pivot_longer(all_of(pair_cols_p), names_to = "pair", values_to = "pval") %>%
    separate(pair, c("A","B"), sep="\\|", fill="right") %>%
    mutate(A = str_squish(A), B = str_squish(B)) %>%
    # 保留 gene_a/gene_b（若不存在也没关系）
    dplyr::select(interacting_pair,
                  gene_a = any_of("gene_a"),
                  gene_b = any_of("gene_b"),
                  A, B, pval) %>%
    group_by(interacting_pair, A, B) %>%
    summarise(
      pval   = suppressWarnings(min(as.numeric(pval), na.rm = TRUE)),
      gene_a = dplyr::first(gene_a),
      gene_b = dplyr::first(gene_b),
      .groups = "drop"
    )
  
  # ---- means: 拉长 + 让键唯一（对重复键取均值） ----
  m_long <- mean %>%
    pivot_longer(all_of(pair_cols_m), names_to = "pair", values_to = "expr_mean") %>%
    separate(pair, c("A","B"), sep="\\|", fill="right") %>%
    mutate(A = str_squish(A), B = str_squish(B)) %>%
    dplyr::select(interacting_pair,
                  gene_a_m = any_of("gene_a"),
                  gene_b_m = any_of("gene_b"),
                  A, B, expr_mean) %>%
    group_by(interacting_pair, A, B) %>%
    summarise(
      expr_mean = suppressWarnings(mean(as.numeric(expr_mean), na.rm = TRUE)),
      gene_a_m  = dplyr::first(gene_a_m),
      gene_b_m  = dplyr::first(gene_b_m),
      .groups = "drop"
    )
  
  # ---- 合并（此时每个键唯一，不再 many-to-many）----
  out <- p_long %>%
    left_join(m_long, by = c("interacting_pair","A","B")) %>%
    mutate(
      organ   = organ_label,
      # 兜底基因名：优先用 pvals 的，若无，用 means 的；再无，用 interacting_pair
      gene_a2 = dplyr::coalesce(gene_a, gene_a_m),
      gene_b2 = dplyr::coalesce(gene_b, gene_b_m),
      LR      = dplyr::if_else(!is.na(gene_a2) & !is.na(gene_b2),
                               paste0(gene_a2, " – ", gene_b2),
                               interacting_pair),
      pval    = as.numeric(pval),
      log10p  = -log10(pmax(pval, 1e-300))
    ) %>%
    dplyr::select(organ, interacting_pair, gene_a = gene_a2, gene_b = gene_b2,
                  A, B, pval, log10p, expr_mean, LR)
  
  out
}

df_all <- purrr::map2_dfr(org_dirs, organ_names, read_one_organ)

# 判定是否为 CellPhoneDB 里那类“复合体”命名
.is_complex <- function(x) {
  !is.na(x) & str_detect(x, "(?i)(^integrin_|_complex\\b)")
  # 说明：
  #  - ^integrin_    : 以 integrin_ 开头（CPDB 常见）
  #  - _complex\b    : 以 _complex 结尾（或单词边界），如 *_complex
}

# 过滤：去掉 gene_a/gene_b 缺失 & 去掉复合体命名
drop_complex_or_na <- function(df) {
  df %>%
    # 1) 基因名必须都有
    filter(!is.na(gene_a), !is.na(gene_b)) %>%
    # 2) 两侧都不能是“复合体”命名
    filter(!.is_complex(gene_a),
           !.is_complex(gene_b)) %>%
    # 3) 兜底：如果某些版本没有 gene_a/gene_b，或者列名不全，
    #    但 interacting_pair 里仍带有 complex，也顺手去掉
    filter(!str_detect(interacting_pair, "(?i)integrin_.*_complex|_complex\\b"))
}
df_all_clean <- df_all %>% drop_complex_or_na()

## 卡显著表达L-R pairs
df_mt <- df_all_clean %>%
  filter(A %in% incl_macros, B == target_immune) %>%  # 你的方向性筛选
  filter(pval < 0.05) 

TopN_per_organ <- 20
Final_cap <- 20

# 每个巨噬细胞A -> target_immune 取各器官top10
top_per_org <- df_mt %>%
  group_by(organ, A) %>%
  arrange(desc(expr_mean), .by_group = TRUE) %>%
  slice_head(n = TopN_per_organ) %>%
  ungroup()

# 各器官top10取并集
# 4) 对每个巨噬亚型 A：汇总“并集”并二次排序（频次→平均表达→统计强度）
rank_union_per_A <- top_per_org %>%
  group_by(A, organ, LR) %>%
  summarise(expr_mean = max(expr_mean, na.rm = TRUE),
            log10p = max(log10p, na.rm = TRUE),
            .groups = "drop") %>%
  group_by(A, LR) %>%
  summarise(
    freq_org   = n_distinct(organ),          # 跨器官出现频次
    mean_expr  = mean(expr_mean, na.rm = TRUE),
    sum_log10p = sum(log10p, na.rm = TRUE),  # 累积统计强度（可换成 mean 或 max）
    .groups = "drop"
  ) %>%
  arrange(A, desc(freq_org), desc(mean_expr), desc(sum_log10p)) %>%
  group_by(A) %>%
  mutate(rank = row_number()) %>%
  ungroup()

# 超过20选top20
keep_LR <- rank_union_per_A %>%
  group_by(A) %>%
  slice_head(n = Final_cap) %>%
  ungroup() %>%
  dplyr::select(A, LR)

# 制成绘图dataframe
plot_df <- df_mt %>%
  semi_join(keep_LR, by = c("A","LR")) %>%
  group_by(organ, A, LR) %>%
  summarise(
    expr_mean = mean(expr_mean, na.rm = TRUE),
    log10p    = max(log10p, na.rm = TRUE),   # 多条同 LR 取更强的统计
    .groups = "drop"
  ) %>%
  complete(organ = organ_names, nesting(A, LR), fill = list(expr_mean = 0, log10p = 0)) %>%
  mutate(
    organ = factor(organ, levels = organ_names),
    ## Row 排序：按 A 分面内，使用 4) 中的 rank
    LR = factor(LR, levels = rank_union_per_A %>%
                  semi_join(keep_LR, by = c("A","LR")) %>%
                  arrange(A, rank) %>%
                  pull(LR) %>% unique())
  )


eps <- 1e-8  # 判定“为 0”的容差

## 先找出“有任何非零表达”的器官
valid_orgs <- plot_df %>%
  dplyr::group_by(organ) %>%
  dplyr::summarise(any_pos = any(expr_mean > eps, na.rm = TRUE), .groups = "drop") %>%
  dplyr::filter(any_pos) %>%
  dplyr::pull(organ)

plot_df2 <- plot_df %>%
  dplyr::filter(organ %in% valid_orgs, expr_mean > eps) %>%
  dplyr::mutate(
    # 重新设定因子顺序，防止掉列后顺序错乱
    organ = factor(organ, levels = intersect(levels(plot_df$organ), valid_orgs))
  ) %>%
  droplevels()  # 清理掉未用到的水平

plot_df2$A <- factor(plot_df2$A, levels = incl_macros)

ggplot(plot_df2, aes(x = organ, y = LR)) +
  geom_point(aes(# size = log10p, 
    colour = expr_mean),
    alpha = 1, na.rm = TRUE) +
  # scale_size_area(name = expression(-log[10](p)), max_size = 6,
  #                 breaks = c(1,2,3,5,10), limits = c(0, NA)) +
  scale_colour_gradientn(colors = c('white','#77C3BC','#4EA9C6','#2B84B2','#0E63A0',
                                    '#073E79','#0C203B'), 
                         name = "Expression (mean)") +
  facet_grid(A ~ ., scales = "free_y", space = "free_y") +
  labs(x = NULL, y = NULL,
       title    = "Macrophage subtypes (sender) → T/NK cell (receiver)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid = element_blank(), 
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text.y  = element_text(face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    axis.text.x   = element_text(angle = 0, hjust = 0.5),
    legend.position = "right"
  ) # size 8*12

#### Macro to B cell L-R pairs分析 ---------------------------------------------------------
## load data
root <- "/storage/data/KAI/Pan_Inflammation/Homo/Mac2025/analysis/CellPhoneDB/result_new"
org_dirs <- fs::dir_ls(root, type = "directory")

incl_macros <- c(
  "IL1B+NLRP3+_Macro",
  "SPP1+CCL2+_Macro",
  "CD1C+_DC-like_Macro",
  "CCL13+_Complement-associated_Macro",
  "MMP3+CXCL8+_Macro"
)

target_immune <- "B cell"

organ_names <- basename(org_dirs) |>
  str_replace("_inflamed$","")

read_one_organ <- function(dir_path, organ_label) {
  # 精确匹配文件名（只取 overall means，不取 significant_means）
  f_pval <- list.files(dir_path, "^statistical_analysis_pvalues_.*\\.txt$", full.names = TRUE)
  f_mean <- list.files(dir_path, "^statistical_analysis_means_.*\\.txt$",  full.names = TRUE)
  
  # 若同类出现多个文件，取最新一个
  pick_latest <- function(paths) if (length(paths) <= 1) paths else paths[which.max(file.info(paths)$mtime)]
  if (length(f_pval) != 1) f_pval <- pick_latest(f_pval)
  if (length(f_mean) != 1) f_mean <- pick_latest(f_mean)
  
  # 若缺文件，返回空 tibble
  if (length(f_pval) == 0 || length(f_mean) == 0) {
    message("Skip (missing files): ", organ_label)
    return(tibble(
      organ=character(), interacting_pair=character(),
      gene_a=character(), gene_b=character(),
      A=character(), B=character(),
      pval=double(), log10p=double(), expr_mean=double(),
      LR=character()
    ))
  }
  
  pval <- read_tsv(f_pval, show_col_types = FALSE)
  mean <- read_tsv(f_mean, show_col_types = FALSE)
  
  # 找到 "clusterA|clusterB" 的列
  pair_cols_p <- grep("\\|", names(pval), value = TRUE)
  pair_cols_m <- grep("\\|", names(mean), value = TRUE)
  if (length(pair_cols_p) == 0 || length(pair_cols_m) == 0) {
    message("No pair columns in: ", organ_label)
    return(tibble(
      organ=character(), interacting_pair=character(),
      gene_a=character(), gene_b=character(),
      A=character(), B=character(),
      pval=double(), log10p=double(), expr_mean=double(),
      LR=character()
    ))
  }
  
  # ---- pvalues: 拉长 + 让键唯一（对重复键取最小 p） ----
  p_long <- pval %>%
    pivot_longer(all_of(pair_cols_p), names_to = "pair", values_to = "pval") %>%
    separate(pair, c("A","B"), sep="\\|", fill="right") %>%
    mutate(A = str_squish(A), B = str_squish(B)) %>%
    # 保留 gene_a/gene_b（若不存在也没关系）
    dplyr::select(interacting_pair,
                  gene_a = any_of("gene_a"),
                  gene_b = any_of("gene_b"),
                  A, B, pval) %>%
    group_by(interacting_pair, A, B) %>%
    summarise(
      pval   = suppressWarnings(min(as.numeric(pval), na.rm = TRUE)),
      gene_a = dplyr::first(gene_a),
      gene_b = dplyr::first(gene_b),
      .groups = "drop"
    )
  
  # ---- means: 拉长 + 让键唯一（对重复键取均值） ----
  m_long <- mean %>%
    pivot_longer(all_of(pair_cols_m), names_to = "pair", values_to = "expr_mean") %>%
    separate(pair, c("A","B"), sep="\\|", fill="right") %>%
    mutate(A = str_squish(A), B = str_squish(B)) %>%
    dplyr::select(interacting_pair,
                  gene_a_m = any_of("gene_a"),
                  gene_b_m = any_of("gene_b"),
                  A, B, expr_mean) %>%
    group_by(interacting_pair, A, B) %>%
    summarise(
      expr_mean = suppressWarnings(mean(as.numeric(expr_mean), na.rm = TRUE)),
      gene_a_m  = dplyr::first(gene_a_m),
      gene_b_m  = dplyr::first(gene_b_m),
      .groups = "drop"
    )
  
  # ---- 合并（此时每个键唯一，不再 many-to-many）----
  out <- p_long %>%
    left_join(m_long, by = c("interacting_pair","A","B")) %>%
    mutate(
      organ   = organ_label,
      # 兜底基因名：优先用 pvals 的，若无，用 means 的；再无，用 interacting_pair
      gene_a2 = dplyr::coalesce(gene_a, gene_a_m),
      gene_b2 = dplyr::coalesce(gene_b, gene_b_m),
      LR      = dplyr::if_else(!is.na(gene_a2) & !is.na(gene_b2),
                               paste0(gene_a2, " – ", gene_b2),
                               interacting_pair),
      pval    = as.numeric(pval),
      log10p  = -log10(pmax(pval, 1e-300))
    ) %>%
    dplyr::select(organ, interacting_pair, gene_a = gene_a2, gene_b = gene_b2,
                  A, B, pval, log10p, expr_mean, LR)
  
  out
}

df_all <- purrr::map2_dfr(org_dirs, organ_names, read_one_organ)

# 判定是否为 CellPhoneDB 里那类“复合体”命名
.is_complex <- function(x) {
  !is.na(x) & str_detect(x, "(?i)(^integrin_|_complex\\b)")
  # 说明：
  #  - ^integrin_    : 以 integrin_ 开头（CPDB 常见）
  #  - _complex\b    : 以 _complex 结尾（或单词边界），如 *_complex
}

# 过滤：去掉 gene_a/gene_b 缺失 & 去掉复合体命名
drop_complex_or_na <- function(df) {
  df %>%
    # 1) 基因名必须都有
    filter(!is.na(gene_a), !is.na(gene_b)) %>%
    # 2) 两侧都不能是“复合体”命名
    filter(!.is_complex(gene_a),
           !.is_complex(gene_b)) %>%
    # 3) 兜底：如果某些版本没有 gene_a/gene_b，或者列名不全，
    #    但 interacting_pair 里仍带有 complex，也顺手去掉
    filter(!str_detect(interacting_pair, "(?i)integrin_.*_complex|_complex\\b"))
}
df_all_clean <- df_all %>% drop_complex_or_na()

## 卡显著表达L-R pairs
df_mt <- df_all_clean %>%
  filter(A %in% incl_macros, B == target_immune) %>%  # 你的方向性筛选
  filter(pval < 0.05) 

TopN_per_organ <- 10
Final_cap <- 20

# 每个巨噬细胞A -> target_immune 取各器官top10
top_per_org <- df_mt %>%
  group_by(organ, A) %>%
  arrange(desc(expr_mean), .by_group = TRUE) %>%
  slice_head(n = TopN_per_organ) %>%
  ungroup()

# 各器官top10取并集
# 4) 对每个巨噬亚型 A：汇总“并集”并二次排序（频次→平均表达→统计强度）
rank_union_per_A <- top_per_org %>%
  group_by(A, organ, LR) %>%
  summarise(expr_mean = max(expr_mean, na.rm = TRUE),
            log10p = max(log10p, na.rm = TRUE),
            .groups = "drop") %>%
  group_by(A, LR) %>%
  summarise(
    freq_org   = n_distinct(organ),          # 跨器官出现频次
    mean_expr  = mean(expr_mean, na.rm = TRUE),
    sum_log10p = sum(log10p, na.rm = TRUE),  # 累积统计强度（可换成 mean 或 max）
    .groups = "drop"
  ) %>%
  arrange(A, desc(freq_org), desc(mean_expr), desc(sum_log10p)) %>%
  group_by(A) %>%
  mutate(rank = row_number()) %>%
  ungroup()

# 超过20选top20
keep_LR <- rank_union_per_A %>%
  group_by(A) %>%
  slice_head(n = Final_cap) %>%
  ungroup() %>%
  dplyr::select(A, LR)

# 制成绘图dataframe
plot_df <- df_mt %>%
  semi_join(keep_LR, by = c("A","LR")) %>%
  group_by(organ, A, LR) %>%
  summarise(
    expr_mean = mean(expr_mean, na.rm = TRUE),
    log10p    = max(log10p, na.rm = TRUE),   # 多条同 LR 取更强的统计
    .groups = "drop"
  ) %>%
  complete(organ = organ_names, nesting(A, LR), fill = list(expr_mean = 0, log10p = 0)) %>%
  mutate(
    organ = factor(organ, levels = organ_names),
    ## Row 排序：按 A 分面内，使用 4) 中的 rank
    LR = factor(LR, levels = rank_union_per_A %>%
                  semi_join(keep_LR, by = c("A","LR")) %>%
                  arrange(A, rank) %>%
                  pull(LR) %>% unique())
  )


eps <- 1e-8  # 判定“为 0”的容差

## 先找出“有任何非零表达”的器官
valid_orgs <- plot_df %>%
  dplyr::group_by(organ) %>%
  dplyr::summarise(any_pos = any(expr_mean > eps, na.rm = TRUE), .groups = "drop") %>%
  dplyr::filter(any_pos) %>%
  dplyr::pull(organ)

plot_df2 <- plot_df %>%
  dplyr::filter(organ %in% valid_orgs, expr_mean > eps) %>%
  dplyr::mutate(
    # 重新设定因子顺序，防止掉列后顺序错乱
    organ = factor(organ, levels = intersect(levels(plot_df$organ), valid_orgs))
  ) %>%
  droplevels()  # 清理掉未用到的水平

plot_df2$A <- factor(plot_df2$A, levels = incl_macros)

ggplot(plot_df2, aes(x = organ, y = LR)) +
  geom_point(aes(# size = log10p, 
    colour = expr_mean),
    alpha = 1, na.rm = TRUE) +
  # scale_size_area(name = expression(-log[10](p)), max_size = 6,
  #                 breaks = c(1,2,3,5,10), limits = c(0, NA)) +
  scale_colour_gradientn(colors = c('white','#77C3BC','#4EA9C6','#2B84B2','#0E63A0',
                                    '#073E79','#0C203B'), 
                         name = "Expression (mean)") +
  facet_grid(A ~ ., scales = "free_y", space = "free_y") +
  labs(x = NULL, y = NULL,
       title    = "Macrophage subtypes (sender) → B cell (receiver)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid = element_blank(), 
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text.y  = element_text(face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    axis.text.x   = element_text(angle = 0, hjust = 0.5),
    legend.position = "right"
  ) # size 8*12



#### Macro to Endothelial cell L-R pairs分析 ---------------------------------------------------------
## load data
root <- "/storage/data/KAI/Pan_Inflammation/Homo/Mac2025/analysis/CellPhoneDB/result_new"
org_dirs <- fs::dir_ls(root, type = "directory")

incl_macros <- c(
  "IL1B+NLRP3+_Macro",
  "SPP1+CCL2+_Macro",
  "CD1C+_DC-like_Macro",
  "CCL13+_Complement-associated_Macro",
  "MMP3+CXCL8+_Macro"
)

target_immune <- "Endothelial cell"

organ_names <- basename(org_dirs) |>
  str_replace("_inflamed$","")

read_one_organ <- function(dir_path, organ_label) {
  # 精确匹配文件名（只取 overall means，不取 significant_means）
  f_pval <- list.files(dir_path, "^statistical_analysis_pvalues_.*\\.txt$", full.names = TRUE)
  f_mean <- list.files(dir_path, "^statistical_analysis_means_.*\\.txt$",  full.names = TRUE)
  
  # 若同类出现多个文件，取最新一个
  pick_latest <- function(paths) if (length(paths) <= 1) paths else paths[which.max(file.info(paths)$mtime)]
  if (length(f_pval) != 1) f_pval <- pick_latest(f_pval)
  if (length(f_mean) != 1) f_mean <- pick_latest(f_mean)
  
  # 若缺文件，返回空 tibble
  if (length(f_pval) == 0 || length(f_mean) == 0) {
    message("Skip (missing files): ", organ_label)
    return(tibble(
      organ=character(), interacting_pair=character(),
      gene_a=character(), gene_b=character(),
      A=character(), B=character(),
      pval=double(), log10p=double(), expr_mean=double(),
      LR=character()
    ))
  }
  
  pval <- read_tsv(f_pval, show_col_types = FALSE)
  mean <- read_tsv(f_mean, show_col_types = FALSE)
  
  # 找到 "clusterA|clusterB" 的列
  pair_cols_p <- grep("\\|", names(pval), value = TRUE)
  pair_cols_m <- grep("\\|", names(mean), value = TRUE)
  if (length(pair_cols_p) == 0 || length(pair_cols_m) == 0) {
    message("No pair columns in: ", organ_label)
    return(tibble(
      organ=character(), interacting_pair=character(),
      gene_a=character(), gene_b=character(),
      A=character(), B=character(),
      pval=double(), log10p=double(), expr_mean=double(),
      LR=character()
    ))
  }
  
  # ---- pvalues: 拉长 + 让键唯一（对重复键取最小 p） ----
  p_long <- pval %>%
    pivot_longer(all_of(pair_cols_p), names_to = "pair", values_to = "pval") %>%
    separate(pair, c("A","B"), sep="\\|", fill="right") %>%
    mutate(A = str_squish(A), B = str_squish(B)) %>%
    # 保留 gene_a/gene_b（若不存在也没关系）
    dplyr::select(interacting_pair,
                  gene_a = any_of("gene_a"),
                  gene_b = any_of("gene_b"),
                  A, B, pval) %>%
    group_by(interacting_pair, A, B) %>%
    summarise(
      pval   = suppressWarnings(min(as.numeric(pval), na.rm = TRUE)),
      gene_a = dplyr::first(gene_a),
      gene_b = dplyr::first(gene_b),
      .groups = "drop"
    )
  
  # ---- means: 拉长 + 让键唯一（对重复键取均值） ----
  m_long <- mean %>%
    pivot_longer(all_of(pair_cols_m), names_to = "pair", values_to = "expr_mean") %>%
    separate(pair, c("A","B"), sep="\\|", fill="right") %>%
    mutate(A = str_squish(A), B = str_squish(B)) %>%
    dplyr::select(interacting_pair,
                  gene_a_m = any_of("gene_a"),
                  gene_b_m = any_of("gene_b"),
                  A, B, expr_mean) %>%
    group_by(interacting_pair, A, B) %>%
    summarise(
      expr_mean = suppressWarnings(mean(as.numeric(expr_mean), na.rm = TRUE)),
      gene_a_m  = dplyr::first(gene_a_m),
      gene_b_m  = dplyr::first(gene_b_m),
      .groups = "drop"
    )
  
  # ---- 合并（此时每个键唯一，不再 many-to-many）----
  out <- p_long %>%
    left_join(m_long, by = c("interacting_pair","A","B")) %>%
    mutate(
      organ   = organ_label,
      # 兜底基因名：优先用 pvals 的，若无，用 means 的；再无，用 interacting_pair
      gene_a2 = dplyr::coalesce(gene_a, gene_a_m),
      gene_b2 = dplyr::coalesce(gene_b, gene_b_m),
      LR      = dplyr::if_else(!is.na(gene_a2) & !is.na(gene_b2),
                               paste0(gene_a2, " – ", gene_b2),
                               interacting_pair),
      pval    = as.numeric(pval),
      log10p  = -log10(pmax(pval, 1e-300))
    ) %>%
    dplyr::select(organ, interacting_pair, gene_a = gene_a2, gene_b = gene_b2,
                  A, B, pval, log10p, expr_mean, LR)
  
  out
}

df_all <- purrr::map2_dfr(org_dirs, organ_names, read_one_organ)

# 判定是否为 CellPhoneDB 里那类“复合体”命名
.is_complex <- function(x) {
  !is.na(x) & str_detect(x, "(?i)(^integrin_|_complex\\b)")
  # 说明：
  #  - ^integrin_    : 以 integrin_ 开头（CPDB 常见）
  #  - _complex\b    : 以 _complex 结尾（或单词边界），如 *_complex
}

# 过滤：去掉 gene_a/gene_b 缺失 & 去掉复合体命名
drop_complex_or_na <- function(df) {
  df %>%
    # 1) 基因名必须都有
    filter(!is.na(gene_a), !is.na(gene_b)) %>%
    # 2) 两侧都不能是“复合体”命名
    filter(!.is_complex(gene_a),
           !.is_complex(gene_b)) %>%
    # 3) 兜底：如果某些版本没有 gene_a/gene_b，或者列名不全，
    #    但 interacting_pair 里仍带有 complex，也顺手去掉
    filter(!str_detect(interacting_pair, "(?i)integrin_.*_complex|_complex\\b"))
}
df_all_clean <- df_all %>% drop_complex_or_na()

## 卡显著表达L-R pairs
df_mt <- df_all_clean %>%
  filter(A %in% incl_macros, B == target_immune) %>%  # 你的方向性筛选
  filter(pval < 0.05) 

TopN_per_organ <- 20
Final_cap <- 20

# 每个巨噬细胞A -> target_immune 取各器官top10
top_per_org <- df_mt %>%
  group_by(organ, A) %>%
  arrange(desc(expr_mean), .by_group = TRUE) %>%
  slice_head(n = TopN_per_organ) %>%
  ungroup()

# 各器官top10取并集
# 4) 对每个巨噬亚型 A：汇总“并集”并二次排序（频次→平均表达→统计强度）
rank_union_per_A <- top_per_org %>%
  group_by(A, organ, LR) %>%
  summarise(expr_mean = max(expr_mean, na.rm = TRUE),
            log10p = max(log10p, na.rm = TRUE),
            .groups = "drop") %>%
  group_by(A, LR) %>%
  summarise(
    freq_org   = n_distinct(organ),          # 跨器官出现频次
    mean_expr  = mean(expr_mean, na.rm = TRUE),
    sum_log10p = sum(log10p, na.rm = TRUE),  # 累积统计强度（可换成 mean 或 max）
    .groups = "drop"
  ) %>%
  arrange(A, desc(freq_org), desc(mean_expr), desc(sum_log10p)) %>%
  group_by(A) %>%
  mutate(rank = row_number()) %>%
  ungroup()

# 超过20选top20
keep_LR <- rank_union_per_A %>%
  group_by(A) %>%
  slice_head(n = Final_cap) %>%
  ungroup() %>%
  dplyr::select(A, LR)

# 制成绘图dataframe
plot_df <- df_mt %>%
  semi_join(keep_LR, by = c("A","LR")) %>%
  group_by(organ, A, LR) %>%
  summarise(
    expr_mean = mean(expr_mean, na.rm = TRUE),
    log10p    = max(log10p, na.rm = TRUE),   # 多条同 LR 取更强的统计
    .groups = "drop"
  ) %>%
  complete(organ = organ_names, nesting(A, LR), fill = list(expr_mean = 0, log10p = 0)) %>%
  mutate(
    organ = factor(organ, levels = organ_names),
    ## Row 排序：按 A 分面内，使用 4) 中的 rank
    LR = factor(LR, levels = rank_union_per_A %>%
                  semi_join(keep_LR, by = c("A","LR")) %>%
                  arrange(A, rank) %>%
                  pull(LR) %>% unique())
  )


eps <- 1e-8  # 判定“为 0”的容差

## 先找出“有任何非零表达”的器官
valid_orgs <- plot_df %>%
  dplyr::group_by(organ) %>%
  dplyr::summarise(any_pos = any(expr_mean > eps, na.rm = TRUE), .groups = "drop") %>%
  dplyr::filter(any_pos) %>%
  dplyr::pull(organ)

plot_df2 <- plot_df %>%
  dplyr::filter(organ %in% valid_orgs, expr_mean > eps) %>%
  dplyr::mutate(
    # 重新设定因子顺序，防止掉列后顺序错乱
    organ = factor(organ, levels = intersect(levels(plot_df$organ), valid_orgs))
  ) %>%
  droplevels()  # 清理掉未用到的水平

plot_df2$A <- factor(plot_df2$A, levels = incl_macros)

ggplot(plot_df2, aes(x = organ, y = LR)) +
  geom_point(aes(# size = log10p, 
    colour = expr_mean),
    alpha = 1, na.rm = TRUE) +
  # scale_size_area(name = expression(-log[10](p)), max_size = 6,
  #                 breaks = c(1,2,3,5,10), limits = c(0, NA)) +
  scale_colour_gradientn(colors = c('white','#77C3BC','#4EA9C6','#2B84B2','#0E63A0',
                                    '#073E79','#0C203B'), 
                         name = "Expression (mean)") +
  facet_grid(A ~ ., scales = "free_y", space = "free_y") +
  labs(x = NULL, y = NULL,
       title    = "Macrophage subtypes (sender) → Endothelial cell (receiver)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid = element_blank(), 
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text.y  = element_text(face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    axis.text.x   = element_text(angle = 0, hjust = 0.5),
    legend.position = "right"
  ) # size 8*12


#### Macro to Fibroblast L-R pairs分析 ---------------------------------------------------------
## load data
root <- "/storage/data/KAI/Pan_Inflammation/Homo/Mac2025/analysis/CellPhoneDB/result_new"
org_dirs <- fs::dir_ls(root, type = "directory")

incl_macros <- c(
  "IL1B+NLRP3+_Macro",
  "SPP1+CCL2+_Macro",
  "CD1C+_DC-like_Macro",
  "CCL13+_Complement-associated_Macro",
  "MMP3+CXCL8+_Macro"
)

target_immune <- "Fibroblast"

organ_names <- basename(org_dirs) |>
  str_replace("_inflamed$","")

read_one_organ <- function(dir_path, organ_label) {
  # 精确匹配文件名（只取 overall means，不取 significant_means）
  f_pval <- list.files(dir_path, "^statistical_analysis_pvalues_.*\\.txt$", full.names = TRUE)
  f_mean <- list.files(dir_path, "^statistical_analysis_means_.*\\.txt$",  full.names = TRUE)
  
  # 若同类出现多个文件，取最新一个
  pick_latest <- function(paths) if (length(paths) <= 1) paths else paths[which.max(file.info(paths)$mtime)]
  if (length(f_pval) != 1) f_pval <- pick_latest(f_pval)
  if (length(f_mean) != 1) f_mean <- pick_latest(f_mean)
  
  # 若缺文件，返回空 tibble
  if (length(f_pval) == 0 || length(f_mean) == 0) {
    message("Skip (missing files): ", organ_label)
    return(tibble(
      organ=character(), interacting_pair=character(),
      gene_a=character(), gene_b=character(),
      A=character(), B=character(),
      pval=double(), log10p=double(), expr_mean=double(),
      LR=character()
    ))
  }
  
  pval <- read_tsv(f_pval, show_col_types = FALSE)
  mean <- read_tsv(f_mean, show_col_types = FALSE)
  
  # 找到 "clusterA|clusterB" 的列
  pair_cols_p <- grep("\\|", names(pval), value = TRUE)
  pair_cols_m <- grep("\\|", names(mean), value = TRUE)
  if (length(pair_cols_p) == 0 || length(pair_cols_m) == 0) {
    message("No pair columns in: ", organ_label)
    return(tibble(
      organ=character(), interacting_pair=character(),
      gene_a=character(), gene_b=character(),
      A=character(), B=character(),
      pval=double(), log10p=double(), expr_mean=double(),
      LR=character()
    ))
  }
  
  # ---- pvalues: 拉长 + 让键唯一（对重复键取最小 p） ----
  p_long <- pval %>%
    pivot_longer(all_of(pair_cols_p), names_to = "pair", values_to = "pval") %>%
    separate(pair, c("A","B"), sep="\\|", fill="right") %>%
    mutate(A = str_squish(A), B = str_squish(B)) %>%
    # 保留 gene_a/gene_b（若不存在也没关系）
    dplyr::select(interacting_pair,
                  gene_a = any_of("gene_a"),
                  gene_b = any_of("gene_b"),
                  A, B, pval) %>%
    group_by(interacting_pair, A, B) %>%
    summarise(
      pval   = suppressWarnings(min(as.numeric(pval), na.rm = TRUE)),
      gene_a = dplyr::first(gene_a),
      gene_b = dplyr::first(gene_b),
      .groups = "drop"
    )
  
  # ---- means: 拉长 + 让键唯一（对重复键取均值） ----
  m_long <- mean %>%
    pivot_longer(all_of(pair_cols_m), names_to = "pair", values_to = "expr_mean") %>%
    separate(pair, c("A","B"), sep="\\|", fill="right") %>%
    mutate(A = str_squish(A), B = str_squish(B)) %>%
    dplyr::select(interacting_pair,
                  gene_a_m = any_of("gene_a"),
                  gene_b_m = any_of("gene_b"),
                  A, B, expr_mean) %>%
    group_by(interacting_pair, A, B) %>%
    summarise(
      expr_mean = suppressWarnings(mean(as.numeric(expr_mean), na.rm = TRUE)),
      gene_a_m  = dplyr::first(gene_a_m),
      gene_b_m  = dplyr::first(gene_b_m),
      .groups = "drop"
    )
  
  # ---- 合并（此时每个键唯一，不再 many-to-many）----
  out <- p_long %>%
    left_join(m_long, by = c("interacting_pair","A","B")) %>%
    mutate(
      organ   = organ_label,
      # 兜底基因名：优先用 pvals 的，若无，用 means 的；再无，用 interacting_pair
      gene_a2 = dplyr::coalesce(gene_a, gene_a_m),
      gene_b2 = dplyr::coalesce(gene_b, gene_b_m),
      LR      = dplyr::if_else(!is.na(gene_a2) & !is.na(gene_b2),
                               paste0(gene_a2, " – ", gene_b2),
                               interacting_pair),
      pval    = as.numeric(pval),
      log10p  = -log10(pmax(pval, 1e-300))
    ) %>%
    dplyr::select(organ, interacting_pair, gene_a = gene_a2, gene_b = gene_b2,
                  A, B, pval, log10p, expr_mean, LR)
  
  out
}

df_all <- purrr::map2_dfr(org_dirs, organ_names, read_one_organ)

# 判定是否为 CellPhoneDB 里那类“复合体”命名
.is_complex <- function(x) {
  !is.na(x) & str_detect(x, "(?i)(^integrin_|_complex\\b)")
  # 说明：
  #  - ^integrin_    : 以 integrin_ 开头（CPDB 常见）
  #  - _complex\b    : 以 _complex 结尾（或单词边界），如 *_complex
}

# 过滤：去掉 gene_a/gene_b 缺失 & 去掉复合体命名
drop_complex_or_na <- function(df) {
  df %>%
    # 1) 基因名必须都有
    filter(!is.na(gene_a), !is.na(gene_b)) %>%
    # 2) 两侧都不能是“复合体”命名
    filter(!.is_complex(gene_a),
           !.is_complex(gene_b)) %>%
    # 3) 兜底：如果某些版本没有 gene_a/gene_b，或者列名不全，
    #    但 interacting_pair 里仍带有 complex，也顺手去掉
    filter(!str_detect(interacting_pair, "(?i)integrin_.*_complex|_complex\\b"))
}
df_all_clean <- df_all %>% drop_complex_or_na()

## 卡显著表达L-R pairs
df_mt <- df_all_clean %>%
  filter(A %in% incl_macros, B == target_immune) %>%  # 你的方向性筛选
  filter(pval < 0.05) 

TopN_per_organ <- 10
Final_cap <- 20

# 每个巨噬细胞A -> target_immune 取各器官top10
top_per_org <- df_mt %>%
  group_by(organ, A) %>%
  arrange(desc(expr_mean), .by_group = TRUE) %>%
  slice_head(n = TopN_per_organ) %>%
  ungroup()

# 各器官top10取并集
# 4) 对每个巨噬亚型 A：汇总“并集”并二次排序（频次→平均表达→统计强度）
rank_union_per_A <- top_per_org %>%
  group_by(A, organ, LR) %>%
  summarise(expr_mean = max(expr_mean, na.rm = TRUE),
            log10p = max(log10p, na.rm = TRUE),
            .groups = "drop") %>%
  group_by(A, LR) %>%
  summarise(
    freq_org   = n_distinct(organ),          # 跨器官出现频次
    mean_expr  = mean(expr_mean, na.rm = TRUE),
    sum_log10p = sum(log10p, na.rm = TRUE),  # 累积统计强度（可换成 mean 或 max）
    .groups = "drop"
  ) %>%
  arrange(A, desc(freq_org), desc(mean_expr), desc(sum_log10p)) %>%
  group_by(A) %>%
  mutate(rank = row_number()) %>%
  ungroup()

# 超过20选top20
keep_LR <- rank_union_per_A %>%
  group_by(A) %>%
  slice_head(n = Final_cap) %>%
  ungroup() %>%
  dplyr::select(A, LR)

# 制成绘图dataframe
plot_df <- df_mt %>%
  semi_join(keep_LR, by = c("A","LR")) %>%
  group_by(organ, A, LR) %>%
  summarise(
    expr_mean = mean(expr_mean, na.rm = TRUE),
    log10p    = max(log10p, na.rm = TRUE),   # 多条同 LR 取更强的统计
    .groups = "drop"
  ) %>%
  complete(organ = organ_names, nesting(A, LR), fill = list(expr_mean = 0, log10p = 0)) %>%
  mutate(
    organ = factor(organ, levels = organ_names),
    ## Row 排序：按 A 分面内，使用 4) 中的 rank
    LR = factor(LR, levels = rank_union_per_A %>%
                  semi_join(keep_LR, by = c("A","LR")) %>%
                  arrange(A, rank) %>%
                  pull(LR) %>% unique())
  )


eps <- 1e-8  # 判定“为 0”的容差

## 先找出“有任何非零表达”的器官
valid_orgs <- plot_df %>%
  dplyr::group_by(organ) %>%
  dplyr::summarise(any_pos = any(expr_mean > eps, na.rm = TRUE), .groups = "drop") %>%
  dplyr::filter(any_pos) %>%
  dplyr::pull(organ)

plot_df2 <- plot_df %>%
  dplyr::filter(organ %in% valid_orgs, expr_mean > eps) %>%
  dplyr::mutate(
    # 重新设定因子顺序，防止掉列后顺序错乱
    organ = factor(organ, levels = intersect(levels(plot_df$organ), valid_orgs))
  ) %>%
  droplevels()  # 清理掉未用到的水平

plot_df2$A <- factor(plot_df2$A, levels = incl_macros)

ggplot(plot_df2, aes(x = organ, y = LR)) +
  geom_point(aes(# size = log10p, 
    colour = expr_mean),
    alpha = 1, na.rm = TRUE) +
  # scale_size_area(name = expression(-log[10](p)), max_size = 6,
  #                 breaks = c(1,2,3,5,10), limits = c(0, NA)) +
  scale_colour_gradientn(colors = c('white','#77C3BC','#4EA9C6','#2B84B2','#0E63A0',
                                    '#073E79','#0C203B'), 
                         name = "Expression (mean)") +
  facet_grid(A ~ ., scales = "free_y", space = "free_y") +
  labs(x = NULL, y = NULL,
       title    = "Macrophage subtypes (sender) → Fibroblast (receiver)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid = element_blank(), 
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text.y  = element_text(face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    axis.text.x   = element_text(angle = 0, hjust = 0.5),
    legend.position = "right"
  ) # size 8*12



