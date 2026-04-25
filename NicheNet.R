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

library(nichenetr)


### load data -------------------------------
moMac <- readRDS("/storage/data/KAI/Pan_Inflammation/Homo/Mac2025/scVI/new/mono_mac.RDS")

# 进行NicheNet分析前进行标准化。
moMac <- NormalizeData(moMac, normalization.method = "LogNormalize", scale.factor = 10000)



# 根据链接下载到本地
lr_network <- readRDS("/storage/data/KAI/Pan_Inflammation/Homo/Infla_R/Macrophage/re-annotation/results/NicheNet/input/lr_network_human_21122021.rds") # ligand和target的对应关系

ligand_target_matrix <- readRDS("/storage/data/KAI/Pan_Inflammation/Homo/Infla_R/Macrophage/re-annotation/results/NicheNet/input/ligand_target_matrix_nsga2r_final.rds") # ligand到target矩阵，行是靶基因，列是ligand，数值为调控的可能性

weighted_networks <- readRDS("/storage/data/KAI/Pan_Inflammation/Homo/Infla_R/Macrophage/re-annotation/results/NicheNet/input/weighted_networks_nsga2r_final.rds") # 在一个由"LR pair和下游信号分子"构成的网络中LR pair/其他pair的互作潜力


# 检查数据
lr_network <- lr_network %>% distinct(from, to)
ligand_target_matrix[1:5,1:5]

weighted_networks_lr <- weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from, to), by = c("from", "to")) #与lr_network的from进行合并
head(weighted_networks$lr_sig)
head(weighted_networks$gr)


################################################ NicheNet IL1B+NLRP3+_Macro ####################################################
# 1.定义潜在ligand集合
receiver <- "IL1B+NLRP3+_Macro"
Idents(moMac) <- moMac$cell_type_level3
expressed_genes_receiver <- get_expressed_genes(receiver, moMac, pct = 0.10)

all_receptors <- unique(lr_network$to)  
expressed_receptors <- intersect(all_receptors, expressed_genes_receiver) # 244个表达受体

potential_ligands <- lr_network %>% filter(to %in% expressed_receptors) %>% pull(from) %>% unique() # 729个潜在配体

background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)] # 6466

# 2.定义感兴趣基因集（炎症样本"IL1B+NLRP3+_Macro"差异基因）
seurat_obj_receiver <- subset(moMac, idents = receiver)
seurat_obj_receiver <- SetIdent(seurat_obj_receiver, value = seurat_obj_receiver$Inflammation)
condition_oi <- "Inflamed"
condition_reference <- "Normal"
DE_table_receiver <- FindMarkers(object = seurat_obj_receiver, ident.1 = condition_oi, ident.2 = condition_reference, min.pct = 0.10) %>% rownames_to_column("gene") # 5639
geneset_oi <- DE_table_receiver %>% filter(p_val_adj <= 0.01 & avg_log2FC >= 0.25) %>% pull(gene)
geneset_oi <- geneset_oi %>% .[. %in% rownames(ligand_target_matrix)] # 6264 



# 3.NicheNet ligand activity
# 一个ligand的靶基因落在"兴趣基因集"里面越多，认为这个ligand的活性越高
ligand_activities <- predict_ligand_activities(geneset = geneset_oi, 
                                               background_expressed_genes = background_expressed_genes,
                                               ligand_target_matrix = ligand_target_matrix,
                                               potential_ligands = potential_ligands
) # 729*5 pearson衡量ligand activity
# pearson： ligand预测的target基因和感兴趣基因的相关性
ligand_activities <- ligand_activities %>% 
  arrange(desc(pearson)) %>%   # 按照 pearson 降序排列
  mutate(rank = rank(-pearson))  # 对 pearson 列降序排列后计算排名

best_upstream_ligands <- ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique() # top20 ligands

# DotPlot(Macro, features = best_upstream_ligands %>% rev(), cols = "RdYlBu") + RotatedAxis()


top20_ligands <- ligand_activities %>% top_n(20, pearson) %>% arrange(pearson)
top20_ligands_matrix <- as.numeric(top20_ligands$pearson) %>% as.matrix()
rownames(top20_ligands_matrix) <- top20_ligands$test_ligand
colnames(top20_ligands_matrix) <- "pearson"
make_heatmap_ggplot(top20_ligands_matrix, "Prioritized ligands", "",
                    color = "orange", legend_title = "Ligand activity") +
  coord_fixed(ratio = 1) +
  scale_fill_gradient2(low = "whitesmoke",  high = "orange")  





# 4.根据best_upstream_ligands，推测receptors和top-predicted target genes
# active target genes inference
active_ligand_target_links_df <- best_upstream_ligands %>% 
  lapply(get_weighted_ligand_target_links,
         geneset = geneset_oi,
         ligand_target_matrix = ligand_target_matrix,
         n = 100) %>%
  bind_rows() %>% drop_na()
nrow(active_ligand_target_links_df) # 1115, 带权重

# 将非常小的权重作为0处理（ligand_target_matrix所有值的0.33分位数）
active_ligand_target_links <- prepare_ligand_target_visualization(
  ligand_target_df = active_ligand_target_links_df,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = 0.33)  #461*16 16ligands


order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links))

# IL1B+ Macro top5 regulon
top5_regulon <- c('CHD1','ELF2','NFIL3','RCOR1','STAT3')

intersect(order_targets, top5_regulon)

vis_ligand_target <- t(active_ligand_target_links[order_targets,order_ligands])

make_heatmap_ggplot(vis_ligand_target, "Prioritized ligands", "Predicted target genes",
                    color = "purple", legend_title = "Regulatory potential") +
  scale_fill_gradient2(low = "whitesmoke",  high = "purple")




intersect(DE_table_receiver$gene, colnames(vis_ligand_target)) #461

length(colnames(vis_ligand_target)) # 461
length(DE_table_receiver$gene) # 6900





# 5.识别receptors of top-ranked ligands
ligand_receptor_links_df <- get_weighted_ligand_receptor_links(
  best_upstream_ligands, expressed_receptors,
  lr_network, weighted_networks$lr_sig
) 

vis_ligand_receptor_network <- prepare_ligand_receptor_visualization(
  ligand_receptor_links_df,
  best_upstream_ligands,
  order_hclust = "both") 

(make_heatmap_ggplot(t(vis_ligand_receptor_network[, order_ligands]), 
                     y_name = "Ligands", x_name = "Receptors",  
                     color = "mediumvioletred", legend_title = "Prior interaction potential"))




lr_network_top <- lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_genes_receiver) %>% distinct(from, to)
best_upstream_receptors <- lr_network_top %>% pull(to) %>% unique()

lr_network_top_df_large <- weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)
lr_network_top_df <- lr_network_top_df_large %>% spread("from", "weight", fill = 0)
lr_network_top_matrix <- lr_network_top_df %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)

dist_receptors <- dist(lr_network_top_matrix, method = "binary")
hcluster_receptors <- hclust(dist_receptors, method = "ward.D2")
order_receptors <- hcluster_receptors$labels[hcluster_receptors$order]

order_ligands <- best_upstream_ligands %>% rev()
order_receptors <- order_receptors %>% intersect(colnames(lr_network_top_matrix)) # 0


### nichenet_seuratobj_aggregate封装函数

IL1B.nichenet <- nichenet_seuratobj_aggregate(
  seurat_obj = moMac, # Seurat对象，其active.ident需设置为细胞类型
  expression_pct = 0.10,  # 界定细胞类型是否表达配/受体的比例阈值，默认为0.1
  # organism = "human", # 交代物种信息，默认为人类 c("human","mouse")
  #Group
  condition_colname = "Inflammation",   # 交代分组的meta名
  condition_oi = "Inflamed", condition_reference = "Normal", # 交代实验组与对照组名 
  # receiver
  receiver = "IL1B+NLRP3+_Macro",  # 交代receiver细胞类型
  geneset = "DE", # 判断特定基因集的方法，默认使用全部差异基因(oi/ref)c("DE","up","down")
  lfc_cutoff = 0.25, # 判断差异基因的阈值
  # sender
  top_n_targets = 200,  #每个ligand最多考虑200个target gene
  top_n_ligands = 20,   #给出最有可能的20个上游ligand
  cutoff_visualization = 0.33,  #设置可视化ligand-target scores的阈值
  # refer data
  ligand_target_matrix = ligand_target_matrix, 
  lr_network = lr_network, 
  weighted_networks = weighted_networks, 
)

saveRDS(IL1B.nichenet,"/storage/data/KAI/Pan_Inflammation/Homo/Mac2025/analysis/NicheNet/data/IL1B_nichenet.RDS")

IL1B.nichenet <- readRDS("/storage/data/KAI/Pan_Inflammation/Homo/Mac2025/analysis/NicheNet/data/IL1B_nichenet.RDS")

IL1B.nichenet$ligand_activities[1:20,]
IL1B.nichenet$top_ligands

IL1B.nichenet$ligand_expression_dotplot
IL1B.nichenet$ligand_differential_expression_heatmap
IL1B.nichenet$ligand_activity_target_heatmap

IL1B.nichenet$ligand_receptor_heatmap$data


heatmap_data_sorted <- IL1B.nichenet$ligand_receptor_heatmap$data %>%
  mutate(y = factor(y, levels = IL1B.nichenet$top_ligands %>% rev())) %>%  # 设置 y 列的因子顺序
  group_by(y) %>%        # 按照配体（y）分组
  arrange(y, desc(score))  # 按照 score 降序排序，并保持 y 列的顺序

heatmap_data_sorted <- IL1B.nichenet$ligand_receptor_heatmap$data %>%
  mutate(y = factor(y, levels = IL1B.nichenet$top_ligands %>% rev())) %>%  # 设置 'y' 列的因子顺序
  group_by(y) %>%        # 按照配体（y）分组
  arrange(y, desc(score)) %>%  # 按照 score 降序排序，并保持 y 列的顺序
  ungroup() %>%         # 解除分组
  mutate(x = factor(x, levels = unique(x) %>% rev()))  # 设置 'x' 列按顺序排列


ggplot(heatmap_data_sorted, aes(x = x, y = y, fill = score)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "mediumvioletred") +
  theme_minimal() +
  labs(x = "Receptors", y = "Ligands", fill = "Prior interaction potential") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))








write.table(as.matrix(moMac@assays$RNA$data),file = "/storage/data/KAI/Pan_Inflammation/Homo/Mac2025/analysis/CytoSig/data/moMac_log.txt",sep = '\t', quote = F)



meta_data <- cbind(rownames(moMac@meta.data),moMac@meta.data[,'cell_type_level3', drop=F],moMac@meta.data[,'Inflammation', drop=F])
meta_data <- as.matrix(meta_data)
colnames(meta_data) <- c('Cell','cell_type_level3','Inflammation')
is.na(meta_data) %>% sum() # 没有NA值

write.table(meta_data, "/storage/data/KAI/Pan_Inflammation/Homo/Mac2025/analysis/CytoSig/data/moMac_meta.txt",sep = '\t', quote = F)













################################################ NicheNet SPP1+CCL2+_Macro ####################################################
SPP1.nichenet <- nichenet_seuratobj_aggregate(
                    seurat_obj = moMac, # Seurat对象，其active.ident需设置为细胞类型
                    expression_pct = 0.10,  # 界定细胞类型是否表达配/受体的比例阈值，默认为0.1
                    # organism = "human", # 交代物种信息，默认为人类 c("human","mouse")
                    # Group
                    condition_colname = "Inflammation",   # 交代分组的meta名
                    condition_oi = "Inflamed", condition_reference = "Normal", # 交代实验组与对照组名 
                    # receiver
                    receiver = "SPP1+CCL2+_Macro",  # 交代receiver细胞类型
                    geneset = "DE", # 判断特定基因集的方法，默认使用全部差异基因(oi/ref)c("DE","up","down")
                    lfc_cutoff = 0.25, # 判断差异基因的阈值
                    # sender
                    top_n_targets = 200,  #每个ligand最多考虑200个target gene
                    top_n_ligands = 20,   #给出最有可能的20个上游ligand
                    cutoff_visualization = 0.33,  #设置可视化ligand-target scores的阈值
                    # refer data
                    ligand_target_matrix = ligand_target_matrix, 
                    lr_network = lr_network, 
                    weighted_networks = weighted_networks, 
                  )
saveRDS(SPP1.nichenet,"/storage/data/KAI/Pan_Inflammation/Homo/Mac2025/analysis/NicheNet/data/SPP1_nichenet.RDS")

SPP1.nichenet <- readRDS("/storage/data/KAI/Pan_Inflammation/Homo/Mac2025/analysis/NicheNet/data/SPP1_nichenet.RDS")

SPP1.nichenet$ligand_activities

SPP1.nichenet$top_ligands






