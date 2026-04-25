### Package ------------------------------
library(ggplot2)
library(readxl)
library(patchwork)
library(tidyr)
library(dplyr)

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
library(boot)

library(monocle)
library(SingleCellExperiment)
library(ggsci)
library(igraph)

### global -----------------------------------
level3_color = c(
  "Monocyte"="#8184e2",
  "mm_Mac_01"="#FFD372",
  "mm_Mac_02"="#FF98A4",
  "mm_Mac_03"="#82CCDD",
  "mm_Mac_04"="#D7BCE7",
  "mm_Mac_05"="#C9C9C9",
  "mm_Mac_06"="#FFA500",
  "mm_Mac_07"="#34C0B8",
  "mm_Mac_08"="#E76F51",
  "mm_Mac_09"="#406196",
  "mm_Mac_10"="#82447F"
)

### load data --------------------------------
mm_mono2mac <- readRDS("/storage/data/KAI/Pan_Inflammation/Mouse/Pneumonia/LC/analysis/monocle2/data/plmouse_mono2mac.RDS")


### 创建CellDataSet对象 ----------------------
ct <- mm_mono2mac@assays$RNA$counts # count 矩阵 
gene_ann <- data.frame(
  gene_short_name = row.names(ct), 
  row.names = row.names(ct)
)
fd <- new("AnnotatedDataFrame",
          data=gene_ann)       # 基因注释

pd <- new("AnnotatedDataFrame",
          data=mm_mono2mac@meta.data) # 临床信息
sc_cds <- newCellDataSet(
  ct, 
  phenoData = pd,
  featureData =fd,
  expressionFamily = negbinomial.size(),
  lowerDetectionLimit=0.5
)

### 用于估算每个细胞的大小因子,目的是归一化基因表达数据以纠正测序深度的差异
sc_cds <- estimateSizeFactors(sc_cds)

### 估算基因表达的离散度，帮助识别高变异的基因，用于差异基因分析和排序
sc_cds <- estimateDispersions(sc_cds)


disp_table <- dispersionTable(sc_cds)

ordering_genes_temp <- subset(disp_table, mean_expression >= 0.6 & dispersion_empirical >= 2*dispersion_fit) 
ordering_genes<-ordering_genes_temp$gene_id
# 设置用于排序的基因
sc_cds <- setOrderingFilter(sc_cds, ordering_genes)

plot_ordering_genes(sc_cds)

# 对细胞进行降维和去批次，常用于将高维的基因表达数据转换为低维（如 2D）空间，便于可视化和后续分析
sc_cds <- reduceDimension(sc_cds, max_components = 2, method = 'DDRTree',norm_method = "log", num_dim = 30)

# 根据降维后的数据和排序基因，进行细胞的拟时序排序，推断细胞在轨迹中的位置
sc_cds <- orderCells(sc_cds)



plot_cell_trajectory(sc_cds, color_by = "cell_type_level3", cell_size = 0.5) + 
  scale_color_manual(values=level3_color) + 
  theme(plot.title = element_blank(),
        legend.position = "none",
        ) + 
  # theme_void() +
  guides(color = "none") 



plot_cell_trajectory(sc_cds, color_by = "Pseudotime",cell_size = 0.5)+
  scale_color_gradientn(
    colours = c("#eeecdf","#eaebea", "#becdd2", "#6f9ad1", "#44679f", "#3f4f71")
  ) +
  # theme_void() +
  theme(legend.position = "none")
  
plot_cell_trajectory(sc_cds, color_by = "Sample")

plot_cell_trajectory(sc_cds, color_by = "State")

plot_cell_trajectory(sc_cds, color_by = 'cell_type_level3', cell_size = 0.5,cell_link_size = 0.5) +facet_wrap(~cell_type_level3, nrow = 3, scales = "free")+scale_color_manual(values=level3_color)+ggtitle('cell_type_level3') + theme(plot.title = element_text(hjust = 0.5),legend.position = "right")

ggplot(as.data.frame(pData(sc_cds)), aes(Pseudotime, fill = cell_type_level3)) +geom_density() +facet_wrap(~cell_type_level3) +theme_bw() +RotatedAxis() +theme(  strip.text = element_blank(),  strip.background = element_rect(color = "white", fill = "white"),  panel.grid = element_blank()) +scale_fill_manual(values = level3_color)


plot_complex_cell_trajectory(sc_cds,color_by = "cell_type_level3") +  
  scale_color_manual(values=level3_color) +
  theme(legend.title = element_blank())


pseudo_sign <- c('Il1b','Spp1')
plot_genes_in_pseudotime(sc_cds[pseudo_sign,],color_by = "cell_type_level3", cell_size = 0.5)+
  scale_color_manual(values=level3_color)











