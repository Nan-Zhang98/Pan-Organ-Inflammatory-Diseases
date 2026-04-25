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

library(scRNAtoolVis)
library(clusterProfiler)
library(org.Hs.eg.db) 
library(AnnotationDbi)

### load data -------------------------------
Mac <- readRDS("/storage/data/KAI/Pan_Inflammation/Homo/Mac2025/scVI/new/pan_Mac_filtered.RDS") 
DEG_Mac <- readRDS("/storage/data/KAI/Pan_Inflammation/Homo/Mac2025/analysis/DEG/data/DEG_Mac.RDS")


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



############################################## Findallmarkers DEG ############################################################
# Idents(Mac) <- Mac$cell_type_level3

Mac2025 <- NormalizeData(Mac2025, normalization.method = "LogNormalize", scale.factor = 10000)


### 细胞类型间差异基因计算
DEG_Mac <- FindAllMarkers(Mac2025, log2FC.threshold = 0, test.use = "wilcox",
                          min.pct = 0.25, min.diff.pct = 0, only.pos = F,
                          assay = "RNA")

saveRDS(DEG_Mac,"/storage/data/KAI/Pan_Inflammation/Homo/Mac2025/analysis/DEG/data/DEG_Mac.RDS")
# 选取top基因展示

### 可视化

DEG_moMac <- DEG_Mac[which(DEG_Mac$cluster %in% c('CCL13+_Complement-associated_Macro','CD1C+_DC-like_Macro','IL1B+NLRP3+_Macro','MMP3+CXCL8+_Macro','SPP1+CCL2+_Macro')),]
DEG_moMac$cluster <- droplevels(DEG_moMac$cluster)

## 环状火山图
jjVolcano(diffData = DEG_moMac,
          log2FC.cutoff = 0.5,
          #myMarkers = c('SPP1','IL1B','NLRP3','CCL2','APOE','FABP2','FOLR2'),
          topGeneN = 5,
          adjustP.cutoff = 0.05,
          tile.col = Mac_color,
          polar = T,
          fontface = 'italic',
          pSize = 1,
          base_size = 20,
          expand = c(-1.4,1.4),
          celltypeSize = 2,
          order.by = 'avg_log2FC',
          legend.position	= c(0.95,0.9)
          ) 


## 多组火山图 (p值=0基因过多分布过于集中，可视化不美观)
top10 <- DEG_moMac %>% 
  group_by(cluster) %>%
  top_n(10, abs(avg_log2FC)) %>%
  ungroup() %>% 
  as.data.frame()

DEG_moMac$log10fdr <- -log10(DEG_moMac$p_val_adj)
DEG_moMac$log10fdr[DEG_moMac$log10fdr=="Inf"] <- max(DEG_moMac$log10fdr[DEG_moMac$log10fdr!=Inf])

DEG_moMac$log10p <- -log10(DEG_moMac$p_val)
DEG_moMac$log10p[DEG_moMac$log10p=="Inf"] <- max(DEG_moMac$log10p[DEG_moMac$log10p!=Inf])


p <- ggplot() +
  geom_point(data = DEG_moMac, aes(x = avg_log2FC, y = log10fdr),size = 0.8, color = 'grey') +
  coord_flip() + # 坐标轴翻转
  facet_grid(. ~ cluster,scales = "free") + # 一行多列;
  geom_point(data = top10, aes(x = avg_log2FC, y = log10fdr,color = cluster)) + # 添加top点颜色
  geom_vline(xintercept = c(-0.5, 0.5), size = 0.5, color = "grey50", lty = 'dashed')+ #添加阈值线
  scale_color_manual(values = Mac_color) + #更改配色
  xlab(label = "avg_log2FC") + 
  ylab(label = "") + 
  theme_bw()+
  theme( legend.position = 'none', #去掉图例
         panel.grid = element_blank(), #去掉背景网格
         axis.text = element_text(size = 10), #坐标轴标签大小
         axis.text.x = element_text(angle = 45, vjust = 0.8), #x轴标签旋转
         strip.text.x = element_text(size = 10, face = 'bold') #加粗分面标题
  )


## 多组火山图2
log2FC.cutoff = 0.5
pvalue.cutoff = 0.05
adjustP.cutoff = 0.05
# tile.col = jjAnno::useMyCol("paired",n = 9)

xnum <- seq_along(levels(DEG_moMac$cluster))

tile.df <- data.frame(
  cluster = levels(DEG_moMac$cluster),
  xmin = xnum - 0.5,
  xmax = xnum + 0.5,
  ymin = -8,
  ymax = -10
)



p1 <- ggplot(DEG_moMac, aes(x = cluster, y = avg_log2FC, fill = cluster)) +
  geom_jitter(aes(color = cluster,alpha = 0.9)) +
  scale_fill_manual(values = Mac_color) +
  scale_color_manual(values = Mac_color) +
  scale_y_continuous(
    limits = c(-10,8.5),   # 上下限
    breaks = c(-8,-6,-4,-2, 0, 2, 4, 6, 8),   # 刻度位置
    labels = c( "-8","-6","-4", "-2", "0", "2", "4", "6","8")  # 刻度标签
  ) +
  geom_hline(yintercept = c(-0.5,0.5),linetype = "dashed", color = "red") +
  geom_rect(data = tile.df, 
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,fill = cluster),
            alpha = 0.8, inherit.aes = FALSE) +
  guides(fill="none",alpha = "none",color="none") +
  theme(
    axis.text.x = element_text(angle = 0,size = 12,hjust = 0.5,colour = "black"),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 10),
    axis.title.y.left = element_text(size = 12,colour = "black"),
    axis.line = element_line(colour = "black"),
    panel.background = element_blank()) 


diff.marker <- DEG_moMac %>%
  dplyr::filter(
    abs(avg_log2FC) >= log2FC.cutoff & p_val < pvalue.cutoff
  ) %>%
  dplyr::mutate(
    type = ifelse(
      avg_log2FC >= log2FC.cutoff,
      "sigUp",
      "sigDown"
    )
  ) %>%
  dplyr::mutate(
    type2 = ifelse(
      p_val_adj < adjustP.cutoff,
      paste("adjust Pvalue < ", adjustP.cutoff, sep = ""),
      paste("adjust Pvalue >= ", adjustP.cutoff, sep = "")
    )
  )


top.marker.tmp <- diff.marker %>% dplyr::group_by(cluster)
top.marker.max <- top.marker.tmp %>% dplyr::slice_max(n = 5, order_by = avg_log2FC)

top.marker.min <- top.marker.tmp %>% dplyr::slice_min(n = 5, order_by = avg_log2FC)
top.marker <- rbind(top.marker.max, top.marker.min)
highlight.genes <- DEG_moMac %>% 
  filter(gene %in% c("IL1B", "SPP1", "NLRP3", "CCL2","MMP3","CXCL8","CCL13","C1QA","C1QB","C1QC"))

p2 <- p1 +
  geom_text_repel(
    data = top.marker,
    aes(x = cluster, y = avg_log2FC, label = gene,color = cluster),
    inherit.aes = FALSE,
    size=4,
    fontface = "italic",
    force = 3,
    nudge_y = 0.3,
    segment.size = 0.3,  # 线段粗细
    segment.color = "grey40",  # 线段颜色
    segment.alpha = 1  # 线段透明度
  ) +
  geom_text_repel(
    data = highlight.genes,
    aes(x = cluster, y = avg_log2FC, label = gene,color = cluster),
    inherit.aes = FALSE,
    # color = "black",           # 可以让它颜色更明显
    fontface = "italic", 
    force = 3,
    size = 4,
    nudge_y = 0.3,
    segment.size = 0.3,  # 线段粗细
    segment.color = "grey40",  # 线段颜色
    segment.alpha = 1  # 线段透明度
  )


