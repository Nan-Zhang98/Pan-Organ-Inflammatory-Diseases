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
library(data.table)

library(diffusionMap)
library(pcaMethods)
library(destiny)
library(scatterplot3d)
library(rgl)
library(ggtern)
library(viridis)
library(AUCell)

library(slingshot)
library(tradeSeq)
library(BiocParallel)

library(monocle3)



### load data -------------------------------
## 同monocle2抽样输入数据
Mac_10 <- readRDS("/storage/data/KAI/Pan_Inflammation/Homo/Mac2025/analysis/diffusion_map/data/Mac_10.RDS")

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



######################################## Diffusion map all moMac downsampled 10% #################################################

scvi <- Embeddings(Mac_10, "scVI")
dmm <- DiffusionMap(scvi)
dpt = DPT(dmm)

fwrite(dmm@transitions %>% as.matrix(),file='/storage/data/KAI/Pan_Inflammation/Homo/Mac2025/analysis/diffusion_map/data/DM_Similarity.txt',sep = '\t',quote = F,row.names = F)

### useful visual ------------------------------

## 3D
plot(dmm)

## 2D
common_cells <- rownames(dmm@eigenvectors)  # 获取 DiffusionMap 仍然保留的细胞
filtered_level3 <- Mac_10$cell_type_level3[common_cells]  # 重新匹配细胞类型

dm_eigen <- data.frame(
  DC1 = dmm@eigenvectors[, 1],
  DC2 = dmm@eigenvectors[, 2],
  cell_type = filtered_level3  # 细胞类型
)

ggplot(dm_eigen, aes(x = DC1, y = DC2, fill = cell_type)) +
  geom_point(shape = 21, size = 3, stroke = 0.5) +  # 使用 shape=21 以便填充颜色
  theme_bw() +
  scale_fill_manual(values = Mac_color) +  # 使用预定义颜色
  xlab("DC 1") +
  ylab("DC 2") +
  ggtitle("Diffusion Map of Macrophage Subtypes") +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    plot.title = element_text(size = 18, hjust = 0.5)
  )
dm_eigen <- data.frame(
  DC1 = dmm@eigenvectors[, 1],  # Diffusion component 1
  DC2 = dmm@eigenvectors[, 2],  # Diffusion component 2
  DC3 = dmm@eigenvectors[, 3],  # Diffusion component 2
  cell_type = filtered_level3,  # 细胞类型
  pseudotime = dpt$dpt          # DPT 计算的伪时间
  
)
write.table(dm_eigen,"/storage/data/KAI/Pan_Inflammation/Homo/Mac2025/analysis/diffusion_map/data/dm_eigen.txt",sep = '\t',row.names = T,quote = F)

write.table(dmm@eigenvalues,"/storage/data/KAI/Pan_Inflammation/Homo/Mac2025/analysis/diffusion_map/data/DM_EigenValues.txt",sep = '\t',row.names = T,quote = F)
write.table(dmm@eigenvectors,"/storage/data/KAI/Pan_Inflammation/Homo/Mac2025/analysis/diffusion_map/data/DM_EigenVectors.txt",sep = '\t',row.names = T,quote = F)


ggplot(dm_eigen, aes(x = DC1, y = DC2, color = pseudotime)) +
  geom_point(size = 1) + 
  scale_color_viridis_c(option = "plasma") +  # 颜色渐变方案
  theme_minimal() +
  labs(title = "Diffusion Map with Pseudotime", 
       x = "DC 1", 
       y = "DC 2", 
       color = "Pseudotime") +
  theme(
    plot.title = element_text(size=16, face="bold", hjust = 0.5),
    axis.text = element_text(size=14),
    axis.title = element_text(size=14, face="bold"),
    legend.text = element_text(size=12)
  )

### 导入DC坐标和dpt返回seurat对象Mac_10 -----------------------------------------
## Diffusion map嵌入回seurat对象
Mac_10[["diffusion_map"]] <- CreateDimReducObject(
  embeddings = dmm@eigenvectors[, 1:2],  # 选择前两个维度
  key = "DC_",
  assay = "RNA"
)

## 去除DC值缺失的细胞
missing_cells <- setdiff(colnames(Mac_10), rownames(dmm@eigenvectors))
Mac_10_DC <- subset(Mac_10, cells = setdiff(colnames(Mac_10), missing_cells))


## Visual
# diffusion map
DimPlot(Mac_10_DC, 
        reduction = "diffusion_map",
        cols = Mac_color,
        pt.size = 1
        ) + 
  coord_cartesian(xlim = c(-0.075, 0.02), ylim = c(-0.05, 0.025)) +
  theme_void() + theme(legend.position = "none")

# pseudotime

palatte=c('#edf1bb',"#fff799", '#f9d3e3',"#ecb0c1", "#b83570")
palatte=c('#7389D5',"#BCC1CF", '#FFFDc5',"#EFB98D", "#E59973")
palatte=c('#82197B',"#B63679", '#E58063',"#FB8760", "#FEC488")
palatte=c('#209A8C',"#28AE7F", '#5FC963',"#ADD830", "#FDE724")


all(rownames(dm_eigen) %in% colnames(Mac_10_DC))  # TRUE
destiny_pseudotime <- dm_eigen$pseudotime
names(destiny_pseudotime) <- rownames(dm_eigen) 

Mac_10_DC$Destiny_pseudotime <- destiny_pseudotime[colnames(Mac_10_DC)]

DimPlot(Mac_10_DC, 
        reduction = "diffusion_map",
        cols = Mac_color,
        pt.size = 1
)+ theme_void() + theme(legend.position = "none")

FeaturePlot(Mac_10_DC, 
            reduction = "diffusion_map",
            features = "Destiny_pseudotime", 
            pt.size = 1) +
  scale_color_gradientn(colours = c("#7F7F85", "#A9ADC3", "#CFCFEC", "#E0B6C6", "#E39D93", "#E26B51")) +
  # scale_color_gradientn(colours = c('#5C569E','#3D87AF','#68BEA2','#A5D5A2',
  #                                   '#EBF0B1','#F6C37E','#F18D5E','#DC5C52')) +
  coord_cartesian(xlim = c(-0.075, 0.02), ylim = c(-0.05, 0.025)) + 
  theme_void() + theme(legend.position = "none",
                       plot.title = element_blank()
                       )



### SCENIC -------------------------------
## 展示STAT3(+)的RSS
Mac_10_DC@meta.data$`STAT3(+)` <- AUC["STAT3(+)",colnames(Mac_10_DC)]

FeaturePlot(Mac_10_DC, features = "STAT3(+)", reduction = "diffusion_map") +
  scale_color_distiller(palette = "Spectral") +
  coord_cartesian(xlim = c(-0.075, 0.02), ylim = c(-0.05, 0.025)) 


## 焦亡打分
# 提取IL1B亚型进行展示
IL1B_DC <- subset(Mac_10_DC,subset = cell_type_level3 %in% c("IL1B+NLRP3+_Macro"))
IL1B_DC$cell_type_level3 <- droplevels(IL1B_DC$cell_type_level3)


Pyroptosis_genes <- c('BAK1','BAX','CASP1','CASP3','CASP4','CASP5','CHMP2A','CHMP2B','CHMP3','CHMP4A','CHMP4B','CHMP4C','CHMP6','CHMP7','CYCS','ELANE','GSDMD','GSDME','GZMB','HMGB1','IL18','IL1A','IL1B','IRF1','IRF2','TP53','TP63')
NFKB_genes <- c('CHUK','FADD','IKBKB','IKBKG','IL1A','IL1R1','MAP3K1','MAP3K14','MAP3K7','MYD88','NFKB1','NFKBIA','RELA','RIPK1','TAB1','TNF','TNFAIP3','TNFRSF1A','TNFRSF1B','TRADD','TRAF6')

geneSets <- list(
  Pyroptosis = Pyroptosis_genes,
  NFkB = NFKB_genes
)


### All Mac AUCell score
Mac10_DC_rankings <- AUCell_buildRankings(GetAssayData(Mac_10_DC, layer = "counts"))

Mac10_DC_AUC <- AUCell_calcAUC(geneSets, Mac10_DC_rankings)

all(colnames(Mac_10_DC) == colnames(Mac10_DC_AUC)) # TRUE

scores <- as.data.frame(t(getAUC(Mac10_DC_AUC)))

Mac_10_DC@meta.data <- cbind(Mac_10_DC@meta.data, scores$NFkB,scores$Pyroptosis)

Mac_10_DC <- NormalizeData(Mac_10_DC, normalization.method = "LogNormalize", scale.factor = 10000)



IL1B_DC <- subset(Mac_10_DC,subset = cell_type_level3 %in% c("IL1B+NLRP3+_Macro"))
IL1B_DC$cell_type_level3 <- droplevels(IL1B_DC$cell_type_level3)

FeaturePlot(IL1B_DC, 
            reduction = "diffusion_map",
            features = 'scores$NFkB', 
            pt.size = 1) +
  # scale_color_gradientn(colours = c("#7F7F85", "#A9ADC3", "#CFCFEC", "#E0B6C6", "#E39D93", "#E26B51")) +
  scale_color_gradientn(colours = c('#5C569E','#3D87AF','#68BEA2','#A5D5A2',
                                    '#EBF0B1','#F6C37E','#F18D5E','#DC5C52')) +
  coord_cartesian(xlim = c(-0.075, 0.02), ylim = c(-0.05, 0.025)) + 
  theme_void() + 
  theme(legend.position = "none",
        plot.title = element_blank()
  )
FeaturePlot(IL1B_DC, 
            reduction = "diffusion_map",
            features = 'scores$Pyroptosis', 
            pt.size = 1) +
  # scale_color_gradientn(colours = c("#7F7F85", "#A9ADC3", "#CFCFEC", "#E0B6C6", "#E39D93", "#E26B51")) +
  scale_color_gradientn(colours = c('#5C569E','#3D87AF','#68BEA2','#A5D5A2',
                                    '#EBF0B1','#F6C37E','#F18D5E','#DC5C52')) +
  coord_cartesian(xlim = c(-0.075, 0.02), ylim = c(-0.05, 0.025)) + 
  theme_void() + 
  theme(legend.position = "none",
        plot.title = element_blank()
  )


FeaturePlot(IL1B_DC, 
            reduction = "diffusion_map",
            features = 'CASP1', 
            pt.size = 1) +
  # scale_color_gradientn(colours = c("#7F7F85", "#A9ADC3", "#CFCFEC", "#E0B6C6", "#E39D93", "#E26B51")) +
  scale_color_gradientn(colours = c('#5C569E','#3D87AF','#68BEA2','#A5D5A2',
                                    '#EBF0B1','#F6C37E','#F18D5E','#DC5C52')) +
  coord_cartesian(xlim = c(-0.075, 0.02), ylim = c(-0.05, 0.025)) + 
  theme_void() + 
  theme(legend.position = "none",
        plot.title = element_blank()
  )










### IL1B Mac AUCell score
IL1B_rankings <- AUCell_buildRankings(GetAssayData(IL1B_DC, layer = "counts"))

IL1B_AUC <- AUCell_calcAUC(geneSets, IL1B_rankings)

all(colnames(IL1B_DC) == colnames(IL1B_AUC)) # TRUE

Pyroptosis_scores <- as.data.frame(t(getAUC(IL1B_AUC)))
NFKB_scores <- as.data.frame(t(getAUC(IL1B_AUC)))

IL1B_DC@meta.data <- cbind(IL1B_DC@meta.data, NFKB_scores$NFkB,Pyroptosis_scores$Pyroptosis)

FeaturePlot(Mac_10_DC, features = 'NFKB_scores$NFkB', reduction = "diffusion_map") +
  scale_color_distiller(palette = "Spectral") +
  coord_cartesian(xlim = c(-0.075, 0.02), ylim = c(-0.05, 0.025)) +
  ggtitle("NFκB signaling pathway")


FeaturePlot(IL1B_DC, 
            reduction = "diffusion_map",
            features = 'NFKB_scores$NFkB', 
            pt.size = 1) +
  scale_color_gradientn(colours = c("#7F7F85", "#A9ADC3", "#CFCFEC", "#E0B6C6", "#E39D93", "#E26B51")) +
  # scale_color_gradientn(colours = c('#5C569E','#3D87AF','#68BEA2','#A5D5A2',
  #                                   '#EBF0B1','#F6C37E','#F18D5E','#DC5C52')) +
  coord_cartesian(xlim = c(-0.075, 0.02), ylim = c(-0.05, 0.025)) + 
  theme_void() + 
  theme(legend.position = "none",
                       plot.title = element_blank()
  )


FeaturePlot(IL1B_DC, 
            reduction = "diffusion_map",
            features = 'Pyroptosis_scores$Pyroptosis', 
            pt.size = 1) +
  # scale_color_gradientn(colours = c("#7F7F85", "#A9ADC3", "#CFCFEC", "#E0B6C6", "#E39D93", "#E26B51")) +
  scale_color_gradientn(colours = c('#5C569E','#3D87AF','#68BEA2','#A5D5A2',
                                    '#EBF0B1','#F6C37E','#F18D5E','#DC5C52')) +
  coord_cartesian(xlim = c(-0.075, 0.02), ylim = c(-0.05, 0.025)) + 
  theme_void() + 
  theme(legend.position = "none",
        plot.title = element_blank()
  )






### slingshot -------------------------------------------------

## Seurat对象转换为sce对象
sce <- as.SingleCellExperiment(Mac_10_DC)

## slingshot
sce <- slingshot(sce, clusterLabels = "cell_type_level3", reducedDim = "DIFFUSION_MAP")

lin1 <- getLineages(sce, 
                    clusterLabels = "cell_type_level3", 
                    start.clus = 'IL1B+NLRP3+_Macro',#可指定起始细胞簇
                    # end.clus=c("IL1B+NLRP3+_Macro","SPP1+CCL2+_Macro"),
                    reducedDim = "DIFFUSION_MAP")
# color
sce$cell_type_level3 <- factor(sce$cell_type_level3, levels = names(momac_color))

# scatter plot (Diffusion map)
plot(reducedDims(sce)$DIFFUSION_MAP, col = momac_color[sce$cell_type_level3], pch=16, asp = 1) 

lines(SlingshotDataSet(lin1), lwd=2, col = 'black',type = 'lineages')

crv1 <- getCurves(lin1)
plot(reducedDims(sce)$DIFFUSION_MAP, col = momac_color[sce$cell_type_level3], pch=16, asp = 1) 
lines(SlingshotDataSet(crv1), lwd = 3, col = 'black')


pseudotime <- slingPseudotime(sce, na = FALSE)
cellWeights <- slingCurveWeights(sce)
crv <- SlingshotDataSet(sce)



### monocle3 -----------------------------------
### 创建CellDataSet对象 (Mac_10_DC)
ct <- GetAssayData(Mac_10_DC, assay = "RNA", layer = 'counts') # count 矩阵 
cell_metadata <- Mac_10_DC@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(ct))
rownames(gene_annotation) <- rownames(ct)

cds <- new_cell_data_set(ct,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation
)

### CellDataSet对象降维
cds <- preprocess_cds(cds, num_dim = 50)
cds <- reduce_dimension(cds, preprocess_method = "PCA")

### 从seurat对象导入整合过的umap坐标
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(Mac_10_DC, reduction = "diffusion_map")
int.embed <- int.embed[colnames(cds),]

cds@int_colData$reducedDims$UMAP <- int.embed

plot_cells(cds, reduction_method="UMAP", color_cells_by="cell_type_level3") + ggtitle('int.umap') # 未加入轨迹

### 计算轨迹
cds <- cluster_cells(cds)
cds <- learn_graph(cds)


cds <- order_cells(cds)
plot_cells(cds, color_cells_by="pseudotime")

















### useless visual -----------------

## 3D
dm_eigen <- data.frame(
  DC1 = dmm@eigenvectors[, 1],
  DC2 = dmm@eigenvectors[, 2],
  DC3 = dmm@eigenvectors[, 3],
  pseudotime = dpt$dpt  # 拟时序数据
)


colors <- viridis(n = nrow(dm_eigen), option = "plasma")[rank(dm_eigen$pseudotime)]

scatterplot3d(dm_eigen$DC1, dm_eigen$DC2, dm_eigen$DC3, 
              color = colors, 
              pch = 16,  # 点的类型
              main = "3D Diffusion Map with Pseudotime",
              xlab = "DC1", ylab = "DC2", zlab = "DC3",
              type = "p",  # 只显示点
              angle = 45,
              theta = 45      # 调整视角
              )  
legend("topright", legend = round(seq(min(dm_eigen$pseudotime), max(dm_eigen$pseudotime), length.out = 5), 2),
       fill = viridis(5, option = "plasma"), title = "Pseudotime")


dm_eigen <- data.frame(
  DC1 = dmm@eigenvectors[, 1],
  DC2 = dmm@eigenvectors[, 2],
  DC3 = dmm@eigenvectors[, 3],
  cell_type = filtered_level3  # 细胞类型
)

scatterplot3d(dm_eigen$DC1, dm_eigen$DC2, dm_eigen$DC3, 
              color = momac_color, 
              pch = 16,  # 点的类型
              main = "3D Diffusion Map with Macrophages",
              xlab = "DC1", ylab = "DC2", zlab = "DC3",
              type = "p",  # 只显示点
              angle = 45,
              theta = 45      # 调整视角
)  





plot_eigenVal <- function(dm=dm){
  linepad <- .5
  plot(
    eigenvalues(dm), 
    ylim = 0:1, 
    pch = 20, 
    xlab ='Diffusion component (DC)', 
    ylab ='Eigenvalue'
  )
}
plot_dm_2D <- function(dm=dm, dc=2, condition=condition, colours=colours){
  DCs <- paste("DC",1:dc, sep="")
  
  dm_eigen <- data.frame(
    dm@eigenvectors
  )
  
  DCs.combinations <- combn(DCs,2)
  g <- apply(
    DCs.combinations,
    2,
    function(combination)
    {
      p1 <- ggplot(dm_eigen, aes_string(x=combination[1], y=combination[2])) +
        geom_point(shape = 21, size = 2.5, stroke=0.5, aes(fill=condition)) +
        theme_bw() +
        scale_fill_manual(
          values=colours,
          name=""
        ) +
        xlab(combination[1])+
        ylab(combination[2])+
        ggtitle("Diffusion Map")+
        theme(
          axis.text=element_text(size=16),
          axis.title=element_text(size=16),
          legend.text = element_text(size =16),
          legend.title = element_text(size =16 ,face="bold"),
          plot.title = element_text(size=18, face="bold", hjust = 0.5),
          aspect.ratio=1
        )
      print(p1)
    }
  )
  
}
plot_dm_3D <- function(dm=dm, dc=c(1:3), condition=condition, colours=colours){
  cond <- factor(condition)
  col <- factor(condition)
  levels(col) <- colours
  col <- as.vector(col)
  DCs <- paste("DC",dc, sep="")
  
  data <- data.frame(
    dm@eigenvectors[,DCs[1]], 
    dm@eigenvectors[,DCs[2]], 
    dm@eigenvectors[,DCs[3]]
  )
  colnames(data) <- DCs
  
  plot3d(
    data,
    col=col,
    size=6.5,
    box = FALSE
  )
  
   legend3d("topright", legend = unique(condition), pch = 21, pt.bg = unique(col), cex=1.5, inset=c(0.02), bty = "n")
  
   plot3d(
    data,
    size=7.5,
    add = TRUE
   )
  
}

condition <- Mac_10@meta.data$cell_type_level3
condition <- as.data.frame(condition)
condition <- cbind(rownames(Mac_10@meta.data),condition)
clustering=condition[,2]
names(clustering)=condition[,1]

plot_eigenVal(dm=dmm)

plot_dm_3D(
  dm=dmm,
  dc=c(1:3),
  condition=clustering,
  colours=momac_color
  )

plot_dm_2D(dm=dmm, dc=2, condition = clustering, colours=)








######################################## Diffusion map IL1B->SPP1 downsampled 10% #################################################

IL1BSPP1_Mac <- subset(Mac_10, subset = cell_type_level3 %in% c("IL1B+NLRP3+_Macro", "SPP1+CCL2+_Macro"))

scvi <- Embeddings(IL1BSPP1_Mac, "scVI")
dmm <- DiffusionMap(scvi)
dpt = DPT(dmm)

## 2D
common_cells <- rownames(dmm@eigenvectors)  # 获取 DiffusionMap 仍然保留的细胞
filtered_level3 <- IL1BSPP1_Mac$cell_type_level3[common_cells]  # 重新匹配细胞类型

dm_eigen <- data.frame(
  DC1 = dmm@eigenvectors[, 1],
  DC2 = dmm@eigenvectors[, 2],
  cell_type = filtered_level3  # 细胞类型
)

ggplot(dm_eigen, aes(x = DC1, y = DC2, fill = cell_type)) +
  geom_point(shape = 21, size = 3, stroke = 0.5) +  # 使用 shape=21 以便填充颜色
  theme_bw() +
  scale_fill_manual(values = momac_color) +  # 使用预定义颜色
  xlab("DC 1") +
  ylab("DC 2") +
  ggtitle("Diffusion Map of Macrophage Subtypes") +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    plot.title = element_text(size = 18, hjust = 0.5)
  )

### DPT pseudotime -------------------------------------------------
dm_eigen <- data.frame(
  DC1 = dmm@eigenvectors[, 1],  # Diffusion component 1
  DC2 = dmm@eigenvectors[, 2],  # Diffusion component 2
  pseudotime = dpt$dpt          # DPT 计算的伪时间
)
ggplot(dm_eigen, aes(x = DC1, y = DC2, color = pseudotime)) +
  geom_point(size = 2) + 
  scale_color_viridis_c(option = "plasma") +  # 颜色渐变方案
  theme_minimal() +
  labs(title = "Diffusion Map with Pseudotime", 
       x = "DC 1", 
       y = "DC 2", 
       color = "Pseudotime") +
  theme(
    plot.title = element_text(size=16, face="bold", hjust = 0.5),
    axis.text = element_text(size=14),
    axis.title = element_text(size=14, face="bold"),
    legend.text = element_text(size=12)
  )

### CytoTRACE -----------------------------------------------------

cytotrace_score <- IL1BSPP1_Mac$CytoTRACE[common_cells]
dm_cytotrace <- data.frame(
  DC1 = dmm@eigenvectors[, 1],  # Diffusion component 1
  DC2 = dmm@eigenvectors[, 2],  # Diffusion component 2
  cytotrace = cytotrace_score # DPT 计算的伪时间
)
ggplot(dm_cytotrace, aes(x = DC1, y = DC2, color = cytotrace)) +
  geom_point(size = 2) + 
  scale_color_viridis_c(option = "plasma") +  # 颜色渐变方案
  theme_minimal() +
  labs(title = "Diffusion Map with CytoTRACE", 
       x = "DC 1", 
       y = "DC 2", 
       color = "CytoTRACE") +
  theme(
    plot.title = element_text(size=16, face="bold", hjust = 0.5),
    axis.text = element_text(size=14),
    axis.title = element_text(size=14, face="bold"),
    legend.text = element_text(size=12)
  )


### Slingshot ----------------------------------------------------
## Diffusion map嵌入回seurat对象
IL1BSPP1_Mac[["diffusion_map"]] <- CreateDimReducObject(
  embeddings = dmm@eigenvectors[, 1:2],  # 选择前两个维度
  key = "DC_",
  assay = "RNA"
)

## 去除DC值缺失的细胞
missing_cells <- setdiff(colnames(IL1BSPP1_Mac), rownames(dmm@eigenvectors))
sce <- subset(IL1BSPP1_Mac, cells = setdiff(colnames(IL1BSPP1_Mac), missing_cells))


## Visual
DimPlot(sce, reduction = "diffusion_map")

## Seurat对象转换为sce对象
sce <- as.SingleCellExperiment(sce)

## slingshot
sce <- slingshot(sce, clusterLabels = "cell_type_level3", reducedDim = "DIFFUSION_MAP")

lin1 <- getLineages(sce, 
                    clusterLabels = "cell_type_level3", 
                    start.clus = 'IL1B+NLRP3+_Macro',#可指定起始细胞簇
                    end.clus=c("IL1B+NLRP3+_Macro","SPP1+CCL2+_Macro"),
                    reducedDim = "DIFFUSION_MAP")
# color
sce$cell_type_level3 <- factor(sce$cell_type_level3, levels = names(momac_color))

# scatter plot (Diffusion map)
plot(reducedDims(sce)$DIFFUSION_MAP, col = momac_color[sce$cell_type_level3], pch=16, asp = 1) 

lines(SlingshotDataSet(lin1), lwd=2, col = 'black',type = 'lineages')

crv1 <- getCurves(lin1)
plot(reducedDims(sce)$DIFFUSION_MAP, col = momac_color[sce$cell_type_level3], pch=16, asp = 1) 
lines(SlingshotDataSet(crv1), lwd = 3, col = 'black')


pseudotime <- slingPseudotime(sce, na = FALSE)
cellWeights <- slingCurveWeights(sce)
crv <- SlingshotDataSet(sce)


### tradeSeq -----------------------------------------------------
counts <- subset(IL1BSPP1_Mac, cells = setdiff(colnames(IL1BSPP1_Mac), missing_cells))@assays$RNA$counts
set.seed(5)
icMat <- evaluateK(counts = counts %>% as.matrix(), sds = crv, k = 3:10, 
                   nGenes = 200, verbose = T)







