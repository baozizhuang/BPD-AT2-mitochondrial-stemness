
# -----------------------------------------------------------------------------#
#               用于获取单细胞机器学习可用的表达矩阵和分组文件
# -----------------------------------------------------------------------------#

# 1.环境变量准备================================================================
# 1.1 清空环境------------------------------------------------------------------
rm(list=ls())
gc()


# 2.R包安装与加载===============================================================
library(ggplot2)
library(ggrepel)
library(dplyr)
library(plyr)
library(magrittr)
library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(tibble)
library(pheatmap)



# 3. 提取需要进行机器学习的单细胞分组文件=======================================
# 读取AT2_CytoTRACE.rds
sc_integer <- readRDS("input/AT2_CytoTRACE.rds")


# 降采样
sc_integer <- SketchData(
         object = sc_integer,
         assay = "RNA",           # 使用RNA assay（或SCT等）
         ncells = 400,            # 采样400个细胞
         method = 'Uniform',      # 均匀采样方法
         sketched.assay = 'sketch_uniform',  # 创建新的sketch assay
         seed = 123               # 设置随机种子确保可重现
  )

# 提取降采样结果并赋值到新变量
downsample_cells <- Cells(sc_integer[["sketch_uniform"]])
sc_integer <- subset(sc_integer, cells = downsample_cells)

# 使用均匀采样后需要重新设置默认assay为RNA
DefaultAssay(sc_integer) <- "RNA"



# 提取高干性和低干性组的细胞用于机器学习
Idents(sc_integer) <- sc_integer$CytoTRACE_group
High_Low_stemness <- subset(sc_integer, idents = c("High_stemness", "Low_stemness"))
DimPlot(High_Low_stemness)


# 从Seurat对象提取表达矩阵和分组标签
exp <- as.matrix(GetAssayData(High_Low_stemness, layer = "data"))  # 基因×细胞矩阵
group <- data.frame(
  cell = colnames(High_Low_stemness),
  group = High_Low_stemness$CytoTRACE_group  # "High_stemness" / "Low_stemness"
)
rownames(group) <- group$cell


# 读取共表达模块基因
module_genes <- read.csv("input/Sig_Modules_Top20_genes.csv")$gene_name
# 取与表达矩阵的交集
log.res <- intersect(module_genes, rownames(exp))



# 对齐细胞顺序（确保exp的列与group的行顺序一致）
exp <- exp[, rownames(group)]

# 导出对齐后的表达矩阵和分组标签
write.csv(exp[log.res, ], file = "output/ML_expression_matrix.csv", quote = FALSE)
write.csv(group, file = "output/ML_group_labels.csv", quote = FALSE, row.names = FALSE)















