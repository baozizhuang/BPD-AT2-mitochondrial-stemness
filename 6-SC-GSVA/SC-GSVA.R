# 0.加载R包=====================================================================

library(GSEABase)
library(clusterProfiler)
library(GSVA)
library(limma)
library(GSEABase)
library(clusterProfiler)
library(msigdbr)
library(Seurat)

# 加载seurat对象
sub_cell <- readRDS("input/AT2_CytoTRACE.rds")
# 设定细胞身份
Idents(sub_cell) <- "CytoTRACE_group"  # 根据不同项目调整
  
# 加载基因集
geneset <- getGmt(file.path("input/mh.all.v2025.1.Mm.symbols.gmt"))

# 计算每个细胞簇的基因的平均表达量并转换为矩阵
exp <- AverageExpression(sub_cell, assays = "RNA", slot = "counts")[[1]]
group_list <- colnames(exp)
exp = exp[rowSums(exp)>0,]   # 过滤低表达基因
exp <- as.matrix(exp)


# 按细胞簇进行gsva分析及可视化
gsva_par <- gsvaParam(exp,
                      geneset,
                      kcdf = "Poisson",
                      maxDiff = FALSE)
es.max <- gsva(gsva_par, verbose = TRUE)  # 显示详细运行信息
gsva.df <- data.frame(Genesets=rownames(es.max),es.max,check.names = F)


# ---------------------------------说明----------------------------------------#
# count，kcdf选择"Poisson"进行计算
# log-normalized kcdf选择 "Gaussian"

# slot = "counts"	原始UMI计数（整数，大量0值）	                         用"Poisson"
# slot = "data"	（连续值，近似正态，适用于经过 log 转换后如log-TPM 等）	 用"Gaussian" 

# 但官方推荐GSVA在单细胞数据中更推荐用Poisson（即使数据已经标准化）
# ---------------------------------说明----------------------------------------#


plot <- pheatmap::pheatmap(es.max,show_colnames = T,
                           width = 13,
                           height = 9,
                           cluster_cols = FALSE,             # 列不聚类
                           scale = "row",
                           angle_col = "0",
                           color = colorRampPalette(c("#4DBBD5","white","#E64B35"))(50),
                           fontsize =12,
                           filename = 'output/5-GSVA.pdf')


