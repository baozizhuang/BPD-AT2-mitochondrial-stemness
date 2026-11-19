# -----------------------------------------------------------------------------#
# 描述: 本脚本用于执行虚拟敲除分析 scTenifoldKnk.R
# R版本：4.3.2
# Seurat (Version:5.0.1)
# 日期: 2026/07/14
# 修改补充：阿伦
# 代码版本：V2
# -----------------------------------------------------------------------------#


# install_github("wjawaid/enrichR")


# 添加bioconductor清华镜像
options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
# 添加bioconductor西湖大学镜像
options(BioC_mirror = "https://mirrors.westlake.edu.cn/bioconductor")

# 添加R清华镜像
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))


# 配置动态R包库,我服务器安装R包没权限，得设置动态R包库--------------------------
# 我的个人R包库在/work/run/projects/bio-35/R_library
# 每次得运行！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！
.libPaths(c("/work/run/projects/bio-35/R_library", .libPaths()))



#==============================================================================#
#                           0. 安装并加载所需包
#==============================================================================#
library(remotes)
# install_github('cailab-tamu/scTenifoldKnk')
# 
# #或者用devtools工具安装
# # install.packages("devtools")
# devtools::install_github("josephzaccaria/scTenifoldKnk")


library(igraph)
library(dplyr)
library(Seurat)
library(scTenifoldKnk)
library(scTenifoldNet)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(Matrix)
library(reshape2)
library(SeuratData)


#----------------------说明  (Ciallo～(∠・ω< )⌒★)----------------------------#

# 1.使用的单细胞数据如果没有说明就使用全部细胞类型的数据，如果说了那就用特定的细胞提取来做
# 2.推荐  抽样细胞/选择高变基因  来做
# 3.可以使用纯WT(野生型)/疾病对照混合的单细胞seurat对象来做

# 本项目聚焦在VSMC细胞

#----------------------说明  (Ciallo～(∠・ω< )⌒★)----------------------------#


#==============================================================================#
#             1. 选择敲除基因，读单细胞数据，进行数据预处理
#==============================================================================#


# 先设定要敲除的目标基因
gKO <- "Rps3"


# 1.1 使用全细胞数据集进行的预处理(选择1）--------------------------------------

# 1.1.1 导入 Seurat 对象
sc_integer <- readRDS("input/AT2_CytoTRACE.rds")

# 1.1.2 查看对象信息
sc_integer

# 1.1.3 归一化 （不需要，但教程有）
if (!"data" %in% Layers(sc_integer)) {
  sc_integer <- NormalizeData(sc_integer, normalization.method = "LogNormalize", scale.factor = 10000)
  message("数据已归一化")
} else {
  message("数据已经归一化，跳过此步骤")
}

# 1.1.4 在提取5000细胞的基础上提取5000高变基因的counts矩阵~~~~~~~~~~~~~~~~~~~~~~
DefaultAssay(sc_integer) <- "RNA"

table(sc_integer$CytoTRACE_group) # 查看细胞类型和数量


# #-------------------随机抽样来抽样5000个细胞(方法1)----------------------------#
# set.seed(123)
# sc_integer_5000cells <- subset(sc_integer,
#                                cells = sample(colnames(sc_integer),
#                                size = 5000))
# #------------------------------------------------------------------------------#


#----------使用 SketchData 均匀采样抽样5000个细胞作为备选方案(方法2)-----------#
sc_integer <- SketchData(
  object = sc_integer,
  assay = "RNA",           # 使用RNA assay（或SCT等）
  ncells = 5000,           # 采样5000个细胞
  method = 'Uniform',      # 均匀采样方法
  sketched.assay = 'sketch_uniform',  # 创建新的sketch assay
  seed = 123               # 设置随机种子确保可重现
)

# 提取降采样结果并赋值到新变量
downsample_cells <- Cells(sc_integer[["sketch_uniform"]])
sc_integer_5000cells <- subset(sc_integer, cells = downsample_cells)

# 使用均匀采样后需要重新设置默认assay为RNA
DefaultAssay(sc_integer_5000cells) <- "RNA"
#------------------------------------------------------------------------------#


# 查找5000高变基因
sc_integer_5000cells <- FindVariableFeatures(sc_integer_5000cells,
                                             selection.method = "vst",
                                             nfeatures = 5000)


# 提取5000细胞的5000高变基因的counts矩阵~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

h5000g <- VariableFeatures(sc_integer_5000cells)

# 强制添加目标基因到高变基因列表（如果不在其中）
if (!gKO %in% h5000g) {
  # 移除一个基因，为目标基因腾出位置
  h5000g <- h5000g[-length(h5000g)]
  # 添加目标基因
  h5000g <- c(gKO, h5000g)
  message("已将 ", gKO, " 添加到高变基因列表")
}

# 提取矩阵(如果没有强行加入目标基因，就直接提取)
countMatrix <- GetAssayData(sc_integer_5000cells[["RNA"]], layer = "counts")[h5000g, ]  # 或layer提取
countMatrix <- as.data.frame(countMatrix)  # 解决维度错误,不是必须

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# 查看维度和前几个基因
dim(countMatrix)
head(rownames(countMatrix))[1:5]

# 检查 gKO 基因的表达 是否全为 0
sum(countMatrix[gKO, ])  # 如果是 0，说明它根本没表达 → 网络肯定没有


# 1.2 提取指定细胞类型数据集进行的预处理(选择2)---------------------------------
# 未写，待以后补充



#==============================================================================#
#                  2. 检查敲除基因并运行scTenifoldKnk
#==============================================================================#

# 检查基因是否在高变基因中
gKO %in% rownames(countMatrix)


#运行 scTenifoldKnk 虚拟敲除分析
result <- scTenifoldKnk(
  countMatrix = countMatrix,
  gKO = gKO,
  qc = TRUE,             # 是否进行质量控制   如果目标敲除基因的表达太低，可以不质控 qc = FALSE   qc = TRUE
  qc_mtThreshold = 0.1,   # 线粒体比例阈值  qc = TRUE才生效
  qc_minLSize = 1000,     # 最小文库大小    qc = TRUE才生效  常规是1000,可以500
  nc_nNet = 10,           # 构建多少个子网络
  nc_nCells = 300,        # 每个子网络包含的细胞数量（建议 300-500）
  nc_nComp = 3,           # PCA 主成分数量
  nCores = 30             # 多线程
)


#==============================================================================#
#                         3. 差异调控结果可视化
#==============================================================================#

# 3.1 Top 20 基因柱状图---------------------------------------------------------
# 根据调整后的 p.adj 值排序并选择前 20 个基因
top_genes <- head(result$diffRegulation[order(result$diffRegulation$p.adj), ], 20)

# 计算 -log10(p.adj) 以便更直观地展示
top_genes$log_pval <- -log10(top_genes$p.adj)

# 绘制柱状图
pdf(paste0("Top20_", gKO, "_Knockout_HVGs.pdf"), width = 8, height = 6)
ggplot(top_genes, aes(x = reorder(gene, log_pval), y = log_pval)) +
  geom_bar(stat = 'identity', fill = 'steelblue') +
  coord_flip() +
  labs(title = paste0("Top 20 Differentially Regulated Genes (KO ", gKO, ", HVGs only)"),
       x = "Gene", y = "-log10(p.adj)") +
  theme_minimal()
dev.off()


# 3.2 Z-score vs -log10(p.adj) 散点图 （选择1）---------------------------------
# 如果筛选后结果中的基因过多，就用选择1

df <- result$diffRegulation
df$Color <- ifelse(abs(df$Z) > 2 & df$p.adj < 0.05, "#9F0000", "gray50")
df$log_pval <- -log10(df$p.adj)
label_genes <- subset(df, abs(Z) > 2 & p.adj < 0.05)

pdf(paste0("Differential Regulation after ", gKO, " Knockout_HVGs.pdf"), width = 8, height = 6)
ggplot(df, aes(x = Z, y = log_pval)) +
  geom_point(aes(color = Color), alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "blue") +
  geom_label_repel(data = label_genes, aes(label = gene),
                   size = 2.5, max.overlaps = 50, direction = "both",
                   box.padding = 0.4, point.padding = 0.3,
                   segment.color = "black", segment.size = 0.3,
                   segment.alpha = 0.6, min.segment.length = 0,
                   force = 2, force_pull = 1) +
  scale_color_identity() +
  coord_cartesian(clip = "off") +
  labs(title = paste0("Differential Regulation after ", gKO, " Knockout"),
       x = "Z-score", y = "-log10(p.adj)") +
  theme_classic()
dev.off()


# 3.2 Z-score vs -log10(p.adj) 散点图（选择2）----------------------------------
# 如果筛选后结果中的基因不多，就用选择2

df <- result$diffRegulation
df$Color <- ifelse(abs(df$Z) > 2 & df$p.adj < 0.05, "#9F0000", "gray50")
df$log_pval <- -log10(df$p.adj)
label_genes <- subset(df, abs(Z) > 2 & p.adj < 0.05)

pdf(paste0("Differential Regulation after ", gKO, " Knockout_HVGs.pdf"), width = 8, height =6)
ggplot(df, aes(x = Z, y = log_pval)) +
  geom_point(aes(color = Color), alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "blue") +
  geom_text_repel(data = label_genes,aes(label = gene),
                  size = 3,max.overlaps = 50, direction = "y",
                  nudge_y = -1) +
  scale_color_identity() +
  coord_cartesian(clip = "off") +
  labs(title = paste0("Differential Regulation after ", gKO, " Knockout"),
       x = "Z-score",y = "-log10(p.adj)")+theme_classic()
dev.off()


# 3.3 柱状图2 Z-score vs -log10(p.adj)（显著差异显著基因过少使用） -------------
# 根据|Z-score|排序取Top 20基因（包含显著和不显著的）
top20_genes <- head(result$diffRegulation[order(-abs(result$diffRegulation$Z)), ], 20)
# 标记显著性：p.adj < 0.05 为显著
top20_genes$significance <- ifelse(top20_genes$p.adj < 0.05, "p.adj < 0.05", "p.adj >= 0.05")

pdf(paste0("Zscore_Barplot_", gKO, "_Knockout_HVGs.pdf"), width = 7, height = 4.3)
ggplot(top20_genes, aes(x = reorder(gene, Z), y = Z, fill = significance)) +
  geom_bar(stat = 'identity', width = 0.7) +
  scale_fill_manual(values = c("p.adj < 0.05" = "#9F0000", "p.adj >= 0.05" = "gray50")) +
  coord_flip() +
  labs(title = paste0("Top 20 Differentially Regulated Genes (KO ", gKO, ")"),
       x = "Gene", y = "Z-score", fill = "Significance") +
  theme_classic()
dev.off()



#==============================================================================#
#                         4. 绘制敲除基因的局部调控网络
#==============================================================================#

# 检查Rps3是否在差异调控结果中
"Rps3" %in% result$diffRegulation$gene

# 查看Rps3的详细调控信息
result$diffRegulation[result$diffRegulation$gene == "Rps3", ]

##打开代理
Sys.setenv(
  http_proxy = "http://localhost:1086", 
  https_proxy = "http://localhost:1086"
)
#加载包和函数
library(enrichR)  # 这个包加载容易超时，可以通过加载上一段的代理解决
source("plotKO.R")

# 函数默认选择了所有满足p.adj < 0.05的扰动基因
pdf(paste0("Egocentric_KO_", gKO, "_HVGs.pdf"), width = 20, height = 10)
plotKO(X = result, gKO = gKO)
dev.off()

##关闭代理
Sys.unsetenv(c("http_proxy", "https_proxy"))



#==============================================================================#
#                             5. 查看差异调控结果
#==============================================================================#
head(result$diffRegulation)

saveRDS(result, paste0("Egocentric_KO_", gKO, "_result.rds"))
result <- readRDS(paste0("Egocentric_KO_", gKO, "_result.rds"))

save.image(file = "KO.RData")
