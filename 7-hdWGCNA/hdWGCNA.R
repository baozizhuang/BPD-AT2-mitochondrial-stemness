# -----------------------------------------------------------------------------#
# 描述: 本脚本用于执行单细胞转录组的hdWGCNA的分析 hdWGCNA.R
# R版本：4.3.2
# Seurat (Version:5.0.1)
# 日期: 2026/07/07
# 创建：阿伦
# 代码版本：V2
#
# 存在问题：会R包加载冲突，导致分析完后的数据再次打开无法跑，得删除所有环境重新加载已经保存的数据
# -----------------------------------------------------------------------------#

# rm(list = ls())

# 添加镜像======================================================================
# 添加bioconductor清华镜像
options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
# 添加bioconductor西湖大学镜像
options(BioC_mirror = "https://mirrors.westlake.edu.cn/bioconductor")
# 添加R清华镜像
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))



# 配置动态R包库,我服务器安装R包没权限，得设置动态R包库(服务器)==================
# 我的个人R包库在/work/run/projects/bio-35/R_library
.libPaths(c("/work/run/projects/bio-35/R_library", .libPaths()))



# =============================================================================#
#                             0. 初始化设置
# =============================================================================#
R.version

options(stringsAsFactors = FALSE) # 禁用字符串自动转为因子

# 从GitHub 安装 hdWGCNA 包的开发版本
# devtools::install_github('smorabit/hdWGCNA', ref='dev',.libPaths = "")
# devtools::install_github('smorabit/hdWGCNA', ref='dev')


# 加载R包
library(dplyr)  
library(ggplot2)
library(Seurat)
library(harmony)
library(WGCNA)
library(hdWGCNA)
library(tidyverse)
library(cowplot)
library(patchwork)
library(qs)
library(GeneOverlap)
library(magrittr)


# 并行计算设置
enableWGCNAThreads(nThreads = 30) # 根据CPU核心数调整
set.seed(12345) # 随机种子

# 增加允许的对象大小限制
options(future.globals.maxSize = 2000 * 1024^2)  # 设置为 2000 MB (2 GB)


# =============================================================================#
#                                1. 数据准备
# =============================================================================#
# 加载Seurat对象
scRNA <- readRDS("input/AT2_CytoTRACE.rds") 
# 设置细胞标签身份
Idents(scRNA) <- "CytoTRACE_group" 

# 数据质量检查
if(!"umap" %in% names(scRNA@reductions)){
  stop("未找到UMAP降维结果，请先运行RunUMAP()")
}

# 可视化细胞聚类
DimPlot(scRNA, reduction = 'umap', label = TRUE) + NoLegend()



# =============================================================================#
#                  2. hdWGCNA初始化(为hdWGCNA 设置 Seurat 对象)
# =============================================================================#
# 选择要用于分析的背景基因
seurat_obj <- SetupForWGCNA(
  scRNA,
  gene_select = "fraction",    # 基因选择策略：fraction(表达比例), variable(高变基因)   custom(自定义)
  fraction = 0.05,             # 基因至少在5%细胞中表达
  wgcna_name = "MyWGCNA"       # 分析项目名称
)



# =============================================================================#
#                          3. 构建MetaCells(核心步骤)
# =============================================================================#

# 先查看下每个分组的具体细胞数量
table(seurat_obj$CytoTRACE_group, seurat_obj$sample)

# 构建MetaCells
seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("CytoTRACE_group", "sample"), # 分组依据(根据数据调整),2个标签分组是为了双重控制,sample改成orig.ident也可以
  reduction = 'harmony',             # 使用harmony降维结果,使用已存在的Harmony降维坐标
  k = 25,                            # KNN邻居数(影响MetaCells数量)
  max_shared = 10,                   # MetaCells间共享细胞上限
  ident.group = 'CytoTRACE_group'    # 细胞标识字段,指定使用Seurat对象中已存在的哪一列元数据来标记每个MetaCell的身份
)

# MetaCells标准化
seurat_obj <- NormalizeMetacells(seurat_obj)


# =============================================================================#
#                              4. 共表达网络构建
# =============================================================================#

# ---------------------------------解读----------------------------------------#
# 为什么hdWGCNA不同于bulk RNAseq的WGCNA,需要指定细胞进行后续分析
# Bulk WGCNA：一个样本=一个混合组织的表达谱,所有样本直接构建一个共表达网络
# scRNA hdWGCNA：数据里混合了 多种细胞类型 (如 T 细胞、上皮细胞、巨噬细胞等)
#                不同细胞类型的基因表达程序完全不同，如果混合所有细胞一起建网 ，
#                模块会被细胞类型之间的差异"主导"，检测到的其实是"哪些基因在A细胞高表达、在B细胞低表达"，
#                这不是真正的共表达关系
# hdWGCNA 的设计理念是为感兴趣的细胞类型单独建网 
# 如果要多个细胞一起构建也不是不行,具体看官方教程
# ---------------------------------解读----------------------------------------#

# 设置目标细胞群(根据方案选择要关注的细胞类型)
target_celltype <- "High_stemness" # 选择要聚焦分析的细胞类型，重要！！！！！！！！！
seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = target_celltype, 
  group.by = 'CytoTRACE_group',    # 选择分组标签,和上述对应,即去 metadata 的哪一列找标签
  assay = 'RNA',
  layer = 'data'                   # 使用标准化后的数据
)

# 软阈值筛选|(关键步骤)|
seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = 'signed' # 网络类型：unsigned/signed/signed hybrid  官方教程推荐使用'signed',可区分基因间的正/负调控关系。
)

power_plot <- PlotSoftPowers(seurat_obj)                           # 创建软阈值评估图
p_power <- cowplot::plot_grid(plotlist = power_plot, ncol = 2)     # 将多个图排列成2列显示
p_power

# 保存软阈值评估图
pdf("output/0-SoftPowers.pdf", width = 8, height = 6)
print(p_power)
dev.off()


# ---------------------------------解读----------------------------------------#
# 软阈值选择解读(PlotSoftPowers输出4张诊断图) (看图1就行)
#
#  ┌─────────────────────────────────────┬─────────────────────────────────────┐
#  │ 图1(左上) Scale Free Topology Fit   │ 图2(右上) Mean Connectivity         │
#  │ 横轴:Soft Power | 纵轴:R²           │ 横轴:Soft Power | 纵轴:平均连接度   │
#  │ 红线:R²=0.8阈值                     │ 随Power增大连接度下降               │
#  │ 选首个R²≥0.8的最小整数作为软阈值    │ 需平衡网络稀疏性与信息保留          │
#  ├─────────────────────────────────────┼─────────────────────────────────────┤
#  │ 图3(左下) Median Connectivity       │ 图4(右下) Max Connectivity          │
#  │ 同图2但展示中位数,作为鲁棒性验证    │ 确保hub gene不会因Soft Power        │
#  │                                     │ 过大而被过度压制                    │
#  └─────────────────────────────────────┴─────────────────────────────────────┘
# 若所有R²均<0.8,则选择曲线趋于平稳的最小power值
# R²阈值默认0.8,0.8是金标准
#
# 出的模块推荐15-20个，如太多建议人工调整
# “灰色”模块由未被分组到任何共表达模块中的基因组成。对于所有下游分析和解释，应忽略灰色模块。
# ---------------------------------解读----------------------------------------#


# (选择1)自动构建共表达网络(自动选择软阈值)(默认参数)
seurat_obj <- ConstructNetwork(
  seurat_obj,
  tom_name = target_celltype,   # TOM 矩阵写入磁盘时的文件名前缀
  overwrite_tom = TRUE          # 允许覆盖已存在的同名文件
)


# (选择2)自动构建共表达网络(自动选择软阈值)(调整后的参数)
# 为了减少模块数,过多模块不方便画图
# seurat_obj <- ConstructNetwork(
#   seurat_obj,
#   tom_name = target_celltype,
#   overwrite_tom = TRUE,
#   minModuleSize = 80,          # ← 增大最小模块大小（默认一般30-50）
#   mergeCutHeight = 0.25        # ← 提高合并阈值，让相似模块自动合并（默认0.15）
# )


# 可视化层次聚类树
pdf("output/0-PlotDendrogram.pdf", width = 10, height = 6)
PlotDendrogram(seurat_obj,
              main = paste(target_celltype, 'Dendrogram'))
dev.off()


# 可选：检查topoligcal重叠矩阵(TOM)
# TOM <- GetTOM(seurat_obj)
# TOM


# 保存seurat_obj
# saveRDS(seurat_obj,"output/seurat_obj.rds")
# 读取seurat_obj
# seurat_obj <- readRDS("output/seurat_obj.rds")


# =============================================================================#
#                    5. 模块分析----模块特征基因和连通性
# =============================================================================#


# -------------------------------概念解读----------------------------------------#
# Module : 模块 = 在干性高的 AT2 细胞中高度协同表达的基因集合，即目标细胞构建的共表达网络中的协同表达的基因集合
# ME : 模块特征基因（Module Eigengenes） = 模块特征向量 = 模块特征超级基因 =  PCA第一主成分,一个细胞在一个模块上的活性值，这个模块在这个细胞中的活跃程度
#      模块特征基因是用于总结整个共表达模块基因表达谱的常用指标。
#      模块特征基因是通过对方差分析中每个模块的基因表达矩阵子集进行主成分分析（PCA）来计算的,这些 PCA 矩阵的第一主成分即为 ME。
# hME : 去批次后的模块特征基因 (harmonized module eigengenes) 
# hMEs : hME 的复数形式,PCA第一主#成分,所有细胞各自在所有模块上的活性值
# 
# kME ： 在模块内,每个基因与模块特征基因（ME）的相关系数 = 连通性 = 连通度
# kMEs ： kME的复数形式，即所有基因各自在所有模块上的kME值矩阵
# -------------------------------概念解读----------------------------------------#

# 计算去批次后的模块特征基因 (hME)
seurat_obj <- ModuleEigengenes(
  seurat_obj,
  group.by.vars = "sample"   # 批次校正字段orig.ident或sample,去批次后再算 hME
)

# 计算基因与模块的连接度(kME)
seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by = 'CytoTRACE_group',    # 选择分组标签,和上述对应,即去 metadata 的哪一列找标签
  group_name = target_celltype
)

# 重命名模块,加前缀(可选)
seurat_obj <- ResetModuleNames(seurat_obj, new_name = "Module")

# 提取模块信息,并去掉 grey 模块
modules <- GetModules(seurat_obj) %>% 
  subset(module != 'grey') # 移除未分配基因

# 可视化模块连接度(kME)
kme_plot <- PlotKMEs(seurat_obj, ncol = 8)  # 可视化每个模块中按 kME 排名的基因
ggsave("output/1-KME_plot.pdf", kme_plot, width = 27, height = 14)



# =============================================================================#
#                 6. 获取每个模块的hub基因与模块活性评分分析
# =============================================================================#

# 获取每个模块的Top hub基因
hub_genes <- GetHubGenes(seurat_obj, n_hubs = 20)   # kME 基因连接度排名前 20 的 hub 基因
  
# 模块基因评分(ModuleExprScore)
seurat_obj <- ModuleExprScore(
  seurat_obj,
  n_genes = 20,    # 和上面一致,Top20基因,用每个模块中TOP20的hub gene 来打分
  method = 'UCell' # 可选：UCell/Seurat/SCTransform
)

# ---------------------------------解读----------------------------------------#
# 至此,seurat_obj中包含了每个模块的两种活性指标:
#   1. hMEs (去批次后模块特征基因的值): PCA第一主成分,所有细胞各自在所有模块上的活性值，适合模块间相关性分析,定位模块在细胞图谱中的空间分布
#   2. scores (模块基因评分): UCell富集评分,适合FeaturePlot可视化,验证核心基因在单细胞水平的实际活性
# ---------------------------------解读----------------------------------------#


# 保存seurat_obj
# saveRDS(seurat_obj,"output/seurat_obj.rds")
# 读取seurat_obj
seurat_obj <- readRDS("output/seurat_obj.rds")



# =============================================================================#
#                            7. 核心可视化
# =============================================================================#

# 7.1 模块特征值图(基于hMEs)----------------------------------------------------
# 通常不用这个
feature_plot1 <- ModuleFeaturePlot(
  seurat_obj,
  features = 'hMEs',
  order = TRUE # 按表达强度排序
)
wrap_plots(feature_plot1, ncol = 5)
ggsave("output/2-Module_hMEs_features.pdf", plot = wrap_plots(feature_plot1, ncol = 5),
       width = 15, height = 16)


# 7.2 模块特征图(基于scores) 基于hub gene signature score-----------------------
# 一般来说通常用这个scores指标的图，而不是hMEs的
feature_plot2 <- ModuleFeaturePlot(
  seurat_obj,
  features = 'scores',
  order='shuffle', 
  ucell = TRUE 
)
wrap_plots(feature_plot2, ncol = 5)
ggsave("output/3-Module_scores_features.pdf", plot = wrap_plots(feature_plot2, ncol = 5),
       width = 15, height = 14)


# 7.3 模块雷达图(多细胞类型比较)------------------------------------------------
# 使用场景1：用亚型某个簇构建了这个簇的共表达网络,然后展示所有亚型的雷达图  思路1：亚型构建网络的，哪些模块是高干性细胞特有的调控程序
# 使用场景2：用某个细胞大类构建了个簇的共表达网络,然后展示所有亚型的雷达图  思路2：大类构建网络的，在大类通用的调控程序中，哪些模块在不同亚型中差异化活跃

radar_plot <- ModuleRadarPlot(
  seurat_obj,
  group.by = 'CytoTRACE_group'
)
radar_plot
ggsave("output/4-Radar_plot.pdf", radar_plot, width = 40, height = 29)



# =============================================================================#
#                       8. 网络高级可视化(按需运行)
# =============================================================================#

# ---------------------------------说明----------------------------------------#
# hdWGCNA提供了丰富的共表达网络可视化
# 
# https://smorabit.github.io/hdWGCNA/articles/network_visualizations.html
# 
# ModuleNetworkPlot： 为每个模块单独绘制网络图，展示按 kME 值排序的Top基因。
# HubGeneNetworkPlot：绘制包含所有模块的网络图，每个模块显示指定数量的枢纽基因。
# ModuleUMAPPlot：    使用UMAP算法同时可视化共表达网络中的所有基因。
# ---------------------------------说明----------------------------------------#


# 8.1 单个模块(module)的网络图--------------------------------------------------

# ---------------------------------说明----------------------------------------#
# 在这个网络中，每个节点代表一个基因，每条边代表网络中两个基因之间的共表达关系。
# 每个模块网络图都根据 hdWGCNA 模块分配表 GetModules(seurat_obj) 中的 color 列进行着色。
# 默认按 kME 排名的前 10 个核心基因位于图的中心，其余 15 个基因位于外圈。
# 如果设置了自定义值那根据自定义值看
# ---------------------------------说明----------------------------------------#

# 绘制所有模块(module)的网络，并将文件导出储存在ModuleNetworks文件夹
ModuleNetworkPlot(                   # 默认可视化每个模块中前 25 个核心基因所对应的网络
  seurat_obj,                   
  n_inner = 6,                       # 内圈基因数
  n_outer = 14,                      # 外圈基因数
  outdir = 'output/ModuleNetworks'   # 输出的文件夹
)


# 8.2 组合hub基因网络图---------------------------------------------------------

# 该函数指定的前n个hub基因以及其他随机选择的基因，使用力导向图绘制算法构建一个联合网络
# 无法编辑图片中的基因位置

HubGeneNetworkPlot(
  seurat_obj,
  n_hubs = 3,        # 每个module需要展示的hub基因数量
  n_other=6,         # 每个module需要展示的非hub基因(随机)数量
  edge_prop = 0.25,  # 对网络中的边数进行降采样，保留网络中n%的边
  mods = "all"       # 展示的模块,all:所有
)


# 8.3 模块UMAP投影--------------------------------------------------------------
# 非常费时,不建议
# seurat_obj <- RunModuleUMAP(
#     seurat_obj,
#     n_hubs =  10,      # 每个模块选择 top 10 个 hub 基因
#     n_neighbors= 15,   # UMAP 的邻居数，控制局部结构
#     min_dist= 0.1,     # UMAP 的最小距离，控制点的聚集程度
#     supervised = TRUE  # 是否使用细胞类型信息
#   )
# 
# 
# umap_plot <- ModuleUMAPPlot(
#   seurat_obj,
#   edge.alpha= 0.25,       # 边的透明度
#   sample_edges= TRUE,     # 是否对边进行降采样
#   edge_prop= 0.1,         # 保留边的比例（10%）
#   label_hubs= 2 ,         # 每个模块标注几个 hub 基因
#   keep_grey_edges= FALSE  # 是否保留灰色边（模块间边）
#   
# )



# =============================================================================#
#                      9 模块-表型（干性状态）关联分析
#                    本代码只适合单个数值类表型,如干性得分
# =============================================================================#


# 9.1 模块与表型相关性分析------------------------------------------------------
cur_traits <- "CytoTRACE"       # 选择的表型,这里是干性得分

# 计算模块与干性状态的相关性
seurat_obj <- ModuleTraitCorrelation(
  seurat_obj,
  traits = cur_traits,
  group.by = NULL,              # 按所有细胞混合展示目标表型的相关性
  # group.by= "CytoTRACE_group" # 按细胞分类展示各个细胞类型各自与目标表型的相关性
  features = "hMEs",            # "hMEs", "MEs", "scores"  使用哪个特性来总结每个模块？
  cor_method = "spearman"       # "pearson", "spearman", "kendall"
)

# ModuleTraitCorrelation 会输出 fdr pval cor
# 我们默认使用pval作为显著性值

# 提取相关性结果
mtc_results <- GetModuleTraitCorrelation(seurat_obj)

cor_df <- data.frame(
  Module = colnames(mtc_results$cor$all_cells),
  Correlation = as.numeric(mtc_results$cor$all_cells),
  Pvalue = as.numeric(mtc_results$pval$all_cells)
) %>%
  mutate(Significance = case_when(
    Pvalue < 0.001 ~ "***",
    Pvalue < 0.01 ~ "**",
    Pvalue < 0.05 ~ "*",
    TRUE ~ ""
  ))

# 可视化
p <- ggplot(cor_df, aes(x = reorder(Module, Correlation), y = Correlation, fill = Correlation)) +
  geom_bar(stat = "identity", color = "black", size = 0.3) +
  geom_text(aes(label = Significance), 
            hjust = ifelse(cor_df$Correlation > 0, -0.2, 1.2),  
            size = 4, fontface = "bold") +
  scale_fill_gradient2(low = "#4DBBD5", mid = "white", high = "#E64B35", midpoint = 0) +
  coord_flip() +
  labs(x = "Module", y = "Correlation with CytoTRACE Score", title = "Module Correlation") +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "none") +
  expand_limits(y = max(abs(cor_df$Correlation)) * 1.15)
ggsave("output/5-Module Correlation.pdf", p, width = 6, height = 5)



# 9.2 筛选显著且强相关模块（p < 0.05 且 |cor| > 0.5）---------------------------
sig_modules <- cor_df %>%
  filter(Pvalue < 0.05, abs(Correlation) > 0.5) %>%
  arrange(desc(abs(Correlation)))
print(sig_modules, row.names = FALSE)


# 9.3 提取显著且强相关模块的 kME TOP20的基因------------------------------------

sig_module_names <- sig_modules$Module

# 从 hub_genes 中筛选显著模块的 Top20 基因
sig_hub_genes <- hub_genes %>%
  filter(module %in% sig_module_names)
print(sig_hub_genes, row.names = FALSE)



# =============================================================================#
#                           10. 结果保存
# =============================================================================#

# 保存显著且强相关模块的kME TOP20的基因（仅第一列基因名）
write.csv(sig_hub_genes[, 1, drop = FALSE], "output/Sig_Modules_Top20_genes.csv", row.names = FALSE)

# 保存hdWGCNA对象
qs::qsave(seurat_obj, file = 'output/hdWGCNA_results.qs')

# 保存关键结果表格
write.csv(modules, "output/unfiltered/All_WGCNA_modules.csv")
write.csv(hub_genes, "output//unfiltered/All_Hub_genes.csv")

# 保存所有模块的每个模块的top基因为列表
write.csv(hub_genes[, 1, drop = FALSE], "output/unfiltered/Hub_genes_list.csv", row.names = FALSE)

# 会话信息记录
capture.output(sessionInfo(), file = "output/session_info.txt")




# =============================================================================#
#                               11. 备用
# =============================================================================#


# 官方绘制热图函数,但效果不理想
# # 可视化模块-性状关联热图
# mtc_plot <- PlotModuleTraitCorrelation(
#   seurat_obj,
#   label = 'pval',          # fdr矫正  pval  fdr
#   label_symbol = 'stars', # 显示p值标签为“星星”或为“数字”
#   text_digits = 2,        # 保留p值的两位小数
#   text_size = 3           # 调整数值字体大小
# )
# 
# ggsave("output/5-Module_Trait_Correlation.pdf", mtc_plot, width = 10, height = 12)

