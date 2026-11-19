# -----------------------------------------------------------------------------#
# 描述: 本脚本用于执行单细胞转录组亚群分析 SubType.R
# R版本：4.3.2
# Seurat (Version:5.0.1)
# 日期: 2026/04/23
# 创建：阿伦
# 代码版本：V2
# -----------------------------------------------------------------------------#


# 配置动态R包库,我服务器安装R包没权限，得设置动态R包库--------------------------
# 我的个人R包库在/work/run/projects/bio-35/R_library
# 每次得运行！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！
.libPaths(c("/work/run/projects/bio-35/R_library", .libPaths()))


# 0.R包的安装与加载=============================================================
library(Seurat)
library(tidyverse)
library(msigdbr)
library(GSVA)
library(Seurat)
library(magrittr)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(GEOquery)
library(harmony)
library(clustree)
library(celldex)
library(SingleR)
library(ExperimentHub)
library(httr)
library(biomaRt)
library(ggsci)
library(ggplot2)
library(GEOquery)
library(data.table)
library(Nebulosa)
library(enrichR)


# 配色 =========================================================================
# 主要用于UMAP/TSNE图,具体根据聚类数灵活判断
  
# 选择1:高对比度颜色------------------------------------------------------------
# 适用细胞聚类数少的情况 如<15
High_contrast_colors <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488",
                          "#F39B7F", "#8491B4", "#91D1C2", "#DC0000",
                          "#7E6148", "#B09C85", "#E18727", "#20854E",
                          "#7876B1", "#6F99AD", "#FFDC91", "#EE4C97",
                          "#CD534C", "#7AA6DC", "#8C6BB1", "#F0A3FF",
                          "#A0D2DB", "#1C7C54", "#E88873", "#F3D266",
                          "#B39EB5", "#59A14F")

# 直接扩展到50个，足够用
High_contrast_colors <- colorRampPalette(High_contrast_colors)(50)



# 选择2:渐变颜色,清新(只有30个色号)---------------------------------------------
# 1.适用细胞聚类数>25
# 2.或者分亚型展示
Gradation_colors <-  c("#CED7E4","#a4cde1", "#67a4cc", "#277fb8", "#549da3","#184879",
                       "#96cb8f","#8bc96d", "#4dae47", "#5c9e43", 
                       "#b79973", "#E39B96","#f38989","#ec5051", "#e32427",
                       "#ef6a45", "#f9b769", "#f9a341","#f48521", "#ee8e46",
                       "#F7DEDB","#d4a6a8", "#af93c4", "#8660a8","#815e99","#6a3d9a",
                       "#c6b598", "#f6f28f", "#d4a55b", "#b05a28"
)



# 选择3:丰富颜色(只有30个色号)--------------------------------------------------
# 所有情况都几乎适用,万金油
More_colors <- c(
  "#f38989", "#549da3", "#b79973", "#f9b769", "#277fb8",
  "#4dae47", "#af93c4", "#96cb8f", "#f48521", "#a4cde1",
  "#ec5051", "#815e99", "#5c9e43", "#d4a6a8", "#f9a341",
  "#e32427", "#8bc96d", "#67a4cc", "#ee8e46", "#8660a8",
  "#f6f28f", "#b05a28", "#ef6a45", "#d4a55b", "#F0A3FF",
  "#A0D2DB", "#4DBBD5", "#EE4C97", "#1C7C54", "#F3D266")


# 1.数据加载和亚群提取==========================================================
# 导入大群rds数据为sc_integer！！！！,注意哈
sc_integer <- readRDS(file = "input/Final_sc_integer_singleR.rds")


# 提取 AT2 亚群
Idents(sc_integer)
AT2 = subset(sc_integer, idents = 'AT2')
DimPlot(AT2)


# 保存提取细胞的seurat对象
saveRDS(AT2, "output/AT2_subtype_no_type.rds")



# # 2.重新聚类分析================================================================
# cat("基因数量（行）:", nrow(AT2), "\n")
# cat("细胞数量（列）:", ncol(AT2), "\n")
# 
# # 最好补上对count的logx+1标准化，保险
# AT2 = NormalizeData(AT2)
# 
# # 寻找高变基因
# AT2 = FindVariableFeatures(AT2, selection.method = "vst", nfeatures = 2000)  
# 
# # 数据标准化和缩放
# # all.genes = rownames(AT2)
# AT2 = ScaleData(AT2, features = VariableFeatures(AT2))  # 第一轮缩放
# # AT2 = ScaleData(AT2, vars.to.regress = "nFeature_RNA") # 第二轮缩放（技术偏差校正）
# 
# # PCA降维和聚类
# AT2 = RunPCA(AT2,
#              features = VariableFeatures(object = AT2))
# 
# 
# # 去批次
# AT2_integer <- AT2 %>%
#   harmony::RunHarmony(group.by.vars = "sample")
# 
# # 可视化比较
# p1 <- DimPlot(AT2, reduction = "pca", pt.size = .1, group.by = "sample") + ggtitle("Before Harmony")
# p2 <- DimPlot(AT2_integer, reduction = "harmony", pt.size = .1, group.by = "sample") + ggtitle("After Harmony")
# p1+p2
# 

# AT2_integer = FindNeighbors(AT2_integer,reduction = "harmony", dims = 1:20)
# AT2_integer = FindClusters(AT2_integer,reduction = "harmony",resolution = 0.2)  # 较低的分辨率产生较少聚类
# 
# 
# # UMAP降维~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# AT2_integer = RunUMAP(AT2_integer,
#                                  dims = 1:30,
#                                  min.dist = 0.3,#控制点的聚合程度，越小点越紧密，越大点越松散,默认0.3,RunUMAP特有参数
#                                  n.neighbors = 30,#UMAP 计算的邻居点数量，决定局部与全局结构,越大越聚焦全局结构,默认30
#                                  reduction = "harmony")
# 
# # UMAP可视化
# UMAP1 <- DimPlot(AT2_integer,
#                  reduction = 'umap',
#                  label = TRUE,
#                  alpha = 0.4, # 点透明度
#                  cols = More_colors,
#                  pt.size = 0.1)+
#   theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),  # 添加边框
#         axis.line = element_blank())  # 去掉默认轴线 
# UMAP1
# ggsave(filename = "output/1-UMAP.pdf", plot = UMAP1, width = 15, height =13, units = "cm")
# 
# 
# ## 4.2细胞注释重命名------------------------------------------------------------
# # 创建重命名映射~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# new_names <- c(
#   "0" = "Epi C01",
#   "1" = "Epi C02",
#   "2" = "Epi C03",
#   "3" = "Epi C04",
#   "4" = "Epi C05",
#   "5" = "Epi C06",
#   "6" = "Epi C07",
#   "7" = "Epi C08",
#   "8" = "Epi C09",
#   "9" = "Epi C10")
# 
# # 应用重命名~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# cluster_to_celltype_map <- as.list(new_names)
# cluster_to_celltype_map
# 
# 
# # 创建新的细胞类型标识符列~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# cell_types <- as.character(Idents(AT2_integer))
# for (cluster_id in names(cluster_to_celltype_map)) {
#   cell_types[cell_types == cluster_id] <- cluster_to_celltype_map[[cluster_id]]
# }
# 
# # 将细胞类型注释添加到Seurat对象~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# AT2_integer$cell_type <- cell_types
# 
# 
# # 可视化~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# AT2_UMAP <- DimPlot(AT2_integer, 
#                              reduction = "umap", 
#                              group.by = "cell_type",
#                              alpha = 0.4, # 点透明度
#                              cols = More_colors,
#                              label = TRUE, 
#                              pt.size = 0.5) +
#   theme(plot.title = element_blank(),#不显示标题
#         panel.border = element_rect(color = "black", fill = NA, linewidth = 1),  # 添加边框
#         axis.line = element_blank())  # 去掉默认轴线
# AT2_UMAP
# ggsave(filename = "output/2-AT2_UMAP.pdf", plot = AT2_UMAP, width = 16, height = 13, units = "cm")
# 
# 
# 
# # 3 将注释好的AT2亚群细胞标签打回总细胞对象=========================
# # 创建总细胞对象的细胞类型注释列（如果不存在则初始化）
# if (!("cell_subtype" %in% colnames(sc_integer@meta.data))) {
#   sc_integer$cell_subtype <- as.character(sc_integer$cell_anno)
# }
# 
# # 将AT2亚群注释结果映射回总细胞对象
# # 获取AT2细胞的barcodes
# AT2_barcodes <- colnames(AT2_integer)
# 
# # 将AT2亚群注释更新到总细胞对象
# sc_integer$cell_subtype[AT2_barcodes] <- as.character(AT2_integer$cell_type)
# 
# # 查看更新后的细胞类型分布
# table(sc_integer$cell_subtype)
# 
# # 保存更新后的总细胞对象
# saveRDS(sc_integer, "output/sc_integer_with_AT2_subtype.rds")
# 
# 
# # 4.绘制打回注释的UMAP制图======================================================
# UMAP2 <- DimPlot(sc_integer,
#                  group.by="cell_subtype",
#                  label=T,
#                  alpha = 0.4, # 点透明度
#                  pt.size = 0.05,
#                  cols = More_colors,
#                  reduction='umap',
#                  ) +  
#   labs(title = NULL)+
#   theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1), # 添加边框
#         axis.line = element_blank(),                                            # 去掉默认轴线
#         axis.ticks = element_blank(),                                           # 去掉刻度线
#         axis.text = element_blank()                                             # 去掉刻度标签
#   )
# UMAP2
# ggsave(filename = "output/3-cell_subtype_Anno.pdf", plot = UMAP2, width = 22, height =14, units = "cm")
# 
# 
# 
# # 5.上皮细胞分面UMAP图==========================================================
# # 按tissue分组显示
# UMAP_split1 <- DimPlot(AT2_integer,
#                        group.by = "cell_type",
#                        split.by = "tissue",
#                        label = TRUE,
#                        cols = More_colors,
#                        # alpha = 0.4,# 点透明度
#                        pt.size = 0.05,
#                        reduction = 'umap'
# ) +
#   labs(title = NULL)+
#   theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1), # 添加边框
#         axis.line = element_blank(),                                            # 去掉默认轴线
#         axis.ticks = element_blank(),                                           # 去掉刻度线
#         axis.text = element_blank()                                             # 去掉刻度标签
#   )
# UMAP_split1
# 
# ggsave(filename = "output/4-cell_subtype_Anno_split.pdf", plot = UMAP_split1, width = 25, height =13, units = "cm")

