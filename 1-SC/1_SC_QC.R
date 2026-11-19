# -----------------------------------------------------------------------------#
# 描述: 本脚本用于执行单细胞转录组分析 1_SC_QC.R
# R版本：4.3.2
# Seurat (Version:5.0.1)
# 日期: 2026/06/25
# 创建：阿伦
# 代码版本：V2
# 
# 更新内容：
# 2026.6.24 更新 单独样本的质控,MAD质控
# 2026.6.25 更新 decontX 去除环境RNA
# 2026.7.4 更新 3种可视化配色
# -----------------------------------------------------------------------------#


# -----------------------------------------------------------------------------#
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<< 此项目为人类 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#
# -----------------------------------------------------------------------------#


# 哈佛单细胞教程
# https://hbctraining.github.io/scRNA-seq/lessons/06_SC_SCT_and_integration.html


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


# 添加镜像======================================================================
# 添加bioconductor清华镜像
options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
# 添加bioconductor西湖大学镜像
options(BioC_mirror = "https://mirrors.westlake.edu.cn/bioconductor")
# 添加R清华镜像
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))


# 加载R包=======================================================================
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
library(homologene)
library(DoubletFinder)
library(scDblFinder)
library(decontX)

# 首先根据自己要做的不同类型原始单细胞数据的格式选择读取的方式 


# # --------------------------------解压参考-------------------------------------# 
# R.utils::gunzip("raw/GSM8209057_extra-pulmonary.csv.gz", 
#                 destname = "raw/GSM8209057_extra-pulmonary.csv", 
#                 overwrite = FALSE,
#                 remove = FALSE)  # 保留原始压缩文件
# # --------------------------------解压参考-------------------------------------# 



# 定义变量======================================================================
# 定义
treat_name = "BPD"       # 疾病组命名
con_name = "Control"     # 对照组命名



# -----------------------------------------------------------------------------#
# <<<<<<<<<<<<<<<<<<<<< 一.标准 10X形式的单细胞文件操作>>>>>>>>>>>>>>>>>>>>>>>>#
# -----------------------------------------------------------------------------#


# 1.读入数据====================================================================
# 1.1 批量10X形式数据读取与构建Seurat对象---------------------------------------

# 先手动给每个样本创建对应GSM编号的文件夹，将10X的3个文件丢进去。
# 别问为什么不用代码批量读取，不许问！！！！

# 例子：单个文件结构
#  - GSM6613031
#      -- barcodes.tsv.gz
#      -- features.tsv.gz     # 有的人命名 genes.tsv.gz 得改成features.tsv.gz,不然读取不了
#      -- matrix.mtx.gz

# Read10X~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
filenames <- dir("raw/GSE217488/", full.names = T) # 10X文件地址
sc_list <- lapply(filenames, function(x) Read10X(x, gene.column = 2))


# 为每个样本构建Seurat对象~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
seuratobject_data <- lapply(sc_list, function(x){
  CreateSeuratObject(counts = x,
                     min.cells = 3,      #保留至少在3个细胞中表达的基因
                     min.features = 200) #保留至少具有200个基因的细胞
})


# 2 给元数据赋标签 GEOID / sample / tissue =====================================
# 2.1 给每个Seurat对象添加样本标签(GEOID)---------------------------------------
GEOID <- stringr::str_extract(filenames, "GSM\\d+")

for(i in 1:length(GEOID)){                    # i为样本的长度，把第 i 个 Seurat 对象的 sample 列设为第 i 个样本名
  seuratobject_data[[i]]$GEOID <- GEOID[i]    # 为Seurat对象添加样本标签
}

# GSM6720374	Hyperoxia   高氧 (Hyperoxia) 疾病
# GSM6720373	Room Air    常氧 (Room Air)  对照



# 2.2 给每个Seurat对象添加标签(sample/tissue)-----------------------------------
# 手动定义每个 GSM 对应的样本名(疾病+对照都可以识别)
sample_mapping <- c(
  # 疾病样本
  GSM6720374 = paste0(treat_name, "1"),
# 对照样本(没有对照就注释掉)
  GSM6720373 = paste0(con_name, "1")
)

# --------------------------------示例-----------------------------------------#
# GSM8029962 = paste0(treat_name, "1"),  看作 GSM8029962 = AS1
# GSM8029962 = paste0(treat_name, "2"),
# GSM8029963 = paste0(con_name, "1"),
# GSM8029963 = paste0(con_name, "2"),
# -----------------------------------------------------------------------------#

# 为每个 Seurat 对象添加 sample 和 tissue 列
for(i in 1:length(GEOID)) {
  gsm_id <- GEOID[i]
  samp_name <- sample_mapping[[gsm_id]]

  seuratobject_data[[i]]$sample <- samp_name

  # 自动判断：如果 sample 名包含 con_name，则为对照，否则为处理组
  if (grepl(con_name, samp_name, fixed = TRUE)) {
    seuratobject_data[[i]]$tissue <- con_name
  } else {
    seuratobject_data[[i]]$tissue <- treat_name
  }
}



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~检查样本是否匹配~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
data.frame(
  filename = basename(filenames),
  GEOID = GEOID,
  sample = sample_mapping[GEOID],
  tissue = sapply(seuratobject_data, function(x) x$tissue[1])
)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~检查样本是否匹配~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




# -----------------------------------------------------------------------------#
# <<<<<<<<<<<<<<<<<<<<<<< 二.5h形式的单细胞文件操作>>>>>>>>>>>>>>>>>>>>>>>>>>>>#
# -----------------------------------------------------------------------------#


# # 1. 构建Seurat对象处理=========================================================
# # 1.1 5h数据构建Seurat对象------------------------------------------------------
# sceList = lapply(samples,function(pro){
#   print(pro)  # 打印当前处理的文件名
# 
#   # 创建Seurat对象
#   sce = CreateSeuratObject(
#     counts = Read10X_h5(file.path(dir,pro)),  # 读取10X Genomics的h5文件
#     project = gsub('_filtered_feature_bc_matrix.h5','',
#                    gsub('^GSM[0-9]*_','',pro)),  # 从文件名提取项目名
#     min.cells = 3,      # 基因至少在3个细胞中表达
#     min.features = 200  # 细胞至少检测到200个基因
#   )
#   return(sce)
# })
# 
# # 如果h5格式的单细胞数据有很多样本用以下的方式合并~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# sce.all=merge(x=sceList[[1]],
# y=sceList[ -1 ],
# add.cell.ids =   gsub('_filtered_feature_bc_matrix.h5','',gsub('^GSM[0-9]*_','',samples)))
# 
# 
# # 从list获得seurat对象：简洁的方式，适用于任意数量的样本
# GSM4998196 <- Reduce(function(x, y) merge(x, y), sceList)
# 
# 
# # 1.2 分别为Seurat对象添加样本标签----------------------------------------------
# # 标注疾病和对照样本(实验分组信息)
# GSM4998196[["tissue"]] <- "SCLC"
# 
# # 标注GEOID
# GSM4998196[["GEOID"]] <- "GSM4998196"



# -----------------------------------------------------------------------------#
# <<<<<<<<<<<<<< 三.txt / .csv /复杂形式矩阵形式的单细胞文件操作>>>>>>>>>>>>>>>#
# -----------------------------------------------------------------------------#


# # 1. 第一种样式 非10X的csv形式数据==============================================
# 
# # -----------------------------------------------------------------------------#
# # 样本
# # 2例饮食诱导造模小鼠
# # 4例正常
# # -----------------------------------------------------------------------------#
# 
# 
# # 读取原始的csv (矩阵)文件~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # 要注意不同人上传的矩阵文件的格式可能会有不同！！！！！！！！，要保证行名为基因名
# 
# 
# # 疾病组~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # 读取矩阵文件(制表符分隔的文件)
# # 先解压文件
# R.utils::gunzip("raw/GSM8209057_extra-pulmonary.csv.gz",
#                 destname = "raw/GSM8209057_extra-pulmonary.csv",
#                 overwrite = FALSE,
#                 remove = FALSE)  # 保留原始压缩文件
# # 再读取
# GSM8209057 <- read.csv("raw/GSM8209057_extra-pulmonary.csv",
#                        header = TRUE,
#                        row.names = 1,
#                        check.names = FALSE)
# 
# # 对照组~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # 读取矩阵文件(制表符分隔的文件)
# # 先解压文件
# R.utils::gunzip("raw/GSM8209058_control.csv.gz",
#                 destname = "raw/GSM8209058_control.csv",
#                 overwrite = FALSE,
#                 remove = FALSE)  # 保留原始压缩文件
# # 再读取
# GSM8209058 <- read.csv("raw/GSM8209058_control.csv",
#                        header = TRUE,
#                        row.names = 1,
#                        check.names = FALSE)



# # 2. 第二种样式 非10X的每个样本单独的.txt形式数据===============================
# # 2.1 数据读取------------------------------------------------------------------
# # 要注意不同人上传的矩阵文件的格式可能会有不同！！！！！！！！，要保证行名为基因名
# # 低速读取法
# # 2疾病组~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# GSM4795788 <- read.delim("raw/GSE158241/GSM4795788_DIO-NASH1.txt",
#                          header = TRUE,
#                          sep = "\t",
#                          check.names = FALSE)
# 
# GSM4795789 <- read.delim("raw/GSE158241/GSM4795789_DIO-NASH2.txt",
#                          header = TRUE,
#                          sep = "\t",
#                          check.names = FALSE)
# 
# 
# # 4对照组~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# GSM4795784 <- read.delim("raw/GSE158241/GSM4795784_Lean-chow1.txt",
#                          header = TRUE,
#                          sep = "\t",
#                          check.names = FALSE)
# 
# GSM4795785 <- read.delim("raw/GSE158241/GSM4795785_Lean-chow2.txt",
#                          header = TRUE,
#                          sep = "\t",
#                          check.names = FALSE)
# 
# GSM4795786 <- read.delim("raw/GSE158241/GSM4795786_Lean-chow3.txt",
#                          header = TRUE,
#                          sep = "\t",
#                          check.names = FALSE)
# 
# GSM4795787 <- read.delim("raw/GSE158241/GSM4795787_Lean-chow4.txt",
#                          header = TRUE,
#                          sep = "\t",
#                          check.names = FALSE)


# # 2.2 构建Seurat对象(分别创建6个样本的Seurat对象)-------------------------------
# # 疾病~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# seurat_GSM4795788 <- CreateSeuratObject(counts = GSM4795788,
#                                         min.cells = 3, # 保留至少在3个细胞中表达的基因
#                                         min.features = 200)# 保留至少具有200个基因的细胞
# 
# seurat_GSM4795789 <- CreateSeuratObject(counts = GSM4795789,
#                                         min.cells = 3, # 保留至少在3个细胞中表达的基因
#                                         min.features = 200)# 保留至少具有200个基因的细胞
# 
# # 对照~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# seurat_GSM4795784 <- CreateSeuratObject(counts = GSM4795784,
#                                         min.cells = 3, # 保留至少在3个细胞中表达的基因
#                                         min.features = 200)# 保留至少具有200个基因的细胞
# 
# seurat_GSM4795785 <- CreateSeuratObject(counts = GSM4795785,
#                                         min.cells = 3, # 保留至少在3个细胞中表达的基因
#                                         min.features = 200)# 保留至少具有200个基因的细胞
# 
# seurat_GSM4795786 <- CreateSeuratObject(counts = GSM4795786,
#                                         min.cells = 3, # 保留至少在3个细胞中表达的基因
#                                         min.features = 200)# 保留至少具有200个基因的细胞
# 
# seurat_GSM4795787 <- CreateSeuratObject(counts = GSM4795787,
#                                         min.cells = 3, # 保留至少在3个细胞中表达的基因
#                                         min.features = 200)# 保留至少具有200个基因的细胞



# # 2.3 分别为Seurat对象添加样本标签sample和tissu和GEOID--------------------------
# # 添加sample~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# seurat_GSM4795788[["sample"]] <- "AS1"
# seurat_GSM4795789[["sample"]] <- "AS2"
# seurat_GSM4795784[["sample"]] <- "Control1"
# seurat_GSM4795785[["sample"]] <- "Control2"
# seurat_GSM4795786[["sample"]] <- "Control3"
# seurat_GSM4795787[["sample"]] <- "Control4"
# 
# 
# 
# # 标注疾病和对照样本(实验分组信息)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# seurat_GSM4795788[["tissue"]] <- "AS"
# seurat_GSM4795789[["tissue"]] <- "AS"
# seurat_GSM4795784[["tissue"]] <- "Control"
# seurat_GSM4795785[["tissue"]] <- "Control"
# seurat_GSM4795786[["tissue"]] <- "Control"
# seurat_GSM4795787[["tissue"]] <- "Control"
# 
# 
# # 标注GEOID~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# seurat_GSM4795788[["GEOID"]] <- "GSM4795788 AS1"
# seurat_GSM4795789[["GEOID"]] <- "GSM4795789 AS2"
# seurat_GSM4795784[["GEOID"]] <- "GSM4795784 Control1"
# seurat_GSM4795785[["GEOID"]] <- "GSM4795785 Control2"
# seurat_GSM4795786[["GEOID"]] <- "GSM4795786 Control3"
# seurat_GSM4795787[["GEOID"]] <- "GSM4795787 Control4"
# 
# 
# # 2.4合并Seurat对象-------------------------------------------------------------
# #创建一个包含所有Seurat对象的列表
# seuratobject_data <- list(seurat_GSM4795788,
#                           seurat_GSM4795789,
#                           seurat_GSM4795784,
#                           seurat_GSM4795785,
#                           seurat_GSM4795786,
#                           seurat_GSM4795787
# )



# # 3. 第三种样式 Matrix Market数据形式数据=======================================
# # 3.1 读取Matrix Market格式的稀疏矩阵文件---------------------------------------
# # 使用Matrix包读取.mtx文件
# courtine_lumbar_EES_rehab <- Matrix::readMM("raw/courtine_lumbar-EES-rehab_UMI.mtx")
# 
# # # 解压并读取基因名和细胞barcode文件
# # # 解压features文件（基因名）
# # R.utils::gunzip("raw/GSE184370_features.txt.gz",
# #                 destname = "raw/GSE184370_features.txt",
# #                 overwrite = TRUE,
# #                 remove = FALSE)
# # # 解压barcodes文件
# # R.utils::gunzip("raw/GSE184370_barcodes.txt.gz",
# #                 destname = "raw/GSE184370_barcodes.txt",
# #                 overwrite = TRUE,
# #                 remove = FALSE)
# 
# # 读取基因名和barcodes
# features <- read.delim("raw/GSE184370_features.txt", header = FALSE, stringsAsFactors = FALSE)
# barcodes <- read.delim("raw/GSE184370_barcodes.txt", header = FALSE, stringsAsFactors = FALSE)
# 
# # 设置行名（基因名）和列名（细胞barcode）
# # 注意：Matrix Market格式的矩阵通常是基因×细胞，需要转置为细胞×基因
# # 但Seurat期望的是基因×细胞的矩阵，所以不需要转置
# rownames(courtine_lumbar_EES_rehab) <- features$V1
# colnames(courtine_lumbar_EES_rehab) <- barcodes$V1



# # 3.2 样本提取功能--------------------------------------------------------------
# # 读取meta.txt文件（细胞注释信息）
# meta_data <- read.delim("raw/meta.txt", header = TRUE, stringsAsFactors = FALSE)
# cat("Meta数据维度:", dim(meta_data), "\n")
# 
# # 读取MAP.csv文件，获取GSM编号与样本标识映射关系
# map_data <- read.csv("raw/MAP.csv", header = TRUE, stringsAsFactors = FALSE)
# cat("MAP.csv文件内容:\n")
# print(map_data)
# 
# # 从barcode中提取样本标识（例如从"A-45:AAAGTAGCAGCTGCACx"中提取"A-45"）
# sample_ids <- sapply(strsplit(colnames(courtine_lumbar_EES_rehab), ":"), "[[", 1)
# 
# # 查看所有样本标识
# cat("所有样本标识:", unique(sample_ids), "\n")
# cat("样本分布:\n")
# print(table(sample_ids))


# # 3.3 构建Seurat对象-------------------------------------------------------------
# # courtine_lumbar-EES-rehab数据构建Seurat对象
# 
# # 为每个GSM样本创建独立的Seurat对象
# cat("=== 为每个GSM样本创建独立Seurat对象 ===\n")
# 
# # 根据MAP.csv文件，为每个GSM样本提取对应的细胞并创建Seurat对象
# for (i in 1:nrow(map_data)) {
#   gsm_id <- map_data$GSM[i]
#   sample_id <- map_data$name[i]
# 
#   # 提取当前样本的细胞
#   sample_cells <- which(sample_ids == sample_id)
# 
#   if (length(sample_cells) > 0) {
#     # 提取子矩阵
#     sample_matrix <- courtine_lumbar_EES_rehab[, sample_cells]
# 
#     # 创建Seurat对象，使用GSM编号作为对象名
#     assign(gsm_id, CreateSeuratObject(counts = sample_matrix,
#                                       min.cells = 3,      # 保留至少在3个细胞中表达的基因
#                                       min.features = 200, # 保留至少具有200个基因的细胞
#                                       project = gsm_id))
# 
#     cat("成功创建", gsm_id, "Seurat对象，细胞数量:", ncol(get(gsm_id)), "\n")
#   } else {
#     cat("警告：样本", sample_id, "未找到对应的细胞\n")
#   }
# }
# 
# cat("所有GSM样本Seurat对象创建完成！\n")


# # 3.4 给元数据赋标签 GEOID / sample / tissue---------------------------------------
# # 给每个Seurat对象添加样本标签(GEOID)
# GSM5585216[["GEOID"]] <- "GSM5585216 SCI"
# GSM5585217[["GEOID"]] <- "GSM5585217 SCI"
# GSM5585218[["GEOID"]] <- "GSM5585218 SCI"
# GSM5585219[["GEOID"]] <- "GSM5585219 Control"
# GSM5585220[["GEOID"]] <- "GSM5585220 Control"
# GSM5585221[["GEOID"]] <- "GSM5585221 Control"
# 
# 
# # 给每个Seurat对象添加标签(sample)
# GSM5585216[["sample"]] <- "SCI 1"
# GSM5585217[["sample"]] <- "SCI 2"
# GSM5585218[["sample"]] <- "SCI 3"
# GSM5585219[["sample"]] <- "Control 1"
# GSM5585220[["sample"]] <- "Control 2"
# GSM5585221[["sample"]] <- "Control 3"
# 
# 
# # 给每个Seurat对象添加标签(tissue)
# GSM5585216[["tissue"]] <- "SCI"
# GSM5585217[["tissue"]] <- "SCI"
# GSM5585218[["tissue"]] <- "SCI"
# GSM5585219[["tissue"]] <- "Control"
# GSM5585220[["tissue"]] <- "Control"
# GSM5585221[["tissue"]] <- "Control"
# 
# 
# # 3.5 合并所有单个矩阵为矩阵list（seuratobject_data）---------------------------
# seuratobject_data <- list(GSM5585216,GSM5585217,GSM5585218,GSM5585219,GSM5585220,GSM5585221)


# # 4. 第四种样式 非10X的总.txt表形式数据=========================================
# # 4.1 读取总表达矩阵文件--------------------------------------------------------
# expr_matrix <- read.table("raw/GSE149614_HCC.scRNAseq.S71915.count.txt.gz", 
#                           header = TRUE, 
#                           row.names = 1,
#                           check.names = FALSE,
#                           stringsAsFactors = FALSE)
# 
# # expr_matrix <- read.delim("raw/GSE149614_HCC.scRNAseq.S71915.count.txt.gz",
# #                            header = TRUE,
# #                            sep = "\t",
# #                            check.names = FALSE)
# 
# 
# 
# # # 转换为稀疏矩阵以节省内存
# # expr_matrix <- as.matrix(expr_matrix)
# # expr_matrix <- Matrix::Matrix(expr_matrix, sparse = TRUE)
# 
# 
# 
# # 4.2 从barcode中提取样本信息并筛选PT和NTL样本--------------------------------
# # 从barcode中提取样本标识（例如从"HCC01T_AAACCTGAGGGCATGT"中提取"HCC01T"）
# sample_ids <- sapply(strsplit(colnames(expr_matrix), "_"), "[[", 1)
# 
# # 需要的样本：PT (HCC01T-HCC10T) 和 NTL (HCC03N-HCC10N)
# selected_samples <- c(
#   "HCC01T", "HCC02T", "HCC03T", "HCC04T", "HCC05T",
#   "HCC06T", "HCC07T", "HCC08T", "HCC09T", "HCC10T",  # PT
#   "HCC03N", "HCC04N", "HCC05N", "HCC06N",             # NTL
#   "HCC07N", "HCC08N", "HCC09N", "HCC10N"              # NTL
# )
# 
# # 筛选只包含PT和NTL的细胞
# cells_to_keep <- which(sample_ids %in% selected_samples)
# expr_matrix <- expr_matrix[, cells_to_keep]
# sample_ids <- sample_ids[cells_to_keep]
# 
# 
# 
# # 4.3 为每个样本创建独立的Seurat对象-------------------------------------------
# # 定义样本分组信息
# # PT (Primary Tumor, 原发肿瘤) -> HCC
# # NTL (Non-Tumor Liver, 非肿瘤肝组织) -> Control
# sample_groups <- c(
#   HCC01T = "HCC", HCC02T = "HCC", HCC03T = "HCC",
#   HCC04T = "HCC", HCC05T = "HCC", HCC06T = "HCC",
#   HCC07T = "HCC", HCC08T = "HCC", HCC09T = "HCC", HCC10T = "HCC",
#   HCC03N = "Control", HCC04N = "Control", HCC05N = "Control",
#   HCC06N = "Control", HCC07N = "Control", HCC08N = "Control",
#   HCC09N = "Control", HCC10N = "Control"
# )
# 
# # 定义sample列的命名（HCC 1-10, Control 1-8）
# sample_names <- c(
#   HCC01T = "HCC 1", HCC02T = "HCC 2", HCC03T = "HCC 3",
#   HCC04T = "HCC 4", HCC05T = "HCC 5", HCC06T = "HCC 6",
#   HCC07T = "HCC 7", HCC08T = "HCC 8", HCC09T = "HCC 9", HCC10T = "HCC 10",
#   HCC03N = "Control 1", HCC04N = "Control 2", HCC05N = "Control 3",
#   HCC06N = "Control 4", HCC07N = "Control 5", HCC08N = "Control 6",
#   HCC09N = "Control 7", HCC10N = "Control 8"
# )
# 
# # 循环为每个样本创建Seurat对象
# seuratobject_data <- list()
# unique_samples <- unique(sample_ids)
# 
# for (sample_id in unique_samples) {
#   cat("处理样本:", sample_id, "\n")
#   
#   # 提取当前样本的细胞
#   sample_cells <- which(sample_ids == sample_id)
#   
#   if (length(sample_cells) > 0) {
#     # 提取子矩阵
#     sample_matrix <- expr_matrix[, sample_cells]
#     
#     # 注意：原始列名已包含样本前缀（如HCC01T_AAACCTGAGGGCATGT），无需修改
#     
#     # 创建Seurat对象
#     seurat_obj <- CreateSeuratObject(
#       counts = sample_matrix,
#       min.cells = 3,      # 保留至少在3个细胞中表达的基因
#       min.features = 200, # 保留至少具有200个基因的细胞
#       project = sample_id
#     )
#     
#     # 添加元数据
#     seurat_obj[["GEOID"]] <- sample_id
#     seurat_obj[["sample"]] <- sample_names[[sample_id]]
#     seurat_obj[["tissue"]] <- sample_groups[[sample_id]]
#     
#     # 存储到列表
#     seuratobject_data[[sample_id]] <- seurat_obj
#     
#     cat("  成功创建", sample_id, "Seurat对象，细胞数量:", ncol(seurat_obj), "\n")
#   }
# }
# 
# cat("\nHCC和Control样本Seurat对象创建完成！共", length(seuratobject_data), "个样本\n")


# # 4. 第四种样式 非10X的总.txt表形式数据=========================================
# # 4.1 读取总表达矩阵文件--------------------------------------------------------
# GSE131907 <- readRDS("raw/GSE131907_Lung_Cancer_raw_UMI_matrix.rds")
# 
# 
# # 我需要的样本信息供查阅--------------------------------------------------------#
# # 全部注释
# # nLung，正常肺；     11个
# # tLung,早期肿瘤肺；  11个
# # tL/B：晚期肿瘤肺；  4个
# # mLN，淋巴结转移；
# # nLN，正常淋巴结；
# # PE，胸腔积液；      5个
# # mBrain，脑转移；
# 
# # 11个tLung早期肿瘤肺,11个nLung正常肺,共22个~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # GEO ID	 SamplesTissue origins
# # 11个nLung正常肺
# # GSM3827114	LUNG_N01  nLung
# # GSM3827115	LUNG_N06  nLung
# # GSM3827116	LUNG_N08  nLung
# # GSM3827117	LUNG_N09  nLung
# # GSM3827118	LUNG_N18  nLung
# # GSM3827119	LUNG_N19  nLung
# # GSM3827120	LUNG_N20  nLung
# # GSM3827121	LUNG_N28  nLung
# # GSM3827122	LUNG_N30  nLung
# # GSM3827123	LUNG_N31  nLung
# # GSM3827124	LUNG_N34  nLung
# # 11个tLung早期肿瘤肺
# # GSM3827125	LUNG_T06	tLung
# # GSM3827126	LUNG_T08	tLung
# # GSM3827127	LUNG_T09	tLung
# # GSM3827128	LUNG_T18	tLung
# # GSM3827129	LUNG_T19	tLung
# # GSM3827130	LUNG_T20	tLung
# # GSM3827131	LUNG_T25	tLung
# # GSM3827132	LUNG_T28	tLung
# # GSM3827133	LUNG_T30	tLung
# # GSM3827134	LUNG_T31	tLung
# # GSM3827135	LUNG_T34	tLung
# #------------------------------------------------------------------------------#
# 
# 
# # 4.2 筛选项目要用的样本--------------------------------------------------------
# # 筛选各个样本的前缀对应的列~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # 11个nLung正常肺
# LUNG_N01_indices <- grep("LUNG_N01$", colnames(GSE131907))
# LUNG_N06_indices <- grep("LUNG_N06$", colnames(GSE131907))
# LUNG_N08_indices <- grep("LUNG_N08$", colnames(GSE131907))
# LUNG_N09_indices <- grep("LUNG_N09$", colnames(GSE131907))
# LUNG_N18_indices <- grep("LUNG_N18$", colnames(GSE131907))
# LUNG_N19_indices <- grep("LUNG_N19$", colnames(GSE131907))
# LUNG_N20_indices <- grep("LUNG_N20$", colnames(GSE131907))
# LUNG_N28_indices <- grep("LUNG_N28$", colnames(GSE131907))
# LUNG_N30_indices <- grep("LUNG_N30$", colnames(GSE131907))
# LUNG_N31_indices <- grep("LUNG_N31$", colnames(GSE131907))
# LUNG_N34_indices <- grep("LUNG_N34$", colnames(GSE131907))
# # 11个tLung早期肿瘤肺
# LUNG_T06_indices <- grep("LUNG_T06$", colnames(GSE131907)) 
# LUNG_T08_indices <- grep("LUNG_T08$", colnames(GSE131907))
# LUNG_T09_indices <- grep("LUNG_T09$", colnames(GSE131907))
# LUNG_T18_indices <- grep("LUNG_T18$", colnames(GSE131907))
# LUNG_T19_indices <- grep("LUNG_T19$", colnames(GSE131907))
# LUNG_T20_indices <- grep("LUNG_T20$", colnames(GSE131907))
# LUNG_T25_indices <- grep("LUNG_T25$", colnames(GSE131907))
# LUNG_T28_indices <- grep("LUNG_T28$", colnames(GSE131907))
# LUNG_T30_indices <- grep("LUNG_T30$", colnames(GSE131907))
# LUNG_T31_indices <- grep("LUNG_T31$", colnames(GSE131907))
# LUNG_T34_indices <- grep("LUNG_T34$", colnames(GSE131907))
# 
# 
# # 提取对应样本的数据并改名~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # 11个nLung正常肺
# GSM3827114 <- GSE131907[, LUNG_N01_indices]
# GSM3827115 <- GSE131907[, LUNG_N06_indices]
# GSM3827116 <- GSE131907[, LUNG_N08_indices]
# GSM3827117 <- GSE131907[, LUNG_N09_indices]
# GSM3827118 <- GSE131907[, LUNG_N18_indices]
# GSM3827119 <- GSE131907[, LUNG_N19_indices]
# GSM3827120 <- GSE131907[, LUNG_N20_indices]
# GSM3827121 <- GSE131907[, LUNG_N28_indices]
# GSM3827122 <- GSE131907[, LUNG_N30_indices]
# GSM3827123 <- GSE131907[, LUNG_N31_indices]
# GSM3827124 <- GSE131907[, LUNG_N34_indices]
# # 11个tLung早期肿瘤肺
# GSM3827125 <- GSE131907[, LUNG_T06_indices]
# GSM3827126 <- GSE131907[, LUNG_T08_indices]
# GSM3827127 <- GSE131907[, LUNG_T09_indices]
# GSM3827128 <- GSE131907[, LUNG_T18_indices]
# GSM3827129 <- GSE131907[, LUNG_T19_indices]
# GSM3827130 <- GSE131907[, LUNG_T20_indices]
# GSM3827131 <- GSE131907[, LUNG_T25_indices]
# GSM3827132 <- GSE131907[, LUNG_T28_indices]
# GSM3827133 <- GSE131907[, LUNG_T30_indices]
# GSM3827134 <- GSE131907[, LUNG_T31_indices]
# GSM3827135 <- GSE131907[, LUNG_T34_indices]
# 
# 
# 
# # 4.3 构建Seurat对象(分别创建22个样本的Seurat对象)------------------------------
# # 样本list
# sample_names <- c("GSM3827114", "GSM3827115", "GSM3827116", "GSM3827117", "GSM3827118",
#                   "GSM3827119", "GSM3827120", "GSM3827121", "GSM3827122", "GSM3827123",
#                   "GSM3827124", "GSM3827125", "GSM3827126", "GSM3827127", "GSM3827128",
#                   "GSM3827129", "GSM3827130", "GSM3827131", "GSM3827132", "GSM3827133",
#                   "GSM3827134", "GSM3827135")
# 
# 
# # 创建Seurat对象列表
# seurat_list <- lapply(sample_names, function(sample) {
#   CreateSeuratObject(counts = get(sample),
#                      min.cells = 3,      #保留至少在3个细胞中表达的基因
#                      min.features = 200) #保留至少具有200个基因的细胞
# })
# 
# 
# # 为列表中的每个对象命名
# names(seurat_list) <- paste0("seurat_", sample_names)
# 
# 
# # 将列表中的对象分配到全局环境
# list2env(seurat_list, envir = .GlobalEnv)
# 
# 
# # 4.4 分别为Seurat对象添加样本标签sample、tissue和GEOID--------------------------
# # 添加sample~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # 11个正常肺（Control组）
# seurat_GSM3827114[["sample"]] <- paste0(con_name, "01")
# seurat_GSM3827115[["sample"]] <- paste0(con_name, "02")
# seurat_GSM3827116[["sample"]] <- paste0(con_name, "03")
# seurat_GSM3827117[["sample"]] <- paste0(con_name, "04")
# seurat_GSM3827118[["sample"]] <- paste0(con_name, "05")
# seurat_GSM3827119[["sample"]] <- paste0(con_name, "06")
# seurat_GSM3827120[["sample"]] <- paste0(con_name, "07")
# seurat_GSM3827121[["sample"]] <- paste0(con_name, "08")
# seurat_GSM3827122[["sample"]] <- paste0(con_name, "09")
# seurat_GSM3827123[["sample"]] <- paste0(con_name, 10)
# seurat_GSM3827124[["sample"]] <- paste0(con_name, 11)
# # 11个肿瘤肺（NSCLC组）
# seurat_GSM3827125[["sample"]] <- paste0(treat_name, "01")
# seurat_GSM3827126[["sample"]] <- paste0(treat_name, "02")
# seurat_GSM3827127[["sample"]] <- paste0(treat_name, "03")
# seurat_GSM3827128[["sample"]] <- paste0(treat_name, "04")
# seurat_GSM3827129[["sample"]] <- paste0(treat_name, "05")
# seurat_GSM3827130[["sample"]] <- paste0(treat_name, "06")
# seurat_GSM3827131[["sample"]] <- paste0(treat_name, "07")
# seurat_GSM3827132[["sample"]] <- paste0(treat_name, "08")
# seurat_GSM3827133[["sample"]] <- paste0(treat_name, "09")
# seurat_GSM3827134[["sample"]] <- paste0(treat_name, 10)
# seurat_GSM3827135[["sample"]] <- paste0(treat_name, 11)
# 
# 
# # 标注疾病和对照样本(实验分组信息)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # 11个正常肺（Control组）
# seurat_GSM3827114[["tissue"]] <- con_name
# seurat_GSM3827115[["tissue"]] <- con_name
# seurat_GSM3827116[["tissue"]] <- con_name
# seurat_GSM3827117[["tissue"]] <- con_name
# seurat_GSM3827118[["tissue"]] <- con_name
# seurat_GSM3827119[["tissue"]] <- con_name
# seurat_GSM3827120[["tissue"]] <- con_name
# seurat_GSM3827121[["tissue"]] <- con_name
# seurat_GSM3827122[["tissue"]] <- con_name
# seurat_GSM3827123[["tissue"]] <- con_name
# seurat_GSM3827124[["tissue"]] <- con_name
# # 11个肿瘤肺（NSCLC组）
# seurat_GSM3827125[["tissue"]] <- treat_name
# seurat_GSM3827126[["tissue"]] <- treat_name
# seurat_GSM3827127[["tissue"]] <- treat_name
# seurat_GSM3827128[["tissue"]] <- treat_name
# seurat_GSM3827129[["tissue"]] <- treat_name
# seurat_GSM3827130[["tissue"]] <- treat_name
# seurat_GSM3827131[["tissue"]] <- treat_name
# seurat_GSM3827132[["tissue"]] <- treat_name
# seurat_GSM3827133[["tissue"]] <- treat_name
# seurat_GSM3827134[["tissue"]] <- treat_name
# seurat_GSM3827135[["tissue"]] <- treat_name
# 
# 
# # 标注GEOID~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # 11个正常肺（Control组）
# seurat_GSM3827114[["GEOID"]] <- paste0("GSM3827114 ", con_name, "01")
# seurat_GSM3827115[["GEOID"]] <- paste0("GSM3827115 ", con_name, "02")
# seurat_GSM3827116[["GEOID"]] <- paste0("GSM3827116 ", con_name, "03")
# seurat_GSM3827117[["GEOID"]] <- paste0("GSM3827117 ", con_name, "04")
# seurat_GSM3827118[["GEOID"]] <- paste0("GSM3827118 ", con_name, "05")
# seurat_GSM3827119[["GEOID"]] <- paste0("GSM3827119 ", con_name, "06")
# seurat_GSM3827120[["GEOID"]] <- paste0("GSM3827120 ", con_name, "07")
# seurat_GSM3827121[["GEOID"]] <- paste0("GSM3827121 ", con_name, "08")
# seurat_GSM3827122[["GEOID"]] <- paste0("GSM3827122 ", con_name, "09")
# seurat_GSM3827123[["GEOID"]] <- paste0("GSM3827123 ", con_name, 10)
# seurat_GSM3827124[["GEOID"]] <- paste0("GSM3827124 ", con_name, 11)
# # 11个肿瘤肺（NSCLC组）
# seurat_GSM3827125[["GEOID"]] <- paste0("GSM3827125 ", treat_name, "01")
# seurat_GSM3827126[["GEOID"]] <- paste0("GSM3827126 ", treat_name, "02")
# seurat_GSM3827127[["GEOID"]] <- paste0("GSM3827127 ", treat_name, "03")
# seurat_GSM3827128[["GEOID"]] <- paste0("GSM3827128 ", treat_name, "04")
# seurat_GSM3827129[["GEOID"]] <- paste0("GSM3827129 ", treat_name, "05")
# seurat_GSM3827130[["GEOID"]] <- paste0("GSM3827130 ", treat_name, "06")
# seurat_GSM3827131[["GEOID"]] <- paste0("GSM3827131 ", treat_name, "07")
# seurat_GSM3827132[["GEOID"]] <- paste0("GSM3827132 ", treat_name, "08")
# seurat_GSM3827133[["GEOID"]] <- paste0("GSM3827133 ", treat_name, "09")
# seurat_GSM3827134[["GEOID"]] <- paste0("GSM3827134 ", treat_name, 10)
# seurat_GSM3827135[["GEOID"]] <- paste0("GSM3827135 ", treat_name, 11)
# 
# 
# 
# # 2.3 合并所有单个矩阵为矩阵list（seuratobject_data）----------------------------
# seuratobject_data <- list(seurat_GSM3827114, seurat_GSM3827115, seurat_GSM3827116,
#                           seurat_GSM3827117, seurat_GSM3827118, seurat_GSM3827119,
#                           seurat_GSM3827120, seurat_GSM3827121, seurat_GSM3827122,
#                           seurat_GSM3827123, seurat_GSM3827124, seurat_GSM3827125,
#                           seurat_GSM3827126, seurat_GSM3827127, seurat_GSM3827128,
#                           seurat_GSM3827129, seurat_GSM3827130, seurat_GSM3827131,
#                           seurat_GSM3827132, seurat_GSM3827133, seurat_GSM3827134,
#                           seurat_GSM3827135)



# -----------------------------------------------------------------------------#
# <<<<<<<<<<<<<<<<<<< 四.自测10x数据读取的单细胞文件操作 >>>>>>>>>>>>>>>>>>>>>>#
# -----------------------------------------------------------------------------#
# 
# # 将文件夹设置为你要的样本名字，这样好对应
# 
# # 例子：单个文件结构
# #  - B1
# #      -- barcodes.tsv.gz
# #      -- features.tsv.gz
# #      -- matrix.mtx.gz
# 
# 
# # Read10X~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # 分别读取样本
# B1 <- Read10X("raw/B1/", gene.column = 2)
# B2 <- Read10X("raw/B2/", gene.column = 2)
# B3 <- Read10X("raw/B3/", gene.column = 2)
# TC <- Read10X("raw/TC/", gene.column = 2)
# Control <- Read10X("raw/Control/", gene.column = 2)
# 
# 
# 
# # 为每个样本构建Seurat对象~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# B1 <- CreateSeuratObject(counts = B1,
#                                 min.cells = 3,      #保留至少在3个细胞中表达的基因
#                                 min.features = 200) #保留至少具有200个基因的细胞
# 
# B2 <- CreateSeuratObject(counts = B2,
#                                 min.cells = 3,      #保留至少在3个细胞中表达的基因
#                                 min.features = 200) #保留至少具有200个基因的细胞
# 
# B3 <- CreateSeuratObject(counts = B3,
#                                 min.cells = 3,      #保留至少在3个细胞中表达的基因
#                                 min.features = 200) #保留至少具有200个基因的细胞
# 
# TC <- CreateSeuratObject(counts = TC,
#                                 min.cells = 3,      #保留至少在3个细胞中表达的基因
#                                 min.features = 200) #保留至少具有200个基因的细胞
# 
# Control <- CreateSeuratObject(counts = Control,
#                               min.cells = 3,      #保留至少在3个细胞中表达的基因
#                               min.features = 200) #保留至少具有200个基因的细胞
# 
# 
# # 分别为Seurat对象添加样本标签sample和tissue------------------------------------
# # 添加sample~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# B1[["sample"]] <- "B1"
# B2[["sample"]] <- "B2"
# B3[["sample"]] <- "B3"
# TC[["sample"]] <- "TC"
# Control[["sample"]] <- "Control"
# 
# 
# # 标注疾病和对照样本(实验分组信息)tissue~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# B1[["tissue"]] <- "Case"
# B2[["tissue"]] <- "Case"
# B3[["tissue"]] <- "Case"
# TC[["tissue"]] <- "Case"
# Control[["tissue"]] <- "Control"
# 
# 
# # 合并所有单个矩阵为矩阵list（seuratobject_data）--------------------------------
# seuratobject_data <- list(B1, B2, B3, TC, Control)





# ==============================================================================
# ----------------所有的不同读取方法读取完后 就来这里开始-----------------------
# ==============================================================================



# 1.合并Seurat对象==============================================================
sc_merge <- merge(seuratobject_data[[1]], seuratobject_data[-1])

# JoinLayers合并单独样本⼩矩阵为大矩阵
sc_merge <- JoinLayers(sc_merge) # 由于seurat5每个样本是单独的小矩阵，所以后面有步骤会无法运行，要用combineAandB <- JoinLayers(combineAandB)转换成大矩阵


#----------------------------保存sc_merge的rds---------------------------------#
# saveRDS(sc_merge, file = "save/sc_merge.rds")
sc_merge <- readRDS("save/sc_merge.rds")
#----------------------------保存sc_merge的rds---------------------------------#




# 2.质控========================================================================
# 2.1 基于多个QC指标初步过滤细胞------------------------------------------------

#------------------------------------------------------------------------------#
# sc_merge$log10GenesPerUMI <- log10(sc_merge$nFeature_RNA) / log10(sc_merge$nCount_RNA)
# 是哈佛教程里的步骤是18年的老步骤,现在可能过时了，但目前我们的项目流程里依旧留用
# 默认是log10GenesPerUMI > 0.80,大部分情况是合理的，但有的项目，过高的阈值可能会筛掉过多
#------------------------------------------------------------------------------#

# 计算线粒体基因百分比
# 人类为MT
# sc_merge[["percent.mt"]] = PercentageFeatureSet(sc_merge, pattern = "^MT-")/100
# 小鼠为mt
sc_merge[["percent.mt"]] = PercentageFeatureSet(sc_merge, pattern = "^mt-")/100


# 初步过滤~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # 初步过滤阈值（标准1）（含log10GenesPerUMI）
# sc_merge <- subset(sc_merge, subset= (nCount_RNA >= 500) &
#                   (nFeature_RNA >= 250) & # 常规是250
#                   (log10GenesPerUMI > 0.80) &   # 常规是0.8
#                   (percent.mt < 0.2)) # 常规是0.2


# 初步过滤阈值（标准2）
sc_merge <- subset(sc_merge,
                   subset= (nFeature_RNA >= 250) & # 常规是250
                     (nCount_RNA >= 500) &
                     (percent.mt < 0.2)) # 常规是0.2


# 从Seurat对象中提取初次质控后原始计数矩阵
counts <- GetAssayData(sc_merge,  layer = "counts")

# 检查 是否缺失分组
unique(sc_merge$sample)



# 2.2 基因质控，过滤基因--------------------------------------------------------
# 2.2.1 过滤低表达基因~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
nonzero <- counts > 0 # 找出哪些基因在哪些细胞中表达（非零）
keep_genes_expr <- Matrix::rowSums(nonzero) >= 10# 保留至少在10个细胞中表达的基因


# # 2.2.2 过滤性别特异性/MT-RNR2-like 核基因家族基因等干扰基因~~~~~~~~~~~~~~~~~~~~
# # 是从根本上从表达矩阵去除
# 
# # 人类原始通用过滤版本~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # ── 1. 性别特异性基因（男女一起）
# # 男性特有（Y 染色体编码，几乎只在男性表达）
# male_genes <- c(
#   "RPS4Y1", "DDX3Y", "KDM5D", "EIF1AY", "USP9Y", "UTY",
#   "NLGN4Y", "CYorf15A", "CYorf15B", "PRKY", "TMSB4Y", "ZFY"
# )
# 
# # 女性特有（XIST 是经典标志；其余为 X 失活相关或逃逸后仍男女差异明显）
# female_genes <- c("XIST", "TSIX", "JPX", "FTX")
# 
# sex_genes <- union(male_genes, female_genes)
# 
# # ── 2. MT-RNR2-like 核基因家族（11 个，无分型价值）
# mtrnr2l_genes <- c(
#   "MTRNR2L1", "MTRNR2L2", "MTRNR2L3", "MTRNR2L4", "MTRNR2L5",
#   "MTRNR2L6", "MTRNR2L7", "MTRNR2L8", "MTRNR2L9", "MTRNR2L10", "MTRNR2L11"
# )
# 
# ## ── 3. 核糖体基因：只去 top10 高表达的(要谨慎,不要乱删,看情况)
# # ribo_all <- grep("^RPL|^RPS", rownames(sc_merge), value = TRUE)
# # ribo_genes <- character(0)
# # if (length(ribo_all) > 0) {
# #   ribo_mean_expr <- rowMeans(counts[ribo_all, ])
# #   ribo_genes <- names(sort(ribo_mean_expr, decreasing = TRUE))[1:min(10, length(ribo_all))]
# # }
# 
# 
# # ── 4. 合并并去重，只保留实际存在于数据中的基因
# # exclude_genes <- Reduce(union, list(sex_genes, mtrnr2l_genes, ribo_genes))# 去除核糖体基因的版本
# exclude_genes <- Reduce(union, list(sex_genes, mtrnr2l_genes))# 不去核糖体基因的版本
# 
# exclude_genes <- intersect(exclude_genes, rownames(sc_merge))
# 
# # ── 5. 最终用于下游分析的干净基因
# highly_expressed_genes <- rownames(counts)[keep_genes_expr]
# keep_genes <- setdiff(highly_expressed_genes, exclude_genes)



# # 人类肌肉相关单细胞分析特化过滤~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # ── 1. 线粒体基因：保留！心肌细胞能量代谢核心，心衰关键指标
# # mt_genes <- grep("^MT-", rownames(sc_merge), value = TRUE)  # 注释掉，不去除
# 
# # ── 2. 性别特异性基因（男女混合样本，去除减少批次效应）
# male_genes <- c(
#   "RPS4Y1", "DDX3Y", "KDM5D", "EIF1AY", "USP9Y", "UTY",
#   "NLGN4Y", "CYorf15A", "CYorf15B", "PRKY", "TMSB4Y", "ZFY"
# )
# female_genes <- c("XIST", "TSIX", "JPX", "FTX")
# sex_genes <- union(male_genes, female_genes)
# 
# # ── 3. MT-RNR2-like 核基因家族（假基因，无分型价值）
# mtrnr2l_genes <- c(
#   "MTRNR2L1", "MTRNR2L2", "MTRNR2L3", "MTRNR2L4", "MTRNR2L5",
#   "MTRNR2L6", "MTRNR2L7", "MTRNR2L8", "MTRNR2L9", "MTRNR2L10", "MTRNR2L11"
# )
# 
# # ── 4. 核糖体蛋白基因：仅去除表达量最高的top10（技术噪音），保留其余
# # 修改：不去除全部，只去除最极端高表达的管家基因
# ribo_all <- grep("^RPL|^RPS", rownames(sc_merge), value = TRUE)
# ribo_mean_expr <- rowMeans(counts[ribo_all, ])
# ribo_genes <- names(sort(ribo_mean_expr, decreasing = TRUE))[1:10]  # 只去top10
# 
# # ── 5. 新增：血红蛋白基因（血液污染，心脏组织常见）
# blood_genes <- c("HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2")
# 
# # ── 6. 合并并去重，只保留实际存在于数据中的基因
# exclude_genes <- Reduce(union, list(sex_genes, mtrnr2l_genes, ribo_genes, blood_genes))
# exclude_genes <- intersect(exclude_genes, rownames(sc_merge))
# 
# # ── 7. 最终用于下游分析的干净基因
# highly_expressed_genes <- rownames(counts)[keep_genes_expr]
# keep_genes <- setdiff(highly_expressed_genes, exclude_genes)


# 小鼠通用过滤版本~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ── 1. 性别特异性基因（小鼠）
# 雄性特有
male_genes_mouse <- c(
  "Eif2s3y", "Kdm5d", "Uty", "Ddx3y", "Usp9y", "Zfy1", "Zfy2",
  "Sry", "Rbmy", "Prky", "Tspyl2"
)

# 雌性特有
female_genes_mouse <- c("Xist", "Tsix", "Jpx", "Ftx")

sex_genes_mouse <- union(male_genes_mouse, female_genes_mouse)

# ── 2. MT-RNR2-like 核基因家族
mtrnr2l_genes_mouse <- c(
  "Mtrnr2l1", "Mtrnr2l2", "Mtrnr2l3", "Mtrnr2l4", "Mtrnr2l5",
  "Mtrnr2l6", "Mtrnr2l7", "Mtrnr2l8", "Mtrnr2l9", "Mtrnr2l10", "Mtrnr2l11"
)

# ── 3. 合并并去重，只保留实际存在于数据中的基因
exclude_genes_mouse <- Reduce(union, list(sex_genes_mouse, mtrnr2l_genes_mouse))
exclude_genes_mouse <- intersect(exclude_genes_mouse, rownames(sc_merge))

# ── 4. 最终用于下游分析的干净基因
highly_expressed_genes <- rownames(counts)[keep_genes_expr]
keep_genes <- setdiff(highly_expressed_genes, exclude_genes_mouse)



# 2.2.3 保留筛选后的矩阵~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
filtered_counts <- counts[keep_genes, ]


# 2.3 创建过滤后的Seurat对象----------------------------------------------------
sc_filter <- CreateSeuratObject(filtered_counts,
                                meta.data = sc_merge@meta.data)#同时将sc_merge对象中的元数据复制到新创建的Seurat对象sc_filter



# 2.4 绘制质控小提琴图初步查看质控结果------------------------------------------
# 所有样本混合的小提琴图(整体质量)
vlnplot <- VlnPlot(sc_filter, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                   pt.size = 0, cols = c('#E8D2B3','#BD9AAD'), ncol = 3, group.by = "tissue")
vlnplot <- vlnplot & 
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1))
vlnplot
unique(sc_filter$sample)


# 每个样本单独的小提琴图(个体质量)
sample_colors <- colorRampPalette(brewer.pal(12, "Set3"))(length(unique(sc_filter$sample)))
vlnplot2 <- VlnPlot(sc_filter, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                    pt.size = 0, cols = sample_colors, ncol = 3, group.by = "sample")
vlnplot2 <- vlnplot2 & 
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1))
vlnplot2



# 2.5 依据小提琴图再进行质控----------------------------------------------------
# 依据小提琴图再质控

# 2.5.1再过滤~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 再过滤阈值
sc_filter <- subset(sc_filter,
                    subset= (nFeature_RNA >= 250) & # 常规是250
                      (nFeature_RNA <= 7500) &
                      (nCount_RNA >= 500) &
                      (nCount_RNA <= 40000) &
                      (percent.mt < 0.2)) # 常规是0.2


# 2.5.2按样本为单位细致质控~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 使用前需要先进行前面2步,初步质控和再质控后再进行！！！！！！！！！！！！！！！
# 只是过滤个别样本的极端值

# 对每个样本自适应计算上限阈值（median + n*MAD），只过滤离群高值
# nFeature_RNA / nCount_RNA：过高可能为双细胞，percent.mt：过高为破损细胞
n_mad <- 5  # 越大越宽松 一般为3,宽松选5

upper_thresholds <- FetchData(sc_filter, vars = c("sample", "nFeature_RNA", "nCount_RNA", "percent.mt")) %>%
  group_by(sample) %>%
  summarise(across(c(nFeature_RNA, nCount_RNA, percent.mt),
                   ~ median(.) + n_mad * mad(.)), .groups = "drop")

keep_cells <- FetchData(sc_filter, vars = "sample") %>%
  left_join(upper_thresholds, by = "sample") %>%
  { sc_filter$nFeature_RNA <= .$nFeature_RNA &
      sc_filter$nCount_RNA  <= .$nCount_RNA &
      sc_filter$percent.mt  <= .$percent.mt }

sc_filter <- sc_filter[, keep_cells]


# 从Seurat对象中提取初次质控后原始计数矩阵
counts <- GetAssayData(sc_filter,  layer = "counts")

# 检查 是否缺失分组
unique(sc_filter$sample)


# 2.5.3绘制最终质控小提琴图~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 所有样本混合的小提琴图(整体质量)
vlnplot <- VlnPlot(sc_filter, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                   pt.size = 0, cols = c('#E8D2B3','#BD9AAD'), ncol = 3, group.by = "tissue")
# 添加边框
vlnplot <- vlnplot & 
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1))

vlnplot
ggsave(filename = "output/1-ViolinPlot1.pdf", plot = vlnplot, width = 16.6, height = 12, units = "cm")


# 每个样本单独的小提琴图(个体质量)
sample_colors <- colorRampPalette(brewer.pal(12, "Set3"))(length(unique(sc_filter$sample)))
vlnplot2 <- VlnPlot(sc_filter, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                    pt.size = 0, cols = sample_colors, ncol = 3, group.by = "sample")
vlnplot2 <- vlnplot2 & 
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1))
vlnplot2
ggsave(filename = "output/1-ViolinPlot2.pdf", plot = vlnplot2, width = 16.6, height = 12, units = "cm")



# # 3.decontX 去除环境RNA=========================================================
# 
# #-----------------------------------说明---------------------------------------#
# # 得再去除双细胞前,未进行降维聚类前进行
# # 如果进行降维聚类可视化后的Umap/tsne图有很多散的毛边，就可以考虑进行去除环境RNA，没有就不用
# # 去除环境RNA后可能会去除一些过渡态的细胞，会影响后续分析，特别是有粘连过渡态的细胞的肿瘤样本，慎用！！！
# # 具体是否要用具体根据项目判断
# # decontX可以输入两类数据：SingleCellExperiment对象，或counts矩阵。
# # decontX的运行有两种方式，一种是直接跑，不提供背景。另外一种是使用包含空液滴的原始/液滴矩阵作为背景。通常是第一种
# #-----------------------------------说明---------------------------------------#
# 
# # 提取表达矩阵
# counts <- GetAssayData(sc_filter, layer = "counts")
# 
# # 运行decontX
# decontX_res <- decontX(counts,
#                        batch = sc_filter$sample) # 按样本分批进行处理,重要参数
# 
# # 存储结果中获取的污染值
# sc_filter$contamination <- decontX_res$contamination
# 
# # 保存校正后的表达矩阵
# sc_filter[["decontX"]] <- CreateAssay5Object(
#   counts = decontX_res$decontXcounts  # 注意这里是校正矩阵
# )
# 
# # 查看污染率分布，决定过滤阈值
# hist(sc_filter$contamination, breaks = 50)
# abline(v = 0.2, col = "red", lwd = 2) # 红线标出0.2阈值
# 
# # 过滤污染率高的细胞(保留污染值<0.2的细胞)(和下面一步2选1)(但不建议选这个步骤)
# # sc_filter <- sc_filter[, sc_filter$contamination < 0.2]
# 
# # 切换默认assay为校正后的数据(重要步骤) (和上一步2选1)
# DefaultAssay(sc_filter) <- "decontX"



# 4. 质控后的标准化/鉴定高变基因/归一化/去细胞周期影响/线性降维/去批次比较======
# 4.1 标准化、鉴定高变基因、归一化----------------------------------------------
sc_filter <- NormalizeData(sc_filter)  #标准化
sc_filter <- FindVariableFeatures(sc_filter,nfeatures = 2000)  #鉴定2000高变基因，常规是2000高变基因
sc_filter <- ScaleData(sc_filter)  #归一化, 如果不进行去除细胞周期直接运行后跳转 4.3 PCA


# 4.2 去细胞周期（可选）（默认不去除）------------------------------------------
# 只是对细胞降维聚类的影响，不影响原始表达的数据

#-----------------------ScaleData()/SCTransform() 区别:------------------------#
# ScaleData():老流程的函数,必须在 NormalizeData() 和 FindVariableFeatures() 之后调用。
# SCTransform():新函数,此函数整合了NormalizeData() 和 FindVariableFeatures(),一步调用即可
# SCTransform()也可以用于去除细胞周期效应
# 用了SCTransform()的话，后续的（如 DotPlot, FindAllMarkers, DimPlot 等）前什么的得先 DefaultAssay(seurat对象) <- "RNA"
# 当你对 Seurat 对象进行下游分析（如 DotPlot, FindAllMarkers, DimPlot 等）时，Seurat 默认会使用当前的 DefaultAssay 所指向的数据层。
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# 注意：如果分析的数据为肿瘤数据则不建议去除细胞周期的影响！！！
# 由于肿瘤的发生与细胞周期相关的基因有一定联系，去除细胞周期可能导致错误的分析结果
# 一般不建议对肿瘤单细胞测序数据进行细胞周期去除
# 但如心肌细胞：终末分化，不增殖，细胞周期影响小
#------------------------------------------------------------------------------#


# 4.2.1 计算细胞周期分数~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 获取G2M期相关基因
g2m_genes <- cc.genes$g2m.genes
g2m_genes <- CaseMatch(search=g2m_genes, match=rownames(sc_filter))

# 获取S期相关基因
s_genes <- cc.genes$s.genes
s_genes <- CaseMatch(search=s_genes, match=rownames(sc_filter))

# 细胞周期评分
sc_filter <- CellCycleScoring(sc_filter,
                              s.features = s_genes,
                              g2m.features = g2m_genes,
                              set.ident = TRUE
)
# 可以对去除细胞周期效应前后进行针对细胞周期基因的PCA查看比较,但此处省略

# 4.2.2 缩放数据并回归协变量,将细胞周期的影响回归掉~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sc_filter <- ScaleData(sc_filter,
                       vars.to.regress =  c("S.Score", "G2M.Score"))

# 查看小提琴图看评分，后接4.3
VlnPlot(sc_filter, features = c("S.Score", "G2M.Score"), group.by = "GEOID")



# 4.3 PCA线性降维 --------------------------------------------------------------
sc_filter <- RunPCA(sc_filter) 
# 如果是去除细胞周期影响后,RunPCA(sc_filter) 默认用的去除细胞周期后的scale.data数据

# 查看前30 PCA结果(仅供查看)
ElbowPlot(sc_filter, ndims = 30)  # 确定重要PCs数量

# 查看去除细胞周期后的效果(如果进行了去除细胞周期效应)
DimPlot(sc_filter, reduction = "pca", group.by = "Phase")  # 按细胞周期阶段着色
DimPlot(sc_filter, reduction = "pca", group.by = "GEOID")  # 按样本分组


# 4.4 scDblFinder 去双细胞（可选）----------------------------------------------
# 必须在Harmony去批次前

# 4.4.1 将seurat对象转换为 SingleCellExperiment SCE对象~~~~~~~~~~~~~~~~~~~~~~~~~
DefaultAssay(sc_filter)  # 如果上面使用了了decontX矫正的矩阵,那应该显示 "decontX" 要确认
library(zellkonverter) 
sc_filter <- as.SingleCellExperiment(sc_filter)  # 默认会提取 counts 层 如果上面使用了decontX矫正的矩阵，那会默认使用decontX矫正的矩阵层


# 4.4.2 基于单样本~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 看具体的教程 https://github.com/plger/scDblFinder 


# 4.4.3 基于多样本~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(BiocParallel) # 多线程并行运行

# 运行scDblFinder
set.seed(123)  # ← 随机种子保证scDblFinder结果可复现,scDblFinder算法有随机性
sc_filter <- scDblFinder(sc_filter,
                         samples="sample",# 按样本，对每个样本进行去双细胞
                         BPPARAM=MulticoreParam(20))

# 统计结果
table(sc_filter$scDblFinder.class) 

# 转为seurat对象
sc_filter <- as.Seurat(sc_filter) # 默认转换,可能不能用，如果前面不用decontX矫正的矩阵就可以运行这个 as.Seurat(sc_filter) 默认会使用 counts 作为counts矩阵。
# sc_filter <- as.Seurat(sc_filter, counts = "counts", data = "counts")
DefaultAssay(sc_filter) <- "RNA"   # 再加这个固定使用标化后的“RNA”数据层的数据进行后续分析，保险！！

# 只提取单细胞的细胞，即去除双细胞
sc_filter <- subset(sc_filter, subset = scDblFinder.class == "singlet")

sc_filter$orig.ident <- NULL # 去除元数据中的orig.ident,这个数据是无用的，必须去除，不去除可能会影响后续的流程的连贯性
sc_filter$ident <- NULL # 去除元数据中的ident,这个数据是无用的，必须去除，不去除可能会影响后续的流程的连贯性

table(sc_filter$scDblFinder.class, useNA = "ifany") # 查看是否去干净



#--------------------特别注意，如果进行了去双细胞，得运行----------------------#
# 去除双细胞分析中的as.SingleCellExperiment会丢失降维信息,得重新如下！！！！！！

# 重新 标准化、鉴定高变基因、归一化、PCA线性降维
sc_filter <- NormalizeData(sc_filter) 
sc_filter <- FindVariableFeatures(sc_filter,nfeatures = 2000) #和上面对应
sc_filter <- ScaleData(sc_filter)  # 不进行去细胞周期的分析就直接运行这个，与上面的NormalizeData/FindVariableFeatures，2步

# 如果进行了去细胞周期的分析，上一步的：sc_filter <- ScaleData(sc_filter)  要改为如下2步！！！！！！！！
# 1.对去除双细胞后的所有细胞的周期效应进行重新打分
sc_filter <- CellCycleScoring(sc_filter,
                              s.features = s_genes,
                              g2m.features = g2m_genes,
                              set.ident = TRUE
)
# 2.重新按照周期评分进行去周期效应
sc_filter <- ScaleData(sc_filter,
                       vars.to.regress =  c("S.Score", "G2M.Score"))


# 然后再重新PCA
sc_filter <- RunPCA(sc_filter)
#-----------------------------特别注意，得运行----------------------------------#



# 4.5 在PCA后进行批次效应校正(Harmony去批次)------------------------------------
# 去批次必须以单个样本为单位进行去批次，千万不能以分组为单位去批次

set.seed(123) # 确保可重复性

# 默认参数Harmony
sc_integer <- sc_filter %>% 
  harmony::RunHarmony(group.by.vars = "sample")

# # 精细Harmony(可选)
# sc_integer <- harmony::RunHarmony(
#   object = sc_integer,
#   group.by.vars = "sample",     # 指定批次变量
#   max.iter = 30,                # 最大迭代次数,默认 max.iter = 10
#   lambda = 3                    # 调整聚类数,默认 lambda = 1,越大矫正越强
#   )

#------------------------------------------------------------------------------#
# 后续分析中：所有基于降维的操作（如RunUMAP、聚类）需指定 reduction = "harmony"
# harmony后是添加了独⽴的嵌⼊矩阵，默认不会运⾏harmony后的效果，要指定
#------------------------------------------------------------------------------#

# 可视化比较
p1 <- DimPlot(sc_filter, reduction = "pca", pt.size = 0.1, group.by = "sample") + ggtitle("Before Harmony")
p2 <- DimPlot(sc_integer, reduction = "harmony", pt.size = 0.1, group.by = "sample") + ggtitle("After Harmony")
p1+p2


# 查看细胞数量
length(colnames(sc_integer))
# 更正后：14260


# 可视化 PCA 结果
pca_plot <- DimPlot(sc_filter, reduction = "pca",
                    pt.size = 0.05,
                    cols = More_colors,
                    # alpha = 0.5, # 点透明度
                    group.by = "sample")+
  labs(x = "PCA 1", y = "PCA 2") +
  guides(color = guide_legend(override.aes = list(size = 4))) +                 # 调大图例圆圈
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1), # 添加边框
        axis.line = element_blank(),                                            # 去掉默认轴线
        axis.ticks = element_blank(),                                           # 去掉刻度线
        axis.text = element_blank(),                                            # 去掉刻度标签
        axis.title = element_text(size = 12, color = "black")  # 控制坐标轴标签字体大小
  )
pca_plot
ggsave("output/2-PCA_plot.pdf", pca_plot, width = 16, height = 13, units = "cm")



## 5 细胞聚类降维分析===========================================================

# 5.1 细胞聚类------------------------------------------------------------------
# 去批次后的数据必须加上reduction = "harmony"！！！
sc_integer <- FindNeighbors(sc_integer,dims = 1:20, verbose = T,reduction = "harmony")  # 构建KNN图，将前20个PC作为输入

# 按照每0.2分辨率为梯度，分别计算从0.1至2.0范围内的细胞聚类信息
sc_integer <- FindClusters(sc_integer,reduction = "harmony", resolution = c(seq(.1,2,.2))) 
# 绘制分辨率聚类树，Seurat对象sc_integer的meta.data中包含以 RNA_snn_res.开头的多列聚类结果,
clustree(sc_integer@meta.data, prefix = "RNA_snn_res.")

sc_integer <- FindClusters(sc_integer,reduction = "harmony", resolution = 0.5) # 聚类，分辨率0.5


# 5.2 非线性降维----------------------------------------------------------------
set.seed(123) # 确保可重复性
# UMAP降维~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# UMAP降维选择1
sc_integer <- RunUMAP(sc_integer,
                      dims = 1:20,
                      min.dist = 0.3,#控制点的聚合程度，越小点越紧密，越大点越松散,默认0.3,RunUMAP特有参数
                      n.neighbors = 30,#UMAP 计算的邻居点数量，决定局部与全局结构,越大越聚焦全局结构,默认30
                      reduction = "harmony")

# # UMAP降维选择2
# sc_integer <- RunUMAP(sc_integer,
#                       dims = 1:30,
#                       reduction = "harmony",
#                       metric = "cosine",        # 距离度量方式
#                       min.dist = 0.3,           # 控制点的聚合程度
#                       n.neighbors = 30          # 邻居点数量
#                      )

# # TSNE降维~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# sc_integer <- RunTSNE(sc_integer,
#                       dims = 1:30,
#                       perplexity = 30,# 决定局部与全局结构,越大越聚焦全局结构,默认是30,类似n.neighbors 
#                       reduction = "harmony")


#---------------------------保存sc_integer的rds--------------------------------#
# saveRDS(sc_integer, "save/sc_integer.rds")
sc_integer <- readRDS("save/sc_integer.rds")
#---------------------------保存sc_integer的rds--------------------------------#



# 5.3 可视化--------------------------------------------------------------------
# UMAP可视化
UMAP1 <- DimPlot(sc_integer,reduction = "umap",
                 label = T,
                 #alpha = 0.4, # 点透明度
                 cols = More_colors,
                 pt.size = 0.05,) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  guides(color = guide_legend(override.aes = list(size = 4))) +                 # 调大图例圆圈
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1), # 添加边框
        axis.line = element_blank(),                                            # 去掉默认轴线
        axis.ticks = element_blank(),                                           # 去掉刻度线
        axis.text = element_blank(),                                            # 去掉刻度标签
        axis.title = element_text(size = 12, color = "black")  # 控制坐标轴标签字体大小
  )
UMAP1
ggsave(filename = "output/3-UMAP.pdf", plot = UMAP1, width = 18, height = 14, units = "cm") # 小尺寸
# ggsave(filename = "output/3-UMAP.pdf", plot = UMAP1, width = 21, height = 16, units = "cm") # 大尺寸


# # TSNE可视化
# TSNE1 <- DimPlot(sc_integer,
#                  reduction = "tsne",
#                  label=T,
#                  # cols = More_colors,
#                  # alpha = 0.4, # 点透明度
#                  pt.size = 0.05)+
#   labs(x = "UMAP 1", y = "UMAP 2") +
#   guides(color = guide_legend(override.aes = list(size = 4))) +                 # 调大图例圆圈
# theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1), # 添加边框
#       axis.line = element_blank(),                                            # 去掉默认轴线
#       axis.ticks = element_blank(),                                           # 去掉刻度线
#       axis.text = element_blank(),                                            # 去掉刻度标签
#       axis.title = element_text(size = 12, color = "black")  # 控制坐标轴标签字体大小
# )
# TSNE1
# ggsave(filename = "output/3-TSNE.pdf", plot = TSNE1, width = 17, height =13, units = "cm")




# 5.3 查看不同样本细胞分布情况--------------------------------------------------
sample_distribution_UMAP <- DimPlot(sc_integer,
                                    reduction = "umap",
                                    label=T,
                                    cols = More_colors,
                                    # alpha = 0.4, # 点透明度
                                    group.by = "sample")+
  labs(x = "UMAP 1", y = "UMAP 2") +
  guides(color = guide_legend(override.aes = list(size = 4))) +                 # 调大图例圆圈
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1), # 添加边框
        axis.line = element_blank(),                                            # 去掉默认轴线
        axis.ticks = element_blank(),                                           # 去掉刻度线
        axis.text = element_blank(),                                            # 去掉刻度标签
        axis.title = element_text(size = 12, color = "black")         # 控制坐标轴标签字体大小
  )
sample_distribution_UMAP
ggsave(filename = "output/4-sample_distribution-UMAP.pdf", plot = sample_distribution_UMAP, width = 16, height = 12, units = "cm")


# sample_distribution_TSNE <- DimPlot(sc_integer,
#                                     reduction = "tsne",
#                                     alpha = 0.4, # 点透明度
#                                     label=T,
#                                   # cols = More_colors,
#                                     group.by = "sample")+
#   labs(x = "UMAP 1", y = "UMAP 2") +
#   guides(color = guide_legend(override.aes = list(size = 4))) +                 # 调大图例圆圈
#   theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1), # 添加边框
#         axis.line = element_blank(),                                            # 去掉默认轴线
#         axis.ticks = element_blank(),                                           # 去掉刻度线
#         axis.text = element_blank(),                                            # 去掉刻度标签
#         axis.title = element_text(size = 12, color = "black")         # 控制坐标轴标签字体大小

#   )
# sample_distribution_TSNE
# ggsave(filename = "output/4-sample_distribution-TSNE.pdf", plot = sample_distribution_TSNE, width = 15, height = 12, units = "cm")




# # 6 检查前面bulk RNA-Seq筛到的基因在scRNA-Seq中是否全表达 ======================
# 
# # -------------------------------话术-------------------------------------------
# # 如果进行了取交集！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！
# 
# # 为确保候选基因在单细胞分辨率下具有生物学可解释性，我们首先构建了 bulk RNA-seq 与 scRNA-seq 
# # 数据集中共同检测到的基因集合（n = X,XXX）。在此交集基因集上，
# # 我们重新进行了差异表达分析与机器学习特征筛选，最终鉴定出 Y 个高置信度候选基因。
# # 这些基因在 scRNA-seq 数据中均显示出明确的细胞类型特异性表达模式（图X），
# # 支持其作为潜在功能基因或生物标志物的可靠性。”
# # -------------------------------话术-------------------------------------------


# # 人类基因转成小鼠基因>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# # 读取人类hubGenes
# model_marker <-  data.table::fread("input/hubGenes.csv",data.table = F)[,1]
#
# # 获取基因列表
# genelist<- model_marker
#
# #使用homologene函进行转换
# #@genelist是要转换的基因列表
# #@inTax是输入的基因列表所属的物种号，10090是小鼠
# #@outTax是要转换成的物种号，9606是人
#
# # 开始转换
# genelistid <- homologene(genelist, inTax = 9606, outTax = 10090)
# # 提取1,2列
# genelistid2 <- genelistid[,c(1,2)]
#
# colnames(genelistid2) <- c("ID","gene")
# genelistid2 <- as.data.frame(genelistid2)
#
# model_marker <- genelistid2$gene

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# 读取模型基因/要聚焦的基因
# model_marker <-  read.table("input/hubgene.txt",header=T)[,1]
# model_marker <- read.csv("input/8-LASSO_hubgenes.csv", header = TRUE)[, 1]
# model_marker <-  c("BCHE","CYP1A2","TTR")


# 6.1 检查模型基因是否在单细胞数据集中表达--------------------------------------
# # 检查基因是否存在于数据
# model_marker_exist <- model_marker[model_marker %in% rownames(sc_integer)]
# model_marker_missing <- model_marker[!model_marker %in% rownames(sc_integer)]
# 
# if(length(model_marker_missing) > 0){
#   cat("以下基因在单细胞数据中未找到:", paste(model_marker_missing, collapse = ", "), "\n")
# }
# cat("成功匹配的基因:", paste(model_marker_exist, collapse = ", "), "\n")


# 6.2 提取sc_integer中所有基因名字到数据----------------------------------------
# sc_genes <- data.frame(gene = rownames(sc_integer))
# cat("单细胞数据集中的基因总数:", nrow(sc_genes), "\n")
#
# # 输出csv文件用于和bulk转录组数据取交集
# write.csv(sc_genes, "output/sc_all_genes.csv", row.names = FALSE)


save.image(file = "save/SCQC.RData")


