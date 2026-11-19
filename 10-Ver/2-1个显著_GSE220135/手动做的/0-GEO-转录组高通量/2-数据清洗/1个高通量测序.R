# 添加bioconductor清华镜像
options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
# 添加bioconductor西湖大学镜像
options(BioC_mirror = "https://mirrors.westlake.edu.cn/bioconductor")

# 添加R清华镜像
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))


# 1 环境变量准备================================================================
# 1.1 清空环境------------------------------------------------------------------
rm(list=ls())
gc()

# # 1.2 R包安装-------------------------------------------------------------------
# # 定义一个函数来检查并安装所需的R包
# install_if_missing <- function(packages) {
#   installed_packages <- rownames(installed.packages())
#   for (pkg in packages) {
#     if (!(pkg %in% installed_packages)) {
#       BiocManager::install(pkg, dependencies = TRUE, update = FALSE)
#     }
#   }
# }
# 
# # 检查安装的包列表
# required_packages <- c("GEOquery","dplyr","magrittr","ggplot2","sva","limma",
#                        "FactoMineR","factoextra")
# 
# 
# # 调用函数进行检查和安装
# install_if_missing(required_packages)
# 
# BiocManager::install("genefilter")
# BiocManager::install("Biostrings")
# BiocManager::install("GEOquery")
# 

# 1.3 加载R包------------------------------------------------------------------
library(GEOquery)
library(dplyr)
library(magrittr)
library(ggplot2)
library(sva)
library(limma)
library(FactoMineR)
library(factoextra)
library(genefilter)
library(Biostrings)
library(magrittr)  # 加载 magrittr 包
library(tidyr)
library(stringr)
library(biomaRt)
library(readxl)
library(homologene)
library(clusterProfiler)
library(org.Hs.eg.db)


# 2.2 定义所需函数-------------------------------------------------------------
# 判断是否需要log2处理
log2_if_need <- function(data)
{
  print(c("before:",range(data)))
  ex <- as.data.frame(data)
  qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = T))
  LogC <- (qx[5] > 100) ||
    (qx[6]-qx[1] > 50 && qx[2] > 0) ||
    (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
  ## 如果LogC=T,就进行循环，如果存在小于0的数，就用非数填入
  if (LogC) {
    for(i in 1:ncol(ex)){
      ex[which(ex[,i] < 0),i] <- NaN
    }
    data <- log2(data + 1)
    print("log2 transform finished")
    print(c("after:",range(data)))
  }else{
    print("log2 transform not needed")
  }
  return(data)
}

# 整理分组信息
geo_group <- function(pdata, sample_id, group_id, treat_name, con_name, 
                      treat_key, con_key, desc){
  geogroup <- pdata %>% 
    dplyr::select(ID = sample_id, group0 = group_id) %>%
    dplyr::mutate(group = dplyr::case_when(
      grepl(treat_key, group0) ~ treat_name,
      grepl(con_key, group0) ~ con_name
    )) %>% 
    dplyr::select(1,3) %>% 
    dplyr::filter(group == treat_name | group == con_name) %>% 
    dplyr::arrange(group)
  if(desc){
    geogroup <- geogroup %>% 
      dplyr::arrange(desc(group))
  }
  message("Please check the order of samples!")
  return(geogroup)
}

# 整理芯片平台
gpl_pro <- function(gpl_name, probe_id, gene_id)
{
  gpl_name$gene_id <- strsplit(gpl_name[[gene_id]], split = " /// ") %>%
    lapply(function(x){
      x1 <- x
      if(length(x1) > 1) x1 <- x1[!grepl("^OTTHU", x1)]
      if(length(x1) == 0) x1 <- x[1]
      vec_gene = strsplit(x, split = " // ") %>% lapply(
        function(x) if(length(x) > 1) 
          return(x[2]) else return(x)) %>%
        unlist() %>% unique()
      vec_gene <- vec_gene[!grepl("^OTTHU", vec_gene)]
      vec_gene <- paste0(vec_gene, collapse = " /// ")
      return(vec_gene)
    }) %>% unlist()
  gpl_name <- gpl_name %>% 
    dplyr::select(probe_id, gene_id) %>% # 选择探针名和基因名列
    dplyr::filter(!grepl("///", .[,2])) %>% # 去除基因名中有///的
    dplyr::filter(.[,2] != "") %>% # 去除基因名为空的
    dplyr::filter(!grepl("---", .[,2])) %>% # 去除基因名中有---的
    dplyr::mutate(ID = as.character(ID)) %>% # 探针列转字符
    dplyr::rename(GENE_SYMBOL = gene_id) %>%  # 重命名基因列为GENE_SYMBOL
    dplyr::select(probe_id, GENE_SYMBOL)
  return(gpl_name)
}

# 探针ID转换，同一基因取均值
geo_matrix <- function(data_matrix, data_group, gpl_name){
  data_matrix <- data_matrix[, data_group$ID] %>% # 按group筛选并排序样本
    data.frame() %>% 
    dplyr::mutate(across(!where(is.numeric), as.numeric)) %>% 
    dplyr::mutate(ID = row.names(.) %>% as.character()) %>% 
    dplyr::left_join(., gpl_name, by = "ID") %>% 
    na.omit()
  gene_row <- data_matrix$GENE_SYMBOL
  data_matrix <- data_matrix %>% 
    dplyr::select(-ID, -GENE_SYMBOL) %>% 
    as.matrix()
  rownames(data_matrix) <- gene_row
  data_matrix <- data_matrix %>%  
    limma::avereps() %>%
    data.frame()
  print(dim(data_matrix))
  return(data_matrix)
}



# 1.4 创建文件夹---------------------------------------------------------------
dir.create("input/")
dir.create("output/")


# 1.5 定义变量-----------------------------------------------------------------
# # 训练集
GSE_2 <- "GSE220135"   # 修改：


treat_name = "BPD"    # 修改：疾病组命名
con_name = "Control"  # 对照组命名



# ========================= 1个高通量验证集的处理 ==============================

# 5.1 读取数据------------------------------------------------------------------
# 5.1.1 标准化后的TPM数据~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 读取 Step 2 生成的 TPM 矩阵
GSE220135 <- read.csv("../1-Count-to-TPM/output/TPM_matrix.csv", row.names = 1,
                      check.names = FALSE)

cat(paste0("TPM 矩阵维度: ", nrow(GSE220135), " 基因 x ",
           ncol(GSE220135), " 样本\n"))


# # 5.1.2 用于DEseq2差异分析的数据~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # 导入原始计数矩阵
# GSE220135_raw_count <- data.table::fread("input/count/GSE220135_raw_counts_GRCh38.p13_NCBI.tsv") 
# 
# # 转换为 data.frame 并设置行名
# GSE220135_raw_count <- as.data.frame(GSE220135_raw_count)
# rownames(GSE220135_raw_count) <- GSE220135_raw_count[, 1]
# GSE220135_raw_count <- GSE220135_raw_count[, -1]  # 移除第一列


## 5.3 数据清洗-----------------------------------------------------------------
# 5.2.1 标准化后的数据（TPM/FPKM）的处理~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
range(GSE220135)
GSE220135 <- na.omit(GSE220135)##检查每一行，如果某一行中至少有一个元素是 NA，则整行都会被移除
GSE220135_adj <- log2_if_need(GSE220135)

# GSE220135_adj <- log2(GSE220135 + 1)#手动log2化


# # 5.3.2 原始count数据的处理~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # 不需要log化！！！！！！！！
# # 去除含有NA的行
# range(GSE220135_raw_count)
# GSE220135_raw_count <- na.omit(GSE220135_raw_count)##检查每一行，如果某一行中至少有一个元素是 NA，则整行都会被移除


# 5.4 整理分组信息------------------------------------------------------------
# 读取本地的 series_matrix.txt 文件获取分组信息
aset_2 <- getGEO(filename = "input/GSE220135_series_matrix.txt/GSE220135_series_matrix.txt",
                 AnnotGPL = F, getGPL = F)

PD_2 <- pData(aset_2)

group_2 <- geo_group(pdata = PD_2,
                     # 样本编号列名；分组信息列名（修改）
                     sample_id = "geo_accession", group_id = "title", 
                     # 疾病组命名；对照组命名；是否降序可修改）
                     treat_name = treat_name, con_name = con_name, desc = T, 
                     # 疾病组关键字符串；对照组关键字符串（修改）
                     treat_key = "_BPD", con_key = "_nonBPD") 



# 5.6 矩阵行基因名转换----------------------------------------------------------
# # <<<<<<<<<<<<<<<<<<<<<<<（entrez → symbol）版本>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# # 读取映射文件
# id_map <- read.csv("input/id_trans.csv", stringsAsFactors = FALSE)
# 
# #  提取有效映射（去掉空 symbol）
# id_map <- id_map %>%
#   filter(!is.na(entrez_id), entrez_id != "", !is.na(gene_name), gene_name != "")
# 
# # 创建命名向量（entrez → symbol）
# entrez_to_symbol <- setNames(id_map$gene_name, id_map$entrez_id)
# 
# # 应用到矩阵（只保留能映射的行）
# 
# # 辅助函数：把行名从 entrez_id 转成 gene_symbol，
# # 若多个 entrez 对应同一 symbol，则合并这些行并取平均表达
# convert_rownames_to_symbol <- function(mat, mapping) {
#   common <- intersect(rownames(mat), names(mapping))
#   mat <- mat[common, ]
#   rownames(mat) <- mapping[rownames(mat)]
#   
#   # 如果有重复 symbol，取平均表达
#   mat <- mat %>%
#     tibble::rownames_to_column("symbol") %>%
#     group_by(symbol) %>%
#     summarise(across(everything(), mean, na.rm = TRUE)) %>%
#     tibble::column_to_rownames("symbol") %>%
#     as.data.frame()
#   
#   return(mat)
# }
# 
# # 应用到标准化后的数据集~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# GSE220135_adj <- convert_rownames_to_symbol(GSE220135_adj, entrez_to_symbol)
# # # 应用到原始count数据集~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # GSE220135_raw_count_adj <- convert_rownames_to_symbol(GSE220135_raw_count_adj, entrez_to_symbol)


# <<<<<<<<<<<<<<<<<<<<<<<（ensembl → symbol）版本>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# # 读取映射文件
# id_map <- read.csv("input/id_trans.csv", stringsAsFactors = FALSE)
# 
# #  提取有效映射（去掉空 symbol）
# id_map <- id_map %>%
#   filter(!is.na(id), entrez_id != "", !is.na(gene_name), gene_name != "")
# 
# # 创建命名向量（ensembl → symbol）
# ensembl_to_symbol <- setNames(id_map$gene_name, id_map$id)
# 
# # 应用到矩阵（只保留能映射的行）
# 
# # 辅助函数：把行名从 ensembl 转成 gene_symbol，
# # 若多个 ensembl 对应同一 symbol，则合并这些行并取平均表达
# convert_rownames_to_symbol <- function(mat, mapping) {
#   common <- intersect(rownames(mat), names(mapping))
#   mat <- mat[common, ]
#   rownames(mat) <- mapping[rownames(mat)]
#   
#   # 如果有重复 symbol，取平均表达
#   mat <- mat %>%
#     tibble::rownames_to_column("symbol") %>%
#     group_by(symbol) %>%
#     summarise(across(everything(), mean, na.rm = TRUE)) %>%
#     tibble::column_to_rownames("symbol") %>%
#     as.data.frame()
#   
#   return(mat)
# }

# # 应用到标准化后的数据集~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# GSE220135_adj <- convert_rownames_to_symbol(GSE220135_adj, ensembl_to_symbol)
# # 应用到原始count数据集~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# GSE220135_raw_count_adj <- convert_rownames_to_symbol(GSE220135_raw_count_adj, ensembl_to_symbol)


# 将GSE220135_adj 矩阵里的小鼠基因转成人源基因~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# （本数据集为人类数据，无需转换）
# # 获取GSE220135_adj数据集的所有小鼠基因并去重
# mouse_genes <- unique(rownames(GSE220135_adj))
# 
# # 进行基因转换（小鼠taxid=10090，人类taxid=9606）
# homolog_genes <- homologene(mouse_genes, inTax = 10090, outTax = 9606)
# 
# # 提取基因映射关系(第1列是小鼠基因,第2列是人类基因)并去重
# homolog_genes <- homolog_genes[, c(1, 2)]
# colnames(homolog_genes) <- c("mouse_gene", "human_gene")
# homolog_genes <- homolog_genes[!duplicated(homolog_genes$mouse_gene), ]
# rownames(homolog_genes) <- homolog_genes$mouse_gene
# 
# # 定义基因转换函数
# convert_mouse_to_human <- function(expr_matrix, mapping) {
#   common_genes <- intersect(rownames(expr_matrix), rownames(mapping))
#   expr_subset <- expr_matrix[common_genes, , drop = FALSE]
#   human_genes <- mapping[common_genes, "human_gene"]
#   
#   # 如果有重复的人类基因，取平均值
#   if (any(duplicated(human_genes))) {
#     # 使用 aggregate 对重复基因取平均值
#     expr_subset$human_gene <- human_genes
#     expr_subset <- aggregate(. ~ human_gene, data = expr_subset, FUN = mean)
#     rownames(expr_subset) <- expr_subset$human_gene
#     expr_subset$human_gene <- NULL
#   } else {
#     rownames(expr_subset) <- human_genes
#   }
#   return(expr_subset)
# }
# 
# # 转换GSE220135_adj数据集的基因名
# GSE220135_adj <- convert_mouse_to_human(GSE220135_adj, homolog_genes)
# 
# # 输出转换结果
# cat("小鼠基因转人类基因转换完成\n")
# cat("转换前基因数量:", length(mouse_genes), "\n")
# cat("转换后基因数量:", nrow(GSE220135_adj), "\n")
# cat("基因转换率:", round(nrow(GSE220135_adj)/length(mouse_genes)*100, 2), "%\n")


# 5.7 使用biomaRt包统一保留mRNA-------------------------------------------------
# 连接Ensembl数据库
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# 获取基因类型信息
gene_attributes <- getBM(attributes = c("hgnc_symbol", "gene_biotype"),
                         mart = ensembl)

# 筛选蛋白编码基因
coding_genes <- gene_attributes$hgnc_symbol[gene_attributes$gene_biotype == "protein_coding"]

# 筛选保留标准化后的表达矩阵中的蛋白编码基因~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
GSE220135_adj <- GSE220135_adj[rownames(GSE220135_adj) %in% coding_genes, ]
# 筛选保留原始count表达矩阵中的蛋白编码基因~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# GSE220135_raw_count_adj <- GSE220135_raw_count_adj[rownames(GSE220135_raw_count_adj) %in% coding_genes, ]


# 5.8 再质控，相当于芯片的标化什么的 -------------------------------------------
# 5.8.1 标准化后的数据集部分~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
GSE220135_adj <- GSE220135_adj

# 负值处理：
# GSE220135_adj <- GSE220135_adj+abs(apply(GSE220135_adj,1,min))
range(GSE220135_adj)

boxplot(GSE220135_adj, outline = F, las = 2) # 再再次查看数据集分布

# 归一化，去除样本间的技术差异
# 去除这一步，高通量测序数据不需要这一步，芯片数据菜肴！！！！！！！
# djusted_counts <- normalizeBetweenArrays(adjusted_counts)

# 去除方差为0的行
GSE220135_adj <- GSE220135_adj[apply(GSE220135_adj, 1, var) != 0,]
GSE220135_adj <- as.data.frame(GSE220135_adj)

# 关闭异常值显示，标签垂直于坐标轴
boxplot(GSE220135_adj, outline = F, las = 2)


# # 5.8.2 原始count的数据部分~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# GSE220135_adj_raw_count <- GSE220135_raw_count_adj
# GSE220135_adj_raw_count <- as.data.frame(GSE220135_adj_raw_count)
# 
# # 去除方差为0的行
# GSE220135_adj_raw_count <- GSE220135_adj_raw_count[apply(GSE220135_adj_raw_count, 1, var) != 0,]
# GSE220135_adj_raw_count <- as.data.frame(GSE220135_adj_raw_count)


## 6.3 整理分组文件-------------------------------------------------------------
# ！！！！！！！！！！一定要检查对照组排序在不在最前面，一定要检查！！！！！
# 这步非常重要，仔细检查！！！！！！！！！
group_2 <- group_2[order(group_2$group,decreasing = T),]#升序排序

# 标准化后的数据部分~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 根据合并数据集的分组信息对合并数据集样本进行重新排列!!!!!!!!!!!!!!!!
GSE220135_adj <- GSE220135_adj[,group_2$ID]
# 再次确认
identical(colnames(GSE220135_adj),group_2$ID)

# # 原始count的数据部分~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # 根据合并数据集的分组信息对合并数据集样本进行重新排列!!!!!!!!!!!!!!!!
# GSE220135_adj_raw_count <- GSE220135_adj_raw_count[,group_2$ID]
# # 再次确认
# identical(colnames(GSE220135_adj_raw_count),group_2$ID)


# 7 输出结果--------------------------------------------------------------------
# 7.1 保存合并数据集结果为csv文件-----------------------------------------------
write.csv(group_2,"output/GSE220135_Group.csv",row.names = F)
write.csv(GSE220135_adj,"output/GSE220135_Matrix_TPM.csv")

# # 7.2 原始count矩阵-------------------------------------------------------------
# write.csv(GSE220135_adj_raw_count,"output/raw_count_GSE220135_Matrix.csv")








