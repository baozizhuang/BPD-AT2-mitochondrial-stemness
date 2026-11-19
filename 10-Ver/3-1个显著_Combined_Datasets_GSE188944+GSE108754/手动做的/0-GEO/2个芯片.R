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

# 1.2 R包安装-------------------------------------------------------------------
# 定义一个函数来检查并安装所需的R包
install_if_missing <- function(packages) {
  installed_packages <- rownames(installed.packages())
  for (pkg in packages) {
    if (!(pkg %in% installed_packages)) {
      BiocManager::install(pkg, dependencies = TRUE, update = FALSE)
    }
  }
}

# 检查安装的包列表
required_packages <- c("GEOquery","dplyr","magrittr","ggplot2","sva","limma",
                       "FactoMineR","factoextra")
install.packages("ggplot2")

# 调用函数进行检查和安装
install_if_missing(required_packages)

BiocManager::install("genefilter")
BiocManager::install("Biostrings")
BiocManager::install("GEOquery")

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
library(homologene)
library(clusterProfiler)
library(org.Hs.eg.db)

# 1.4 创建文件夹---------------------------------------------------------------
dir.create("input/")
dir.create("output/")


# 1.5 定义变量-----------------------------------------------------------------
# 测试集
GSE_1 <- "GSE188944"   # 修改：GEO数据集
GSE_2 <- "GSE108754"

treat_name = "BPD"    # 修改：疾病组命名
con_name = "Control"  # 对照组命名


# 2 获取数据--------------------------------------------------------------------
# 2.1 下载原始矩阵数据----------------------------------------------------------
# 先手动下载到input
# 然后本地获取测试集芯片数据
# 测试集
aset_1 <- getGEO(GSE_1,destdir = "input", AnnotGPL = F, getGPL = F)
aset_2 <- getGEO(GSE_2,destdir = "input", AnnotGPL = F, getGPL = F)


# 2.2 定义所需函数--------------------------------------------------------------
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


## 2.3 数据清洗-----------------------------------------------------------------
# 测试集的2个芯片数据处理
# 第1个芯片表达谱的处理GSE188944
GSE_matrix_1 <- as.data.frame(exprs(aset_1$GSE188944_series_matrix.txt.gz))
range(GSE_matrix_1)
#移除所有包含至少一个NA值的行
GSE_matrix_1 <- na.omit(GSE_matrix_1)
#处理负值数据
GSE_matrix_1 <- GSE_matrix_1 + abs(apply(GSE_matrix_1,1,min))
range(GSE_matrix_1)
#log2化
GSE_matrix_1 <- log2_if_need(GSE_matrix_1)

# 第2个芯片表达谱的处理GSE108754
GSE_matrix_2 <- as.data.frame(exprs(aset_2$GSE108754_series_matrix.txt.gz))
range(GSE_matrix_2)
#移除所有包含至少一个NA值的行
GSE_matrix_2 <- na.omit(GSE_matrix_2)
#处理负值数据
# GSE_matrix_2 <- GSE_matrix_2 + abs(apply(GSE_matrix_2,1,min))
range(GSE_matrix_2)
GSE_matrix_2 <- log2_if_need(GSE_matrix_2)


# 2.4 整理分组信息------------------------------------------------------------
# 测试集的2个芯片数据分组信息处理
PD_1 <- pData(aset_1[["GSE188944_series_matrix.txt.gz"]])
  
# 如果有研究疾病之外的其他疾病样本，运行下面一行代码去除(此处有)
# PD_1 <- PD_1[!grepl("liver_SS_CL", PD_1$title), ]  # 修改：去除无关样本
# PD_1 <- PD_1[!grepl("liver_SS_HLD", PD_1$title), ]  # 修改：去除无关样本


group_1 <- geo_group(pdata = PD_1,
                     # 样本编号列名；分组信息列名（修改）
                     sample_id = "geo_accession", group_id = "title", 
                     # 疾病组命名；对照组命名；是否降序可修改）
                     treat_name = treat_name, con_name = con_name, desc = T, 
                     # 分组信息列名的疾病组关键字符串；对照组关键字符串（修改）
                     treat_key = "^BPD", con_key = "^nonBPD") 


PD_2 <- pData(aset_2$GSE108754_series_matrix.txt.gz)
# PD_2 <- PD_2[!grepl("Healthy obese|Steatosis", PD_2$source_name_ch1), ] # 修改：去除无关样本
# PD_2 <- PD_2[!grepl("surgery", PD_2$source_name_ch1), ] # 修改：去除无关样本

group_2 <- geo_group(pdata = PD_2,
                     # 样本编号列名；分组信息列名（修改）
                     sample_id = "geo_accession", group_id = "title", 
                     # 疾病组命名；对照组命名；是否降序（修改）
                     treat_name = treat_name, con_name = con_name, desc = T, 
                     # 分组信息列名的疾病组关键字符串；对照组关键字符串（修改）
                     treat_key = "BPD", con_key = "Normal") 


# 2.5 探针转换,同一基因取均值--------------------------------------------------
# 测试集探针替换
GPL_1 <- Table(getGEO(unique(PD_1$"platform_id"), 
                      destdir = "input", AnnotGPL = F))

# 探针和基因ID视情况修改为相应名称
gpl_1 <- gpl_pro(gpl_name = GPL_1, probe_id = "ID", gene_id = "Symbol")#修改基因信息
GSE_matrix_processed_1 <- geo_matrix(data_matrix = GSE_matrix_1, 
                                     data_group = group_1, gpl_name = gpl_1)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
GPL_2 <- Table(getGEO(unique(PD_2$"platform_id"), 
                      destdir = "input", AnnotGPL = F))

gpl_2 <- gpl_pro(gpl_name = GPL_2, probe_id = "ID", gene_id = "GENE_SYMBOL")
GSE_matrix_processed_2 <- geo_matrix(data_matrix = GSE_matrix_2, 
                                     data_group = group_2, gpl_name = gpl_2)


# 2.6 将训练集表达矩阵的基因转换------------------------------------------------
# # 2.6.1 将训练集表达矩阵的大鼠基因转成人类基因~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # 使用homologene包进行基因转换
# # inTax = 10116: 大鼠, outTax = 9606: 人类
# 
# # 获取两个数据集的所有基因并去重
# rat_genes <- unique(c(rownames(GSE_matrix_processed_1),
#                       rownames(GSE_matrix_processed_2)))
# 
# # 进行基因转换
# homolog_genes <- homologene(rat_genes, inTax = 10116, outTax = 9606)
# 
# # 提取基因映射关系(第1列是大鼠基因,第2列是人类基因)并去重
# homolog_genes <- homolog_genes[, c(1, 2)]
# colnames(homolog_genes) <- c("rat_gene", "human_gene")
# homolog_genes <- homolog_genes[!duplicated(homolog_genes$rat_gene), ]
# rownames(homolog_genes) <- homolog_genes$rat_gene
# 
# # 定义基因转换函数
# convert_rat_to_human <- function(expr_matrix, mapping) {
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
# # 转换三个数据集的基因名
# GSE_matrix_processed_1 <- convert_rat_to_human(GSE_matrix_processed_1, homolog_genes)
# GSE_matrix_processed_2 <- convert_rat_to_human(GSE_matrix_processed_2, homolog_genes)


# # 2.6.2 将训练集的表达矩阵的大鼠基因转成小鼠基因~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # 使用homologene包进行基因转换
# # inTax = 10116: 大鼠, outTax = 10090: 小鼠
# 
# # 获取两个数据集的所有基因并去重
# rat_genes <- unique(c(rownames(GSE_matrix_processed_1),
#                       rownames(GSE_matrix_processed_2)))
# 
# # 进行基因转换
# homolog_genes <- homologene(rat_genes, inTax = 10116, outTax = 10090)
# 
# # 提取基因映射关系(第1列是大鼠基因,第2列是小鼠基因)并去重
# homolog_genes <- homolog_genes[, c(1, 2)]
# colnames(homolog_genes) <- c("rat_gene", "mouse_gene")
# homolog_genes <- homolog_genes[!duplicated(homolog_genes$rat_gene), ]
# rownames(homolog_genes) <- homolog_genes$rat_gene
# 
# # 定义基因转换函数
# convert_rat_to_mouse <- function(expr_matrix, mapping) {
#   common_genes <- intersect(rownames(expr_matrix), rownames(mapping))
#   expr_subset <- expr_matrix[common_genes, , drop = FALSE]
#   mouse_genes <- mapping[common_genes, "mouse_gene"]
# 
#   # 如果有重复的小鼠基因，取平均值
#   if (any(duplicated(mouse_genes))) {
#     # 使用 aggregate 对重复基因取平均值
#     expr_subset$mouse_gene <- mouse_genes
#     expr_subset <- aggregate(. ~ mouse_gene, data = expr_subset, FUN = mean)
#     rownames(expr_subset) <- expr_subset$mouse_gene
#     expr_subset$mouse_gene <- NULL
#   } else {
#     rownames(expr_subset) <- mouse_genes
#   }
#   return(expr_subset)
# }
# 
# # 转换三个数据集的基因名
# GSE_matrix_processed_1 <- convert_rat_to_mouse(GSE_matrix_processed_1, homolog_genes)
# GSE_matrix_processed_2 <- convert_rat_to_mouse(GSE_matrix_processed_2, homolog_genes)


# # 2.6.3 将训练集的表达矩阵的人类基因转成小鼠基因~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # 使用homologene包进行基因转换
# # inTax = 9606: 人类, outTax = 10090: 小鼠
# 
# # 获取两个数据集的所有基因并去重
# human_genes <- unique(c(rownames(GSE_matrix_processed_1),
#                       rownames(GSE_matrix_processed_2)))
# 
# # 进行基因转换
# homolog_genes_human_mouse <- homologene(human_genes, inTax = 9606, outTax = 10090)
# 
# # 提取基因映射关系(第1列是人类基因,第2列是小鼠基因)并去重
# homolog_genes_human_mouse <- homolog_genes_human_mouse[, c(1, 2)]
# colnames(homolog_genes_human_mouse) <- c("human_gene", "mouse_gene")
# homolog_genes_human_mouse <- homolog_genes_human_mouse[!duplicated(homolog_genes_human_mouse$human_gene), ]
# rownames(homolog_genes_human_mouse) <- homolog_genes_human_mouse$human_gene
# 
# # 定义基因转换函数
# convert_human_to_mouse <- function(expr_matrix, mapping) {
#   common_genes <- intersect(rownames(expr_matrix), rownames(mapping))
#   expr_subset <- expr_matrix[common_genes, , drop = FALSE]
#   mouse_genes <- mapping[common_genes, "mouse_gene"]
# 
#   # 如果有重复的小鼠基因，取平均值
#   if (any(duplicated(mouse_genes))) {
#     # 使用 aggregate 对重复基因取平均值
#     expr_subset$mouse_gene <- mouse_genes
#     expr_subset <- aggregate(. ~ mouse_gene, data = expr_subset, FUN = mean)
#     rownames(expr_subset) <- expr_subset$mouse_gene
#     expr_subset$mouse_gene <- NULL
#   } else {
#     rownames(expr_subset) <- mouse_genes
#   }
#   return(expr_subset)
# }
# 
# # 转换两个数据集的基因名
# GSE_matrix_processed_1 <- convert_human_to_mouse(GSE_matrix_processed_1, homolog_genes_human_mouse)
# GSE_matrix_processed_2 <- convert_human_to_mouse(GSE_matrix_processed_2, homolog_genes_human_mouse)


# 2.7 使用biomaRt包统一保留mRNA-------------------------------------------------
# 人类版本~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 连接Ensembl数据库
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# 获取基因类型信息
gene_attributes <- getBM(attributes = c("hgnc_symbol", "gene_biotype"),
                         mart = ensembl)

# 筛选蛋白编码基因
coding_genes <- gene_attributes$hgnc_symbol[gene_attributes$gene_biotype == "protein_coding"]

# 筛选保留表达矩阵中的蛋白编码基因(3组数据)
GSE_matrix_processed_1 <- GSE_matrix_processed_1[rownames(GSE_matrix_processed_1) %in% coding_genes, ]
GSE_matrix_processed_2 <- GSE_matrix_processed_2[rownames(GSE_matrix_processed_2) %in% coding_genes, ]


# # 小鼠版本~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # 连接Ensembl数据库（小鼠基因组）
# ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
# 
# # 获取基因类型信息（使用小鼠基因符号）
# gene_attributes <- getBM(attributes = c("mgi_symbol", "gene_biotype"),
#                          mart = ensembl)
# 
# # 筛选蛋白编码基因
# coding_genes <- gene_attributes$mgi_symbol[gene_attributes$gene_biotype == "protein_coding"]
# 
# # 筛选保留表达矩阵中的蛋白编码基因(2组数据)
# GSE_matrix_processed_1 <- GSE_matrix_processed_1[rownames(GSE_matrix_processed_1) %in% coding_genes, ]
# GSE_matrix_processed_2 <- GSE_matrix_processed_2[rownames(GSE_matrix_processed_2) %in% coding_genes, ]



# 3 数据集合并（2个）-----------------------------------------------------------
# 3.1 保留数据集中共有基因的数据部分--------------------------------------------
# 找到2个芯片表达谱中共同的基因
same_gene <- intersect(rownames(GSE_matrix_processed_1), 
                       rownames(GSE_matrix_processed_2))
comGSE_1 <- GSE_matrix_processed_1[same_gene,]
comGSE_2 <- GSE_matrix_processed_2[same_gene,]

# 检查行名是否完全相同
identical(rownames(comGSE_1), rownames(comGSE_2))

# 合并矩阵数据
count_matrix <- cbind(comGSE_1, comGSE_2) %>% na.omit()

# 查看数据集间是否有批次
boxplot(count_matrix, outline = F, las = 2) 


## 3.2 去除批次效应-------------------------------------------------------------
## "1"代表 data1,"2"代表 data2
batch <- c(rep("1", length(group_1$group)),
           rep("2", length(group_2$group)))

adjusted_counts <- ComBat(count_matrix, batch = batch) 
boxplot(adjusted_counts, outline = F, las = 2) # 再次查看数据集间是否有批次
range(adjusted_counts)


##负值处理：
# adjusted_counts <- adjusted_counts+abs(apply(adjusted_counts,1,min))
range(adjusted_counts)
boxplot(adjusted_counts, outline = F, las = 2)

# 标准化
adjusted_counts <- normalizeBetweenArrays(adjusted_counts)

# 去除方差为0的行
adjusted_counts <- adjusted_counts[apply(adjusted_counts, 1, var) != 0,]
adjusted_counts <- as.data.frame(adjusted_counts)

# 关闭异常值显示，标签垂直于坐标轴
boxplot(adjusted_counts, outline = F, las = 2)



## 3.3 整理合并后的分组文件-----------------------------------------------------
combined_pd <- rbind(group_1,group_2)

# ！！！！！！！！！！一定要检查对照组排序,确保在最前面，一定要检查！！！！！
# 这步非常重要，仔细检查！！！！！！！！！
# adjusted_pd <- combined_pd %>%
#   arrange(group)

# 或者(第二种排序方法)
adjusted_pd <- combined_pd[order(combined_pd$group, decreasing = T),]


# 根据合并数据集的分组信息对合并数据集样本进行重新排列!!!!!!!!!!!!!!!!
adjusted_counts <- adjusted_counts[,adjusted_pd$ID]
# 再次确认
identical(colnames(adjusted_counts),adjusted_pd$ID)


# 3 绘制校正前后boxplot,PCA图---------------------------------------------------
## 3.1 定义绘图函数-------------------------------------------------------------
# mytheme,义了一个自定义的 ggplot2 主题 mytheme，用于设置图表的外观和样式。
mytheme <- 
  theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
        plot.title = element_text(size = 7)) + 
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = NA, size = 0.75 * 0.47)) +
  theme(axis.line = element_line(size = 0.75 * 0.47),
        axis.text = element_text(size = 6, color = "black"), 
        axis.title = element_text(size = 6),
        axis.ticks = element_line(size = 0.75 * 0.47)) +
  theme(legend.key = element_rect(color = NA,fill = NA),
        legend.key.size = unit(c(0.3, 0.3), "cm"),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.margin = margin(),
        legend.box.margin = margin(),
        legend.box.spacing = unit(0, "cm"),
        legend.background = element_blank(), legend.spacing = unit(0, "cm"),
        legend.box.background = element_blank())


# 转长数据,自定义函数
dat_to_long <- function(raw_dat,combined_pd,cla){
  raw_dat <- cbind(rownames(raw_dat),raw_dat)
  rownames(raw_dat) <- NULL
  colnames(raw_dat)[1] <- "Gene"
  rownames(combined_pd) <- NULL
  raw_dat<-raw_dat[,-1]
  
  dat_long <- raw_dat %>%
    t() %>%
    data.frame() %>%
    dplyr::mutate(sample = rownames(.)) %>%
    dplyr::mutate(group = cla) %>%
    dplyr::mutate(sample = factor(.$sample, levels = .$sample)) %>%
    tidyr::gather(key = geneid, value, - c(sample, group))
  dat_long$sample <- factor(dat_long$sample,
                            levels = as.vector(unique(dat_long$sample)))
  return(dat_long)
}

# Boxplot绘图，自定义函数
dat_boxplot <- function(dat_long,GSEID,desc,group_col,boxplot_width,boxplot_high)
{
  boxplot <- 
    ggplot(dat_long, aes(sample, value)) + 
    geom_boxplot(aes(fill = group), lwd = 0.1 *0.47, outlier.shape = NA) + 
    scale_fill_manual(values = group_col) +
    labs(title = paste0(GSEID," ",desc," Normalization")) + 
    mytheme + # 自己的主题
    theme(legend.position = 'top', 
          legend.direction = "horizontal", 
          legend.title = element_blank(), 
          axis.title.x = element_blank(), 
          axis.title.y = element_blank(), 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank()) + 
    scale_y_continuous(limits = boxplot.stats(dat_long$value)$stats[c(1, 5)]*1.05,
                       expand = c(0,0)) +
    coord_cartesian(clip = "off")
  boxplot
  ggsave(paste0("output/", GSEID, " Boxplot ", desc, ".pdf"), 
         plot = boxplot, units = "cm", 
         width = boxplot_width, height = boxplot_high)
}     

# PCA绘图
pca_plot <- function(data,GSEID,desc,group_col,pca_width,pca_height,cla){
  data <- as.matrix(t(data))
  iris.pca<- PCA(data,graph =F,scale.unit = TRUE)
  pdf(paste0("output/", GSEID," PCA ", desc, ".pdf"), 
      width = pca_width, height = pca_height)
  
  ind.p <- fviz_pca_ind(iris.pca, geom.ind = "point", col.ind = cla,
                        palette = group_col,
                        addEllipses = TRUE,
                        legend.title = "Group",
                        ggtheme = theme_bw() +  # 使用theme_bw为基础
                          theme(panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank()))
  
  ind.p <- ggpubr::ggpar(ind.p, title = paste(GSEID,desc, "Normalization"))
  print(ind.p)
  dev.off()
}

# 3.2 绘图----------------------------------------------------------------------
# 画合并数据集校正前后的boxplot图-----------------------------------------------
group_col_1 = c("#FDD692","#A593E0")  # 修改：自定义绘图分组颜色
cla <- c(rep(GSE_1, length(group_1$group)),
         rep(GSE_2, length(group_2$group)))


dat_before_long <- dat_to_long(count_matrix,combined_pd,cla = cla)
boxplot_before <- dat_boxplot(dat_before_long,GSEID = "Combined Datasets",
                              desc = "Before", group_col_1,
                              boxplot_width=14.5, boxplot_high=5)

dat_after_long <- dat_to_long(adjusted_counts,adjusted_pd,cla)
boxplot_after <- dat_boxplot(dat_after_long,GSEID = "Combined Datasets",
                             desc = "After",group_col_1,
                             boxplot_width=14.5, boxplot_high=5)


## 画合并数据集校正前后的PCA图--------------------------------------------------
pca_before <- pca_plot(count_matrix,GSEID = "Combined Datasets",
                       desc = "Before",group_col=group_col_1,
                       pca_width = 6, pca_height = 4.5,cla=cla)
pca_after <- pca_plot(adjusted_counts,GSEID = "Combined Datasets",
                      desc = "After",group_col=group_col_1,
                      pca_width = 6, pca_height = 4.5,cla=cla)



# 4 输出结果--------------------------------------------------------------------
# 4.1 保存为csv文件------------------------------------------------------------
write.csv(adjusted_pd,"output/Combined_Datasets_Group.csv",row.names = F)
write.csv(adjusted_counts,"output/Combined_Datasets_Matrix.csv")


