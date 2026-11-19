# -----------------------------------------------------------------------------#
# 描述: 本脚本用于执行单细胞scPagwas分析 scPagwas.R
#       此代码的单细胞数据是小鼠的 GWAS数据是人类的,跨物种版（单细胞小鼠 + GWAS人类)
# R版本：4.3.2
# Seurat (Version:5.4.1)
# 日期: 2026/07/07
# 优化修改：阿伦
# 代码版本：V1 
# 
# 更新内容：
# 2026.7.7 优化了注释结构
# 2026.7.20 添加了降采样
# 
# 存在问题：会R包加载冲突，导致分析完后的数据再次打开无法跑，得删除所有环境重新加载
# -----------------------------------------------------------------------------#


# 配置动态R包库,我服务器安装R包没权限，得设置动态R包库(服务器）=================
# 我的个人R包库在/work/run/projects/bio-35/R_library
.libPaths(c("/work/run/projects/bio-35/R_library", .libPaths()))



# 0.准备工作====================================================================
library(devtools)
# install_git("https://github.com/WMU-SuLab/scPagwas.git", ref = "main")

# 加载包
library(scPagwas)
library(Seurat)
library(ggplot2)
library(tidyverse)
library(reshape2)
library(gridExtra)
library(magrittr)
library(data.table)
library(dplyr)
library(KEGGREST)
library(purrr)
library(cli)
library(gwasglue)
# library(VariantAnnotation)


# 1. 加载数据===================================================================
# 1.1 单细胞数据----------------------------------------------------------------
# 含有AT2亚群的总的单细胞数据
scdata <- readRDS("input/sc_integer_all_with_AT2_CytoTRACE.rds")


# ---------------------------可选步骤（降采样）--------------------------------#

scdata <- SketchData(
  object = scdata,
  assay = "RNA",           # 使用RNA assay（或SCT等）
  ncells = 5000,           # 采样5000个细胞
  method = 'Uniform',      # 均匀采样方法
  sketched.assay = 'sketch_uniform',  # 创建新的sketch assay
  seed = 123               # 设置随机种子确保可重现
)

# 提取降采样结果并赋值到新变量
downsample_cells <- Cells(scdata[["sketch_uniform"]])
scdata <- subset(scdata, cells = downsample_cells)

# 使用均匀采样后需要重新设置默认assay为RNA
DefaultAssay(scdata) <- "RNA"

# ---------------------------可选步骤（降采样）--------------------------------#



SeuratObject::Version(scdata) # 确认版本>=4.0.0
Idents(scdata)[1:10] # 确认为细胞类型作标识
Idents(scdata) <- "cell_anno_CytoTRACE"   # 定义后续分析的细胞身份，仔细注意这个


# 1.2 GWAS数据------------------------------------------------------------------
# 读取教程里方法1下载的数据
# 要看教程，不同疾病对应不同GWAS数据
gwas_data <- data.table::fread("input/finngen_R12_P16_RESPI_DISTR_NEWBO.gz", header = T)
head(gwas_data)

# GWAS数据格式必须包含的内容
# chrom | pos | rsid | beta | se | maf


# 数据格式整理
gwas_data <- gwas_data[,c('#chrom', 'pos', 'rsids', 'beta', 'sebeta', 'af_alt')]
gwas_data[, `:=` (chrom = `#chrom`,
                  rsid = rsids,
                  se = sebeta,
                  maf = pmin(af_alt, 1-af_alt))]
gwas_data <- gwas_data[,.(chrom, pos, rsid, beta, se, maf)] %>% as.data.frame()
head(gwas_data)

# 保存处理好格式后的GWAS数据
# saveRDS(gwas_data,"input/gwas_data.RDS")
gwas_data <- readRDS("input/gwas_data.RDS")



# 2 KEGG通路基因集==============================================================
# 请选择以下方法的1种

# 2.1 scPagwas 内置KEGG通路数据-------------------------------------------------
# data(Genes_by_pathway_kegg)


# 2.2 最新/指定版本KEGG通路数据集-----------------------------------------------
# 即自定义KEGG文件,根据物种不同选择下载的脚本（人类）

# # ---------------获取人类 homo sapiens 的最新 KEGG 通路列表--------------------#
# # 获取通路名称
# 
# dir.create("KEGG_pathway_names", showWarnings = FALSE)
# # 设置代理
# Sys.setenv(
#   http_proxy = "http://localhost:1086",
#   https_proxy = "http://localhost:1086"
# )
# 
# pathways.list <- keggList("pathway", "hsa")
# cli::cli_alert_info("找到 {length(pathways.list)} 个通路")
# 
# # 分批处理函数（避免请求过多）
# source("utils/get_all_pathway_genes_safely.R")
# 
# all_genes <- get_all_pathway_genes_safely("hsa", batch_size = 5, delay = 0.3)
# 
# # 将 ID 转换为 Symbol
# all_genes %<>%
#   mutate(gene_symbol = str_extract(description, "^[^;]+"))
# 
# # 生成一个列表，其中包含通路名称及其对应的基因符号
# Genes_by_pathway_kegg <- all_genes %>%
#   group_by(pathway_name) %>%
#   summarise(genes = list(gene_symbol)) %>%
#   deframe()
# 
# # 创建一个表格，包含通路 ID（pwy id）、通路名称（name），以及用分号（;）拼接在一起的基因符号
# collap <-  all_genes %>%
#   group_by(pathway_id, pathway_name) %>%
#   summarise(genes = paste(gene_symbol, collapse = ";")) %>%
#   ungroup() %>%
#   as.data.frame()
# 
# data.table::fwrite(collap,
#   file = str_glue("KEGG_pathway_names/Homo_KEGG_pathway_names_{Sys.Date()}.tsv"),
#   sep = "\t")
# 
# # 最新的KEGG通路基因集以当天日期为后缀
# # 关闭代理
# Sys.unsetenv(c("http_proxy", "https_proxy"))
# # ---------------获取人类 homo sapiens 的最新 KEGG 通路列表--------------------#



# 2.2.1 处理下载好的KEGG通路文件~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 导入人类版
kegg_pathway_file <- "KEGG_pathway_names/Homo_KEGG_pathway_names_2026-07-07.tsv" # 按下载好的文件日期修改名字
kegg_pathway_data <- data.table::fread(kegg_pathway_file, header = T)

# 转换格式的KEGG通路文件
Genes_by_pathway_kegg <- kegg_pathway_data %>%
  mutate(genes = strsplit(genes, ";")) %>%
  dplyr::select(pathway_id, genes) %>%     # 或者pathway_name
  deframe()


# -------------------------小鼠基因→人类同源基因转换---------------------------#

# 由于 scPagwas 的基因组/LD/通路数据均为人类，单细胞数据为小鼠，
# 需将小鼠基因名转换为其人类同源基因，使 GWAS 信号能正确匹配
library(homologene)

# 获取 Seurat 对象中的小鼠基因
mouse_genes <- rownames(scdata)

# 使用 homologene 进行鼠→人同源转换（10090=小鼠，9606=人类）
homologene_map <- homologene(mouse_genes, inTax = 10090, outTax = 9606)

# 整理映射表：小鼠基因 → 人类基因，保留 1:1 映射
orthologs <- homologene_map[, c(1, 2)] %>%
  setNames(c("mouse_gene", "human_gene")) %>%
  filter(human_gene != "") %>%
  filter(!duplicated(mouse_gene)) %>% # 一个小鼠基因只对应一个人类基因
  filter(!duplicated(human_gene))     # 一个人类基因只来自一个小鼠基因

cat(sprintf("小鼠基因总数: %d → 成功映射的人类同源基因: %d\n",
            length(mouse_genes), nrow(orthologs)))

# 步骤1：过滤——只保留有 1:1 人类同源基因的基因（此时基因名仍为小鼠名）
scdata <- scdata[orthologs$mouse_gene, ]

# 步骤2：重命名——将 counts 矩阵中的小鼠基因名改为人类基因名
counts <- GetAssayData(scdata, assay = "RNA", layer = "counts")
rownames(counts) <- orthologs$human_gene

# 重建 RNA assay（Assay5 需通过重建才能正确更新 feature 名）
scdata[["RNA"]] <- CreateAssay5Object(counts = counts)

# 重新归一化，生成 data 层（重建后 data 不存在，必须重新计算）
scdata <- NormalizeData(scdata, verbose = FALSE)

rownames(scdata)

# -------------------------小鼠基因→人类同源基因转换---------------------------#



# 3 运行scPagwas================================================================
# 3.1 跨物种版（单细胞小鼠 + GWAS人类）-----------------------------------------

# --------------------------------注意事项-------------------------------------#
# 你可小心点,小心踩坑
# block_annotation = block_annotation_hg37 要根据GWAS的数据选择对应的参考基因组数据
# 如：finngen_R12_P16_RESPI_DISTR_NEWBO.gz 是hg38,不能用hg37

# block_annotation 参数解释
# block_annotation      = hg38 基因组
# block_annotation_hg37 = hg37 基因组
# --------------------------------注意事项-------------------------------------#

# 查看scPagwas包有什么内置数据集
data(package = "scPagwas")

# 设置随机种子
set.seed(2026)
Pagwas_data <- scPagwas_main(Pagwas = NULL,
                             gwas_data = gwas_data ,         # GWAS 汇总统计数据（人类）
                             Single_data = scdata,           # 已标准化的 Seurat 单细胞对象（基因已转为人类同源基因）
                             output.prefix = "real",         # 输出文件前缀名
                             output.dirs = "scPagwastest_output", # 结果输出文件夹
                             block_annotation = block_annotation, # 基因组注释信息（根据GWAS版本选择hg37或hg38）
                             assay = "RNA",                  # 使用的表达矩阵类型
                             Pathway_list = Genes_by_pathway_kegg, # 通路基因列表（人类 KEGG）
                             n.cores = 50,                   # 并行计算核心数
                             iters_singlecell = 100,         # 单细胞分析的bootstrap迭代次数
                             chrom_ld = chrom_ld,            # 连锁不平衡（LD）参考数据
                             singlecell = TRUE,              # 是否进行单细胞层面的分析
                             celltype = TRUE,                # 是否进行细胞类型层面的分析
                             n_topgenes = 1000,              # 用于计算TRS的高变基因数量
                             min.pathway.size = 5,           # 通路最少基因数
                             max.pathway.size = 300,         # 通路最多基因数
                             maf_filter = 0.01               # SNP的最小等位基因频率过滤阈值
)


saveRDS(Pagwas_data, file = "Pagwas_data_OS.rds")
# Pagwas_data <- readRDS("Pagwas_data_OS.rds")



# 4. 结果可视化=================================================================
Idents(Pagwas_data) <- "cell_anno_CytoTRACE"   # 定义后续分析的细胞身份

# 4.1 绘制得分图----------------------------------------------------------------
scPagwas_Visualization(
  Single_data = Pagwas_data,
  p_thre = 0.05,       # p阈值
  FigureType = "umap",
  width = 6,
  height = 6,
  lowColor = "white",
  highColor = "red",
  output.dirs = "figure",
  size = 0.5,          # 点大小参数
  do_plot = TRUE
)

# 4.2 绘制Bootstrap柱状结果图---------------------------------------------------
cairo_pdf("figure/2-bootstrap_barplot.pdf",
          width = 6,
          height = 4)
Bootstrap_P_Barplot(Pagwas_data@misc$bootstrap_results$bp_value[-1], 
                    p_names = rownames(Pagwas_data@misc$bootstrap_results)[-1])
dev.off()


# 4.3 绘制Bootstrap结果森林图---------------------------------------------------
cairo_pdf("figure/3-bootstrap_forestplot.pdf", width = 7, height = 5.5)

Bootstrap_estimate_Plot(
  Pagwas_data@misc$bootstrap_results,
  figurenames = NULL,
  width = 9,
  height = 7,
  do_plot = TRUE
)

dev.off()


# 4.4 绘制基因遗传力gPAS与基因表达相关性散点图----------------------------------
cairo_pdf("figure/4-gene_heritability_scatter.pdf", width = 5, height = 5)

heritability_cor_scatterplot(
  gene_heri_cor = Pagwas_data@misc$PCC,
  topn_genes_label = 10,
  color_low = "#035397",
  color_high = "#F32424",
  color_mid = "white",
  text_size = 3.5,
  max.overlaps = 20,
  do_plot = TRUE
)

dev.off()


# 4.5 TRS 按细胞类型分布图------------------------------------------------------
trs_data <- data.frame(
  TRS = Pagwas_data$scPagwas.TRS.Score1,
  celltype = Idents(Pagwas_data)
)

p_trs <- ggplot(trs_data, aes(x = celltype, y = TRS, fill = celltype)) +
  geom_boxplot(outlier.size = 0.5) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(fill = NA, color = "black")) +
  labs(title = "TRS Score Distribution Across Cell Types",
       x = "Cell Type",
       y = "TRS Score") +
  guides(fill = "none")

cairo_pdf("figure/5-TRS_by_celltype.pdf", width = 5, height = 4)
print(p_trs)
dev.off()


# 4.6 通路贡献 Dot 图----------------------------------------------------------
library(tidyverse)
library("rhdf5")
library(grDevices)
library(stats)
library(FactoMineR)
library(scales)
library(reshape2)
library(ggdendro)
library(grImport2)
library(gridExtra)
library(grid)
library(sisal)

source(system.file("extdata", 
                   "plot_scpathway_contri_dot.R", 
                   package = "scPagwas"))

dcast <- reshape2::dcast  # 强制使用 reshape2::dcast
select <- dplyr::select

cairo_pdf("figure/6-pathway_contribution_dotplot.pdf", width = 9, height = 6)

plot_scpathway_dot(
  Pagwas = Pagwas_data,
  celltypes = unique(Idents(Pagwas_data)),
  topn_path_celltype = 3,  # 每个细胞类型最多展示前 5 个最显著的通路
  filter_p = 0.05,
  max_logp = 15,
  display_max_sizes = FALSE,
  size_var = "logrankPvalue",
  col_var = "proportion",
  shape.scale = 8,
  cols.use = c("lightgrey", "#E45826"),
  dend_x_var = "logrankPvalue",
  dist_method = "euclidean",
  hclust_method = "ward.D",
  do_plot = TRUE
)

dev.off()


# # TSMR分析细胞亚群特征基因与性状风险之间的因果关系==============================
# # 提取目的细胞群的特征基因
# # 从TRS得分中提取重要基因
# trs_genes <- Pagwas_data@misc$PCC %>%
#   as.data.frame() %>%
#   mutate(abs_cor = abs(.[[1]])) %>%  # 使用 [[1]] 提取第一列的数值向量
#   arrange(desc(abs_cor)) %>%
#   filter(abs_cor >= quantile(abs_cor, 0.9)) %>% # 取前10 % 重要性的基因
#   rownames()
# 
# unique(Idents(scdata))
# 
# CD8Tem_markers <- FindMarkers(
#   scdata,
#   ident.1 = "CD8_T:Tem",        # 目标细胞类型：CD8+ Tem
#   ident.2 = NULL,               # 对比所有其他细胞类型
#   assay = "RNA",
#   slot = "data",
#   only.pos = TRUE,              # 只保留上调基因
#   min.pct = 0.25,               # 至少在25%的CD8_Tem细胞中表达
#   logfc.threshold = 1,          # 适当降低阈值，保留更多相关基因
#   test.use = "wilcox"
# )
# 
# IPFR_CD8Tem_Gs <- intersect(trs_genes,rownames(CD8Tem_markers))
# 
# saveRDS(IPFR_CD8Tem_Gs,"output/IPFR_CD8Tem_Gs.RDS")
# write.csv(IPFR_CD8Tem_Gs,"output/IPFR_CD8Tem_Gs.csv")
# 
# IPFR_CD8Tem_Gs <- read.csv("output/IPFR_CD8Tem_Gs.csv")

# 保存工作空间
save.image(file = "scPagwas.RData")

