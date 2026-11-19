
# 添加bioconductor清华镜像
options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
# 添加bioconductor西湖大学镜像
options(BioC_mirror = "https://mirrors.westlake.edu.cn/bioconductor")

# 添加R清华镜像
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))


# 2.R包安装与加载===============================================================
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("AUCell")

library(AUCell)#AUcell是分析感兴趣的基因集在所有细胞是否存在富集
library(Seurat)
library(ggpubr)
library(dplyr)


# 3.AUCell======================================================================
# 3.1加载数据-------------------------------------------------------------------
# 加载Final_sc_integer_singleR
Final_sc_integer_singleR <- readRDS("input/Final_sc_integer_singleR.rds")

# 加载基因集
phenogene <- read.csv("input/MEMRGs.csv")  # 修改为自己的数据集
# phenogene <- read.table("input/hubgenes.txt", header = TRUE, sep = "\t")


gene_set <- phenogene$gene



# 3.2分析-----------------------------------------------------------------------
#细胞内基因排序
data <- GetAssayData(Final_sc_integer_singleR[['RNA']], layer = "counts")   #5.0版本需要先把数据提出来
cells_rankings <- AUCell_buildRankings(data,splitByBlocks=TRUE)  #基因排序，使用10个核，加速计算


# 计算每个细胞特定基因集的富集分数
cells_AUC <- AUCell_calcAUC(gene_set, cells_rankings,nCores = 10)


#自动确定阈值,得到特定基因集富集较高的细胞数,无需再手动修改阈值，自动选择Global_k1阈值
set.seed(55)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist = TRUE, assign = TRUE)
auc_threshold <- cells_assignment$geneSet$aucThr$thresholds %>%
  as.data.frame() %>%
  tibble::rownames_to_column("threshold_name") %>%  # 将行名转为列
  filter(threshold_name == "Global_k1") %>%  
  pull(threshold)  # 提取threshold列的值


#根据AUC阈值对细胞进行分类
Final_sc_integer_singleR$AUC <- as.numeric(getAUC(cells_AUC))
Final_sc_integer_singleR$AUCell_Group <- ifelse(Final_sc_integer_singleR$AUC > auc_threshold, 'AUCell_High', 'AUCell_Low')



# 3.3绘制细胞AUC值的UMAP分布图（不分面）----------------------------------------
# 绘制细胞AUC值的UMAP分布图
if ("umap" %in% names(Final_sc_integer_singleR@reductions)) {
  plot.df <- data.frame(Final_sc_integer_singleR@meta.data, Final_sc_integer_singleR@reductions$umap@cell.embeddings)
  colnames(plot.df)[(ncol(plot.df)-1):ncol(plot.df)] <- c("UMAP_1", "UMAP_2")
  reduction <- "UMAP"
  cat("使用UMAP坐标\n")
} else if ("tsne" %in% names(Final_sc_integer_singleR@reductions)) {
  plot.df <- data.frame(Final_sc_integer_singleR@meta.data, Final_sc_integer_singleR@reductions$tsne@cell.embeddings) 
  colnames(plot.df)[(ncol(plot.df)-1):ncol(plot.df)] <- c("TSNE_1", "TSNE_2")
  reduction <- "TSNE"
  cat("使用tSNE坐标\n")
} else {
  stop("未找到UMAP或tSNE降维结果")
}

# 设置坐标列名
x_col <- paste0(reduction, "_1")
y_col <- paste0(reduction, "_2")

# 计算细胞类型中心
class_avg <- plot.df %>%
  dplyr::group_by(cell_anno) %>%
  dplyr::summarise(
    x_center = median(.data[[x_col]]),
    y_center = median(.data[[y_col]])
  )

# 绘图
aucell_plot <- ggplot() +
  geom_point(data = plot.df, 
             aes(x = .data[[x_col]], 
                 y = .data[[y_col]], 
                 colour = AUC), 
             size = 0.6) +
  viridis::scale_color_viridis(option = "C") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text = element_text(colour = 'black', size = 14),
        axis.title = element_text(size = 14),
        axis.line = element_line(colour = 'black'),
        panel.border = element_rect(size = 0.5, linetype = "solid", colour = "black"),
        aspect.ratio = 1) +
  labs(x = paste(reduction, "1"), y = paste(reduction, "2")) +
  ggrepel::geom_label_repel(
    aes(x = x_center, y = y_center, label = cell_anno),
    data = class_avg,
    size = 4,
    label.size = 0,
    segment.color = NA
  )
aucell_plot
ggsave("output/3-AUCell_UMAP.pdf", plot = aucell_plot, scale = 1.2, width = 16.6, height = 12, units =c("cm"))


# 3.3 绘制分面UMAP图（按疾病和对照分面）--------------------------------------
# 计算分面所需的标签位置（按细胞类型和分组同时分组）
class_avg_facet <- plot.df %>%
  dplyr::group_by(cell_anno, tissue) %>%
  dplyr::summarise(
    x_center = median(.data[[x_col]]),
    y_center = median(.data[[y_col]]),
    .groups = "drop"
  )

# 绘制分面UMAP图
aucell_plot_facet <- ggplot() +
  geom_point(data = plot.df, 
             aes(x = .data[[x_col]], 
                 y = .data[[y_col]], 
                 colour = AUC), 
             size = 0.6) +
  viridis::scale_color_viridis(option = "C") +
  # 按分组分面
  facet_wrap(~ tissue) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text = element_text(colour = 'black', size = 14),
        axis.title = element_text(size = 14),
        axis.line = element_line(colour = 'black'),
        panel.border = element_rect(size = 0.5, linetype = "solid", colour = "black"),
        aspect.ratio = 1,
        # 分面标题样式
        strip.background = element_rect(fill = "white", colour = "black"),
        strip.text = element_text(size = 14, colour = "black")) +
  labs(x = paste(reduction, "1"), y = paste(reduction, "2")) +
  ggrepel::geom_label_repel(
    aes(x = x_center, y = y_center, label = cell_anno),
    data = class_avg_facet,
    size = 4,
    label.size = 0,
    segment.color = NA
  )
aucell_plot_facet
ggsave("output/4-AUCell_UMAP_facet.pdf", plot = aucell_plot_facet, scale = 1.2, width = 23, height = 10, units =c("cm"))



# 3.4 绘制细胞的AUC值箱线图（不分面）-------------------------------------------
aucellboxplot <- ggboxplot(plot.df,x = 'cell_anno',y = 'AUC',color = 'cell_anno') +
  stat_compare_means(method = "anova")+
  scale_color_manual(values = c("#C76753", "#CC8593", "#D6B48B", "#A19881",
                                "#718592", "#7989AB", "#B193B7", "#69B3A2",
                                "#ECA869", "#4E79A7", "#F28E2B", "#59A14F",
                                "#8CD17D", "#B07AA1", "#FF9DA7", "#76B7B2"
  )) + 
  xlab('') +
  theme_classic2()+
  theme(legend.position = 'none',
        panel.border = element_rect(size = 0.5, linetype = "solid", colour = "black", fill = NA)) +
  rotate_x_text(45)+
  theme(axis.text.x = element_text(size = 12, colour = "black"))+
  theme(axis.text.y = element_text(size = 12, colour = "black"))
aucellboxplot
ggsave("output/5-AUCell_Boxplot.pdf", plot = aucellboxplot, scale = 1, width = 18, height =13, units =c("cm"))


# 3.4 绘制分面箱线图（按疾病和对照分面）--------------------------------------
# 按分组和细胞类型过滤，仅保留实际存在的组合
plot_filtered <- plot.df %>%
  group_by(tissue, cell_anno) %>%
  filter(n() > 0) %>%
  ungroup()

# 绘制分面箱线图
aucellboxplot_facet <- ggboxplot(
  plot_filtered,
  x = 'cell_anno',
  y = 'AUC',
  color = 'cell_anno',
  facet.by = "tissue",
  palette = c("#C76753", "#CC8593", "#D6B48B", "#A19881",
              "#718592", "#7989AB", "#B193B7", "#69B3A2",
              "#ECA869", "#4E79A7", "#F28E2B", "#59A14F",
              "#8CD17D", "#B07AA1", "#FF9DA7", "#76B7B2"),
  scales = "free_x"
) +
  stat_compare_means(method = "anova") +
  xlab('') +
  theme_classic2() +
  theme(legend.position = 'none',
        panel.border = element_rect(size = 0.5, linetype = "solid", colour = "black", fill = NA)) +
  rotate_x_text(45) +
  theme(axis.text.x = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"),
        strip.background = element_rect(fill = "white", color = "black"),
        strip.text = element_text(size = 14, colour = "black"))
aucellboxplot_facet
ggsave("output/5-AUCell_Boxplot_facet.pdf", plot = aucellboxplot_facet, scale = 1, width = 25, height = 13, units =c("cm"))




# 3.5 总体 AUCell得分较高的细胞类型选取-----------------------------------------
# 查看各个细胞类型的数量
cat(paste0(names(table(Final_sc_integer_singleR$cell_anno)),
           ": ", table(Final_sc_integer_singleR$cell_anno),
           " cells", collapse = "\n"))


# 计算每个细胞类型的AUCell得分（用中位数作为主要指标）
auc_by_celltype <- plot.df %>%
  dplyr::group_by(cell_anno) %>%
  dplyr::summarise(
    mean_AUC = mean(AUC, na.rm = TRUE),
    median_AUC = median(AUC, na.rm = TRUE),
    sd_AUC = sd(AUC, na.rm = TRUE),
    cell_count = n()
  ) %>%
  dplyr::arrange(desc(median_AUC))  # 改为按中位数降序排列

# 输出AUCell评分排名
cat("\n========== AUCell评分排名（按中位数降序） ==========\n")
for (i in 1:nrow(auc_by_celltype)) {
  cat(sprintf("Rank %d: %s | Median AUC = %.4f | Mean AUC = %.4f | Cell Count = %d\n",
              i, 
              auc_by_celltype$cell_anno[i],
              auc_by_celltype$median_AUC[i],
              auc_by_celltype$mean_AUC[i],
              auc_by_celltype$cell_count[i]))
}

# 输出评分最高的细胞簇（基于中位数）
top_cluster <- auc_by_celltype$cell_anno[1]
cat("\n========== AUCell评分最高的细胞簇（基于中位数） ==========\n")
cat(sprintf("Top 1: %s (Median AUC = %.4f)\n", top_cluster, auc_by_celltype$median_AUC[1]))


# 得出评分最高的：
# Top 1: AT2 (Median AUC = 0.0928)




