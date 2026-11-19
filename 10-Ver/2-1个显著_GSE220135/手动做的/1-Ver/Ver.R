# -----------------------------------------------------------------------------#
#                            GSE220135 数据分析
# -----------------------------------------------------------------------------#


# 1 加载R包=====================================================================
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(magrittr)
library(ggpubr)
library(ggthemes)
library(ggcorrplot)


# 2 读取数据====================================================================
# 1.1读取基因-------------------------------------------------------------------
gene <- read.csv("input/human_hub_genes.csv")
colnames(gene)[1] <- "gene"


# 1.2 读取分组------------------------------------------------------------------
group <- data.table::fread("input/GSE220135_Group.csv", data.table = F) %>% 
  dplyr::rename(Sample = ID, Group = group) %>%
  dplyr::mutate(Group = trimws(Group)) %>%            ## 去除组名两端的空格
  dplyr::arrange(desc(Group))                         ## 排序样本，疾病(BPD)在前


# 1.3 读取矩阵------------------------------------------------------------------
Matrix <- data.table::fread("input/GSE220135_Matrix_TPM.csv", data.table = F) %>% 
  tibble::column_to_rownames("V1") %>% 
  dplyr::select(group$Sample)                          ## 矩阵样本与group样本顺序一致


# 4 绘制Hub基因分组比较图=======================================================

# 提取Hub基因表达矩阵（取交集，避免基因名不匹配报错）
genes_common <- intersect(gene$gene, rownames(Matrix))
message(paste0("匹配到的Hub基因数: ", length(genes_common)))
Matrix_Hub_genes <- Matrix[genes_common, , drop = FALSE]


# 宽形数据转长格式
box_long <- Matrix_Hub_genes %>%
  t() %>%
  data.frame() %>%
  dplyr::mutate(Sample = rownames(.)) %>%
  merge(group, by.x = "Sample", by.y = "Sample") %>%
  tidyr::gather(key = geneid, value, -c(Sample, Group))

# 修改：疾病组和对照组名称！！！！！！！！！！！！！！！！！！！！！！！
box_long$Group <- factor(box_long$Group, levels = c("Control", "BPD"))#左对照，右实验疾病组

# 绘制分组箱线图
barplot <- ggplot(box_long, aes(x = geneid, y = value, fill = Group)) +
  geom_boxplot(width = 0.7, size = 0.5, outlier.color = "black", outlier.size = 1, fatten = 1,
               position = position_dodge(width = 0.8)) +
  # 修改：分组颜色
  theme_bw() + scale_fill_manual(values = c("Control" = "#8cc3dd", "BPD" = "#f4c0c6")) +
  # 去除背景网格线
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),# 去除X轴标题
        axis.text.x  = element_text(angle = 60, hjust = 1, colour = "black",
                                    size = 14),  # ← X 轴刻度字号
        axis.text.y  = element_text(size = 14),  # ← Y 轴刻度字号
        axis.title.y = element_text(size = 16)   # ← Y 轴标题字号
  ) +
  # 添加组间差异显著性标记（n=3小样本，用t.test提高统计功效） 
  stat_compare_means(method = "t.test",
                     symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                        symbols = c("***", "**", "*", "ns")),
                     label = "p.signif") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, color = "black")) +
  # 设置y轴标题
  ylab('Gene Expression') +
  # 重命名图例标题
  labs(fill = 'Group') +
  # 所有文本元素（轴标签、图例、标题）设为14号字
  theme(text = element_text(size = 14))

barplot

ggsave(barplot, file = "output/1-Boxplot_GSE220135_HubGenes.pdf", width = 9, height = 5)


# # 添加组间差异显著性标记
# stat_compare_means(method = "wilcox.test",
