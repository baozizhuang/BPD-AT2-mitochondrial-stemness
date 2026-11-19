# -----------------------------------------------------------------------------#
# 描述: 本脚本用于执行单细胞转录组分析 3_Analysis.R
# R版本：4.3.2
# Seurat (Version:5.0.1)
# 日期: 2026/06/23
# 创建：阿伦
# 代码版本：V2
# -----------------------------------------------------------------------------#



# 加载R包=======================================================================
library(clusterProfiler)
library(org.Hs.eg.db)
library(dittoSeq)
library(RColorBrewer)
library(tidyr)
library(stringr)



# 9 重新对细胞簇命名（根据最终注释结果）========================================

# 设置当前细胞身份标识
Idents(sc_integer_singleR) <- "seurat_clusters"

# 根据 celltype 数据框中 ClusterID 列的值来筛选 sc_integer_singleR对象中的细胞
sc_integer_singleR <- subset(sc_integer_singleR, idents =  celltype$ClusterID)

# 重新赋予细胞注释
sc_integer_singleR$cell_anno <- factor(sc_integer_singleR$seurat_clusters, labels = celltype2$cell)

# 设置当前细胞身份标识
Idents(sc_integer_singleR) <- "cell_anno"


# # 9.1去除各样本中细胞数量<10的细胞类型(技术与生物学背景噪音) (可选)-------------
# # 这种为技术与生物学背景噪音,在统计上没有任何意义,可以忽略的噪音
# cell_counts <- table(sc_integer_singleR$cell_anno, sc_integer_singleR$sample)
# idx <- which(cell_counts < 10, arr.ind = TRUE)
# 
# if (nrow(idx) > 0) {
#   combo <- paste(sc_integer_singleR$cell_anno, sc_integer_singleR$sample, sep = "_")
#   rare <- paste(rownames(cell_counts)[idx[,1]], colnames(cell_counts)[idx[,2]], sep = "_")
#   keep <- !combo %in% rare
#   sc_integer_singleR <- subset(sc_integer_singleR, cells = colnames(sc_integer_singleR)[keep])
#   cat(">>> 去除", sum(!keep), "个细胞（细胞数<10）\n")
# }


# 9.2可视化完全注释好后的UMAP图-------------------------------------------------
UMAP2 <- DimPlot(sc_integer_singleR,
                 group.by="cell_anno",
                 label=T,
                 cols = More_colors,
                 # alpha = 0.4, # 点透明度
                 pt.size = 0.05,
                 reduction='umap')+
  labs(x = "UMAP 1", y = "UMAP 2") +
  guides(color = guide_legend(override.aes = list(size = 4))) +                 # 调大图例圆圈
  labs(title = NULL)+
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1), # 添加边框
        axis.line = element_blank(),                                            # 去掉默认轴线
        axis.ticks = element_blank(),                                           # 去掉刻度线
        axis.text = element_blank(),                                            # 去掉刻度标签
        axis.title = element_text(size = 12, color = "black")  # 控制坐标轴标签字体大小
  )
UMAP2
ggsave(filename = "output/6-NEW-UMAP_Anno.pdf", plot = UMAP2, width = 18, height =14, units = "cm")


# # 9.2统计各细胞类型数量与比例并可视化
# # 计算各细胞类型的数量和比例
# cell_count <- table(sc_integer_singleR$cell_anno)
# cell_pct  <- round(prop.table(cell_count) * 100, 2)
# cell_stats <- data.frame(CellType = names(cell_count),
#                          Count = as.integer(cell_count),
#                          Proportion = paste0(cell_pct, "%"))
# print(cell_stats, row.names = FALSE)
# 
# # 创建带数量和比例的标签并绘制UMAP
# cell_label <- paste0(names(cell_count), " (n=", cell_count, ", ", cell_pct, "%)")
# sc_integer_singleR$cell_anno_label <- factor(sc_integer_singleR$cell_anno,
#                                              levels = names(cell_count), labels = cell_label)
# 
# UMAP2_stats <- DimPlot(sc_integer_singleR,
#                        group.by = "cell_anno_label", label = TRUE, pt.size = 0.05, reduction = "umap") +
#   labs(title = NULL) +
#   theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
#         axis.line = element_blank(), axis.ticks = element_blank(),
#         axis.text = element_blank()) +
#   NoLegend()
# ggsave("output/6-NEW-UMAP_Anno_Stats.pdf", plot = UMAP2_stats, width = 19, height = 17, units = "cm")


# # 9.3 可视化完全注释好后的TSNE图------------------------------------------------
# TSNE2 <- DimPlot(sc_integer_singleR,
#                  group.by="cell_anno",
#                  label=T,
#                 #cols = More_colors,
#                  alpha = 0.4, # 点透明度
#                  pt.size = 0.05,
#                  reduction='tsne')+
#   labs(x = "UMAP 1", y = "UMAP 2") +
#   guides(color = guide_legend(override.aes = list(size = 4))) +                 # 调大图例圆圈
#   labs(title = NULL)+
#   theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1), # 添加边框
#         axis.line = element_blank(),                                            # 去掉默认轴线
#         axis.ticks = element_blank(),                                           # 去掉刻度线
#         axis.text = element_blank(),                                            # 去掉刻度标签
#         axis.title = element_text(size = 12, color = "black")  # 控制坐标轴标签字体大小
#   )
# TSNE2
# ggsave(filename = "output/6-NEW-TSNE_Anno.pdf", plot = TSNE2, width = 19, height =13.5, units = "cm")



# 9.3正常样本的细胞聚类+患病样本的细胞聚类--------------------------------------
# UMAP版~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 按tissue分组显示
UMAP_split1 <- DimPlot(sc_integer_singleR,
                       group.by = "cell_anno",
                       split.by = "tissue",
                       label = TRUE,
                       cols = More_colors,
                       # alpha = 0.4,# 点透明度
                       pt.size = 0.05,
                       reduction = 'umap') +
  labs(x = "UMAP 1", y = "UMAP 2") +
  guides(color = guide_legend(override.aes = list(size = 4))) +                 # 调大图例圆圈
  labs(title = NULL)+
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1), # 添加边框
        axis.line = element_blank(),                                            # 去掉默认轴线
        axis.ticks = element_blank(),                                           # 去掉刻度线
        axis.text = element_blank(),                                            # 去掉刻度标签
        axis.title = element_text(size = 12, color = "black")  # 控制坐标轴标签字体大小
  )
UMAP_split1
ggsave(filename = "output/6-NEW-split-UMAP_Anno1.pdf", plot = UMAP_split1, width = 27, height =13, units = "cm")


# # 按样本个体sample分组显示
# UMAP_split2 <- DimPlot(sc_integer_singleR,
#                       group.by = "cell_anno",
#                       split.by = "sample",
#                       label = TRUE,
#                     # cols = More_colors,
#                       pt.size = 0.05,
#                       reduction = 'umap',
#                       alpha = 0.4) + # 点透明度
#   labs(x = "UMAP 1", y = "UMAP 2") +
#   guides(color = guide_legend(override.aes = list(size = 4))) +                 # 调大图例圆圈
#   labs(title = NULL)+
#   theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1), # 添加边框
#         axis.line = element_blank(),                                            # 去掉默认轴线
#         axis.ticks = element_blank(),                                           # 去掉刻度线
#         axis.text = element_blank(),                                            # 去掉刻度标签
#         axis.title = element_text(size = 12, color = "black")  # 控制坐标轴标签字体大小
#         )
# UMAP_split2
# # ggsave(filename = "output/6-NEW-split-UMAP_Anno2.pdf", plot = UMAP_split2, width = 27, height =13, units = "cm") # 短款
# ggsave(filename = "output/6-NEW-split-UMAP_Anno2.pdf", plot = UMAP_split2, width = 50, height =13, units = "cm") # 5连图



# # TSNE版~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# TSNE_split <- DimPlot(sc_integer_singleR,
#                       group.by = "cell_anno",
#                       split.by = "tissue",
#                       label = TRUE,
#                      # cols = More_colors,
#                       alpha = 0.4, # 点透明度
#                       pt.size = 0.05,
#                       reduction = 'tsne') +
#   labs(x = "UMAP 1", y = "UMAP 2") +
#   guides(color = guide_legend(override.aes = list(size = 4))) +                 # 调大图例圆圈
#   labs(title = NULL)+
#   theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1), # 添加边框
#         axis.line = element_blank(),                                            # 去掉默认轴线
#         axis.ticks = element_blank(),                                           # 去掉刻度线
#         axis.text = element_blank(),                                            # 去掉刻度标签
#         axis.title = element_text(size = 12, color = "black")  # 控制坐标轴标签字体大小
#   )
# TSNE_split
# ggsave(filename = "output/6-NEW-split-TSNE_Anno.pdf", plot = TSNE_split, width = 27, height =14, units = "cm")



# 10 绘制诊断标志物marker的气泡图与热图与小提琴图===============================
# 10.1找差异表达基因------------------------------------------------------------
# 对已经注释好的细胞类型进行差异表达基因（DEGs）的验证和筛选
# 检查通过FindAllMarkers找到的每个细胞类型的差异基因,最终为每个细胞类型筛选出3个最具代表性的基因

# FindAllMarkers寻找每个细胞簇的差异表达基因
dat_degs <- FindAllMarkers(sc_integer_singleR,
                           min.pct = 0.2,  # 可以是0.3
                           logfc.threshold = 0.25)

# 要注意啊，我这里的差异基因是筛选的满足|logFC| >0.25的，
# 后续的AUcell前的差异分析火山图筛选要|logFC| >2的，后续的差异基因得按照
# 阈值|logFC| >2再筛选，但此处的不需要，容易混淆

# 用与后续的AUcell之前的差异分析，后续需要再筛选
saveRDS(dat_degs,"output2/dat_degs.rds")
dat_degs <- readRDS("output2/dat_degs.rds")

# 将dat_degs转换成一个数据框，赋值给DEG
DEG <- dat_degs %>%as.data.frame() 



# # 10.2 top3的诊断标志物marker数据筛选-------------------------------------------
# 用于给客户展示诊断标志物
# 从差异表达分析结果中，为每个细胞类型提取前3个表达量最显著上调的基因
# 输出到gene3

# ----------------新增：定义非编码/低质量基因过滤函数--------------------------#

# # 人类数据专用：排除常见非编码/预测/未知功能基因~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# is_likely_protein_coding <- function(gene) {
#   g <- toupper(gene)  # 人类基因通常用大写表示，统一转为大写
#   
#   # 排除线粒体基因（人类通常是MT-前缀）
#   is_mt <- grepl("^MT-", g)
#   
#   # 排除核糖体蛋白（人类也是 RPL*, RPS*）
#   is_rp <- grepl("^RP[L|S][0-9]", g)
#   
#   # 排除假基因（通常以P结尾或有P前缀）
#   is_pseudogene <- grepl("P[0-9]+$", g) | grepl("-P[0-9]+$", g) | grepl("^P[0-9]+-", g)
#   
#   # 排除LOC基因（未定位基因）
#   is_loc <- grepl("^LOC[0-9]", g)
#   
#   # 排除非编码RNA（常见的人类非编码RNA模式）
#   is_ncrna <- grepl("^(MIR|SNORD|SNHG|RNU|NEAT|MALAT|XIST|HOTAIR|DANCR|TUG1|H19|TERC|TERRA)", g)
#   
#   # 排除预测基因（常见模式）
#   is_predicted <- grepl("^OR[0-9]", g) |  # 嗅觉受体
#     grepl("^C[0-9]+ORF", g) |  # 开放阅读框
#     grepl("^LINC[0-9]", g) |  # 长链非编码RNA
#     grepl("^AC[0-9]", g) |  # 序列片段
#     grepl("^AL[0-9]", g)  # 序列片段
#   
#   # 排除免疫相关假基因（人类常见）
#   is_ig_pseudo <- grepl("^IG[H|K|L][V|D|J|C]P", g) | grepl("^TR[G|D|B|A][V|D|J|C]P", g)
#   
#   # 合并所有排除条件
#   exclude <- is_mt | is_rp | is_pseudogene | is_loc | is_ncrna | is_predicted | is_ig_pseudo
#   
#   # 返回：不是排除项的基因
#   !exclude
# }

# 小鼠数据专用：排除常见非编码/预测/未知功能基因~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
is_likely_protein_coding <- function(gene) {
  g <- gene  # 小鼠基因通常保持原始大小写（首字母大写）

  # 排除线粒体基因（小鼠是mt-前缀）
  is_mt <- grepl("^mt-", g)

  # 排除核糖体蛋白（小鼠 Rpl*, Rps*）
  is_rp <- grepl("^Rp[l|s][0-9]", g)

  # 排除假基因
  is_pseudogene <- grepl("ps$", g) |  # 如 Rpl23ps
    grepl("-ps[0-9]*$", g) |
    grepl("^Gm[0-9]+-ps", g)  # Gm假基因

  # 排除LOC基因（未定位基因）
  is_loc <- grepl("^LOC[0-9]", g)

  # 排除Gm基因（Gene model，多为预测或未知功能）
  # 注意：少数Gm基因有功能，可根据需要注释掉这行
  is_gm <- grepl("^Gm[0-9]", g)

  # 排除RIKEN cDNA（典型格式：数字+字母+Rik）
  # 如 9230105E05Rik, 4930519G02Rik 等
  is_riken <- grepl("Rik$", g)

  # 排除非编码RNA
  is_ncrna <- grepl("^(Mir|Snord|Snhg|Rnu|Neat|Malat|Xist|Hotair|Dancr|Tug1|H19|Terc|Gm[0-9]+$)", g) |
    grepl("^Linc[0-9]", g) |
    grepl("-as$", g)  # 反义链RNA

  # 排除预测基因（小鼠特有模式）
  is_predicted <- grepl("^AU[0-9]", g) |  # AU开头的序列标签
    grepl("^AV[0-9]", g) |  # AV开头的序列标签
    grepl("^BB[0-9]", g) |  # BB开头的序列标签
    grepl("^AK[0-9]", g) |  # AK开头的cDNA
    grepl("^BC[0-9]", g)    # BC开头的cDNA

  # 排除免疫相关假基因
  is_ig_pseudo <- grepl("^Igh[V|D|J|C]p", g) | grepl("^Tr[G|D|B|A][V|D|J|C]p", g)

  # 合并所有排除条件
  exclude <- is_mt | is_rp | is_pseudogene | is_loc | is_gm |
    is_riken | is_ncrna | is_predicted | is_ig_pseudo

  !exclude
}
# -----------------------------------------------------------------------------#

DEG_celltype <- unique(DEG$cluster)
testseu <- sc_integer_singleR
gene3 <- data.frame()

for(ct in DEG_celltype){
  DEG_gene <- DEG %>% 
    filter(cluster == ct) %>%
    filter(p_val_adj < 0.05) %>%          # 显著性筛选
    mutate(pct = pct.1 - pct.2) %>%       # 计算比例差
    arrange(desc(avg_log2FC), desc(pct))  # 双重排序
  
  # 非编码基因过滤
  DEG_gene <- DEG_gene %>%
    filter(is_likely_protein_coding(gene))
  
  # 取最多前3个
  n_take <- min(3, nrow(DEG_gene))
  if (n_take > 0) {
    gene3 <- rbind(gene3, DEG_gene[1:n_take, c("cluster", "gene")])
  }
}


# 10.3 top3的分组气泡图与小提琴图可视化-----------------------------------------
# 10.3.1 气泡图1 (自动top3的 图片类型1 )~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 必要的输入：
# seurat对象，meta数据框中有表示细胞类型的列和表示分组的列；
# 提供细胞类型的顺序；
# 表示marker基因的数据框，至少有两列，一列为cluster，一列为gene，基因名无重复。
marker <- unique(gene3$gene)
pdf("output2/7-cellmarker_bubble_1.pdf", height = 5, width = 17)
dotplot <- DotPlot(testseu,features = marker,group.by = ,dot.scale = 8)+
  ggpubr::rotate_x_text(90)+
  scale_color_gradientn(colours = c("grey","red")) #颜色
dotplot
dev.off()


# 10.3.2 气泡图2 (自动top3的 图片类型2 )~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
marker <- unique(gene3$gene)

pdf("output2/7-cellmarker_bubble_2.pdf", height = 4.5, width = 17)
dotplot_enhanced <- DotPlot(testseu, features = marker, group.by = "cell_anno") +
  scale_color_gradientn(
    colors = rev(brewer.pal(7, "RdBu")),
    guide = guide_colorbar(title = "Mean expression", barwidth = 0.8, barheight = 4)
  ) +
  scale_size_continuous(
    range = c(0, 8), # 👈 原来是 c(0,5)，现在改成 c(0,8) 让气泡更大！
    breaks = c(25, 50, 75, 100), # 手动指定显示的数字，跳过 0
    labels = c("25%", "50%", "75%", "100%")
  ) +
  guides(size = guide_legend(
    title = "Percent Expressed",
    override.aes = list(shape = 21, colour = "black", fill = NA, stroke = 0.5)
  )) +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    strip.text = element_blank(),
    panel.spacing = unit(0.1, "lines"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid.major = element_line(color = "white", linewidth = 0.2),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, color = "black"),
    axis.text.y = element_text(size = 13, color = "black"),
    axis.title = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  )

print(dotplot_enhanced)
dev.off()


# FeaturePlot(sc_integer_singleR,
#             features = "S100A9",
#             pt.size = 0.05,
#             label.size = 3
# )

# 先看看
print(gene3)


# 10.3.3 手动赋值版 气泡图1~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 可根据gene3修改
# 绘制报告与文章要用的marker图,避免掉非编码基因
# 完整细胞类型标记基因词典(气泡图2也沿用这个列表的)
Display_genes_list <- list(
  "B cell" = c("Cd79a", "Ms4a1", "Pax5"),
  "gCAP" = c("Gpihbp1", "Aplnr", "Kit"),
  "Fibroblast" = c("Col1a1", "Dcn", "Lum"),
  "SMC" = c("Tagln", "Cnn1","Myh11"),
  "aCAP" = c("Ednrb", "Apln", "Car4"),
  "T cell" = c("Cd3d", "Cd3e"),
  "AT2" = c("Sftpc", "Sftpa1", "Sftpb"),
  "Macrophage" = c("Cd68","Csf1r","Cd163"),
  "Club cell" = c("Scgb1a1", "Scgb3a2"), 
  "Neutrophil" = c("S100a8", "S100a9", "Csf3r"), 
  "NK cell " = c("Ncr1","Klrb1c"),
  "Lymph" = c("Prox1", "Lyve1", "Flt4"),
  "Pericyte" = c("Rgs5", "Pdgfrb"),
  "AT1" = c("Ager", "Hopx", "Pdpn"),
  "Ciliated cell" = c("Foxj1", "Tuba1a", "Cfap53"),
  "Platelet" = c("Treml1", "Ppbp", "Gp1ba")
)

# 设置细胞簇标识
Idents(sc_integer_singleR) <- "cell_anno"


# 气泡图可视化~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Display_dot_plot <- DotPlot(sc_integer_singleR,
                            features = Display_genes_list,
                            group.by = "cell_anno", # 确保按细胞类型分组
                            dot.scale = 8) +  # 增大气泡大小
  # 🔴🔵 关键：用 scale_color_gradientn 控制红蓝配色
  scale_color_gradientn(
    colors = rev(brewer.pal(n = 7, name = "RdBu")),  # 红-白-蓝（高表达为红）
    name = "Mean\nexpression"
  ) +
  scale_size_continuous(
    range = c(0, 8), # 增大气泡尺寸
    breaks = c(25, 50, 75, 100), # 手动指定显示的百分比点
    labels = c("25%", "50%", "75%", "100%"),
    name = "Percent Expressed"
  ) +
  guides(size = guide_legend(
    title = "Percent Expressed",
    override.aes = list(shape = 21, colour = "black", fill = NA, stroke = 0.5)
  )) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey98", color = "black", linewidth = 0.3), # 显示strip背景
    strip.text = element_text(face = "bold", size = 11, color = "black"), # 显示strip文本
    panel.spacing.x = unit(0.3, "lines"), # 增加水平间距
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid.major = element_line(color = "grey98", linewidth = 0.2),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, color = "black"),
    axis.text.y = element_text(size = 13, color = "black"),
    axis.title = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  labs(x = "Markers", y = "Cell Clusters", title = "") + # 标题留空
  theme(plot.title = element_text(hjust = 0.5, size = 13, face = "bold"))

pdf("output2/7-cellmarker_bubble_3.pdf", height = 5, width = 19) # 大尺寸
# pdf("output2/7-cellmarker_bubble_3.pdf", height = 3.8, width = 14) # 小尺寸

print(Display_dot_plot)
dev.off()


# 10.3.4 手动赋值版 气泡图2~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 加载scCustomize包
library(scCustomize)

# 绘图颜色设置
my36colors <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488",
                "#F39B7F", "#8491B4", "#91D1C2", "#DC0000",
                "#7E6148", "#B09C85", "#E18727", "#20854E",
                "#7876B1", "#6F99AD", "#FFDC91", "#EE4C97",
                "#CD534C", "#7AA6DC", "#8C6BB1", "#F0A3FF",
                "#A0D2DB", "#1C7C54", "#E88873", "#F3D266",
                "#B39EB5", "#59A14F")

# 聚类点图（使用手动赋值版marker列表）
pdf("output2/7-cellmarker_bubble_4.pdf", width = 13,height = 3.6,pointsize = 8,onefile = FALSE) # 大尺寸
# pdf("output2/7-cellmarker_bubble_4.pdf", width = 11,height = 3.3,pointsize = 7,onefile = FALSE) # 小尺寸

bubble_plot <- Clustered_DotPlot(seurat_object = sc_integer_singleR,
                                 colors_use_exp = rev(brewer.pal(n = 7, name = "RdBu")), # 蓝-白-红（高表达为红），与气泡图1统一
                                 colors_use_idents = my36colors,
                                 features = unname(unlist(Display_genes_list)), # 展开为向量
                                 plot_km_elbow = FALSE,
                                 cluster_feature = FALSE,  # 关闭基因聚类
                                 cluster_ident = FALSE,    # 关闭细胞类型聚类
                                 flip = TRUE, # 翻转坐标轴
                                 x_lab_rotate = 60)  # 横坐标角度为60
print(bubble_plot)
dev.off()



# 12 展示细胞组成比例===========================================================
# 12.1 计算并可视化样本间的细胞比例---------------------------------------------
# 计算每个样本的细胞比例
# 计算各组样本不同细胞群比例: 1 = row, 2 = column, default is NULL
Cellratio_sample <- prop.table(table(Idents(sc_integer_singleR), gsub("/","",sc_integer_singleR$sample)), margin = 2)
Cellratio_sample <- as.data.frame(Cellratio_sample)
colnames(Cellratio_sample)[1] <- "Cell"

# (选择1)在这里添加：强制 Var2 按原始样本出现顺序排列
# original_sample_order <- rev(unique(gsub("/", "", sc_integer_singleR$GEOID))) 

# (选择2) 按字母顺序排列样本
original_sample_order <- str_sort(unique(gsub("/", "", sc_integer_singleR$sample)), 
                                  numeric = TRUE, 
                                  decreasing = TRUE)


Cellratio_sample$Var2 <- factor(Cellratio_sample$Var2, levels = original_sample_order)


# 绘制样本细胞比例堆叠柱状图
ratioplot_sample <- ggplot(Cellratio_sample) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Cell),
           stat = "identity", width = 0.7, size = 0.4, colour = NA)+ 
  theme_classic() +
  labs(x='Sample', y = 'Ratio')+
  scale_fill_manual(values= colorRampPalette(c('#c76753', '#cd8a9d', '#dac386', '#7f7f7f', '#658ba5', '#9b86b6'))(NROW(table(Cellratio_sample$Cell))))+
  # scale_fill_futurama()+
  # 确保轴标签和文本为黑色
  theme(
    # 添加边框：fill=NA 表示边框内部透明（不遮挡背景），color 为黑色，size 控制粗细
    panel.border = element_rect(fill = NA, color = "black", size = 0.6, linetype = "solid"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.text = element_text(color = "black"),
    legend.title = element_text(color = "black")
  ) +
  coord_flip()
#theme(panel.border = element_rect(fill=NA, color="black", size=0.5, linetype="solid"),
#legend.position = "bottom")
ratioplot_sample
ggsave("output2/8-Cellratio_Sample.pdf", plot = ratioplot_sample, scale = 1, width = 14, height = 9, units =c("cm"))


# 12.2 计算并可视化组间的细胞比例-----------------------------------------------
#计算各组样本不同细胞群比例: 1 = row, 2 = column, default is NULL
Cellratio_group <- prop.table(table(Idents(sc_integer_singleR), sc_integer_singleR$tissue), margin = 2)
Cellratio_group

Cellratio_group <- as.data.frame(Cellratio_group)
colnames(Cellratio_group)[1] <- "Cell"

ratioplot_group <- ggplot(Cellratio_group) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Cell),
           stat = "identity", width = 0.7, size = 0.4, colour = NA)+ 
  theme_classic() +
  labs(x='Group', y = 'Ratio')+
  scale_fill_manual(values= colorRampPalette(c('#c76753', '#cd8a9d', '#dac386', '#7f7f7f', '#658ba5', '#9b86b6'))(NROW(table(Cellratio_group$Cell))))+
  #scale_fill_futurama()+
  # 确保轴标签和文本为黑色
  theme(
    # 添加边框：fill=NA 表示边框内部透明（不遮挡背景），color 为黑色，size 控制粗细
    panel.border = element_rect(fill = NA, color = "black", size = 0.6, linetype = "solid"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.text = element_text(color = "black"),
    legend.title = element_text(color = "black")
  ) +
  coord_flip()
#theme(panel.border = element_rect(fill=NA, color="black", size=0.5, linetype="solid"),
#legend.position = "bottom")
ratioplot_group
ggsave("output2/9-Cellratio_Group.pdf", plot = ratioplot_group, scale = 1, width = 16, height = 9, units =c("cm"))



# # 13 可视化模型关键基因表达点图和表达特征图=====================================
# # 13.1点图----------------------------------------------------------------------
# # 13.1.1不分面~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# dotplot <- DotPlot(sc_integer_singleR,features = model_marker,dot.scale = 5)+
#   ggpubr::rotate_x_text(90)+
#   scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33'))+ #颜色
#   theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8))# 加边框
# dotplot
# ggsave("output2/10-Dotplot.pdf", plot = dotplot, width = 17, height = 12, units = "cm")
# 
# 
# # 13.1.2分面(z score)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # 检查基因是否在数据中
# key_genes <- model_marker[model_marker %in% rownames(sc_integer_singleR)]
# 
# # 提取metadata和表达数据
# plot_data <- cbind(sc_integer_singleR@meta.data, FetchData(sc_integer_singleR, vars = key_genes))
# 
# # 计算每个细胞类型中基因的平均表达
# plot_data_long <- plot_data %>%
#   group_by(cell_anno, tissue) %>%
#   summarise(across(all_of(key_genes), mean), .groups = "drop") %>%
#   pivot_longer(cols = all_of(key_genes),
#                names_to = "Gene",
#                values_to = "Expression")
# 
# # 按基因进行z-score标准化（跨细胞类型），使表达模式可比
# plot_data_long <- plot_data_long %>%
#   group_by(Gene) %>%
#   mutate(Expression = scale(Expression, center = TRUE, scale = TRUE)[, 1]) %>%
#   ungroup()
# 
# # 计算细胞类型在每组中的比例
# cell_prop <- plot_data %>%
#   group_by(tissue, cell_anno) %>%
#   summarise(Count = n(), .groups = "drop") %>%
#   group_by(tissue) %>%
#   mutate(Proportion = Count / sum(Count))
# 
# # 合并表达数据和比例数据
# plot_data_final <- left_join(plot_data_long, cell_prop, by = c("cell_anno", "tissue"))
# 
# # 绘制气泡图
# dotplot_facet <- ggplot(plot_data_final, aes(x = cell_anno, y = Gene)) +
#   geom_point(aes(size = Proportion, color = Expression)) +
#   scale_color_gradientn(colours = rev(brewer.pal(n = 7, name = "RdBu"))) +
#   facet_grid(. ~ tissue) +
#   theme_classic() +
#   theme(
#     axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
#     axis.line = element_line(linewidth = 0.5),
#     axis.ticks = element_line(linewidth = 0.5),
#     panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8) # 加边框
#   ) +
#   labs(x = "Cell Type", y = "Hub Gene")
# 
# dotplot_facet
# ggsave("output2/10-Dotplot_facet.pdf", plot = dotplot_facet, width = 20, height = 10, units = "cm")
# 
# 
# # 13.2 特征图-------------------------------------------------------------------
# # 13.2.1 uamp~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# featureplot <- FeaturePlot(sc_integer_singleR,
#                            features = model_marker,
#                            pt.size = 0.1,
#                            label.size = 3,
#                            ncol = 3) + #一行显示3个
#     theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1), # 添加边框
#           axis.line = element_blank(),                                            # 去掉默认轴线
#           axis.ticks = element_blank(),                                           # 去掉刻度线
#           axis.text = element_blank()                                             # 去掉刻度标签
# )
# featureplot
# ggsave("output2/11-Featureplot.pdf", plot = featureplot, width = 32, height = 18, units = "cm")
# 
# 
# # 单个基因展示
# featureplot <- FeaturePlot(sc_integer_singleR,
#                            features = "CHAF1A",
#                            pt.size = 0.5,
#                            label.size = 3,
#                            ncol = 1)+  #一行显示1个
#   theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1), # 添加边框
#         axis.line = element_blank(),                                            # 去掉默认轴线
#         axis.ticks = element_blank(),                                           # 去掉刻度线
#         axis.text = element_blank()                                             # 去掉刻度标签
#   )
# featureplot
# ggsave("output2/11-CHAF1A-Featureplot.pdf", plot = featureplot, width = 15, height = 14, units = "cm")
# 
# 
# # 13.2.2 tsne~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# featureplot <- FeaturePlot(sc_integer_singleR,
#                            reduction = "tsne",
#                            features = model_marker,
#                            pt.size = 0.05,
#                            label.size = 3,
#                            ncol = 3)+ #一行显示3个
#   theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1), # 添加边框
#         axis.line = element_blank(),                                            # 去掉默认轴线
#         axis.ticks = element_blank(),                                           # 去掉刻度线
#         axis.text = element_blank()                                             # 去掉刻度标签
#   )
# featureplot
# ggsave("output2/11-Featureplot.pdf", plot = featureplot, width = 23, height = 10, units = "cm")



# 保存最终的seurat对象,注意啊,我保存为Final_sc_integer_singleR
saveRDS(sc_integer_singleR,file = "save/Final_sc_integer_singleR.rds")


# 读取为Final_sc_integer_singleR,用于后续分析
# sc_integer_singleR <- readRDS("save/Final_sc_integer_singleR.rds")



# # ---------------------------可选步骤（降采样）--------------------------------#
# 
# sc_integer_singleR <- SketchData(
#          object = sc_integer_singleR,
#          assay = "RNA",           # 使用RNA assay（或SCT等）
#          ncells = 20000,           # 采样20000个细胞
#          method = 'Uniform',      # 均匀采样方法
#          sketched.assay = 'sketch_uniform',  # 创建新的sketch assay
#          seed = 123               # 设置随机种子确保可重现
#   )
# 
# # 提取降采样结果并赋值到新变量
# downsample_cells <- Cells(sc_integer_singleR[["sketch_uniform"]])
# sc_integer_singleR_sampled <- subset(sc_integer_singleR, cells = downsample_cells)
# 
# # 使用均匀采样后需要重新设置默认assay为RNA
# DefaultAssay(sc_integer_singleR_sampled) <- "RNA"
# 
# # 保存降采样后的数据
# saveRDS(sc_integer_singleR_sampled,file = "save/sc_integer_singleR_sampled.rds")
# 
# # ---------------------------可选步骤（降采样）--------------------------------#


save.image(file = "save/SC.RData.gz", compress = "gzip")

# 读取整个工作空间
load("save/SC.RData")
