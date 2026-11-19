# -----------------------------------------------------------------------------#
# 描述: 本脚本用于执行单细胞转录组分析 2_Cell_Annotation.R
# R版本：4.3.2
# Seurat (Version:5.0.1)
# 日期: 2026/06/23
# 创建：阿伦
# 代码版本：V2
# -----------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
# <<<<<<<<<<<<<< 本注释流程采用SingleR细胞注释定义seurat 主文件>>>>>>>>>>>>>>>>>
# <<<<<<<<<<<<<<<<<<< mLLMCelltype与DeepCellSeek只是辅助 >>>>>>>>>>>>>>>>>>>>>>>
#------------------------------------------------------------------------------#


# 添加镜像======================================================================
# 添加bioconductor清华镜像
options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
# 添加bioconductor西湖大学镜像
options(BioC_mirror = "https://mirrors.westlake.edu.cn/bioconductor")
# 添加R清华镜像
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))


# 配置动态R包库,我服务器安装R包没权限，得设置动态R包库(服务器）=================
# 我的个人R包库在/work/run/projects/bio-35/R_library
.libPaths(c("/work/run/projects/bio-35/R_library", .libPaths()))



# 下载R包=======================================================================
# install.packages("mLLMCelltype")#mLLMCelltype包
# install.packages("SCpubr")
# install.packages("showtext")
# install.packages("assertthat")
# install.packages("stringr")


# 加载R包=======================================================================
library(mLLMCelltype)
library(SCpubr)
library(viridis)  # 用于调色板
library(showtext) # 帮助加载中文字体，从而正确显示中文字符
library(assertthat)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichR)

showtext_auto()#自动加载系统中的中文字体，从而在R图形中正确显示中文字符



# 6 SingleR 细胞注释============================================================

# 6.1.加载参考数据集------------------------------------------------------------
# 加载参考数据集(人类)
# hpca <- celldex::HumanPrimaryCellAtlasData()#人类原代细胞的转录组数据
# saveRDS(hpca,file = "input/hpca.rds")
# hpca <- readRDS(file = "input/hpca.rds")

# 加载参考数据集(小鼠)
# mrsd <- celldex::MouseRNAseqData()
# saveRDS(mrsd,file = "input/mrsd.rds")
mrsd <- readRDS(file = "input/mrsd.rds")

# 其他数据集可供选择
# hpca <- celldex::HumanPrimaryCellAtlasData() #人类原代细胞的转录组数据
# bep <- celldex::BlueprintEncodeData()
# igd <- celldex::ImmGenData()
# ice <- celldex::DatabaseImmuneCellExpressionData()


# 获取标准化后矩阵
testdata <- GetAssayData(sc_integer, layer = "data") 


# 物种间的基因名转换(用不到)
# library(homologene)
# genelist <- homologene(testdata@Dimnames[[1]],inTax = 10116,outTax = 10090)
# genelist <- genelist[!duplicated(genelist$`10090`),]
# testdata <- testdata[genelist$`10116`,]
# dim(testdata)
# rownames(testdata) <- genelist$`10090`


# 6.2运行SingleR注释------------------------------------------------------------
# 运行SingleR进行细胞类型预测(小鼠)
cellpred <- SingleR::SingleR(test = testdata,
                             ref = mrsd,              # 小鼠参考注释数据集
                             labels = mrsd$label.main,# 参考集的细胞类型标签
                             de.method = "classic",   # 指定用于寻找标记基因的差异表达分析方法，默认也为classic        # 修改
                             clusters = sc_integer@meta.data$seurat_clusters) # 包含了每个细胞所属的簇（cluster）的编号

# # 运行SingleR进行细胞类型预测(人)
# cellpred <- SingleR::SingleR(test = testdata,
#                              ref = hpca,                # 人类参考注释数据集
#                              labels = hpca$label.main,  # 参考集的细胞类型标签
#                              de.method = "classic",     # 指定用于寻找标记基因的差异表达分析方法，默认也为classic        # 修改
#                              clusters = sc_integer@meta.data$seurat_clusters)# 包含了每个细胞所属的簇（cluster）的编号

# 整理预测结果
celltype = data.frame(ClusterID=rownames(cellpred),
                      cell=cellpred$labels,
                      stringsAsFactors = F)

# 得分热图查看
# print(plotScoreHeatmap(cellpred))


# 6.3 查看这些细胞类型对应的Cluster（判断后续是否要人工注释）-------------------
# 提取singleR注释标签
# labels 是 SingleR 方法给每个细胞群或每个细胞分配的标签，表示该细胞的分类结果。这个标签是基于参考数据（如 MouseRNAseqData）进行注释的。例如，“T_cells”、“B_cell”等代表细胞类型。
# pruned.labels 是经过修剪的标签结果。当 delta.next 的值过小时，意味着该细胞的分类不够可靠，SingleR 会将这些分类结果设置为 NA，即未分类，以提高分类的准确性。因此，pruned.labels 是去除低置信度分类后的标签，保留了分类结果较为可靠的细胞标签。
# 但我用labels
new.cluster.ids <- cellpred$pruned.labels
# 确保每个细胞的聚类标签与对应的细胞标识符（idents）对齐
names(new.cluster.ids) <- levels(sc_integer)

# 构建新对象进行重命名聚类标识符
sc_integer_singleR <- RenameIdents(sc_integer,new.cluster.ids)

# 将 SingleR 注释保存到元数据列 cell_anno 中
sc_integer_singleR$cell_anno <- Idents(sc_integer_singleR)


# 绘制新对象的umap图
singleR_UMAP <- DimPlot(sc_integer_singleR,
                        reduction = "umap",
                        label = T,
                        pt.size = 0.5)+
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1), # 添加边框
        axis.line = element_blank(),                                            # 去掉默认轴线
        axis.ticks = element_blank(),                                           # 去掉刻度线
        axis.text = element_blank()                                             # 去掉刻度标签
  )
singleR_UMAP
# ggsave(filename = "output/5-singleR_UMAP.pdf", plot = singleR_UMAP, width = 16, height = 12, units = "cm")
ggsave(filename = "output/5-singleR_UMAP.pdf", plot = singleR_UMAP, width = 21, height = 16, units = "cm")


# # 绘制新对象的tsne图
# singleR_TSNE <- DimPlot(sc_integer_singleR, reduction = "tsne",label = T,pt.size = 0.5)
# singleR_TSNE
# ggsave(filename = "output/5-singleR_TSNE.pdf", plot = singleR_TSNE, width = 17, height = 13, units = "cm")



#------------------------保存sc_integer_singleR的rds----------------------------#
# saveRDS(sc_integer_singleR, "save/sc_integer_singleR.rds")
# sc_integer_singleR <- readRDS("save/sc_integer_singleR.rds")
#------------------------保存sc_integer_singleR的rds----------------------------#



# ### 7.mLLMCelltype细胞注释 =====================================================
# # github教程：https://github.com/cafferychen777/mLLMCelltype/blob/main/README_CN.md
# 
# # 7.1 支持模型------------------------------------------------------------------
# #OpenAI: GPT-4.1/GPT-4.5/GPT-4o
# #Anthropic: Claude-3.7-Sonnet/Claude-3.5-Haiku
# #Google: Gemini-2.0-Pro/Gemini-2.0-Flash
# #Alibaba: Qwen2.5-Max Qwen3-Max
# #DeepSeek: DeepSeek-V3/DeepSeek-R1
# #Minimax: MiniMax-Text-01
# #Stepfun: Step-2-16K
# #Zhipu: GLM-4
# #X.AI: Grok-3/Grok-3-mini
# #OpenRouter: 通过单一API访问多种模型
# #支持来自OpenAI、Anthropic、Meta、Google、Mistral等多家提供商的模型
# #格式: 'provider/model-name'（例如：'openai/gpt-4o'、'anthropic/claude-3-opus'）
# #提供免费模型，使用:free后缀（例如：'microsoft/mai-ds-r1:free'、'deepseek/deepseek-chat:free'）
# 
# #------------------------------------------------------------------------------#
# # 对于中国大陆用户，强烈推荐直接使用网页版：https://mllmcelltype.com!!!!!!
# #------------------------------------------------------------------------------#
# 
# 
# # 7.2 使用国内模型进行多种大模型共识注释----------------------------------------
# # 怎么调都只能注释的很细分，无法进行大类注释，但可用于参考
# # 还可能会出现同一种细胞不同命名写法的尴尬情况，得仔细查看，如同时存在B cell 和B cells
# 
# # 从CSV配置文件读取API密钥
# Credential-related configuration omitted from publication copy; configure locally.
# Credential-related configuration omitted from publication copy; configure locally.
# 
# # 创建专属于用于mLLMCelltype注释的seurat对象
# sc_integer_mLLMCelltype <- copy(sc_integer)
# 
# # 为每个聚类寻找标记基因(使用FindAllMarkers)
# sc_integer_markers_mLLMCelltype <- FindAllMarkers(sc_integer_mLLMCelltype, 
#                                      only.pos = TRUE,           # 仅保留上调基因
#                                      min.pct = 0.2,             # 基因在至少20%的细胞中表达
#                                      logfc.threshold = 0.25)    # 表达量差异≥0.25倍
# 
# # 设置缓存目录以加快处理速度
# cache_dir <- "./mllmcelltype_cache"
# dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
# 
# 
# # 支持的模型(请查看作者的github网页)
# 
# 
# # 使用多个LLM模型运行LLMCelltype注释
# # 讨论round数不确定，整体运行时间很长10-30分钟，需耐心等待
# consensus_results <- interactive_consensus_annotation(
#   input = sc_integer_markers_mLLMCelltype,
#   tissue_name = "人类的骨肉瘤（osteosarcoma）组织",#组织的描述或者名称或组织的上下文 promart
#   models = c(
#     "deepseek-chat",        # DeepSeek
#     "qwen-max-2025-01-25"   # 阿里云千问
#     #"glm-4-plus"           # 智谱 GLM
#   ),
# Credential-related configuration omitted from publication copy; configure locally.
#   top_gene_count = 10,         #每个细胞簇中用于注释的标记基因数量，最显著的前10个差异表达基因
#   controversy_threshold = 0.7, # 值:(0-1) 设置争议检测的阈值衡量不同模型对同一细胞簇注释的一致性程度，用于识别模型间分歧大的注释，值越低表示要求越严格
#   entropy_threshold = 1.0,     #值:(0-2) 设置熵阈值，用于检测单个模型输出的不确定性,1为中等程度的不确定性，熵值高表示模型"犹豫不决"，给出了多个可能选项，熵值低表示模型很确信自己的单一答案
#   cache_dir = cache_dir#指定缓存目录路径
# )
# 
# 
# # 保存结果
# saveRDS(consensus_results, "save/mLLMCelltype_results.rds")
# consensus_results <- readRDS("save/mLLMCelltype_results.rds")
# 
# 
# # 7.3将注释添加到Seurat对象-----------------------------------------------------
# # 从 consensus_results$final_annotations 获取细胞类型注释
# cluster_to_celltype_map <- consensus_results$final_annotations
# 
# # 创建新的细胞类型标识符列
# cell_types <- as.character(Idents(sc_integer_mLLMCelltype))
# for (cluster_id in names(cluster_to_celltype_map)) {
#   cell_types[cell_types == cluster_id] <- cluster_to_celltype_map[[cluster_id]]
# }
# 
# # 将细胞类型注释添加到Seurat对象
# sc_integer_mLLMCelltype$cell_type <- cell_types
# 
# # 提取包含指标的详细共识结果
# consensus_details <- consensus_results$initial_results$consensus_results
# 
# # 为每个聚类创建包含指标的数据框
# # 可以用于进一步分析聚类结果的稳定性和不确定性
# uncertainty_metrics <- data.frame(
#   cluster_id = names(consensus_details),
#   consensus_proportion = sapply(consensus_details, function(res) res$consensus_proportion),
#   entropy = sapply(consensus_details, function(res) res$entropy)
# )
# 
# # 获取每个细胞的聚类ID号
# current_clusters <- sc_integer_mLLMCelltype$seurat_clusters
# 
# # 将每个细胞的聚类ID与uncertainty_metrics中的相应指标匹配，并添加到元数据
# sc_integer_mLLMCelltype$consensus_proportion <- uncertainty_metrics$consensus_proportion[match(current_clusters, uncertainty_metrics$cluster_id)]
# sc_integer_mLLMCelltype$entropy <- uncertainty_metrics$entropy[match(current_clusters, uncertainty_metrics$cluster_id)]
# 
# 
# # 7.4 可视化--------------------------------------------------------------------
# # 7.4.1 UMAP 可视化
# mLLMCelltype_UMAP <- DimPlot(sc_integer_mLLMCelltype, 
#                              reduction = "umap", 
#                              group.by = "cell_type", 
#                              label = TRUE, 
#                              pt.size = 0.05) +
#   theme(plot.title = element_blank())#不显示标题
# mLLMCelltype_UMAP
# ggsave(filename = "output/5-mLLMCelltype-UMAP.pdf", plot = mLLMCelltype_UMAP, width = 17, height =12, units = "cm")
# 
# # TSNE可视化
# mLLMCelltype_TSNE <- DimPlot(sc_integer_mLLMCelltype,
#                              reduction = "tsne",
#                              group.by = "cell_type",
#                              label = TRUE,
#                              pt.size = 0.05) +
#   theme(plot.title = element_blank())#不显示标题
# mLLMCelltype_TSNE
# ggsave(filename = "output/5-mLLMCelltype-TSNE.pdf", plot = mLLMCelltype_TSNE, width = 25, height =18, units = "cm")



# 8 DeepCellSeek 细胞注释=======================================================
# 输出每个细胞簇的TOP20差异基因用于DeepCellSeek

# https://deepcellseek.kiz.ac.cn/

# 为每个聚类寻找标记基因
sc_integer_markers <- FindAllMarkers(sc_integer, 
                                     only.pos = TRUE,           # 仅保留上调基因
                                     min.pct = 0.3,             # 基因在至少30%的细胞中表达
                                     logfc.threshold = 0.25)    # 表达量差异≥0.25倍

# 筛选显著基因
markers_for_DeepCellSeek <- sc_integer_markers %>%
  dplyr::filter(p_val_adj < 0.05) %>%        # 显著
  dplyr::arrange(cluster, desc(avg_log2FC))  # 按cluster分组，按表达差异降序排

# 提取每个 cluster 的 Top 20 基因
top_genes_for_DeepCellSeek <- markers_for_DeepCellSeek %>%
  group_by(cluster) %>%
  slice_head(n = 20) %>%    
  ungroup()

# 生成目标格式：clusterX: gene1,gene2,gene3...
cluster_groups_for_DeepCellSeek <- split(top_genes_for_DeepCellSeek$gene, top_genes_for_DeepCellSeek$cluster)

# 构建输出行
result_lines_for_DeepCellSeek <- sapply(names(cluster_groups_for_DeepCellSeek), function(clust) {
  genes_for_DeepCellSeek <- paste(cluster_groups_for_DeepCellSeek[[clust]], collapse = ",")
  paste0("cluster", clust, ": ", genes_for_DeepCellSeek)
})

# 按数字顺序排序（避免 cluster10 排在 cluster2 前面）
cluster_nums_for_DeepCellSeek <- as.numeric(names(cluster_groups_for_DeepCellSeek))
result_lines_for_DeepCellSeek <- result_lines_for_DeepCellSeek[order(cluster_nums_for_DeepCellSeek)]

# 输出
writeLines(result_lines_for_DeepCellSeek, "output/top_genes_for_DeepCellSeek.txt")

# 在网页中操作选择物种与组织，最多每个簇上传20个基因，具体看网页里的示例，然后将
# top_genes_for_DeepCellSeek.txt 中的结果粘进去提交，比对各个模型的注释结果，选择注释频率最高的细胞类型



#------------------------------------------------------------------------------#
# 承接singleR注释的结果
# 对singleR注释的结果根据mLLM与DeepCellSeek与文章与经典marker进行再修正注释！！！！！！！！
# 为了避免seurat对象冗余，后续使用的seurat对象都是基于singleR的seurat对象
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
# 如果前面用了SCT的标化步骤，先对于后续要用的seurat对象进行如下，避免使用错误的表达数据
# DefaultAssay(seurat对象) <- "RNA"
#------------------------------------------------------------------------------#


#----------------------------------说明----------------------------------------#
# 当你对 Seurat 对象进行下游分析（如 DotPlot, FindAllMarkers, DimPlot 等）时，Seurat 
# 默认会使用当前的 DefaultAssay 所指向的数据层。
# DefaultAssay(seurat对象) <- "RNA" 来避免使用错误的 assay 数据
# 如果你之前运行过 SCTransform()，Seurat 会自动将默认 assay 切换为 "SCT"。
# 而 DotPlot、VlnPlot 等可视化函数会从默认 assay 中提取表达值。
#----------------------------------说明----------------------------------------#



## 8 人工注释===================================================================
# 承接singleR注释,在singleR注释基础上继续编辑

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<特别注意>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#
# 如果注释到Proliferating cell 增殖期的细胞簇，要特别注意，不要单独分为一类，
# 将这个Proliferating cell簇要回归到具体的细胞类型中，一般它只来自一个主要的细胞类型簇，
# 也有特殊情况，但不多，我们统一认为来自一种具体的细胞类型，看它靠近什么细胞类型来具体判断，
# 原因是只是这个簇的top基因高表达了细胞周期的基因，它本质还是具体的细胞类型
# 
# 要对Proliferating cell单独用去除细胞周期基因后的基因合集重新进行FindMarkers，
# 找特异性基因来注释回去原来的细胞类型，之后如果要亚型分析可以作为亚型！
# 具体应对流程模块在下面有👇
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<特别注意>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#

#------------------------------------------------------------------------------#
# 另外可以根据https://cell-ontology.github.io/ 标准化细胞命名
#------------------------------------------------------------------------------#

# 8.1 准备工作------------------------------------------------------------------
# 创建注释副本并标记待验证的簇
celltype2 <- celltype #将celltype 数据框的内容复制到一个新的数据框 celltype2

# 找出 celltype2 数据框中 ClusterID 列值为xxx的所有行，将这些行的 cell 列的值修改为 "clusterxxx"
celltype2$cell[which(celltype2$ClusterID == 0)] <- "cluster0"
celltype2$cell[which(celltype2$ClusterID == 1)] <- "cluster1"
celltype2$cell[which(celltype2$ClusterID == 2)] <- "cluster2"
celltype2$cell[which(celltype2$ClusterID == 3)] <- "cluster3"
celltype2$cell[which(celltype2$ClusterID == 4)] <- "cluster4"
celltype2$cell[which(celltype2$ClusterID == 5)] <- "cluster5"
celltype2$cell[which(celltype2$ClusterID == 6)] <- "cluster6"
celltype2$cell[which(celltype2$ClusterID == 7)] <- "cluster7"
celltype2$cell[which(celltype2$ClusterID == 8)] <- "cluster8"
celltype2$cell[which(celltype2$ClusterID == 9)] <- "cluster9"
celltype2$cell[which(celltype2$ClusterID == 10)] <- "cluster10"
celltype2$cell[which(celltype2$ClusterID == 11)] <- "cluster11"
celltype2$cell[which(celltype2$ClusterID == 12)] <- "cluster12"
celltype2$cell[which(celltype2$ClusterID == 13)] <- "cluster13"
celltype2$cell[which(celltype2$ClusterID == 14)] <- "cluster14"
celltype2$cell[which(celltype2$ClusterID == 15)] <- "cluster15"
celltype2$cell[which(celltype2$ClusterID == 16)] <- "cluster16"
celltype2$cell[which(celltype2$ClusterID == 17)] <- "cluster17"
celltype2$cell[which(celltype2$ClusterID == 18)] <- "cluster18"
celltype2$cell[which(celltype2$ClusterID == 19)] <- "cluster19"
celltype2$cell[which(celltype2$ClusterID == 20)] <- "cluster20"
celltype2$cell[which(celltype2$ClusterID == 21)] <- "cluster21"
celltype2$cell[which(celltype2$ClusterID == 22)] <- "cluster22"
# celltype2$cell[which(celltype2$ClusterID == 23)] <- "cluster23"
# celltype2$cell[which(celltype2$ClusterID == 24)] <- "cluster24"
# celltype2$cell[which(celltype2$ClusterID == 25)] <- "cluster25"
# celltype2$cell[which(celltype2$ClusterID == 26)] <- "cluster26"
# celltype2$cell[which(celltype2$ClusterID == 27)] <- "cluster27"


#------------------------------------------------------------------------------#
#             之后的分析都是针对sc_integer_singleR seurat主体的
#------------------------------------------------------------------------------#


# -------------------------普惠marker注释知识点--------------------------------#
# 上皮细胞（Epithelial Cells）
# ├── 鳞状上皮细胞（Squamous Epithelial Cells）
# │       ├── 角质形成细胞（Keratinocytes）← 复层鳞状上皮，有角化能力
# │       ├── 非角化鳞状上皮细胞（Non-keratinized）← 如口腔、食管
# │       └── 单层鳞状上皮细胞（Simple squamous）← 如血管内皮、肺泡I型
# ├── 柱状上皮细胞（Columnar Epithelial Cells）
# Credential-related configuration omitted from publication copy; configure locally.
# ├── 纤毛上皮细胞（Ciliated Epithelial Cells）
# └── 基底/储备细胞（Basal/Reserve Cells）

# 成熟T 细胞（T cells）
# ├── αβ T 细胞（主流 T 细胞，占 95% 以上）
# │       ├── CD4⁺ T 细胞（辅助性 T 细胞）
# │       │       └── Th1, Th2, Th17, Treg, Tfh 等亚群
# │       └── CD8⁺ T 细胞（细胞毒性 T 细胞）
# │               └── 细胞毒性效应 T 细胞、耗竭 T 细胞等
# └── γδ T 细胞（非主流 T 细胞，占 1-5%）  # γδ T 细胞通常不表达 CD4 和 CD8A/CD8B。
# ├── Vδ1 亚群（主要存在于黏膜组织）
# ├── Vδ2 亚群（主要存在于外周血）
# └── 其他 Vδ 亚群


# 未成熟与成熟细胞，发育谱系注释知识点
# # DN（双阴性）（未成熟）—— CD4- CD8- CD1A 通常在 DN3 阶段才开始表达,但此处不区分，就直接纳入
# # 核心表达模式：CD4 –, CD8A –, CD8B –, CD1A +, CCR7 –
# 
# # DP（双阳性）（未成熟）—— CD4+ CD8+   
# # 核心表达模式：CD4 +, CD8A +, CD8B +, CD1A +, CCR7 –, RAG1 +, DNTT +
# 
# # CD4+ T（CD4单阳性）（成熟）—— CD4+ CD8-
# # 核心表达模式：CD4 +, CD8A –, CD8B –, CD1A –, CCR7 +, IL7R +, SELL +, RAG1 –, DNTT –
# 
# # CD8+ SP（CD8单阳性）（成熟）—— CD4- CD8+
# # 核心表达模式：CD4 –, CD8A +, CD8B +, CD1A –, CCR7 +, IL7R +, SELL +, RAG1 –, DNTT –


# 间充质干细胞
# │
# ┌──────────────────┼──────────────────┐
# │                  │                  │
# 成纤维细胞       周细胞              SMC
# (Fibroblast)    (Pericyte)     (Smooth Muscle)
# │                  │                  │
# │◄────────────────►│◄────────────────►│
# │    互相转化      │    互相转化      │
# │                  │                  │
# └──────────────────┴──────────────────┘
# │
# 肌成纤维细胞 (Myofibroblast)
# (三者均可转化而来)


# -------------------------普惠marker注释知识点--------------------------------#


# ## 5大细胞谱系及对应细胞类型（Table 1）
# ### 1. Epithelial cells（上皮细胞大类）
# AT1         Alveolar type I cells（肺泡I型上皮细胞）    Ager, Hopx, Pdpn 
# AT2         Alveolar type II cells（肺泡II型上皮细胞）  Sftpc, Sftpa1, Sftpb, Lamp3 
# Ciliated    Ciliated cells（纤毛细胞）                  Foxj1, Tuba1a, Cfap53 
# Club        Club cells（克拉拉细胞）                    Scgb1a1, Scgb3a2, Cyp2f2
 
# ### 2. Endothelial cells（内皮细胞）
# aCap     Aerocyte capillary endothelial cells（气-血屏障毛细血管内皮）     Ednrb, Apln, Car4 
# gCap     General capillary endothelial cells（通用毛细血管内皮）           Gpihbp1, Aplnr, Kit 
# Artery   Arterial endothelial cells（动脉内皮）                            Bmx, Dkk2, Gja5 
# Vein     Venous endothelial cells（静脉内皮）                              Vwf, Nr2f2, Ackr1 
# Lymph    Lymphatic endothelial cells（淋巴管内皮）                         Prox1, Lyve1, Flt4
# 
# ### 3. Stromal cells（间质细胞）
# Fib             Fibroblasts（成纤维细胞）                                    Col1a1, Dcn, Lum, Pdgfra 
# Myofib          Myofibroblasts（肌成纤维细胞）                               Acta2, Tagln, Myh11 
# Fibromyo./SMCs  Fibromyocytes / Smooth muscle cells（纤维肌细胞/平滑肌细胞） Acta2, Myh11, Cnn1 
# Pericyte        Pericytes（周细胞）                                          Rgs5, Pdgfrb, Abcc9
# 
# ### 4. Myeloid cells（髓系细胞）
# MΦ       Macrophages（巨噬细胞）        Cd68, Cd163, Csf1r, Lyz2 
# Mono.    Monocytes（单核细胞）          Ly6c2, Fcn1, S100a8, S100a9 
# Neu      Neutrophils（中性粒细胞）      S100a8, S100a9, Csf3r 
# DC       Dendritic cells（树突状细胞）  Flt3, Xcr1, Cd209a, H2-Ab1
# 
# ### 5. Lymphocytes（淋巴细胞）
# T cells     T cells（T细胞）     Cd3d, Cd3e, Cd4/Cd8a 
# B cells     B cells（B细胞）     Cd79a, Ms4a1 (Cd20), Pax5


# 8.2 最终的手动赋值注释----------------------------------------------------
celltype2$cell[which(celltype2$ClusterID == 0)] <- "B cell"
celltype2$cell[which(celltype2$ClusterID == 1)] <- "gCAP"
celltype2$cell[which(celltype2$ClusterID == 2)] <- "Fibroblast" 
celltype2$cell[which(celltype2$ClusterID == 3)] <- "SMC"  # fibroblast
celltype2$cell[which(celltype2$ClusterID == 4)] <- "aCAP"
celltype2$cell[which(celltype2$ClusterID == 5)] <- "T cell" 
celltype2$cell[which(celltype2$ClusterID == 6)] <- "AT2"
celltype2$cell[which(celltype2$ClusterID == 7)] <- "aCAP" 
celltype2$cell[which(celltype2$ClusterID == 8)] <- "Macrophage"
celltype2$cell[which(celltype2$ClusterID == 9)] <- "Club cell"
celltype2$cell[which(celltype2$ClusterID == 10)] <- "Neutrophil"
celltype2$cell[which(celltype2$ClusterID == 11)] <- "Macrophage"
celltype2$cell[which(celltype2$ClusterID == 12)] <- "Macrophage"
celltype2$cell[which(celltype2$ClusterID == 13)] <- "NK cell"
celltype2$cell[which(celltype2$ClusterID == 14)] <- "Lymph"   # Endothelial cell
celltype2$cell[which(celltype2$ClusterID == 15)] <- "Pericyte"
celltype2$cell[which(celltype2$ClusterID == 16)] <- "Fibroblast"
celltype2$cell[which(celltype2$ClusterID == 17)] <- "SMC"
celltype2$cell[which(celltype2$ClusterID == 18)] <- "T cell"
celltype2$cell[which(celltype2$ClusterID == 19)] <- "AT1"
celltype2$cell[which(celltype2$ClusterID == 20)] <- "gCAP"
celltype2$cell[which(celltype2$ClusterID == 21)] <- "Ciliated cell"
celltype2$cell[which(celltype2$ClusterID == 22)] <- "Platelet"
# celltype2$cell[which(celltype2$ClusterID == 23)] <- "CD8+ T"
# celltype2$cell[which(celltype2$ClusterID == 24)] <- "NK"
# celltype2$cell[which(celltype2$ClusterID == 25)] <- "pDC"
# celltype2$cell[which(celltype2$ClusterID == 26)] <- "Cycling T" # Cycling CD8+ T
# celltype2$cell[which(celltype2$ClusterID == 27)] <- "Endothelial cell"


FeaturePlot(sc_integer_singleR,
            features = "S100a8",
            pt.size = 0.05,
            label.size = 3
)



# # # 根据胸腺T细胞发育阶段marker查看注释~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # 胸腺内T细胞注释判断
# FeaturePlot(sc_integer_singleR,   # CD1A 是成熟度的核心分界线。CD1A+ = 未成熟，CD1A- = 已成熟
#             features = "MYH11",
#             pt.size = 0.05,
#             label.size = 3
# )
#
# FeaturePlot(sc_integer_singleR,
#             features = "CD8A",
#             pt.size = 0.05,
#             label.size = 3
# )
# # 胸腺T细胞发育相关基因表达气泡图
# thymus_genes <- c("CD4", "CD8A", "CD8B", "CD1A", "CCR7", "RAG1", "DNTT", "IL7R", "SELL")
# p <- DotPlot(sc_integer_singleR,
#              features = thymus_genes,
#              group.by = "seurat_clusters",
#              dot.scale = 8,
#              cols = c("grey", "red")) +
#   RotatedAxis() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
#   scale_color_gradientn(colours = c("grey","red"))
# print(p)



# 根据fphar-13-1039103文献marker查看如何注释~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 新生小鼠肺脏(BPD)5大细胞谱系标志物列表
Paper_marker_list <- list(
  # ==================== 1. Epithelial cells（上皮细胞） ====================
  "AT1" = c("Ager", "Hopx", "Pdpn"),                 # 肺泡I型上皮细胞
  "AT2" = c("Sftpc", "Sftpa1", "Sftpb", "Lamp3"),     # 肺泡II型上皮细胞
  "Ciliated" = c("Foxj1", "Tuba1a", "Cfap53"),         # 纤毛细胞
  "Club" = c("Scgb1a1", "Scgb3a2", "Cyp2f2"),          # 克拉拉细胞(Club cell)
  # ==================== 2. Endothelial cells（内皮细胞） ====================
  "aCap" = c("Ednrb", "Apln", "Car4"),                 # 气-血屏障毛细血管内皮
  "gCap" = c("Gpihbp1", "Aplnr", "Kit"),                # 通用毛细血管内皮
  "Artery" = c("Bmx", "Dkk2", "Gja5"),                  # 动脉内皮
  "Vein" = c("Vwf", "Nr2f2", "Ackr1"),                  # 静脉内皮
  "Lymphatic" = c("Prox1", "Lyve1", "Flt4"),            # 淋巴管内皮
  # ==================== 3. Stromal cells（间质细胞） ====================
  "Fibroblast" = c("Col1a1", "Dcn", "Lum", "Pdgfra"),   # 成纤维细胞
  "Myofibroblast" = c("Acta2", "Tagln", "Myh11"),       # 肌成纤维细胞
  "Pericyte" = c("Rgs5", "Pdgfrb", "Abcc9"),            # 周细胞
  # ==================== 4. Myeloid cells（髓系细胞） ====================
  "Macrophage" = c("Adgre1", "Cd68", "Csf1r"),          # 巨噬细胞
  "Monocyte" = c("Ly6c2", "Ccr2", "Plac8"),              # 单核细胞
  "Neutrophil" = c("S100a8", "S100a9", "Csf3r"),        # 中性粒细胞
  "DC" = c("Flt3", "Xcr1", "Cd209a", "H2-Ab1"),         # 树突状细胞
  # ==================== 5. Lymphocytes（淋巴细胞） ====================
  "T cell" = c("Cd3d", "Cd3e", "Cd8a", "Cd4"),          # T细胞
  "B cell" = c("Cd79a", "Ms4a1", "Cd19", "Pax5"),       # B细胞
  "NK cell" = c("Klrd1", "Nkg7", "Gzma")                # NK细胞
)
p <- DotPlot(sc_integer_singleR,
             features = Paper_marker_list,  # 直接传入列表，实现分组
             group.by = "seurat_clusters",         # 【重要】指定细胞分组的元数据列名
             dot.scale = 8,                 # 点的大小
             cols = c("grey", "red")) +     # 设置颜色渐变：低表达灰色，高表达红色
  RotatedAxis() +                           # 使用 Seurat 自带的旋转函数，等同于 ggpubr::rotate_x_text
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + # 进一步调整X轴文字位置防止遮挡
  scale_color_gradientn(colours = c("grey","red")) # 确保颜色映射正确
print(p)


# 查看对象当前的默认Assay
DefaultAssay(sc_integer_singleR)


# #  cluster0---------------------------------------------------------------------
# # 寻找簇0的差异表达基因
# markers_cluster0 <- FindMarkers(sc_integer_singleR,
#                                 min.pct = 0.3, # 基因在至少30%的细胞中表达,提高筛选的参数
#                                 ident.1 = "cluster0"
#                                 )
# 
# # 筛选和排序标记基因
# markers_cluster0_1 <- markers_cluster0 %>%
#   dplyr::filter(p_val_adj < 0.05) %>%          # 筛选显著基因
#   dplyr::mutate(pct = pct.1 - pct.2) %>%       # 计算表达比例差异
#   dplyr::arrange(desc(avg_log2FC), desc(pct))  # 按表达差异排序
# 
# # 可视化前10个标记基因
# png("output/C0.png")
# DotPlot(sc_integer_singleR, assay = "RNA",
#         features = rownames(markers_cluster0_1)[1:10],
#         group.by = "seurat_clusters") +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# dev.off()
# 
# # 输出前10个标记基因名字(用于Cellmarker2)
# cluster0_top_10_genes <- rownames(markers_cluster0_1)[1:10]
# cat(paste(cluster0_top_10_genes, collapse = ","), "\n")
# 
# # FOLR2,C1QB,C1QC,GPR34,C3,SELENOP,C1QA,SLCO2B1,MS4A7,FCGBP
# # 我在对我在进行人类颈动脉粥样硬化斑块组织(carotid athersclerotic plaques)的单细胞转录组细胞注释，我有以下10个top marker 基因：(FOLR2,C1QB,C1QC,GPR34,C3,SELENOP,C1QA,SLCO2B1,MS4A7,FCGBP)
# 
# 
# FeaturePlot(sc_integer_singleR,
#             features = "FOLR2",
#             pt.size = 0.05,
#             label.size = 3
# )
# 
# # 最终定夺注释
# celltype2$cell[which(celltype2$ClusterID == 0)] <- "Macrophage"
# # Macrophage （巨噬细胞）


# # cluster1----------------------------------------------------------------------
# # 寻找簇1的差异表达基因
# markers_cluster1 <- FindMarkers(sc_integer_singleR,
#                                 min.pct = 0.3, # 基因在至少30%的细胞中表达,提高筛选的参数
#                                 ident.1 = "cluster1"
#                                 )
# 
# # 筛选和排序标记基因
# markers_cluster1_1 <- markers_cluster1 %>%
#   dplyr::filter(p_val_adj < 0.05) %>%          # 筛选显著基因
#   dplyr::mutate(pct = pct.1 - pct.2) %>%       # 计算表达比例差异
#   dplyr::arrange(desc(avg_log2FC), desc(pct))  # 按表达差异排序
# 
# # 可视化前10个标记基因
# png("output/C1.png")
# DotPlot(sc_integer_singleR, assay = "RNA",
#         features = rownames(markers_cluster1_1)[1:10],
#         group.by = "seurat_clusters") +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# dev.off()
# 
# # 输出前10个标记基因名字(用于Cellmarker2)
# cluster1_top_10_genes <- rownames(markers_cluster1_1)[1:10]
# cat(paste(cluster1_top_10_genes, collapse = ","), "\n")
# 
# # CD40LG,LTB,KLRB1,IL7R,PRKCQ-AS1,SPOCK2,CD6,TRAC,TRBC2,IL32 
# # 我在对我在进行人类颈动脉粥样硬化斑块组织(carotid athersclerotic plaques)的单细胞转录组细胞注释，我有以下10个top marker 基因：(CD40LG,LTB,KLRB1,IL7R,PRKCQ-AS1,SPOCK2,CD6,TRAC,TRBC2,IL32)
# 
# 
# FeaturePlot(sc_integer_singleR,
#             features = "CD40LG",
#             pt.size = 0.05,
#             label.size = 3
# )
# 
# # 最终定夺注释
# celltype2$cell[which(celltype2$ClusterID == 1)] <- "CD4+ T cell"
# 
# 
# # cluster2----------------------------------------------------------------------
# # 寻找簇2的差异表达基因
# markers_cluster2 <- FindMarkers(sc_integer_singleR,
#                                 min.pct = 0.3, # 基因在至少30%的细胞中表达,提高筛选的参数
#                                 ident.1 = "cluster2"
#                                 )
# 
# # 筛选和排序标记基因
# markers_cluster2_1 <- markers_cluster2 %>%
#   dplyr::filter(p_val_adj < 0.05) %>%          # 筛选显著基因
#   dplyr::mutate(pct = pct.1 - pct.2) %>%       # 计算表达比例差异
#   dplyr::arrange(desc(avg_log2FC), desc(pct))  # 按表达差异排序
# 
# # 可视化前10个标记基因
# png("output/C2.png")
# DotPlot(sc_integer_singleR, assay = "RNA",
#         features = rownames(markers_cluster2_1)[1:10],
#         group.by = "seurat_clusters") +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# dev.off()
# 
# # 输出前10个标记基因名字(用于Cellmarker2)
# cluster2_top_10_genes <- rownames(markers_cluster2_1)[1:10]
# cat(paste(cluster2_top_10_genes, collapse = ","), "\n")
# 
# # KLRD1,NKG7,CD8B,CD8A,XCL2,PRF1,KLRK1,GZMH,CTSW,CCL5 
# # 我在对我在进行人类颈动脉粥样硬化斑块组织(carotid athersclerotic plaques)的单细胞转录组细胞注释，我有以下10个top marker 基因：(KLRD1,NKG7,CD8B,CD8A,XCL2,PRF1,KLRK1,GZMH,CTSW,CCL5)
# 
# FeaturePlot(sc_integer_singleR,
#             features = "KLRD1",
#             pt.size = 0.05,
#             label.size = 3
# )
# 
# # 最终定夺注释
# celltype2$cell[which(celltype2$ClusterID == 2)] <- "CD8+ T cell"


# # cluster3----------------------------------------------------------------------
# # 寻找簇3的差异表达基因
# markers_cluster3 <- FindMarkers(sc_integer_singleR,
#                                 min.pct = 0.3, # 基因在至少30%的细胞中表达,提高筛选的参数
#                                 ident.1 = "cluster3",
#                                 )
# 
# # 筛选和排序标记基因
# markers_cluster3_1 <- markers_cluster3 %>%
#   dplyr::filter(p_val_adj < 0.05) %>%          # 筛选显著基因
#   dplyr::mutate(pct = pct.1 - pct.2) %>%       # 计算表达比例差异
#   dplyr::arrange(desc(avg_log2FC), desc(pct))  # 按表达差异排序
# 
# # 可视化前10个标记基因
# png("output/C3.png")
# DotPlot(sc_integer_singleR, assay = "RNA",
#         features = rownames(markers_cluster3_1)[1:10],
#         group.by = "seurat_clusters") +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# dev.off()
# 
# # 输出前10个标记基因名字(用于Cellmarker2)
# cluster3_top_10_genes <- rownames(markers_cluster3_1)[1:10]
# cat(paste(cluster3_top_10_genes, collapse = ","), "\n")
# 
# # GNG4,HSPB6,RBP1,HAND2,TAGLN,MYL9,TWIST2,ULK4,S100A16,RORB 
# # 我在对我在进行人类蜕膜组织(decidua)的单细胞转录组细胞注释，我有以下10个top marker 基因：(GNG4,HSPB6,RBP1,HAND2,TAGLN,MYL9,TWIST2,ULK4,S100A16,RORB)
# 
# FeaturePlot(sc_integer_singleR,
#             features = "LUM",
#             pt.size = 0.05,
#             label.size = 3
# )
# 
# # 最终定夺注
# celltype2$cell[which(celltype2$ClusterID == 3)] <- "DSC"



# # cluster4--------------------------------------------------------------------
# # 寻找簇4的差异表达基因
# markers_cluster4 <- FindMarkers(sc_integer_singleR,
#                                 min.pct = 0.3, # 基因在至少30%的细胞中表达,提高筛选的参数
#                                 ident.1 = "cluster4",
#                                 )
# 
# # 筛选和排序标记基因
# markers_cluster4_1 <- markers_cluster4 %>%
#   dplyr::filter(p_val_adj < 0.05) %>%          # 筛选显著基因
#   dplyr::mutate(pct = pct.1 - pct.2) %>%       # 计算表达比例差异
#   dplyr::arrange(desc(avg_log2FC), desc(pct))  # 按表达差异排序
# 
# # 可视化前10个标记基因
# png("output/C4.png")
# DotPlot(sc_integer_singleR, assay = "RNA",
#         features = rownames(markers_cluster4_1)[1:10],
#         group.by = "seurat_clusters") +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# dev.off()
# 
# # 输出前10个标记基因名字(用于Cellmarker2)
# cluster4_top_10_genes <- rownames(markers_cluster4_1)[1:10]
# cat(paste(cluster4_top_10_genes, collapse = ","), "\n")
# 
# # Mob3b,Anln,Klk6,Enpp6,Cyp27a1,Cntn3,Mapre2,A330015K06Rik,Cpox,Elovl7 
# 
# # 我在进行小鼠脊髓组织spinal的单细胞转录组细胞注释，我有以下10个top marker 基因：(Mob3b,Anln,Klk6,Enpp6,Cyp27a1,Cntn3,Mapre2,A330015K06Rik,Cpox,Elovl7 )
# 
# FeaturePlot(sc_integer_singleR,
#             features = "Mob3b",
#             pt.size = 0.05,
#             label.size = 3
# )
# 
# 
# # 最终定夺注释
# celltype2$cell[which(celltype2$ClusterID == 4)] <- "Kupffer cell"


# # cluster5----------------------------------------------------------------------
# # 寻找簇5的差异表达基因
# markers_cluster5 <- FindMarkers(sc_integer_singleR,
#                                 min.pct = 0.3, # 基因在至少30%的细胞中表达,提高筛选的参数
#                                 ident.1 = "cluster5",
#                                 )
# 
# # 筛选和排序标记基因
# markers_cluster5_1 <- markers_cluster5 %>%
#   dplyr::filter(p_val_adj < 0.05) %>%          # 筛选显著基因
#   dplyr::mutate(pct = pct.1 - pct.2) %>%       # 计算表达比例差异
#   dplyr::arrange(desc(avg_log2FC), desc(pct))  # 按表达差异排序
# 
# # 可视化前10个标记基因
# png("output/C5.png")
# DotPlot(sc_integer_singleR, assay = "RNA",
#         features = rownames(markers_cluster5_1)[1:10],
#         group.by = "seurat_clusters") +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# dev.off()
# 
# # 输出前10个标记基因名字(用于Cellmarker2)
# cluster5_top_10_genes <- rownames(markers_cluster5_1)[1:10]
# cat(paste(cluster5_top_10_genes, collapse = ","), "\n")
# 
# # Fermt1,Sftpd,Pkhd1,Dcdc2a,Onecut1,Scara3,Plet1,Tff2,Adam33,Sox9 
# # 我在对我在进行人类的颈动脉carotid artery组织的单细胞转录组细胞注释，我有以下10个top marker 基因：(S100A12,S100A8,APOBEC3A,MCEMP1,FCN1,CD300E,EREG,S100A9,AQP9,G0S2 )
# 
# FeaturePlot(sc_integer_singleR,
#             features = "CD163",
#             pt.size = 0.05,
#             label.size = 3
# )
# 
# # cluster5 KEGG富集~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# cluster5_top_50_genes <- rownames(markers_cluster5_1)[1:50]
# cluster5_enriched <- enrichr(cluster5_top_50_genes, "KEGG_2026") # 或者 KEGG_2021_Human/KEGG_2019_Mouse    可选用的数据库参考https://maayanlab.cloud/Enrichr/#libraries
# cluster5_kegg_sig <- cluster5_enriched[["KEGG_2026"]] %>% filter(Adjusted.P.value < 0.05)
# 
# if (nrow(cluster5_kegg_sig) > 0) {
#   cluster5_top_10_KEGG <- head(cluster5_kegg_sig[order(cluster5_kegg_sig$Adjusted.P.value), ], 10)
#   cat("Cluster5 前10个KEGG通路:\n")
#   # 直接打印序号和通路名称
#   for (i in 1:nrow(cluster5_top_10_KEGG)) {
#     cat(i, ". ", cluster5_top_10_KEGG$Term[i], "\n", sep="")
#   }
# } else {
#   cat("未发现显著富集的KEGG通路\n")
# }
# 
# # 然后把富集的通路与marker一并丢给大模型进行辅助注释
# 
# # 我在对我在进行人类的颈动脉carotid artery组织的单细胞转录组细胞注释，我有以下10个top marker 基因：(S100A12,S100A8,APOBEC3A,MCEMP1,FCN1,CD300E,EREG,S100A9,AQP9,G0S2 )
# # 并且50个top marker KEGG富集到...这个细胞类型是什么
# 
# 
# # 最终定夺注释
# celltype2$cell[which(celltype2$ClusterID == 5)] <- "Macrophage"


# # cluster6----------------------------------------------------------------------
# # 寻找簇6的差异表达基因
# markers_cluster6 <- FindMarkers(sc_integer_singleR,
#                                 min.pct = 0.3, # 基因在至少30%的细胞中表达,提高筛选的参数
#                                 ident.1 = "cluster6",
#                                )
# 
# # 筛选和排序标记基因
# markers_cluster6_1 <- markers_cluster6 %>%
#   dplyr::filter(p_val_adj < 0.05) %>%          # 筛选显著基因
#   dplyr::mutate(pct = pct.1 - pct.2) %>%       # 计算表达比例差异
#   dplyr::arrange(desc(avg_log2FC), desc(pct))  # 按表达差异排序
# 
# # 可视化前10个标记基因
# png("output/C6.png")
# DotPlot(sc_integer_singleR, assay = "RNA",
#         features = rownames(markers_cluster6_1)[1:10],
#         group.by = "seurat_clusters") +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# dev.off()
# 
# # 输出前10个标记基因名字(用于Cellmarker2)
# cluster6_top_10_genes <- rownames(markers_cluster6_1)[1:10]
# cat(paste(cluster6_top_10_genes, collapse = ","), "\n")
# 
# # MFAP5,CHRDL1,FBLN1,SFRP2,APOD,PI16,DCN,SFRP4,GAS1,C3 
# 
# # 我在对我在进行人类的颈动脉carotid artery组织的单细胞转录组细胞注释，我有以下10个top marker 基因：(MFAP5,CHRDL1,FBLN1,SFRP2,APOD,PI16,DCN,SFRP4,GAS1,C3 )
# 
# FeaturePlot(sc_integer_singleR,
#             features = "PI16",
#             pt.size = 0.05,
#             label.size = 3
# )
# 
# 
# # 最终定夺注
# celltype2$cell[which(celltype2$ClusterID == 6)] <- "Fibroblast"  #


# # cluster7----------------------------------------------------------------------
# # 寻找簇7的差异表达基因
# markers_cluster7 <- FindMarkers(sc_integer_singleR,
#                                 min.pct = 0.3, # 基因在至少30%的细胞中表达,提高筛选的参数
#                                 ident.1 = "cluster7",
# )
# 
# # 筛选和排序标记基因
# markers_cluster7_1 <- markers_cluster7 %>%
#   dplyr::filter(p_val_adj < 0.05) %>%          # 筛选显著基因
#   dplyr::mutate(pct = pct.1 - pct.2) %>%       # 计算表达比例差异
#   dplyr::arrange(desc(avg_log2FC), desc(pct))  # 按表达差异排序
# 
# # 可视化前10个标记基因
# png("output/C7.png")
# DotPlot(sc_integer_singleR, assay = "RNA",
#         features = rownames(markers_cluster7_1)[1:10],
#         group.by = "seurat_clusters") +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# dev.off()
# 
# # 输出前10个标记基因名字
# cluster7_top_10_genes <- rownames(markers_cluster7_1)[1:10]
# cat(paste(cluster7_top_10_genes, collapse = ","), "\n")
# 
# # ITGAX,ARAP2,ACAP1,ARHGEF1,ITGAD,SAMD3,TRG-AS1,NKTR,EVL,RUNX3
# # 我在进行人类蜕膜组织(decidua)的单细胞转录组细胞注释，我有以下10个top marker 基因：(ITGAX,ARAP2,ACAP1,ARHGEF1,ITGAD,SAMD3,TRG-AS1,NKTR,EVL,RUNX3 )
# 
# 
# celltype2$cell[which(celltype2$ClusterID == 7)] <- "dNK"    # Decidual NK cell       蜕膜自然杀伤细胞


# # cluster8----------------------------------------------------------------------
# # 寻找簇8的差异表达基因
# markers_cluster8 <- FindMarkers(sc_integer_singleR,
#                                 min.pct = 0.3, # 基因在至少30%的细胞中表达,提高筛选的参数
#                                 ident.1 = "cluster8",
# )
# 
# # 筛选和排序标记基因
# markers_cluster8_1 <- markers_cluster8 %>%
#   dplyr::filter(p_val_adj < 0.05) %>%          # 筛选显著基因
#   dplyr::mutate(pct = pct.1 - pct.2) %>%       # 计算表达比例差异
#   dplyr::arrange(desc(avg_log2FC), desc(pct))  # 按表达差异排序
# 
# # 可视化前10个标记基因
# png("output/C8.png")
# DotPlot(sc_integer_singleR, assay = "RNA",
#         features = rownames(markers_cluster8_1)[1:10],
#         group.by = "seurat_clusters") +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# dev.off()
# 
# # 输出前10个标记基因名字
# cluster8_top_10_genes <- rownames(markers_cluster8_1)[1:10]
# cat(paste(cluster8_top_10_genes, collapse = ","), "\n")
# 
# # Vcan,Pdgfra,Megf11,3110035E14Rik,C1ql1,Ppfibp1,Itga9,Gm4876,Tnr,Sema3d 
# # 我在进行小鼠脊髓组织spinal的单细胞转录组细胞注释，我有以下10个top marker 基因：(Vcan,Pdgfra,Megf11,3110035E14Rik,C1ql1,Ppfibp1,Itga9,Gm4876,Tnr,Sema3d )
# 
# 
# celltype2$cell[which(celltype2$ClusterID == 8)] <- "OPC" 



# # cluster8----------------------------------------------------------------------
# # 寻找簇8的差异表达基因
# markers_cluster8 <- FindMarkers(sc_integer_singleR,
#                                 min.pct = 0.3, # 基因在至少30%的细胞中表达,提高筛选的参数
#                                 ident.1 = "cluster8",
#                               )
# 
# # 筛选和排序标记基因
# markers_cluster8_1 <- markers_cluster8 %>%
#   dplyr::filter(p_val_adj < 0.05) %>%          # 筛选显著基因
#   dplyr::mutate(pct = pct.1 - pct.2) %>%       # 计算表达比例差异
#   dplyr::arrange(desc(avg_log2FC), desc(pct))  # 按表达差异排序
# 
# # 可视化前10个标记基因
# png("output/C8.png")
# DotPlot(sc_integer_singleR, assay = "RNA",
#         features = rownames(markers_cluster8_1)[1:10],
#         group.by = "seurat_clusters") +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# dev.off()
# 
# # 输出前10个标记基因名字
# cluster8_top_10_genes <- rownames(markers_cluster8_1)[1:10]
# cat(paste(cluster8_top_10_genes, collapse = ","), "\n")
# 
# # DLGAP5,SPC25,HIST1H3B,UBE2C,ASPM,KIF15,KIF2C,CDCA8,TYMS,HIST1H1B 
# # 我在进行人类蜕膜组织(decidua)的单细胞转录组细胞注释，我有以下10个top marker 基因：(DLGAP5,SPC25,HIST1H3B,UBE2C,ASPM,KIF15,KIF2C,CDCA8,TYMS,HIST1H1B )
# 
# # celltype2$cell[which(celltype2$ClusterID == 8)] <- ""
# 
# # 判断出是增殖期细胞（proliferating cells）   cycling cell 
# 
# 
# #-------------------------------------------------------------------------------#
# # <<<<<<<<<<<<<<<<<<<<<<<<<增殖期细胞来了，特殊流程>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# 
# # 1. 定义标准细胞周期基因列表（人类版）无G1期 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # a) Seurat 默认
# s_g2m <- c(cc.genes$s.genes, cc.genes$g2m.genes)
# 
# # b) GO:0000278 (mitotic cell cycle)
# go_cc <- bitr("GO:0000278", "GO", "SYMBOL", org.Hs.eg.db)$SYMBOL
# 
# # c) MSigDB Hallmark: HALLMARK_G2M_CHECKPOINT + HALLMARK_E2F_TARGETS
# #    （这是目前最全面、经实验验证的周期基因集）
# #    来源: https://www.gsea-msigdb.org/gsea/msigdb/cards/HALLMARK_G2M_CHECKPOINT.html
# #         https://www.gsea-msigdb.org/gsea/msigdb/cards/HALLMARK_E2F_TARGETS.html
# 
# hallmark_g2m <- c(
#   "AURKA", "AURKB", "BUB1", "BUB1B", "CCNA2", "CCNB1", "CCNB2", "CDC20", "CDC25C",
#   "CDK1", "CENPA", "CENPF", "DLGAP5", "ESPL1", "KIF11", "KIF20A", "KIF23", "MAD2L1",
#   "MKI67", "NCAPD2", "NCAPG", "NDC80", "NUF2", "PBK", "PLK1", "PKMYT1", "RAD51",
#   "RRM2", "TOP2A", "TPX2", "TROAP", "TYMS", "UBE2C", "ZWINT"
# )
# 
# hallmark_e2f <- c(
#   "ASPM", "BIRC5", "BRIP1", "CCNE1", "CDC6", "CDCA2", "CDCA3", "CDCA5", "CDCA8",
#   "CEP55", "CKS1B", "CKS2", "DEPDC1", "E2F1", "E2F2", "E2F7", "E2F8", "ECT2",
#   "FEN1", "FOXM1", "GINS2", "GMNN", "HJURP", "HMGB2", "KIF2C", "MCM2", "MCM3",
#   "MCM4", "MCM5", "MCM6", "MCM7", "MELK", "MND1", "ORC6", "PCNA", "PCLAF", "PKMYT1",
#   "RAD51", "RFC2", "RFC3", "RFC4", "RFC5", "RRM2", "RSRC2", "RPA2", "SMC4", "SPC25",
#   "TIMELESS", "TK1", "TOP2A", "TYMS", "UBE2C", "UNG", "USP1", "WDR76"
# )
# 
# # d) 手动补充常见但可能遗漏的周期/有丝分裂基因（包括组蛋白等）
# manual_cc <- c(
#   "SPC25", "PCLAF", "ASPM", "BUB1B", "FAM111B",
#   "HIST1H1B", "HIST1H3B", "NCAPG", "MND1", "CEP55",
#   "RRM2", "MKI67", "TOP2A", "UBE2C", "DLGAP5", "CDK1",
#   "KIF15", "KIF20A", "AURKB", "CCNB1", "CDC20",
#   "KIFC1", "TROAP", "KNL1", "CCNA2", "ESCO2", "TK1",
#   "FOXM1", "DIAPH3", "MELK", "ZWINT",
#   "HIST1H2AA", "HIST1H2AB", "HIST1H2AC", "HIST1H2BA", "HIST1H2BB",
#   "HIST1H3A", "HIST1H3C", "HIST1H3D", "HIST1H4A", "HIST1H4B",
#   "PLK1", "PBK", "NUF2", "CENPF", "HMGB2", "CKS2", "GINS2",
#   "CDKN3", "SHCBP1", "PRC1", "POC1A", "CENPM", "NCAPH", "UBE2T",
#   "ASF1B", "SGO1", "CHEK1", "FANCI", "CENPN", "TCF19",
#   "BRCA1", "SGO2", "ARHGAP11A", "STMN1", "NCAPG2", "CENPU", "PARPBP",
#   "HMGB3", "RACGAP1", "GGH", "DHFR", "C21orf58", "CENPK", "CIP2A",
#   "PTTG1", "FBP1", "CENPH", "EZH2", "ATAD5", "HIST1H2AG", "FBXO5",
#   "CENPS", "LMNB1", "KNSTRN", "PAQR4", "SMC2", "HIST1H4C",
#   "HIST1H2AH", "HIST1H1D", "LIG1", "SPDL1", "BRCA2", "CCDC34", 
#   "KPNA2", "CHAF1A", "UBE2S", "DNAJC9",
#   "H2AFX", "HIST1H1E", "PHF19", "TMEM106C", "HES6", "HIRIP3", 
#   "DUT", "VRK1", "ITGB3BP", "HIST1H1C", "H2AFZ", "MIS18A", 
#   "DTYMK", "LRR1", "ACAT2", "NCAPD3", "NSD2", "CDCA4", 
#   "BARD1", "BCL2L12", "NUDT8", "AP000547.3", "TTF2", "EBP", 
#   "ACY3", "ACOT7", "HMGN2", "WDR34", "RPA3", "SASS6", 
#   "C19orf48", "LMNB2", "CDKN2A", "ITM2A", "MTHFD1", 
#   "HIST2H2AC", "TESC", "CENPJ", "NUCB2", "H2AFV", "SMC1A", 
#   "DBF4", "MSH6", "IFI27L1"
# )
# 
# # 合并所有来源并去重
# cell_cycle_genes <- unique(c(s_g2m, go_cc, hallmark_g2m, hallmark_e2f, manual_cc))
# 
# # 只保留存在于你数据中的
# cell_cycle_genes <- cell_cycle_genes[cell_cycle_genes %in% rownames(sc_integer_singleR)]
# 
# 
# # 2. 对 cluster8 重新找 marker，排除细胞周期基因~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# markers_cluster8_clean <- FindMarkers(sc_integer_singleR,
#                                        ident.1 = "cluster8",
#                                        min.pct = 0.3
# )
# 
# # 过滤掉细胞周期基因
# markers_cluster8_clean <- markers_cluster8_clean[!rownames(markers_cluster8_clean) %in% cell_cycle_genes, ]
# 
# 
# # 3. 筛选显著且高表达的非周期 marker ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# markers_cluster8_final <- markers_cluster8_clean %>%
#   dplyr::filter(p_val_adj < 0.05) %>%          # 筛选显著基因
#   dplyr::mutate(pct = pct.1 - pct.2) %>%       # 计算表达比例差异
#   dplyr::arrange(desc(avg_log2FC), desc(pct))  # 按表达差异排序
# 
# # 查看 top 10 非周期 marker
# cat("Cluster8 non-cell-cycle top markers:\n")
# print(rownames(markers_cluster8_final)[1:10])
# 
# # "MZB1"   "GNLY"   "GZMA"   "RANBP1" "TEX30"  "RAB33A" "TOPBP1" "DDX39A" "OXCT1"  "MZT1"  
# # 我在进行人类蜕膜组织(decidua)的单细胞转录组细胞注释，我有以下10个top marker 基因：("MZB1"   "GNLY"   "GZMA"   "RANBP1" "TEX30"  "RAB33A" "TOPBP1" "DDX39A" "OXCT1"  "MZT1"  ) 
# 
# 
# FeaturePlot(sc_integer_singleR,
#             features = "GNLY",
#             pt.size = 0.05,
#             label.size = 3
# )
# 
# 
# # 最后定义细胞类型
# celltype2$cell[which(celltype2$ClusterID == 8)] <- "Cycling dNK"

#-------------------------------------------------------------------------------#

  
# # cluster9----------------------------------------------------------------------
# # 寻找簇9的差异表达基因
# markers_cluster9 <- FindMarkers(sc_integer_singleR,
#                                 min.pct = 0.3, # 基因在至少30%的细胞中表达,提高筛选的参数
#                                 ident.1 = "cluster9"
#                                )
# 
# # 筛选和排序标记基因
# markers_cluster9_1 <- markers_cluster9 %>%
#   dplyr::filter(p_val_adj < 0.05) %>%          # 筛选显著基因
#   dplyr::mutate(pct = pct.1 - pct.2) %>%       # 计算表达比例差异
#   dplyr::arrange(desc(avg_log2FC), desc(pct))  # 按表达差异排序
# 
# # 可视化前10个标记基因
# png("output/C9.png")
# DotPlot(sc_integer_singleR, assay = "RNA",
#         features = rownames(markers_cluster9_1)[1:10],
#         group.by = "seurat_clusters") +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# dev.off()
# 
# # 输出前10个标记基因名字
# cluster9_top_10_genes <- rownames(markers_cluster9_1)[1:10]
# cat(paste(cluster9_top_10_genes, collapse = ","), "\n")
# 
# # CCND1,ANTXR1,SYNPO,MYH10,NDUFA4L2,COPZ2,MYH11,PALLD,SUGCT,ITGBL1
# # 我在对我在进行人类肝细胞癌组织(hepatocellular carcinoma)的单细胞转录组细胞注释，我有以下10个top marker 基因：(MKI67,RRM2,SHCBP1,GTSE1,AURKB,TYMS,ASF1B,ASPM,KIAA0101,CCNA2)
# 
# celltype2$cell[which(celltype2$ClusterID == 9)] <- "Cycling hepatocyte"
# 
# 
# FeaturePlot(sc_integer_singleR,
#             features = "APOA1",
#             pt.size = 0.05,
#             label.size = 3
# )


# # cluster10---------------------------------------------------------------------
# # 寻找簇10的差异表达基因
# markers_cluster10 <- FindMarkers(sc_integer_singleR,
#                                  min.pct = 0.3, # 基因在至少30%的细胞中表达,提高筛选的参数
#                                  ident.1 = "cluster10",
#                                  )
# 
# # 筛选和排序标记基因
# markers_cluster10_1 <- markers_cluster10 %>%
#   dplyr::filter(p_val_adj < 0.05) %>%          # 筛选显著基因
#   dplyr::mutate(pct = pct.1 - pct.2) %>%       # 计算表达比例差异
#   dplyr::arrange(desc(avg_log2FC), desc(pct))  # 按表达差异排序
# 
# # 可视化前10个标记基因
# png("output/C10.png")
# DotPlot(sc_integer_singleR, assay = "RNA",
#         features = rownames(markers_cluster10_1)[1:10],
#         group.by = "seurat_clusters") +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# dev.off()
# 
# # 输出前10个标记基因名字(用于Cellmarker2)
# cluster10_top_10_genes <- rownames(markers_cluster10_1)[1:10]
# cat(paste(cluster10_top_10_genes, collapse = ","), "\n")
# 
# 
# # CD1E,CD1C,FCER1A,LGALS2,CLEC10A,PKIB,IL1R2,HLA-DQB2,CPVL,HLA-DQA1 
# # 我在对我在进行人类的颈动脉carotid artery组织的单细胞转录组细胞注释，我有以下10个top marker 基因，我有以下10个top marker 基因：(CD1E,CD1C,FCER1A,LGALS2,CLEC10A,PKIB,IL1R2,HLA-DQB2,CPVL,HLA-DQA1)
# 
# 
# FeaturePlot(sc_integer_singleR,
#             features = "Ly6c1",
#             pt.size = 0.05,
#             label.size = 3
# )
# 
# celltype2$cell[which(celltype2$ClusterID == 10)] <- "cDC"


# # cluster11---------------------------------------------------------------------
# # 寻找簇11的差异表达基因
# markers_cluster11 <- FindMarkers(sc_integer_singleR,
#                                  min.pct = 0.3, # 基因在至少30%的细胞中表达,提高筛选的参数
#                                  ident.1 = "cluster11",
# )
# 
# # 筛选和排序标记基因
# markers_cluster11_1 <- markers_cluster11 %>%
#   dplyr::filter(p_val_adj < 0.05) %>%          # 筛选显著基因
#   dplyr::mutate(pct = pct.1 - pct.2) %>%       # 计算表达比例差异
#   dplyr::arrange(desc(avg_log2FC), desc(pct))  # 按表达差异排序
# 
# 
# # 可视化前10个标记基因
# png("output/C11.png")
# DotPlot(sc_integer_singleR, assay = "RNA",
#         features = rownames(markers_cluster11_1)[1:10],
#         group.by = "seurat_clusters") +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# dev.off()
# 
# # 输出前10个标记基因名字(用于Cellmarker2)
# cluster11_top_10_genes <- rownames(markers_cluster11_1)[1:10]
# cat(paste(cluster11_top_10_genes, collapse = ","), "\n")
# 
# # SPC25,UBE2C,DLGAP5,PKMYT1,MKI67,RRM2,TOP2A,PCLAF,KIF15,ASPM
# # 我在进行人类蜕膜组织(decidua)的单细胞转录组细胞注释，我有以下10个top marker 基因：(AC211486.5,NOTUM,CORO6,PAPPA2,FAT2,CLDN19,GPR78,KISS1R,AC211433.1,COL27A1 )
# 
# # 最后定义细胞类型
# celltype2$cell[which(celltype2$ClusterID == 11)] <- "EVT"
# 
# FeaturePlot(sc_integer_singleR,
#             features = "DKK1",
#             pt.size = 0.05,
#             label.size = 3
# )
# 
# # Mmarker :DKK1

# # cluster12---------------------------------------------------------------------
# # 寻找簇12的差异表达基因
# markers_cluster12 <- FindMarkers(sc_integer_singleR,
#                                 min.pct = 0.3, # 基因在至少30%的细胞中表达,提高筛选的参数
#                                 ident.1 = "cluster12",
#                                 )
# 
# # 筛选和排序标记基因
# markers_cluster12_1 <- markers_cluster12 %>%
#   dplyr::filter(p_val_adj < 0.05) %>%          # 筛选显著基因
#   dplyr::mutate(pct = pct.1 - pct.2) %>%       # 计算表达比例差异
#   dplyr::arrange(desc(avg_log2FC), desc(pct))  # 按表达差异排序
# 
# # 可视化前10个标记基因
# png("output/C12.png")
# DotPlot(sc_integer_singleR, assay = "RNA",
#         features = rownames(markers_cluster12_1)[1:10],
#         group.by = "seurat_clusters") +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# dev.off()
# 
# # 输出前10个标记基因名字(用于Cellmarker2)
# cluster12_top_10_genes <- rownames(markers_cluster12_1)[1:10]
# cat(paste(cluster12_top_10_genes, collapse = ","), "\n")
# 
# # Scn3a,Slc38a5,Slc36a2,Col13a1,Slit2,Limch1,Npnt,Inmt,Robo2,Ptprd  
# # 我在进行小鼠肺组织(Lung)的单细胞转录组细胞注释，我有以下10个top marker 基因：(Scn3a,Slc38a5,Slc36a2,Col13a1,Slit2,Limch1,Npnt,Inmt,Robo2,Ptprd )
# 
# celltype2$cell[which(celltype2$ClusterID == 12)] <- "dM"



# # cluster13---------------------------------------------------------------------
# # 寻找簇13的差异表达基因
# markers_cluster13 <- FindMarkers(sc_integer_singleR,
#                                  min.pct = 0.3, # 基因在至少30%的细胞中表达,提高筛选的参数
#                                  ident.1 = "cluster13",
#                                  )
# 
# # 筛选和排序标记基因
# markers_cluster13_1 <- markers_cluster13 %>%
#   dplyr::filter(p_val_adj < 0.05) %>%          # 筛选显著基因
#   dplyr::mutate(pct = pct.1 - pct.2) %>%       # 计算表达比例差异
#   dplyr::arrange(desc(avg_log2FC), desc(pct))  # 按表达差异排序
# 
# # 可视化前10个标记基因
# png("output/C13.png")
# DotPlot(sc_integer_singleR, assay = "RNA",
#         features = rownames(markers_cluster13_1)[1:10],
#         group.by = "seurat_clusters") +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# dev.off()
# 
# # 输出前10个标记基因名字(用于Cellmarker2)
# cluster13_top_10_genes <- rownames(markers_cluster13_1)[1:10]
# cat(paste(cluster13_top_10_genes, collapse = ","), "\n")
# 
# # MKI67,RRM2,SHCBP1,AURKB,GTSE1,ASPM,ASF1B,NCAPG,KIAA0101,TYMS 
# # 我在进行人类肝组织(liver)的单细胞转录组细胞注释，我有以下10个top marker 基因：(MKI67,RRM2,SHCBP1,AURKB,GTSE1,ASPM,ASF1B,NCAPG,KIAA0101,TYMS )
# 
# FeaturePlot(sc_integer_singleR,
#             features = "KLRB1",
#             pt.size = 0.05,
#             label.size = 3
# )
# 
# celltype2$cell[which(celltype2$ClusterID == 13)] <- ""


# # cluster14---------------------------------------------------------------------
# # 寻找簇14的差异表达基因
# markers_cluster14 <- FindMarkers(sc_integer_singleR,
#                                  min.pct = 0.3, # 基因在至少30%的细胞中表达,提高筛选的参数
#                                  ident.1 = "cluster14",
#                                  )
# 
# # 筛选和排序标记基因
# markers_cluster14_1 <- markers_cluster14 %>%
#   dplyr::filter(p_val_adj < 0.05) %>%          # 筛选显著基因
#   dplyr::mutate(pct = pct.1 - pct.2) %>%       # 计算表达比例差异
#   dplyr::arrange(desc(avg_log2FC), desc(pct))  # 按表达差异排序
# 
# # 可视化前10个标记基因
# png("output/C14.png")
# DotPlot(sc_integer_singleR, assay = "RNA",
#         features = rownames(markers_cluster14_1)[1:10],
#         group.by = "seurat_clusters") +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# dev.off()
# 
# # 输出前10个标记基因名字(用于Cellmarker2)
# cluster14_top_10_genes <- rownames(markers_cluster14_1)[1:10]
# cat(paste(cluster14_top_10_genes, collapse = ","), "\n")
# 
# # GOLGA8A,VEGFA,WWP2,KCNQ1OT1,COL11A1,CCNL2,PTH1R,RASD1,ACAN,RNMT 
# # 我在对我在进行人类骨肉瘤osteosarcoma组织的单细胞转录组细胞注释，我有以下10个top marker 基因：(GOLGA8A,VEGFA,WWP2,KCNQ1OT1,COL11A1,CCNL2,PTH1R,RASD1,ACAN,RNMT)
# 
# celltype2$cell[which(celltype2$ClusterID == 14)] <- "Unidentified"


# # cluster15---------------------------------------------------------------------
# # 寻找簇15的差异表达基因
# markers_cluster15 <- FindMarkers(sc_integer_singleR,
#                                  min.pct = 0.3, # 基因在至少30%的细胞中表达,提高筛选的参数
#                                  ident.1 = "cluster15",
#                                  )
# 
# # 筛选和排序标记基因
# markers_cluster15_1 <- markers_cluster15 %>%
#   dplyr::filter(p_val_adj < 0.05) %>%          # 筛选显著基因
#   dplyr::mutate(pct = pct.1 - pct.2) %>%       # 计算表达比例差异
#   dplyr::arrange(desc(avg_log2FC), desc(pct))  # 按表达差异排序
# 
# # 可视化前10个标记基因
# png("output/C15.png")
# DotPlot(sc_integer_singleR, assay = "RNA", 
#         features = rownames(markers_cluster15_1)[1:10],
#         group.by = "seurat_clusters") +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# dev.off()
# 
# # 输出前10个标记基因名字(用于Cellmarker2)
# cluster15_top_10_genes <- rownames(markers_cluster15_1)[1:10]
# cat(paste(cluster15_top_10_genes, collapse = ","), "\n")
# 
# # Mfap4,Col6a2,Gdf10,Fmod,Adamts2,Tbx20,Col6a1,Loxl1,Cpxm1,Col1a1 
# 
# # 我在进行小鼠肝组织的单细胞转录组细胞注释，以下是top10 marker：Mfap4,Col6a2,Gdf10,Fmod,Adamts2,Tbx20,Col6a1,Loxl1,Cpxm1,Col1a1  ，这是什么细胞类型？最好都要有相应的文献支持，准确率得高
# 
# celltype2$cell[which(celltype2$ClusterID == 15)] <- "HSC"



# # cluster16---------------------------------------------------------------------
# # 寻找簇16的差异表达基因
# markers_cluster16 <- FindMarkers(sc_integer_singleR,
#                                  min.pct = 0.3, # 基因在至少30%的细胞中表达,提高筛选的参数
#                                  ident.1 = "cluster16",
#                                  )
# 
# # 筛选和排序标记基因
# markers_cluster16_1 <- markers_cluster16 %>%
#   dplyr::filter(p_val_adj < 0.05) %>%          # 筛选显著基因
#   dplyr::mutate(pct = pct.1 - pct.2) %>%       # 计算表达比例差异
#   dplyr::arrange(desc(avg_log2FC), desc(pct))  # 按表达差异排序
# 
# # 可视化前10个标记基因
# png("output/C16.png")
# DotPlot(sc_integer_singleR, assay = "RNA",
#         features = rownames(markers_cluster16_1)[1:10],
#         group.by = "seurat_clusters") +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# dev.off()
# 
# # 输出前10个标记基因名字(用于Cellmarker2)
# cluster16_top_10_genes <- rownames(markers_cluster16_1)[1:10]
# cat(paste(cluster16_top_10_genes, collapse = ","), "\n")
# 
# # Dnah12,Tmem212,1500015O10Rik,Odf3b,1700007G11Rik,Cfap43,Spag16,Ak7,3300002A11Rik,Rarres2
# # 我在进行小鼠脊髓组织spinal的单细胞转录组细胞注释，我有以下10个top marker 基因：(Dnah12,Tmem212,1500015O10Rik,Odf3b,1700007G11Rik,Cfap43,Spag16,Ak7,3300002A11Rik,Rarres2 )
# 
# FeaturePlot(sc_integer_singleR,
#             features = "SiglecF",
#             pt.size = 0.05,
#             label.size = 3
# )
# 
# celltype2$cell[which(celltype2$ClusterID == 16)] <- "Ependymal cell"


# # cluster17---------------------------------------------------------------------
# # 寻找簇17的差异表达基因
# markers_cluster17 <- FindMarkers(sc_integer_singleR,
#                                  min.pct = 0.3, # 基因在至少30%的细胞中表达,提高筛选的参数
#                                  ident.1 = "cluster17")
# 
# # 筛选和排序标记基因
# markers_cluster17_1 <- markers_cluster17 %>%
#   dplyr::filter(p_val_adj < 0.05) %>%          # 筛选显著基因
#   dplyr::mutate(pct = pct.1 - pct.2) %>%       # 计算表达比例差异
#   dplyr::arrange(desc(avg_log2FC), desc(pct))  # 按表达差异排序
# 
# # 可视化前10个标记基因
# png("output/C17.png")
# DotPlot(sc_integer_singleR, assay = "RNA",
#         features = rownames(markers_cluster17_1)[1:10],
#         group.by = "seurat_clusters") +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# dev.off()
# 
# # 输出前10个标记基因名字(用于Cellmarker2)
# cluster17_top_10_genes <- rownames(markers_cluster17_1)[1:10]
# cat(paste(cluster17_top_10_genes, collapse = ","), "\n")
# 
# # RSAD2,OASL,IFIT2,IFIT3,IFIT1,ISG15,CMPK2,TNFSF10,HERC5,MX1
# # 我在进行人类骨肉瘤osteosarcoma组织的单细胞转录组细胞注释，我有以下10个top marker 基因：(RSAD2,OASL,IFIT2,IFIT3,IFIT1,ISG15,CMPK2,TNFSF10,HERC5,MX1)
# 
# FeaturePlot(sc_integer_singleR,
#             features = "RSAD2",
#             pt.size = 0.05,
#             label.size = 3
# )
# 
# celltype2$cell[which(celltype2$ClusterID == 17)] <- "Fibroblast"


# #  cluster18--------------------------------------------------------------------
# # 寻找簇18的差异表达基因
# markers_cluster18 <- FindMarkers(sc_integer_singleR,
#                                  min.pct = 0.3, # 基因在至少30%的细胞中表达,提高筛选的参数
#                                  ident.1 = "cluster18")
# 
# # 筛选和排序标记基因
# markers_cluster18_1 <- markers_cluster18 %>%
#   dplyr::filter(p_val_adj < 0.05) %>%          # 筛选显著基因
#   dplyr::mutate(pct = pct.1 - pct.2) %>%       # 计算表达比例差异
#   dplyr::arrange(desc(avg_log2FC), desc(pct))  # 按表达差异排序
# 
# # 可视化前10个标记基因
# png("output/C18.png")
# DotPlot(sc_integer_singleR, assay = "RNA",
#         features = rownames(markers_cluster18_1)[1:10],
#         group.by = "seurat_clusters") +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# dev.off()
# 
# # 输出前10个标记基因名字(用于Cellmarker2)
# cluster18_top_10_genes <- rownames(markers_cluster18_1)[1:10]
# cat(paste(cluster18_top_10_genes, collapse = ","), "\n")
# 
# # KRTAP13-1,ABCA12,COCH,EPHA8,TAC3,PLA2G4F,LINC01836,FABP7,FAM166B,TTC29 
# # 我在进行人类骨肉瘤osteosarcoma组织的单细胞转录组细胞注释，我有以下10个top marker 基因：(KRTAP13-1,ABCA12,COCH,EPHA8,TAC3,PLA2G4F,LINC01836,FABP7,FAM166B,TTC29 )
# 
# FeaturePlot(sc_integer_singleR,
#             features = c("CYP19A1"))
# 
# celltype2$cell[which(celltype2$ClusterID == 18)] <- "Unidentified"




# # cluster19---------------------------------------------------------------------
# # 寻找簇19的差异表达基因
# markers_cluster19 <- FindMarkers(sc_integer_singleR,
#                                  min.pct = 0.3, # 基因在至少30%的细胞中表达,提高筛选的参数
#                                  ident.1 = "cluster19")
# 
# # 筛选和排序标记基因
# markers_cluster19_1 <- markers_cluster19 %>%
#   dplyr::filter(p_val_adj < 0.05) %>%          # 筛选显著基因
#   dplyr::mutate(pct = pct.1 - pct.2) %>%       # 计算表达比例差异
#   dplyr::arrange(desc(avg_log2FC), desc(pct))  # 按表达差异排序
# 
# # 可视化前10个标记基因
# png("output/C19.png")
# DotPlot(sc_integer_singleR, assay = "RNA",
#         features = rownames(markers_cluster19_1)[1:10],
#         group.by = "seurat_clusters") +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# dev.off()
# 
# # 输出前10个标记基因名字(用于Cellmarker2)
# cluster19_top_10_genes <- rownames(markers_cluster19_1)[1:10]
# cat(paste(cluster19_top_10_genes, collapse = ","), "\n")
# 
# # Hhip,Csmd1,Cntfr,Cyp2e1,Cdh4,Lgr6,Enpp2,Aspn,Thbs2,Pde10a 
# # 我在进行小鼠肺组织(Lung)的单细胞转录组细胞注释，我有以下10个top marker 基因：(Hhip,Csmd1,Cntfr,Cyp2e1,Cdh4,Lgr6,Enpp2,Aspn,Thbs2,Pde10a )
# 
# FeaturePlot(sc_integer_singleR,
#             features = "Hhip",
#             pt.size = 0.05,
#             label.size = 3
# )
# 
# celltype2$cell[which(celltype2$ClusterID == 19)] <- "Myofibroblast"


# # cluster20---------------------------------------------------------------------
# # 寻找簇20的差异表达基因
# markers_cluster20 <- FindMarkers(sc_integer_singleR,
#                                  min.pct = 0.3, # 基因在至少30%的细胞中表达,提高筛选的参数
#                                  ident.1 = "cluster20")
# 
# # 筛选和排序标记基因
# markers_cluster20_1 <- markers_cluster20 %>%
#   dplyr::filter(p_val_adj < 0.05) %>%          # 筛选显著基因
#   dplyr::mutate(pct = pct.1 - pct.2) %>%       # 计算表达比例差异
#   dplyr::arrange(desc(avg_log2FC), desc(pct))  # 按表达差异排序
# 
# # 可视化前10个标记基因
# png("output/C20.png")
# DotPlot(sc_integer_singleR, assay = "RNA",
#         features = rownames(markers_cluster20_1)[1:10],
#         group.by = "seurat_clusters") +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# dev.off()
# 
# # 输出前10个标记基因名字(用于Cellmarker2)
# cluster20_top_10_genes <- rownames(markers_cluster20_1)[1:10]
# cat(paste(cluster20_top_10_genes, collapse = ","), "\n")
# 
# # LRP2,ITGB6,F5,DSC3,ITGB4,FBN2,GPR37,PRSS12,CMYA5,TENM3 
# # 我在进行人类蜕膜组织(decidua)的单细胞转录组细胞注释，我有以下10个top marker 基因：(LRP2,ITGB6,F5,DSC3,ITGB4,FBN2,GPR37,PRSS12,CMYA5,TENM3 )
# 
# FeaturePlot(sc_integer_singleR,
#             features = "HIGD1B",
#             pt.size = 0.05,
#             label.size = 3
# )
# 
# celltype2$cell[which(celltype2$ClusterID == 20)] <- "pDC"

# cluster21---------------------------------------------------------------------
# 寻找簇21的差异表达基因
markers_cluster21 <- FindMarkers(sc_integer_singleR,
                                 min.pct = 0.3, # 基因在至少30%的细胞中表达,提高筛选的参数
                                 ident.1 = "cluster21")

# 筛选和排序标记基因
markers_cluster21_1 <- markers_cluster21 %>%
  dplyr::filter(p_val_adj < 0.05) %>%          # 筛选显著基因
  dplyr::mutate(pct = pct.1 - pct.2) %>%       # 计算表达比例差异
  dplyr::arrange(desc(avg_log2FC), desc(pct))  # 按表达差异排序

# 可视化前10个标记基因
png("output/C21.png")
DotPlot(sc_integer_singleR, assay = "RNA",
        features = rownames(markers_cluster21_1)[1:10],
        group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

# 输出前10个标记基因名字(用于Cellmarker2)
cluster21_top_10_genes <- rownames(markers_cluster21_1)[1:10]
cat(paste(cluster21_top_10_genes, collapse = ","), "\n")

# CEP55,CDK1,HMMR,GTSE1,SHCBP1,RAD51AP1,MKI67,TOP2A,KIAA0101,TK1 
# 我在进行人类肝组织(Liver)的单细胞转录组细胞注释，我有以下10个top marker 基因：(CEP55,CDK1,HMMR,GTSE1,SHCBP1,RAD51AP1,MKI67,TOP2A,KIAA0101,TK1  )

FeaturePlot(sc_integer_singleR,
            features = "CD14",
            pt.size = 0.05,
            label.size = 3
)

celltype2$cell[which(celltype2$ClusterID == 21)] <- "pDC"


# # cluster23---------------------------------------------------------------------
# # 寻找簇23的差异表达基因
# markers_cluster23 <- FindMarkers(sc_integer_singleR,
#                                  min.pct = 0.3, # 基因在至少30%的细胞中表达,提高筛选的参数
#                                  ident.1 = "cluster23")
# 
# # 筛选和排序标记基因
# markers_cluster23_1 <- markers_cluster23 %>%
#   dplyr::filter(p_val_adj < 0.05) %>%          # 筛选显著基因
#   dplyr::mutate(pct = pct.1 - pct.2) %>%       # 计算表达比例差异
#   dplyr::arrange(desc(avg_log2FC), desc(pct))  # 按表达差异排序
# 
# # 可视化前10个标记基因
# png("output/C23.png")
# DotPlot(sc_integer_singleR, assay = "RNA",
#         features = rownames(markers_cluster23_1)[1:10],
#         group.by = "seurat_clusters") +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# dev.off()
# 
# # 输出前10个标记基因名字(用于Cellmarker2)
# cluster23_top_10_genes <- rownames(markers_cluster23_1)[1:10]
# cat(paste(cluster23_top_10_genes, collapse = ","), "\n")
# 
# # Adcy8,Clstn2,Cstdc2,Agtr1a,Higd1b,Trpc6,Fam162b,Vtn,Tafa1,Tmem200a 
# # 我在进行小鼠肺组织(Lung)的单细胞转录组细胞注释，我有以下10个top marker 基因：(Adcy8,Clstn2,Cstdc2,Agtr1a,Higd1b,Trpc6,Fam162b,Vtn,Tafa1,Tmem200a )
# 
# FeaturePlot(sc_integer_singleR,
#             features = "Ager",
#             pt.size = 0.05,
#             label.size = 3
# )
# 
# celltype2$cell[which(celltype2$ClusterID == 23)] <- "Pericyte"


# # cluster24---------------------------------------------------------------------
# # 寻找簇24的差异表达基因
# markers_cluster24 <- FindMarkers(sc_integer_singleR,
#                                  min.pct = 0.3, # 基因在至少30%的细胞中表达,提高筛选的参数
#                                  ident.1 = "cluster24")
# 
# # 筛选和排序标记基因
# markers_cluster24_1 <- markers_cluster24 %>%
#   dplyr::filter(p_val_adj < 0.05) %>%          # 筛选显著基因
#   dplyr::mutate(pct = pct.1 - pct.2) %>%       # 计算表达比例差异
#   dplyr::arrange(desc(avg_log2FC), desc(pct))  # 按表达差异排序
# 
# # 可视化前10个标记基因
# png("output/C24.png")
# DotPlot(sc_integer_singleR, assay = "RNA",
#         features = rownames(markers_cluster24_1)[1:10],
#         group.by = "seurat_clusters") +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# dev.off()
# 
# # 输出前10个标记基因名字(用于Cellmarker2)
# cluster24_top_10_genes <- rownames(markers_cluster24_1)[1:10]
# cat(paste(cluster24_top_10_genes, collapse = ","), "\n")
# 
# # NKG7,CCL5,CD2,VSIG4,FYB1,GZMA,FOLR2,GIMAP7,CD74,MS4A6A 
# # 我在进行人类骨肉瘤osteosarcoma组织的单细胞转录组细胞注释，我有以下10个top marker 基因：(NKG7,CCL5,CD2,VSIG4,FYB1,GZMA,FOLR2,GIMAP7,CD74,MS4A6A )
# 
# FeaturePlot(sc_integer_singleR,
#             features = "CD74",
#             pt.size = 0.05,
#             label.size = 3
# )
# 
# celltype2$cell[which(celltype2$ClusterID == 24)] <- "Pericyte"


# # cluster26---------------------------------------------------------------------
# # 寻找簇26的差异表达基因
# markers_cluster26 <- FindMarkers(sc_integer_singleR,
#                                  min.pct = 0.3, # 基因在至少30%的细胞中表达,提高筛选的参数
#                                  ident.1 = "cluster26")
# 
# # 筛选和排序标记基因
# markers_cluster26_1 <- markers_cluster26 %>%
#   dplyr::filter(p_val_adj < 0.05) %>%          # 筛选显著基因
#   dplyr::mutate(pct = pct.1 - pct.2) %>%       # 计算表达比例差异
#   dplyr::arrange(desc(avg_log2FC), desc(pct))  # 按表达差异排序
# 
# # 可视化前10个标记基因
# png("output/C26.png")
# DotPlot(sc_integer_singleR, assay = "RNA",
#         features = rownames(markers_cluster26_1)[1:10],
#         group.by = "seurat_clusters") +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# dev.off()
# 
# # 输出前10个标记基因名字(用于Cellmarker2)
# cluster26_top_10_genes <- rownames(markers_cluster26_1)[1:10]
# cat(paste(cluster26_top_10_genes, collapse = ","), "\n")
# 
# # NKG7,CCL5,CD2,VSIG4,FYB1,GZMA,FOLR2,GIMAP7,CD74,MS4A6A
# # 我在进行人类骨肉瘤osteosarcoma组织的单细胞转录组细胞注释，我有以下10个top marker 基因：(NKG7,CCL5,CD2,VSIG4,FYB1,GZMA,FOLR2,GIMAP7,CD74,MS4A6A )
# 
# FeaturePlot(sc_integer_singleR,
#             features = "CD74",
#             pt.size = 0.05,
#             label.size = 3
# )
# 
# celltype2$cell[which(celltype2$ClusterID == 26)] <- "Pericyte"



#------------------------------------------------------------------------------#
# ！！！！！！！！！！！恭喜你 最繁琐的细胞注释终于完成了！！！！！！！！！！！！
#------------------------------------------------------------------------------#



