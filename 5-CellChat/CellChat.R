# 0.加载R包=====================================================================
# devtools::install_github("sqjin/CellChat")

library(CellChat)
library(Seurat)
library(magrittr)
library(ggplot2)


# 1.分析========================================================================
# 读取sc_integer.rds数据
sc_integer <- readRDS("input/AT2_CytoTRACE.rds")

# # 去除未知和增殖期的细胞
# sc_integer = subset(sc_integer, idents = c("Mesothelial cell", "APC", "M2 Macrophage", "Stressed mesothelial cells",
#                                                                        "T cell", "B cell", "NK", "Pericyte", "Mast cell"))

Idents(sc_integer) <- "CytoTRACE_group"

# 粗略查看
DimPlot(sc_integer)


# 导出normalized data和细胞类型注释信息
data.input <- GetAssayData(sc_integer[["RNA"]], layer = "data") # normalized data matrix
meta <- data.frame(row.names = rownames(sc_integer@meta.data),
                   CellType = sc_integer$CytoTRACE_group)   # 修改标签

# 创建一个Cell Chat对象
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "CellType")

# 设置cellType为默认的细胞分类方式
cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "CellType")
levels(cellchat@idents)

# 导入人类配体受体数据库
# CellChatDB <- CellChatDB.human
# 导入小鼠配体受体数据库
CellChatDB <- CellChatDB.mouse


#使用整个配体库来识别相互作用
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use

#提取细胞通讯相关基因
cellchat <- subsetData(cellchat)
future::plan("default")

#鉴定高表达的配体和受体
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

#储存到cellchat的cellchat@LR$LRsig中，将配体受体相关基因表达值投射到互作网络上并进行校正，时间较久
# cellchat <- projectData(cellchat,
#                         PPI.human)
cellchat <- projectData(cellchat,
                        PPI.mouse)

# 如果使用基因表达数据投射到蛋白质互作（PPI）网络上,后续要raw.use = FALSE,以使用投射到PPI上的数据

#------------------------------------------------------------------------------#
# （删除不要的细胞簇后的特别操作）
# 清理元数据和 CellChat 内部标识符的因子水平
# cellchat@meta$CellType <- droplevels(cellchat@meta$CellType)
# cellchat@idents <- droplevels(cellchat@idents)
#------------------------------------------------------------------------------#


#推断细胞通讯网络，为每个通讯网络计算一个概率值和对应的p值
cellchat <- computeCommunProb(cellchat, raw.use = TRUE) #raw.use表示使用未校正的数据
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)

#统计细胞通讯数量和强度（概率）
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))



# 2.可视化======================================================================
## 2.1第一部分：圆形网络图（权重）----------------------------------------------
pdf("output/2-cellchat_circle_weight.pdf", height = 6, width = 6)
netVisual_circle(cellchat@net$weight,
                 vertex.weight = groupSize,
                 weight.scale = T,
                 label.edge= F,
                 title.name = "Weight of interactions")
dev.off()


# 2.2第二部分：热图-------------------------------------------------------------
# pdf("output/1-cellchat_heatmap.pdf", height = 10/2.2, width = 12/2.2) # 大尺寸
pdf("output/1-cellchat_heatmap.pdf", height = 4, width = 4.5) # 小尺寸
netVisual_heatmap(cellchat, color.heatmap = "Reds")
dev.off()


# 2.3第三部分：气泡图-----------------------------------------------------------
#（未提取重要细胞通讯通路）（未细分细胞）的粗略查看
pdf("output/4-cellchat_buubble_all.pdf", height = 8, width = 6)
netVisual_bubble(cellchat,
                 remove.isolate = FALSE)
dev.off()

# 仅供参考
#-------------------------------------------------------------------------------#
# 查看发送接收数量
hub <- netVisual_bubble(cellchat, 
                        sources.use = NULL, 
                        targets.use = NULL, 
                        return.data = T,
                        remove.isolate = T)#隐藏无相互作用的项
hub
table(hub$communication$source)#信号发送方统计
table(hub$communication$target)#信号接收方统计
#-------------------------------------------------------------------------------#


# 提取重要细胞通讯通路~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
df.net <- subsetCommunication(cellchat)  # 获取所有互作数据框
top_interactions <- df.net %>%
  arrange(desc(prob)) %>%  # 按概率降序排列
  slice_head(n = 30)       # 取前30行


# 提取这些个互作涉及的通路名称
top_pathways <- unique(top_interactions$pathway_name)


# 可视化
# 筛选细胞类型进行分析

# 全部细胞（提取重要细胞通讯通路）气泡图----------------------------------------
pdf("output/3-CellChat_Bubble_all_cells.pdf", height = 35/5, width = 26/2.5)
netVisual_bubble(cellchat,
                 signaling = top_pathways,
                 remove.isolate = T ) #不展示孤立的细胞类型，都会在可视化中显示出来
dev.off()


# # 针对目标细胞（提取重要细胞通讯通路）气泡图------------------------------------
# # 统计每个 cell_subtype 的细胞数量
# table(Idents(sc_integer))
# 
# key_cells <- "Foam-like VSMC"   # 高代谢VSMC亚群
# 
# # 高代谢VSMC亚群做发送~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# pdf("output/4-CellChat_Bubble_immune_cells高代谢VSMC亚群做发送.pdf", height = 35/2, width = 26/3)
# cellchat@idents
# netVisual_bubble(cellchat,
#                  # signaling = top_pathways, # 针对Top通路
#                  sources.use =  key_cells,
#                  # targets.use = in_cells,
#                  remove.isolate = T ) #不展示孤立的细胞类型，都会在可视化中显示出来
# dev.off()
# 
# 
# pdf("output/5-CellChat_Bubble_immune_cells高代谢VSMC亚群做发送.pdf", height = 35/5, width = 26/4)
# cellchat@idents
# netVisual_bubble(cellchat,
#                  signaling = top_pathways, # 针对Top通路
#                  sources.use =  key_cells,
#                  # targets.use = in_cells,
#                  remove.isolate = T ) #不展示孤立的细胞类型，都会在可视化中显示出来
# dev.off()
# 
# 
# 
# # 高代谢VSMC亚群做接收~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# pdf("output/6-CellChat_Bubble_immune_cells高代谢VSMC亚群做接收.pdf", height = 35/3, width = 26/3)
# cellchat@idents
# netVisual_bubble(cellchat,
#                  # signaling = top_pathways, # 针对Top通路
#                  # sources.use =  ,
#                  targets.use = key_cells,
#                  remove.isolate = T ) #不展示孤立的细胞类型，都会在可视化中显示出来
# dev.off()
# 
# 
# pdf("output/7-CellChat_Bubble_immune_cells高代谢VSMC亚群做接收.pdf", height = 35/6, width = 26/4)
# cellchat@idents
# netVisual_bubble(cellchat,
#                  signaling = top_pathways, # 针对Top通路
#                  # sources.use =  ,
#                  targets.use = key_cells,
#                  remove.isolate = T ) #不展示孤立的细胞类型，都会在可视化中显示出来
# dev.off()
# 
# 
# # 指定细胞间通讯强度的网络图~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# pdf("output/8-CellChat_CircleWeight_高代谢VSMC亚群与巨噬细胞.pdf", height = 8, width = 8)
# netVisual_circle(cellchat@net$weight,
#                  vertex.weight = groupSize,
#                  sources.use = c("Foam-like VSMC","Macrophage"),
#                  targets.use = c("Foam-like VSMC","Macrophage"),
#                  weight.scale = T,
#                  label.edge= F,
#                  title.name = "Weight of interactions")
# dev.off()


# 细胞通讯信号角色热图~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# pdf("output/9-CellChat_Signel.pdf", height = 15/1.6, width = 13/2.3)
pdf("output/9-CellChat_Signel.pdf", height = 5, width = 5) # 小尺寸

cellchat <- netAnalysis_computeCentrality(cellchat,
                                          slot.name = "netP")
netAnalysis_signalingRole_heatmap(cellchat,
                                  width = 6,
                                  height = 8,
                                  pattern = "all") 
dev.off()




