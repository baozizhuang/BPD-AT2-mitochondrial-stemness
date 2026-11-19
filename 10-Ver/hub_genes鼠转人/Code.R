
# 1.加载R包---------------------------------------------------------------------
library(dplyr)
library(openxlsx)
library(rvest)
library(msig)
library(stringr)
library(readxl)
library(readr)
library(homologene)

# 2.创建文件夹------------------------------------------------------------------
dir.create("input/")
dir.create("output/")


# 3.将小鼠的基因转人------------------------------------------------------------

# 读取小鼠hub基因数据
mus_genes <- read_csv("input/hub_genes.csv")
# 获取基因列表
genelist <- mus_genes$shared_genes

#使用homologene函进行转换
#@genelist是要转换的基因列表
#@inTax是输入的基因列表所属的物种号，10090是小鼠
#@outTax是要转换成的物种号，9606是人

# 开始转换
genelistid <- homologene(genelist, inTax = 10090, outTax = 9606)
# 提取1,2列
genelistid2 <- genelistid[,c(1,2)]

colnames(genelistid2) <- c("ID","gene")
genelistid2 <- as.data.frame(genelistid2)

# 提取genelistid2的第二列，把genelistid2赋值为RGs
RGs <- genelistid2[, 2, drop = FALSE]



# 导出--------------------------------------------------------------------------
write_csv(RGs, "output/human_hub_genes.csv")




