
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


# 3定义变量---------------------------------------------------------------------
# 表型 线粒体能量代谢（mitochondrial energy metabolism） (MEMRGs)
pheno1 <- "mitochondrial energy metabolism"
pheno1_file1 <- paste0("input/GeneCards-SearchResults_",pheno1,".csv") #genecards里筛选的表型基因,注意输入的后缀
pheno1_file2 <- "input/mitochondrial energy metabolism-inPaperGenes.xlsx" #文献里筛选的表型基因,,列名改为gene


#定义表型基因的简写
pheno1_abbr <- "MEM"

# 3 表型合并、去重=============================================================
# 3.1 GeneCards来源------------------------------------------------------------
pheno1_RGs_genecards <- data.table::fread(file = pheno1_file1, data.table = F)

## 选择蛋白编码基因并依据相关性分数筛选基因（可修改）（选择所有蛋白编码基因）
pheno1_RGs_genecards_filter <- pheno1_RGs_genecards %>%
  dplyr::filter(Category == 'Protein Coding' & `Relevance score` > 1)
## 统一命名基因列
pheno1_RGs_genecards_filter <- pheno1_RGs_genecards_filter %>%
  dplyr::select("Gene Symbol") %>%
  dplyr::rename("gene" = "Gene Symbol")


## 3.2 文献来源的基因-----------------------------------------------------------
pheno1_RGs_paper <- read_xlsx("input/mitochondrial energy metabolism-inPaperGenes.xlsx",sheet = 1)   #文献里筛选的表型基因
colnames(pheno1_RGs_paper)[1] <- "gene"


# 表型基因的合并、去重（取并集）------------------------------------------------
pheno1_RGs <- union(pheno1_RGs_genecards_filter$gene, pheno1_RGs_paper$gene) %>%
  as.data.frame()
colnames(pheno1_RGs) <- "gene"



# 将筛选到的人的表型基因转成小鼠的基因再分析，这样就有小鼠的表型基因了----------

# 获取基因列表
genelist<- pheno1_RGs$gene

#使用homologene函进行转换
#@genelist是要转换的基因列表
#@inTax是输入的基因列表所属的物种号，10090是小鼠
#@outTax是要转换成的物种号，9606是人

# 开始转换
genelistid <- homologene(genelist, inTax = 9606, outTax = 10090)
# 提取1,2列
genelistid2 <- genelistid[,c(1,2)]

colnames(genelistid2) <- c("ID","gene")
genelistid2 <- as.data.frame(genelistid2)

# 提取genelistid2的第二列，把genelistid2赋值为RGs
RGs <- genelistid2[, 2, drop = FALSE]



# 导出表型1的表型基因(GeneCards与文献取并集)------------------------------------
write_csv(RGs, "output/MEMRGs.csv")




