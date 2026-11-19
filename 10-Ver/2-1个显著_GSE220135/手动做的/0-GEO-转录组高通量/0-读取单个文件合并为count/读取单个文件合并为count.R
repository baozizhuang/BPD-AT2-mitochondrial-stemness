# -----------------------------------------------------------------------------#
#  读取 GEO 单个样本 count 文件，合并为总 count 矩阵
#  输入: input/GSE220135_RAW/ 下所有 *_counts.txt.gz 文件
#  输出: output/GSE220135_raw_count_matrix.csv
#  说明: 每个文件两列（gene_name, count），按基因名合并
# -----------------------------------------------------------------------------#

# 0 环境变量准备================================================================
rm(list = ls())
gc()

# 1 配置（用户修改这里）=========================================================
input_dir  <- "input/GSE220135_RAW"             # 原始 count 文件目录
output_dir <- "output"
output_file <- file.path(output_dir, "GSE220135_raw_count_matrix.csv")


# 2 读取所有 count 文件并合并====================================================
# 获取所有 count 文件列表
file_list <- list.files(input_dir, pattern = "_counts\\.txt\\.gz$", full.names = TRUE)
cat("找到", length(file_list), "个 count 文件\n")

# 提取样本 ID（GSM号，即文件名第一个下划线前的部分）
get_sample_id <- function(fname) {
  sub("_.*", "", basename(fname))
}

# 读取第一个文件作为基础矩阵骨架
cat("正在读取:", basename(file_list[1]), "\n")
first <- read.table(file_list[1], header = FALSE, sep = "\t",
                    stringsAsFactors = FALSE, comment.char = "")
colnames(first) <- c("gene_name", get_sample_id(file_list[1]))
count_mat <- first
cat("初始矩阵维度:", nrow(count_mat), "基因 x 1 样本\n")

# 遍历其余文件，按基因名列合并
for (f in file_list[-1]) {
  sid <- get_sample_id(f)
  df <- read.table(f, header = FALSE, sep = "\t",
                   stringsAsFactors = FALSE, comment.char = "")
  colnames(df) <- c("gene_name", sid)
  
  # 按 gene_name 合并
  count_mat <- merge(count_mat, df, by = "gene_name", all = TRUE)
}

cat("合并后矩阵维度:", nrow(count_mat), "基因 x", ncol(count_mat) - 1, "样本\n")


# 3 数据清洗====================================================================
# 设置基因名为行名
rownames(count_mat) <- count_mat$gene_name
count_mat$gene_name <- NULL

# NA 替换为 0（某个基因在某样本中未被检测到）
count_mat[is.na(count_mat)] <- 0

# 确保为数值矩阵（保存基因名后再转换，因 as.data.frame 会重置行名）
gene_names <- rownames(count_mat)
count_mat <- as.data.frame(lapply(count_mat, as.numeric))
rownames(count_mat) <- gene_names

# 去除在所有样本中 count 都为 0 的基因（未表达基因）
count_mat <- count_mat[rowSums(count_mat) > 0, ]
cat("去除全零基因后:", nrow(count_mat), "基因 x", ncol(count_mat), "样本\n")


# 4 输出结果====================================================================
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(count_mat, output_file)

