# -----------------------------------------------------------------------------#
#                  Step 1: 从 GTF 文件获取基因长度（人/鼠通用）
#  输入: 参考基因组 GTF 文件（Ensembl 格式）
#  输出: output/{species}_gene_length.csv（gene_id, gene_name, length）
# -----------------------------------------------------------------------------#


# 0 配置（用户修改这里）=========================================================
# 指定 GTF 文件路径（支持人和小鼠）
gtf_file <- "input/Homo_sapiens.GRCh38.114.gtf"

# 指定物种标签（用于输出文件名，如 "Human" 或 "Mouse"）
species <- "Human"



# 1 安装/加载R包================================================================
required_pkgs <- c("GenomicFeatures", "rtracklayer")
for (pkg in required_pkgs) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE)
    library(pkg, character.only = TRUE)
  }
}

# 检查必要的基础包
if (!require("GenomicFeatures") || !require("rtracklayer")) {
  stop("GenomicFeatures 和 rtracklayer 包加载失败，请检查安装。")
}


# 2 从 GTF 计算基因长度=========================================================
cat("正在读取 GTF 文件...\n")
txdb <- makeTxDbFromGFF(gtf_file, format = "gtf")

cat("正在计算基因长度（外显子合并后总长度）...\n")
exons_gene <- exonsBy(txdb, by = "gene")
exons_gene_lens <- lapply(exons_gene, function(x) {
  sum(width(reduce(x)))
})
gene_length_df <- data.frame(
  gene_id   = names(exons_gene_lens),
  length    = unlist(exons_gene_lens),
  stringsAsFactors = FALSE
)
cat(paste0("基因长度计算完成，共 ", nrow(gene_length_df), " 个基因\n"))


# 3 获取 gene_id → gene_name 映射===============================================
cat("正在提取 gene_id → gene_name 映射...\n")
gtf_gr <- rtracklayer::import(gtf_file)
gtf_df <- as.data.frame(gtf_gr)

# 提取 gene_id 和 gene_name，去重
id_map <- unique(gtf_df[, c("gene_id", "gene_name")])
id_map <- id_map[id_map$gene_id != "" & !is.na(id_map$gene_id), ]
id_map <- id_map[!duplicated(id_map$gene_id), ]
rownames(id_map) <- id_map$gene_id

cat(paste0("基因名映射完成，共 ", nrow(id_map), " 个基因\n"))


# 4 合并基因长度 + 基因名=======================================================
# 取 gene_length 和 id_map 的交集
common_genes <- intersect(gene_length_df$gene_id, id_map$gene_id)
gene_length_df <- gene_length_df[match(common_genes, gene_length_df$gene_id), ]
id_map <- id_map[common_genes, ]

result <- data.frame(
  gene_id   = common_genes,
  gene_name = id_map$gene_name,
  length    = gene_length_df$length,
  stringsAsFactors = FALSE
)

# 去除 gene_name 为空的
result <- result[result$gene_name != "" & !is.na(result$gene_name), ]

cat(paste0("最终基因长度文件: ", nrow(result), " 个基因\n"))


# 5 输出结果====================================================================
dir.create("output", showWarnings = FALSE)
out_file <- paste0("output/", species, "_gene_length.csv")
write.csv(result, out_file, row.names = FALSE)


