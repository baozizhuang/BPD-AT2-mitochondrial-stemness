# -----------------------------------------------------------------------------#
#                        2.机器学习可解释性: SHAP
#                       适用H2O / GBM / XGBoost
# -----------------------------------------------------------------------------#

# SHAP解释通常是在 测试集 上进行的
# 在模型表现良好且没有过拟合或欠拟合的前提下，SHAP也可以在 训练集 或 整体数据集 上进行

# 一键清理
rm(list = ls())
gc()

# 安装包
# install.packages("devtools")
# devtools::install_github("ModelOriented/shapviz")
# devtools::install_github("futureverse/doFuture")   # shap值计算
# devtools::install_github("ModelOriented/kernelshap")   # shap值计算

# 加载包
library(dplyr)
library(mlr3)
library(mlr3learners)
library(kernelshap)
library(shapviz)
library(ggplot2)


# (手动更改)选择用于SHAP分析的模型 (从之前训练好的多个模型中，选一个来做 SHAP 解释)
best_model <- "XGBoost"   # 可选：H2O / GBM / XGBoost（与1-Build_ML_Model+CI.R中的模型名一致）哪个AUC值最高选哪个

# (手动更改)定义分组变量名(实验和对照组)
positive_class <- "High_stemness"  # 正类：实验组
negative_class <- "Low_stemness"   # 负类：对照组


# 输入文件路径
input_learner <- "output/learners.rds"             # 训练好的模型列表（含多个模型）
input_param <- "output/best_param_list.rds"        # 各模型的最佳超参数
input_train <- "output/model_data_train.rds"       # 训练集（用于训练SHAP解释模型）


# 1 读入数据====================================================================
# 训练集
data_t <- readRDS(input_train)

# 模型列表
learner_list <- readRDS(input_learner)
learner <- learner_list[[best_model]]$learner

# 最佳超参数
best_param_list <- readRDS(input_param)
best_params <- best_param_list[[best_model]] %>% 
  select(1:(which(names(.) == "classif.auc") - 1)) %>% 
  as.list()
learner$param_set$values <- modifyList(
  learner$param_set$values,  # 默认参数
  best_params                # 最优参数（覆盖同名默认值）
)


# 2 创建任务====================================================================
# 分类任务
task_t <- as_task_classif(data_t, target = "group", positive = positive_class) # positive填实验组 # 把训练数据包装成 mlr3 能识别的分类任务格式
task_t


# 3 SHAP 值计算（训练集）=======================================================
# 训练模型
learner$train(task_t)

# 提取特征矩阵（不含目标变量）
X_explain <- data_t %>% select(-group)

# 计算 SHAP 值
set.seed(123)
shap_values <- kernelshap(
  object = learner$model,     # 模型
  X = X_explain,              # 特征数据框
  bg_X = X_explain,           # 背景数据（可用训练集子集）
  predict_type = "response"   # 概率尺度（SHAP值=对概率的贡献，f(x)在0~1之间）
)

# 转换为 shapviz 对象
shap <- shapviz(shap_values, X = X_explain)

# S：SHAP 值矩阵。每一行对应一个样本，每一列对应一个特征。每个元素表示该特征对该样本的预测贡献。
# X：用于可视化的特征值数据框。包含与 X_pred 相同的列名。它通常包含解释数据集的特征值。
# baseline：基线值。表示在 SHAP 值尺度上的平均预测值。
# S_inter（可选）：SHAP 交互值（如果计算了交互值）。


# 4 SHAP 可视化=================================================================
# 4.1 特征重要性图(全局解释,条形图)---------------------------------------------
p_bar <- sv_importance(shap, 
                       kind = "bar",        # c("bar", "beeswarm", "both", "no")
                       max_display = 18L, #根据变量数修改
                       fill = "#ffd460", 
                       number_size = 2.5,
                       show_numbers = F)+
  theme_classic() +
  theme(axis.text = element_text(colour = "black", size = 11),
        axis.title = element_text(colour = "black", size = 13),
        plot.margin = unit(c(0.1, 0.01, 0.1, 0.01), "cm")     # 上右下左
  )
p_bar
dev.off()

if (!dir.exists("output")) {
  dir.create("output")
}
ggsave("output/1-shap_bar.pdf", plot = p_bar, width = 11, height = 10, units = "cm", dpi = 300)



# 4.2 汇总图(summary plot)/蜂群图（beeswarm plot）------------------------------
shap2 <- shap
shap2$S[shap2$S > 10] <- 11
p_bee <- sv_importance(shap2, 
                       kind = "bee",
                       max_display = 18L,    # 显示前多少个
                       size = 0.66,          # 点的大小
                       # viridis_args = list(begin = 0.2, end = 0.9, option = "plasma"),  # "viridis"、"magma"、"plasma"、"inferno" 或 "cividis"
                       color_bar_title = "Feature value")+
  scale_colour_gradient(low = "#6dd5fa", high = "#fe0c56")+
  # scale_colour_gradient(low = "#1c89e4", high = "#ffd460")+
  # scale_colour_gradient(low = "#5c715e", high = "#fccf4d")+
  theme_classic() +
  theme(axis.text = element_text(colour = "black", size = 8),
        axis.title = element_text(colour = "black", size = 10),
        plot.margin = unit(c(0.1, 0.01, 0.1, 0.01), "cm"),     # 上右下左
        legend.text = element_blank(),   # 隐藏图例文本
        legend.ticks = element_blank(),
        legend.title = element_text(angle = 90, vjust = 0.5, size = 8),
        legend.key.size = unit(1.5, "lines"),  # 调整图例键的高度
        legend.key.width = unit(0.6, "lines")   # 调整图例键的宽度
  )
p_bee
dev.off()

ggsave("output/2-shap_bee.pdf", plot = p_bee, width = 12, height = 8, units = "cm", dpi = 300)


# "viridis"：黄绿蓝渐变。
# "magma"：红黄渐变。
# "plasma"：黄紫渐变。
# "inferno"：黄黑渐变。
# "cividis"：从蓝色到黄色的渐变。


# 4.3 瀑布图(最常规的SHAP累加归因图，单样本解释)--------------------------------
p_wf <- sv_waterfall(shap,
                     fill_colors = c("#fe0c56","#6dd5fa"),
                     row_id = 1L,                  # 默认使用第一个样本
                     size = 8,
                     annotation_size = 4.5    # f(x)字体大小
)+
  theme_classic() +
  theme(axis.text = element_text(colour = "black", size = 8),
        axis.title = element_text(colour = "black", size = 10),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),     # 上右下左
        axis.ticks = element_blank(),
        axis.line.y = element_blank()  # 隐藏坐标轴的线
  )
p_wf 
dev.off()

ggsave("output/3-shap_waterfall.pdf", plot = p_wf, width = 15, height = 10, units = "cm", dpi = 300)


# 4.4 力图(瀑布图的一种变体，单样本解释)----------------------------------------
# 负类中 f(x) 最低的 10 个样本
tibble(
  row_id = 1:nrow(X_explain),
  sample = rownames(X_explain),
  outcome = data_t$group,
  f_x = shap$baseline + rowSums(shap$S)
) %>%
  filter(outcome == negative_class) %>%
  arrange(f_x) %>%
  head(10)

# 正类中 f(x) 最高的 10 个样本
tibble(
  row_id = 1:nrow(X_explain),
  sample = rownames(X_explain),
  outcome = data_t$group,
  f_x = shap$baseline + rowSums(shap$S)
) %>%
  filter(outcome == positive_class) %>%
  arrange(desc(f_x)) %>%
  head(10)


# 分别挑了一个负类和一个正类的样本，展示各自的特征贡献
# 红色：推高预测值（朝正类方向），蓝色：拉低预测值（朝负类方向）
# 条块长度代表特征影响大小

# 负类，对照样本，f(x) < baseline（蓝色占主导，预测为负类）
p_fc <- sv_force(shap,
                 fill_colors = c("#fe0c56","#6dd5fa"),
                 row_id = 3L,                  # 第3个样本
                 size = 8,                     # 文本大小
                 annotation_size = 4.5
)+
  theme_classic() +
  theme(axis.text = element_text(colour = "black", size = 8),
        axis.title = element_text(colour = "black", size = 10),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),     # 上右下左
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        axis.line.y = element_blank()  # 隐藏坐标轴的线
  )
p_fc
ggsave("output/4-shap_force-0.pdf", plot = p_fc, width = 16, height = 6, units = "cm", dpi = 300)


# 正类，实验样本，f(x) > baseline（红色占主导，预测为正类）
p_fc <- sv_force(shap,
                 fill_colors = c("#fe0c56","#6dd5fa"),
                 row_id = 59L,                  # 第6个样本
                 size = 8,                      # 文本大小
                 annotation_size = 4.5
)+
  theme_classic() +
  theme(axis.text = element_text(colour = "black", size = 8),
        axis.title = element_text(colour = "black", size = 10),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),     # 上右下左
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        axis.line.y = element_blank()  # 隐藏坐标轴的线
  )
p_fc
ggsave("output/5-shap_force-1.pdf", plot = p_fc, width = 16, height = 6, units = "cm", dpi = 300)


