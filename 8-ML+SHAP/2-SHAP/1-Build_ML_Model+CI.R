# -----------------------------------------------------------------------------#
#                         1.机器学习建模，评估性能
#                 此流程仅用一个训练集来做SHAP,未融合验证集
#                          由于R包问题建议本地跑
#    脚本解释：
#    1-Build_ML_Model+CI.R  — 多模型（XGBoost/RF/SVM/GBM）超参调优 + 10折CV评估 + ROC，选出AUC最高模型保存
#    2-SHAP_Explanation1.R  — 对最佳模型（H2O / GBM / XGBoost）做SHAP解释
#    3-SHAP_Explanation2.R  — 对最佳模型（任意模型）做SHAP解释
# 
#                    流程：先跑1选最佳模型 → 再跑2或3
# -----------------------------------------------------------------------------#


# ---------------------Build_ML_Model+CI.R 详细说明----------------------------#
# 用5个Hub基因作为特征，分别训练 XGBoost/RF/SVM/GBM 多种模型，先各自调优超参数，
# 再用 10折交叉验证 全面评估性能（AUC + 置信区间 + ROC图），评估这5个基因作为特征的分类能力
# 最终目的是选出AUC最高的模型 ，保存结果给后续SHAP分析用。
# -----------------------------------------------------------------------------#


# 一键清理
rm(list = ls())
gc()

# 安装包
# if (!require(mlr3verse)) install.packages("mlr3verse")
# remotes::install_github("mlr-org/mlr3")
# if (!require(mlr3extralearners)) remotes::install_github("mlr-org/mlr3extralearners")
# remotes::install_github("mlr-org/mlr3learners")
# pak::pak("mlr-org/mlr3torch")
# install.packages("mlr3torch")
# install.packages("precrec")
# install.packages("gbm")
# install.packages("cvAUC")
# install.packages("xgboost")

# 加载包和函数
library(dplyr)
library(mlr3)
library(mlr3verse)
library(mlr3learners)
library(mlr3extralearners)
library(mlr3torch)
library(mlr3tuning)
library(paradox)
library(precrec)
library(ggplot2)
library(tibble)
library(cvAUC)
library(xgboost)
source("source/my_fun.R")


# 批量创建文件夹
# dir <- c("input","output","source")
# purrr::map(dir,dir.create)


# 1 读入数据====================================================================
# 1.1 读取训练集----------------------------------------------------------------
# 读取表达矩阵（行=基因，列=细胞）(所有未机器学习的基因矩阵)
expr <- read.csv("../1-ML-Pre/output/ML_expression_matrix.csv", row.names = 1, check.names = FALSE)

# 转置为 细胞 × 基因
expr_t <- as.data.frame(t(expr))
expr_t <- expr_t %>%
  tibble::rownames_to_column(var = "cell")

# 读取分组信息
labels <- read.csv("../1-ML-Pre/output/ML_group_labels.csv")

# 合并
data <- labels %>%
  left_join(expr_t, by = "cell") %>%
  tibble::column_to_rownames(var = "cell")

# 将 group 转为因子
data$group <- factor(data$group, levels = c("Low_stemness", "High_stemness"))


# 1.2 筛选Hub Genes（多种机器学习交集）-----------------------------------------
hub_genes <- read.csv("../2-ML/20260710_ML/Analysis/3-TRAIN/output/hub_genes.csv")
data <- data[, c("group", hub_genes$shared_genes)]


# 1.3 构建 Task 和 Learners ----------------------------------------------------
# 创建一个 mlr3 的 分类任务 （把数据转换成 mlr3 能理解的格式）
# 分类任务
task <- as_task_classif(data, target = "group", positive = "High_stemness")  # positive填实验组
task


# 2 超参数调整==================================================================

#-----------------------------------说明---------------------------------------#
# 优先调参顺序 
# 第一梯队（高敏）：XGBoost/LightGBM/GBM/RF/SVM
# 第二梯队（中敏）：AdaBoost/ET/ANN/DT
# 最后（通常不需要）：LR/KNN
#-----------------------------------说明---------------------------------------#


# 查看学习器的超参数,列出所有可用的模型
# param_set_all <- list(
#   XGBoost = lrn("classif.xgboost")$param_set,
#   LightGBM = lrn("classif.lightgbm")$param_set,
#   GBM = lrn("classif.gbm")$param_set,
#   RF = lrn("classif.randomForest")$param_set,
#   SVM = lrn("classif.svm")$param_set,
#   AdaBoost = lrn("classif.AdaBoostM1")$param_set,
#   ET = lrn("classif.ranger")$param_set,
#   DT = lrn("classif.reptree")$param_set,
#   LR = lrn("classif.log_reg")$param_set,
#   KNN = lrn("classif.kknn")$param_set
# )


# 定义学习器和对应的超参数搜索空间
learners <- list(
  # 树模型内核
  "XGBoost" = list(
    learner = lrn("classif.xgboost", predict_type = "prob"),
    param_set = ps(
      nrounds = p_int(lower = 100, upper = 500),         # 迭代次数
      eta = p_dbl(lower = 0.01, upper = 0.1),            # 学习率
      max_depth = p_int(lower = 3, upper = 10),          # 树的最大深度
      min_child_weight = p_dbl(lower = 1, upper = 10),
      subsample = p_dbl(lower = 0.5, upper = 1.0),
      gamma = p_dbl(lower = 0, upper = 5)
    )
  ),
  # # # 线性内核
  # # "XGBoost" = list(
  # #   learner = mlr3extralearners::lrn("classif.xgboost",
  # #                 predict_type = "prob",
  # #                 booster = "gblinear"),
  # #   param_set = ps(
  # #     nrounds = p_int(lower = 50, upper = 300), # 迭代次数（相当于训练轮数）
  # #     eta = p_dbl(lower = 0.01, upper = 0.3),   # 学习率
  # #     alpha = p_dbl(lower = 0, upper = 10),     # L1 正则化系数
  # #     lambda = p_dbl(lower = 0, upper = 10)     # L2 正则化系数
  # #     )
  # #   ),
  # # "LightGBM" = list(
  # #   learner = lrn("classif.lightgbm", predict_type = "prob"),
  # #   param_set = ps(
  # #     num_leaves = p_int(lower = 30, upper = 100),
  # #     max_depth = p_int(lower = 3, upper = 10),
  # #     learning_rate = p_dbl(lower = 0.01, upper = 0.1),
  # #     min_data_in_leaf = p_int(lower = 10, upper = 50),
  # #     feature_fraction = p_dbl(lower = 0.5, upper = 1.0)
  # #   )
  # # ),
  # # "GBM" = list(
  # #   learner = lrn("classif.gbm", predict_type = "prob"),
  # #   param_set = ps(
  # #     n.trees = p_int(lower = 100, upper = 500),           # 树的数量
  # #     interaction.depth = p_int(lower = 1, upper = 10),    # 树的深度
  # #     shrinkage = p_dbl(lower = 0.01, upper = 0.1),        # 学习率
  # #     n.minobsinnode = p_int(lower = 10, upper = 50)
  # #   )
  # # ),
  "RF" = list(
    learner = lrn("classif.randomForest", predict_type = "prob"),
    param_set = ps(
      ntree = p_int(lower = 300, upper = 500),        # 树的数量
      mtry = p_int(lower = 6, upper = 12),
      nodesize = p_int(lower = 2, upper = 5)
    )
  ),
  "SVM" = list(
    # 创建SVM分类器，设置预测类型为概率和分类类型
    learner = lrn("classif.svm",predict_type = "prob"),
    # 定义参数搜索空间
    param_set = ps(
      type = p_fct(levels = c("C-classification")),
      kernel = p_fct(
        levels = c("radial", "linear", "polynomial", "sigmoid")
      ),
      cost = p_dbl(lower = 0.1, upper = 10),
      gamma = p_dbl(lower = 0.01, upper = 1,
                    depends = quote(kernel %in% c("radial", "polynomial", "sigmoid"))),
      degree = p_int(lower = 2, upper = 5,
                     depends = quote(kernel == "polynomial")),  # 仅多项式核需要
      coef0 = p_dbl(lower = 0, upper = 1,        # 核函数独立项
                    depends = quote(kernel %in% c("polynomial", "sigmoid")))  # 仅某些核需要
    )
  ),
  "GBM" = list(
    learner = lrn("classif.gbm", predict_type = "prob"),
    param_set = ps(
      n.trees = p_int(lower = 100, upper = 500),           # 树的数量
      interaction.depth = p_int(lower = 1, upper = 10),    # 树的深度
      shrinkage = p_dbl(lower = 0.01, upper = 0.1),        # 学习率
      n.minobsinnode = p_int(lower = 1, upper = 10)        # 默认lower = 10, upper = 50
    )
  )
  # # "AdaBoost" = list(
  # #   learner = lrn("classif.AdaBoostM1", predict_type = "prob"),
  # #   param_set = ps(
  # #     I = p_int(lower = 50, upper = 200)         # 迭代次数
  # #   )
  # # ),
  # # "ET" = list(
  # #   learner = lrn("classif.ranger",
  # #                 predict_type = "prob",
  # #                 splitrule = "extratrees"),
  # #   param_set = ps(
  # #     num.trees = p_int(lower = 100, upper = 500),     # 树的数量
  # #     mtry = p_int(lower = 6, upper = 12),
  # #     min.node.size = p_int(lower = 1, upper = 10)
  # #   )
  # # ),
  # # "DT" = list(
  # #   learner = lrn("classif.reptree", predict_type = "prob"),
  # #   param_set = ps(
  # #     M = p_int(lower = 2, upper = 10),           # 树的最大深度
  # #     V = p_dbl(lower = 0.0001, upper = 0.01),
  # #     N = p_int(lower = 2, upper = 5)
  # #   )
  # # ),
  # # "LR" = list(
  # #   learner = lrn("classif.log_reg", predict_type = "prob"),
  # #   param_set = ps(
  # #     maxit = p_int(lower = 10, upper = 100),         # 最大迭代次数
  # #     epsilon = p_dbl(lower = 1e-10, upper = 1e-5)    # 收敛阈值
  # #   )
  # # ),
  # # "KNN" = list(
  # #   learner = lrn("classif.kknn", predict_type = "prob"),
  # #   param_set = ps(
  # #     k = p_int(lower = 1, upper = 20),
  # #     distance = p_dbl(lower = 1, upper = 5)
  # #   )
  # # ),
  # # 弹性网分类器
  # "ENet" = list(
  #   learner = lrn("classif.glmnet", predict_type = "prob"),
  #   param_set = ps(
  #     alpha = p_dbl(lower = 0.1, upper = 0.9),         # L1/L2 混合比例
  #     lambda = p_dbl(lower = 0.0001, upper = 1)    # 正则化强度
  #   )
  # )
)
# install.packages("RWeka")
# 创建组合终止器,告诉自动调优器"什么时候可以停了"。三个条件任意一个满足就停止
terminator <- trm("combo", list(
  trm("evals", n_evals = 50),                      # 最多评估50次
  trm("stagnation", threshold = 0.1, iters = 5),   # 基于性能提升
  trm("perf_reached", level = 0.95)                # 基于性能目标
))

# 创建自动调优器并执行调优
best_param_list <- list()
for (learner_name in names(learners)) {
  learner_info <- learners[[learner_name]]
  learner <- learner_info$learner
  param_set <- learner_info$param_set
  
  # 创建自动调优器
  at <- AutoTuner$new(
    learner = learner,
    resampling = rsmp("cv", folds = 5),  
    measure = msr("classif.auc"),         # 使用分类 AUC 作为评估指标
    search_space = param_set,             # 超参数搜索空间
    terminator = terminator,              # 终止条件
    tuner = tnr("grid_search", resolution = 5)  # 使用网格搜索
  )
  
  at$train(task)
  best_param_list[[learner_name]] <- at$archive$best() 
}

# 输出是一个 best_param_list ，里面存了每个模型的最佳参数组合，后续评估阶段直接用这个最优参数去训练和评估模型


# 3 使用最佳超参数重新训练模型==================================================
# 用最佳参数做10折交叉验证评估,产出每个机器学习模型的性能指标和ROC图

# 3.1 交叉验证用这个------------------------------------------------------------
# 评估结果
meas_rt <- data.frame(
  Model = character(),
  AUC = numeric(),
  Accuracy = numeric(),
  Precision = numeric(),
  Recall = numeric(),
  F1 = numeric(),
  AUC_CI_L   = numeric(),
  AUC_CI_U   = numeric(),
  stringsAsFactors = FALSE
)

# 交叉验证结果
comb_predict_list <- list()
# predict_list <- list()
comb_predic_rt_list <- list()
p_roc_list <- list()

# 对每个模型进行训练和评估
for (learner_name in names(learners)) {
  cat("正在评估模型：", learner_name, "\n")
  
  # 获取最佳超参数
  best_params <- best_param_list[[learner_name]] %>% 
    select(1:(which(names(.) == "classif.auc") - 1)) %>% 
    as.list()
  
  # 更新学习器的超参数
  learner <- learners[[learner_name]]$learner
  # learner$param_set$values <- best_params
  learner$param_set$values <- modifyList(
    learner$param_set$values,  # 默认参数
    best_params                # 最优参数（覆盖默认值）
  )
  
  # 设置重采样策略为交叉验证（小样本量用10折）
  resampling <- rsmp("cv", folds = 10)
  
  # 执行评估
  resampling_rt <- resample(task, learner, resampling)
  comb_predict_list[[learner_name]] <- resampling_rt
   
  # 设置评估指标
  measures <- msrs(c("classif.auc", "classif.acc", "classif.precision", "classif.recall"))
  fbeta_measure <- msr("classif.fbeta")
  fbeta_measure$param_set$values$beta <- 1
  measures <- c(measures, fbeta_measure)
  
  # 计算评估指标
  meas_scores <- resampling_rt$score(measures)
  
  
  # --- 交叉验证 AUC + bootstrap 95% CI ---
  preds_cv <- resampling_rt$prediction()
  truth_cv <- preds_cv$truth
  prob_cv  <- preds_cv$prob[, task$positive]
  
  # 构造 fold 向量（与样本一一对应）
  folds_list <- resampling_rt$resampling$instance
  nfolds <- length(folds_list)
  fold_vector <- rep(NA_integer_, nrow(data))
  for (k in seq_len(nfolds)) {
    fold_vector[folds_list[[k]]] <- k
  }
  
  set.seed(123)
  ci_obj <- cvAUC::ci.cvAUC(
    predictions  = prob_cv,
    labels       = as.character(truth_cv),   # 关键：字符型，顺序与 factor 一致
    label.ordering = c("Low_stemness", "High_stemness"),  # 负类在前，正类在后
    folds        = fold_vector
  )
  
  auc_ci_l <- round(ci_obj$ci[1], 4)
  auc_ci_u <- round(ci_obj$ci[2], 4)
  
  
  # 将结果添加到数据框中
  meas_rt <- rbind(meas_rt, data.frame(
    Model = learner_name,
    AUC = round(ci_obj$cvAUC, 4),
    Accuracy = mean(meas_scores$"classif.acc",na.rm = TRUE),
    Precision = mean(meas_scores$"classif.precision",na.rm = TRUE),
    Recall = mean(meas_scores$"classif.recall",na.rm = TRUE),
    F1 = mean(meas_scores$"classif.fbeta",na.rm = TRUE),
    AUC_CI_L  = auc_ci_l,
    AUC_CI_U  = auc_ci_u,
    stringsAsFactors = FALSE
  ))
  
  
  # 计算auc
  if (exists("auc")) {
    rm(auc)
  }
  auc <- paste0(round(ci_obj$cvAUC, 4),"(",
                auc_ci_l ,"-",auc_ci_u,")" 
                )
  
  # 提取预测结果
  comb_predic_rt_list[[learner_name]] <- resampling_rt$prediction()  # 汇总的所有折
  # predict_list[[learner_name]] <- resampling_rt$predictions()    # 单折预测结果
  
  # 绘制 ROC 曲线
  p_roc_list[[learner_name]] <- plot_roc(data = comb_predic_rt_list[[learner_name]], 
                                         learner = learner_name,
                                         curve_color = "#ff8c31",
                                         lw = 0.75,
                                         x_site = 0.65,
                                         y_site = 0.1,
                                         auc = auc )
}

# p_roc_list[[1]]

meas_rt <- meas_rt %>%
  mutate(across(everything(), ~ifelse(is.na(.), 0, .))) %>%  # 将NA替换为0
  mutate(across(-Model, ~round(., 4)))  # 保留4位小数（除了Model列）
write.csv(meas_rt,file = "output/train_measures_result.csv",row.names = F)
saveRDS(meas_rt,file = "output/train_measures_result.rds")


# 保存每个ROC曲线为单独的PDF
for (i in seq_along(names(p_roc_list))) {
  if (exists("name")) {
    rm(name)
  }
  name <- names(p_roc_list)[i]

  # 生成文件名，包含数字编号
  if (!dir.exists("output/roc_plots")) {
    dir.create("output/roc_plots")
  }
  filename <- paste0("output/roc_plots/",i,"-train_", name, "_roc.pdf")
  
  # 保存图像
  ggsave(
    filename = filename,
    plot = p_roc_list[[name]],
    device = "pdf",
    width = 9,   
    height = 9, 
    units = "cm", 
    dpi = 300  
  )
}


# 备份--------------------------------------------------------------------------
# 清理临时变量
envir_var <- c("best_params","learner_info","meas_scores","measures",
               "auc","auc1","learner_name","name") 
rm(list = envir_var)

# 保存给SHAP用的文件
save.image(file = "output/ml.RData")
saveRDS(data,file = "output/model_data_train.rds")
saveRDS(learners,file = "output/learners.rds")
saveRDS(best_param_list,file = "output/best_param_list.rds")




