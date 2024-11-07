
devtools::install("/home/pyk/stochtree")

library(stochtree)


# 加载数据集
data("mtcars")

# 设置随机种子
set.seed(123)

# 定义训练集比例
train_ratio <- 0.8

# 获取行数
n <- nrow(mtcars)

# 随机选择行索引，生成训练集和测试集
train_indices <- sample(1:n, size = train_ratio * n)

# 分割数据集
train_set <- mtcars[train_indices, ]
test_set <- mtcars[-train_indices, ]

# 将最后一列作为响应变量
X_train <- train_set[, -ncol(mtcars)]  # 训练集特征
train_y <- train_set[, ncol(mtcars)]   # 训练集响应变量

test_x <- test_set[, -ncol(mtcars)]    # 测试集特征
test_y <- test_set[, ncol(mtcars)]     # 测试集响应变量

# 查看分割结果
cat("训练集大小:", nrow(train_set), "\n")
cat("测试集大小:", nrow(test_set), "\n")

# 调用 R 函数
bart(X_train,train_y)

train_cov_preprocess_list <- preprocessTrainData(X_train)
X_train_metadata <- train_cov_preprocess_list$metadata
X_train <- train_cov_preprocess_list$data
original_var_indices <- X_train_metadata$original_var_indices
feature_types <- X_train_metadata$feature_types
feature_types <- as.integer(feature_types)
forest_dataset_train <- createForestDataset(X_train)
forest_model_mean <- createForestModel(forest_dataset_train, feature_types, 200, nrow(X_train), 0.95, 2, 5, 10)

forest_model_mean
bart_new(X_train,train_y)
bart_test(X_train,train_y)
