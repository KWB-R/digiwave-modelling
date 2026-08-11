# ============================================================================
# 02_train_models.R — train the E. coli inlet-concentration model.
#   1) build the training frame from ModelInput_CLEAN_final.Rdata
#      (3h rolling-mean features joined to lab E. coli labels)
#   2) Quantile-Random-Forest variable selection + fit the deployed `models`
#   INPUT : data/raw/ModelInput_CLEAN_final.Rdata (config: MODEL_INPUT_RDATA)
#   OUTPUT: models/models.RData                 (list: full / sparse / bestR2)
#           data/derived/modelling_objects.rds  (train/test matrices for analyses)
#   Figure: plots/Per_vs_Vars.png
# Note: the stepwise/Lasso/NN model *comparison* now lives in
#       analyses/model_comparison.R (report figure, not the deployed model).
# ============================================================================

source("R/config.R")
source("R/functions.R")
library(tidyverse)
library(zoo)
library(quantregForest)
library(Metrics)

set.seed(42)

load(MODEL_INPUT_RDATA)   # provides df_KA

# ---- build training features (3h rolling means) ----------------------------
d <- df_KA %>% select(1:22, E_Coli)

window_size <- 12  # 12 x 15 min = 3 h
df_roll <- d %>% select(-E_Coli) %>%
  arrange(DateTime) %>%
  mutate(across(where(is.numeric),
                ~ rollapply(.x, width = window_size, FUN = mean, align = "right", fill = NA),
                .names = "{.col}_roll3h"))
df_roll3h_only <- df_roll %>% select(DateTime, contains("_roll3h"))

# ---- E. coli labels --------------------------------------------------------
ecoli    <- d %>% select(DateTime, E_Coli) %>% filter(!is.na(E_Coli))
ec_clean <- ecoli %>% filter(grepl("\\d{2}:\\d{2}:\\d{2}$", DateTime))

# ---- training frame --------------------------------------------------------
mod <- df_roll3h_only %>% inner_join(ec_clean)
mod$log_ec <- log10(mod$E_Coli)
mod$E_Coli <- NULL
mod <- mod %>% select(-contains("SAK"))

# ---- train/test split (data from 2024-04-01 onwards) -----------------------
mod <- mod %>% filter(DateTime > "2024-04-01")
n <- nrow(mod)
train_idx <- sample(seq_len(n), size = 0.7 * n)
train <- mod[train_idx, ]
test  <- mod[-train_idx, ]

predictors <- setdiff(names(train), c("DateTime", "log_ec"))
X_train <- as.matrix(train[, predictors]); y_train <- train$log_ec
X_test  <- as.matrix(test[, predictors]);  y_test  <- test$log_ec

# ============================================================================
# Quantile-Random-Forest variable selection + final models (from AccuracyVSnVars.R)
# ============================================================================
# 1. Fit initial model on full feature set
qrf_model <- quantregForest(
  x = X_train,
  y = y_train,
  ntree = 500,
  nodesize = 5
)

# 2. Get variable importance
var_imp <- importance(qrf_model)  # this returns %IncMSE-style importance
var_imp <- sort(var_imp[,1], decreasing = TRUE)

# We'll create an evaluation function that:
# - trains a model on a given subset of variables
# - predicts on test set
# - returns RMSE, MAE, R2

eval_subset <- function(vars) {
  # train with only 'vars'
  model <- quantregForest(
    x = X_train[, vars, drop = FALSE],
    y = y_train,
    ntree = 500,
    nodesize = 5
  )
  
  # predict median (50% quantile) on test set
  preds <- predict(model, newdata = X_test[, vars, drop = FALSE], what = 0.5)
  
  # metrics
  rmse_val <- rmse(y_test, preds)
  mae_val  <- mae(y_test, preds)
  
  # R² (coefficient of determination)
  ss_res <- sum((y_test - preds)^2)
  ss_tot <- sum((y_test - mean(y_test))^2)
  r2_val <- 1 - ss_res/ss_tot
  
  data.frame(
    n_vars = length(vars),
    rmse = rmse_val,
    mae  = mae_val,
    r2   = r2_val
  )
}

# 3. Iteratively drop the least important variable each round
#
# Approach:
# - Start with all variables, ordered by importance (most → least)
# - Keep top k variables for k = p, p-1, ..., 1
#
# This gives you a "how well do we do with the top k features?" curve.

all_vars_ordered <- names(var_imp)  # already sorted dec = most important first
results_list <- list()

for (k in seq_along(all_vars_ordered)) {
  # keep the top N variables
  vars_keep <- all_vars_ordered[1:(length(all_vars_ordered) - k + 1)]
  
  message("Fitting with ", length(vars_keep), " vars ...")
  res <- eval_subset(vars_keep)
  results_list[[k]] <- res
}

results_df <- bind_rows(results_list)

# 4. Plot performance vs number of variables
# we'll do 3 panels or 3 lines. I'll show both options.

## Option A: 3 separate y-axes using facets (clean)
results_long <- results_df %>%
  tidyr::pivot_longer(cols = c(rmse, mae, r2),
                      names_to = "metric",
                      values_to = "value")

results_long$metric <- factor(results_long$metric,
                              levels = c("rmse", "mae", "r2"),
                              labels = c("RMSE", "MAE", "R2"))


p <- ggplot(results_long, aes(x = n_vars, y = value)) +
  geom_line(lwd =1.3, col = "blue") +
  geom_point(size = 3, col = "blue") +
  scale_x_continuous(breaks = unique(results_df$n_vars)) +
  scale_y_continuous(breaks = seq(0,1, 0.05), limits = c(0,1)) +
  labs(
    x = "Anzahl Variablen im Modell",
    y = "Wert",
    title = "Modellperformance vs Anzahl Vorhersagevariablen"
  ) +
  facet_wrap(~ metric) +
  theme_grey(base_size = 18, base_line_size = 1)
p  
ggsave("plots/Per_vs_Vars.png", plot = p, 
       width = 16, 
       height = 9, dpi = "retina")




#all_vars_ordered <- names(var_imp)  # from importance(), high -> low importance
#results_list <- list()

#for (k in seq_along(all_vars_ordered)) {
#  vars_keep <- all_vars_ordered[1:(length(all_vars_ordered) - k + 1)]
#  res <- eval_subset(vars_keep)
#  results_list[[k]] <- res
#}

#results_df <- dplyr::bind_rows(results_list)



library(dplyr)

results_df2 <- results_df %>%
  arrange(n_vars) %>%
  mutate(
    r2_lag = dplyr::lag(r2),
    r2_gain = r2 - r2_lag
  )

improving_steps <- results_df2 %>%
  filter(!is.na(r2_gain), r2_gain > 0.1)
improving_steps


if (nrow(improving_steps) > 0) {
  vars_to_keep_idx <- max(improving_steps$n_vars)
  # Explanation:
  # If 3 vars gave big gain and 5 vars also gave big gain,
  # we'll keep the top 5, because that includes the top 3 anyway.
} else {
  # fallback: no single step had +0.1 R² improvement
  # we'll pick the number of vars that gave the best overall R²
  best_row <- results_df %>% arrange(desc(r2)) %>% slice(1)
  vars_to_keep_idx <- best_row$n_vars
}

vars_final <- all_vars_ordered[1:vars_to_keep_idx]
vars_final <- setdiff(vars_final, "KA_Inf_T_roll3h")  # drop inlet temp (not available at scoring)


final_model <- quantregForest(
  x = X_train[, vars_final, drop = FALSE],
  y = y_train,
  ntree = 500,
  nodesize = 5
)

# Evaluate on test set
final_preds <- predict(final_model,
                       newdata = X_test[, vars_final, drop = FALSE],
                       what = 0.5)

final_rmse <- Metrics::rmse(y_test, final_preds)
final_mae  <- Metrics::mae(y_test, final_preds)

ss_res <- sum((y_test - final_preds)^2)
ss_tot <- sum((y_test - mean(y_test))^2)
final_r2 <- 1 - ss_res/ss_tot

cat("Final vars used:\n")
print(vars_final)

cat("\nFinal performance:\n")
cat("RMSE:", final_rmse, "\n")
cat("MAE :", final_mae,  "\n")
cat("R2  :", final_r2,   "\n")


library(dplyr)

# make sure results_df is sorted so that row i = top i variables
results_df2 <- results_df %>%
  arrange(n_vars) %>%
  mutate(
    r2_prev = lag(r2),
    marginal_gain = r2 - r2_prev
  )

var_gains <- data.frame(
  var = all_vars_ordered,
  n_vars = seq_along(all_vars_ordered)  # 1,2,3,...
) %>%
  left_join(results_df2 %>%
              select(n_vars, marginal_gain, r2),
            by = "n_vars")


selected_vars <- var_gains %>%
  filter(
    n_vars == 1 |        # keep the #1 variable no matter what
      marginal_gain > 0.1  # keep this variable if it helped a lot when added
  ) %>%
  pull(var)

selected_vars 
selected_vars <- setdiff(selected_vars, "KA_Inf_T_roll3h")  # drop inlet temp (not available at scoring)

final_model_sparse <- quantregForest(
  x = X_train[, selected_vars, drop = FALSE],
  y = y_train,
  ntree = 500,
  nodesize = 5
)

final_preds <- predict(
  final_model_sparse,
  newdata = X_test[, selected_vars, drop = FALSE],
  what = 0.5
)

# metrics
library(Metrics)

final_rmse <- rmse(y_test, final_preds)
final_mae  <- mae(y_test, final_preds)

ss_res <- sum((y_test - final_preds)^2)
ss_tot <- sum((y_test - mean(y_test))^2)
final_r2 <- 1 - ss_res/ss_tot

cat("Variables kept:\n"); print(selected_vars)
cat("\nPerformance of sparse model:\n")
cat("RMSE:", final_rmse, "\n")
cat("MAE :", final_mae,  "\n")
cat("R2  :", final_r2,   "\n")


library(dplyr)

best_row <- results_df %>%
  arrange(desc(r2)) %>%
  slice(1)

best_n_vars <- best_row$n_vars
best_n_vars

vars_best_prefix <- all_vars_ordered[1:best_n_vars]

vars_best_prefix <- setdiff(vars_best_prefix, "KA_Inf_T_roll3h")  # drop inlet temp (not available at scoring)

library(quantregForest)
library(Metrics)

final_model_bestR2 <- quantregForest(
  x = X_train[, vars_best_prefix, drop = FALSE],
  y = y_train,
  ntree = 500,
  nodesize = 5
)

# predict median (50% quantile) on test set
best_preds <- predict(
  final_model_bestR2,
  newdata = X_test[, vars_best_prefix, drop = FALSE],
  what = 0.5
)

# evaluate
best_rmse <- rmse(y_test, best_preds)
best_mae  <- mae(y_test, best_preds)

ss_res <- sum((y_test - best_preds)^2)
ss_tot <- sum((y_test - mean(y_test))^2)
best_r2 <- 1 - ss_res/ss_tot

cat("Variables in best-R²-prefix model:\n")
print(vars_best_prefix)

cat("\nPerformance of best-R²-prefix model:\n")
cat("RMSE:", best_rmse, "\n")
cat("MAE :", best_mae,  "\n")
cat("R2  :", best_r2,   "\n")


# Named + ordered so downstream colours map: bestR2 (11 vars)=blue, best7 (6)=purple, best4of7 (3)=green
models <- list(
  bestR2   = final_model_bestR2,   # 11 vars — best cross-val R2 prefix
  best7    = final_model,          #  6 vars — top-7 importance prefix
  best4of7 = final_model_sparse    #  3 vars — high marginal-gain subset
)
save(models, file = "models/models.RData")

# ---- STANDALONE OUTPUT — train/test objects for analyses/model_comparison.R -
save_derived(list(
  mod = mod, train = train, test = test, predictors = predictors,
  X_train = X_train, y_train = y_train, X_test = X_test, y_test = y_test,
  df_roll3h_only = df_roll3h_only, ec_clean = ec_clean
), "modelling_objects")
