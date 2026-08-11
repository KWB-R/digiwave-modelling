# ============================================================================
# analyses/model_comparison.R — REPORT FIGURE, not part of the deployed pipeline.
#   Compares Stepwise LM / Lasso / Quantile RF / Neural Net for the inlet
#   E. coli model (Modelling.R) and their rolling-origin CV (Modelling2.R).
#   Self-contained: rebuilds the training frame from ModelInput and produces
#   plots/modellvergleich.png. The DEPLOYED model is fit in R/02_train_models.R.
# ============================================================================


source("R/config.R")
library(tidyverse)

load(MODEL_INPUT_RDATA)   # provides df_KA (see EXTERNAL_INPUTS.md)




df_KA %>% filter(is.na(E_Coli) == F)

vector <- which(colnames(df_KA) %>% 
                  str_detect(pattern = "linInterpoliert"))


d <- df_KA %>% select(1:22, E_Coli)

#d %>% gather(Parameter, Value, -DateTime) %>% 
#  ggplot(aes(x = DateTime, y = Value)) + 
#  geom_point()+ facet_wrap(.~Parameter)

library(zoo)

# Example: df has a datetime column and multiple numeric columns
# Assume your datetime column is called "datetime"

window_size <- 12  # 12 x 15min = 3h

df_roll <- d %>% select(-E_Coli) %>%
  arrange(DateTime) %>%
  mutate(across(
    where(is.numeric),
    ~ rollapply(.x, width = window_size, FUN = mean, align = "right", fill = NA),
    .names = "{.col}_roll3h"
  ))

df_roll3h_only <- df_roll %>%
  select(DateTime, contains("_roll3h"))



ecoli <- d %>% select(DateTime, E_Coli) %>% filter(is.na(E_Coli) == F)
ecoli %>% group_by(DateTime) %>% reframe(n = n()) %>% filter(n >1)

ec_clean <- ecoli %>%
  filter(grepl("\\d{2}:\\d{2}:\\d{2}$", DateTime))

mod <- df_roll3h_only %>% inner_join(ec_clean)

mod$log_ec <- log10(mod$E_Coli)
mod$E_Coli <- NULL



full_model <- lm(log_ec ~ . - DateTime, data = mod)
step_model_both <- step(full_model, direction = "both", trace = TRUE)
summary(step_model_both)

step_model_back <- step(full_model, direction = "backward", trace = TRUE)
summary(step_model_back)

step_model_for <- step(full_model, direction = "forward", trace = TRUE)
summary(step_model_for)


library(glmnet)

mod <- mod %>%
  select(-contains("SAK"))


# Select predictors (all except datetime and E_Coli)
predictors <- setdiff(names(mod), c("DateTime", "log_ec"))

X <- as.matrix(mod[, predictors])
y <- mod$log_ec


set.seed(123)  # for reproducibility

cv_lasso <- cv.glmnet(X, y, alpha = 1, nfolds = 10)  # alpha=1 → Lasso


plot(cv_lasso)   # shows CV error vs lambda

best_lambda <- cv_lasso$lambda.min   # lambda with lowest CV error
best_lambda_1se <- cv_lasso$lambda.1se  # simpler model within 1 SE rule
best_lambda
best_lambda_1se


coef(cv_lasso, s = "lambda.min")   # coefficients at best lambda
coef(cv_lasso, s = "lambda.1se")   # more regularized model



set.seed(42)

mod1 <- mod

mod <- mod %>% filter(DateTime > "2024-04-01")
n <- nrow(mod)

train_idx <- sample(seq_len(n), size = 0.7 * n)  # 70% training

train <- mod[train_idx, ]
test  <- mod[-train_idx, ]


#train_size <- floor(0.7 * n)
#train <- mod[1:train_size, ]
#test  <- mod[(train_size + 1):n, ]


# Build formula dynamically (exclude datetime and log_ec)
predictors <- setdiff(names(train), c("DateTime", "log_ec"))
formula <- as.formula(paste("log_ec ~", paste(predictors, collapse = " + ")))

# Full model
full_model <- lm(formula, data = train)

# Stepwise selection
step_model <- step(full_model, direction = "both", trace = 0)

# Predictions on test set
pred_step <- predict(step_model, newdata = test)


library(glmnet)

# Prepare matrices
X_train <- as.matrix(train[, predictors])
y_train <- train$log_ec
X_test  <- as.matrix(test[, predictors])
y_test  <- test$log_ec

# Fit lasso with CV
set.seed(42)
cv_lasso <- cv.glmnet(X_train, y_train, alpha = 1, nfolds = 10)

# Predict on test
pred_lasso <- predict(cv_lasso, s = "lambda.min", newx = X_test)



# Define metrics
rmse <- function(actual, predicted) {
  sqrt(mean((actual - predicted)^2))
}

rsq <- function(actual, predicted) {
  1 - sum((actual - predicted)^2) / sum((actual - mean(actual))^2)
}

# Calculate
rmse_step  <- rmse(y_test, pred_step)
rmse_lasso <- rmse(y_test, pred_lasso)

rsq_step  <- rsq(y_test, pred_step)
rsq_lasso <- rsq(y_test, pred_lasso)

cat("Stepwise: RMSE =", rmse_step, " R² =", rsq_step, "\n")
cat("Lasso:    RMSE =", rmse_lasso, " R² =", rsq_lasso, "\n")


library(ggplot2)

# Combine results into one data frame
results <- data.frame(
  actual = y_test,
  stepwise = pred_step,
  lasso = as.numeric(pred_lasso)
)

# Stepwise plot
p1 <- ggplot(results, aes(x = actual, y = stepwise)) +
  geom_point(color = "steelblue", alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Stepwise Regression",
       x = "Actual log_ec",
       y = "Predicted log_ec") +
  theme_minimal()

# Lasso plot
p2 <- ggplot(results, aes(x = actual, y = lasso)) +
  geom_point(color = "darkgreen", alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Lasso Regression",
       x = "Actual log_ec",
       y = "Predicted log_ec") +
  theme_minimal()

p1
p2

# Summary of the model
summary(step_model)

# Standardized coefficients (beta weights)
library(lm.beta)
lm.beta(step_model)
coef(cv_lasso, s = "lambda.min")

library(caret)

# For stepwise model
varImp(step_model)

# Extract coefficients as sparse matrix
coef_lasso <- coef(cv_lasso, s = "lambda.min")

# Convert to data frame
coef_df <- data.frame(
  term = rownames(coef_lasso),
  estimate = as.numeric(coef_lasso)
)

# Drop intercept
coef_df <- coef_df[coef_df$term != "(Intercept)", ]

# Rank by absolute size
coef_df <- coef_df %>%
  dplyr::mutate(abs_estimate = abs(estimate)) %>%
  dplyr::arrange(desc(abs_estimate))

head(coef_df, 10)   # top 10 predictors

head(coef_df, 10)   # top 10 predictors

library(ggplot2)

ggplot(coef_df %>% head(15), aes(x = reorder(term, abs_estimate), y = abs_estimate)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(title = "Lasso Variable Importance",
       x = "Predictor",
       y = "Absolute Coefficient") +
  theme_minimal()


library(quantregForest)


# Fit Quantile Random Forest
set.seed(42)
qrf_model <- quantregForest(X_train, y_train, ntree = 500, nodesize = 5)

# Predict median (0.5 quantile)
pred_qrf <- predict(qrf_model, newdata = X_test, what = 0.5)

# You can also get prediction intervals, e.g. 10% and 90% quantiles
pred_qrf_interval <- predict(qrf_model, newdata = X_test, what = c(0.1, 0.9))

newdata <- na.omit(df_roll3h_only)
pred <- predict(qrf_model, newdata = newdata, what = c(0.1, 0.5, 0.9))

pred_df <- as.data.frame(pred)

colnames(pred_df) <- c("q10", "q50", "q90" )

newdata <- newdata %>% left_join(ec_clean)

pred_df <- pred_df %>% bind_cols(newdata)


pred_df$Date <- lubridate::date(pred_df$DateTime)


daily <- pred_df %>% group_by(Date) %>% reframe(m = mean(q50),
                                                l = mean(q10),
                                                u = mean(q90),
                                                ec = mean(E_Coli,na.rm = T ))


ggplot(daily, aes(Date, y = m, ymin = l, ymax = u))+ 
  #geom_linerange(col = "steelblue") +
  geom_line(col = "blue", lwd = 1.3)+ 
  geom_ribbon(alpha = 0.3, fill = "blue")+ theme_gray(base_size = 16) +
  scale_x_date(date_breaks = "2 months")+ geom_point(mapping = aes(Date, log10(ec)), 
                                                     bg = "white", 
                                                     pch = 22,
                                                     col = "black", 
                                                     size = 3, 
                                                     alpha = 0.9)+
  ggtitle("Zeitreihe E.coli Vorhersage vs. Messwerte") + labs(x = "Datum", y = "E.coli [log10]")
  

library(nnet)

# Scale predictors (important for neural nets)
X_train_scaled <- scale(X_train)
X_test_scaled  <- scale(X_test, center = attr(X_train_scaled, "scaled:center"),
                        scale  = attr(X_train_scaled, "scaled:scale"))

# Fit a small neural network
set.seed(42)
nn_model <- nnet(X_train_scaled, y_train, size = 5, linout = TRUE, maxit = 500, decay = 0.01)

# Predict
pred_nn <- predict(nn_model, X_test_scaled)




rmse_qrf <- rmse(y_test, pred_qrf)
rsq_qrf  <- rsq(y_test, pred_qrf)

rmse_nn <- rmse(y_test, pred_nn)
rsq_nn  <- rsq(y_test, pred_nn)

cat("QRF: RMSE =", rmse_qrf, " R² =", rsq_qrf, "\n")
cat("NN:  RMSE =", rmse_nn, " R² =", rsq_nn, "\n")


results_qrf_nn <- data.frame(
  actual = y_test,
  qrf = pred_qrf,
  nn  = pred_nn
)


results_all <- results %>% select(-actual)%>% bind_cols(results_qrf_nn)
res_long <- results_all %>% gather(model, value, -actual)
res_long$model <- factor(res_long$model, ordered = T,
                         levels = c("stepwise", "lasso", "qrf", "nn"),
                         labels = c("Stepwise LR", "Lasso LR", "QRF", "Neuronales Netz"))

modelvergleich <- ggplot(res_long, aes(x = actual, y = value, col = model)) +
  facet_wrap(.~model) +
  geom_point( alpha = 0.6, size =4) +
  geom_abline(slope = 1, intercept = 0, lwd = 1.2, linetype = "dashed", color = "black") +
  theme_minimal(base_size = 20)+ 
  labs(x="gemessen E.coli [lg] MPN/100mL",
       y="vorhergesagt E.coli [lg] MPN/100mL")+
  scale_x_continuous(limits = c(0,6), breaks = seq(0, 6, 0.5))+ 
  scale_y_continuous(limits = c(0,6), breaks = seq(0, 6, 0.5))+
  theme(legend.position = "bottom") + 
  ggtitle("Vorhersagen gegen Testdaten")


ggsave(filename = "plots/modellvergleich.png", plot = modelvergleich, width = 18, height = 9, dpi = "retina")

ggplot(results_qrf_nn, aes(x = actual, y = qrf)) +
  geom_point(color = "purple", alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Quantile Random Forest (Median)",
       x = "Actual log_ec", y = "Predicted log_ec") +
  theme_minimal() + 
  scale_x_continuous(limits = c(0,6))+ 
  scale_y_continuous(limits = c(0,6))



ggplot(results_qrf_nn, aes(x = actual, y = nn)) +
  geom_point(color = "orange", alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Neural Network",
       x = "Actual log_ec", y = "Predicted log_ec") +
  theme_minimal()+ 
  scale_x_continuous(limits = c(0,6))+ 
  scale_y_continuous(limits = c(0,6))



# ============================================================================
# STANDALONE OUTPUT — persist objects consumed by other scripts
#   AccuracyVSnVars.R needs: X_train, y_train, X_test, y_test, df_roll3h_only, ec_clean
#   Modelling2.R      needs: mod
# ============================================================================
save_derived(list(
  mod             = mod,
  train           = train,
  test            = test,
  predictors      = predictors,
  X_train         = X_train,
  y_train         = y_train,
  X_test          = X_test,
  y_test          = y_test,
  df_roll3h_only  = df_roll3h_only,
  ec_clean        = ec_clean
), "modelling_objects")

# ============================ rolling-origin CV (was Modelling2.R) ===========

source("R/config.R")
library(caret)
library(glmnet)
library(quantregForest)
library(nnet)
library(doParallel)
library(ggplot2)
library(dplyr)
library(tidyr)

# STANDALONE INPUT — needs `mod` from Modelling.R
if (!exists("mod")) mod <- load_derived("modelling_objects")$mod

set.seed(42)

#------------------------------
# 1. Define predictors
#------------------------------
predictors <- setdiff(names(mod), c("DateTime", "log_ec"))

X <- mod[, predictors]
y <- mod$log_ec

#------------------------------
# 2. Define rolling CV scheme
#------------------------------
n <- nrow(mod)
initial_window <- floor(0.7 * n)
horizon <- floor(0.1 * n)

train_control <- trainControl(
  method = "timeslice",
  initialWindow = initial_window,
  horizon = horizon,
  fixedWindow = TRUE,
  savePredictions = "final"
)

#------------------------------
# 3. Set up parallel backend
#------------------------------
cl <- makeCluster(parallel::detectCores() - 1)  # leave one core free
registerDoParallel(cl)

#------------------------------
# 4. Fit models in parallel
#------------------------------
lm_model <- train(x = X, y = y, method = "lm", trControl = train_control)

lasso_model <- train(
  x = as.matrix(X), y = y,
  method = "glmnet",
  trControl = train_control,
  tuneLength = 10
)

qrf_model <- train(
  x = X, y = y,
  method = "qrf",
  trControl = train_control,
  tuneLength = 5,
  ntree = 300
)

nn_model <- train(
  x = scale(X), y = y,
  method = "nnet",
  trControl = train_control,
  tuneLength = 5,
  linout = TRUE,
  trace = FALSE,
  maxit = 300
)

#------------------------------
# 5. Stop parallel backend
#------------------------------
stopCluster(cl)
registerDoSEQ()

#------------------------------
# 6. Compare results
#------------------------------
results <- resamples(list(
  Stepwise = lm_model,
  Lasso    = lasso_model,
  QRF      = qrf_model,
  NeuralNet= nn_model
))

summary(results)
bwplot(results, metric = "RMSE")
dotplot(results, metric = "Rsquared")

#------------------------------
# 7. Extract predictions for plotting
#------------------------------
preds_all <- bind_rows(
  lm_model$pred    %>% mutate(Model = "Stepwise", pred = as.numeric(pred)),
  lasso_model$pred %>% mutate(Model = "Lasso",    pred = as.numeric(pred)),
  qrf_model$pred   %>% mutate(Model = "QRF",      pred = as.numeric(pred)),
  nn_model$pred    %>% mutate(Model = "NeuralNet",pred = as.numeric(pred))
)


# Align with DateTime (rownames in pred = row indices in mod)
preds_all <- preds_all %>%
  mutate(RowIndex = as.integer(rowIndex)) %>%
  mutate(DateTime = mod$DateTime[RowIndex])

#------------------------------
# 8. Plot: actual vs predicted over time
#------------------------------
ggplot(preds_all, aes(x = DateTime, y = pred)) +
  geom_line(color = "steelblue", alpha = 0.7) +
  geom_point(aes(y = obs), color = "black", size = 1) +
  facet_wrap(~ Model, scales = "free_y") +
  labs(title = "Predicted vs Actual log_ec (Rolling CV)",
       x = "DateTime", y = "log_ec") +
  theme_minimal()

ggplot(preds_all, aes(y = pred, x = obs, col = Model))+ 
  geom_point(alpha = 0.4) + 
  facet_wrap(.~Model) + geom_abline(slope = 1, intercept = 0, col = "red")
