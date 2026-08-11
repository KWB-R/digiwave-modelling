# ============================================================================
# R/functions.R — shared helpers for the DigiWave modelling pipeline.
# Sourced by every 0x_*.R stage (after config.R). Keeps the repeated logic
# (metrics, channel names, quantile prediction, 30-min binning, dose<->LRV)
# in one place so there is a single definition of each concept.
# Assumes dplyr / lubridate are loaded by the calling stage.
# ============================================================================

# ---- metrics ---------------------------------------------------------------
rmse <- function(actual, predicted) sqrt(mean((actual - predicted)^2))
rsq  <- function(actual, predicted) 1 - sum((actual - predicted)^2) /
                                        sum((actual - mean(actual))^2)

# ---- canonical sensor-channel names (post 3h rolling) ----------------------
# Order matches the plant CSV after selecting the *_roll3h columns.
CHANNELS_ROLL3H <- c(
  "DateTime",
  "KA_Eff_NO3_roll3h",  "KA_Eff_NH4_roll3h",  "KA_Eff_Turb_roll3h", "KA_Eff_Q_roll3h",
  "KA_Inf_LF_roll3h",   "KA_Inf_T_roll3h",    "KA_Inf_Q_roll3h",    "Bio_Eff_T_roll3h",
  "Bio_Eff1_O2_roll3h", "Bio_Eff2_O2_roll3h", "Bio_Eff1_NH4_roll3h","Bio_Eff2_NH4_roll3h",
  "Bio_Eff1_NO3_roll3h","Bio_Eff2_NO3_roll3h","Bio_Eff1_TR_roll3h", "Bio_Eff2_TR_roll3h",
  "GAK_Inf_Q_roll3h",   "GAK_Eff_SAK_roll3h", "SF_Eff_SAK_roll3h",  "SF_Inf_LF_roll3h",
  "SF_Inf_TURB_roll3h"
)
# Channels dropped before modelling (redundant / not predictive).
CHANNELS_DROP <- c("KA_Inf_T_roll3h", "SF_Eff_SAK_roll3h", "GAK_Eff_SAK_roll3h")

# ---- quantile prediction ---------------------------------------------------
# Predict a set of quantiles from a quantregForest model and return a tidy
# data.frame with columns q10, q50, ... bound to the newdata rows.
predict_quantiles <- function(model, newdata, percentiles = c(0.1, 0.5, 0.9)) {
  p <- predict(model, newdata = newdata, what = percentiles)
  p <- as.data.frame(p)
  colnames(p) <- paste0("q", percentiles * 100)
  dplyr::bind_cols(p, newdata)
}

# ---- 30-minute binning (mean of numeric columns) ---------------------------
bin_30min <- function(df, time_col = "DateTime") {
  df %>%
    dplyr::mutate(bin_30min = lubridate::floor_date(.data[[time_col]], "30 minutes")) %>%
    dplyr::group_by(bin_30min) %>%
    dplyr::summarize(dplyr::across(dplyr::where(is.numeric),
                                   ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
    dplyr::mutate(date = lubridate::date(bin_30min))
}

# ---- dose <-> LRV relationship ---------------------------------------------
# Empirical fit used across the reports: log(UV dose) = 4.1 + 1.7967 * log(LRV)
UV_A <- 4.1
UV_B <- 1.7967
# Required UV dose [J/m2] to achieve a target log-reduction value.
required_dose <- function(lrv) exp(UV_A + UV_B * log(lrv))
# Inverse: expected LRV for a given applied dose.
expected_lrv  <- function(dose) exp((log(dose) - UV_A) / UV_B)
