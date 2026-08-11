# ============================================================================
# 05_predict_validate.R — score plant data with the models and validate
#   predictions against lab (uv_zulauf) + AquaBio measurements, overlaying
#   the applied UV-dose schedule. Canonical validator (replaces the duplicate
#   validation blocks formerly in ValidationPlots.R, CleanNewKAData.R, Database.R).
#   INPUTS (auto-loaded): models, stanDR, v_roll3h_only, dose_schedule
#   Figures: plots/model_*_3panel.png, plots/all_models_3panel.pdf
# ============================================================================
################################################################################
# 3-panel validation plots for 3 models
# Top:    time series (q50 + q10/q90 ribbon + lab points) + daily RMSE/R²
# Middle: applied UV dose schedule (from interval table `df`)
# Bottom: UV validation (ribbon + AquaBio points)
#
# Assumptions (must exist in your workspace OR you load/create them below):
#   - models       : list with at least 3 fitted models that support predict(..., what=percentiles)
#   - v_roll3h_only: data.frame used as newdata for prediction (has DateTime column)
#   - daily        : data.frame with columns Datum (Date) and dose (numeric) used with stanDR
#   - stanDR       : rstanarm model used for posterior_epred()
#
# Files used:
#   - ./uv_zulauf.csv (lab measurements; must have columns date and ecoli)
#   - Historical_ecoli.csv (AquaBio; first col datetime, second col ecoli)
################################################################################

# --------------------------- Libraries ---------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(data.table)
  library(ggplot2)
  library(gridExtra)
  library(quantregForest)  # for predict.quantregForest (models are QRF)
  library(rstanarm)   # for posterior_epred()
})

# --------------------------- User settings -----------------------------------
source("R/config.R")
TZ_USE <- "Europe/Berlin"

# Output directory for the generated figures
OUT_DIR <- PLOTS_DIR

# --------------------------- STANDALONE INPUT --------------------------------
# Objects produced by upstream scripts (auto-loaded from disk if not in session)
if (!exists("models"))        load(model_path("models.RData"))                # 02_train_models.R
if (!exists("v_roll3h_only")) v_roll3h_only <- load_derived("v_roll3h_only")  # 01_clean_plant_data.R
if (!exists("stanDR"))        stanDR <- readRDS(model_path("stanDR.rds"))     # 03_dose_response.R
# `daily` is overloaded (predictions vs dose schedule) — need the DOSE SCHEDULE
if (!exists("daily") || !all(c("Datum", "dose") %in% names(daily)))
  daily <- load_derived("dose_schedule")$daily                               # 04_dose_schedule.R (needs Datum + dose)

# --------------------------- UV dose interval table --------------------------
txt <- "
31.12.2025\t4\t54,5\t432,6\t80
01.12.2025\t4\t54,5\t432,6\t80
20.11.2025\t3,2\t61,8\t613,9\t81,6
30.10.2025\t4,6\t74,8\t526,4\t83
27.10.2025\t4,6\t71,9\t505,8\t82,2
24.09.2025\t3,5\t65,3\t604,4\t78,6
11.09.2025\t2,5\t58,4\t756,8\t77,9
05.09.2025\t3,9\t57,8\t479,6\t
20.08.2025\t4,8\t58,6\t395,4\t78,5
12.08.2025\t4\t67,2\t665,5\t83,2
"

df_intervals <- readr::read_tsv(
  I(txt),
  col_names = c("Datum", "m3_h", "W_m2", "J_m2", "pct"),
  col_types = cols(.default = col_character()),
  locale = locale(decimal_mark = ",")
) %>%
  mutate(
    Datum = dmy(Datum),
    across(-Datum, ~ suppressWarnings(parse_number(.x, locale = locale(decimal_mark = ","))))
  ) %>%
  arrange(Datum)

# --------------------------- Helper: expand interval table to daily series ----
make_dose_ts <- function(start_dates_df, date_from, date_to, dose_col = "J_m2") {
  day_seq <- seq.Date(as.Date(date_from), as.Date(date_to), by = "day")
  
  idx <- findInterval(day_seq, start_dates_df$Datum)
  idx[idx == 0] <- 1
  
  tibble(
    day  = day_seq,
    dose = start_dates_df[[dose_col]][idx]
  )
}

# --------------------------- Data prep: newdata -------------------------------
nd <- na.omit(v_roll3h_only)

# Quantiles to predict
percentiles <- c(0.1, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1)
q_names <- paste0("q", percentiles * 100)

# --------------------------- Read lab measurements ----------------------------
# Expect columns: date, ecoli
ec <- fread("data/uv_zulauf.csv")
if (!all(c("date", "ecoli") %in% names(ec))) {
  stop("uv_zulauf.csv must contain columns `date` and `ecoli`.")
}
ec$date <- dmy(ec$date)
ec$date <- as.POSIXct(ec$date, tz = TZ_USE)

# --------------------------- Read AquaBio ------------------------------------
aquabio <- fread(AQUABIO_ECOLI_HIST)
ablauf_uv <- aquabio %>% select(1, 2)
colnames(ablauf_uv) <- c("date", "ecoli")

# Try to coerce to POSIXct safely
# If your first column is already POSIXct, this will keep it; otherwise it will try parsing
ablauf_uv$date <- as.POSIXct(ablauf_uv$date, tz = TZ_USE)
if (all(is.na(ablauf_uv$date))) {
  # fallback parse common formats (edit if needed)
  ablauf_uv$date <- ymd_hms(aquabio[[1]], tz = TZ_USE, quiet = TRUE)
}
ablauf_uv <- ablauf_uv %>% filter(date > as.POSIXct("2025-08-19 10:00", tz = TZ_USE))

# --------------------------- Model keys + colors ------------------------------
# Use first 3 models in list if you don't have exactly model1/model2/model3 names
if (is.null(names(models)) || any(names(models)[1:3] == "")) {
  model_keys <- as.character(seq_len(min(3, length(models))))
  names(models) <- names(models) %||% model_keys
} else {
  model_keys <- names(models)[seq_len(min(3, length(models)))]
}

if (length(model_keys) < 3) stop("`models` must contain at least 3 models.")

model_cols <- c("blue", "purple", "forestgreen")  # model 1/2/3

# --------------------------- Validation window (zoom) ------------------------
# Focus every panel on the span where validation data actually exist (lab +
# AquaBio) so the three stacked panels are directly comparable on time.
.val_dates <- as.Date(c(ec$date, ablauf_uv$date))
val_min <- min(.val_dates, na.rm = TRUE) - 3
val_max <- max(.val_dates, na.rm = TRUE) + 3
xlim_d  <- c(val_min, val_max)                                   # Date axes
xlim_t  <- as.POSIXct(paste(xlim_d, "00:00:00"), tz = TZ_USE)    # datetime axis

# --------------------------- Loop over models --------------------------------
combo_per_model <- vector("list", length(model_keys))
names(combo_per_model) <- model_keys
model_metrics <- list()   # accumulate the comparison table

for (i in seq_along(model_keys)) {
  
  mname <- model_keys[i]
  mcol  <- model_cols[i]
  mobj  <- models[[mname]]
  
  # ---- Predict quantiles (log-space) ----
  log_pred <- predict(mobj, newdata = nd, what = percentiles)
  pred_logdf <- as.data.frame(log_pred)
  colnames(pred_logdf) <- q_names
  
  pred_logdf <- pred_logdf %>%
    bind_cols(nd) %>%
    mutate(
      # as_datetime handles both POSIXct (the .rds case) and character without
      # dropping midnight rows. ymd_hms() turned "2025-07-01 00:00:00" into the
      # bare string "2025-07-01" and failed to parse it, silently losing one row/day.
      DateTime = as_datetime(DateTime, tz = TZ_USE),
      Date     = as.Date(DateTime)
    )
  
  # ---- 30-min aggregation ----
  df_30min <- pred_logdf %>%
    mutate(bin_30min = floor_date(DateTime, unit = "30 minutes")) %>%
    group_by(bin_30min) %>%
    summarize(
      across(.cols = where(is.numeric), .fns = ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    mutate(
      bin_30min = as.POSIXct(bin_30min, tz = TZ_USE),
      date      = as.Date(bin_30min)
    )
  
  # ============================================================================
  # TOP PANEL: time series + DAILY metrics (daily mean lab vs daily mean model)
  # ============================================================================
  pred_daily <- df_30min %>%
    group_by(date) %>%
    summarize(q50_day = mean(q50, na.rm = TRUE), .groups = "drop")
  
  ec_daily <- ec %>%
    mutate(date = as.Date(date),
           logec = log10(ecoli)) %>%
    group_by(date) %>%
    summarize(logec_day = mean(logec, na.rm = TRUE), .groups = "drop")
  
  ts_eval_daily <- pred_daily %>%
    inner_join(ec_daily, by = "date") %>%
    filter(is.finite(q50_day), is.finite(logec_day))
  
  rmse_val <- sqrt(mean((ts_eval_daily$q50_day - ts_eval_daily$logec_day)^2, na.rm = TRUE))
  
  sse <- sum((ts_eval_daily$logec_day - ts_eval_daily$q50_day)^2, na.rm = TRUE)
  sst <- sum((ts_eval_daily$logec_day - mean(ts_eval_daily$logec_day, na.rm = TRUE))^2, na.rm = TRUE)
  r2_val <- if (isTRUE(all.equal(sst, 0))) NA_real_ else 1 - (sse / sst)
  r_ts   <- if (nrow(ts_eval_daily) >= 3)
              cor(ts_eval_daily$logec_day, ts_eval_daily$q50_day) else NA_real_

  lab_txt <- paste0(
    "Daily mean\n",
    "Pearson r = ", ifelse(is.na(r_ts), "NA", sprintf("%.2f", r_ts)), "\n",
    "RMSE = ", sprintf("%.2f", rmse_val), "\n",
    "R2 = ", ifelse(is.na(r2_val), "NA", sprintf("%.2f", r2_val)), "\n",
    "N days = ", nrow(ts_eval_daily)
  )

  # annotation placed at the right edge of the zoom window
  x_anno <- xlim_t[2] - days(1)
  y_anno <- max(df_30min$q90, na.rm = TRUE)
  
  p_ts <- ggplot(df_30min, aes(x = bin_30min, y = q50, ymin = q10, ymax = q90)) +
    geom_line(col = mcol, alpha = 1) +
    geom_ribbon(fill = mcol, alpha = 0.12) +
    geom_point(
      data = ec, inherit.aes = FALSE,
      mapping = aes(x = date, y = log10(ecoli)),
      col = "red", size = 2, alpha = 1
    ) +
    scale_x_datetime(limits = xlim_t) +
    theme_grey(base_size = 18) +
    labs(x = NULL, y = "E.coli/100mL [lg]") +
    ggtitle(paste0("Labormessung vs Modellvorhersage (", mname, ")")) +
    annotate(
      "label",
      x = x_anno, y = y_anno,
      label = lab_txt,
      hjust = 1, vjust = 1,
      size = 4.5,
      label.size = 0.25,
      colour = mcol,
      fill = "white",
      alpha = 0.9
    )
  
  # ============================================================================
  # MIDDLE PANEL: applied dose schedule (from df_intervals)
  # ============================================================================
  dose_ts <- make_dose_ts(
    start_dates_df = df_intervals,
    date_from = min(df_30min$date, na.rm = TRUE),
    date_to   = max(df_30min$date, na.rm = TRUE),
    dose_col  = "J_m2"     # change to "W_m2" or "m3_h" or "pct" if desired
  )
  
  p_dose <- ggplot(dose_ts, aes(x = day, y = dose)) +
    geom_step(colour = mcol, linewidth = 1) +
    scale_x_date(limits = xlim_d) +
    theme_grey(base_size = 16) +
    labs(x = NULL, y = "UV Dose [J/m2]") +
    ggtitle("Applied UV dose (interval schedule)")
  
  # ============================================================================
  # BOTTOM PANEL: UV validation (your existing logic, with model-specific color)
  # ============================================================================
  # Note: daily must contain `Datum` (Date) and `dose` (numeric)
  if (!all(c("Datum", "dose") %in% names(daily))) {
    stop("`daily` must contain columns `Datum` (Date) and `dose` (numeric).")
  }
  daily$Datum <- as.Date(daily$Datum)
  
  lrv <- apply(
    exp(posterior_epred(stanDR, newdata = data.frame(logdose = log(daily$dose)))),
    2, quantile, probs = c(0.025, 0.5, 0.95)
  ) %>%
    t() %>%
    as.data.frame()
  
  daily$LRV <- lrv$`2.5%`
  
  rmse_prep <- df_30min %>%
    group_by(date) %>%
    reframe(
      q10  = mean(q10),  q50  = mean(q50),  q60  = mean(q60),  q70  = mean(q70),
      q80  = mean(q80),  q90  = mean(q90),  q95  = mean(q95),  q100 = mean(q100)
    )
  
  eff_df <- rmse_prep %>% inner_join(daily, by = c("date" = "Datum"))
  # gather the percentile columns (q*) to long; keep everything else fixed.
  # (Robust to column count — the original hardcoded indices c(1,10:17) broke
  #  when `daily` changed width.)
  remain <- setdiff(colnames(eff_df), q_names)
  eff_long <- eff_df %>% gather(percentile, value, -all_of(remain))
  eff_long$ec_eff <- eff_long$value - eff_long$LRV
  
  eff_long2 <- eff_long %>%
    mutate(y = 10^ec_eff) %>%
    arrange(percentile, date) %>%
    filter(!(percentile %in% c("q100")))

  # ---- Pearson correlation: predicted post-UV (q50) vs measured AquaBio ----
  uv_pred  <- eff_df %>% transmute(date = as.Date(date), pred_uv = 10^(q50 - LRV))
  # as.numeric guards against fread's integer64 (mean() on it yields garbage).
  ab_daily <- ablauf_uv %>% mutate(d = as.Date(date), ecoli = as.numeric(ecoli)) %>%
    group_by(d) %>% summarize(meas = mean(ecoli, na.rm = TRUE), .groups = "drop")
  uv_cmp <- inner_join(uv_pred, ab_daily, by = c("date" = "d")) %>%
    filter(is.finite(pred_uv), is.finite(meas))
  uv_r <- if (nrow(uv_cmp) >= 3) cor(uv_cmp$meas, uv_cmp$pred_uv) else NA_real_
  uv_lab <- paste0("Pearson r = ", ifelse(is.na(uv_r), "NA", sprintf("%.2f", uv_r)),
                   "\nN = ", nrow(uv_cmp))

  p_uv <- ggplot(eff_long2, aes(x = date, y = 10^ec_eff)) +
    geom_ribbon(
      aes(ymin = 0, ymax = 10^ec_eff, group = percentile),
      alpha = 0.2, fill = mcol, colour = NA
    ) +
    geom_point(
      data = ablauf_uv,
      mapping = aes(x = as.Date(date), y = ecoli),
      colour = "black", pch = 21, bg = "white", size = 3,
      inherit.aes = FALSE
    ) +
    annotate("label", x = xlim_d[1], y = 290, hjust = 0, vjust = 1,
             label = uv_lab, size = 4.5, colour = mcol, fill = "white",
             alpha = 0.9, label.size = 0.25) +
    scale_x_date(limits = xlim_d) +
    scale_y_continuous(limits = c(0, 300)) +
    theme_gray(base_size = 16) +
    labs(x = "Datum", y = "Vorhergesagte E.coli Konzentration") +
    ggtitle("UV validation")

  # ============================================================================
  # EXTRA standalone figure: UV validation with dose as grey background
  # (validation_UV_i.png) — dose shown behind the ribbon + AquaBio points.
  # ============================================================================
  dose_bg <- daily %>% transmute(xmin = as.Date(Datum), xmax = as.Date(Datum) + 1, dose = dose)

  p_uv_grey <- ggplot(eff_long2, aes(x = date, y = 10^ec_eff)) +
    geom_rect(data = dose_bg,
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = dose),
              inherit.aes = FALSE, alpha = 0.35) +
    scale_fill_gradient(name = "UV dose", low = "grey95", high = "grey20") +
    geom_ribbon(aes(ymin = 0, ymax = 10^ec_eff, group = percentile),
                alpha = 0.20, fill = mcol, show.legend = FALSE, colour = NA) +
    geom_point(data = ablauf_uv, aes(x = as.Date(date), y = ecoli),
               colour = "black", pch = 21, bg = "white", size = 3, inherit.aes = FALSE) +
    annotate("label", x = xlim_d[1], y = 290, hjust = 0, vjust = 1,
             label = uv_lab, size = 5, colour = mcol, fill = "white",
             alpha = 0.9, label.size = 0.25) +
    scale_x_date(limits = xlim_d, date_breaks = "month", date_labels = "%b %Y") +
    scale_y_continuous(breaks = seq(0, 6000, 50), limits = c(0, 300)) +
    theme_gray(base_size = 20) +
    labs(x = "Datum", y = "Vorhergesagte E.coli Konzentration") +
    ggtitle(paste0("UV / E.coli Vorhersage (", mname, ")"))

  ggsave(file.path(OUT_DIR, paste0("validation_UV_", i, ".png")),
         p_uv_grey, width = 16, height = 7, dpi = "retina")
  
  # ============================================================================
  # Combine 3 panels — width-align the plot areas so the shared time axis lines
  # up vertically across panels (gtable::unit.pmax on the grob widths).
  # ============================================================================
  g_ts <- ggplotGrob(p_ts); g_dose <- ggplotGrob(p_dose); g_uv <- ggplotGrob(p_uv)
  max_w <- grid::unit.pmax(g_ts$widths, g_dose$widths, g_uv$widths)
  g_ts$widths <- max_w; g_dose$widths <- max_w; g_uv$widths <- max_w

  combo <- gridExtra::arrangeGrob(
    g_ts, g_dose, g_uv,
    ncol = 1,
    heights = c(2.2, 1.0, 2.0)
  )
  
  combo_per_model[[mname]] <- combo

  # ---- collect metrics for the comparison table ----
  model_metrics[[mname]] <- data.frame(
    model         = mname,
    n_vars        = length(rownames(mobj$importance)),
    RMSE_inlet    = round(rmse_val, 3),
    R2_inlet      = round(r2_val, 3),
    Pearson_inlet = round(r_ts, 3),
    N_days_inlet  = nrow(ts_eval_daily),
    Pearson_UV    = round(uv_r, 3),
    N_UV          = nrow(uv_cmp),
    stringsAsFactors = FALSE
  )
  
  ggsave(
    filename = file.path(OUT_DIR, paste0("model_", i, "_3panel.png")),
    plot = combo,
    width = 18, height = 18, dpi = "retina"
  )
}

# Optional: also export all models into a single multi-page PDF (one page per model)
pdf(file.path(OUT_DIR, "all_models_3panel.pdf"), width = 18, height = 18)
for (i in seq_along(combo_per_model)) {
  if (i > 1) grid::grid.newpage()          # one page per model (was overlapping on 1 page)
  grid::grid.draw(combo_per_model[[i]])
}
dev.off()

# ---- model comparison table -------------------------------------------------
model_comparison <- do.call(rbind, model_metrics)
rownames(model_comparison) <- NULL
write.csv(model_comparison, file.path(OUT_DIR, "model_comparison.csv"), row.names = FALSE)
cat("\n================ MODEL COMPARISON ================\n")
print(model_comparison, row.names = FALSE)
cat("(Pearson_inlet/RMSE_inlet/R2_inlet: daily lab vs prediction; Pearson_UV: post-UV pred vs AquaBio)\n")

message("Done. Created per-model PNGs, all_models_3panel.pdf, model_comparison.csv in: ", normalizePath(OUT_DIR))