# ============================================================================
# 01_clean_plant_data.R — clean the raw plant CSV into modelling features.
#   QC report + timestamp parsing + duplicate handling + linear interpolation
#   + 3h trailing rolling means + canonical channel names.
#   INPUT : data/klaerwerkNOVD.csv (config: KLAERWERK_NOVD_CSV / in_file)
#   OUTPUT: data/KlaerwerkS2025_clean_interpolated*.csv  and
#           data/derived/v_roll3h_only.rds  (scoring features for stage 05)
# ============================================================================


  source("R/config.R")
  source("R/functions.R")
  library(readr)
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(stringr)
  library(zoo)
  library(slider)

# -----------------------------
# FILES (edit if desired)
# -----------------------------
in_file  <- "data/klaerwerkNOVD.csv"
out_clean_interp <- "data/KlaerwerkS2025_clean_interpolated.csv"
out_roll3h       <- "data/KlaerwerkS2025_clean_interpolated_roll3h.csv"

tz_used <- "Europe/Berlin"

# Rolling window settings
roll_window <- hours(3)   # 3-hour trailing mean

# Interpolation settings
# If you want to limit interpolation over long gaps, set e.g. max_gap_hours <- 6
# and the script will leave NAs for gaps larger than that (based on expected step).
max_gap_hours <- NA_real_   # NA => no max-gap limit; interpolate all internal gaps

# -----------------------------
# Helpers
# -----------------------------
mean_na_safe <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

to_num <- function(x) {
  # Convert German decimal comma to dot; handle blanks; keep NA if not numeric
  x <- str_trim(as.character(x))
  x[x == ""] <- NA_character_
  suppressWarnings(as.numeric(str_replace_all(x, ",", ".")))
}

longest_na_run <- function(v) {
  r <- rle(is.na(v))
  if (!any(r$values)) return(0L)
  max(r$lengths[r$values])
}

# -----------------------------
# 1) Read raw (all as character), drop empty rows
# -----------------------------
raw <- readr::read_delim(
  in_file,
  delim = ";",
  col_types = cols(.default = col_character()),
  locale = locale(encoding = "UTF-8")
)

# Keep only rows with Datum and Zeit present (your file contains many trailing empty rows)
raw <- raw %>%
  filter(!is.na(Datum), !is.na(Zeit), str_trim(Datum) != "", str_trim(Zeit) != "")

# -----------------------------
# 2) Build timestamp and convert numerics
# -----------------------------
df <- raw %>%
  mutate(
    datetime = dmy(Datum, tz = tz_used) + hm(Zeit),
    datetime = as.POSIXct(datetime, tz = tz_used)
  )

# Identify numeric candidate columns (everything except Datum/Zeit/datetime)
value_cols <- setdiff(names(df), c("Datum", "Zeit", "datetime"))

df <- df %>%
  mutate(across(all_of(value_cols), to_num)) %>%
  select(datetime, all_of(value_cols)) %>%
  arrange(datetime)

# -----------------------------
# 3) Validation checks
# -----------------------------
cat("\n================ DATA QUALITY REPORT ================\n")

# Timestamp parsing failures
bad_time <- sum(is.na(df$datetime))
cat(sprintf("Rows after dropping empty lines: %d\n", nrow(df)))
cat(sprintf("Unparsed datetime rows: %d\n", bad_time))

if (bad_time > 0) {
  cat("Dropping rows with NA datetime.\n")
  df <- df %>% filter(!is.na(datetime))
}

# Handle duplicate timestamps by averaging numeric values (robust default)
dup_n <- sum(duplicated(df$datetime))
cat(sprintf("Duplicate timestamps: %d\n", dup_n))

if (dup_n > 0) {
  df <- df %>%
    group_by(datetime) %>%
    summarise(across(where(is.numeric), ~ if (all(is.na(.x))) NA_real_ else mean(.x, na.rm = TRUE)),
              .groups = "drop") %>%
    arrange(datetime)
  cat("Duplicates aggregated by mean (numeric columns).\n")
}

# Step size summary (your data is mostly 15-min = 900 seconds)
dt <- diff(df$datetime)
if (length(dt) > 0) {
  dt_sec <- as.numeric(dt)
  cat(sprintf("Time span: %s to %s\n",
              format(min(df$datetime), "%Y-%m-%d %H:%M:%S %Z"),
              format(max(df$datetime), "%Y-%m-%d %H:%M:%S %Z")))
  cat(sprintf("Median step (sec): %.0f\n", median(dt_sec, na.rm = TRUE)))
  cat(sprintf("Max step (sec): %.0f\n", max(dt_sec, na.rm = TRUE)))
  cat(sprintf("Non-15min steps: %d\n", sum(dt_sec != 900, na.rm = TRUE)))
}

# Missingness and longest missing run per variable
miss_tbl <- df %>%
  summarise(across(where(is.numeric), ~ mean(is.na(.x)))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "missing_frac") %>%
  arrange(desc(missing_frac))

run_tbl <- tibble(
  variable = names(df)[sapply(df, is.numeric)],
  longest_na_run_obs = sapply(df[names(df)[sapply(df, is.numeric)]], longest_na_run)
) %>%
  # Convert run length to hours using median step (fallback 15min)
  mutate(step_sec = if (length(dt) > 0) median(as.numeric(dt), na.rm = TRUE) else 900,
         longest_na_run_hours = (longest_na_run_obs * step_sec) / 3600) %>%
  arrange(desc(longest_na_run_hours))

cat("\nTop missing fractions (numeric):\n")
print(head(miss_tbl, 15))

cat("\nLongest NA runs (hours, numeric):\n")
print(head(run_tbl, 15))

# Generic negative-value check (common sanity check for sensors/flows)
neg_tbl <- df %>%
  summarise(across(where(is.numeric), ~ sum(.x < 0, na.rm = TRUE))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "n_negative") %>%
  filter(n_negative > 0) %>%
  arrange(desc(n_negative))

if (nrow(neg_tbl) == 0) {
  cat("\nNegative values: none detected.\n")
} else {
  cat("\nNegative values detected:\n")
  print(neg_tbl)
}

# -----------------------------
# 4) Linear interpolation (time-based)
# -----------------------------
# Optional max-gap limit: only interpolate gaps up to max_gap_hours.
# zoo::na.approx maxgap is in *number of observations*, so we estimate from median step.
step_sec <- if (length(dt) > 0) median(as.numeric(dt), na.rm = TRUE) else 900
maxgap_obs <- if (is.na(max_gap_hours)) Inf else ceiling((max_gap_hours * 3600) / step_sec)

df_interp <- df %>%
  mutate(across(
    where(is.numeric),
    ~ zoo::na.approx(.x,
                     x = as.numeric(datetime),
                     na.rm = FALSE,
                     rule = 1,
                     maxgap = maxgap_obs)
  ))

# Save interpolated-clean dataset
write_csv2(df_interp, out_clean_interp)
cat(sprintf("\nWrote interpolated dataset: %s\n", out_clean_interp))

# -----------------------------
# 5) 3-hour rolling averages (time-based trailing window)
# -----------------------------
df_roll <- df_interp %>%
  mutate(across(
    where(is.numeric),
    ~ slider::slide_index_dbl(
      .x,
      .i = datetime,
      .f = mean_na_safe,
      .before = roll_window,
      .complete = FALSE
    ),
    .names = "{.col}_roll3h"
  ))

write_csv2(df_roll, out_roll3h)
cat(sprintf("Wrote dataset with 3h rolling means: %s\n", out_roll3h))

cat("\nDone.\n")


v_roll3h_only <- df_roll %>%
  select(datetime, contains("_roll3h"))

colnames(v_roll3h_only) <- CHANNELS_ROLL3H    # canonical names (functions.R)
v_roll3h_only[CHANNELS_DROP] <- NULL          # drop redundant/non-predictive channels

# STANDALONE OUTPUT — prediction newdata consumed by
# ValidationPlots.R and Modellvergleich.R
save_derived(v_roll3h_only, "v_roll3h_only")
