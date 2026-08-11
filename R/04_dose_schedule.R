# ============================================================================
# 04_dose_schedule.R — reconstruct the applied UV-dose schedule at the plant.
#   Parses the logged UV interval settings (m³/h, W/m²) and derives the daily
#   applied dose (J/m²). Consumed by 05/06 to validate predictions vs. real UV runs.
#   OUTPUT: data/derived/dose_schedule.rds  (list: daily [Datum, dose, ...], df)
# ============================================================================

source("R/config.R")
source("R/functions.R")

library(dplyr)
library(tidyr)
library(lubridate)
library(zoo)
library(readr)

# UV interval settings as logged at the plant (date, m³/h, W/m², J/m², %).
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

df <- read_tsv(
  I(txt),
  col_names = c("Datum", "m3_h", "W_m2", "J_m2", "pct"),  # ASCII names (multibyte names break across())
  col_types = cols(.default = col_character()),
  locale = locale(decimal_mark = ",")
) %>%
  mutate(
    Datum = dmy(Datum),
    across(-Datum, ~ suppressWarnings(parse_number(.x,
                                                   locale = locale(decimal_mark = ","))))
  ) %>%
  arrange(Datum)

# Expand to one row per day, carrying each setting forward until it next changes.
daily <- df %>%
  complete(Datum = seq(min(Datum), max(Datum), by = "day")) %>%
  arrange(Datum) %>%
  mutate(across(-Datum, ~ na.locf(.x, na.rm = FALSE)))

# Applied dose (J/m²) from lamp irradiance and contact time: 0.008 / (Q/3600) * W.
colnames(daily)[4] <- "dose"
daily$dose <- 0.008 / (daily$m3_h / 3600) * daily$W_m2

# ============================================================================
# STANDALONE OUTPUT — applied UV-dose schedule consumed by
# 05_predict_validate.R and 06_energy_savings.R (they reference `daily` and `df`)
#   daily : one row per day with columns Datum (Date) and dose (numeric)
#   df    : parsed UV interval table (Datum, m³/h, W/m², J/m², %)
# ============================================================================
save_derived(list(daily = daily, df = df), "dose_schedule")
