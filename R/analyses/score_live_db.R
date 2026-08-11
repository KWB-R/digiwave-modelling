source("R/config.R")
library(tidyverse)
library(DBI)
library(RSQLite)
library(data.table)

# STANDALONE INPUT — needs the fitted `models` from AccuracyVSnVars.R
if (!exists("models")) load(model_path("models.RData"))

con <- dbConnect(SQLite(), SINK_DB)


dbListTables(con)

dbListFields(con, "input")
log_df <- dbReadTable(con, "input")



log_df$DateTime <- as.POSIXct(
  substring(log_df$sampled_at,1, last = 19),
  format = "%Y-%m-%d %H:%M:%S",  # %OS6 = seconds with 6 digits after decimal
  tz = "Europe/Berlin"            # pick your timezone
)
colnames(log_df)

log_d <- fread("data/KlaerwerkS2025.csv")

log_d$DateTime <- dmy_hm(paste(log_d$Datum, log_d$Zeit))
log_d$Datum <- NULL
log_d$Zeit <- NULL
log_d$`Filtration, Leitfähigkeit` <- as.numeric(log_d$`Filtration, Leitfähigkeit`)
library(zoo)

# Example: df has a datetime column and multiple numeric columns
# Assume your datetime column is called "datetime"

window_size <- 12  # 12 x 15min = 3h

logdf_roll <- log_df %>% 
  arrange(DateTime) %>%
  mutate(across(
    where(is.numeric),
    ~ rollapply(.x, width = window_size, FUN = mean, align = "right", fill = NA),
    .names = "{.col}_roll3h"
  ))

logdf_roll3h_only <- logdf_roll %>%
  select(DateTime, contains("_roll3h"), -contains("Ecotouch"), -id_roll3h, -Filtration...GAK.Summe3_5_7_roll3h)


logdf_roll3h_only <- logdf_roll %>%
  select(DateTime, contains("_roll3h"))

colnames(logdf_roll3h_only)

colnames(logdf_roll3h_only) <- c("DateTime",
  "KA_Eff_NO3_roll3h",
  "KA_Eff_NH4_roll3h",   
  "KA_Eff_Turb_roll3h",  
  "KA_Eff_Q_roll3h",     
  "KA_Inf_LF_roll3h",    
  "KA_Inf_T_roll3h",  
  "KA_Inf_Q_roll3h",
  "Bio_Eff_T_roll3h",
  "Bio_Eff1_O2_roll3h",
  "Bio_Eff2_O2_roll3h",
  "Bio_Eff1_NH4_roll3h",
  "Bio_Eff2_NH4_roll3h",
  "Bio_Eff1_NO3_roll3h",
  "Bio_Eff2_NO3_roll3h",
  "Bio_Eff1_TR_roll3h",
  "Bio_Eff2_TR_roll3h",
  "GAK_Inf_Q_roll3h",
  "GAK_Eff_SAK_roll3h",
  "SF_Eff_SAK_roll3h",  
  "SF_Inf_LF_roll3h",
  "SF_Inf_TURB_roll3h" )


logdf_roll3h_only$KA_Inf_T_roll3h <- NULL
logdf_roll3h_only$SF_Eff_SAK_roll3h <- NULL
logdf_roll3h_only$GAK_Eff_SAK_roll3h <- NULL

nd <- na.omit(logdf_roll3h_only)

nd <- logdf_roll3h_only

check <- logdf_roll3h_only %>% gather(param, value, -DateTime)

ggplot(check, aes(x=DateTime, y = value)) + 
  geom_line()+ 
  facet_wrap(.~param, scale = "free_y")


log_pred <- predict(models$full, newdata = nd, what = c(0.1, 0.5, 0.9))

pred_logdf <- as.data.frame(log_pred)

colnames(pred_logdf) <- c("q10", "q50", "q90" )

#newdata <- newdata %>% left_join(ec_clean)

pred_logdf <- pred_logdf %>% bind_cols(nd)

pred_logdf$DateTime <- lubridate::ymd_hms(pred_logdf$DateTime)


pred_logdf$Date <- lubridate::date(pred_logdf$DateTime)


daily <- pred_logdf %>% group_by(Date) %>% reframe(m = mean(q50),
                                                l = mean(q10),
                                                u = mean(q90))

df_30min <- pred_logdf %>%
  mutate(bin_30min = floor_date(DateTime, unit = "30 minutes")) %>%
  group_by(bin_30min) %>%
  summarize(
    across(
      .cols = where(is.numeric),
      .fns  = ~ mean(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )

df_30min$date <- lubridate::date(df_30min$bin_30min)

library(data.table)
ec <- fread("data/uv_zulauf.csv")
ec$date <- lubridate::dmy(ec$date)
ec$date <- as.POSIXct(ec$date)

x <- ggplot(df_30min, aes(bin_30min, y = q50, ymin = q10, ymax = q90))+ 
  geom_point(col = "forestgreen", alpha = 1)+ 
  geom_linerange(col = "forestgreen", alpha = 0.02)+
  #geom_line(col = "blue", lwd = 1.3) +
  #geom_ribbon(alpha = 0.3, fill = "blue")+ theme_gray(base_size = 16) +
  scale_x_datetime() + 
  geom_point(data = ec, inherit.aes = F,
             mapping =  aes(date, y = log10(ecoli)), 
             col = "red", size = 5, alpha = 1)+
  labs(x = "Zeit", y = "E.coli/100mL [lg]") +
  theme_grey(base_size = 22)+
  ggtitle("Labormessung vs Modellvorhersage (Model 3)")

x
ggsave(filename = "plots/validation_model_3.png", plot = x, 
       width = 18, height = 9, dpi = "retina")  
  

