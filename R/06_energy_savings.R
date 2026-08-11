# ============================================================================
# 06_energy_savings.R — recommended UV dose per prediction + energy savings.
#   Turns the model's predicted inlet E. coli into the UV dose required to meet
#   the target, then quantifies the saving vs. a constant high-dose regime.
#   INPUTS (auto-loaded): models (02), df_roll3h_only/ec_clean (02)
#   Figures: plots/timeplot.png, plots/UVdose.png
#   Also includes the standalone parametric saving model (was Estimate_enegysavings.R).
# ============================================================================

source("R/config.R")
source("R/functions.R")
library(tidyverse)
library(quantregForest)
library(lubridate)
library(purrr)

# ---- STANDALONE INPUT ------------------------------------------------------
if (!exists("models")) load(model_path("models.RData"))                 # 02_train_models.R
if (!exists("df_roll3h_only")) {
  .mo <- load_derived("modelling_objects")                              # 02_train_models.R
  df_roll3h_only <- .mo$df_roll3h_only
  ec_clean       <- .mo$ec_clean
}

library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(purrr)

## 0. Prep newdata (same for all models)
##    You were doing: newdata <- na.omit(df_roll3h_only)
##    We'll keep that and also pre-join ec_clean once so we don't redo it per model.

newdata_base <- df_roll3h_only %>%
  na.omit() 

r <- list()
for(model in names(models)){
  pred <- predict(models[[model]], newdata = newdata_base, what = c(0.1, 0.5, 0.9))
  r[[model]] <- as.data.frame(pred)
  colnames(r[[model]]) <- c("q10", "q50", "q90") 
}

nd <- newdata_base %>% left_join(ec_clean)

for(model in names(r)){
  r[[model]] <- r[[model]]%>% bind_cols(nd)
  r[[model]]$model <- model
  
}


daily <- r

for(model in names(r)){
  daily[[model]]$date <- lubridate::date(r[[model]]$DateTime)
  daily[[model]] <- daily[[model]]%>% group_by(date) %>% reframe(m = mean(q50),
                                                                 l = mean(q10),
                                                                 u = mean(q90),
                                                                 ec = mean(E_Coli,na.rm = T ))
  daily[[model]]$model <- model
}

predictions <- bind_rows(daily)

predictions$model <- factor(predictions$model, ordered = T,
                            levels = c("bestR2", "best7", "best4of7"),
                            labels = c("bestR2 (11 Var.)", "best7 (6 Var.)", "best4of7 (3 Var.)"))

timeplot <- ggplot(predictions, aes(date, y = m, ymin = l, ymax = u, col = model, fill = model))+ 
  #geom_linerange(col = "steelblue") +
  geom_line(lwd = 1.3)+ 
  geom_ribbon(alpha = 0.3, color = NA)+ 
  theme_gray(base_size = 16) +
  scale_color_manual(values = c("blue", "purple", "forestgreen"))+
  scale_fill_manual(values = c("blue", "purple", "forestgreen"))+
  facet_grid(model~.)+
  scale_x_date(date_breaks = "2 months")+ geom_point(mapping = aes(date, log10(ec)), 
                                                     bg = "white", 
                                                     pch = 22,
                                                     col = "black", 
                                                     size = 3, 
                                                     alpha = 0.9)+
  ggtitle("Zeitreihe E.coli Vorhersage vs. Messwerte") + labs(x = "Datum", y = "E.coli [log10]")+
  theme(legend.position = "bottom")

timeplot
ggsave("plots/timeplot.png", plot = timeplot, 
       width = 20, 
       height = 9, dpi = "retina")

ggplot(predictions, aes(date, y = m, ymin = l, ymax = u, col = model, fill = model))+ 
  #geom_linerange(col = "steelblue") +
  geom_line(lwd = 1.3)+ 
#  geom_ribbon(alpha = 0.3, color = NA)+ 
  theme_gray(base_size = 16) +
#  facet_grid(model~.)+
  scale_x_date(date_breaks = "2 months")+ geom_point(mapping = aes(date, log10(ec)), 
                                                     bg = "white", 
                                                     pch = 22,
                                                     col = "black", 
                                                     size = 3, 
                                                     alpha = 0.9)+
  ggtitle("Zeitreihe E.coli Vorhersage vs. Messwerte") + labs(x = "Datum", y = "E.coli [log10]")




predictions$zielwert <- 50
predictions$rLRVm <- predictions$m-log10(predictions$zielwert)
predictions$rLRVu <- predictions$u-log10(predictions$zielwert)

predictions$rUVm <- exp(4.1+1.7967*log(predictions$rLRVm))
predictions$rUVu <- exp(4.1+1.7967*log(predictions$rLRVu))



uv_dose <- ggplot(predictions, aes(date, ymin = rUVm, ymax = rUVu, col = model, fill = model))+ 
  #geom_linerange(col = "steelblue") +
  #geom_line(lwd = 1.3)+ 
  geom_ribbon(alpha = 0.6, color = NA)+ 
  theme_gray(base_size = 16) +
  scale_fill_manual(values = c("blue", "purple", "forestgreen"))+
  facet_grid(model~.)+
  scale_x_date(date_breaks = "2 months")+ 
  ggtitle("Zeitreihe empfohlene UV Dosis") + 
  labs(x = "Datum", y = "UV Dosis PSS [J/m²]")+
  theme(legend.position = "bottom")

ggsave("plots/UVdose.png", plot = uv_dose, 
       width = 20, 
       height = 9, dpi = "retina")

energie <- predictions %>% group_by(model) %>% reframe(avg_UV = mean(rUVm),
                                                       upperUV = mean(rUVu))

energie$ec_q95<- quantile(log10(predictions$ec), na.rm = T, probs = .95)

energie$constant <- exp(4.1 + 1.79*log(energie$ec_q95 - log10(50)))

energie$reduction_m <- 1-energie$avg_UV/energie$constant
energie$reduction_u <- 1-energie$upperUV/energie$constant



# ============================================================================
# Parametric energy-saving model (dose reduction -> energy saving)
# (was Estimate_enegysavings.R — independent of the fitted models)
# ============================================================================


# parameters you can tune
f <- 0.30  # fixed fraction of UV energy demand
v <- 0.70  # variable fraction (1 - f)

# sequence of relative dose levels (r)
r <- seq(1, 0.4, by = -0.1)  # from 100% dose down to 40%

# compute relative energy use and savings
E_rel <- f + v * r                # fraction of baseline energy
saving <- 1 - E_rel               # fraction saved

result <- data.frame(
  dose_reduction_percent   = (1 - r) * 100,
  relative_energy_use      = E_rel * 100,
  energy_saving_percent    = saving * 100
)

print(result)

# plot: dose reduction (%) on x, energy saving (%) on y
plot(result$dose_reduction_percent,
     result$energy_saving_percent,
     xlab = "UV dose reduction (%)",
     ylab = "Energy saving (%)",
     main = "Expected Energy Saving vs. UV Dose Reduction\n(fixed = 30%, variable = 70%)",
     pch = 19)
lines(result$dose_reduction_percent,
      result$energy_saving_percent,
      lwd = 2)
