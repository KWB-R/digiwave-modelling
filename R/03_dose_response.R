# ============================================================================
# 03_dose_response.R — UV dose -> log-reduction (LRV) relationship.
#   Fits the bench-scale dose-response used to translate a predicted inlet
#   E. coli load into the UV dose required to meet the target.
#   INPUT : data/raw/plot-data.csv  (config: UV_DOSE_RESPONSE_CSV; cols dose, lrv)
#   OUTPUT: models/stanDR.rds       (rstanarm dose-response model)
#   Figure: plots/uv_dose_response.png
# (Extracted from the original Laboratory.R, which also held unrelated lab plots.)
# ============================================================================

source("R/config.R")
source("R/functions.R")
library(tidyverse)
library(data.table)
library(rstanarm)

uv <- fread(UV_DOSE_RESPONSE_CSV)
names(uv) <- c("dose", "lrv")

m2 <- lm(lrv ~ log(dose), data = uv)

library(rstanarm)

uvlog <- data.frame(logdose = log(uv$dose), loglrv = log(uv$lrv))

m3 <- lm(logdose ~ loglrv, data = uvlog)
m4 <- lm(loglrv ~ logdose, data = uvlog)

stanDR <- rstanarm::stan_glm(loglrv ~ logdose, data = uvlog)

# STANDALONE OUTPUT — dose-response model consumed by
# CleanNewKAData.R, ValidationPlots.R, Modellvergleich.R
saveRDS(stanDR, model_path("stanDR.rds"))


hist(bayes_R2(stanDR))

summary(m4)

newdata <- data.frame(dose = c(725,395,238,
452,
478,
446,
409,
725,
588,
484,
473))

predict(m2, newdata = newdata)


summary(m3)

summary(m2)

residuals(m2)

predictions <- predict(m2, newdata = data.frame(dose = 100:1000))


dplot <- data.frame(dose = 100:1000, predictions = predictions)


uv_plot <- ggplot(dplot, aes(x = dose, y = predictions)) + 
  geom_line(col = "blue", size = 1)  + 
  theme_minimal(base_size = 33)+
  geom_point(data = uv, 
             col = "blue", 
             size = 4, 
             mapping = aes(x = dose, y = lrv), 
             inherit.aes = F) + 
  labs(x = "UV Dose [J/m²]", 
       y = "Log10-Reduktion") + 
  #theme_gray(base_size = 22) + 
  scale_y_continuous(breaks = seq(0, 5.5, 1), 
                     limits = c(0, 5.5))+
  scale_x_continuous(breaks = seq(0, 1000, 100), 
                     limits = c(0, 1000))


uv_plot
ggsave(plot = uv_plot, filename = "plots/uv_dose_response.png", 
       width = 13, 
       height = 9, 
       dpi = "retina")
