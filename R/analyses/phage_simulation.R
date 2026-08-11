library(brms)     # For Bayesian modeling
library(tidyverse)  # For data wrangling

# Original data: 13 samples of 100 mL with observed phage counts


phage_counts <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 7, 10, 0)
data <- data.frame(
  count = phage_counts,
  volume = rep(100, length(phage_counts))  # all 100 mL samples
)

# Step 2: Fit a negative binomial model using brms
# We use a simple intercept-only model since no covariates are provided

fit <- brm(
  formula = count ~ 1,
  family = "negbinomial2",
  data = data,
  chains = 4,
  iter = 4000,
  seed = 123
)

# Step 3: Draw posterior predictive samples
# Simulate a large number of 100 mL samples from the posterior predictive distribution

n_fake_samples <- 10000  # More samples = smoother aggregation distribution

pp_samples <- posterior_predict(fit, draws = n_fake_samples)
# Each row = one posterior draw of the full dataset. We want a flat vector of fake 100 mL counts:
simulated_100ml <- as.vector(pp_samples)

# Step 4: Simulate 270 mL samples by aggregating 2.7 × 100 mL units
set.seed(42)

simulated_270ml <- replicate(10000, {
  full_units <- sample(simulated_100ml, 2, replace = TRUE)
  partial_unit <- sample(simulated_100ml, 1, replace = TRUE) * 0.7
  round(sum(full_units) + partial_unit)
})

# Step 5: Plot the distribution
hist(log10(simulated_270ml), breaks = 30, main = "Posterior Predictive Distribution (270 mL)",
     xlab = "Phage Count per 270 mL", col = "lightgreen", border = "gray")



l$Somatische_Coliphagen$ablauf <- ifelse(l$Somatische_Coliphagen$ablauf>0, NA, l$Somatische_Coliphagen$ablauf)

datentabelle <- data.frame(zulaufwerte = l$Somatische_Coliphagen$zulauf,
                           
                           ablaufwerte = l$Somatische_Coliphagen$ablauf)



control <- list(adapt_delta = 0.99)

MCMC_iterations <- 4000



Zulaufverteilung <- brms::brm(zulaufwerte ~ 1, 
                              
                              data = datentabelle, 
                              
                              family = "lognormal",
                              
                              iter = MCMC_iterations,
                              
                              control = control)


Ablaufverteilung <- brms::brm(ablaufwerte ~ 1, 
                              
                              data = datentabelle,
                              
                              family = "negbinomial2",
                              
                              iter = MCMC_iterations, 
                              
                              control = control)


# Extrahieren der MCMC Werte von den jeweiligen Zu-und Ablaufverteilungen

MCMC_zulauf <-  brms::as_draws_df(Zulaufverteilung)

MCMC_ablauf <- brms::as_draws_df(Ablaufverteilung)

Perzentilverteilung <- c(1:100)#nrow(MCMC_zulauf))



for(i in 1:100){
  
  # Generieren von 10000 Zufallswerten auf Basis der i-ten Zeile der MCMC-Werte
  
  
  
  # Zulauf
  
  in_sim <- rlnorm(n = 10000, 
                   
                   meanlog = MCMC_zulauf$b_Intercept[i],
                   
                   sdlog = MCMC_zulauf$sigma[i])*4.7
  
  # Ablauf
  
  out_sim <- rnbinom(10000,
                     
                     mu = exp(MCMC_ablauf$b_Intercept[i]),
                     
                     size = MCMC_ablauf$sigma[i])
  
  
  
  simulated_270ml <- replicate(10000, {
    full_units <- sample(out_sim, 4, replace = TRUE)
    partial_unit <- sample(out_sim, 1, replace = TRUE) * 0.7
    round(sum(full_units) + partial_unit)
  })
  
  
  
  # Berechnugn der Log10-Reduktionsverteilung aus 10000 Simulationswerten 
  
  lrv <- log10(in_sim/(simulated_270ml+1))
  
  
  
  # Berechnung des 10. und 90. Perzentils
  
  Perzentilverteilung[i] <- quantile(lrv, probs = c(0.1))
  
}

Perzentilverteilungx <- Perzentilverteilung
library(ggplot2)





alpha_value <- 0.05

df <- data.frame(p = Perzentilverteilung)

lower <- quantile(df$p, probs = c(alpha_value))



p <-  ggplot(data = df, 
             
             aes(x = p)) + 
  
  geom_histogram(
    
    fill = "steelblue") + 
  
  geom_vline(xintercept = quantile(df$p, 
                                   
                                   probs = c(alpha_value, 1-alpha_value)), 
             
             lwd = 1, 
             
             lty = 2, 
             
             col = "red3") + 
  
  annotate("text", x = quantile(df$p, 
                                
                                probs = c(alpha_value))-0.3,
           
           col = "red3",
           
           size = 4,
           
           #fontface="bold",
           
           y = 700, 
           
           label = bquote("TI"[.("U")]~"(P = 0,9"~alpha~"= 0,05)"~ "="*.(round(lower,1)))
           
  ) +
  
  labs(x = bquote("Log" [.("10")] ~ "Reduktion"), y= "Anzahl Simulationen") +
  
  ggtitle("Verteilung des 10. Perzentils als Histogram", 
          subtitle = paste("Indikator =", "Phages_sim"))

p
ggsave(plot = p, width = 12, height = 9, dpi = "retina", filename = paste0("plots/preVal_", "phages_sim", ".png")) 


