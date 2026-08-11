source("R/config.R")
library(tidyverse)
library(readxl)
library(data.table)
library(brms)

zu <- data.table::fread(ZULAUF_KA_CSV)

ab <- data.table::fread(ABLAUF_UV_CSV)



ab$E_Coli <- ifelse(ab$E_Coli == "< 1", 
                    0, 
                    ab$E_Coli)

ab$E_Coli <- ifelse(ab$E_Coli == "< 15", 
                    NA, 
                    ab$E_Coli)

ab$E_Coli <- as.numeric(ab$E_Coli)

ab$C_Perfringens_Sporen <- ifelse(ab$C_Perfringens_Sporen == "< 1", 
                                  0, 
                                  ab$C_Perfringens_Sporen) %>% as.numeric()

ab$Somatische_Coliphagen <- str_replace_all(ab$Somatische_Coliphagen, 
                                            "< ", 
                                            "") %>% as.numeric()

ab$Somatische_Coliphagen <- ifelse(ab$Somatische_Coliphagen == 1, 
                                   0, 
                                   ab$Somatische_Coliphagen)

#ab$Somatische_Coliphagen <- ifelse(ab$Somatische_Coliphagen > 0, 
 #                                  NA, 
  #                                 ab$Somatische_Coliphagen)


ab$Datum_Start_Probenahme <- dmy(ab$Datum_Start_Probenahme)
zu$Datum_Start_Probenahme <- dmy(zu$Datum_Start_Probenahme)
zu$E_Coli <- ifelse(str_detect(zu$E_Coli, pattern = ">"), NA, zu$E_Coli) %>% as.numeric()



l <- list()

for(i in c("E_Coli", "C_Perfringens_Sporen", "Somatische_Coliphagen")){
  
  param_ab <- ab %>% select(any_of(c("Datum_Start_Probenahme", i))) %>% rename(ablauf = i)
  param_zu <- zu %>% select(any_of(c("Datum_Start_Probenahme", i)))%>% rename(zulauf = i)

  l[[i]] <- param_zu %>% full_join(param_ab, by = "Datum_Start_Probenahme")

}



for(j in c("E_Coli", "C_Perfringens_Sporen", "Somatische_Coliphagen")){
  
  datentabelle <- data.frame(zulaufwerte = l[[j]]$zulauf,
                             
                             ablaufwerte = l[[j]]$ablauf)
  
  
  
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
  
  Perzentilverteilung <- c(1:nrow(MCMC_zulauf))
  
  
  
  for(i in 1:nrow(MCMC_zulauf)){
    
    # Generieren von 10000 Zufallswerten auf Basis der i-ten Zeile der MCMC-Werte
    
    
    
    # Zulauf
    
    in_sim <- rlnorm(n = 10000, 
                     
                     meanlog = MCMC_zulauf$b_Intercept[i],
                     
                     sdlog = MCMC_zulauf$sigma[i])
    
    # Ablauf
    
    out_sim <- rnbinom(10000,
                       
                       mu = exp(MCMC_ablauf$b_Intercept[i]),
                       
                       size = MCMC_ablauf$sigma[i])+1
    
    
    
    # Berechnugn der Log10-Reduktionsverteilung aus 10000 Simulationswerten 
    
    lrv <- log10(in_sim/out_sim)
    
    
    
    # Berechnung des 10. und 90. Perzentils
    
    Perzentilverteilung[i] <- quantile(lrv, probs = c(0.1))
    
  }
  
  
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
            subtitle = paste("Indikator =", j))
 
    ggsave(plot = p, width = 12, height = 9, dpi = "retina", filename = paste0("plots/preVal_", j, ".png")) 

}

hist((posterior_predict(Ablaufverteilung)))
