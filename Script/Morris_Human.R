  # --------------------------------------------------------------------------- #
  # Script for running the Morris Test
  # CP, 30-06-2025
  # --------------------------------------------------------------------------- #
  
  rm(list=ls()) # to clear out the global environment
  
  # Load packages
  library(here)
  library(ggplot2)
  library(deSolve)
  library(tidyverse)
  library(sensitivity)
  library(ggrepel)
  library(scales)
  library(patchwork)
  library(pracma)
  
  
  # For plotting
  library(showtext)
  font_add(family = "Garamond", regular = "GARA.TTF")
  showtext_auto()

  
  CP_theme <- theme_minimal() +
    theme(
      plot.background = element_rect(fill = 'white'), #element_blank(), #
      plot.margin = margin(0.3, 0.3, 0.3, 0.3, "cm"), 
      plot.title = element_text(hjust = 0, size = 48, face = "bold"),
      plot.title.position = "plot", 
      panel.background = element_blank(), #rect(fill = 'white', color = NA),
      panel.grid = element_blank(), #element_line(size = 0.3),  
      axis.title = element_text(size = 44, face = "bold"),
      axis.text = element_text(size = 44),
      axis.line = element_line(size = 0.5),
      legend.text = element_text(size = 44),
      legend.title = element_text(size = 44), 
      legend.key.size = unit(0.5, "lines"),
      legend.spacing.y = unit(0.5, "cm"), 
      legend.margin = margin(0,0,0,0),
      legend.position = "right"
    )
  
  # Set output storage directory
  OUTPUT <- here("Output", format(Sys.Date(), "%Y-%m-%d"), format(Sys.time(), "%H-%M-%S"))
  dir.create(OUTPUT, recursive = TRUE)
  
  # Model Code
  source(here("Script", "GSA","PBK_model_human.R")) #
  

  # SENSITIVITY FUNCTION ####
  # ---------------------------------------------------------------------------- #
  
  ## Morris test ####
  SENSI_model <- function(P){
    
    
    
    A_init <- c(
      AOD_ZEN = 0.55,
      
      ADuoL_ZEN = 0, AJejL_ZEN = 0, AIleL_ZEN = 0, 
      ADuo_ZEN = 0, AJej_ZEN = 0, AIle_ZEN = 0, 
      AL_ZEN = 0, 
      AR_ZEN = 0,
      AS_ZEN = 0,
      AF_ZEN = 0, 
      AB_ZEN = 0,
      AUR_ZEN = 0, 
      AFE_ZEN = 0, 
      
      A_bZEL = 0,
      
      ADuoL_ZEN_gluc = 0, AJejL_ZEN_gluc = 0, AIleL_ZEN_gluc = 0, 
      AJej_IC_ZEN_gluc = 0, AJej_EC_ZEN_gluc = 0, 
      AL_IC_ZEN_gluc = 0, AL_EC_ZEN_gluc = 0,
      ABile_ZEN_gluc = 0, 
      AR_ZEN_gluc = 0,
      AS_ZEN_gluc = 0,
      AF_ZEN_gluc = 0, 
      AB_ZEN_gluc = 0, 
      AK_ZEN_gluc = 0,
      AUR_ZEN_gluc = 0, 
      AFE_ZEN_gluc = 0,
      
      AIleL_aZEL = 0, 
      AJej_aZEL = 0, AIle_aZEL = 0,
      AL_aZEL = 0, 
      AR_aZEL = 0,
      AS_aZEL = 0,
      AF_aZEL = 0,
      AB_aZEL = 0,
      AUR_aZEL = 0,
      AFE_aZEL = 0,
      
      ADuoL_aZEL_gluc = 0, AJejL_aZEL_gluc = 0, AIleL_aZEL_gluc = 0, 
      AJej_IC_aZEL_gluc = 0, AJej_EC_aZEL_gluc = 0,
      AL_IC_aZEL_gluc = 0, AL_EC_aZEL_gluc = 0,
      ABile_aZEL_gluc = 0, 
      AR_aZEL_gluc = 0, 
      AS_aZEL_gluc = 0,
      AF_aZEL_gluc = 0,
      AB_aZEL_gluc = 0,
      AK_aZEL_gluc = 0,
      AUR_aZEL_gluc = 0,
      AFE_aZEL_gluc = 0
      
      
    )
    
    
    PBK.out <- lsoda(y = A_init,
                     times = seq(0, 24,by=1/60),
                     func = H.PBK.model, 
                     parms = P)
    PBK_OUTPUT <- as.data.frame(PBK.out)
    
    CmaxB_ZEN <- max(PBK.out[,"CB_ZEN"])
    AUC_ZEN <- trapz(PBK_OUTPUT[ , "time"], PBK_OUTPUT[ , "CB_ZEN"])
    
    CmaxB_ZEN_gluc <- max(PBK.out[,"CB_ZEN_gluc"])
    AUC_ZEN_gluc <- trapz(PBK_OUTPUT[ , "time"], PBK_OUTPUT[ , "CB_ZEN_gluc"])
    
    CmaxB_aZEL <- max(PBK.out[,"CB_aZEL"])
    AUC_aZEL <- trapz(PBK_OUTPUT[ , "time"], PBK_OUTPUT[ , "CB_aZEL"])
    
    CmaxB_aZEL_gluc <- max(PBK.out[,"CB_aZEL_gluc"])
    AUC_aZEL_gluc <- trapz(PBK_OUTPUT[ , "time"], PBK_OUTPUT[ , "CB_aZEL_gluc"])
    
    
    return( c(CmaxB_ZEN = CmaxB_ZEN, 
              AUC_ZEN = AUC_ZEN ,

              CmaxB_ZEN_gluc = CmaxB_ZEN_gluc,
              AUC_ZEN_gluc = AUC_ZEN_gluc,

              CmaxB_aZEL = CmaxB_aZEL,
              AUC_aZEL = AUC_aZEL,

              CmaxB_aZEL_gluc = CmaxB_aZEL_gluc,
              AUC_aZEL_gluc = AUC_aZEL_gluc
              ))
    
  }
  
  # SET MORRIS PARAMETERS ####
  # ---------------------------------------------------------------------------- #
  
  # Parameter upper and lower bounds
  
  # Parameters
  Parameters_df <- read.csv(here("Input", "Parameters_Human.csv")) 
  
  # Dosing
  PhysChem <- read_delim("Input/PhysChem.csv", 
                         delim = ";", escape_double = FALSE, trim_ws = TRUE)

  ParametersMorris <- Parameters_df 
  
  ParametersMorris$Value <- as.numeric(ParametersMorris$Value)

  P <- ParametersMorris$Value
  P <- ParametersMorris %>% column_to_rownames("Parameter")
  
  binf <- P$Value*0.9
  bsup <- P$Value*1.1
  
  Morris.factors <- ParametersMorris$Parameter  
  N.par <- nrow(ParametersMorris)
  
  
  # SET MORRIS EXPERIMENT DESIGN ####
  # ---------------------------------------------------------------------------- #
  
  
  # Experiment design
  Morris.design <- list(type = "oat", levels = 6, grid.jump = 3)
  
  # Repetitions
  Morris.r <- 5
  
  # Morris test function
  Experience <- morris(model =  NULL, 
                   factors = Morris.factors, # number of parameters
                   r = Morris.r,  # number of repetitions
                   design = Morris.design,  
                   binf = binf,  
                   bsup = bsup, 
                   scale = TRUE)
  
  # Save Morris test design
  design <- Experience$X
  colnames(design) <- row.names(P)
  write.csv(design, "MorrisExpDesign.csv", row.names = FALSE)
  
  pdf("Distr_ParFull.pdf")
  par(mfrow=c(4,3))
  for( i in 1:ncol(design)){ hist(design[,i], breaks=100, col="purple", main=colnames(design)[i] ) }
  dev.off()
  save(Experience, file = here(OUTPUT, "ExperienceFull.RData"))
  
  # load(here(OUTPUT,"ExperienceFull.RData"))
  
  # RUN MORRIS DESIGN ####
  # ---------------------------------------------------------------------------- #

  # Runs r (n of repetitions) * (param+1) simulations * (n) model outputs * (n) model outputs per time
  results <- apply(design, 1,  SENSI_model) 
  
  y <- results
  y <- as.data.frame(y) 
  save(y, file = here(OUTPUT, "y.RData"))
  
  # load(here(OUTPUT, "y.RData"))
  
  # CALCULATE MORRIS INDICES ####
  # ---------------------------------------------------------------------------- #
  
  sim.results.morris <- as.matrix(y)
  
  ZEN.sim.results.morris <- sim.results.morris[str_detect(rownames(sim.results.morris), "ZEN") &
                                              !str_detect(rownames(sim.results.morris), "ZEN_gluc"), ]
  
  aZEL.sim.results.morris <- sim.results.morris[str_detect(rownames(sim.results.morris), "aZEL") &
                                                    str_detect(rownames(sim.results.morris), "aZEL_gluc"), ]
  
  ZEN_gluc.sim.results.morris <- sim.results.morris[str_detect(rownames(sim.results.morris), "ZEN_gluc"), ]
  
  aZEL_gluc.sim.results.morris <- sim.results.morris[str_detect(rownames(sim.results.morris), "aZEL_gluc"), ]
  
  
  
  mu = mu.star = sigma = NULL

  # Tell morris y for each model output
  for (i in (1:length(sim.results.morris[,1])))
  {
    tell(Experience, sim.results.morris[i,])
    mu = rbind( mu, apply(Experience$ee, 2, mean, na.rm = TRUE))
    mu.star = rbind( mu.star, apply(Experience$ee, 2, function(x) mean(abs(x), na.rm = TRUE)))
    sigma = rbind( sigma, apply(Experience$ee, 2, sd, na.rm = TRUE))

  }
  
  
  plot(Experience)
  
  # Save mu.star, mu, sigma, scaled indices and global indices
  rownames( mu.star ) = rownames( mu ) = rownames( sigma ) = row.names(sim.results.morris)
  
  mu.star_scale = mu.star/apply(mu.star, 1, mean, na.rm=TRUE) 
  sigma_scale = sigma/apply(sigma, 1, mean, na.rm=TRUE)
  
  # Calculate global indices
  IGlobal_mu.star = apply(mu.star_scale, 2, mean, na.rm=TRUE)
  IGlobal_sigma = apply(sigma_scale, 2, mean, na.rm = TRUE )
  
  Morris_IGlobal.df <- data.frame(
    IGlobal_mu.star = apply(mu.star_scale, 2, mean, na.rm=TRUE),
    IGlobal_sigma = apply(sigma_scale, 2, mean, na.rm = TRUE )
  )
  
  IGlobal <- data.frame(IGlobal_sigma, IGlobal_mu.star) 
  IGlobal <- IGlobal %>% mutate(parm_names = rownames(IGlobal))

  Morris.indices <- list(
    mu = mu,
    mu_star = mu.star,
    sigma = sigma,
    mu_star_scale = mu.star_scale,
    sigma_scale = sigma_scale,
    IGlobal = IGlobal
  )
  
  save(Morris.indices, file = here(OUTPUT,"Morris.indices.RData"))
  
  
  ## ZEN ----
  
  mu = mu.star = sigma = NULL
  
  # Tell morris y for each model output
  for (i in (1:length(ZEN.sim.results.morris[,1])))
  {
    tell(Experience, ZEN.sim.results.morris[i,])
    mu = rbind( mu, apply(Experience$ee, 2, mean, na.rm = TRUE))
    mu.star = rbind( mu.star, apply(Experience$ee, 2, function(x) mean(abs(x), na.rm = TRUE)))
    sigma = rbind( sigma, apply(Experience$ee, 2, sd, na.rm = TRUE))
    
  }
  plot(Experience)
  
  # Save mu.star, mu, sigma, scaled indices and global indices
  rownames( mu.star ) = rownames( mu ) = rownames( sigma ) = row.names(ZEN.sim.results.morris)
  
  mu.star_scale = mu.star/apply(mu.star, 1, mean, na.rm=TRUE) 
  sigma_scale = sigma/apply(sigma, 1, mean, na.rm=TRUE)
  
  # Calculate global indices
  IGlobal_mu.star = apply(mu.star_scale, 2, mean, na.rm=TRUE)
  IGlobal_sigma = apply(sigma_scale, 2, mean, na.rm = TRUE )
  
  Morris_IGlobal.df <- data.frame(
    IGlobal_mu.star = apply(mu.star_scale, 2, mean, na.rm=TRUE),
    IGlobal_sigma = apply(sigma_scale, 2, mean, na.rm = TRUE )
  )
  
  IGlobal <- data.frame(IGlobal_sigma, IGlobal_mu.star) 
  IGlobal <- IGlobal %>% mutate(parm_names = rownames(IGlobal))
  
  Morris.indices <- list(
    mu = mu,
    mu_star = mu.star,
    sigma = sigma,
    mu_star_scale = mu.star_scale,
    sigma_scale = sigma_scale,
    IGlobal = IGlobal
  )
  
  save(Morris.indices, file = here(OUTPUT,"Morris_ZEN.indices.RData"))
  
  ## ZEN glucuronide ----
  
  mu = mu.star = sigma = NULL
  
  # Tell morris y for each model output
  for (i in (1:length(ZEN_gluc.sim.results.morris[,1])))
  {
    tell(Experience, ZEN_gluc.sim.results.morris[i,])
    mu = rbind( mu, apply(Experience$ee, 2, mean, na.rm = TRUE))
    mu.star = rbind( mu.star, apply(Experience$ee, 2, function(x) mean(abs(x), na.rm = TRUE)))
    sigma = rbind( sigma, apply(Experience$ee, 2, sd, na.rm = TRUE))
    
  }
  plot(Experience)
  
  # Save mu.star, mu, sigma, scaled indices and global indices
  rownames( mu.star ) = rownames( mu ) = rownames( sigma ) = row.names(ZEN_gluc.sim.results.morris)
  
  mu.star_scale = mu.star/apply(mu.star, 1, mean, na.rm=TRUE) 
  sigma_scale = sigma/apply(sigma, 1, mean, na.rm=TRUE)
  
  # Calculate global indices
  IGlobal_mu.star = apply(mu.star_scale, 2, mean, na.rm=TRUE)
  IGlobal_sigma = apply(sigma_scale, 2, mean, na.rm = TRUE )
  
  Morris_IGlobal.df <- data.frame(
    IGlobal_mu.star = apply(mu.star_scale, 2, mean, na.rm=TRUE),
    IGlobal_sigma = apply(sigma_scale, 2, mean, na.rm = TRUE )
  )
  
  IGlobal <- data.frame(IGlobal_sigma, IGlobal_mu.star) 
  IGlobal <- IGlobal %>% mutate(parm_names = rownames(IGlobal))
  
  Morris.indices <- list(
    mu = mu,
    mu_star = mu.star,
    sigma = sigma,
    mu_star_scale = mu.star_scale,
    sigma_scale = sigma_scale,
    IGlobal = IGlobal
  )
  
  save(Morris.indices, file = here(OUTPUT,"Morris_ZEN_gluc.indices.RData"))
  
  ## aZEL ----
  
  mu = mu.star = sigma = NULL
  
  # Tell morris y for each model output
  for (i in (1:length(aZEL.sim.results.morris[,1])))
  {
    tell(Experience, aZEL.sim.results.morris[i,])
    mu = rbind( mu, apply(Experience$ee, 2, mean, na.rm = TRUE))
    mu.star = rbind( mu.star, apply(Experience$ee, 2, function(x) mean(abs(x), na.rm = TRUE)))
    sigma = rbind( sigma, apply(Experience$ee, 2, sd, na.rm = TRUE))
    
  }
  plot(Experience)
  
  # Save mu.star, mu, sigma, scaled indices and global indices
  rownames( mu.star ) = rownames( mu ) = rownames( sigma ) = row.names(aZEL.sim.results.morris)
  
  mu.star_scale = mu.star/apply(mu.star, 1, mean, na.rm=TRUE) 
  sigma_scale = sigma/apply(sigma, 1, mean, na.rm=TRUE)
  
  # Calculate global indices
  IGlobal_mu.star = apply(mu.star_scale, 2, mean, na.rm=TRUE)
  IGlobal_sigma = apply(sigma_scale, 2, mean, na.rm = TRUE )
  
  Morris_IGlobal.df <- data.frame(
    IGlobal_mu.star = apply(mu.star_scale, 2, mean, na.rm=TRUE),
    IGlobal_sigma = apply(sigma_scale, 2, mean, na.rm = TRUE )
  )
  
  IGlobal <- data.frame(IGlobal_sigma, IGlobal_mu.star) 
  IGlobal <- IGlobal %>% mutate(parm_names = rownames(IGlobal))
  
  Morris.indices <- list(
    mu = mu,
    mu_star = mu.star,
    sigma = sigma,
    mu_star_scale = mu.star_scale,
    sigma_scale = sigma_scale,
    IGlobal = IGlobal
  )
  
  save(Morris.indices, file = here(OUTPUT,"Morris_aZEL.indices.RData"))
  
  ## aZEL glucuronide ----
  
  mu = mu.star = sigma = NULL
  
  # Tell morris y for each model output
  for (i in (1:length(aZEL_gluc.sim.results.morris[,1])))
  {
    tell(Experience, aZEL_gluc.sim.results.morris[i,])
    mu = rbind( mu, apply(Experience$ee, 2, mean, na.rm = TRUE))
    mu.star = rbind( mu.star, apply(Experience$ee, 2, function(x) mean(abs(x), na.rm = TRUE)))
    sigma = rbind( sigma, apply(Experience$ee, 2, sd, na.rm = TRUE))
    
  }
  plot(Experience)
  
  # Save mu.star, mu, sigma, scaled indices and global indices
  rownames( mu.star ) = rownames( mu ) = rownames( sigma ) = row.names(aZEL_gluc.sim.results.morris)
  
  mu.star_scale = mu.star/apply(mu.star, 1, mean, na.rm=TRUE) 
  sigma_scale = sigma/apply(sigma, 1, mean, na.rm=TRUE)
  
  # Calculate global indices
  IGlobal_mu.star = apply(mu.star_scale, 2, mean, na.rm=TRUE)
  IGlobal_sigma = apply(sigma_scale, 2, mean, na.rm = TRUE )
  
  Morris_IGlobal.df <- data.frame(
    IGlobal_mu.star = apply(mu.star_scale, 2, mean, na.rm=TRUE),
    IGlobal_sigma = apply(sigma_scale, 2, mean, na.rm = TRUE )
  )
  
  IGlobal <- data.frame(IGlobal_sigma, IGlobal_mu.star) 
  IGlobal <- IGlobal %>% mutate(parm_names = rownames(IGlobal))
  
  Morris.indices <- list(
    mu = mu,
    mu_star = mu.star,
    sigma = sigma,
    mu_star_scale = mu.star_scale,
    sigma_scale = sigma_scale,
    IGlobal = IGlobal
  )
  
  save(Morris.indices, file = here(OUTPUT,"Morris_aZEL_gluc.indices.RData"))
  