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
  
  # For plotting
  library(showtext)
  font_add(family = "Garamond", regular = "GARA.TTF")
  showtext_auto()
  CP_theme <- theme_minimal() +
    theme(
      axis.text = element_text(size = 56),
      axis.title = element_text(size = 62),
      plot.title = element_text(size = 65, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 56, hjust = 0.5),
      legend.position = "none",
      plot.margin = margin(0.3, 0.3, 0.3, 0.3, "cm") 
    )
  
  # Set output storage directory
  # OUTPUT <- here("Output", format(Sys.Date(), "%Y-%m-%d"), format(Sys.time(), "%H-%M-%S"))
  # dir.create(OUTPUT, recursive = TRUE)

  # Model Code
  source(here("Script", "PBK_model_rat.R"))
  

  # SENSITIVITY FUNCTION ####
  # ---------------------------------------------------------------------------- #
  
  ## Morris test ####
  SENSI_model <- function(P){
    
    
    
    A_init <- c(
      AOD_ZEN = 0, 
      
      ADuoL_ZEN = 0, AJejL_ZEN = 0, AIleL_ZEN = 0, ACOL_ZEN = 0,
      ADuo_ZEN = 0, AJej_ZEN = 0, AIle_ZEN = 0, 
      AL_ZEN = 0, 
      AR_ZEN = 0,
      AS_ZEN = 0,
      AF_ZEN = 0, 
      AB_ZEN = 0, 
      AUR_ZEN = 0, 
      AFE_ZEN = 0, 
      
      A_bZEL = 0,
      
      ADuoL_ZEN_gluc = 0, AJejL_ZEN_gluc = 0, AIleL_ZEN_gluc = 0, ACOL_ZEN_gluc = 0,
      AJej_IC_ZEN_gluc = 0, AJej_EC_ZEN_gluc = 0, 
      AL_IC_ZEN_gluc = 0, AL_EC_ZEN_gluc = 0,
      ABile_ZEN_gluc = 0, 
      AR_ZEN_gluc = 0,
      AS_ZEN_gluc = 0,
      AF_ZEN_gluc = 0, 
      AB_ZEN_gluc = 0, 
      AUR_ZEN_gluc = 0, 
      AFE_ZEN_gluc = 0,
      
      AIleL_aZEL = 0, ACOL_aZEL = 0,
      AJej_aZEL = 0, AIle_aZEL = 0,
      AL_aZEL = 0, 
      AR_aZEL = 0,
      AS_aZEL = 0,
      AF_aZEL = 0,
      AB_aZEL = 0,
      AUR_aZEL = 0,
      AFE_aZEL = 0,
      
      ADuoL_aZEL_gluc = 0, AJejL_aZEL_gluc = 0, AIleL_aZEL_gluc = 0, ACOL_aZEL_gluc = 0,
      AJej_IC_aZEL_gluc = 0, AJej_EC_aZEL_gluc = 0,
      AL_IC_aZEL_gluc = 0, AL_EC_aZEL_gluc = 0,
      ABile_aZEL_gluc = 0, 
      AR_aZEL_gluc = 0, 
      AS_aZEL_gluc = 0,
      AF_aZEL_gluc = 0,
      AB_aZEL_gluc = 0,
      AUR_aZEL_gluc = 0,
      AFE_aZEL_gluc = 0
      
    )
    
    
    PBK.out <- lsoda(y = A_init,
                     times = seq(0, 24,by=1/10),
                     func =  R.PBK.model,
                     parms = P)
    
    CB_ZEN <- PBK.out[,"CB_ZEN"]
    CB_ZEN_gluc <- PBK.out[,"CB_ZEN_gluc"]
    CB_aZEL <- PBK.out[,"CB_aZEL"]
    CB_aZEL_gluc <- PBK.out[,"CB_aZEL_gluc"]
    
    return( c(CB_ZEN=CB_ZEN, 
              CB_ZEN_gluc = CB_ZEN_gluc, 
              CB_aZEL = CB_aZEL, 
              CB_aZEL_gluc = CB_aZEL_gluc))
    
  }
  
  # SET MORRIS PARAMETERS ####
  # ---------------------------------------------------------------------------- #
  
  # Parameter upper and lower bounds
  
  # Parameters
  Parameters_df <- read.csv(here("Input", "Parameters_rat.csv"))
  # Parameters <- as.list(Parameters_df$Value)
  # names(Parameters) <- Parameters_df$Parameter

  # Dosing
  PhysChem <- read_delim("Input/PhysChem.csv", 
                         delim = ";", escape_double = FALSE, trim_ws = TRUE)
 
  # oral dose
  MW_ZEN <- pull(PhysChem[PhysChem$Substance == "Zearalenone", 
  "MolecularWeight"])
  BW_rat <- 0.26 #kg
  ODOSEmg <- 1 #{mg/kg bw} # oral dose, variable
  ODOSEumol <- ODOSEmg*1000/MW_ZEN*0.26 # (μmol) unit change to μmol
  
  Dosing <- c(Paramter = "ODOSEumol",
              Value = ODOSEumol)
  ParametersMorris <- rbind(Parameters_df, Dosing)
  
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
  Morris.r <- 3 #50
  
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
  save(Experience, file = "ExperienceFull.RData")
  
  load(here("ExperienceFull.RData"))
  
  # RUN MORRIS DESIGN ####
  # ---------------------------------------------------------------------------- #

  # Runs r (n of repetitions) * (param+1) simulations * (n) model outputs * (n) model outputs per time
  results <- apply(design, 1,  SENSI_model) 
  
  y <- results
  y <- as.data.frame(y) #%>% mutate(Tout = "times = seq(0, 5*365,by=1/10)")
  save(y, file = "y.RData")
  
  load(here("y.RData"))
  
  # CALCULATE MORRIS INDICES ####
  # ---------------------------------------------------------------------------- #
  
  sim.results.morris <- as.matrix(y)
  
  CB_ZEN.sim.results.morris <- sim.results.morris[str_detect(rownames(sim.results.morris), "CB_ZEN") &
                                              !str_detect(rownames(sim.results.morris), "CB_ZEN_gluc"), ]
  
  CB_aZEL.sim.results.morris <- sim.results.morris[str_detect(rownames(sim.results.morris), "CB_aZEL") &
                                                    str_detect(rownames(sim.results.morris), "CB_aZEL_gluc"), ]
  
  CB_ZEN_gluc.sim.results.morris <- sim.results.morris[str_detect(rownames(sim.results.morris), "CB_ZEN_gluc"), ]
  
  CB_aZEL_gluc.sim.results.morris <- sim.results.morris[str_detect(rownames(sim.results.morris), "CB_aZEL_gluc"), ]
  
  
  
  mu = mu.star = sigma = NULL

  # Tell morris y for each model output
  for (i in (1:length(CB_ZEN.sim.results.morris[,1])))
  {
    tell(Experience, CB_ZEN.sim.results.morris[i,])
    mu = rbind( mu, apply(Experience$ee, 2, mean, na.rm = TRUE))
    mu.star = rbind( mu.star, apply(Experience$ee, 2, function(x) mean(abs(x), na.rm = TRUE)))
    sigma = rbind( sigma, apply(Experience$ee, 2, sd, na.rm = TRUE))

  }
  plot(Experience)
  
  # Save mu.star, mu, sigma, scaled indices and global indices
  rownames( mu.star ) = rownames( mu ) = rownames( sigma ) = row.names(CB_ZEN.sim.results.morris)
  
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
  
  save(Morris.indices, file = here("CB_ZEN.Morris.indices.RData"))
  
  
  # PLOT MORRIS RESULTS ####
  # ---------------------------------------------------------------------------- #
  
  ggplot(IGlobal) +
    aes(x = IGlobal_mu.star, y = IGlobal_sigma, label = rownames(IGlobal)) +
    # xlim (min,2.5) + #to zoom in/out, determine max value
    # ylim (min,2.5) + #to zoom in/out, determine max value
    geom_point()
