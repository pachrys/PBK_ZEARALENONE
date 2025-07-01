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
  source(here("Script", "PBK_code_Mendez.R"))
  

  # SENSITIVITY FUNCTION ####
  # ---------------------------------------------------------------------------- #
  
  ## Morris test ####
  SENSI_model <- function(parm.c){
    
    
    A_init <- c(
      Ast = 0,
      Ain1 = 0,
      Ain2 = 0,
      Ain3 = 0,
      Ain4 = 0,
      Ain5 = 0,
      Ain6 = 0,
      Ain7 = 0,
      ASIZEN = 0,
      ASIZENGLU = 0,
      ALIZEN = 0,
      ALIaZEL = 0,
      ALIbZEL = 0,
      AfeZ = 0,
      AfeaZ = 0,
      ALZEN = 0,
      ALaZEL = 0,
      ALbZEL = 0,
      ALZENGLU = 0,
      ALaZELGLU = 0,
      AF = 0,
      AR = 0,
      AS = 0,
      AB = 0,
      AZur = 0,
      AFaZEL = 0,
      ARaZEL = 0,
      ASaZEL = 0,
      ABaZEL = 0,
      AaZur = 0,
      ATotal = 0
    )
    
    
    PBK.out <- lsoda(y = A_init,
                     times = seq(0, 24,by=1/10),
                     func =  PBK.model,
                     parms = parm.c)
    
    CBaZEL <- PBK.out[,"CBaZEL"]
    CB <- PBK.out[,"CB"]
    CLZEN <- PBK.out[,"CLZEN"]
    CLZENGLU <- PBK.out[,"CLZENGLU"]
    
    return( c(CBaZEL=CBaZEL, 
           CB = CB, 
           CLZEN = CLZEN, 
           CLZENGLU = CLZENGLU))
    
  }
  
  # SET MORRIS PARAMETERS ####
  # ---------------------------------------------------------------------------- #
  
  # Parameter upper and lower bounds
  
  # Parameters
  Parameters <- read.csv(here("Input", "Parm.rat.csv"))
  
  # Dosing
  # molecular weight
  MWZEN <- 318.37 # molecular weight ZEN
  
  # oral dose
  ODOSEmg <- 0.5 #{mg/kg bw} # oral dose, variable
  ODOSEumol <- ODOSEmg*1000/MWZEN*Parameters[1, "Value"] # (μmol) unit change to μmol
  
  # IV dose
  IVDOSEmg1 <- 0.5 #{mg/kg bw} # IV dose, variable
  IVDOSEumol <-IVDOSEmg1*1000/MWZEN*Parameters[1, "Value"] # unit change to μmol
  # time
  
  Dosing <- data.frame(
    Parameter = c("ODOSEumol", "IVDOSEumol"),
    Value = c(ODOSEumol, IVDOSEumol)
  )
  
  
  ParametersMorris <- rbind(Parameters, Dosing)
  P <- ParametersMorris$Value
  P <- ParametersMorris %>% tibble::column_to_rownames("Parameter")
  
  binf <- P$Value*0.9
  bsup <- P$Value*1.1
  
  Morris.factors <- ParametersMorris$Parameter  
  N.par <- nrow(ParametersMorris)
  
  
  # SET MORRIS EXPERIMENT DESIGN ####
  # ---------------------------------------------------------------------------- #
  
  
  # Experiment design
  Morris.design <- list(type = "oat", levels = 6, grid.jump = 3)
  
  # Repetitions
  Morris.r <- 50
  
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
  rownames( mu.star ) = rownames( mu ) = rownames( sigma ) = row.names(y)
  
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
  
  save(Morris.indices, file = "Morris.indices.RData")
  
  
  # PLOT MORRIS RESULTS ####
  # ---------------------------------------------------------------------------- #
  
  ggplot(IGlobal) +
    aes(x = IGlobal_mu.star, y = IGlobal_sigma, label = rownames(IGlobal)) +
    # xlim (min,2.5) + #to zoom in/out, determine max value
    # ylim (min,2.5) + #to zoom in/out, determine max value
    geom_point()

  
  morris_plot <- IGlobal %>% 
    # filter(mu.star_scale >= 0.01) %>%
    ggplot(aes(x = IGlobal_mu.star, y = IGlobal_sigma, label = parm_names)) +
    geom_point(alpha = 0.5, color = "darkblue") +
    geom_text_repel(
      size = 7,
      angle = 30,
      max.overlaps = Inf,
      segment.color = "darkblue",
      segment.size = 0.1,
      segment.alpha = 0.5,
      segment.linetype = "dotted",
      direction = "x",
      ylim = c(7, 7),
      point.padding = 0.5,
      box.padding = 0.2,
      set.seed(123)
    ) +
    coord_cartesian(clip = "off") +
    labs(
      title = "Morris Sensitivity Analysis",
      x = "μ* (Mean Effect)",
      y = "σ (Standard Deviation)"
      # subtitle = "Parameters with relative μ* ≥ 5%"
    ) +
    # theme_minimal() +
    CP_theme+
    theme(plot.margin = margin(0, 0.15, 0, 0, "cm"))
  morris_plot
  ggsave(
    filename = here("morrisGI_plot.png"),
    plot = morris_plot,
    dpi = 1000, 
    width = 12, height = 9, units = "cm"
  )

  #Zoom out
  morris_plot <- IGlobal %>% 
    filter(IGlobal_mu.star >= 0.15) %>%
    ggplot(aes(x = IGlobal_mu.star, y = IGlobal_sigma, label = parm_names)) +
    geom_point(alpha = 0.5, color = "darkblue") +
    geom_text_repel(
      size = 9,
      angle = 30,
      max.overlaps = Inf,
      segment.color = "darkblue",
      segment.size = 0.1,
      segment.alpha = 0.5,
      segment.linetype = "dotted",
      direction = "x",
      ylim = c(7, 7),
      point.padding = 0.4,
      box.padding = 0.1,
      set.seed(123)
    ) +
    coord_cartesian(clip = "off") +
    labs(
      title = "Morris Sensitivity Analysis",
      x = "μ* (Mean Effect)",
      y = "σ (Standard Deviation)",
      subtitle = "Parameters with μ* ≥ 0.15"
    ) +
    CP_theme+
    # theme_minimal() +
    theme(plot.margin = margin(0, 0.15, 0, 0, "cm"))
  morris_plot
  ggsave(
    filename = here("morrisZoomOut_plot.png"),
    plot = morris_plot,
    dpi = 1000, 
    width = 12, height = 9, units = "cm"
  )
  
  #Zoom in
  morris_plot <- IGlobal %>% 
    filter(IGlobal_mu.star <=0.15) %>%
    ggplot(aes(x = IGlobal_mu.star, y = IGlobal_sigma, label = parm_names)) +
    geom_point(alpha = 0.5, color = "darkblue") +
    geom_text_repel(
      size = 9,
      angle = 30,
      max.overlaps = Inf,
      segment.color = "darkblue",
      segment.size = 0.1,
      segment.alpha = 0.5,
      segment.linetype = "dotted",
      direction = "x",
      ylim = c(0.1, 0.1),
      point.padding = 0.4,
      box.padding = 0.3,
      set.seed(123)
    ) +
    coord_cartesian(clip = "off") +
    labs(
      title = "Morris Sensitivity Analysis",
      x = "μ* (Mean Effect)",
      y = "σ (Standard Deviation)",
      subtitle = "Parameters with μ* <= 0.15"
    ) +
    CP_theme+
    # theme_minimal() +
    theme(plot.margin = margin(0, 0.15, 0, 0, "cm"))
  morris_plot
  ggsave(
    filename = here("morrisZoomIn_plot.png"),
    plot = morris_plot,
    dpi = 1000, 
    width = 12, height = 9, units = "cm"
  )
  # D <- ggplotly(Plot)
  