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
  Parameters_df <- read.csv(here("Input", "Parameters_human.csv")) #Parameters_rat.csv Parameters_human.csv
  # Parameters <- as.list(Parameters_df$Value)
  # names(Parameters) <- Parameters_df$Parameter

  # Dosing
  PhysChem <- read_delim("Input/PhysChem.csv", 
                         delim = ";", escape_double = FALSE, trim_ws = TRUE)
 
  # oral dose
  # MW_ZEN <- pull(PhysChem[PhysChem$Substance == "Zearalenone", 
  # "MolecularWeight"])
  # BW_rat <- 0.26 #kg
  # ODOSEmg <- 1 #{mg/kg bw} # oral dose, variable
  # ODOSEumol <- ODOSEmg*1000/MW_ZEN*0.26 # (μmol) unit change to μmol
  # 
  # Dosing <- c(Paramter = "ODOSEumol",
  #             Value = ODOSEumol)
  ParametersMorris <- Parameters_df #rbind(Parameters_df, Dosing)
  
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
  y <- as.data.frame(y) #%>% mutate(Tout = "times = seq(0, 5*365,by=1/10)")
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
  
  
#   # Plotting ----
#   # 
#   # load(here(OUTPUT,"Morris_ZEN.indices.RData"))
#   # morris_data <- Morris.indices$IGlobal
#   # ALL <- morris_data %>% arrange(desc(IGlobal_mu.star)) %>% 
#   #   mutate(n = seq(from=1, to=113, by=1)) %>% 
#   #   filter(IGlobal_sigma > 0)
#   # 
#   # ALL <- ALL %>%  mutate(
#   #   "ParameterType" = case_when(
#   #     str_detect(parm_names, regex("MRP2|BCRP")) ~ "Glucuronide excretion",
#   #     str_detect(parm_names, regex("MRP3")) ~ "Glucuronide absorption",
#   #     str_detect(parm_names, regex("Hydro")) ~ "Hydrolysis",
#   #     str_detect(parm_names, regex("L_ZEN_gluc|SI_ZEN_gluc|L_aZEL_gluc|SI_aZEL_gluc")) ~ "Glucuronidation",
#   #     str_detect(parm_names, regex("L_bZEL|SI_bZEL|L_aZEL|SI_aZEL")) ~ "Reduction",
#   #     str_detect(parm_names, regex("PSI|PL|PR|PS|PR_gluc|PS_gluc|PF|PK|BPR")) ~ "Partitioning",
#   #     str_detect(parm_names, regex("ODOSEumol")) ~ "Dose",
#   #     str_detect(parm_names, regex("Papp")) ~ "Aglycone absorption",
#   #     str_detect(parm_names, regex("fup|fut")) ~ "Fraction unbound",
#   #     TRUE ~ "Physiological"))      
#   # 
#   # colouring <- c(
#   #   "Glucuronide excretion" = "#FF8C00",
#   #   "Glucuronide absorption" = "#FF6D8D",
#   #   "Hydrolysis" = "#1B60FF",
#   #   "Glucuronidation" = "#1B85FF",
#   #   "Reduction" = "#1B35FF",
#   #   "Aglycone absorption" = "#FF3D8D",
#   #   "Partitioning" = "#8631A7",
#   #   "Fraction unbound" = "#4631A7",
#   #   "Physiological" = "#008475",
#   #   "Dose" = "black")
#   # 
#   # set.seed(1234)
#   # ALL_morris_plot <- ALL %>% filter(n <= 20) %>% 
#   #   ggplot(aes(x = IGlobal_mu.star, y = IGlobal_sigma, 
#   #              label = parm_names, colour = ParameterType)) +
#   #   geom_jitter(width = .5, alpha = 0.7, size = 1)  +
#   #   geom_text_repel(aes(label = as.character(parm_names)), 
#   #                   size=15, 
#   #                   max.overlaps = Inf,
#   #                   # min.segment.length = 0.05,
#   #                   segment.size = 0.1
#   #   ) +
#   #   scale_color_manual(name = "Parameter category",
#   #                      values = colouring, 
#   #                      breaks = c(
#   #                        "Dose" = "Dose",
#   #                        "Physiological" = "Physiological",
#   #                        "Partitioning" = "Partitioning",
#   #                        "Fraction unbound" = "Fraction unbound",
#   #                        "Reduction" = "Reduction",
#   #                        "Glucuronidation" = "Glucuronidation",
#   #                        "Hydrolysis" = "Hydrolysis",
#   #                        "Aglycone absorption" = "Aglycone absorption",
#   #                        "Glucuronide absorption" = "Glucuronide absorption",
#   #                        "Glucuronide excretion" = "Glucuronide excretion"
#   #                      )) +
#   #   labs(
#   #     # title = "A",
#   #     x = expression(paste("Mean Effect (μ*)")),
#   #     y = expression(paste("Standard Deviation (σ)"))
#   #   ) +
#   #   CP_theme +
#   #   theme(legend.position = "none")
#   # 
#   # ggsave(
#   #   filename = here(OUTPUT,"Morris_ALL.png"),
#   #   plot = ALL_morris_plot,
#   #   dpi = 600,  
#   #   width = 24, 
#   #   height = 18, 
#   #   units = "cm",
#   #   limitsize = FALSE
#   # )
#   # 
#   load(here(OUTPUT,"Morris_ZEN.indices.RData"))
#   ZEN_morris_data <- Morris.indices$IGlobal
#   ZEN <- ZEN_morris_data %>% arrange(desc(IGlobal_mu.star)) %>% 
#     mutate(n = seq(from=1, to=113, by=1)) %>% 
#     filter(IGlobal_sigma > 0)
#   write.csv(ZEN, here(OUTPUT,"ZEN_Morris_results.csv"))
#   ZEN <- ZEN %>%  mutate(
#     "ParameterType" = case_when(
#       str_detect(parm_names, regex("MRP2|BCRP")) ~ "Glucuronide excretion",
#       str_detect(parm_names, regex("MRP3")) ~ "Glucuronide absorption",
#       str_detect(parm_names, regex("Hydro")) ~ "Hydrolysis",
#       str_detect(parm_names, regex("L_ZEN_gluc|SI_ZEN_gluc|L_aZEL_gluc|SI_aZEL_gluc")) ~ "Glucuronidation",
#       str_detect(parm_names, regex("L_bZEL|SI_bZEL|L_aZEL|SI_aZEL")) ~ "Reduction",
#       str_detect(parm_names, regex("PSI|PL|PR|PS|PR_gluc|PS_gluc|PF|PK|BPR")) ~ "Partitioning",
#       str_detect(parm_names, regex("ODOSEumol")) ~ "Dose",
#       str_detect(parm_names, regex("Papp")) ~ "Aglycone absorption",
#       str_detect(parm_names, regex("fup|fut")) ~ "Fraction unbound",
#       TRUE ~ "Physiological"))      
#   
#   colouring <- c(
#     "Glucuronide excretion" = "#FF8C00",
#     "Glucuronide absorption" = "#FF6D8D",
#     "Hydrolysis" = "#1B60FF",
#     "Glucuronidation" = "#1B85FF",
#     "Reduction" = "#1B35FF",
#     "Aglycone absorption" = "#FF3D8D",
#     "Partitioning" = "#8631A7",
#     "Fraction unbound" = "#4631A7",
#     "Physiological" = "#008475",
#     "Dose" = "black")
#   
#   set.seed(1234)
#   ZEN_morris_plot <- ZEN %>% filter(n <= 20) %>% 
#   ggplot(aes(x = IGlobal_mu.star, y = IGlobal_sigma, 
#              label = parm_names, colour = ParameterType)) +
#     geom_jitter(width = .5, alpha = 0.7, size = 1)  +
#     geom_text_repel(aes(label = as.character(parm_names)), 
#                     size=15, 
#                     max.overlaps = Inf,
#                     # min.segment.length = 0.05,
#                     segment.size = 0.1
#                     ) +
#     scale_color_manual(name = "Parameter category",
#                        values = colouring, 
#                        breaks = c(
#                          "Dose" = "Dose",
#                          "Physiological" = "Physiological",
#                          "Partitioning" = "Partitioning",
#                          "Fraction unbound" = "Fraction unbound",
#                          "Reduction" = "Reduction",
#                          "Glucuronidation" = "Glucuronidation",
#                          "Hydrolysis" = "Hydrolysis",
#                          "Aglycone absorption" = "Aglycone absorption",
#                          "Glucuronide absorption" = "Glucuronide absorption",
#                          "Glucuronide excretion" = "Glucuronide excretion"
#                        )) +
#     labs(
#       title = "A",
#       x = expression(paste("Mean Effect (μ*)")),
#       y = expression(paste("Standard Deviation (σ)"))
#     ) +
#     CP_theme +
#     theme(legend.position = "none")
#     
#   
# 
#   load(here(OUTPUT,"Morris_ZEN_gluc.indices.RData"))
#   ZEN_glucmorris_data <- Morris.indices$IGlobal
#   ZEN_gluc <- ZEN_glucmorris_data %>% arrange(desc(IGlobal_mu.star)) %>% 
#     mutate(n = seq(from=1, to=113, by=1)) %>% 
#     filter(IGlobal_sigma > 0)
#   ZEN_gluc <- ZEN_gluc %>%  mutate(
#     "ParameterType" = case_when(
#       str_detect(parm_names, regex("MRP2|BCRP")) ~ "Glucuronide excretion",
#       str_detect(parm_names, regex("MRP3")) ~ "Glucuronide absorption",
#       str_detect(parm_names, regex("Hydro")) ~ "Hydrolysis",
#       str_detect(parm_names, regex("L_ZEN_gluc|SI_ZEN_gluc|L_aZEL_gluc|SI_aZEL_gluc")) ~ "Glucuronidation",
#       str_detect(parm_names, regex("L_bZEL|SI_bZEL|L_aZEL|SI_aZEL")) ~ "Reduction",
#       str_detect(parm_names, regex("PSI|PL|PR|PS|PR_gluc|PS_gluc|PF|PK|BPR")) ~ "Partitioning",
#       str_detect(parm_names, regex("ODOSEumol")) ~ "Dose",
#       str_detect(parm_names, regex("Papp")) ~ "Aglycone absorption",
#       str_detect(parm_names, regex("fup|fut")) ~ "Fraction unbound",
#       TRUE ~ "Physiological"))      
#   write.csv(ZEN_gluc, here(OUTPUT,"ZEN_gluc_Morris_results.csv"))
#   
#   ZEN_gluc_morris_plot <- ZEN_gluc %>% filter(n <= 20) %>% 
#     ggplot(aes(x = IGlobal_mu.star, y = IGlobal_sigma, 
#                label = parm_names, colour = ParameterType)) +
#     geom_jitter(width = .5, alpha = 0.7, size = 1)  +
#     geom_text_repel(aes(label = as.character(parm_names)), 
#                     size=15, 
#                     max.overlaps = Inf,
#                     # min.segment.length = 0.05,
#                     segment.size = 0.1
#     ) +
#     scale_color_manual(name = "Parameter category",
#                        values = colouring, 
#                        breaks = c(
#                          "Dose" = "Dose",
#                          "Physiological" = "Physiological",
#                          "Partitioning" = "Partitioning",
#                          "Fraction unbound" = "Fraction unbound",
#                          "Reduction" = "Reduction",
#                          "Glucuronidation" = "Glucuronidation",
#                          "Hydrolysis" = "Hydrolysis",
#                          "Aglycone absorption" = "Aglycone absorption",
#                          "Glucuronide absorption" = "Glucuronide absorption",
#                          "Glucuronide excretion" = "Glucuronide excretion"
#                        )) +
#     labs(
#       title = "B",
#       x = expression(paste("Mean Effect (μ*)")),
#       y = expression(paste("Standard Deviation (σ)"))
#     ) +
#     CP_theme +
#     theme(legend.position = "none")
#   
#   load(here(OUTPUT,"Morris_aZEL.indices.RData"))
#   aZEL_morris_data <- Morris.indices$IGlobal
#   aZEL <- aZEL_morris_data %>% arrange(desc(IGlobal_mu.star)) %>% 
#     mutate(n = seq(from=1, to=113, by=1)) %>% 
#     filter(IGlobal_sigma > 0)
#   aZEL <- aZEL %>%  mutate(
#     "ParameterType" = case_when(
#       str_detect(parm_names, regex("MRP2|BCRP")) ~ "Glucuronide excretion",
#       str_detect(parm_names, regex("MRP3")) ~ "Glucuronide absorption",
#       str_detect(parm_names, regex("Hydro")) ~ "Hydrolysis",
#       str_detect(parm_names, regex("L_ZEN_gluc|SI_ZEN_gluc|L_aZEL_gluc|SI_aZEL_gluc")) ~ "Glucuronidation",
#       str_detect(parm_names, regex("L_bZEL|SI_bZEL|L_aZEL|SI_aZEL")) ~ "Reduction",
#       str_detect(parm_names, regex("PSI|PL|PR|PS|PR_gluc|PS_gluc|PF|PK|BPR")) ~ "Partitioning",
#       str_detect(parm_names, regex("ODOSEumol")) ~ "Dose",
#       str_detect(parm_names, regex("Papp")) ~ "Aglycone absorption",
#       str_detect(parm_names, regex("fup|fut")) ~ "Fraction unbound",
#       TRUE ~ "Physiological"))      
#   write.csv(aZEL, here(OUTPUT,"aZEL_Morris_results.csv"))
#   
#   aZEL_morris_plot <- aZEL %>% filter(n <= 20) %>% 
#     ggplot(aes(x = IGlobal_mu.star, y = IGlobal_sigma, 
#                label = parm_names, colour = ParameterType)) +
#     geom_jitter(width = .5, alpha = 0.7, size = 1)  +
#     geom_text_repel(aes(label = as.character(parm_names)), 
#                     size=15, 
#                     max.overlaps = Inf,
#                     # min.segment.length = 0.05,
#                     segment.size = 0.1
#     ) +
#     scale_color_manual(name = "Parameter category",
#                        values = colouring, 
#                        breaks = c(
#                          "Dose" = "Dose",
#                          "Physiological" = "Physiological",
#                          "Partitioning" = "Partitioning",
#                          "Fraction unbound" = "Fraction unbound",
#                          "Reduction" = "Reduction",
#                          "Glucuronidation" = "Glucuronidation",
#                          "Hydrolysis" = "Hydrolysis",
#                          "Aglycone absorption" = "Aglycone absorption",
#                          "Glucuronide absorption" = "Glucuronide absorption",
#                          "Glucuronide excretion" = "Glucuronide excretion"
#                        )) +
#     labs(
#       title = "C",
#       x = expression(paste("Mean Effect (μ*)")),
#       y = expression(paste("Standard Deviation (σ)"))
#     ) +
#     CP_theme +
#     theme(legend.position = "none")
#   
#   
#   
#   
#   load(here(OUTPUT,"Morris_aZEL_gluc.indices.RData"))
#   aZEL_gluc_morris_data <- Morris.indices$IGlobal
#   aZEL_gluc <- aZEL_gluc_morris_data %>% arrange(desc(IGlobal_mu.star)) %>% 
#     mutate(n = seq(from=1, to=113, by=1)) %>% 
#     filter(IGlobal_sigma > 0)
#   
#   aZEL_gluc <- aZEL_gluc %>%  mutate(
#     "ParameterType" = case_when(
#       str_detect(parm_names, regex("MRP2|BCRP")) ~ "Glucuronide excretion",
#       str_detect(parm_names, regex("MRP3")) ~ "Glucuronide absorption",
#       str_detect(parm_names, regex("Hydro")) ~ "Hydrolysis",
#       str_detect(parm_names, regex("L_ZEN_gluc|SI_ZEN_gluc|L_aZEL_gluc|SI_aZEL_gluc")) ~ "Glucuronidation",
#       str_detect(parm_names, regex("L_bZEL|SI_bZEL|L_aZEL|SI_aZEL")) ~ "Reduction",
#       str_detect(parm_names, regex("PSI|PL|PR|PS|PR_gluc|PS_gluc|PF|PK|BPR")) ~ "Partitioning",
#       str_detect(parm_names, regex("ODOSEumol")) ~ "Dose",
#       str_detect(parm_names, regex("Papp")) ~ "Aglycone absorption",
#       str_detect(parm_names, regex("fup|fut")) ~ "Fraction unbound",
#       TRUE ~ "Physiological"))      
#   write.csv(aZEL_gluc, here(OUTPUT,"aZEL_gluc_Morris_results.csv"))
#   
#   aZEL_gluc_morris_plot <- aZEL_gluc %>% filter(n <= 20) %>% 
#     ggplot(aes(x = IGlobal_mu.star, y = IGlobal_sigma, 
#                label = parm_names, colour = ParameterType)) +
#     geom_jitter(width = .5, alpha = 0.7, size = 1)  +
#     geom_text_repel(aes(label = as.character(parm_names)), 
#                     size=15, 
#                     max.overlaps = Inf,
#                     # min.segment.length = 0.05,
#                     segment.size = 0.1
#     ) +
#     scale_color_manual(name = "Parameter category",
#                        values = colouring, 
#                        breaks = c(
#                          "Dose" = "Dose",
#                          "Physiological" = "Physiological",
#                          "Partitioning" = "Partitioning",
#                          "Fraction unbound" = "Fraction unbound",
#                          "Reduction" = "Reduction",
#                          "Glucuronidation" = "Glucuronidation",
#                          "Hydrolysis" = "Hydrolysis",
#                          "Aglycone absorption" = "Aglycone absorption",
#                          "Glucuronide absorption" = "Glucuronide absorption",
#                          "Glucuronide excretion" = "Glucuronide excretion"
#                        )) +
#     labs(
#       title = "D",
#       x = expression(paste("Mean Effect (μ*)")),
#       y = expression(paste("Standard Deviation (σ)"))
#     ) +
#     CP_theme +
#     theme(legend.position = "none")
#   
#   Morris_separate <- ZEN_morris_plot + 
#                      ZEN_gluc_morris_plot + 
#                      aZEL_morris_plot + 
#                      aZEL_gluc_morris_plot
#   ggsave(
#     filename = here(OUTPUT,"Morris_separate.png"),
#     plot = Morris_separate,
#     dpi = 600,  
#     width = 24, 
#     height = 18, 
#     units = "cm",
#     limitsize = FALSE
#   )
#   
#   ZEN_eFAST <- ZEN %>% filter(n <= 40) %>% mutate(Mycotoxin = "ZEN")
#   aZEL_eFAST <- aZEL %>% filter(n <= 40) %>% mutate(Mycotoxin = "aZEL")
#   ZEN_gluc_eFAST <- ZEN_gluc %>% filter(n <= 40) %>% mutate(Mycotoxin = "ZEN_gluc")
#   aZEL_gluc_eFAST <- aZEL_gluc %>% filter(n <= 40) %>% mutate(Mycotoxin = "aZEL_gluc")
#   eFASTparams <- rbind(ZEN_eFAST, aZEL_eFAST, ZEN_gluc_eFAST, aZEL_gluc_eFAST) 
#   
#   feFASTparams<- eFASTparams %>% 
#     group_by(parm_names)  %>% summarise(occurrences = n()) %>% 
#     arrange(-occurrences) 
#   
#   # PLOT MORRIS RESULTS ####
#   # ---------------------------------------------------------------------------- #
# 
#   ggplot(IGlobal, aes(x = IGlobal_mu.star, y = IGlobal_sigma, label = parm_names)) +
#     geom_jitter(width = .5, alpha = 0.7)  +
#     geom_text(aes(label=ifelse(IGlobal_mu.star>0.02,as.character(parm_names),'')),hjust=0,vjust=0) 
# 
#   Global_eFAST <- morris_data %>% filter(n<=40) %>% mutate(Mycotoxin = "Global")
# 
#   Occurence <- rbind(ZEN_eFAST, aZEL_eFAST, ZEN_gluc_eFAST, aZEL_gluc_eFAST) %>% 
#   group_by(parm_names)  %>% summarise(occurrences = n()) 
#   
#   Global_mu.star <- rbind(ZEN_eFAST, aZEL_eFAST, ZEN_gluc_eFAST, aZEL_gluc_eFAST) %>% 
#     group_by(parm_names)  %>% 
#     slice_max(IGlobal_mu.star, with_ties = FALSE)
#   
#   eFASTparamsN <- data.frame(
#     Occurence = Occurence$occurrences,
#     parm_names = Occurence$parm_names, 
#     parm_namesT = Global_mu.star$parm_names, 
#     parm_type = Global_mu.star$ParameterType,
#     IGlobal_mu.star = Global_mu.star$IGlobal_mu.star)
# 
#    toeFAST<-  eFASTparamsN %>% #mutate(parm_names = reorder(parm_names, +Occurence)) %>% 
#     ggplot()+
#     geom_point(aes(parm_names, IGlobal_mu.star, size = Occurence, colour = parm_type), alpha = 0.9) +
#     geom_segment(aes(parm_names, xend=parm_names, y=0, yend = IGlobal_mu.star)) +
#     scale_color_manual(name = "Parameter category",
#                        values = colouring, 
#                        breaks = c(
#                          "Dose" = "Dose",
#                          "Physiological" = "Physiological",
#                          "Partitioning" = "Partitioning",
#                          "Fraction unbound" = "Fraction unbound",
#                          "Reduction" = "Reduction",
#                          "Glucuronidation" = "Glucuronidation",
#                          "Hydrolysis" = "Hydrolysis",
#                          "Aglycone absorption" = "Aglycone absorption",
#                          "Glucuronide absorption" = "Glucuronide absorption",
#                          "Glucuronide excretion" = "Glucuronide excretion"
#                        )) +
#     coord_flip() +
#     labs(
#       x = "",
#       y = expression(paste("Mean Effect (μ*)"))) +
#     CP_theme +
#     theme(legend.spacing.y = unit(1, "cm"))
#   
#   ggsave(
#     filename = here(OUTPUT,"toeFAST.png"),
#     plot = toeFAST,
#     dpi = 600,  
#     width = 15, 
#     height = 18, 
#     units = "cm",
#     limitsize = FALSE
#   )
#   
#   
#   
#   morris_data <- Morris.indices$IGlobal
#   morris_data <- morris_data %>% arrange(desc(IGlobal_mu.star)) %>% 
#     mutate(n = seq(from=1, to=113, by=1)) %>% 
#     filter(IGlobal_sigma > 0)
#   morris_data <- morris_data %>%  mutate(
#     "ParameterType" = case_when(
#       str_detect(parm_names, regex("MRP2|BCRP")) ~ "Glucuronide excretion",
#       str_detect(parm_names, regex("MRP3")) ~ "Glucuronide absorption",
#       str_detect(parm_names, regex("Hydro")) ~ "Hydrolysis",
#       str_detect(parm_names, regex("L_ZEN_gluc|SI_ZEN_gluc|L_aZEL_gluc|SI_aZEL_gluc")) ~ "Glucuronidation",
#       str_detect(parm_names, regex("L_bZEL|SI_bZEL|L_aZEL|SI_aZEL")) ~ "Reduction",
#       str_detect(parm_names, regex("PSI|PL|PR|PS|PR_gluc|PS_gluc|PF|PK|BPR")) ~ "Partitioning",
#       str_detect(parm_names, regex("ODOSEumol")) ~ "Dose",
#       str_detect(parm_names, regex("Papp")) ~ "Aglycone absorption",
#       str_detect(parm_names, regex("fup|fut")) ~ "Fraction unbound",
#       TRUE ~ "Physiological"))      
# 
# 
# colouring <- c(
#   "Glucuronide excretion" = "#FF8C00",
#   "Glucuronide absorption" = "#FF6D8D",
#   "Hydrolysis" = "#1B60FF",
#   "Glucuronidation" = "#1B85FF",
#   "Reduction" = "#1B35FF",
#   "Aglycone absorption" = "#FF3D8D",
#   "Partitioning" = "#8631A7",
#   "Fraction unbound" = "#4631A7",
#   "Physiological" = "#008475",
#   "Dose" = "black")
# 
# morris_plotMax <- morris_data %>% filter(n <= 20) %>% 
#   ggplot(aes(x = IGlobal_mu.star, y = IGlobal_sigma, 
#              label = parm_names, colour = ParameterType)) +
#   geom_point(alpha = 0.7, size=1) +
#   # geom_jitter(width = .5, alpha = 0.7, size = 1)  +
#   geom_text_repel(aes(label = as.character(parm_names)), 
#                   size=15, 
#                   max.overlaps = Inf,
#                   # min.segment.length = 0.05,
#                   segment.size = 0.1
#   ) +
#   scale_color_manual(name = "Parameter category",
#                      values = colouring, 
#                      breaks = c(
#                        "Dose" = "Dose",
#                        "Physiological" = "Physiological",
#                        "Partitioning" = "Partitioning",
#                        "Fraction unbound" = "Fraction unbound",
#                        "Reduction" = "Reduction",
#                        "Glucuronidation" = "Glucuronidation",
#                        "Hydrolysis" = "Hydrolysis",
#                        "Aglycone absorption" = "Aglycone absorption",
#                        "Glucuronide absorption" = "Glucuronide absorption",
#                        "Glucuronide excretion" = "Glucuronide excretion"
#                      )) +
#   labs(
#     title = "A",
#     x = expression(paste("Mean Effect (μ*)")),
#     y = expression(paste("Standard Deviation (σ)"))
#   ) +
#   CP_theme +
#   theme(legend.position = "none")
# 
# morris_plotQMax <- morris_data %>% filter(n > 20) %>%  filter(n < 40) %>% 
#   ggplot(aes(x = IGlobal_mu.star, y = IGlobal_sigma, 
#              label = parm_names, colour = ParameterType)) +
#   geom_point(alpha = 0.7, size = 1) +
#   # geom_jitter(width = .5, alpha = 0.7, size = 1)  +
#   geom_text_repel(aes(label = as.character(parm_names)), 
#                   size=15, 
#                   max.overlaps = Inf,
#                   # min.segment.length = 0.05,
#                   segment.size = 0.1
#   ) +
#   scale_color_manual(name = "Parameter category",
#                      values = colouring, 
#                      breaks = c(
#                        "Dose" = "Dose",
#                        "Physiological" = "Physiological",
#                        "Partitioning" = "Partitioning",
#                        "Fraction unbound" = "Fraction unbound",
#                        "Reduction" = "Reduction",
#                        "Glucuronidation" = "Glucuronidation",
#                        "Hydrolysis" = "Hydrolysis",
#                        "Aglycone absorption" = "Aglycone absorption",
#                        "Glucuronide absorption" = "Glucuronide absorption",
#                        "Glucuronide excretion" = "Glucuronide excretion"
#                      )) +
#   labs(
#     title = "B",
#     x = expression(paste("Mean Effect (μ*)")),
#     y = expression(paste("Standard Deviation (σ)"))
#   ) +
#   CP_theme +
#   theme(legend.position = "none")
# 
# 
# morris_plotQMin <- morris_data %>% filter(n > 72) %>% filter(n < 92) %>% 
#   ggplot(aes(x = IGlobal_mu.star, y = IGlobal_sigma, 
#              label = parm_names, colour = ParameterType)) +
#   geom_point(alpha = 0.7, size = 1)+
#   # geom_jitter(width = .5, alpha = 0.7, size = 1)  +
#   geom_text_repel(aes(label = as.character(parm_names)), 
#                   size=15, 
#                   max.overlaps = Inf,
#                   # min.segment.length = 0.05,
#                   segment.size = 0.1
#   ) +
#   scale_color_manual(name = "Parameter category",
#                      values = colouring, 
#                      breaks = c(
#                        "Dose" = "Dose",
#                        "Physiological" = "Physiological",
#                        "Partitioning" = "Partitioning",
#                        "Fraction unbound" = "Fraction unbound",
#                        "Reduction" = "Reduction",
#                        "Glucuronidation" = "Glucuronidation",
#                        "Hydrolysis" = "Hydrolysis",
#                        "Aglycone absorption" = "Aglycone absorption",
#                        "Glucuronide absorption" = "Glucuronide absorption",
#                        "Glucuronide excretion" = "Glucuronide excretion"
#                      )) +
#   labs(
#     title = "C",
#     x = expression(paste("Mean Effect (μ*)")),
#     y = expression(paste("Standard Deviation (σ)"))
#   ) +
#   CP_theme +
#   theme(legend.position = "none")
# 
# 
# morris_plotMin <- morris_data %>% filter(n >= 92) %>% 
#   ggplot(aes(x = IGlobal_mu.star, y = IGlobal_sigma, 
#              label = parm_names, colour = ParameterType)) +
#   geom_point(alpha = 0.7, size = 1)+
#   # geom_jitter(width = .5, alpha = 0.7, size = 1)  +
#   geom_text_repel(aes(x = IGlobal_mu.star, y = IGlobal_sigma, 
#                       label = as.character(parm_names)), 
#                   size=15, 
#                   max.overlaps = Inf,
#                   min.segment.length = 0.01,
#                   segment.size = 0.1
#   ) +
#   scale_color_manual(name = "Parameter category",
#                      values = colouring, 
#                      breaks = c(
#                        "Dose" = "Dose",
#                        "Physiological" = "Physiological",
#                        "Partitioning" = "Partitioning",
#                        "Fraction unbound" = "Fraction unbound",
#                        "Reduction" = "Reduction",
#                        "Glucuronidation" = "Glucuronidation",
#                        "Hydrolysis" = "Hydrolysis",
#                        "Aglycone absorption" = "Aglycone absorption",
#                        "Glucuronide absorption" = "Glucuronide absorption",
#                        "Glucuronide excretion" = "Glucuronide excretion"
#                      )) +
#   labs(
#     title = "D",
#     x = expression(paste("Mean Effect (μ*)")),
#     y = expression(paste("Standard Deviation (σ)"))
#   ) +
#   CP_theme +
#   theme(legend.position = "none")
# 
# 
# Morris_Global <- morris_plotMax + 
#   morris_plotQMax + 
#   morris_plotQMin + 
#   morris_plotMin
# 
# ggsave(
#   filename = here(OUTPUT,"Morris_Global.png"),
#   plot = Morris_Global,
#   dpi = 600,  
#   width = 24, 
#   height = 18, 
#   units = "cm",
#   limitsize = FALSE
# )
# 
# 
#   forlabel <- morris_data %>% 
#     filter(n <= 20) %>%
#     ggplot(aes(x = IGlobal_mu.star, y = IGlobal_sigma, 
#                colour = ParameterType)) +
#     # geom_abline(slope = 1, intercept = 0, color = "grey", linewidth = 0.05) +
#     geom_jitter(width = .5, alpha = 0.7) +
#     scale_color_manual(name = "Parameter category",
#                        values = colouring, 
#                        breaks = c(
#                          "Dose" = "Dose",
#                          "Physiological" = "Physiological",
#                          "Partitioning" = "Partitioning",
#                          "Fraction unbound" = "Fraction unbound",
#                          "Reduction" = "Reduction",
#                          "Glucuronidation" = "Glucuronidation",
#                          "Hydrolysis" = "Hydrolysis",
#                          "Aglycone absorption" = "Aglycone absorption",
#                          "Glucuronide absorption" = "Glucuronide absorption",
#                          "Glucuronide excretion" = "Glucuronide excretion"
#                          )) +
#     # geom_text_repel(
#     #   size = 10,
#     #   angle = 0,
#     #   max.overlaps = Inf,
#     #   # segment.color = "#8631A7",
#     #   segment.size = 0.1,
#     #   segment.linetype = "dotted",
#     #   # min.segment.length = 0.01,
#     #   # point.padding = 0.5,
#     #   box.padding = 0.3,
#     #   set.seed(1234)
#     # )   +
#     # geom_label(data = legend_data, aes(log10(x), log10(y), label = label),
#     #            size = 8, label.size = NA, alpha = 0, hjust = 0) +
#     # coord_cartesian(clip = 'off')  +
#     labs(
#       x = expression(paste("Mean Effect (log"[10] ~ "μ*)")),
#       y = expression(paste("Standard Deviation (log"[10] ~ "σ)"))
#     ) +
#     CP_theme +
#     theme(legend.position = "bottom")
#     
#   
#   Plot.morris <- morris_plot #+ ggtitle("Morris Sensitivity Analysis")
#   ggsave(
#     filename = here(OUTPUT,"forlabel.png"),
#     plot = forlabel,
#     dpi = 600,  
#     width = 24, 
#     height = 18, 
#     units = "cm",
#     limitsize = FALSE
#   )