  # --------------------------------------------------------------------------- #
  # Script for running the eFAST Test
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
  library(readxl)
  library(purrr)

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
  
  # Model Code
  source(here("Script", "PBK_code_Mendez.R"))
  
  
  # SENSITIVITY FUNCTION ####
  # ---------------------------------------------------------------------------- #

  ## eFAST ####
  
  GSA_model <- function(parm.GSA, parm.c) {
    
    # Final parameters per test
    names(parm.GSA) <- GSA.parms
    parm.final <- parm.c
    parm.final[names(parm.GSA)] <- parm.GSA
    
    # (Optional) Print parameters for debugging
    print(parm.final[names(parm.GSA)])
    
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
                 times = seq(0, 10, 1),
                 func = PBK.model,
                 parms = parm.final)
    
    # Return selected outputs
    CB <- PBK.out[,"CB"]
    AUC.CB <- sum(diff(PBK.out[,"time"]) * (PBK.out[-nrow(PBK.out),"CB"] + PBK.out[-1,"CB"]))/2
    
    CLZEN <- PBK.out[,"CLZEN"]
    AUC.CLZEN <- sum(diff(PBK.out[,"time"]) * (PBK.out[-nrow(PBK.out),"CLZEN"] + PBK.out[-1,"CLZEN"]))/2
    
    CLZENGLU <- PBK.out[,"CLZENGLU"]
    AUC.CLZEN <- sum(diff(PBK.out[,"time"]) * (PBK.out[-nrow(PBK.out),"CLZENGLU"] + PBK.out[-1,"CLZENGLU"]))/2
    
    CBaZEL <- PBK.out[,"CBaZEL"]
    AUC.CBaZEL <- sum(diff(PBK.out[,"time"]) * (PBK.out[-nrow(PBK.out),"CBaZEL"] + PBK.out[-1,"CBaZEL"]))/2
    
    
    return(c(CB = CB, AUC.CB = AUC.CB,
             CLZEN = CLZEN, AUC.CLZEN = AUC.CLZEN, 
             CLZENGLU = CLZENGLU, AUC.CLZEN = AUC.CLZEN,
             CBaZEL = CBaZEL, AUC.CBaZEL = AUC.CBaZEL)) 
  }
  
  
  # SET eFAST PARAMETERS ####
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
  
  Parameters <- rbind(Parameters, Dosing)
  
  P <- Parameters %>% 
    filter(Parameter %in% c(
      "BW", 
      "VSI", "VLc", "VFc",
      "QSIc", "QLc", "QFc",
      "PLZEN", "PFZEN", 
      "S9SI", 
      "VmaxSIZENGLUc", "KmSIZENGLU", 
      "VLS9", 
      "VmaxLaZELc", "KmLaZEL", 
      "VmaxLZENGLUc", "KmLZENGLU", 
      "VmazLaZELGLUc", "KmLaZELGLU",
      "ODOSEumol", "IVDOSEumol"
    )) %>% 
    mutate(Binf = Value*0.9, 
           Bsup = Value*1.1, 
           Distribution = "qunif")
  
  # P <- Parameters %>% select(Abbreviation, Distribution, Initial_value, Binf, Bsup, Mean, sdlog, zscore) %>%
  #   filter(Abbreviation %in% c("VILc", "fup", "COral", "R_PTL", "QT", "SF_OAT", "Km_OAT4c", "BW", "GFRc", "QC")) # "SA_SI", "QKc", "SF_OATP1B1", "Vmax_OAT4c", "Km_OATP1B1c", "tco", "Papp_SI", "VPc", "VKc", "VPTc", "Hct", "PAc", "Vmax_OATP1B1c", "Km_OATP1B3c"))

  parm.c <- as.list(setNames(Parameters$Value, Parameters$Parameter))

  
  # SET eFAST EXPERIMENT DESIGN ####
  # ---------------------------------------------------------------------------- #
  
  GSA.parms <- P$Parameter

  eFAST.factors <- nrow(P)

  # Define q: q needs to be a list of character strings, giving the names of the quantile functions
  q <- ifelse(P$Distribution == "LogNormal", "qlnorm", "qunif") 
  
  # Define q.arg: q.arg needs to be a list of lists
  q.arg <- lapply(1:nrow(P), function(i) { # apply the function to eachone of the rows of P
    if (P$Distribution[i] == "qunif") { #if the distribution in uniform
      list(min = P$Binf[i], max = P$Bsup[i]) #create a list containing min (Bing of row i) and max (Bsup of row i) per row
    } else {
      list(meanlog = log(P$Mean[i]), sdlog = P$sdlog[i]) #if not then create a list containing mean and sdlog
    }
  })
 
  names(q.arg) <- GSA.parms
  
  # eFAST test function
  set.seed(1234)

  length(q) == length(q.arg)

  Experience <- fast99(
    model = NULL, 
    factors = GSA.parms, # These are the parameters that will be varying in the model
    n = 1000,            # Integer giving the sample size, i.e. the length of the discretization of the s-space
    q = q,
    q.arg = q.arg
  )

  dim(Experience$X) # dataframe of: n factors(number of parameters) * n n(1000) number of observations, of n factors(number of parameters) number of variables
  design <- Experience$X
  write.csv(design, "eFAST.ExpDesign.csv", row.names = FALSE)

  pdf("eFAST.Distr_ParFull.pdf")
  par(mfrow=c(4,3))
  for( i in 1:ncol(design)){ hist(design[,i], 
                                  breaks=100, 
                                  col="purple", 
                                  main=colnames(design)[i] ) }
  dev.off()
  save(Experience, file = "eFAST.ExperienceFull.RData")
  
  # for (i in 1:nrow(P)) {
  #   if (P$Distribution[i] != "Uniform") {
  #     par_name <- GSA.parms[i]
  #     min_val <- exp(log(P$Mean[i]) - P$zscore[i] * P$sdlog[i])
  #     max_val <- exp(log(P$Mean[i]) + P$zscore[i] * P$sdlog[i])
  #
  #     Experience$X[, par_name] <- pmax(Experience$X[, par_name], min_val)
  #     Experience$X[, par_name] <- pmin(Experience$X[, par_name], max_val)
  #   }
  # }
  # 
  # design <- Experience$X
  # write.csv(design, "eFAST.ExpDesign.csv", row.names = FALSE)
  #
  # pdf("Distr_ParFull.pdf")
  # par(mfrow=c(4,3))
  # for( i in 1:ncol(design)){ hist(design[,i], breaks=100, col="purple", main=colnames(design)[i] ) }
  # dev.off()
  # save(Experience, file = "eFAST.ExperienceFull.RData")

  
  # RUN eFAST DESIGN ####
  # ---------------------------------------------------------------------------- #

  # Run model row-wise from design (i.e each row is an eFAST test). parm.GSA = row
  results <- apply(design, 1, function(row) {
    GSA_model(row, parm.c = parm.c)  
  })
  
  y <- results
  y <- as.data.frame(y) 
  save(y, file = "eFAST.y.RData")
  
  # CALCULATE eFAST INDICES ####
  # ---------------------------------------------------------------------------- #
  
  sim.results.eFAST <- as.matrix(y)

  for (i in (1:length(sim.results.eFAST[,1])))
  {
    tell(Experience, sim.results.eFAST[i,])
  }
  # AUC <- results["AUC",]
  # CP <- results[1:11,]
  # Experience$y <- CP
  # tell(Experience)
  # plot(Experience)
  # 
  plot(Experience)

  plot(Experience)  

  Variance <- Experience$V # Total variance
  names(Variance) <- GSA.parms
  
  Done <- Experience$D1
  names(Done) <- GSA.parms
  
  Dt <- Experience$Dt
  names(Dt) <- GSA.parms
  
  first_order <- Done / Variance  # D1 is the estimated Variance of the Conditional Expectation (VCE) with respect to each factor, normalized by total variance
  names(first_order) <- GSA.parms
  
  total_order <- Dt / Experience$V  # Dt is the estimated VCE with respect to each factor complementary set of factors ("all but Xi"), normalized by total variance
  names(total_order) <- GSA.parms
  
  lowry_data <- data.frame(
    Parameter = GSA.parms,
    Main.Effect = first_order,
    Interaction = total_order) %>% 
    mutate(
      Main.Effect = Main.Effect/sum(Main.Effect), # to normalise
      Interaction = Interaction/sum(Interaction)  # to normalise
    ) 
  
  write.csv(lowry_data, file = here("eFASTresults.csv"), row.names = FALSE)
  
  lowry_data <- read.csv(here("eFASTresults.csv"))
  
  ordered_data <- lowry_data %>%
    arrange(desc(Main.Effect)) %>%
    mutate(
      Parameter = factor(Parameter, levels = Parameter), 
      Total.effect = Main.Effect + Interaction,
      Cumulative.main = cumsum(Main.Effect) 
    )
  long_data <- ordered_data %>% pivot_longer(c(Main.Effect, Interaction))
  long_data$name <- factor(long_data$name, levels = c("Main.Effect", "Interaction"))
  
  
  lowry_plot <- 
    ggplot(long_data) +
    geom_col(aes(x = Parameter, y = value, fill = name),
             position = position_stack(reverse = FALSE),
             width = 0.8,
             alpha = 0.8) +
    geom_ribbon(
      data = ordered_data,
      aes(x = as.numeric(Parameter),
          ymin = lag(Cumulative.main, default = 0),
          ymax = Cumulative.main),
      fill = "grey50", alpha = 0.3, color = "grey50", linewidth = 0.2
    ) +
    scale_fill_manual(
      values = c("Main.Effect" = "darkblue", "Interaction" = "blueviolet"),
      labels = c("Main Effect", "Interaction")
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      expand = expansion(mult = c(0, 0.05)),
      labels = percent_format(),
      name = "Sensitivity Index"
    ) +
    labs(
      x = "Parameter",
      title = "eFAST Sensitivity Analysis"
    ) +
    CP_theme +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = 42),
      legend.key.size = unit(0.5, "lines"),
      legend.spacing.y = unit(0, "cm"),
      legend.margin = margin(0, 0, 0, 0),
      axis.text.x = element_text(angle = 45, hjust = 1,),
      plot.margin = margin(0, 0, 0, 0.15, "cm")
    )
  
  lowry_plot
  ggsave(
    filename = here("lowry_plot.png"),
    plot = lowry_plot,
    dpi = 1000, 
    width = 12, height = 9, units = "cm"
  )