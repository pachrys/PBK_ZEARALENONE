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
  source(here("Script", "PBK_model_4GSA_Rat.R"))
  
  
  
  # SENSITIVITY FUNCTION ####
  # ---------------------------------------------------------------------------- #

  ## eFAST ####
  
  GSA_model <- function(parm.GSA, parm.c) {

# Final parameters per test
names(parm.GSA) <- GSA.parms
parm.final <- parm.c
parm.final[names(parm.GSA)] <- parm.GSA


A_init <- c(
  AOD_ZEN = 0.82,
  
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
 func =  R.PBK.model,
 parms = parm.final)

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
  
  
  # SET eFAST PARAMETERS ####
  # ---------------------------------------------------------------------------- #
  
  # Parameter upper and lower bounds

  # Parameters
  Parameters <- read.csv(here("Input", "Parameters_Rat.csv"))
  eFASTParameters <- read.csv(here("Input", "Parameters_toeFAST_Rat.csv"))$parm_names
 
  P <- Parameters %>% filter(Parameter %in% eFASTParameters)
  P$Value <- as.numeric(P$Value)
  P <- P %>% mutate(Binf = Value*0.9, 
   Bsup = Value*1.1, 
   Distribution = "qunif")
  

  parm.c <- as.list(setNames(Parameters$Value, Parameters$Parameter))

  
  # SET eFAST EXPERIMENT DESIGN ####
  # ---------------------------------------------------------------------------- #
  
  GSA.parms <- P$Parameter

  eFAST.factors <- nrow(P)

  # Define q: q needs to be a list of character strings, giving the names of the quantile functions
  q <- ifelse(P$Distribution == "LogNormal", "qlnorm", "qunif") 
  
  # Define q.arg: q.arg needs to be a list of lists
  q.arg <- lapply(1:nrow(P), function(i) { 
    if (P$Distribution[i] == "qunif") { 
      list(min = P$Binf[i], max = P$Bsup[i]) 
      } else {
        list(meanlog = log(P$Mean[i]), sdlog = P$sdlog[i]) 
        }
  })
 
  names(q.arg) <- GSA.parms
  
  # eFAST test function
  set.seed(1234)

  length(q) == length(q.arg)

  Experience <- fast99(
    model = NULL, 
    factors = GSA.parms, # These are the parameters that will be varying in the model
    n = 1000,# Integer giving the sample size, i.e. the length of the discretization of the s-space
    q = q,
    q.arg = q.arg
    )

  dim(Experience$X) # dataframe of: n factors(number of parameters) * n n(1000) number of observations, of n factors(number of parameters) number of variables
  design <- Experience$X
  write.csv(design, file = here(OUTPUT, "eFAST.ExpDesign.csv"), row.names = FALSE)

  pdf("eFAST.Distr_ParFull.pdf")
  par(mfrow=c(4,3))
  for( i in 1:ncol(design)){ hist(design[,i], 
  breaks=100, 
  col="purple", 
  main=colnames(design)[i] ) }
  dev.off()
  save(Experience, file = here(OUTPUT,"eFAST.ExperienceFull.RData"))
  
  # RUN eFAST DESIGN ####
  # ---------------------------------------------------------------------------- #

  # Run model row-wise from design (i.e each row is an eFAST test). parm.GSA = row
  results <- apply(design, 1, function(row) {
    GSA_model(row, parm.c = parm.c)  
  })
  
  y <- results
  y <- as.data.frame(y) 
  save(y, file = here(OUTPUT, "eFAST.y.RData"))
  
  # CALCULATE eFAST INDICES ####
  # ---------------------------------------------------------------------------- #
  
  sim.results.eFAST <- as.matrix(y)

  for (i in (1:length(sim.results.eFAST[,1])))
    {
    tell(Experience, sim.results.eFAST[i,])
    }

  plot(Experience)
  print(Experience)

  plot(Experience)  

  Variance <- Experience$V # Total variance
  names(Variance) <- GSA.parms
  
  Done <- Experience$D1
  names(Done) <- GSA.parms
  
  Dt <- Experience$Dt
  names(Dt) <- GSA.parms
  
  first_order <- Done / Variance  # D1 is the estimated Variance of the Conditional Expectation (VCE) with respect to each factor, normalized by total variance
  names(first_order) <- GSA.parms
  
  total_order <- 1 - (Dt / Variance)  # Dt is the estimated VCE with respect to each factor complementary set of factors ("all but Xi"), normalized by total variance
  names(total_order) <- GSA.parms
  
  lowry_data <- data.frame(
    Parameter = GSA.parms,
    Main.Effect = first_order,
    Interaction = total_order)
  
  write.csv(lowry_data, file = here(OUTPUT, "eFASTresults.csv"), row.names = FALSE)
  
  lowry_data <- read.csv(here(OUTPUT, "eFASTresults.csv"))
  
  ordered_data <- lowry_data %>%
    arrange(desc(Main.Effect)) %>% 
    mutate(
      Parameter = factor(Parameter, levels = Parameter),
      Total.effect = Main.Effect + Interaction,
      Cumulative.main = cumsum(Main.Effect),
      Cumulative.total = cumsum(Total.effect),
      ymin = Cumulative.main, 
      ymax = pmin(Cumulative.total, 1))
  
  long_data <- ordered_data %>% 
    pivot_longer(c(Main.Effect, Interaction))
  long_data$name <- factor(long_data$name, levels = c("Main.Effect", "Interaction"))
  

  lowry_plot <- 
    ggplot(long_data) +
    geom_col(aes(x = Parameter, y = value, fill = name),
             position = position_stack(reverse = TRUE),
             width = 0.8,
             alpha = 0.8) +
    geom_ribbon(
      data = ordered_data,
      aes(x = as.numeric(Parameter),
          ymin = ymin, 
          ymax = ymax),
      fill = "grey50", alpha = 0.3, color = "grey50", linewidth = 0.2) +
    scale_fill_manual(
      values = c("Main.Effect" = "grey5", "Interaction" = "grey50"),
      labels = c("Main Effect", "Interaction")) +
    scale_y_continuous(
      limits = c(0, 1),
      expand = expansion(mult = c(0, 0.05)),
      labels = percent_format(),
      name = "Sensitivity Index") +
    labs(
      x = "Parameter",
      title = "eFAST Sensitivity Analysis") +
    CP_theme +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      legend.spacing.y = unit(0, "cm"),
      legend.margin = margin(0, 0, 0, 0),
      axis.text.x = element_text(size = 30, angle = 45, hjust = 1,),
      axis.line = element_line(size = 0.1)
      )
  
  lowry_plot
  ggsave(
    filename = here(OUTPUT, "lowry_plot.png"),
    plot = lowry_plot,
    dpi = 600, 
    width = 16, 
    height = 9, 
    units = "cm")
  
  
  
  CB_ZEN.sim.results.eFAST <- sim.results.eFAST[str_detect(rownames(sim.results.eFAST), "ZEN") &
                                                    !str_detect(rownames(sim.results.eFAST), "ZEN_gluc"), ]
  

  for (i in (1:length(CB_ZEN.sim.results.eFAST[,1])))
  {
    tell(Experience, CB_ZEN.sim.results.eFAST[i,])
  }

  plot(Experience)
  print(Experience)

  Variance <- Experience$V # Total variance
  names(Variance) <- GSA.parms
  
  Done <- Experience$D1
  names(Done) <- GSA.parms
  
  Dt <- Experience$Dt
  names(Dt) <- GSA.parms
  
  first_order <- Done / Variance  # D1 is the estimated Variance of the Conditional Expectation (VCE) with respect to each factor, normalized by total variance
  names(first_order) <- GSA.parms
  
  total_order <- 1 - (Dt / Variance)  # Dt is the estimated VCE with respect to each factor complementary set of factors ("all but Xi"), normalized by total variance
  names(total_order) <- GSA.parms
  
  CB_ZEN.lowry_data <- data.frame(
    Parameter = GSA.parms,
    Main.Effect = first_order,
    Interaction = total_order)
  
  write.csv(CB_ZEN.lowry_data, file = here(OUTPUT,"CB_ZEN.eFASTresults.csv"), row.names = FALSE)
  
  CB_ZEN.lowry_data <- read.csv(here(OUTPUT,"CB_ZEN.eFASTresults.csv"))
  
  CB_ZEN.ordered_data <- CB_ZEN.lowry_data %>%
    arrange(desc(Main.Effect)) %>% 
    mutate(
      Parameter = factor(Parameter, levels = Parameter),
      Total.effect = Main.Effect + Interaction,
      Cumulative.main = cumsum(Main.Effect),
      Cumulative.total = cumsum(Total.effect),
      ymin = Cumulative.main, 
      ymax = pmin(Cumulative.total, 1))
  
  CB_ZEN.long_data <- CB_ZEN.ordered_data %>% 
    pivot_longer(c(Main.Effect, Interaction))
  CB_ZEN.long_data$name <- factor(CB_ZEN.long_data$name, levels = c("Main.Effect", "Interaction"))
  
  
  CB_ZEN.lowry_plot <- 
    ggplot(CB_ZEN.long_data) +
    geom_col(aes(x = Parameter, y = value, fill = name),
             position = position_stack(reverse = TRUE),
             width = 0.8,
             alpha = 0.8) +
    geom_ribbon(
      data = CB_ZEN.ordered_data,
      aes(x = as.numeric(Parameter),
          ymin = ymin, 
          ymax = ymax),
      fill = "grey50", alpha = 0.3, color = "grey50", linewidth = 0.2) +
    scale_fill_manual(
      values = c("Main.Effect" = "grey5", "Interaction" = "grey50"),
      labels = c("Main Effect", "Interaction")) +
    scale_y_continuous(
      limits = c(0, 1),
      expand = expansion(mult = c(0, 0.05)),
      labels = percent_format(),
      name = "Sensitivity Index") +
    labs(
      x = "Parameter",
      title = "eFAST Sensitivity Analysis") +
    CP_theme +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      legend.spacing.y = unit(0, "cm"),
      legend.margin = margin(0, 0, 0, 0),
      axis.text.x = element_text(size = 30, angle = 45, hjust = 1,),
      axis.line = element_line(size = 0.1)
    )
  
  CB_ZEN.lowry_plot
  ggsave(
    filename = here(OUTPUT,"CB_ZEN.lowry_plot.png"),
    plot = CB_ZEN.lowry_plot,
    dpi = 600, 
    width = 16, 
    height = 9, 
    units = "cm")
  
  
  
  
  
  CB_aZEL.sim.results.eFAST <- sim.results.eFAST[str_detect(rownames(sim.results.eFAST), "aZEL") &
                                                     str_detect(rownames(sim.results.eFAST), "aZEL_gluc"), ]

  
  for (i in (1:length(CB_aZEL.sim.results.eFAST[,1])))
  {
    tell(Experience, CB_aZEL.sim.results.eFAST[i,])
  }
  
  plot(Experience)
  print(Experience)
  
  Variance <- Experience$V # Total variance
  names(Variance) <- GSA.parms
  
  Done <- Experience$D1
  names(Done) <- GSA.parms
  
  Dt <- Experience$Dt
  names(Dt) <- GSA.parms
  
  first_order <- Done / Variance  # D1 is the estimated Variance of the Conditional Expectation (VCE) with respect to each factor, normalized by total variance
  names(first_order) <- GSA.parms
  
  total_order <- 1 - (Dt / Variance)  # Dt is the estimated VCE with respect to each factor complementary set of factors ("all but Xi"), normalized by total variance
  names(total_order) <- GSA.parms
  
  CB_aZEL.lowry_data <- data.frame(
    Parameter = GSA.parms,
    Main.Effect = first_order,
    Interaction = total_order)
  
  write.csv(CB_aZEL.lowry_data, file = here(OUTPUT,"CB_aZEL.eFASTresults.csv"), row.names = FALSE)
  
  CB_aZEL.lowry_data <- read.csv(here(OUTPUT,"CB_aZEL.eFASTresults.csv"))
  
  CB_aZEL.ordered_data <- CB_aZEL.lowry_data %>%
    arrange(desc(Main.Effect)) %>% 
    mutate(
      Parameter = factor(Parameter, levels = Parameter),
      Total.effect = Main.Effect + Interaction,
      Cumulative.main = cumsum(Main.Effect),
      Cumulative.total = cumsum(Total.effect),
      ymin = Cumulative.main, 
      ymax = pmin(Cumulative.total, 1))
  
  CB_aZEL.long_data <- CB_aZEL.ordered_data %>% 
    pivot_longer(c(Main.Effect, Interaction))
  CB_aZEL.long_data$name <- factor(CB_aZEL.long_data$name, levels = c("Main.Effect", "Interaction"))
  
  
  CB_aZEL.lowry_plot <- 
    ggplot(CB_aZEL.long_data) +
    geom_col(aes(x = Parameter, y = value, fill = name),
             position = position_stack(reverse = TRUE),
             width = 0.8,
             alpha = 0.8) +
    geom_ribbon(
      data = CB_aZEL.ordered_data,
      aes(x = as.numeric(Parameter),
          ymin = ymin, 
          ymax = ymax),
      fill = "grey50", alpha = 0.3, color = "grey50", linewidth = 0.2) +
    scale_fill_manual(
      values = c("Main.Effect" = "grey5", "Interaction" = "grey50"),
      labels = c("Main Effect", "Interaction")) +
    scale_y_continuous(
      limits = c(0, 1),
      expand = expansion(mult = c(0, 0.05)),
      labels = percent_format(),
      name = "Sensitivity Index") +
    labs(
      x = "Parameter",
      title = "eFAST Sensitivity Analysis") +
    CP_theme +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      legend.spacing.y = unit(0, "cm"),
      legend.margin = margin(0, 0, 0, 0),
      axis.text.x = element_text(size = 30, angle = 45, hjust = 1,),
      axis.line = element_line(size = 0.1)
    )
  
  CB_aZEL.lowry_plot
  ggsave(
    filename = here(OUTPUT,"CB_aZEL.lowry_plot.png"),
    plot = CB_aZEL.lowry_plot,
    dpi = 600, 
    width = 16, 
    height = 9, 
    units = "cm")
  
  
  
  
  CB_ZEN_gluc.sim.results.eFAST  <- sim.results.eFAST [str_detect(rownames(sim.results.eFAST), "ZEN_gluc"), ]
  
  
  
  for (i in (1:length(CB_ZEN_gluc.sim.results.eFAST[,1])))
  {
    tell(Experience, CB_ZEN_gluc.sim.results.eFAST[i,])
  }
  
  plot(Experience)
  print(Experience)
  
  Variance <- Experience$V # Total variance
  names(Variance) <- GSA.parms
  
  Done <- Experience$D1
  names(Done) <- GSA.parms
  
  Dt <- Experience$Dt
  names(Dt) <- GSA.parms
  
  first_order <- Done / Variance  # D1 is the estimated Variance of the Conditional Expectation (VCE) with respect to each factor, normalized by total variance
  names(first_order) <- GSA.parms
  
  total_order <- 1 - (Dt / Variance)  # Dt is the estimated VCE with respect to each factor complementary set of factors ("all but Xi"), normalized by total variance
  names(total_order) <- GSA.parms
  
  CB_ZEN_gluc.lowry_data <- data.frame(
    Parameter = GSA.parms,
    Main.Effect = first_order,
    Interaction = total_order)
  
  write.csv(CB_ZEN_gluc.lowry_data, file = here(OUTPUT,"CB_ZEN_gluc.eFASTresults.csv"), row.names = FALSE)
  
  CB_ZEN_gluc.lowry_data <- read.csv(here(OUTPUT,"CB_ZEN_gluc.eFASTresults.csv"))
  
  CB_ZEN_gluc.ordered_data <- CB_ZEN_gluc.lowry_data %>%
    arrange(desc(Main.Effect)) %>% 
    mutate(
      Parameter = factor(Parameter, levels = Parameter),
      Total.effect = Main.Effect + Interaction,
      Cumulative.main = cumsum(Main.Effect),
      Cumulative.total = cumsum(Total.effect),
      ymin = Cumulative.main, 
      ymax = pmin(Cumulative.total, 1))
  
  CB_ZEN_gluc.long_data <- CB_ZEN_gluc.ordered_data %>% 
    pivot_longer(c(Main.Effect, Interaction))
  CB_ZEN_gluc.long_data$name <- factor(CB_ZEN_gluc.long_data$name, levels = c("Main.Effect", "Interaction"))
  
  
  CB_ZEN_gluc.lowry_plot <- 
    ggplot(CB_ZEN_gluc.long_data) +
    geom_col(aes(x = Parameter, y = value, fill = name),
             position = position_stack(reverse = TRUE),
             width = 0.8,
             alpha = 0.8) +
    geom_ribbon(
      data = CB_ZEN_gluc.ordered_data,
      aes(x = as.numeric(Parameter),
          ymin = ymin, 
          ymax = ymax),
      fill = "grey50", alpha = 0.3, color = "grey50", linewidth = 0.2) +
    scale_fill_manual(
      values = c("Main.Effect" = "grey5", "Interaction" = "grey50"),
      labels = c("Main Effect", "Interaction")) +
    scale_y_continuous(
      limits = c(0, 1),
      expand = expansion(mult = c(0, 0.05)),
      labels = percent_format(),
      name = "Sensitivity Index") +
    labs(
      x = "Parameter",
      title = "eFAST Sensitivity Analysis") +
    CP_theme +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      legend.spacing.y = unit(0, "cm"),
      legend.margin = margin(0, 0, 0, 0),
      axis.text.x = element_text(size = 30, angle = 45, hjust = 1,),
      axis.line = element_line(size = 0.1)
    )
  
  CB_ZEN_gluc.lowry_plot
  ggsave(
    filename = here(OUTPUT,"CB_ZEN_gluc.lowry_plot.png"),
    plot = CB_ZEN_gluc.lowry_plot,
    dpi = 600, 
    width = 16, 
    height = 9, 
    units = "cm")
  
  CB_aZEL_gluc.sim.results.eFAST  <- sim.results.eFAST [str_detect(rownames(sim.results.eFAST), "aZEL_gluc"), ]

  for (i in (1:length(CB_aZEL_gluc.sim.results.eFAST[,1])))
  {
    tell(Experience, CB_aZEL_gluc.sim.results.eFAST[i,])
  }
  
  plot(Experience)
  print(Experience)
  
  Variance <- Experience$V # Total variance
  names(Variance) <- GSA.parms
  
  Done <- Experience$D1
  names(Done) <- GSA.parms
  
  Dt <- Experience$Dt
  names(Dt) <- GSA.parms
  
  first_order <- Done / Variance  # D1 is the estimated Variance of the Conditional Expectation (VCE) with respect to each factor, normalized by total variance
  names(first_order) <- GSA.parms
  
  total_order <- 1 - (Dt / Variance)  # Dt is the estimated VCE with respect to each factor complementary set of factors ("all but Xi"), normalized by total variance
  names(total_order) <- GSA.parms
  
  CB_aZEL_gluc.lowry_data <- data.frame(
    Parameter = GSA.parms,
    Main.Effect = first_order,
    Interaction = total_order)
  
  write.csv(CB_aZEL_gluc.lowry_data, file = here(OUTPUT,"CB_aZEL_gluc.eFASTresults.csv"), row.names = FALSE)
  
  CB_aZEL_gluc.lowry_data <- read.csv(here(OUTPUT,"CB_aZEL_gluc.eFASTresults.csv"))
  
  CB_aZEL_gluc.ordered_data <- CB_aZEL_gluc.lowry_data %>%
    arrange(desc(Main.Effect)) %>% 
    mutate(
      Parameter = factor(Parameter, levels = Parameter),
      Total.effect = Main.Effect + Interaction,
      Cumulative.main = cumsum(Main.Effect),
      Cumulative.total = cumsum(Total.effect),
      ymin = Cumulative.main, 
      ymax = pmin(Cumulative.total, 1))
  
  CB_aZEL_gluc.long_data <- CB_aZEL_gluc.ordered_data %>% 
    pivot_longer(c(Main.Effect, Interaction))
  CB_aZEL_gluc.long_data$name <- factor(CB_aZEL_gluc.long_data$name, levels = c("Main.Effect", "Interaction"))
  
  
  CB_aZEL_gluc.lowry_plot <- 
    ggplot(CB_aZEL_gluc.long_data) +
    geom_col(aes(x = Parameter, y = value, fill = name),
             position = position_stack(reverse = TRUE),
             width = 0.8,
             alpha = 0.8) +
    geom_ribbon(
      data = CB_aZEL_gluc.ordered_data,
      aes(x = as.numeric(Parameter),
          ymin = ymin, 
          ymax = ymax),
      fill = "grey50", alpha = 0.3, color = "grey50", linewidth = 0.2) +
    scale_fill_manual(
      values = c("Main.Effect" = "grey5", "Interaction" = "grey50"),
      labels = c("Main Effect", "Interaction")) +
    scale_y_continuous(
      limits = c(0, 1),
      expand = expansion(mult = c(0, 0.05)),
      labels = percent_format(),
      name = "Sensitivity Index") +
    labs(
      x = "Parameter",
      title = "eFAST Sensitivity Analysis") +
    CP_theme +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      legend.spacing.y = unit(0, "cm"),
      legend.margin = margin(0, 0, 0, 0),
      axis.text.x = element_text(size = 30, angle = 45, hjust = 1,),
      axis.line = element_line(size = 0.1)
    )
  
  CB_aZEL_gluc.lowry_plot
  ggsave(
    filename = here(OUTPUT,"CB_aZEL_gluc.lowry_plot.png"),
    plot = CB_aZEL_gluc.lowry_plot,
    dpi = 600, 
    width = 16, 
    height = 9, 
    units = "cm")
  
  