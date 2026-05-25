  # --------------------------------------------------------------------------- #
  # Script for running the PBK model and obtaining the results
  # By: Chrysanthi Pachoulide
  # Date: 30-01-2026
  # --------------------------------------------------------------------------- #
  
  rm(list=ls()) 
  
  
  # Packages
  library(here)
  library(tidyverse)
  library(deSolve)
  library(showtext)
  library(plotly)
  library(readr)
  font_add(family = "Garamond", regular = "GARA.TTF")
  showtext_auto()
  
  # For plotting
  
  CP_theme <- theme_minimal() +
    theme(
      plot.background = element_rect(fill = 'white'),
      plot.margin = margin(0.3, 0.3, 0.3, 0.3, "cm"), 
      plot.title = element_text(hjust = 0, size = 48, face = "bold"),
      plot.title.position = "plot", 
      panel.background = element_blank(), 
      panel.grid = element_line(linewidth = 0.1),  
      axis.title = element_text(size = 44, face = "bold"),
      axis.text = element_text(size = 44),
      legend.text = element_text(size = 44),
      legend.title = element_text(size = 44), 
      legend.key.spacing = unit(0.1, "cm"),
      legend.key.size = unit(0.5, "lines"),
      legend.spacing.y = unit(0.2, "cm"), 
      legend.spacing.x = unit(0.2, "cm"),
      legend.margin = margin(0,0,0,0),
      legend.position = "right"
    )
  
  test = "Ado_M"
  expAGE <- 15 #30 #30                                                          # years, old age at exposure if not provided then age argument is not used physiology is based on BW
  expBW <- NA #                                                                 # kg, if not provided then the BW of the corresponding age and sex is taken; if both BW and Age are not given then a default BW = 70 is taken; if BW is higher than the BW from the lifestage equations then the actual BW overwrites the calculated one
  sex  <-  "F" #"F"                                                             # sex either "F" or "M" if none then default is "M"
  
  # Set output storage directory
  OUTPUT <- here("Output", format(Sys.Date(), "%Y-%m-%d")) 
  dir.create(OUTPUT, recursive = TRUE)
  
  # Choose if physiology should change with age ("Yes" to include physiological changes due to age and "No" to assume the same physiology over time)
  Lifestage = "Yes" 
  
  # Load files
  Physio.c <- read_csv(here("Input", "PhysioVariables.csv"))
  Tissue.c <- read_csv(here("Input", "TissueComposition.csv"))
  PhysChem <- read_delim("Input/PhysChem.csv", 
                         delim = ";", escape_double = FALSE, trim_ws = TRUE)
  source(here("Script", "RUN_and_OUTPUT.R"))
  
  MW_ZEN <- pull(PhysChem[PhysChem$Substance == "Zearalenone", 
                               "MolecularWeight"])
  
  # Oral dose
  ODOSEmg.BW <- 0.25/1000
    # TDI 0.25/1000, Hokkaido  0.016/1000, Tokai 0.020/1000, McKeon 0.0331/1000, Mirocha 100/70, # Warth 0.25/1000                                         # mg/kg bw 
  ODOSEumol.BW <- ((ODOSEmg.BW*1000)/MW_ZEN)                                    # μmol/kg bw
 
  # Load input ----

  Tstart  <-  0                                                                 # days, start of the simulation
  Tstop  <- 10
  # McKeon, Mirocha 1, Warth 4                                                          # days, stop of the simulation
  Dt  <- 1/(24)
  # McKeon, Mirocha 1/(24*60), Warth 1/24
  
  exp <- ODOSEumol.BW                                                           # umol/kg/day
  Tinput <- 1                                                                   # for repeated exposure with in a day (let's say at every meal or so), default = 1
  tinterval <- 1                                                                # for repeated exposure or so default = 1
  expSTART <- 0                                                                 # time in days of the begining of exposure
  expSTOP <- 1 
  # McKeon,Mirocha 1, Warth 4                                                          # time in days after which the exposure stopped
  dosingT <- c(0:1) #, 2, 4, 5) 
  # McKeon,Mirocha c(0:Dt), Warth c(0.4)
 
  # Subject-relevant information
  # expAGE <- 1 #30 #30                                                              # years, old age at exposure if not provided then age argument is not used physiology is based on BW
  expBW <- NA #                                                                 # kg, if not provided then the BW of the corresponding age and sex is taken; if both BW and Age are not given then a default BW = 70 is taken; if BW is higher than the BW from the lifestage equations then the actual BW overwrites the calculated one
  # sex  <-  "F" #"F"                                                             # sex either "F" or "M" if none then default is "M"


  if(Lifestage == "Yes"){
    expBW <- NA                                                                 # forcing BW to be NA as BW is taken based on age in the lifestage model
    
    if(is.na(expAGE)){
        stop("expAGE required to run the lifestage option")
      }
    }

 
 # RUN Model ----
  MODEL_OUTPUT <- if(Lifestage == "Yes"){ 
      
      RUNandOUT_lifestage(
        
        # Exposure-relevant information
        exp = exp, 
        expSTOP = expSTOP, 
        dosingT = dosingT,

        # Subject-relevant information
        expAGE = expAGE, 
        expBW = expBW, 
        sex = sex, 
        
        # Simulation relevant information
        Tstart = Tstart, 
        Tstop = Tstop, 
        Dt = Dt 
        
        )
      
    } else { # Individual, no lifestage
      
      RUNandOUT(
        
        # Exposure-relevant information
        exp = exp, 
        expSTOP = expSTOP, 
        dosingT = dosingT,
        
        # Subject-relevant information
        expAGE = expAGE, 
        expBW = expBW, 
        sex = sex, 
        
        # Simulation relevant information
        Tstart = Tstart, 
        Tstop = Tstop, 
        Dt = Dt 
        
      )
      
    }
  
  
 # Save Results ----

 # RawData <- MODEL_OUTPUT
 save(MODEL_OUTPUT, file = here(OUTPUT, paste0(test, ".RData")))
  
 # print(MODEL_OUTPUT$Cmaxs)
 # 
 # Out <- MODEL_OUTPUT$PBK_OUTPUT
 # Urine <- Out %>%  select(c(time, starts_with("AUR_")))
 # write.csv(Urine, here(OUTPUT, "Tokai.csv"), row.names = FALSE)
 # Dosing <- MODEL_OUTPUT$Dosing