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
  
  
  
  # Set output storage directory
  OUTPUT <- here("Output", format(Sys.Date(), "%Y-%m-%d")) 
  dir.create(OUTPUT, recursive = TRUE)
  
  # Choose if physiology should change with age ("Yes" to include physiological changes due to age and "No" to assume the same physiology over time)
  Lifestage = "No" 
  
  # Load files
  Physio.c <- read_csv(here("Input", "PhysioVariables.csv"))
  Tissue.c <- read_csv(here("Input", "TissueComposition.csv"))
  PhysChem <- read_delim("Input/PhysChem.csv", 
                         delim = ";", escape_double = FALSE, trim_ws = TRUE)
  source(here("Script", "RUN_and_OUTPUT.R"))
  
  MW_ZEN <- pull(PhysChem[PhysChem$Substance == "Zearalenone", 
                               "MolecularWeight"])
  
  # Oral dose
  ODOSEmg.BW <- 45e-3 #0.2/1000  #100/70                                        # mg/kg bw 
  ODOSEumol.BW <- ((ODOSEmg.BW*1000)/MW_ZEN)                                    # μmol/kg bw
 
  # Load input ----

  Tstart  <-  0                                                                 # days, start of the simulation
  Tstop  <- 1                                                                   # days, stop of the simulation
  Dt  <- 1/(24*60)# round(1/10,3) #(24*60)                                      # days, iteration steps (decrease/increase depending on run time)
  
  
  exp <- ODOSEumol.BW                                                           # umol/kg/day
  Tinput <- 1                                                                   # for repeated exposure with in a day (let's say at every meal or so), default = 1
  tinterval <- 1                                                                # for repeated exposure or so default = 1
  expSTART <- 0                                                                 # time in days of the begining of exposure
  expSTOP <- 1#13*Dt                                                            # time in days after which the exposure stopped
  dosingT <-  c(0:Dt)#c(expSTOP,1+expSTOP, 2+expSTOP)
  # exp1STOP <- 1+(13*Dt) 
  # exp2STOP <- 2+(13*Dt) 
  # exp3STOP <- 3+(13*Dt) 
 
  # Subject-relevant information
  expAGE <- 30                                                                  # years, old age at exposure if not provided then age argument is not used physiology is based on BW
  expBW <- NA    #60 #                                                          # kg, if not provided then the BW of the corresponding age and sex is taken; if both BW and Age are not given then a default BW = 70 is taken; if BW is higher than the BW from the lifestage equations then the actual BW overwrites the calculated one
  sex  <-  "F"                                                                  # sex either "F" or "M" if none then default is "M"


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
 save(MODEL_OUTPUT, file = here(OUTPUT, "EC50.ZEN.0.06.RawData_ZEN_PBK.RData"))
  
 print(MODEL_OUTPUT$Cmaxs)
  
 # Out <- MODEL_OUTPUT$PBK_OUTPUT
 # Urine <- Out %>%  select(c(time, starts_with("AUR_")))
 # write.csv(Urine, here(OUTPUT, "UrineMirocha.csv"), row.names = FALSE)
 # Dosing <- MODEL_OUTPUT$Dosing