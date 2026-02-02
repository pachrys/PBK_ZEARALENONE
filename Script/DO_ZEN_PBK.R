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
  ODOSEmg <- 2.1                                                                # mg/kg bw 
  ODOSEumol <- ((ODOSEmg*1000)/MW_ZEN)                                          # μmol/kg bw
 
  # Load input ----
  
 exp <- ODOSEumol                                                               # umol/kg/day
 Tinput <- 1                                                                    # for repeated exposure or so default = 1
 tinterval <- 1                                                                 # for repeated exposure or so default = 1
 expSTOP <- 7                                                                   # time in days after which the exposure stopped
  
 # Subject-relevant information
 expAGE <- 30                                                                   # years, old age at exposure if not provided then age argument is not used physiology is based on BW
 expBW <- NA                                                                    # kg, if not provided then the BW of the corresponding age and sex is taken; if both BW and Age are not given then a default BW = 70 is taken; if BW is higher than the BW from the lifestage equations then the actual BW overwrites the calculated one
 sex  <-  "F"                                                                   # sex either "F" or "M" if none then default is "M"

 # Simulation relevant information
 Tstart  <-  0                                                                  # days, start of the simulation
 Tstop  <-  2                                                                   # days, stop of the simulation
 Dt  <-  0.01                                                                   # days, iteration steps (decrease/increase depending on run time)

 RawData <- list(
   CALC_Parameters = vector("list", 1), 
   PBK_OUTPUT = vector("list", 1)
   )
    
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
        Tinput = Tinput, 
        tinterval = tinterval, 
        expSTOP = expSTOP, 
        
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
        Tinput = Tinput, 
        tinterval = tinterval, 
        expSTOP = expSTOP, 
        
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

  RawData <- MODEL_OUTPUT
  save(RawData, file = here(OUTPUT, "RawData_ZEN_PBK.RData"))
  
  
  PBK_OUTPUT <- RawData$PBK_OUTPUT
  
  
  MWZEN <- 318.37
  MWZENGLUC <- 494.493
  MWaZEL <- 320.4
  BP <- 0.67 #assuming the same for all zearalenone, a-zearalenol and zearalenone glucuronide
  
  BW <- RawData$CALC_Parameters$BW  ##{Kg} #  body weight rat
  VST <-RawData$CALC_Parameters$FV_ST*BW  # fraction volume of small intestine
  VSIL <- RawData$CALC_Parameters$FV_SIL*BW # fractional volume of small intestinal lumen
  # VCOL <- (RawData$CALC_Parameters$FV_LIL*BW)-VREL
  VSI <- RawData$CALC_Parameters$FV_SI*BW # fractional volume of small intestinal tissue
  VDuoL <- VSIL*RawData$CALC_Parameters$FV_duol 
  VJejL <- VSIL*RawData$CALC_Parameters$FV_jejl 
  VIleL <- VSIL*RawData$CALC_Parameters$FV_ilel 
  VCOL <- RawData$CALC_Parameters$FV_COL*BW 
  
  VDuo <- VSI*RawData$CALC_Parameters$FV_Duo
  VJej <- VSI*RawData$CALC_Parameters$FV_Jej
  VJej_IC <- RawData$CALC_Parameters$FV_I_IC*VJej
  VJej_EC <- RawData$CALC_Parameters$FV_I_EC*VJej
  VIle <- VSI*RawData$CALC_Parameters$FV_Ile
  
  # Create df with organ concentrations
  
  C.Pre <- PBK_OUTPUT %>%
    transmute(Time = time/365,
              DuodenumT_ZEN = CDuo_ZEN, DuodenumL_ZEN = CDuoL_ZEN,
              JejunumT_ZEN = CJej_ZEN, JejunumL_ZEN = CJejL_ZEN,
              IleumT_ZEN = CIle_ZEN, IleumL_ZEN = CIleL_ZEN,
              LIntestineL_ZEN = CCOL_ZEN,
              Liver_ZEN = CL_ZEN, 
              Fat_ZEN = CF_ZEN, Rapidly_ZEN = CR_ZEN, Slowly_ZEN = CS_ZEN,
              Blood_ZEN = CB_ZEN, DuodenumT_aZEL = 0, DuodenumL_aZEL = 0,
              JejunumT_aZEL = CJej_aZEL, JejunumL_aZEL = 0,
              IleumT_aZEL = 0, IleumL_aZEL = 0,
              LIntestineL_aZEL = 0, Liver_aZEL = CL_aZEL,
              Fat_aZEL = CF_aZEL, Rapidly_aZEL = CR_aZEL, Slowly_aZEL = CS_aZEL,
              Blood_aZEL = CB_aZEL,
              DuodenumT_ZENGLU = 0, DuodenumL_ZENGLU = CDuoL_ZENGLU,
              JejunumT_ZENGLU = CJej_ZENGLU, JejunumL_ZENGLU = CJejL_ZENGLU,
              IleumT_ZENGLU = 0, IleumL_ZENGLU = CIleL_ZENGLU, LIntestineL_ZENGLU = CCOL_ZENGLU, JejunumT_ZENGLU = CJej_ZENGLU,
              Liver_ZENGLU = CL_ZENGLU, 
              Fat_ZENGLU = CF_ZENGLU, Rapidly_ZENGLU = CR_ZENGLU, Slowly_ZENGLU = CS_ZENGLU,
              Blood_ZENGLU = CB_ZENGLU)
  
  
  C_ZEN <- last(C.Pre$Blood_ZEN)* 1000 #umol to nmol/L* MWZEN # * 1000 #* 1000
  print(C_ZEN)
  C_aZEL <- last(C.Pre$Blood_aZEL) * 1000 #umol to nmol/L* MWaZEL #* 1000 #* 1000
  print(C_aZEL)
  C_ZENGluc <- last(C.Pre$Blood_ZENGLU) * 1000 #umol to nmol/L* MWZENGLUC
  print(C_ZENGluc)
  # 
  # Cmax_ZEN <- max(C.Pre$Blood_ZEN) * 1000 #umol to nmol/L* MWZEN # * 1000 #* 1000
  # print(Cmax_ZEN)
  # Cmax_aZEL <- max(C.Pre$Blood_aZEL) * 1000 #umol to nmol/L* MWaZEL #* 1000 #* 1000
  # print(Cmax_aZEL)
  # Cmax_ZENGluc <- max(C.Pre$Blood_ZENGLU) * 1000 #umol to nmol/L* MWZENGLUC
  # print(C_ZENGluc)
  