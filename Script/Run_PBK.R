# --------------------------------------------------------------------------- #
# Script for the PBK Model of DIANA et al. 
# By: Chrysanthi Pachoulide
# Date: 23-06-2025
# --------------------------------------------------------------------------- #

rm(list=ls()) 

# Packages
library(here)
library(tidyverse)
library(deSolve)
library(PKNCA)
library(pracma)
library(showtext)
library(plotly)
font_add(family = "Garamond", regular = "GARA.TTF")
showtext_auto()
 
CP_theme <- theme_minimal() +
  theme(
    axis.title = element_text(size = 44, face = "bold"),
    axis.text = element_text(size = 40),
    plot.title = element_text(hjust = 0, size = 48, face = "bold"),
    plot.margin = margin(0.3, 0.3, 0.3, 0.3, "cm"), 
    axis.line = element_line(size = 0.1),
    panel.grid = element_line(size = 0.3),
    plot.caption = element_text(hjust = 0, size = 40, face= "italic"),
    plot.title.position = "plot", 
    plot.caption.position = "plot",
    legend.position = "none"
  ) 


# Set output storage directory
OUTPUT <- here("Output", format(Sys.Date(), "%Y-%m-%d"), 
               format(Sys.time(), "%H-%M-%S"))
dir.create(OUTPUT, recursive = TRUE)

# Species
Select_Species = "Human" #Human
Select_Author = "Videmann" #Shin, Q.Lu, Videmann
Select_AdminType = "PO" #IV, PO

# Parameters
Parameters <- read.csv(here("Input", "Parms.csv"))%>% 
  pivot_longer(Rat:Human, names_to = "Species", values_to = "Value")

Parameters <- if(Select_Species == "Rat"){
    
    Parameters <- Parameters %>% filter(Species == "Rat")
    
    }else{
      
      Parameters <- Parameters %>% filter(Species == "Human")
      
} %>% select(!Species)

# Dosing
MWZEN <- 318.37 # molecular weight ZEN

# oral dose
ODOSEmg <-  0.25/1000 #29/1000000# {mg/kg bw} # oral dose, variable Shin is 8, Videmann = 1, Q.Lu = 6.65
ODOSEumol <- ((ODOSEmg*1000)/MWZEN)*Parameters[1, "Value"][[1]] # (μmol) unit change to μmol

# IV dose
IVDOSEmg1 <- 0 #{mg/kg bw} # IV dose, variable Shin is 8, Q.Lu = 6.65
IVDOSEumol <-((IVDOSEmg1*1000)/MWZEN)*Parameters[1, "Value"][[1]] # unit change to μmol
# time

Dosing <- data.frame(Parameter = c("ODOSEumol", "IVDOSEumol"),
                     Value = c(ODOSEumol, IVDOSEumol))

Species <- data.frame(Parameter = "Species", 
                      Value = Select_Species)

Parameters <- rbind(Parameters, Dosing, Species)
parm.c <- as.list(setNames(Parameters$Value, Parameters$Parameter)) #%>% unlist()
parm.c[-length(parm.c)] <- lapply(parm.c[-length(parm.c)], as.numeric)

# Model Code
source(here("Script", "PBKi_code.R"))
# source(here("Script", "PBK_code_Mendez.R"))

# Initials
if(Select_AdminType == "PO"){
  
  A_init <- c(
    
    AST_ZEN = ODOSEumol, 
    
    ADuoL_ZEN = 0, AJejL_ZEN = 0, AIleL_ZEN = 0, 
    ADuoL_ZENGLU = 0, AJejL_ZENGLU = 0, AIleL_ZENGLU = 0, 
    
    ACOL_ZEN = 0, ACOL_ZENGLU = 0, ACOL_aZEL = 0, ACOL_bZEL = 0, 
    AFE_ZEN = 0, AFE_aZEL = 0, AFE_ZENGLU = 0, 
    
    ADuo_ZEN = 0, AJej_ZEN = 0, AIle_ZEN = 0, 
    
    AJej_IC_ZENGLU = 0, AJej_EC_ZENGLU = 0, #AIle_IC_ZENGLU = 0, AIle_EC_ZENGLU = 0,
    AJej_aZEL = 0, AJej_bZEL = 0, AJej_aZELGLU = 0, #AJej_ZENGLU = 0, 
    
    AL_ZEN = 0, AL_IC_ZENGLU = 0, AL_EC_ZENGLU = 0, AL_bZEL = 0, AL_aZEL = 0, AL_aZELGLU = 0, #AL_ZENGLU = 0
    ABile_ZENGLU = 0,
    
    AB_ZEN = 0, AB_aZEL = 0, AB_ZENGLU = 0, 
    
    AUR_ZEN = 0, AUR_aZEL = 0, AUR_ZENGLU = 0, 
    AF_ZEN = 0, AF_aZEL = 0, AF_ZENGLU = 0, 
    AR_ZEN = 0, AR_aZEL = 0, AR_ZENGLU = 0, 
    AS_ZEN = 0, AS_aZEL = 0, AS_ZENGLU = 0
    
  )
  }else{
    
    A_init <- c(
      AST_ZEN = 0, 
      
      ADuoL_ZEN = 0, AJejL_ZEN = 0, AIleL_ZEN = 0, 
      ADuoL_ZENGLU = 0, AJejL_ZENGLU = 0, AIleL_ZENGLU = 0, 
      
      ACOL_ZEN = 0, ACOL_ZENGLU = 0, ACOL_aZEL = 0, ACOL_bZEL = 0, 
      AFE_ZEN = 0, AFE_aZEL = 0, AFE_ZENGLU = 0, 
      
      ADuo_ZEN = 0, AJej_ZEN = 0, AIle_ZEN = 0, 
      
      AJej_IC_ZENGLU = 0, AJej_EC_ZENGLU = 0, #AIle_IC_ZENGLU = 0, AIle_EC_ZENGLU = 0,
      AJej_aZEL = 0, AJej_bZEL = 0, AJej_aZELGLU = 0, #AJej_ZENGLU = 0, 
      
      AL_ZEN = 0, AL_IC_ZENGLU = 0, AL_EC_ZENGLU = 0, AL_bZEL = 0, AL_aZEL = 0, AL_aZELGLU = 0, #AL_ZENGLU = 0
      ABile_ZENGLU = 0,
      
      AB_ZEN = IVDOSEumol, AB_aZEL = 0, AB_ZENGLU = 0, 
      
      AUR_ZEN = 0, AUR_aZEL = 0, AUR_ZENGLU = 0, 
      AF_ZEN = 0, AF_aZEL = 0, AF_ZENGLU = 0, 
      AR_ZEN = 0, AR_aZEL = 0, AR_ZENGLU = 0, 
      AS_ZEN = 0, AS_aZEL = 0, AS_ZENGLU = 0
    ) 
    }



PBK.out <- lsoda(y = A_init,
                 times = seq(0, 24, by=1/60),#per minute
                 func =  H.PBK.model,
                 parms = parm.c
)

  

# PBK.out <-  if(Species == "Rat"){
#   
#   lsoda(y = A_init,
#         times = seq(0, 24, by=1/60),#per minute
#         func =  R.PBK.model,
#         parms = parm.c
#         )
#   
#   }else{
#     
#     lsoda(y = A_init,
#         times = seq(0, 24, by=1/60),#per minute
#         func =  H.PBK.model,
#         parms = parm.c
#         )
#     
#     }

PBK_OUTPUT <- as.data.frame(PBK.out)

save(PBK_OUTPUT, file = here(OUTPUT, "RawData_PBK.RData"))

source(here("Script", "Results_PBK.R"))