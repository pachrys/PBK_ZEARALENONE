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
    axis.text = element_text(size = 56),
    axis.title = element_text(size = 62),
    plot.title = element_text(size = 65, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 56, hjust = 0.5),
    legend.position = "none",
    plot.margin = margin(0.3, 0.3, 0.3, 0.3, "cm") 
  )

# Set output storage directory
OUTPUT <- here("Output", format(Sys.Date(), "%Y-%m-%d"), format(Sys.time(), "%H-%M-%S"))
dir.create(OUTPUT, recursive = TRUE)

# Species
Species = "Rat"
Select_Author = "Shin" #Shin, Q.Lu
Select_AdminType = "IV" #IV, PO

# Parameters
Parameters <- if(Species == "Rat"){
  read.csv(here("Input", "tbdParm.rat.csv"))
}else{
  read.csv(here("Input", "Parm.human.csv"))
}

# Dosing
# molecular weight
MWZEN <- 318.37 # molecular weight ZEN

# oral dose
ODOSEmg <- 1 #{mg/kg bw} # oral dose, variable Shin is 8, Videmann, Q.Lu = 1
ODOSEumol <- ODOSEmg*1000/MWZEN*Parameters[1, "Value"] # (μmol) unit change to μmol

# IV dose
IVDOSEmg1 <- 8 #{mg/kg bw} # IV dose, variable Shin is 8, Q.Lu = 6.63
IVDOSEumol <-IVDOSEmg1*1000/MWZEN*Parameters[1, "Value"] # unit change to μmol
# time

Dosing <- data.frame(
  Parameter = c("ODOSEumol", "IVDOSEumol"),
  Value = c(ODOSEumol, IVDOSEumol)
)


Parameters <- rbind(Parameters, Dosing)
parm.c <- as.list(setNames(Parameters$Value, Parameters$Parameter))


# Model Code
source(here("Script", "PBKi_code.R"))
# source(here("Script", "PBK_code_Mendez.R"))

# Initials
if(Select_AdminType == "PO"){
  A_init <- c(
    AST_ZEN = ODOSEumol, 
    ADuoL_ZEN = 0, AJejL_ZEN = 0, AIleL_ZEN = 0, 
    ADuoL_ZENGLU = 0, AJejL_ZENGLU = 0, AIleL_ZENGLU = 0, 
    ALIL_ZEN = 0, ALIL_ZENGLU = 0, ALIL_aZEL = 0, ALIL_bZEL = 0, 
    AREL_ZEN = 0, 
    AFE_ZEN = 0, AFE_aZEL = 0, AFE_ZENGLU = 0, 
    ADuo_ZEN = 0, AJej_ZEN = 0, AIle_ZEN = 0, 
    AJej_ZENGLU = 0, AJej_aZEL = 0, AJej_bZEL = 0, AJej_aZELGLU = 0,
    AL_ZEN = 0, AL_ZENGLU = 0, AL_bZEL = 0, AL_aZEL = 0, AL_aZELGLU = 0,
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
    ALIL_ZEN = 0, ALIL_ZENGLU = 0, ALIL_aZEL = 0, ALIL_bZEL = 0, 
    AREL_ZEN = 0, 
    AFE_ZEN = 0, AFE_aZEL = 0, AFE_ZENGLU = 0, 
    ADuo_ZEN = 0, AJej_ZEN = 0, AIle_ZEN = 0, 
    AJej_ZENGLU = 0, AJej_aZEL = 0, AJej_bZEL = 0, AJej_aZELGLU = 0,
    AL_ZEN = 0, AL_ZENGLU = 0, AL_bZEL = 0, AL_aZEL = 0, AL_aZELGLU = 0,
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
                 func =  PBK.model,
                 parms = parm.c)

PBK_OUTPUT <- as.data.frame(PBK.out)

save(PBK_OUTPUT, file = here(OUTPUT, "RawData_PBK.RData"))

source(here("Script", "Results_PBK.R"))