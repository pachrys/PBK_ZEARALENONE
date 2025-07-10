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
ODOSEmg <- 1 #{mg/kg bw} # oral dose, variable Shin is 8, Videmann = 1
ODOSEumol <- ODOSEmg*1000/MWZEN*Parameters[1, "Value"] # (μmol) unit change to μmol

# IV dose
IVDOSEmg1 <- 0 #{mg/kg bw} # IV dose, variable Shin is 8
IVDOSEumol <-IVDOSEmg1*1000/MWZEN*Parameters[1, "Value"] # unit change to μmol
# time

Dosing <- data.frame(
  Parameter = c("ODOSEumol", "IVDOSEumol"),
  Value = c(ODOSEumol, IVDOSEumol)
)


Parameters <- rbind(Parameters, Dosing)
parm.c <- as.list(setNames(Parameters$Value, Parameters$Parameter))


# Model Code
source(here("Script", "tbdPBK_code.R"))
# source(here("Script", "PBK_code_Mendez.R"))


# Tell Initials 
 
# A_init <- c(
#   AST_ZEN = ODOSEumol,
#  ASIL_ZEN1 = 0,
#  ASIL_ZEN2 = 0,
#  ASIL_ZEN3 = 0,
#  ASIL_ZEN4 = 0,
#  ASIL_ZEN5 = 0,
#  ASIL_ZEN6 = 0,
#  ASIL_ZEN7 = 0,
#  ASI_ZEN = 0,
#   ASI_ZENGLU = 0,
#   ALIL_ZEN = 0,
#   ALIL_aZEL = 0,
#   ALIL_bZEL = 0,
#  AFE_ZEN = 0,
#  AFE_aZEL = 0,
#   AL_ZEN = 0,
#   ALaZEL = 0,
#   ALbZEL = 0,
#   AL_ZENGLU = 0,
#   AB_ZENGLU = 0, #new
#   AUR_ZENGLU = 0, #new
#   AL_aZELGLU = 0,
#   AF = 0,
#   AR = 0,
#   AS = 0,
#   AB = 0,
#   AZur = 0,
#   AFaZEL = 0,
#   ARaZEL = 0,
#   ASaZEL = 0,
#   ABaZEL = 0,
#   AaZur = 0,
#   ATotal = 0
# ) 

A_init <- c(
  AST_ZEN = ODOSEumol,
  ASIL_ZEN = 0, ASIL_ZENGLU = 0,
  ALIL_ZEN = 0, ALIL_aZEL = 0,ALIL_bZEL = 0,
  AREL_ZEN = 0, ARE_ZEN = 0,
  AFE_ZEN = 0, AFE_aZEL = 0, AFE_ZENGLU = 0,
  ASI_ZEN = 0, ASI_ZENGLU = 0, ASI_bZEL = 0, ASI_aZEL = 0, ASI_aZELGLU = 0,
  AL_ZEN = 0, AL_ZENGLU = 0, AL_bZEL = 0, AL_aZEL = 0, AL_aZELGLU = 0,
  AB_ZEN = 0, AB_aZEL = 0, AB_ZENGLU = 0,
  AUR_ZEN = 0, AUR_aZEL = 0, AUR_ZENGLU = 0,
  AF_ZEN = 0, AF_aZEL = 0, AF_ZENGLU = 0,
  AR_ZEN = 0, AR_aZEL = 0, AR_ZENGLU = 0,
  AS_ZEN = 0, AS_aZEL = 0, AS_ZENGLU = 0#,
  # ATotal = 0
  )


PBK.out <- lsoda(y = A_init,
                 times = seq(0, 24, by=1/60),#per minute
                 func =  PBK.model,
                 parms = parm.c)

PBK_OUTPUT <- as.data.frame(PBK.out)

save(PBK_OUTPUT, file = here(OUTPUT, "RawData_PBK.RData"))

source(here("Script", "Results_PBK.R"))