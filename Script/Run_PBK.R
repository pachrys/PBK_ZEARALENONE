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
  read.csv(here("Input", "Parm.rat.csv"))
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
source(here("Script", "PBK_code.R"))
# source(here("Script", "PBK_code_Mendez.R"))


# Tell Initials 
A_init <- c(
  Ast = ODOSEumol,
  Ain = 0,
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
  ABZENGLU = 0,
  AZENGLUur = 0,
  ALaZELGLU = 0,
  AF = 0,
  AR = 0,
  AS = 0,
  AB = IVDOSEumol,
  AZur = 0,
  AFaZEL = 0,
  ARaZEL = 0,
  ASaZEL = 0,
  ABaZEL = 0,
  AaZur = 0,
  ATotal = 0
  )


PBK.out <- lsoda(y = A_init,
                 times = seq(0, 24, by=1/60),#per minute
                 func =  PBK.model,
                 parms = parm.c)

PBK_OUTPUT <- as.data.frame(PBK.out)

save(PBK_OUTPUT, file = here(OUTPUT, "RawData_PBK.RData"))

source(here("Script", "Results_PBK.R"))