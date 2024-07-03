# Libraries-------------------------------------------------------------------#

library(deSolve)
library(tidyverse)
library(ggplot2)
library(here)
library(showtext)
font_add(family = "Garamond", regular = "GARA.TTF")
showtext_auto()

# Working directory-----------------------------------------------------------#
WhatAmITesting <- "DetailedMetabolism"
TodaysOuput <- here("Output/Data", Sys.Date(), WhatAmITesting)      #folder for saving results
dir.create(TodaysOuput, recursive = TRUE)


#*****************************************************************************#
#---Run script----------------------------------------------------------------#
#*****************************************************************************#

#----Bring PBK function-------------------------------------------------------#
source(here("Script/ZENPBK_Function.R"))


#----Extracting input---------------------------------------------------------#
# Script for calculated input
INPUT <- read.csv(here("InputData/CALCULATED.INPUT.csv"), header=1) %>%  #The function for calculating input can be called here: source(here("Script/ZENPBK_Input"))                                         #call input script
  pivot_wider(names_from = "Parameter", 
              values_from = "Value")

#----Set Parameters for the simulation----------------------------------------#

# Tout
STARTtime <- 0 # start of simulation (h)
STOPtime <- 24 # end of simmulation (h) (can be from 24h to years)
Dt <- 1/60 # steps of simulation (now it per minute) (could be 1/100, 1/200 or 1 for simulations over years)
times <- seq(STARTtime,STOPtime,by=Dt)

# Exposure
PO_CONC = 1 #mg/kg day, Videmann
PO_DOSE = PO_CONC*INPUT$BW/INPUT$MW*1000 #umol

dosing = data.frame(var = c("AInL"),
                   time = c(seq(0, 24, by = 24)), 
                   value = c(PO_DOSE),
                   method = c("add"))


# Model parameters 
PARAMETERS=c(PO_DOSE = PO_DOSE, 
             VInL = INPUT$VInL,
             VBi = INPUT$VBi,
             VIn = INPUT$VIn,
             VLi = INPUT$VLi, 
             VKi = INPUT$VKi, 
             VAd = INPUT$VAd, 
             VSl = INPUT$VSl, 
             VRa = INPUT$VRa, 
             VBl = INPUT$VBl,
             PInb = INPUT$PIn,
             PLib = INPUT$PLi, 
             PKib = INPUT$PKi,
             PAdb = INPUT$PAd,
             PSlb = INPUT$PSl,
             PRab = INPUT$PRa, 
             fub = INPUT$fub,
             PInb.glu = INPUT$PIn, 
             PLib.glu = INPUT$PLi,
             PKib.glu = INPUT$PKi, 
             PAdb.glu = INPUT$PAd,
             PSlb.glu = INPUT$PSl,
             PRab.glu = INPUT$PRa,
             QIn = INPUT$QIn,
             QPo = INPUT$QPo,
             QLi = INPUT$QLi,
             QHa = INPUT$QLi - INPUT$QPo, #hepatic artery (Qliver is sum of both Qhepatic artery and Qportal vein)
             QKi = INPUT$QKi,
             QAd = INPUT$QAd, 
             QSl = INPUT$QSl,  
             QRa = INPUT$QRa,
             ka = INPUT$ka,
             kTIn = INPUT$kTIn,
             kTge = INPUT$kTge, #gastric emptying rate, assumed to mimic bile secretion rate
             CL_Re = INPUT$CL_Re,
             CL_Bi = INPUT$CL_Bi, #biliary clearance
             VmaxMet.In = INPUT$VmaxMet.In, #intestinal glucuronidation
             KmGluc.In = INPUT$KmGluc.In,  #affinity constant zearalenone to UGT
             VmaxMet.Li = INPUT$VmaxMet.Li, #liver glucuronidation
             KmGluc.Li = INPUT$KmGluc.Li,
             VmaxExc.In.glu = INPUT$VmaxExc.In.glu, #active transport from intestinal tissue to lumen
             KmGluc.Exc = INPUT$KmGluc.Exc, #affinity constant zearalenone glucuronide to efflux transporter
             VmaxAbs.In.glu = INPUT$VmaxAbs.In.glu, #active transport from lumen to tissue
             KmGluc.Abs = INPUT$KmGluc.Abs #affinity constant zearalenone glucuronide to influx transporer
             ) 


#Initial state
INITIALstate=c(
  #D = PO_DOSE,
  AInL = 0.0, 
  AIn = 0.0, 
  AFe = 0.0, 
  ALi = 0.0,
  # ABi = 0.0,
  AKi = 0.0, 
  AUr = 0.0, 
  AAd = 0.0, 
  ASl = 0.0, 
  ARa = 0.0, 
  ABl = 0.0,
  AIn.glu = 0.0, 
  AInL.glu = 0.0, 
  AFe.glu = 0.0, 
  ALi.glu = 0.0, 
  ABi.glu = 0.0, 
  AKi.glu = 0.0, 
  AUr.glu = 0.0, 
  AAd.glu = 0.0, 
  ASl.glu = 0.0, 
  ARa.glu = 0.0, 
  ABl.glu = 0.0)


#----Simulation---------------------------------------------------------------#
PBK <- ode(y = INITIALstate, 
           times = times,  
           func = PBK_ZENandGLUC,
           parms = PARAMETERS, 
           atol = 1e-50,
           events = list(data = dosing),
           method = "lsoda")
head(PBK)

#----Results------------------------------------------------------------------#
plot(PBK)

OUTPUT <- as.data.frame(PBK) #(umol)

# Mass Balance
Calculated <- c( #!Change accordingly depending on the compartments in the function#
  "AInL", 
  "AIn", 
  "AFe", 
  "ALi",
  "AKi", 
  "AUr", 
  "AAd", 
  "ASl", 
  "ARa", 
  "ABl",
  "AIn.glu", 
  "AInL.glu", 
  "AFe.glu", 
  "ALi.glu", 
  "ABi.glu", 
  "AKi.glu", 
  "AUr.glu", 
  "AAd.glu", 
  "ASl.glu", 
  "ARa.glu", 
  "ABl.glu")

# Plotting function
source(here("Script/Functions/PlottingFunction.R"))                               #call plotting function


#-----Mass Balance Calculations----------------------------------------------#
Calculated <- rowSums(OUTPUT[ ,Calculated])
Total <- PO_DOSE

Massbalance <- Total - Calculated + 1
ErrorAll <- (Total - Calculated) / (Total + 1*10^(-30)) * 100                     #mass error (% of mass lost)

MassBalance_Error <- cbind("time" = OUTPUT[, "time"],                       #data frame
                           Massbalance, ErrorAll) %>%
  as.data.frame()

# Mass Balance and Error plot
MassBalance_Error %>%
  filter(time >= 0.01) %>% #as PO administration time is 0,01h
  ggplot(aes(time, Massbalance)) +
  geom_line(color = "grey") +
  geom_line(aes(time, ErrorAll), color = "grey6") +
  labs(x = "Time (h)", y = "Value",
       caption = "Massbalance in grey, Error in black") +
  CP_theme()
ggsave(here(TodaysOuput, "MassbalanceError.png"), dpi = 300)             #saving the plot



