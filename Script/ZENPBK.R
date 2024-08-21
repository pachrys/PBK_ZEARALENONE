# Working Document
rm(list=ls()) # to clear out the global environment


# Libraries-------------------------------------------------------------------#

library(deSolve)
library(tidyverse)
library(ggplot2)
library(here)
library(showtext)
library(FME)
font_add(family = "Garamond", regular = "GARA.TTF")
showtext_auto()

# Working directory-----------------------------------------------------------#
WhatAmITesting <- "BackIntoIt"
TodaysOuput <- here("Output/Data", Sys.Date(), WhatAmITesting)
dir.create(TodaysOuput, WhatAmITesting,recursive = TRUE)
PlotFolderIs <- here("Output/Plots", Sys.Date(), WhatAmITesting)
dir.create(PlotFolderIs, WhatAmITesting,recursive = TRUE)

source(here("Script/Functions/PlottingFunction.R"))                               #call plotting function


# ---------------------------------------------------------------------------- #
# PBK MODEL INITIATION ####
# ---------------------------------------------------------------------------- #

## Extracting input parameters ####
## ------------------------------------------------------ #

# Script for CALCULATED input
INPUT <- read.csv(here("InputData/CALCULATED.INPUT.csv"), header=1) %>%  #The function for calculating input can be called here: source(here("Script/ZENPBK_Input"))                                         #call input script
  pivot_wider(names_from = "Parameter", 
              values_from = "Value")
ZEN.dat <- read.csv(here("InputData/Zen.input.csv"), header=1)     #from {Agahi, 2020 #456} and chemicalize, BP from Shin, Hong et al.2009



## Settings ####
## ------------------------------------------------------ #

# Tout
STARTtime <- 0 # start of simulation (h)
STOPtime <- 24 # end of simulation (h) (for repeated pregnancy in Videmann daily exposure time from D7 to D20 (13days)

Dt <- 1/ 100 # steps of simulation (now it per minute) (could be 1/100, 1/200 or 1 for simulations over years)
Tout <- seq(STARTtime,STOPtime,by=Dt)

# Exposure
PO_CONC = 1 #mg/kg day, Videmann
PO_DOSE = PO_CONC*INPUT$BW / INPUT$MW * 1000 #umol
#PO_REGIMEN = seq(0, 13 * 24, by = 24) # (h) dosing regimen

## Parameters ####
## ------------------------------------------------------ #

PARAMETERS=c(PO_DOSE = PO_DOSE, 
             VInL1 = INPUT$VInL/3,
             VInL2 = INPUT$VInL/3,
             VInL3 = INPUT$VInL/3,
             VBi = INPUT$VBi,
             VIn1 = INPUT$VIn/3,
             VIn2 = INPUT$VIn/3,
             VIn3 = INPUT$VIn/3,
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
             PInb.glu = INPUT$PIn.glu, 
             PLib.glu = INPUT$PLi.glu,
             PKib.glu = INPUT$PKi.glu, 
             PAdb.glu = INPUT$PAd.glu,
             PSlb.glu = INPUT$PSl.glu,
             PRab.glu = INPUT$PRa.glu,
             fub.glu = INPUT$fub.glu,
             QIn = INPUT$QIn,
             QHa = INPUT$QHa,
             QLi = INPUT$QLi,
             QKi = INPUT$QKi,
             QAd = INPUT$QAd, 
             QSl = INPUT$QSl,  
             QRa = INPUT$QRa,
             ka = INPUT$ka.is,
             # ka.is = INPUT$ka.is,
             # ka.sl = INPUT$ka.sl, 
             # ka.fa = INPUT$ka.fa, 
             # ka.fit = INPUT$ka.fit, 
             # kex = INPUT$kex, 
             # kfil = INPUT$kfil, 
             # kbil = INPUT$kbil, 
             # ka.sl.glu = INPUT$ka.sl.glu, 
             # ka.fa.glu = INPUT$ka.fa.glu, 
             # ka.fit.glu = INPUT$ka.fit.glu, 
             # kex.glu = INPUT$kex.glu, 
             # kfil.glu = INPUT$kfil.glu, 
             # kbil.glu = INPUT$kbil.glu, 
             # VmaxAbs = INPUT$VmaxAbs, 
             # KmAbs = INPUT$KmAbs, 
             # VmaxSec = INPUT$VmaxSec, 
             # KmSec = INPUT$KmSec, 
             kTIn = INPUT$kTIn,#total intestinal transit time
             kTIn.3 = INPUT$kTIn/3, #intestine is splited in 3
             # kTsi = 1 / ZEN.dat$kTsi, 
             kTCo = 1 / ZEN.dat$kTco,
             kTge = INPUT$kTge, #gastric emptying rate, "Assumed to mimic bile secretion rate
             CL_Re = INPUT$CL_Re,
             CL_Bi = INPUT$CL_Bi, #biliary clearance
             VmaxMet.In = INPUT$VmaxMet.In, #intestinal glucuronidation
             KmGluc.In = INPUT$KmGluc.In,  #affinity constant zearalenone to UGT
             VmaxMet.Li = INPUT$VmaxMet.Li, #liver glucuronidation
             KmGluc.Li = INPUT$KmGluc.Li,
             # VmaxExc.In.glu = INPUT$VmaxExc.In.glu, #active transPOrt from intestinal tissue to lumen
             # KmGluc.Exc = INPUT$KmGluc.Exc, #affinity constant zearalenone glucuronide to efflux transPOrter
             # VmaxAbs.In.glu = INPUT$VmaxAbs.In.glu, #active transPOrt from lumen to tissue
             # KmGluc.Abs = INPUT$KmGluc.Abs, #affinity constant zearalenone glucuronide to influx transPOrer
             kHydro = 1000 # /hydrolysis constant; assumption as in Veronique
) 


## State ####
## ------------------------------------------------------ #

INITIALstate = c(
  AInL1 = PO_DOSE, 
  AInL2 = 0,
  AInL3 = 0,
  AIn1 = 0, 
  AIn2 = 0,
  AIn3 = 0,
  AFe = 0,
  ALi = 0,
  AKi = 0,
  AUr = 0,
  AAd = 0,
  ASl = 0,
  ARa = 0,
  ABl = 0,
  AIn1.glu = 0,
  AIn2.glu = 0,
  AIn3.glu = 0,
  AInL1.glu = 0,
  AInL2.glu = 0,
  AInL3.glu = 0,
  AFe.glu = 0,
  ALi.glu = 0,
  ABi.glu = 0,
  AKi.glu = 0,
  AUr.glu = 0,
  AAd.glu = 0,
  ASl.glu = 0,
  ARa.glu = 0,
  ABl.glu = 0
)




# ---------------------------------------------------------------------------- #
# PBPK MODEL ####
# ---------------------------------------------------------------------------- #
# source(here("Script/ZENPBK_Function.R"))


## Model Function ####
## ------------------------------------------------------ #

PBK_ZENandGLUC <- function(t,state,parameters){
  with(as.list(c(state,parameters)),{
    
    #----Concentrations-------------------------------------------------------#
    
    # Zearalenone
    
    ## Non perfused organ Concentrations (umol/L)
    CInL1 <- AInL1 / VInL1 #intestinal lumen
    CInL2 <- AInL2 / VInL2 #intestinal lumen
    CInL3 <- AInL3 / VInL3 #intestinal lumen
    # CBi <- ABi / VBi #bile
    
    ## Tissue cellular Concentrations (umol/L)                                    #Abbreviation (C for Concentration and then 2 first letters of tissue name)
    CIn1 <- AIn1 / VIn1 #intestine
    CIn2 <- AIn2 / VIn2 #intestine
    CIn3 <- AIn3 / VIn3 #intestine
    CLi <- ALi / VLi #liver
    CKi <- AKi / VKi #kidney
    CAd <- AAd / VAd #adipose
    CSl <- ASl / VSl #slowly perfused
    CRa <- ARa / VRa #rapidly perfused
    CBl <- ABl / VBl #blood
    
    ## Tissue venous Concentrations (umol/L)                                      #Only for perfused organs Equals tissue Concentration/tissue:blood partition coeefficient
    CvIn1 <- CIn1 / PInb #intestinal tissue
    CvIn2 <- CIn2 / PInb #intestinal tissue
    CvIn3 <- CIn3 / PInb #intestinal tissue
    CvLi <- CLi / PLib #liver
    CvKi <- CKi / PKib #kidney
    CvAd <- CAd / PAdb #adipose
    CvSl <- CSl / PSlb #slowly perfused
    CvRa <- CRa / PRab #rapidly perfused
    
    ## Tissue free Concentrations (umol/L), used for clearance/transport
    CKi_free <- CvKi*fub
    CLi_free <- CvLi*fub
    
    
    # Zearalenone glucuronide
    
    ## Non perfused organ concentrations (umol/L)
    CInL1.glu <- AInL1.glu / VInL1 #intestinal lumen glucuronide
    CInL2.glu <- AInL2.glu / VInL2 #intestinal lumen glucuronide
    CInL3.glu <- AInL3.glu / VInL3 #intestinal lumen glucuronide
    CBi.glu <- ABi.glu / VBi #bile
    
    ## Tissue cellular Concentrations (umol/L)                                    #Abbreviation (C for Concentration and then 2 first letters of tissue name)
    CIn1.glu <- AIn1.glu / VIn1 #intestinal tissue
    CIn2.glu <- AIn2.glu / VIn2 #intestinal tissue
    CIn3.glu <- AIn3.glu / VIn3 #intestinal tissue
    CLi.glu <- ALi.glu / VLi #liver
    CKi.glu <- AKi.glu / VKi #kidney
    CAd.glu <- AAd.glu / VAd #adipose tissue
    CSl.glu <- ASl.glu / VSl #slowly perfused tissue
    CRa.glu <- ARa.glu / VRa #rapidly perfused tissue
    CBl.glu <- ABl.glu / VBl #blood
    
    ## Tissue venous Concentrations (umol/L)                                      #Only for perfused organs Equals tissue Concentration/tissue:blood partition LInefficient
    CvIn1.glu <- CIn1.glu / PInb.glu #Concentration leaving intestinal tissue
    CvIn2.glu <- CIn2.glu / PInb.glu #Concentration leaving intestinal tissue
    CvIn3.glu <- CIn3.glu / PInb.glu
    CvLi.glu <- CLi.glu / PLib.glu #Concentration leaving liver
    CvKi.glu <- CKi.glu / PKib.glu #Concentration leaving kidney
    CvAd.glu <- CAd.glu / PAdb.glu #Concentration leaving adipose
    CvSl.glu <- CSl.glu / PSlb.glu #Concentration leaving slowly perfused tissue
    CvRa.glu <- CRa.glu / PRab.glu #Concentration leaving rapidly perfused tissue
    
    
    #----Differential equations-----------------------------------------------#
    # Rate of change in amounts in each tissue (umol/h)
    
    # Zearalenone-------------------- ##
    # Intestinal Lumen
    dAInL1 <- - ka*AInL1 + kHydro*AInL1.glu - kTIn.3*AInL1 # Intestinal lumen 1
    dAInL2 <- - ka*AInL2 + kHydro*AInL2.glu + kTIn.3*AInL1 - kTIn.3*AInL2 #+ kHydro*AInL1.glu    # Intestinal lumen 2
    dAInL3 <- - ka*AInL3 + kHydro*AInL3.glu + kTIn.3*AInL2 - kTCo*AInL3 #+ kHydro*AInL.glu3  # Intestinal lumen 3
    
    # Feces
    dAFe <- kTCo*AInL3
    
    # Intestinal tissue
    dAIn1 <-  + ka*AInL1 - VmaxMet.In*CIn1/(KmGluc.In+CIn1) + QIn*(CBl - CvIn1) - kTIn.3*AIn1  # Intestine 1
    dAIn2 <-  + ka*AInL2 - VmaxMet.In*CIn2/(KmGluc.In+CIn2) + QIn*(CBl - CvIn2) + kTIn.3*AIn1 - kTIn.3*AIn2     # Intestine 2
    dAIn3 <-  + ka*AInL3 - VmaxMet.In*CIn3/(KmGluc.In+CIn3) + QIn*(CBl - CvIn3) + kTIn.3*AIn2    # Intestine 3    
    
    # Liver
    dALi <- QIn*(CvIn1 + CvIn2 + CvIn3) + QHa*CBl - QLi*CvLi - VmaxMet.Li*CLi/(KmGluc.Li+CLi)
    #       fromIntest fromBlo   toBlood    toglucuronide
    # Kidney
    dAKi <- QKi*(CBl - CvKi) - CL_Re*CKi_free #CL_Re is equal to rat GER
    #       fromtoblood        toUrine
    
    dAUr <- CL_Re*CKi_free
    #       fromKidney
    # Adipose
    dAAd <- QAd*(CBl - CvAd)
    #       fromtoblood
    
    # Slowly Perfused Tissue
    dASl <- QSl*(CBl - CvSl)
    #       fromtoslowly
    
    # Rapidly Perfused Tissue
    dARa <- QRa*(CBl - CvRa)
    #       fromtorapidly
    
    dABl <- QLi*CvLi + QKi*CvKi + QAd*CvAd + QSl*CvSl + QRa*CvRa +
      # from  liver      kidney     adipose    slowly     rapidly
      - ((3*QIn) + QHa + QKi + QAd + QSl + QRa)*CBl
    # to int  hepA  kid   adi   slo   rap
    
    # ## Glucuronide---------------------- ##
    
    # Intestinal tissue
    dAIn1.glu <- VmaxMet.In*CIn1/(KmGluc.In+CIn1) - AIn1.glu + QIn*(CBl.glu - CvIn1.glu) - kTIn.3*AIn1.glu # Intestinal tissue 1
    dAIn2.glu <- VmaxMet.In*CIn2/(KmGluc.In+CIn2) - AIn2.glu + QIn*(CBl.glu - CvIn2.glu) + kTIn.3*AIn1.glu - kTIn.3*AIn2.glu# Intestinal tissue 2
    dAIn3.glu <- VmaxMet.In*CIn3/(KmGluc.In+CIn3) - AIn3.glu + QIn*(CBl.glu - CvIn3.glu) + kTIn.3*AIn2.glu# Intestinal tissue 2
    
    # Intestinal lumen
    dAInL1.glu <- AIn1.glu + kTge*ABi.glu - kTIn.3*AInL1.glu - kHydro*AInL1.glu # Intestinal Lumen 1
    dAInL2.glu <- AIn2.glu + kTIn.3*AInL1.glu - kTIn.3*AInL2.glu - kHydro*AInL2.glu # Intestinal Lumen 2
    dAInL3.glu <- AIn3.glu + kTIn.3*AInL2.glu - kTCo*AInL3.glu - kHydro*AInL3.glu # Intestinal Lumen 2
    
    
    # Feces
    dAFe.glu <- kTCo*AInL3.glu
    
    # Liver
    dALi.glu <- VmaxMet.Li*CLi/(KmGluc.Li+CLi) + QIn*(CvIn1.glu+CvIn2.glu+CvIn3.glu) + QHa*CBl.glu - QLi*CvLi.glu - CL_Bi*CLi.glu
    
    # Bile
    dABi.glu <- CL_Bi*CLi.glu - kTge*ABi.glu
    
    # Kidney
    dAKi.glu <- QKi*(CBl.glu-CvKi.glu) - CL_Re*CKi.glu
    
    # Urine
    dAUr.glu <- CL_Re*CKi.glu
    
    # Adipose
    dAAd.glu <- QAd*(CBl.glu - CvAd.glu)
    
    # Slowly Perfused Tissue
    dASl.glu <- QSl*(CBl.glu - CvSl.glu)
    
    # Rapidly Perfused Tissue
    dARa.glu <- QRa*(CBl.glu - CvRa.glu)
    
    # Blood
    dABl.glu <- QLi*CvLi.glu + QKi*CvKi.glu + QAd*CvAd.glu + QSl*CvSl.glu + QRa*CvRa.glu +
      - (((3*QIn) + QHa + QKi + QAd + QSl + QRa) * CBl.glu)
    
    
    #----Output---------------------------------------------------------------#
    
    return(list(c(dAInL1 = dAInL1, 
                  dAInL2 = dAInL2, 
                  dAInL3 = dAInL3,
                  dAIn1 = dAIn1, 
                  dAIn2 = dAIn2, 
                  dAIn3 = dAIn3, 
                  dAFe = dAFe,
                  dALi = dALi, 
                  dAKi = dAKi, 
                  dAUr = dAUr, 
                  dAAd = dAAd, 
                  dASl = dASl, 
                  dARa = dARa, 
                  dABl = dABl,
                  dAIn1.glu = dAIn1.glu,
                  dAIn2.glu = dAIn2.glu,
                  dAIn3.glu = dAIn3.glu,
                  dAInL1.glu = dAInL1.glu, 
                  dAInL2.glu = dAInL2.glu,
                  dAInL3.glu = dAInL3.glu,
                  dAFe.glu = dAFe.glu,
                  dALi.glu = dALi.glu, 
                  dABi.glu = dABi.glu, 
                  dAKi.glu = dAKi.glu, 
                  dAUr.glu = dAUr.glu,
                  dAAd.glu = dAAd.glu, 
                  dASl.glu = dASl.glu, 
                  dARa.glu = dARa.glu, 
                  dABl.glu = dABl.glu),
                list(c(CInL1 = CInL1, 
                       CInL2 = CInL2, 
                       CInL3 = CInL3,
                       CIn1 = CIn1, 
                       CIn2 = CIn2,
                       CIn3 = CIn3,
                       CLi = CLi ,
                       CKi = CKi, 
                       CAd = CAd, 
                       CSl = CSl, 
                       CRa = CRa, 
                       CBl = CBl, 
                       CIn1.glu = CIn1.glu,
                       CIn2.glu = CIn2.glu,
                       CIn3.glu = CIn3.glu,
                       CInL1.glu = CInL1.glu, 
                       CInL2.glu = CInL2.glu, 
                       CInL3.glu = CInL3.glu,
                       CLi.glu = CLi.glu, 
                       CBi.glu = CBi.glu, 
                       CKi.glu = CKi.glu, 
                       CAd.glu = CAd.glu, 
                       CSl.glu = CSl.glu, 
                       CRa.glu = CRa.glu, 
                       CBl.glu = CBl.glu 
                       ))))
  })
  }


## RUN THE MODEL ####
## ------------------------------------------------------ #

PBK <- ode(y = INITIALstate,
           times = Tout,  #used to be times
           func = PBK_ZENandGLUC,
           parms = PARAMETERS,
           atol = 1e-50,
           #events = list(data = events), #for repeated exposure
           method = "lsoda")
head(PBK)



# ---------------------------------------------------------------------------- #
# RESULTS ####
# ---------------------------------------------------------------------------- #
# plot(PBK)
# ggsave(here(PlotFolderIs, "Overall.png"), dpi = 300)

# Output in amounts (umol)
OUTPUT.A <- as.data.frame(PBK) %>% select("time", "AInL1" : "ABl.glu")
write.csv(OUTPUT.A, file = here(TodaysOuput, "OUTPUT.A.csv"), row.names=FALSE)


## Mass Balance / Error ####
## ------------------------------------------------------ #

# CALCULATED <- c( #!Change accordingly depending on the compartments in the function#
#   "AInL1", "AInL2", "AInL3", "AIn1", "AIn2", "AIn3", "AFe",
#   "ALi", "AKi", "AUr", "AAd", "ASl", "ARa", "ABl",
#   "AIn1.glu", "AIn2.glu", "AIn3.glu", "AInL1.glu", "AInL2.glu", "AInL3.glu","AFe.glu",
#   "ALi.glu", "ABi.glu", "AKi.glu", "AUr.glu",
#   "AAd.glu", "ASl.glu", "ARa.glu", "ABl.glu"
# )

CALCULATED <- names(OUTPUT.A)[-1]
CALCULATED <- rowSums(OUTPUT.A[ ,CALCULATED])
TOTAL <- PO_DOSE

MASSBALANCE <- TOTAL - CALCULATED + 1 #better be 1
ERROR <- (TOTAL - CALCULATED) / (TOTAL) * 100 #better be 0

MASSBALANCE_Error <- cbind("time" = OUTPUT.A[, "time"],
                           MASSBALANCE, ERROR) %>%
  as.data.frame()

# Mass Balance and Error plot
MASSBALANCE_Error %>%
  filter(time >= 0.01) %>% #as PO administration time is 0,01h
  ggplot(aes(time, MASSBALANCE)) +
  geom_line(color = "grey") +
  geom_line(aes(time, ERROR), color = "grey6") +
  labs(x = "Time (h)", y = "Value",
       caption = "Mass balance in grey, Error in black") +
  CP_theme()
ggsave(here(TodaysOuput, "MB_ERROR2=.png"), dpi = 300)


## Plotting Simulated Output ####
## ------------------------------------------------------ #

# Concentrations
OUTPUT.C <- as.data.frame(PBK) %>% select(!c("AInL1" : "ABl.glu")) 

write.csv(OUTPUT.C, file = here(TodaysOuput, "OUTPUT.C.csv"), row.names=FALSE)


# To create df of simulated data for plotting
CONCENTRATIONS <- OUTPUT.C %>% #!change accordingly depending on the compartments in the model
  rename(#!change accordingly depending on the compartments in the model
    "Time" = "time", 
    "Intestine_Lumen1" = "CInL1",
    "Intestine_Lumen2" = "CInL2",
    "Intestine_Lumen3" = "CInL3",
    "Intestine1" = "CIn1",
    "Intestine2" = "CIn2",
    "Intestine3" = "CIn3",
    "Liver" = "CLi",
    "Kidney" = "CKi",
    "Adipose" = "CAd",
    "Slowly_perfused" = "CSl",
    "Rapidly_perfused" = "CRa",
    "Blood" = "CBl", 
    "Intestine1_glucuronide" = "CIn1.glu",
    "Intestine2_glucuronide" = "CIn2.glu",
    "Intestine3_glucuronide" = "CIn3.glu",
    "Intestine_Lumen1_glucuronide" = "CInL1.glu",
    "Intestine_Lumen2_glucuronide" = "CInL2.glu",
    "Intestine_Lumen3_glucuronide" = "CInL3.glu",
    "Liver_glucuronide" = "CLi.glu",
    "Bile_glucuronide" = "CBi.glu",
    "Kidney_glucuronide" = "CKi.glu",
    "Adipose_glucuronide" = "CAd.glu",
    "Slowly_perfused_glucuronide" = "CSl.glu",
    "Rapidly_perfused_glucuronide" = "CRa.glu",
    "Blood_glucuronide" = "CBl.glu" 
  ) %>% pivot_longer(names_to = "Tissue",
               values_to = "Concentration",
               Intestine_Lumen1:Blood_glucuronide) 

# Simulated data for Zearalenone
SIM_CONC.ZEN <- CONCENTRATIONS %>%
  filter(!str_ends(Tissue, "_glucuronide")) %>% 
  mutate(Concentration = Concentration * INPUT$MW) %>%  #converting umol/L to ng/ml
  pivot_wider(names_from = Tissue, values_from = Concentration) %>% 
  mutate(Intestine = Intestine_Lumen1 + Intestine_Lumen2 + Intestine_Lumen3 + Intestine1 + Intestine2 + Intestine3) %>% 
  select(!c("Intestine_Lumen1", "Intestine_Lumen2", "Intestine_Lumen3", "Intestine1", "Intestine2", "Intestine3"))%>%
  pivot_longer(names_to = "Tissue", values_to = "Concentration", Liver:Intestine)

# Simulated data for Zearalenone glucuronide
SIM_CONC.ZENGluc <- CONCENTRATIONS %>%
  filter(str_ends(Tissue, "_glucuronide")) %>%
  mutate(Concentration = Concentration * INPUT$MW.glu) %>% #converting umol/L to ng/ml 
  mutate(Tissue = str_replace(Tissue, "_glucuronide$", "")) %>% 
  pivot_wider(names_from = Tissue, values_from = Concentration) %>% 
  mutate(Intestine = Intestine_Lumen1 + Intestine_Lumen2 + Intestine_Lumen3 + Intestine1 + Intestine2 + Intestine3) %>% 
  select(!c("Intestine_Lumen1", "Intestine_Lumen2", "Intestine_Lumen3", "Intestine1", "Intestine2", "Intestine3")) %>%
  pivot_longer(names_to = "Tissue", values_to = "Concentration", Liver:Intestine)



# To plot Zearalenone blood concentration
SIM_CONC.ZEN %>%
  filter(Tissue == "Blood") %>%
  ggplot(aes(Time, Concentration)) +
  facet_wrap(~Tissue) +
  geom_path(linewidth = 0.5) +
  labs(title = "Zearalenone blood concentration",
       x = "Time (h)",
       y = "Concentration (ng/ml)") +
  CP_theme()
ggsave(here(PlotFolderIs, "ZEN.blood.simulated.png"), dpi = 300)


# To plot Zearalenone tissue concentrations
SIM_CONC.ZEN %>%
 # filter(Tissue != "Intestine") %>% 
  ggplot(aes(Time, Concentration, colour = Tissue)) +
  facet_wrap(~Tissue) +
  geom_path(linewidth = 0.5) +
  labs(title = "Zearalenone tissue concentrations",
       x = "Time (h)",
       y = "Concentration (ng/ml)") +
  CP_theme()
ggsave(here(PlotFolderIs, "ZEN.simulated.png"), dpi = 300)



# To plot Zearalenone and Zearalenone glucuronide blood concentrations
SIM_BLOOD.ZEN <- SIM_CONC.ZEN %>% filter(Tissue == "Blood") %>% mutate(Compound = "Zearalenone") 
SIM_BLOOD.ZENGluc <- SIM_CONC.ZENGluc %>% filter(Tissue == "Blood") %>% mutate(Compound = "Z.Glucuronide")
SIM_BLOOD <- rbind(SIM_BLOOD.ZEN, SIM_BLOOD.ZENGluc)

SIM_BLOOD %>%
  ggplot(aes(Time, Concentration, colour = Compound)) +
  geom_path(linewidth = 0.5) +
  labs(title = "Zearalenone and glucuronide blood concentrations",
       x = "Time (h)",
       y = "Concentration (ng/ml)") +
  CP_theme()
ggsave(here(PlotFolderIs, "ZEN&Gluc.blood.simulated.png"), dpi = 300)


# To plot Zearalenone glucuronide tissue concentrations
SIM_CONC.ZENGluc %>%
  ggplot(aes(Time, Concentration, colour = Tissue)) +
  geom_path(linewidth = 0.5) +
  facet_wrap(~Tissue) +
  labs(title = "Zearalenone glucuronide tissue concentrations",
       x = "Time (h)",
       y = "Concentration (ng/ml)") +
  CP_theme()
ggsave(here(PlotFolderIs, "ZEN_gluc.simulated.png"), dpi = 300)



## Plotting Experimental data ####
## ------------------------------------------------------ #

# Input in vivo data
EXPERIMENTAL <- read.csv(here("InputData/InVivoData.csv"), header=1) %>%
  filter(Author == "Videmann", #I chose Videmann as it had the best recoveries, also Videmann separates Zearalenone from it's metabolites
         Gender != "Pregnant.F", #pregnancy not included yet
         TestCompound == "ZEN", #other options are ZEN_glu and a_ZEL
         AdminType == "PO", #"IV"), #Videmann doesn't have IV, but subcutaneous!!
         Tissue  %in% c("Plasma", "Intestine", "Kidney", "Liver")) %>%  #"Urine", "Colon"
  select(c("Gender", "SamplingTime", "Tissue", "Concentration", "SampleCompound"))%>%
  rename(Time = SamplingTime) %>% #(h)
  mutate(Concentration = Concentration/0.76)  #correcting for the recovery (recovery was 72-80% => averaged that to 76% and then corrected the concentrations for that)


## Plotting Simulated Vs Experimental data for Zearalenone ####
## ------------------------------------------------------ #

# Experimental data for only Zearalenone
EXP.PO_CONC.ZEN <- EXPERIMENTAL %>%
  filter(SampleCompound == "ZEN")

# Blood
EXP_BLOOD.ZEN <- EXP.PO_CONC.ZEN %>% filter(Tissue == "Plasma") %>%
  mutate(Concentration = Concentration*ZEN.dat$BP) %>% #(ng/ml) blood concentration: changing from plasma concentrations to blood concentrations by multiplying with BP
  mutate(Tissue = recode(Tissue, "Plasma" = "Blood")) #renaming plasma to blood


SIM_BLOOD.ZEN %>% ggplot() +
  geom_line(aes(Time, Concentration)) +
  geom_point(data = EXP_BLOOD.ZEN, aes(Time, Concentration, color = Gender)) +
  geom_line(data = EXP_BLOOD.ZEN, aes(Time, Concentration, color = Gender),
            linewidth = 0.1) +
  labs(title = "Zearalenone blood concentrations",
       x = "Time (h)",
       y = "Concentration (ng/ml)",
       color = "Gender") +
  CP_theme()
ggsave(here(PlotFolderIs, "BloodSimulatedInvivoVidemann.png"), dpi = 300)

# Intestine
EXP_INTESTINE <- EXP.PO_CONC.ZEN %>% filter(Tissue =="Intestine")
SIM_INTESTINE <- SIM_CONC.ZEN %>% filter(Tissue == "Intestine") %>% 
  filter(Time>=0.1)

SIM_INTESTINE %>% ggplot() +
  geom_line(aes(Time, Concentration)) +
  geom_point(data = EXP_INTESTINE, aes(Time, Concentration, color = Gender)) +
  geom_line(data = EXP_INTESTINE, aes(Time, Concentration, color = Gender),
            linewidth = 0.1) +
  labs(title = "Zearalenone Intestine concentrations",
       x = "Time (h)",
       y = "Concentration (ng/ml)",
       color = "Gender") +
  CP_theme()
ggsave(here(PlotFolderIs, "IntestineSimulatedInvivo.png"), dpi = 300)

# Liver
EXP_LIVER <- EXP.PO_CONC.ZEN %>% filter(Tissue =="Liver")
SIM_LIVER <- SIM_CONC.ZEN %>% filter(Tissue =="Liver")

SIM_LIVER %>% ggplot() +
  geom_line(aes(Time, Concentration)) +
  geom_point(data = EXP_LIVER, aes(Time, Concentration, color = Gender)) +
  geom_line(data = EXP_LIVER, aes(Time, Concentration, color = Gender),
            linewidth = 0.1) +
  labs(title = "Zearalenone Liver concentrations",
       x = "Time (h)",
       y = "Concentration (ng/ml)",
       color = "Gender") +
  CP_theme()
ggsave(here(PlotFolderIs, "LiverSimulatedInvivo.png"), dpi = 300)

# Kidney
EXP_KIDNEY <- EXP.PO_CONC.ZEN %>% filter(Tissue =="Kidney")
SIM_KIDNEY <- SIM_CONC.ZEN %>% filter(Tissue =="Kidney")

EXP_KIDNEY %>% ggplot() +
  geom_line(aes(Time, Concentration)) +
  geom_point(data = EXP_KIDNEY, aes(Time, Concentration, color = Gender)) +
  geom_line(data = EXP_KIDNEY, aes(Time, Concentration, color = Gender),
            linewidth = 0.1) +
  labs(title = "Zearalenone Kidney concentrations",
       x = "Time (h)",
       y = "Concentration (ng/ml)",
       color = "Gender") +
  CP_theme()
ggsave(here(PlotFolderIs, "KidneySimulatedInvivo.png"), dpi = 300)



## Plotting Simulated Vs Experimental data for Zearalenone Glucuronide ####
## ------------------------------------------------------

# Experimental data for only Zearalenone glucuronide
EXP.PO_CONC.ZENGluc <- EXPERIMENTAL %>%
  filter(SampleCompound == "ZEN_glu")

# Blood
EXP_BLOOD.ZEN <- EXP.PO_CONC.ZENGluc %>% filter(Tissue == "Plasma") %>%
  #mutate(Concentration = Concentration*ZEN.dat$BP) %>% # assuming it's one as I don't have an experimental value
  mutate(Tissue = recode(Tissue, "Plasma" = "Blood")) #renaming plasma to blood


SIM_BLOOD.ZENGluc %>% ggplot() +
  geom_line(aes(Time, Concentration)) +
  geom_point(data = EXP_BLOOD.ZEN, aes(Time, Concentration, color = Gender)) +
  geom_line(data = EXP_BLOOD.ZEN, aes(Time, Concentration, color = Gender),
            linewidth = 0.1) +
  labs(title = "Zearalenone glucuronide blood concentrations",
       x = "Time (h)",
       y = "Concentration (ng/ml)",
       color = "Gender") +
  CP_theme()
ggsave(here(PlotFolderIs, "GlucBloodSimulatedInvivo.png"), dpi = 300)

# Intestine
EXP_INTESTINE <- EXP.PO_CONC.ZENGluc %>% filter(Tissue =="Intestine")
SIM_INTESTINE <- SIM_CONC.ZENGluc %>% filter(Tissue =="Intestine")

SIM_INTESTINE %>% ggplot() +
  geom_line(aes(Time, Concentration)) +
  geom_point(data = EXP_INTESTINE, aes(Time, Concentration, color = Gender)) +
  geom_line(data = EXP_INTESTINE, aes(Time, Concentration, color = Gender),
            linewidth = 0.1) +
  labs(title = "Zearalenone glucuronide Intestine concentrations",
       x = "Time (h)",
       y = "Concentration (ng/ml)",
       color = "Gender") +
  CP_theme()
ggsave(here(PlotFolderIs, "GlucIntestineSimulatedInvivo.png"), dpi = 300)

# Liver
EXP_LIVER <- EXP.PO_CONC.ZENGluc %>% filter(Tissue =="Liver")
SIM_LIVER <- SIM_CONC.ZENGluc %>% filter(Tissue =="Liver")

SIM_LIVER %>% ggplot() +
  geom_line(aes(Time, Concentration)) +
  geom_point(data = EXP_LIVER, aes(Time, Concentration, color = Gender)) +
  geom_line(data = EXP_LIVER, aes(Time, Concentration, color = Gender),
            linewidth = 0.1) +
  labs(title = "Zearalenone glucuronide Liver concentrations",
       x = "Time (h)",
       y = "Concentration (ng/ml)",
       color = "Gender") +
  CP_theme()
ggsave(here(PlotFolderIs, "GlucLiverSimulatedInvivo.png"), dpi = 300)

# Kidney
EXP_KIDNEY <- EXP.PO_CONC.ZENGluc %>% filter(Tissue =="Kidney")
SIM_KIDNEY <- SIM_CONC.ZENGluc %>% filter(Tissue =="Kidney")

SIM_KIDNEY %>% ggplot() +
  geom_line(aes(Time, Concentration)) +
  geom_point(data = EXP_KIDNEY, aes(Time, Concentration, color = Gender)) +
  geom_line(data = EXP_KIDNEY, aes(Time, Concentration, color = Gender),
            linewidth = 0.1) +
  labs(title = "Zearalenone glucuronide Kidney concentrations",
       x = "Time (h)",
       y = "Concentration (ng/ml)",
       color = "Gender") +
  CP_theme()
ggsave(here(PlotFolderIs, "GlucKidneySimulatedInvivo.png"), dpi = 300)




#*****************************************************************************#
#------------------------------Data analysis----------------------------------#
# 
# Rearranging the data--------------------------------------------------------#

# Videmann
VidemannSimple.F.data <- EXP_BLOOD.ZEN %>% #female
  filter(Gender == "Female") %>% 
  select(Time, Concentration) 

VidemannSimple.M.data <- EXP_BLOOD.ZEN %>% #male
  filter(Gender == "Male") %>%
  select(Time, Concentration) 

Videmann.TimeP <- VidemannSimple.F.data$Time

closest_times <- sapply(Videmann.TimeP, function(x) {
  SIM_BLOOD.ZEN$Time[which.min(abs(SIM_BLOOD.ZEN$Time - x))]
})

#%>% filter(Time != 0.705)#time points are the same for both male and female
SimVidemann.data <- SIM_BLOOD.ZEN %>% 
  select(c(Time, Concentration)) %>% 
  filter(Time %in% closest_times) #simulated data is the same for both genders as gender differences are not included

#Calculating AUCs------------------------------------------------------------#

gender <- "female" #choose gender "female" | "male" 

if(gender == "female"){ #selecting gender specific data from the invivo dataframe
  VidemannSimple.data <- VidemannSimple.F.data
} else if(gender == "male"){
  VidemannSimple.data <- VidemannSimple.M.data
}

Videmann.AUC <- 0 #creating empty vector for AUC values from experimental data; unit: ng.h/ml
for (i in 2:nrow(VidemannSimple.data)) {
  Videmann.AUC <- Videmann.AUC +
    ((VidemannSimple.data$Concentration[i-1] + VidemannSimple.data$Concentration[i]) / 2) *
    (VidemannSimple.data$Time[i] - VidemannSimple.data$Time[i-1])
}

SimVidemann.AUC <- 0 #creating empty vector for AUC values from simulated data; unit: ng.h/ml
for (i in 2:nrow(SimVidemann.data)) {
  SimVidemann.AUC <- SimVidemann.AUC +
    ((SimVidemann.data$Concentration[i-1] + SimVidemann.data$Concentration[i]) / 2) *
    (SimVidemann.data$Time[i] - SimVidemann.data$Time[i-1])
} 

if(gender == "female"){ #renaming depending on gender
  Videmann.F.AUC <- Videmann.AUC
} else if(gender == "male"){
  Videmann.M.AUC <- Videmann.AUC
}

Combined.df <- bind_cols( #creating df with combined data
  VidemannSimple.F.data, SimVidemann.data$Concentration) %>%
  rename(Experimental = Concentration, ,
         Simulated = ...3)


# Simulated/Experimental linear regression------------------------------------#
Regression <- lm(Experimental~Simulated, data=Combined.df)
summary(Regression)

Regression %>%
  ggplot(aes(Simulated, Experimental)) +
  geom_point(aes(Simulated, Experimental)) +
  geom_smooth(method='lm', color = "black") +
  CP_theme() +
  theme(legend.position = "bottom")
ggsave(here(PlotFolderIs, "ExpVsSimRegression.png"), dpi = 300)


#---End-----------------------------------------------------------------------#
#*****************************************************************************#
