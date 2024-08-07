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
WhatAmITesting <- "Dividingto3Intestines"
TodaysOuput <- here("Output/Data", Sys.Date(), WhatAmITesting)
dir.create(TodaysOuput, WhatAmITesting,recursive = TRUE)
PlotFolderIs <- here("Output/Plots", Sys.Date(), WhatAmITesting)
dir.create(PlotFolderIs, WhatAmITesting,recursive = TRUE)

source(here("Script/Functions/PlottingFunction.R"))                               #call plotting function


# -----------------------------------------------------------------------------
# PBK MODEL INITIATION ####
# -----------------------------------------------------------------------------

## Extracting input parameters ####
## ------------------------------------------------------

# Script for CALCULATED input
INPUT <- read.csv(here("InputData/CALCULATED.INPUT.csv"), header=1) %>%  #The function for calculating input can be called here: source(here("Script/ZENPBK_Input"))                                         #call input script
  pivot_wider(names_from = "Parameter", 
              values_from = "Value")
ZEN.dat <- read.csv(here("InputData/Zen.input.csv"), header=1)     #from {Agahi, 2020 #456} and chemicalize, BP from Shin, Hong et al.2009



## Settings ####
## ------------------------------------------------------

# Tout
STARTtime <- 0 # start of simulation (h)
STOPtime <- 24 #end of simulation (h) (for repeated pregnancy in Videmann daily exposure time from D7 to D20 (13days)

Dt <- 1/ 100 # steps of simulation (now it per minute) (could be 1/100, 1/200 or 1 for simulations over years)
Tout <- seq(STARTtime,STOPtime,by=Dt)

# Exposure
PO_CONC = 1 #mg/kg day, Videmann
PO_DOSE = PO_CONC*INPUT$BW / INPUT$MW * 1000 #umol
#PO_REGIMEN = seq(0, 13 * 24, by = 24) # (h) dosing regimen

## Parameters ####
## ------------------------------------------------------

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
             kTIn.3 = INPUT$kTIn/3, #total intestinal transit time
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
## ------------------------------------------------------

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
  AIn.glu1 = 0,
  AIn.glu2 = 0,
  AIn.glu3 = 0,
  AInL.glu1 = 0,
  AInL.glu2 = 0,
  AInL.glu3 = 0,
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




# -----------------------------------------------------------------------------
# PBPK MODEL ####
# -----------------------------------------------------------------------------
# source(here("Script/ZENPBK_Function.R"))

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
    
    
    # Zearalenone glucuronide
    
    ## Non perfused organ concentrations (umol/L)
    CInL.glu1 <- AInL.glu1 / VInL1 #intestinal lumen glucuronide
    CInL.glu2 <- AInL.glu2 / VInL2 #intestinal lumen glucuronide
    CInL.glu3 <- AInL.glu3 / VInL3 #intestinal lumen glucuronide
    
    CBi.glu <- ABi.glu / VBi #bile
    
    ## Tissue cellular Concentrations (umol/L)                                    #Abbreviation (C for Concentration and then 2 first letters of tissue name)
    CIn.glu1 <- AIn.glu1 / VIn1 #intestinal tissue
    CIn.glu2 <- AIn.glu2 / VIn2 #intestinal tissue
    CIn.glu3 <- AIn.glu3 / VIn3 #intestinal tissue
    
    CLi.glu <- ALi.glu / VLi #liver
    CKi.glu <- AKi.glu / VKi #kidney
    CAd.glu <- AAd.glu / VAd #adipose tissue
    CSl.glu <- ASl.glu / VSl #slowly perfused tissue
    CRa.glu <- ARa.glu / VRa #rapidly perfused tissue
    CBl.glu <- ABl.glu / VBl #blood
    
    ## Tissue venous Concentrations (umol/L)                                      #Only for perfused organs Equals tissue Concentration/tissue:blood partition LInefficient
    CvIn.glu1 <- CIn.glu1 / PInb.glu #Concentration leaving intestinal tissue
    CvIn.glu2 <- CIn.glu2 / PInb.glu #Concentration leaving intestinal tissue
    CvIn.glu3 <- CIn.glu3 / PInb.glu #Concentration leaving intestinal tissue
    
    CvLi.glu <- CLi.glu / PLib.glu #Concentration leaving liver
    CvKi.glu <- CKi.glu / PKib.glu #Concentration leaving kidney
    CvAd.glu <- CAd.glu / PAdb.glu #Concentration leaving adipose
    CvSl.glu <- CSl.glu / PSlb.glu #Concentration leaving slowly perfused tissue
    CvRa.glu <- CRa.glu / PRab.glu #Concentration leaving rapidly perfused tissue
    
    
    #----Differential equations-----------------------------------------------#
    # Rate of change in amounts in each tissue (umol/h)
    
    # Zearalenone-------------------- ##
    # Intestinal Lumen 
    dAInL1 <- - ka*AInL1 + kHydro*AInL.glu1 - kTIn.3*AInL1 # Intestinal lumen 1
    dAInL2 <- - ka*AInL2 + kHydro*AInL.glu2 + kTIn.3*AInL1 - kTIn.3*AInL2   # Intestinal lumen 2
    dAInL3 <- - ka*AInL3 + kHydro*AInL.glu3 + kTIn.3*AInL2 - kTCo*AInL3   # Intestinal lumen 3
    
    # Intestinal tissue
    dAIn1 <-  + ka*AInL1 + QIn*(CBl-CvIn1) - VmaxMet.In*CIn1/(KmGluc.In+CIn1) - kTIn.3*AIn1   # Intestinal tissue 1
    dAIn2 <-  + ka*AInL2 + QIn*(CBl-CvIn2) - VmaxMet.In*CIn2/(KmGluc.In+CIn2) + kTIn.3*AIn1 - kTIn.3*AIn2 # Intestinal tissue 2
    dAIn3 <-  + ka*AInL3 + QIn*(CBl-CvIn3) - VmaxMet.In*CIn3/(KmGluc.In+CIn3) + kTIn.3*AIn2  # Intestinal tissue 3
    
    # Feces
    dAFe <- kTCo*AInL3
    
    # Liver
    dALi <- QIn*(CvIn1 + CvIn2 + CvIn3) + QHa*CBl - QLi*CvLi - VmaxMet.Li*CLi/(KmGluc.Li+CLi)
    
    # Kidney
    dAKi <- QKi*(CBl - CvKi) - CL_Re*CvKi*fub 
    
    # Urine 
    dAUr <- CL_Re*CvKi*fub 
    
    # Adipose
    dAAd <- QAd*(CBl - CvAd)
    
    # Slowly Perfused Tissue
    dASl <- QSl*(CBl - CvSl)
    
    # Rapidly Perfused Tissue
    dARa <- QRa*(CBl - CvRa)
    
    dABl <- QLi*CvLi + QKi*CvKi + QAd*CvAd + QSl*CvSl + QRa*CvRa - (QIn + QHa + QKi + QAd + QSl + QRa)*CBl
    
    
    ## Glucuronide---------------------- ##
    
    # Intestinal tissue 1
    dAIn.glu1 <- VmaxMet.In*CIn1/(KmGluc.In+CIn1) + QIn*(CBl.glu-CvIn.glu1) - AIn.glu1  - kTIn.3*AIn.glu1 # Intestinal tissue 1
    dAIn.glu2 <- VmaxMet.In*CIn2/(KmGluc.In+CIn2) + QIn*(CBl.glu-CvIn.glu2) - AIn.glu2  + kTIn.3*AIn.glu1 - kTIn.3*AIn.glu2 # Intestinal tissue 2
    dAIn.glu3 <- VmaxMet.In*CIn3/(KmGluc.In+CIn3) + QIn*(CBl.glu-CvIn.glu3) - AIn.glu3  + kTIn.3*AIn.glu2 # Intestinal tissue 3
    
    
    # Intestinal lumen 1
    dAInL.glu1 <- AIn.glu1 - kHydro*AInL.glu1 - kTIn.3*AInL.glu1 + kTge*ABi.glu # Intestinal Lumen 1
    dAInL.glu2 <- AIn.glu2 - kHydro*AInL.glu2 + kTIn.3*AInL.glu1 - kTIn.3*AInL.glu2 # Intestinal Lumen 2
    dAInL.glu3 <- AIn.glu3 - kHydro*AInL.glu3 + kTIn.3*AInL.glu2 - kTCo*AInL.glu3  # Intestinal Lumen 3
    
    # Feces
    dAFe.glu <- kTCo*AInL.glu3
    #           fromInLumen
    
    # Liver
    dALi.glu <- VmaxMet.Li*CLi/(KmGluc.Li+CLi) + QIn*(CvIn.glu1 + CvIn.glu2 + CvIn.glu3)  + QHa*CBl.glu - QLi*CvLi.glu - CL_Bi*CLi.glu
    #           fromZen                          fromPortal                                 fromHepaticA  toBlood        toBile
    
    # Bile
    dABi.glu <- CL_Bi*CLi.glu - kTge*ABi.glu
    #           fromLiver       toIntestLumen
    
    # Kidney
    dAKi.glu <- QKi*(CBl.glu-CvKi.glu) - CL_Re*CKi.glu
    #           fromToBlood              toUrine
    
    # Urine
    dAUr.glu <- CL_Re*CKi.glu
    #           fromKidney
    
    # Adipose
    dAAd.glu <- QAd*(CBl.glu - CvAd.glu)
    
    # Slowly Perfused Tissue
    dASl.glu <- QSl*(CBl.glu - CvSl.glu)
    
    # Rapidly Perfused Tissue
    dARa.glu <- QRa*(CBl.glu - CvRa.glu)
    
    # Blood
    dABl.glu <- QLi*CvLi.glu + QKi*CvKi.glu + QAd*CvAd.glu + QSl*CvSl.glu + QRa*CvRa.glu +
      - ((QIn + QHa + QKi + QAd + QSl + QRa) * CBl.glu)
    
    #----Output---------------------------------------------------------------#
    
    return(list(c(dAInL1, dAInL2, dAInL3, dAIn1, dAIn2, dAIn3, dAFe,
                  dALi, dAKi, dAUr, dAAd, dASl, dARa, dABl,
                  dAIn.glu1, dAIn.glu2, dAIn.glu3, 
                  dAInL.glu1, dAInL.glu2, dAInL.glu3, dAFe.glu,
                  dALi.glu, dABi.glu, dAKi.glu, dAUr.glu,
                  dAAd.glu, dASl.glu, dARa.glu, dABl.glu
    ), c(CInL1, #"30"
         CInL2, #"31"
         CInL3, #"32"
         CIn1, #"33"
         CIn2, #"34"
         CIn3, #"35"
         CLi, #"36"
         CKi, #"37"
         CAd, #"38"
         CSl, #"39"
         CRa, #"40"
         CBl, #"41"
         CIn.glu1, #"42"
         CIn.glu2, #"43"
         CIn.glu3, #"44"
         CInL.glu1, #"45"
         CInL.glu2, #"46"
         CInL.glu3, #"47"
         CLi.glu, #"48"
         CBi.glu, #"49"
         CKi.glu, #"50"
         CAd.glu, #"51"
         CSl.glu, #"52"
         CRa.glu, #"53"
         CBl.glu #"54"
    )))
  })
}



# -----------------------------------------------------------------------------
# RUN THE MODEL ####
# -----------------------------------------------------------------------------

PBK <- ode(y = INITIALstate,
           times = Tout,  #used to be times
           func = PBK_ZENandGLUC,
           parms = PARAMETERS,
           atol = 1e-50,
           #events = list(data = events), #for repeated exposure
           method = "lsoda")
head(PBK)



# -----------------------------------------------------------------------------
# RESULTS ####
# -----------------------------------------------------------------------------
# plot(PBK)
# ggsave(here(PlotFolderIs, "Overall.png"), dpi = 300)

# Output in amounts (umol)
OUTPUT.A <- as.data.frame(PBK) %>% select("time", "AInL1" : "ABl.glu")
write.csv(OUTPUT.A, file = here(TodaysOuput, "OUTPUT.A.csv"), row.names=FALSE)


## Mass Balance / Error ####
## ------------------------------------------------------

CALCULATED <- c( #!Change accordingly depending on the compartments in the function#
  "AInL1", "AInL2", "AInL3", "AIn1", "AIn2", "AIn3", "AFe",
  "ALi", "AKi", "AUr", "AAd", "ASl", "ARa", "ABl",
  "AIn.glu1", "AIn.glu2", "AIn.glu3", 
  "AInL.glu1", "AInL.glu2", "AInL.glu3", "AFe.glu",
  "ALi.glu", "ABi.glu", "AKi.glu", "AUr.glu",
  "AAd.glu", "ASl.glu", "ARa.glu", "ABl.glu"
  )

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
## ------------------------------------------------------

# Concentrations
OUTPUT.C <- as.data.frame(PBK) %>% select(!c("AInL1" : "ABl.glu")) %>%
  rename(
    "Time" = "time",
    "Intestinal_Lumen1" = "30",
    "Intestinal_Lumen2" = "31",
    "Intestinal_Lumen3" = "32",
    "Intestinal_Tissue1" = "33",
    "Intestinal_Tissue2" = "34",
    "Intestinal_Tissue3" = "35",
    "Liver" = "36",
    "Kidney" = "37",
    "Adipose" = "38",
    "Slowly_perfused" = "39",
    "Rapidly_perfused" = "40",
    "Blood" = "41",
    "Intestinal_Tissue_glucuronide1" = "42",
    "Intestinal_Tissue_glucuronide2" = "43",
    "Intestinal_Tissue_glucuronide3" = "44",
    "Intestinal_Lumen_glucuronide1" = "45",
    "Intestinal_Lumen_glucuronide2" = "46",
    "Intestinal_Lumen_glucuronide3" = "47",
    "Liver_glucuronide" = "48",
    "Bile_glucuronide" = "49",
    "Kidney_glucuronide" = "50",
    "Adipose_glucuronide" = "51",
    "Slowly_perfused_glucuronide" = "52",
    "Rapidly_perfused_glucuronide" = "53",
    "Blood_glucuronide" = "54"
  )

write.csv(OUTPUT.C, file = here(TodaysOuput, "OUTPUT.C.csv"), row.names=FALSE)


# To create df of simulated data for plotting
CONCENTRATIONS <- OUTPUT.C %>% #!change accordingly depending on the compartments in the model
  # rename(#!change accordingly depending on the compartments in the model
  #   "Time" = "time",
  #   "Intestine" = "CIn",
  #   "Liver" = "CLi",
  #   "Kidney" = "CKi",
  #   "Adipose" = "CAd",
  #   "Slowly_perfused" = "CSl",
  #   "Rapidly_perfused" = "CRa",
  #   "Blood" = "CBl"
  # ) %>%
  pivot_longer(names_to = "Tissue",
               values_to = "Concentration",
               Intestinal_Lumen1:Blood_glucuronide,
  ) %>%
  mutate(Concentration = Concentration * INPUT$MW)  #converting umol/L to ng/ml


# To plot Zearalenone blood concentration
CONCENTRATIONS %>%
  filter(str_starts(Tissue, "Blood")) %>%
 #        !str_ends(Tissue, "_glucuronide")) %>%
  ggplot(aes(Time, Concentration)) +
  facet_wrap(~Tissue) +
  geom_path(linewidth = 0.5) +
  labs(title = "Zearalenone blood concentration",
       x = "Time (h)",
       y = "Concentration (ng/ml)") +
  CP_theme()
ggsave(here(PlotFolderIs, "ZEN.blood.simulated.png"), dpi = 300)


# To plot Zearalenone and Zearalenone glucuronide blood concentrations
CONCENTRATIONS %>%
  filter(str_starts(Tissue, "Blood")) %>%
  ggplot(aes(Time, Concentration, colour = Tissue)) +
  facet_wrap(~Tissue) +
  geom_path(linewidth = 0.5) +
  labs(title = "Zearalenone and glucuronide tissue concentrations",
       x = "Time (h)",
       y = "Concentration (ng/ml)") +
  CP_theme()
ggsave(here(PlotFolderIs, "ZEN&Gluc.blood.simulated.png"), dpi = 300)


# To plot Zearalenone tissue concentrations
CONCENTRATIONS %>%
  filter(!str_ends(Tissue, "_glucuronide"),
         Tissue != "Intestinal_lumen",
         Tissue != "Intestinal_tissue"
         ) %>%
  ggplot(aes(Time, Concentration, colour = Tissue)) +
  facet_wrap(~Tissue) +
  geom_path(linewidth = 0.5) +
  labs(title = "Zearalenone tissue concentrations",
       x = "Time (h)",
       y = "Concentration (ng/ml)") +
  CP_theme()
ggsave(here(PlotFolderIs, "ZEN.simulated.png"), dpi = 300)

# To plot Zearalenone glucuronide tissue concentrations
CONCENTRATIONS %>%
  filter(str_ends(Tissue, "_glucuronide"),
         Tissue != "Intestinal_lumen_glucuronide",
         Tissue != "Intestinal_tissue_glucuronide"
  ) %>%
  ggplot(aes(Time, Concentration, colour = Tissue)) +
  geom_path(linewidth = 0.5) +
  facet_wrap(~Tissue) +
  labs(title = "Zearalenone glucuronide tissue concentrations",
       x = "Time (h)",
       y = "Concentration (ng/ml)") +
  CP_theme()
ggsave(here(PlotFolderIs, "ZEN_gluc.simulated.png"), dpi = 300)



## Plotting Experimental data ####
## ------------------------------------------------------

# Input in vivo data
DataInvivo <- read.csv(here("InputData/InVivoData.csv"), header=1) %>%  
  filter(Author == "Videmann", #I chose Videmann as it had the best recoveries, also Videmann separates Zearalenone from it's metabolites
         Gender != "Pregnant.F", #pregnancy not included yet
         TestCompound == "ZEN", #other options are ZEN_glu and a_ZEL
         AdminType == "PO", #"IV"), #Videmann doesn't have IV, but subcutaneous!!
         Tissue  %in% c("Plasma", "Intestine", "Kidney", "Liver")) %>%  #"Urine", "Colon"
  select(c("Gender", "SamplingTime", "Tissue", "Concentration", "SampleCompound"))%>%
  rename(Time = SamplingTime) %>% #(h)
  mutate(Concentration = Concentration/0.76)  #correcting for the recovery (recovery was 72-80% => averaged that to 76% and then corrected the concentrations for that)


## Plotting Simulated Vs Experimental data for Zearalenone ####
## ------------------------------------------------------

# Experimental data for only Zearalenone
POData.ZEN <- DataInvivo %>%
  filter(SampleCompound == "ZEN")

# Simulated data for Zearalenone
CONCENTRATIONS.ZEN <- CONCENTRATIONS %>%
  filter(!str_ends(Tissue, "_glucuronide"))

# Blood
InvivoBlood <- POData.ZEN %>% filter(Tissue == "Plasma") %>%
  mutate(Concentration = Concentration*ZEN.dat$BP) %>% #(ng/ml) blood concentration: changing from plasma concentrations to blood concentrations by multiplying with BP
  mutate(Tissue = recode(Tissue, "Plasma" = "Blood")) #renaming plasma to blood


SimulatedBlood <- CONCENTRATIONS.ZEN %>% filter(Tissue =="Blood")

SimulatedBlood %>% ggplot() +
  geom_line(aes(Time, Concentration)) +
  geom_point(data = InvivoBlood, aes(Time, Concentration, color = Gender)) +
  geom_line(data = InvivoBlood, aes(Time, Concentration, color = Gender),
            linewidth = 0.1) +
  labs(title = "Zearalenone blood concentrations",
       x = "Time (h)",
       y = "Concentration (ng/ml)",
       color = "Gender") +
  CP_theme()
ggsave(here(PlotFolderIs, "BloodSimulatedInvivoVidemann.png"), dpi = 300)

# Intestine
InvivoIntestine <- POData.ZEN %>% filter(Tissue =="Intestine")
SimulatedIntestine <- CONCENTRATIONS.ZEN %>% filter(Tissue =="Intestine")

SimulatedIntestine %>% ggplot() +
  geom_line(aes(Time, Concentration)) +
  geom_point(data = InvivoIntestine, aes(Time, Concentration, color = Gender)) +
  geom_line(data = InvivoIntestine, aes(Time, Concentration, color = Gender),
            linewidth = 0.1) +
  labs(title = "Zearalenone Intestine concentrations",
       x = "Time (h)",
       y = "Concentration (ng/ml)",
       color = "Gender") +
  CP_theme()
ggsave(here(PlotFolderIs, "IntestineSimulatedInvivo.png"), dpi = 300)

# Liver
InvivoLiver <- POData.ZEN %>% filter(Tissue =="Liver")
SimulatedLiver <- CONCENTRATIONS.ZEN %>% filter(Tissue =="Liver")

SimulatedLiver %>% ggplot() +
  geom_line(aes(Time, Concentration)) +
  geom_point(data = InvivoLiver, aes(Time, Concentration, color = Gender)) +
  geom_line(data = InvivoLiver, aes(Time, Concentration, color = Gender),
            linewidth = 0.1) +
  labs(title = "Zearalenone Liver concentrations",
       x = "Time (h)",
       y = "Concentration (ng/ml)",
       color = "Gender") +
  CP_theme()
ggsave(here(PlotFolderIs, "LiverSimulatedInvivo.png"), dpi = 300)

# Kidney
InvivoKidney <- POData.ZEN %>% filter(Tissue =="Kidney")
SimulatedKidney <- CONCENTRATIONS.ZEN %>% filter(Tissue =="Kidney")

SimulatedKidney %>% ggplot() +
  geom_line(aes(Time, Concentration)) +
  geom_point(data = InvivoKidney, aes(Time, Concentration, color = Gender)) +
  geom_line(data = InvivoKidney, aes(Time, Concentration, color = Gender),
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
POData.ZENGluc <- DataInvivo %>%
  filter(SampleCompound == "ZEN_glu")

# Simulated data for Zearalenone glucuronide
CONCENTRATIONS.ZENGluc <- CONCENTRATIONS %>%
  filter(str_ends(Tissue, "_glucuronide")) %>%
  pivot_wider(names_from = Tissue, values_from = Concentration) %>%
  rename("Intestine" = "Intestine_glucuronide",
         "Intestinal-lumen" = "Intestinal_lumen_glucuronide",
         "Liver" = "Liver_glucuronide",
         "Bile" = "Bile_glucuronide",
         "Kidney" = "Kidney_glucuronide",
         "Adipose" = "Adipose_glucuronide",
         "Slowly-Perfused" = "Slowly_perfused_glucuronide",
         "Rapidly-Perfused" = "Rapidly_perfused_glucuronide",
         "Blood" = "Blood_glucuronide") %>%
  pivot_longer(names_to = "Tissue", values_to = "Concentration", Intestine:Blood)

# Blood
InvivoBlood <- POData.ZENGluc %>% filter(Tissue == "Plasma") %>%
  #mutate(Concentration = Concentration*ZEN.dat$BP) %>% # assuming it's one as I don't have an experimental value
  mutate(Tissue = recode(Tissue, "Plasma" = "Blood")) #renaming plasma to blood


SimulatedBlood <- CONCENTRATIONS.ZENGluc %>% filter(Tissue =="Blood")

SimulatedBlood %>% ggplot() +
  geom_line(aes(Time, Concentration)) +
  geom_point(data = InvivoBlood, aes(Time, Concentration, color = Gender)) +
  geom_line(data = InvivoBlood, aes(Time, Concentration, color = Gender),
            linewidth = 0.1) +
  labs(title = "Zearalenone glucuronide blood concentrations",
       x = "Time (h)",
       y = "Concentration (ng/ml)",
       color = "Gender") +
  CP_theme()
ggsave(here(PlotFolderIs, "GlucBloodSimulatedInvivo.png"), dpi = 300)

# Intestine
InvivoIntestine <- POData.ZENGluc %>% filter(Tissue =="Intestine")
SimulatedIntestine <- CONCENTRATIONS.ZENGluc %>% filter(Tissue =="Intestine")

SimulatedIntestine %>% ggplot() +
  geom_line(aes(Time, Concentration)) +
  geom_point(data = InvivoIntestine, aes(Time, Concentration, color = Gender)) +
  geom_line(data = InvivoIntestine, aes(Time, Concentration, color = Gender),
            linewidth = 0.1) +
  labs(title = "Zearalenone glucuronide Intestine concentrations",
       x = "Time (h)",
       y = "Concentration (ng/ml)",
       color = "Gender") +
  CP_theme()
ggsave(here(PlotFolderIs, "GlucIntestineSimulatedInvivo.png"), dpi = 300)

# Liver
InvivoLiver <- POData.ZENGluc %>% filter(Tissue =="Liver")
SimulatedLiver <- CONCENTRATIONS.ZENGluc %>% filter(Tissue =="Liver")

SimulatedLiver %>% ggplot() +
  geom_line(aes(Time, Concentration)) +
  geom_point(data = InvivoLiver, aes(Time, Concentration, color = Gender)) +
  geom_line(data = InvivoLiver, aes(Time, Concentration, color = Gender),
            linewidth = 0.1) +
  labs(title = "Zearalenone glucuronide Liver concentrations",
       x = "Time (h)",
       y = "Concentration (ng/ml)",
       color = "Gender") +
  CP_theme()
ggsave(here(PlotFolderIs, "GlucLiverSimulatedInvivo.png"), dpi = 300)

# Kidney
InvivoKidney <- POData.ZENGluc %>% filter(Tissue =="Kidney")
SimulatedKidney <- CONCENTRATIONS.ZENGluc %>% filter(Tissue =="Kidney")

SimulatedKidney %>% ggplot() +
  geom_line(aes(Time, Concentration)) +
  geom_point(data = InvivoKidney, aes(Time, Concentration, color = Gender)) +
  geom_line(data = InvivoKidney, aes(Time, Concentration, color = Gender),
            linewidth = 0.1) +
  labs(title = "Zearalenone glucuronide Kidney concentrations",
       x = "Time (h)",
       y = "Concentration (ng/ml)",
       color = "Gender") +
  CP_theme()
ggsave(here(PlotFolderIs, "GlucKidneySimulatedInvivo.png"), dpi = 300)






