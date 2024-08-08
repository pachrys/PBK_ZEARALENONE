# SENSITIVITY ANALYSIS SCRIPT
# BY: CHRYSANTHI PACHOULIDE
# DATE: AUGUST 2024
# For background on the FME package: https://www.rdocumentation.org/packages/FME/versions/1.3.6.3/topics/sensFun

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
WhatAmITesting <- "AddingLSA3"
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


# -----------------------------------------------------------------------------
# LOCAL SENSITIVITY ANALYSIS ####
# -----------------------------------------------------------------------------

pars <- PARAMETERS %>%  as.list()

SENSI_SOLVE <- function(pars) {
  
  PBK_ZENandGLUC
  
  state   <- INITIALstate
  tout    <- Tout
  
  out <- ode(y = state, parms = pars, times = Tout, func = PBK_ZENandGLUC)
  
  return(as.data.frame(out))
}

SENSI_OUT <- SENSI_SOLVE(pars) #saving the results of the model used for Sensitivity Analysis

#All parameters
SENSI_all <- sensFun(
  func = SENSI_SOLVE,
  parms = pars, 
  tiny = 0.1)
head(SENSI_all)

## Analysing the sensitivity on all the variables ####
## ------------------------------------------------------

summary(SENSI_all)

# value : the nominal value of the parameter
# scale : the scale of the parameter
# L1 : the L1 norm of each parameter's sensitivity function
# ---> it's the average of the absolute values of the sensitivities; represents the sensitivity of the x parameter in the y timepoint
# L2 : the L2 norm of each parameter's sensitivity function
# ---> it's the square root of the average of the squares of the sensitivities; represents the sensitivity of the x parameter at the y timepoint
# Mean : the mean value of the sensitivity function
# Min : the minimum value of the sensitivity function
# Max : the maximum value of the sensitivity function
# N : the number of points used in the sensitivity analysis

# L1 and L2 norms provide a measure of the overall sensitivity of the model output to changes in each parameter

DF_summary_SENSI_all <- summary(SENSI_all) %>% as.data.frame() %>% 
  select(Mean) %>%
  mutate(Parameter = rownames(DF_summary_SENSI_all))
  
DF_summary_SENSI_all %>% 
  ggplot(aes(x = Mean, y = Parameter)) +
  geom_bar(stat = "identity") +
  coord_flip() + # Flip coordinates for better readability if there are many parameters
  labs(x = "Mean Sensitivity coefficient", y = "Parameter") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis text 45 degrees
  )
ggsave(here(TodaysOuput, "SENSI_all.png"), dpi = 300)


## Analysing the sensitivity on only the Blood concentration of Zearalenone ####
## ------------------------------------------------------

SENSI_CBl <- SENSI_all %>% filter(var == "CBl")
summary(SENSI_CBl)

DF_summary_SENSI_CBl <- summary(SENSI_CBl) %>% as.data.frame() %>% 
  select(Mean) %>%
  mutate(Parameter = rownames(DF_summary_SENSI_CBl))

DF_summary_SENSI_CBl %>% 
  ggplot(aes(x = Mean, y = Parameter)) +
  geom_bar(stat = "identity") +
  coord_flip() + # Flip coordinates for better readability if there are many parameters
  labs(x = "Mean Sensitivity coefficient", y = "Parameter") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis text 45 degrees
  ) 
ggsave(here(TodaysOuput, "SENSI_CBl.png"), dpi = 300)



## Analysing the sensitivity on only the Blood concentration of Zearalenone glucuronide ####
## ------------------------------------------------------

SENSI_CBl.glu <- SENSI_all %>% filter(var == "CBl.glu")
summary(SENSI_CBl.glu)

DF_summary_SENSI_CBl.glu <- summary(SENSI_CBl.glu) %>% as.data.frame() %>% 
  select(Mean) %>%
  mutate(Parameter = rownames(DF_summary_SENSI_CBl.glu))

DF_summary_SENSI_CBl.glu %>% 
  ggplot(aes(x = Mean, y = Parameter)) +
  geom_bar(stat = "identity") +
  coord_flip() + # Flip coordinates for better readability if there are many parameters
  labs(x = "Mean Sensitivity coefficient", y = "Parameter") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis text 45 degrees
  )
ggsave(here(TodaysOuput, "SENSI_CBl.gluc.png"), dpi = 300)

