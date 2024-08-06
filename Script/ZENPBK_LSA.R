# SENSITIVITY ANALYSIS SCRIPT
# BY: CHRYSANTHI PACHOULIDE
# DATE: AUGUST 2024
# For background on the FME package: https://www.rdocumentation.org/packages/FME/versions/1.3.6.3/topics/sensFun

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
             # kTsi = 1 / ZEN.dat$kTsi, 
             # kTco = 1 / ZEN.dat$kTco,
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
  AInL = PO_DOSE, 
  AIn = 0, 
  AFe = 0,
  ALi = 0,
  AKi = 0,
  AUr = 0,
  AAd = 0,
  ASl = 0,
  ARa = 0,
  ABl = 0,
  AIn.glu = 0,
  AInL.glu = 0,
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
    CInL <- AInL / VInL #intestinal lumen
    # CBi <- ABi / VBi #bile

    ## Tissue cellular Concentrations (umol/L)                                    #Abbreviation (C for Concentration and then 2 first letters of tissue name)
    CIn <- AIn / VIn #intestine
    CLi <- ALi / VLi #liver
    CKi <- AKi / VKi #kidney
    CAd <- AAd / VAd #adipose
    CSl <- ASl / VSl #slowly perfused
    CRa <- ARa / VRa #rapidly perfused
    CBl <- ABl / VBl #blood

    ## Tissue venous Concentrations (umol/L)                                      #Only for perfused organs Equals tissue Concentration/tissue:blood partition coeefficient
    CvIn <- CIn / PInb #intestinal tissue
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
    CInL.glu <- AInL.glu / VInL #intestinal lumen glucuronide
    CBi.glu <- ABi.glu / VBi #bile

    ## Tissue cellular Concentrations (umol/L)                                    #Abbreviation (C for Concentration and then 2 first letters of tissue name)
    CIn.glu <- AIn.glu / VIn #intestinal tissue
    CLi.glu <- ALi.glu / VLi #liver
    CKi.glu <- AKi.glu / VKi #kidney
    CAd.glu <- AAd.glu / VAd #adipose tissue
    CSl.glu <- ASl.glu / VSl #slowly perfused tissue
    CRa.glu <- ARa.glu / VRa #rapidly perfused tissue
    CBl.glu <- ABl.glu / VBl #blood

    ## Tissue venous Concentrations (umol/L)                                      #Only for perfused organs Equals tissue Concentration/tissue:blood partition LInefficient
    CvIn.glu <- CIn.glu / PInb.glu #Concentration leaving intestinal tissue
    CvLi.glu <- CLi.glu / PLib.glu #Concentration leaving liver
    CvKi.glu <- CKi.glu / PKib.glu #Concentration leaving kidney
    CvAd.glu <- CAd.glu / PAdb.glu #Concentration leaving adipose
    CvSl.glu <- CSl.glu / PSlb.glu #Concentration leaving slowly perfused tissue
    CvRa.glu <- CRa.glu / PRab.glu #Concentration leaving rapidly perfused tissue


    #----Differential equations-----------------------------------------------#
    # Rate of change in amounts in each tissue (umol/h)

    # Zearalenone-------------------- ##
    # Intestinal Lumen
    dAInL <- - ka*AInL - kTIn*AInL + kHydro*AInL.glu #- VmaxAbs*CInL/(KmAbs+CInL) + VmaxSec*CIn/(KmSec+CIn)
    #          toInt

    # Intestinal tissue
    dAIn <-  + ka*AInL - VmaxMet.In*CIn/(KmGluc.In+CIn) + QIn*(CBl - CvIn)   #+ VmaxAbs*CInL/(KmAbs+CInL) - VmaxSec*CIn/(KmSec+CIn)
    #        fromLumenActive      toglucuronide           fromBlood toLiver             toLumen                   fromLum

    # # Feces
    dAFe <- kTIn*AInL
    #       fromIntestine

    # Liver
    dALi <- QIn*CvIn + QHa*CBl - QLi*CvLi - VmaxMet.Li*CLi/(KmGluc.Li+CLi)
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
      - (QIn + QHa + QKi + QAd + QSl + QRa)*CBl
    # to int  hepA  kid   adi   slo   rap

    # ## Glucuronide---------------------- ##

    # Intestine
    dAIn.glu <- VmaxMet.In*CIn/(KmGluc.In+CIn) - AIn.glu + QIn*CBl.glu - QIn*CvIn.glu # -VmaxExc.In.glu*CIn.glu/(KmGluc.Exc+CIn.glu) + VmaxAbs.In.glu*CInL.glu/(KmGluc.Abs+CInL.glu)  #(umol/h)
    #           fromIntZen                       toInL     fromblood     toportal                                        fromInL

    # Intestinal lumen
    dAInL.glu <- AIn.glu - kHydro*AInL.glu - kTIn*AInL.glu + kTge*ABi.glu #+ VmaxExc.In.glu*CIn.glu/(KmGluc.Exc+CIn.glu) - VmaxAbs.In.glu*CInL.glu/(KmGluc.Abs+CInL.glu)  #(umol/h)
    #            fromInTis toZen             toFeces         fromBile

    # Feces
    dAFe.glu <- kTIn*AInL.glu
    #           fromInLumen

    # Liver
    dALi.glu <- VmaxMet.Li*CLi/(KmGluc.Li+CLi) + QIn*CvIn.glu + QHa*CBl.glu - QLi*CvLi.glu - CL_Bi*CLi.glu
    #           fromZen                          fromPortal     fromHepaticA  toBlood        toBile

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

    return(list(c(dAInL, dAIn, dAFe,
                  dALi, dAKi, dAUr, dAAd, dASl, dARa, dABl,
                  dAIn.glu, dAInL.glu, dAFe.glu,
                  dALi.glu, dABi.glu, dAKi.glu, dAUr.glu,
                  dAAd.glu, dASl.glu, dARa.glu, dABl.glu
    ), c(CInL, #"22"
         CIn, #"23"
         CLi, #"24"
         CKi, #"25"
         CAd, #"26"
         CSl, #"27"
         CRa, #"28"
         CBl, #"29"
         CIn.glu, #"30"
         CInL.glu, #"31"
         CLi.glu, #"32"
         CBi.glu, #"33"
         CKi.glu, #"34"
         CAd.glu, #"35"
         CSl.glu, #"36"
         CRa.glu, #"37"
         CBl.glu #"38"
    )))
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

DF_summary_SENSI_all <- summary(SENSI_all) %>% as.data.frame()
Plot_DF_summary_SENSI_all <- DF_summary_SENSI_all %>% select(Mean) %>%
  mutate(Parameter = rownames(DF_summary_SENSI_all)) %>%
  ggplot(aes(x = Mean, y = Parameter)) +
  geom_bar(stat = "identity") +
  coord_flip() + # Flip coordinates for better readability if there are many parameters
  labs(x = "Mean Sensitivity coefficient", y = "Parameter") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis text 45 degrees
  )

Plot_DF_summary_SENSI_all


## Analysing the sensitivity on only the Blood concentration of Zearalenone ####
## ------------------------------------------------------

SENSI_CBl <- SENSI_all %>% filter(var == "29")
summary(SENSI_CBl)

DF_summary_SENSI_CBl <- summary(SENSI_CBl) %>% as.data.frame()
Plot_DF_summary_SENSI_CBl <- DF_summary_SENSI_CBl %>% select(Mean) %>%
  mutate(Parameter = rownames(DF_summary_SENSI_CBl)) %>%
  ggplot(aes(x = Mean, y = Parameter)) +
  geom_bar(stat = "identity") +
  coord_flip() + # Flip coordinates for better readability if there are many parameters
  labs(x = "Mean Sensitivity coefficient", y = "Parameter") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis text 45 degrees
  )

Plot_DF_summary_SENSI_CBl



## Analysing the sensitivity on only the Blood concentration of Zearalenone glucuronide ####
## ------------------------------------------------------

SENSI_CBl.glu <- SENSI_all %>% filter(var == "38")
summary(SENSI_CBl.glu)

DF_summary_SENSI_CBl.glu <- summary(SENSI_CBl.glu) %>% as.data.frame()
Plot_DF_summary_SENSI_CBl.glu <- DF_summary_SENSI_CBl.glu %>% select(Mean) %>%
  mutate(Parameter = rownames(DF_summary_SENSI_CBl.glu)) %>%
  ggplot(aes(x = Mean, y = Parameter)) +
  geom_bar(stat = "identity") +
  coord_flip() + # Flip coordinates for better readability if there are many parameters
  labs(x = "Mean Sensitivity coefficient", y = "Parameter") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis text 45 degrees
  )

Plot_DF_summary_SENSI_CBl.glu

