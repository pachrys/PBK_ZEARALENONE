
# STARTING UP ####

# Working Document
rm(list=ls()) # to clear out the global environment


# Libraries-------------------------------------------------------------------#

library(deSolve)
library(tidyverse)
library(ggplot2)
library(here)
library(showtext)
library(FME)
library(purrr)
font_add(family = "Garamond", regular = "GARA.TTF")
showtext_auto()

# Working directory-----------------------------------------------------------#
WhatAmITesting <- "TEST-Hypothesis"
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
  ABl.glu = 0,
  AInL1i = 0, 
  AInL2i = 0,
  AInL3i = 0,
  AIn1i = 0, 
  AIn2i = 0,
  AIn3i = 0,
  AFei = 0,
  ALii = 0,
  AKii = 0,
  AUri = 0,
  AAdi = 0,
  ASli = 0,
  ARai = 0,
  ABli = 0
)




# ---------------------------------------------------------------------------- #
# PBPK MODEL ####
# ---------------------------------------------------------------------------- #


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
    CvIn3.glu <- CIn3.glu / PInb.glu #Concentration leaving intestinal tissue
    
    CvLi.glu <- CLi.glu / PLib.glu #Concentration leaving liver
    CvKi.glu <- CKi.glu / PKib.glu #Concentration leaving kidney
    CvAd.glu <- CAd.glu / PAdb.glu #Concentration leaving adipose
    CvSl.glu <- CSl.glu / PSlb.glu #Concentration leaving slowly perfused tissue
    CvRa.glu <- CRa.glu / PRab.glu #Concentration leaving rapidly perfused tissue
    
    
    ### Zearalenone 2 ####
    
    ## Non perfused organ Concentrations (umol/L)
    CInL1i <- AInL1i / VInL1 #intestinal lumen
    CInL2i <- AInL2i / VInL2 #intestinal lumen
    CInL3i <- AInL3i / VInL3 #intestinal lumen
    # CBi <- ABi / VBi #bile
    
    ## Tissue cellular Concentrations (umol/L)                                    #Abbreviation (C for Concentration and then 2 first letters of tissue name)
    CIn1i <- AIn1i / VIn1 #intestine
    CIn2i <- AIn2i / VIn2 #intestine
    CIn3i <- AIn3i / VIn3 #intestine
    CLii <- ALii / VLi #liver
    CKii <- AKii / VKi #kidney
    CAdi <- AAdi / VAd #adipose
    CSli <- ASli / VSl #slowly perfused
    CRai <- ARai / VRa #rapidly perfused
    CBli <- ABli / VBl #blood
    
    ## Tissue venous Concentrations (umol/L)                                      #Only for perfused organs Equals tissue Concentration/tissue:blood partition coeefficient
    CvIn1i <- CIn1i / PInb #intestinal tissue
    CvIn2i <- CIn2i / PInb #intestinal tissue
    CvIn3i <- CIn3i / PInb #intestinal tissue
    CvLii <- CLii / PLib #liver
    CvKii <- CKii / PKib #kidney
    CvAdi <- CAdi / PAdb #adipose
    CvSli <- CSli / PSlb #slowly perfused
    CvRai <- CRai / PRab #rapidly perfused
    
    ## Tissue free Concentrations (umol/L), used for clearance/transport
    CKi_freei <- CvKii*fub
    CLi_freei <- CvLii*fub
    
    
    
    
    
    #----Differential equations-----------------------------------------------#
    # Rate of change in amounts in each tissue (umol/h)
    
    # Zearalenone-------------------- ##
    # Intestinal Lumen
    dAInL1 <- - ka*AInL1 - kTIn.3*AInL1 # Intestinal lumen 1 + kHydro*AInL1.glu
    dAInL2 <- - ka*AInL2 + kTIn.3*AInL1 - kTIn.3*AInL2 #+ kHydro*AInL1.glu    # Intestinal lumen 2 + kHydro*AInL2.glu
    dAInL3 <- - ka*AInL3 + kTIn.3*AInL2 - kTCo*AInL3 #+ kHydro*AInL.glu3  # Intestinal lumen 3 + kHydro*AInL3.glu
    
    # Feces
    dAFe <- kTCo*AInL3
    
    # Intestinal tissue
    dAIn1 <-  + ka*AInL1 + QIn*(CBl - CvIn1) - kTIn.3*AIn1  # Intestine 1 - VmaxMet.In*CIn1/(KmGluc.In+CIn1)
    dAIn2 <-  + ka*AInL2 + QIn*(CBl - CvIn2) + kTIn.3*AIn1 - kTIn.3*AIn2     # Intestine 2 - VmaxMet.In*CIn2/(KmGluc.In+CIn2)
    dAIn3 <-  + ka*AInL3 - VmaxMet.In*CIn3/(KmGluc.In+CIn3) + QIn*(CBl - CvIn3) + kTIn.3*AIn2    # Intestine 3  # Intestine 3
    
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
    dAIn1.glu <- - AIn1.glu + QIn*(CBl.glu - CvIn1.glu) - kTIn.3*AIn1.glu # Intestinal tissue 1 VmaxMet.In*CIn1/(KmGluc.In+CIn1)
    dAIn2.glu <- - AIn2.glu + QIn*(CBl.glu - CvIn2.glu) + kTIn.3*AIn1.glu - kTIn.3*AIn2.glu# Intestinal tissue 2 VmaxMet.In*CIn2/(KmGluc.In+CIn2)
    dAIn3.glu <- VmaxMet.In*CIn3/(KmGluc.In+CIn3) + VmaxMet.In*CIn3i/(KmGluc.In+CIn3i) - AIn3.glu + QIn*(CBl.glu - CvIn3.glu) + kTIn.3*AIn2.glu# Intestinal tissue 2  
    
    # Intestinal lumen
    dAInL1.glu <- AIn1.glu + kTge*ABi.glu - kTIn.3*AInL1.glu - kHydro*AInL1.glu # Intestinal Lumen 1
    dAInL2.glu <- AIn2.glu + kTIn.3*AInL1.glu - kTIn.3*AInL2.glu - kHydro*AInL2.glu # Intestinal Lumen 2
    dAInL3.glu <- AIn3.glu + kTIn.3*AInL2.glu - kTCo*AInL3.glu - kHydro*AInL3.glu # Intestinal Lumen 2
    
    
    # Feces
    dAFe.glu <- kTCo*AInL3.glu
    
    # Liver
    dALi.glu <- VmaxMet.Li*CLi/(KmGluc.Li+CLi) + VmaxMet.Li*CLii/(KmGluc.Li+CLii) + QIn*(CvIn1.glu+CvIn2.glu+CvIn3.glu) + QHa*CBl.glu - QLi*CvLi.glu - CL_Bi*CLi.glu
    
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
    
    
    
    
    ### Zearalenone 2-------------------- ####
    # Intestinal Lumen
    dAInL1i <- kHydro*AInL1.glu - kTIn.3*AInL1i # Intestinal lumen 1
    dAInL2i <- kHydro*AInL2.glu + kTIn.3*AInL1i - kTIn.3*AInL2i #+ kHydro*AInL1.glu    # Intestinal lumen 2
    dAInL3i <- kHydro*AInL3.glu + kTIn.3*AInL2i - kTCo*AInL3i #+ kHydro*AInL.glu3  # Intestinal lumen 3
    
    # Feces
    dAFei <- kTCo*AInL3i
    
    # Intestinal tissue
    dAIn1i <-  QIn*(CBli - CvIn1i) - kTIn.3*AIn1i  # Intestine 1 - VmaxMet.In*CIn1/(KmGluc.In+CIn1)
    dAIn2i <-  QIn*(CBli - CvIn2i) + kTIn.3*AIn1i - kTIn.3*AIn2i     # Intestine 2 - VmaxMet.In*CIn2/(KmGluc.In+CIn2)
    dAIn3i <-  - VmaxMet.In*CIn3i/(KmGluc.In+CIn3i) + QIn*(CBli - CvIn3i) + kTIn.3*AIn2i    # Intestine 3  # Intestine 3
    
    # Liver
    dALii <- QIn*(CvIn1i + CvIn2i + CvIn3i) + QHa*CBli - QLi*CvLii - VmaxMet.Li*CLii/(KmGluc.Li+CLii)
    #       fromIntest fromBlo   toBlood    toglucuronide
    # Kidney
    dAKii <- QKi*(CBli - CvKii) - CL_Re*CKi_freei #CL_Re is equal to rat GER
    #       fromtoblood        toUrine
    
    dAUri <- CL_Re*CKi_freei
    #       fromKidney
    # Adipose
    dAAdi <- QAd*(CBli - CvAdi)
    #       fromtoblood
    
    # Slowly Perfused Tissue
    dASli <- QSl*(CBli - CvSli)
    #       fromtoslowly
    
    # Rapidly Perfused Tissue
    dARai <- QRa*(CBli - CvRai)
    #       fromtorapidly
    
    dABli <- QLi*CvLii + QKi*CvKii + QAd*CvAdi + QSl*CvSli + QRa*CvRai +
      # from  liver      kidney     adipose    slowly     rapidly
      - ((3*QIn) + QHa + QKi + QAd + QSl + QRa)*CBli
    # to int  hepA  kid   adi   slo   rap
    
    
    
    
    
    
    
    
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
                  dABl.glu = dABl.glu,
                  dAInL1i = dAInL1, 
                  dAInL2i = dAInL2, 
                  dAInL3i = dAInL3,
                  dAIn1i = dAIn1, 
                  dAIn2i = dAIn2, 
                  dAIn3i = dAIn3, 
                  dAFei = dAFei,
                  dALii = dALii, 
                  dAKii = dAKii, 
                  dAUri = dAUri, 
                  dAAdi = dAAdi, 
                  dASli = dASli, 
                  dARai = dARai, 
                  dABli = dABli),
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
                       CBl.glu = CBl.glu,
                       CInL1i = CInL1i, 
                       CInL2i = CInL2i, 
                       CInL3i = CInL3i,
                       CIn1i = CIn1i, 
                       CIn2i = CIn2i,
                       CIn3i = CIn3i,
                       CLii = CLii,
                       CKii = CKii, 
                       CAdi = CAdi, 
                       CSli = CSli, 
                       CRai = CRai, 
                       CBli = CBli 
                       ))))
  })
  }


## SOLVE THE MODEL ####
## ------------------------------------------------------ #

# PBK <- ode(y = INITIALstate,
#              times = Tout,  #used to be times
#              func = PBK_ZENandGLUC,
#              parms = PARAMETERS,
#              atol = 1e-50,
#              #events = list(data = events), #for repeated exposure
#              method = "lsoda")
# head(PBK)


### RUN THE MODEL FOR DIFFERENT SETS OF PARAMETER CHANGES / OPTIMISATION ####
## ------------------------------------------------------ #

DEFAULT_PARAMS <- as.list(PARAMETERS)

VARYING.PARAMETERS <- tibble(
  PInb = c(INPUT$PIn, INPUT$PIn/10, INPUT$PIn/10, INPUT$PIn/10), #good prediction if deviding them /10
  PLib = c(INPUT$PLi, INPUT$PLi/10, INPUT$PLi/10, INPUT$PLi/10),
  PKib = c(INPUT$PKi, INPUT$PKi/10, INPUT$PKi/10, INPUT$PKi/10),
  PAdb = c(INPUT$PAd, INPUT$PAd/10, INPUT$PAd/10, INPUT$PAd/10),
  PSlb = c(INPUT$PSl, INPUT$PSl/10, INPUT$PSl/10, INPUT$PSl/10),
  PRab = c(INPUT$PRa, INPUT$PRa/10, INPUT$PRa/10, INPUT$PRa/10),
  PInb.glu = c(INPUT$PIn.glu , INPUT$PIn.glu*10, INPUT$PIn.glu*10, INPUT$PIn.glu), #good prediction if multiplying them *10; increasing the partition coeffs only doesn't help as it decreases ZEN's coeffs; increasing them while decreasing ZENs doesn't actually do anything
  PLib.glu = c(INPUT$PLi.glu, INPUT$PLi.glu*10, INPUT$PLi.glu*10, INPUT$PLi.glu),
  PKib.glu = c(INPUT$PKi.glu, INPUT$PKi.glu*10, INPUT$PKi.glu*10, INPUT$PKi.glu),
  PAdb.glu = c(INPUT$PAd.glu, INPUT$PAd.glu*10, INPUT$PAd.glu*10, INPUT$PAd.glu),
  PSlb.glu = c(INPUT$PSl.glu, INPUT$PSl.glu*10, INPUT$PSl.glu*10, INPUT$PSl.glu),
  PRab.glu = c(INPUT$PRa.glu, INPUT$PRa.glu*10, INPUT$PRa.glu*10, INPUT$PRa.glu),
  fub = c(INPUT$fub, INPUT$fub, INPUT$fub/10, INPUT$fub),
  ka = c(INPUT$ka, INPUT$ka, INPUT$ka, INPUT$ka), #same ka gives closer concentrations to the in vivo compared to increasing it
  kTge = c(INPUT$kTge, INPUT$kTge, INPUT$kTge, INPUT$kTge),
  CL_Bi = c(INPUT$CL_Bi, INPUT$CL_Bi/10, INPUT$CL_Bi/10, INPUT$CL_Bi/10), # divided /10 gives good prediction; shouldn't be increased as the concentration in the bile is increased
  KmGluc.In = c(INPUT$KmGluc.In, INPUT$KmGluc.In, INPUT$KmGluc.In, INPUT$KmGluc.In), #changing KmGluc.In doesn't make any significant difference
  VmaxMet.In = c(INPUT$VmaxMet.In, INPUT$VmaxMet.In, INPUT$VmaxMet.In, INPUT$VmaxMet.In), #changing VmaxMet.In doesn't make any significant difference
  KmGluc.Li = c(INPUT$KmGluc.Li, INPUT$KmGluc.Li, INPUT$KmGluc.Li, INPUT$KmGluc.Li), #decreasing doesn't help as it decreases the blood concentration too much; increasing it doesn't help as it brings the blood concentration to0 high for sure shouldn't be to increased
  VmaxMet.Li = c(INPUT$VmaxMet.Li, INPUT$VmaxMet.Li, INPUT$VmaxMet.Li, INPUT$VmaxMet.Li), #increasing it doesn't help as it decreases the blood concentration; decreasing it doesn't help as it brings the blood concentration to0 high
  kHydro = c(1000, 1000, 1000, 1000) #changing it doesn't do anything actually
)

combine_parameters <- function(varying_row) {
  modifyList(DEFAULT_PARAMS, as.list(varying_row))
}

# Create a list of parameter sets
parameter_sets <- VARYING.PARAMETERS %>%
  mutate(version = c("Original", "Test1", "Test2", "Test3")) %>%
  split(.$version) %>%
  map(~combine_parameters(.x %>% select(-version)))

print(parameter_sets)

# Function to run the model
run_model_for_parameters <- function(params) {
  PBK <- ode(y = INITIALstate, 
             times = Tout, 
             func = PBK_ZENandGLUC, 
             parms = params)
  return(PBK)
}

# Run the model for each set of parameters
model_results <- map(parameter_sets, run_model_for_parameters)
str(model_results)

print(model_results)

extract_results <- function(ode_result) {
  if (inherits(ode_result, "deSolve")) {
    # Convert the matrix to a data frame
    result_matrix <- ode_result
    result_df <- as.data.frame(result_matrix)
    
    # Check for column names in the attributes
    col_names <- attr(ode_result, "dimnames")[[2]]
    if (!is.null(col_names)) {
      colnames(result_df) <- col_names
    }
    
    return(result_df)
  } else {
    return(NULL)
  }
}

# Apply the function to extract results from all model results
extracted_resultsT <- map(model_results, extract_results)
str(extracted_resultsT)
Original.df <- extracted_resultsT[["Original"]] %>% mutate(Version = "Original")
Test1.df <- extracted_resultsT[["Test1"]] %>% mutate(Version = "Test1")
Test2.df <- extracted_resultsT[["Test2"]] %>% mutate(Version = "Test2")
Test3.df <- extracted_resultsT[["Test3"]] %>% mutate(Version = "Test3")

model_results_df <- rbind(Original.df, Test1.df, Test2.df, Test3.df)

# Save the combined results to a CSV file
write.csv(model_results_df, file = here(TodaysOuput, "OUTPUT_PARAMETERRANGE.csv"), row.names=FALSE)


# ---------------------------------------------------------------------------- #
# RESULTS ####
# ---------------------------------------------------------------------------- #
# plot(PBK)
# ggsave(here(PlotFolderIs, "Overall.png"), dpi = 300)

### Output in amounts (umol) ####
# OUTPUT.A <- as.data.frame(PBK) %>% select("time", "AInL1" : "ABl.glu")
# write.csv(OUTPUT.A, file = here(TodaysOuput, "OUTPUT.A.csv"), row.names=FALSE)

OUTPUT.A <- model_results_df %>% select("time", "AInL1" : "ABli", "Version")
write.csv(OUTPUT.A, file = here(TodaysOuput, "OUTPUT.A.csv"), row.names=FALSE)


## Mass Balance / Error ####
## ------------------------------------------------------ #

#CALCULATED <- rowSums(OUTPUT.A[ ,CALCULATED])
CALCULATEDnames <- names(OUTPUT.A)[-c(1, ncol(OUTPUT.A))]
TOTAL <- PO_DOSE
MB_ER <- OUTPUT.A %>%
  group_by(Version) %>%
  rowwise() %>%
  mutate(CALCULATED_SUM = sum(c_across(all_of(CALCULATEDnames)), na.rm = TRUE)) %>%
  mutate(MB = TOTAL - CALCULATED_SUM + 1) %>% 
  mutate(ER = (TOTAL-CALCULATED_SUM)/(TOTAL*100)) %>% 
  ungroup()

#MASSBALANCE <- TOTAL - CALCULATED + 1 #better be 1
#ERROR <- (TOTAL - CALCULATED) / (TOTAL) * 100 #better be 0
# MASSBALANCE_Error <- cbind("time" = OUTPUT.A[, "time"],
#                            MASSBALANCE, ERROR) %>%
#   as.data.frame()

# Mass Balance and Error plot
# MASSBALANCE_Error %>%
#   filter(time >= 0.01) %>% #as PO administration time is 0,01h
#   ggplot(aes(time, MASSBALANCE)) +
#   geom_line(color = "grey") +
#   geom_line(aes(time, ERROR), color = "grey6") +
#   labs(x = "Time (h)", y = "Value",
#        caption = "Mass balance in grey, Error in black") +
#   theme_CP()

MB_ER %>%
  filter(time >= 0.01) %>% #as PO administration time is 0,01h
  ggplot(aes(time, MB)) +
  geom_line(color = "grey") +
  geom_line(aes(time, ER), color = "grey6") +
  facet_wrap(~Version) +
  labs(x = "Time (h)", y = "Value",
       caption = "Mass balance in grey, Error in black") +
  theme_CP()
ggsave(here(TodaysOuput, "MB_ERROR2=.png"), dpi = 300)


## Plotting Simulated Output ####
## ------------------------------------------------------ #

# OUTPUT.C <- as.data.frame(PBK) %>% select(!c("AInL1" : "ABl.glu")) 
OUTPUT.C <- model_results_df %>% select(!c("AInL1" : "ABli")) 

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
    "Blood_glucuronide" = "CBl.glu",
    "Intestine_Lumen1i" = "CInL1i",
    "Intestine_Lumen2i" = "CInL2i",
    "Intestine_Lumen3i" = "CInL3i",
    "Intestine1i" = "CIn1i",
    "Intestine2i" = "CIn2i",
    "Intestine3i" = "CIn3i",
    "Liveri" = "CLii",
    "Kidneyi" = "CKii",
    "Adiposei" = "CAdi",
    "Slowly_perfusedi" = "CSli",
    "Rapidly_perfusedi" = "CRai",
    "Bloodi" = "CBli") %>% pivot_longer(names_to = "Tissue",
               values_to = "Concentration",
               Intestine_Lumen1:Bloodi) 


SIM_CONC.ZEN <- CONCENTRATIONS %>%
  filter(Tissue %in% c("Intestine_Lumen1",
                       "Intestine_Lumen2",
                       "Intestine_Lumen3",
                       "Intestine1",
                       "Intestine2",
                       "Intestine3",
                       "Liver",
                       "Kidney",
                       "Adipose",
                       "Slowly_perfused",
                       "Rapidly_perfused",
                       "Blood")) %>% 
 # filter(!str_ends(Tissue, "_glucuronide")) %>% 
  mutate(Concentration = Concentration * INPUT$MW) %>%  #converting umol/L to ng/ml
  pivot_wider(names_from = Tissue, values_from = Concentration) %>% 
  mutate(Intestine = Intestine_Lumen1 + Intestine_Lumen2 + Intestine_Lumen3 + Intestine1 + Intestine2 + Intestine3) %>% 
  select(!c("Intestine_Lumen1", "Intestine_Lumen2", "Intestine_Lumen3", "Intestine1", "Intestine2", "Intestine3"))%>%
  pivot_longer(names_to = "Tissue", values_to = "Concentration", Liver:Intestine)


SIM_CONC.ZENi <- CONCENTRATIONS %>%
  filter(Tissue %in% c("Intestine_Lumen1i",
                       "Intestine_Lumen2i",
                       "Intestine_Lumen3i",
                       "Intestine1i",
                       "Intestine2i",
                       "Intestine3i",
                       "Liveri",
                       "Kidneyi",
                       "Adiposei",
                       "Slowly_perfusedi",
                       "Rapidly_perfusedi",
                       "Bloodi")) %>% 
  # filter(!str_ends(Tissue, "_glucuronide")) %>% 
  mutate(Concentration = Concentration * INPUT$MW) %>%  #converting umol/L to ng/ml
  pivot_wider(names_from = Tissue, values_from = Concentration) %>% 
  mutate(Intestinei = Intestine_Lumen1i + Intestine_Lumen2i + Intestine_Lumen3i + Intestine1i + Intestine2i + Intestine3i) %>% 
  select(!c("Intestine_Lumen1i", "Intestine_Lumen2i", "Intestine_Lumen3i", "Intestine1i", "Intestine2i", "Intestine3i"))%>%
  pivot_longer(names_to = "Tissue", values_to = "Concentration", Liveri:Intestinei)


SIM_CONC.ZENGluc <- CONCENTRATIONS %>%
  filter(str_ends(Tissue, "_glucuronide")) %>%
  mutate(Concentration = Concentration * INPUT$MW.glu) %>% #converting umol/L to ng/ml 
  mutate(Tissue = str_replace(Tissue, "_glucuronide$", "")) %>% 
  pivot_wider(names_from = Tissue, values_from = Concentration) %>% 
  mutate(Intestine = Intestine_Lumen1 + Intestine_Lumen2 + Intestine_Lumen3 + Intestine1 + Intestine2 + Intestine3) %>% 
  select(!c("Intestine_Lumen1", "Intestine_Lumen2", "Intestine_Lumen3", "Intestine1", "Intestine2", "Intestine3")) %>%
  pivot_longer(names_to = "Tissue", values_to = "Concentration", Liver:Intestine)


# Data frame for Zearalenone glucuronide blood concentrations
SIM_BLOOD.ZEN <- SIM_CONC.ZEN %>% filter(Tissue == "Blood") %>% mutate(Compound = "Zearalenone") 
SIM_BLOOD.ZENi <- SIM_CONC.ZENi %>% filter(Tissue == "Bloodi") %>% mutate(Compound = "Zearalenonei") 
SIM_BLOOD.ZENGluc <- SIM_CONC.ZENGluc %>% filter(Tissue == "Blood") %>% mutate(Compound = "Z.Glucuronide")
SIM_BLOOD <- rbind(SIM_BLOOD.ZEN, SIM_BLOOD.ZENGluc, SIM_BLOOD.ZENi)

### Simulated data for Zearalenone ####
### ----------------------------------------------------- #

SIM_BLOOD %>%
  #filter(Version == "Original") %>% 
  ggplot(aes(Time, Concentration, colour = Compound)) +
  facet_wrap(~Version) +
  geom_path(linewidth = 0.5) +
  labs(title = "Zearalenone and glucuronide blood concentrations",
       x = "Time (h)",
       y = "Concentration (ng/ml)") +
  theme_CP()
ggsave(here(PlotFolderIs, "ZEN&Gluc.blood.simulated.png"), dpi = 300)


# To plot Zearalenone blood concentration
SIM_BLOOD.ZEN %>%
  ggplot(aes(Time, Concentration)) +
  facet_wrap(~Version) +
  geom_path(linewidth = 0.5) +
  labs(title = "Zearalenone blood concentration",
       x = "Time (h)",
       y = "Concentration (ng/ml)") +
  theme_CP()
ggsave(here(PlotFolderIs, "ZEN.blood.simulated.png"), dpi = 300)


# To plot Zearalenone i blood concentration
SIM_BLOOD.ZENi %>%
  ggplot(aes(Time, Concentration)) +
  facet_wrap(~Version) +
  geom_path(linewidth = 0.5) +
  labs(title = "Zearalenone i blood concentration",
       x = "Time (h)",
       y = "Concentration (ng/ml)") +
  theme_CP()
ggsave(here(PlotFolderIs, "ZENi.blood.simulated.png"), dpi = 300)



# To plot Zearalenone tissue concentrations
SIM_CONC.ZEN %>%
  # filter(Tissue != "Intestine") %>% 
  ggplot(aes(Time, Concentration, colour = Tissue)) +
  facet_wrap(~Version) +
  geom_path(linewidth = 0.5) +
  labs(title = "Zearalenone tissue concentrations",
       x = "Time (h)",
       y = "Concentration (ng/ml)") +
  theme_CP()
ggsave(here(PlotFolderIs, "ZEN.simulated.png"), dpi = 300)



### Simulated data for Zearalenone glucuronide ####
### ----------------------------------------------------- #

# To plot Zearalenone glucuronide tissue concentrations
SIM_CONC.ZENGluc %>%
  ggplot(aes(Time, Concentration, colour = Tissue)) +
  geom_path(linewidth = 0.5) +
  facet_wrap(~Version) +
  labs(title = "Zearalenone glucuronide tissue concentrations",
       x = "Time (h)",
       y = "Concentration (ng/ml)") +
  theme_CP()
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


### Simulated Vs Experimental data for Zearalenone ####
### ----------------------------------------------------- #

# Experimental data for only Zearalenone
EXP.PO_CONC.ZEN <- EXPERIMENTAL %>%
  filter(SampleCompound == "ZEN")

# Blood
EXP_BLOOD.ZEN <- EXP.PO_CONC.ZEN %>% filter(Tissue == "Plasma") %>%
  mutate(Concentration = Concentration*ZEN.dat$BP) %>% #(ng/ml) blood concentration: changing from plasma concentrations to blood concentrations by multiplying with BP
  mutate(Tissue = recode(Tissue, "Plasma" = "Blood")) #renaming plasma to blood


SIM_BLOOD.ZEN %>% ggplot() +
  geom_line(aes(Time, Concentration)) +
  facet_wrap(~Version) +
  geom_point(data = EXP_BLOOD.ZEN, aes(Time, Concentration, color = Gender)) +
  geom_line(data = EXP_BLOOD.ZEN, aes(Time, Concentration, color = Gender),
            linewidth = 0.1) +
  labs(title = "Zearalenone blood concentrations",
       x = "Time (h)",
       y = "Concentration (ng/ml)",
       color = "Gender") +
  theme_CP()
ggsave(here(PlotFolderIs, "BloodSimulatedInvivoVidemann.png"), dpi = 300)

# Intestine
EXP_INTESTINE <- EXP.PO_CONC.ZEN %>% filter(Tissue =="Intestine")
SIM_INTESTINE <- SIM_CONC.ZEN %>% filter(Tissue == "Intestine") %>% 
  filter(Time>=0.1) #%>% filter(Version != "Original")

SIM_INTESTINE %>% ggplot() +
  geom_line(aes(Time, Concentration)) +  
  facet_wrap(~Version) +
  geom_point(data = EXP_INTESTINE, aes(Time, Concentration, color = Gender)) +
  geom_line(data = EXP_INTESTINE, aes(Time, Concentration, color = Gender),
            linewidth = 0.1) +
  labs(title = "Zearalenone Intestine concentrations",
       x = "Time (h)",
       y = "Concentration (ng/ml)",
       color = "Gender") +
  theme_CP()
ggsave(here(PlotFolderIs, "IntestineSimulatedInvivo.png"), dpi = 300)

# Liver
EXP_LIVER <- EXP.PO_CONC.ZEN %>% filter(Tissue =="Liver")
SIM_LIVER <- SIM_CONC.ZEN %>% filter(Tissue =="Liver")

SIM_LIVER %>% ggplot() +
  geom_line(aes(Time, Concentration)) + 
  facet_wrap(~Version) +
  geom_point(data = EXP_LIVER, aes(Time, Concentration, color = Gender)) +
  geom_line(data = EXP_LIVER, aes(Time, Concentration, color = Gender),
            linewidth = 0.1) +
  labs(title = "Zearalenone Liver concentrations",
       x = "Time (h)",
       y = "Concentration (ng/ml)",
       color = "Gender") +
  theme_CP()
ggsave(here(PlotFolderIs, "LiverSimulatedInvivo.png"), dpi = 300)

# Kidney
EXP_KIDNEY <- EXP.PO_CONC.ZEN %>% filter(Tissue =="Kidney")
SIM_KIDNEY <- SIM_CONC.ZEN %>% filter(Tissue =="Kidney")

SIM_KIDNEY %>% ggplot() +
  geom_line(aes(Time, Concentration)) +
  facet_wrap(~Version) +
  geom_point(data = EXP_KIDNEY, aes(Time, Concentration, color = Gender)) +
  geom_line(data = EXP_KIDNEY, aes(Time, Concentration, color = Gender),
            linewidth = 0.1) +
  labs(title = "Zearalenone Kidney concentrations",
       x = "Time (h)",
       y = "Concentration (ng/ml)",
       color = "Gender") +
  theme_CP()
ggsave(here(PlotFolderIs, "KidneySimulatedInvivo.png"), dpi = 300)



### Simulated Vs Experimental data for Zearalenone Glucuronide ####
## ------------------------------------------------------

# Experimental data for only Zearalenone glucuronide
EXP.PO_CONC.ZENGluc <- EXPERIMENTAL %>%
  filter(SampleCompound == "ZEN_glu")

# Blood
EXP_BLOOD.ZENGluc <- EXP.PO_CONC.ZENGluc %>% filter(Tissue == "Plasma") %>%
  #mutate(Concentration = Concentration*ZEN.dat$BP) %>% # assuming it's one as I don't have an experimental value
  mutate(Tissue = recode(Tissue, "Plasma" = "Blood")) #renaming plasma to blood

SIM_BLOOD.ZENGluc %>% ggplot() +
  geom_line(aes(Time, Concentration)) +
  facet_wrap(~Version) +
  geom_point(data = EXP_BLOOD.ZENGluc, aes(Time, Concentration, color = Gender)) +
  geom_line(data = EXP_BLOOD.ZENGluc, aes(Time, Concentration, color = Gender),
            linewidth = 0.1) +
  labs(title = "Zearalenone glucuronide blood concentrations",
       x = "Time (h)",
       y = "Concentration (ng/ml)",
       color = "Gender") +
  theme_CP()
ggsave(here(PlotFolderIs, "GlucBloodSimulatedInvivo.png"), dpi = 300)

# Intestine
EXP_INTESTINE.ZENGluc <- EXP.PO_CONC.ZENGluc %>% filter(Tissue =="Intestine")
SIM_INTESTINE.ZENGluc <- SIM_CONC.ZENGluc %>% filter(Tissue =="Intestine")

SIM_INTESTINE.ZENGluc %>% ggplot() +
  geom_line(aes(Time, Concentration)) +
  facet_wrap(~Version) +
  geom_point(data = EXP_INTESTINE.ZENGluc, aes(Time, Concentration, color = Gender)) +
  geom_line(data = EXP_INTESTINE.ZENGluc, aes(Time, Concentration, color = Gender),
            linewidth = 0.1) +
  labs(title = "Zearalenone glucuronide Intestine concentrations",
       x = "Time (h)",
       y = "Concentration (ng/ml)",
       color = "Gender") +
  theme_CP()
ggsave(here(PlotFolderIs, "GlucIntestineSimulatedInvivo.png"), dpi = 300)

# Liver
EXP_LIVER.ZENGluc <- EXP.PO_CONC.ZENGluc %>% filter(Tissue =="Liver")
SIM_LIVER.ZENGluc <- SIM_CONC.ZENGluc %>% filter(Tissue =="Liver")

SIM_LIVER.ZENGluc %>% ggplot() +
  geom_line(aes(Time, Concentration)) +
  facet_wrap(~Version) +
  geom_point(data = EXP_LIVER.ZENGluc, aes(Time, Concentration, color = Gender)) +
  geom_line(data = EXP_LIVER.ZENGluc, aes(Time, Concentration, color = Gender),
            linewidth = 0.1) +
  labs(title = "Zearalenone glucuronide Liver concentrations",
       x = "Time (h)",
       y = "Concentration (ng/ml)",
       color = "Gender") +
  theme_CP()
ggsave(here(PlotFolderIs, "GlucLiverSimulatedInvivo.png"), dpi = 300)

# Kidney
EXP_KIDNEY.ZENGluc <- EXP.PO_CONC.ZENGluc %>% filter(Tissue =="Kidney")
SIM_KIDNEY.ZENGluc <- SIM_CONC.ZENGluc %>% filter(Tissue =="Kidney")

SIM_KIDNEY.ZENGluc %>% ggplot() +
  geom_line(aes(Time, Concentration)) +
  facet_wrap(~Version) +
  geom_point(data = EXP_KIDNEY.ZENGluc, aes(Time, Concentration, color = Gender)) +
  geom_line(data = EXP_KIDNEY.ZENGluc, aes(Time, Concentration, color = Gender),
            linewidth = 0.1) +
  labs(title = "Zearalenone glucuronide Kidney concentrations",
       x = "Time (h)",
       y = "Concentration (ng/ml)",
       color = "Gender") +
  theme_CP()
ggsave(here(PlotFolderIs, "GlucKidneySimulatedInvivo.png"), dpi = 300)



## Predicted vs Observed Regression ####

# Rearranging the data--------------------------------------------------------#
Version_range <- unique(OUTPUT.A$Version)  # Sequence where each value is 10 times higher than before Adjust the number of points as needed

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
}) %>% as.numeric()

SimVidemann.data <- SIM_BLOOD.ZEN %>% 
  filter(Time %in% closest_times) %>%  #simulated data is the same for both genders as gender differences are not included
  select(c(Time, Concentration, Version)) %>% 
  mutate(Time = as.numeric(Time), 
         Time = round(Time, 1))

ExpVidemann.data <- bind_rows(
  VidemannSimple.F.data %>% mutate(Gender = "Female"), 
  VidemannSimple.M.data %>% mutate(Gender = "Male")) %>% 
  mutate(Time = as.numeric(Time), 
         Time = round(Time, 1))

# A dataframe with all the possible combinations for Timepoints, gender and Version
Expanded.df <- expand.grid(
  Time = closest_times, 
  Gender = c("Female", "Male"), 
  Version = Version_range
) %>% mutate(Time = as.numeric(Time), Time = round(Time, 1))

Combined.df <- Expanded.df %>%
  left_join(ExpVidemann.data, by = c("Time", "Gender")) %>%
  left_join(SimVidemann.data, by = c("Time", "Version")) %>% 
  rename(Experimental = Concentration.x,
         Simulated = Concentration.y) %>% 
  arrange(Version, Gender, Time)


Results <- Combined.df %>%
  group_by(Version, Gender) %>%
  mutate(Time_Diff = c(NA, diff(Time))) %>%  # To calculate AUC; difference between consecutive "Time" values representing the width of each trapezoid, the sum sums up all the trapezoidal areas
  summarise(
    Videmann.AUC = sum((lead(Experimental) + Experimental) / 2 * Time_Diff, na.rm = TRUE), # To calculate AUC; lead() shifts the values forward by one time step every time; allowing to pair each value with the next one for the trapezoidal calculation
    SimVidemann.AUC = sum((lead(Simulated) + Simulated) / 2 * Time_Diff, na.rm = TRUE), # To calculate AUC; lead() shifts the values forward by one time step every time; allowing to pair each value with the next one for the trapezoidal calculation
    Regression = list(lm(Experimental ~ Simulated)) # To model experimental vs simulated regression for each Version
  )

# To view the results
Results %>%
  rowwise() %>%
  mutate(
    Model_Summary = list(summary(Regression[[1]]))
  ) %>%
  unnest_wider(Model_Summary)


# Plotting and further analysis
ggplot(Combined.df, aes(x = Simulated, y = Experimental, colour = Gender)) +
  geom_point() +
  geom_smooth(method = 'lm', color = "black") +
  facet_wrap(~ Version, scales = "free") +
  theme_CP()
  ggsave(here(PlotFolderIs, "ExpVsSimRegression.png"), dpi = 300)

write.csv(VARYING.PARAMETERS, file = here(TodaysOuput, "VARYING.PARAMETERS.csv"), row.names=FALSE)

# END ####
  