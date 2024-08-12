# ORIGINAL WITH ONE INTESTINE SCRIPT

# Libraries-------------------------------------------------------------------#

library(deSolve)
library(tidyverse)
library(ggplot2)
library(here)
library(showtext)
library(FME)
font_add(family = "Garamond", regular = "GARA.TTF")
showtext_auto()

rm(list=ls()) # to clear out the global environment

# Working directory-----------------------------------------------------------#
WhatAmITesting <- "SimplestModel"
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
             kTsi = 1 / ZEN.dat$kTsi, 
             kTco = 1 / ZEN.dat$kTco,
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
  ASt = PO_DOSE, 
  AIn = 0, 
  ACo = 0,
  ALi = 0,
  AKi = 0,
  AAd = 0,
  ASl = 0,
  ARa = 0,
  ABl = 0,
  AFe = 0,
  AUr = 0,
  ALi.glu = 0
)


# -----------------------------------------------------------------------------
# PBPK MODEL ####
# -----------------------------------------------------------------------------
# source(here("Script/ZENPBK_Function.R"))

PBK_ZEN <- function(t,state,parameters){
  with(as.list(c(state,parameters)),{
    
    #----Concentrations-------------------------------------------------------#
    
    # Zearalenone
    
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
    
    ## Tissue cellular Concentrations (umol/L)                                    #Abbreviation (C for Concentration and then 2 first letters of tissue name)
    CLi.glu <- ALi.glu / VLi #liver
    
    ## Tissue venous Concentrations (umol/L)                                      #Only for perfused organs Equals tissue Concentration/tissue:blood partition LInefficient
    CvLi.glu <- CLi.glu / PLib.glu #Concentration leaving liver
    
    
    #----Differential equations-----------------------------------------------#
    # Rate of change in amounts in each tissue (umol/h)
    
    # Zearalenone-------------------- ##
    # Gastro-intestinal tract - not perfused 
    dASt <- - kTge*ASt                                     # stomach 
    dAIn <-  + kTge*ASt + CL_Bi*CLi.glu - ka*AIn - kTsi*AIn  # intestine
    dACo <- + kTsi*AIn - kTco*ACo                                # colon
    
    # Perfused organs
    dALi <- + ka*AIn + QLi*(CBl - CvLi) - VmaxMet.Li*CLi/(KmGluc.Li+CLi) # liver
    dAKi <- QKi*(CBl - CvKi) - CL_Re*CKi_free              # kidney 
    dAAd <- QAd*(CBl - CvAd)                               # adipose
    dASl <- QSl*(CBl - CvSl)                               # slowly perfused
    dARa <- QRa*(CBl - CvRa)                               # rapidly perfused
    dABl <- QLi*CvLi + QKi*CvKi + QAd*CvAd + QSl*CvSl + QRa*CvRa - (QLi + QKi + QAd + QSl + QRa)*CBl # blood
    
    # Excretion compartments
    dAFe <- kTco*ACo                                       # feces
    dAUr <- CL_Re*CKi_free                                 # urine
    
    
    # Glucuronide---------------------- ##
    dALi.glu <- VmaxMet.Li*CLi/(KmGluc.Li+CLi) - CL_Bi*CLi.glu # liver
    
    
    #----Output---------------------------------------------------------------#
    
    return(list(c(dASt = dASt, 
                  dAIn = dAIn, 
                  dACo = dACo,
                  dALi = dALi, 
                  dAKi = dAKi, 
                  dAAd = dAAd, 
                  dASl = dASl, 
                  dARa = dARa, 
                  dABl = dABl,
                  dAFe = dAFe,
                  dAUr = dAUr, 
                  dALi.glu = dALi.glu),
                list(c(CIn = CIn, 
                       CLi = CLi ,
                       CKi = CKi, 
                       CAd = CAd, 
                       CSl = CSl, 
                       CRa = CRa, 
                       CBl = CBl, 
                       CLi.glu = CLi.glu))
    ))
  })
}


# -----------------------------------------------------------------------------
# RUN THE MODEL ####
# -----------------------------------------------------------------------------


PBK <- ode(y = INITIALstate,
           times = Tout,  #used to be times
           func = PBK_ZEN,
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
OUTPUT.A <- as.data.frame(PBK) %>% select("time", "ASt" : "ALi.glu")
write.csv(OUTPUT.A, file = here(TodaysOuput, "OUTPUT.A.csv"), row.names=FALSE)


## Mass Balance / Error ####
## ------------------------------------------------------

CALCULATED <- c( #!Change accordingly depending on the compartments in the function#
  "ASt", "AIn", "ACo",
  "ALi", "AKi", "AAd", "ASl", "ARa", "ABl",
  "AFe", "AUr",
  "ALi.glu")


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
OUTPUT.C <- as.data.frame(PBK) %>% select(!c("ASt" : "ALi.glu"))

write.csv(OUTPUT.C, file = here(TodaysOuput, "OUTPUT.C.csv"), row.names=FALSE)


# To create df of simulated data for plotting
CONCENTRATIONS <- OUTPUT.C %>% #!change accordingly depending on the compartments in the model
  rename(#!change accordingly depending on the compartments in the model
    "Time" = "time",
    "Intestine" = "CIn",
    "Liver" = "CLi",
    "Kidney" = "CKi",
    "Adipose" = "CAd",
    "Slowly_perfused" = "CSl",
    "Rapidly_perfused" = "CRa",
    "Blood" = "CBl",
    "Liver_glucuronide" = "CLi.glu"
  ) %>% pivot_longer(names_to = "Tissue",
                     values_to = "Concentration",
                     Intestine:Liver_glucuronide)

# Simulated data for Zearalenone
CONCENTRATIONS.ZEN <- CONCENTRATIONS %>%
  filter(!str_ends(Tissue, "_glucuronide")) %>%
  mutate(Concentration = Concentration * INPUT$MW)  #converting umol/L to ng/ml

# Simulated data for Zearalenone glucuronide
CONCENTRATIONS.ZENGluc <- CONCENTRATIONS %>%
  filter(str_ends(Tissue, "_glucuronide")) %>%
  mutate(Concentration = Concentration * INPUT$MW.glu) %>% #converting umol/L to ng/ml
  mutate(Tissue = str_replace(Tissue, "_glucuronide$", ""))



# To plot Zearalenone blood concentration
CONCENTRATIONS.ZEN %>%
  filter(Tissue == "Blood") %>%
  ggplot(aes(Time, Concentration)) +
  geom_path(linewidth = 0.5) +
  labs(title = "Zearalenone blood concentration",
       x = "Time (h)",
       y = "Concentration (ng/ml)") +
  CP_theme()
ggsave(here(PlotFolderIs, "ZEN.blood.simulated.png"), dpi = 300)


# To plot Zearalenone tissue concentrations
CONCENTRATIONS.ZEN %>%
  # filter(Tissue != "Intestine") %>%
  ggplot(aes(Time, Concentration, colour = Tissue)) +
  facet_wrap(~Tissue) +
  geom_path(linewidth = 0.5) +
  labs(title = "Zearalenone tissue concentrations",
       x = "Time (h)",
       y = "Concentration (ng/ml)") +
  CP_theme()
ggsave(here(PlotFolderIs, "ZEN.simulatedtissue.png"), dpi = 300)

# To plot Zearalenone tissue concentrations excluding intestine
CONCENTRATIONS.ZEN %>%
  filter(Tissue != "Intestine") %>%
  ggplot(aes(Time, Concentration, colour = Tissue)) +
  facet_wrap(~Tissue) +
  geom_path(linewidth = 0.5) +
  labs(title = "Zearalenone tissue concentrations",
       x = "Time (h)",
       y = "Concentration (ng/ml)") +
  CP_theme()
ggsave(here(PlotFolderIs, "ZEN.simulatedtissueNoInt.png"), dpi = 300)


# To plot Zearalenone and Zearalenone glucuronide blood concentrations
BLOOD_CONCENTRATIONS.ZEN <- CONCENTRATIONS.ZEN %>% filter(Tissue == "Blood") %>% mutate(Compound = "Zearalenone")

BLOOD_CONCENTRATIONS.ZEN %>%
  ggplot(aes(Time, Concentration)) +
  geom_path(linewidth = 0.5) +
  labs(title = "Zearalenone blood concentrations",
       x = "Time (h)",
       y = "Concentration (ng/ml)") +
  CP_theme()
ggsave(here(PlotFolderIs, "ZEN.blood.simulated.png"), dpi = 300)


# To plot Zearalenone glucuronide tissue concentrations
CONCENTRATIONS.ZENGluc %>%
  # filter(Tissue != "Bile") %>%
  ggplot(aes(Time, Concentration)) + #, colour = Tissue)) +
  geom_path(linewidth = 0.5) +
  labs(title = "ZENGlucuronide concentrations",
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
  select(c("Gender", "SamplingTime", "Tissue", "Concentration", "SampleCompound")) %>%
  rename(Time = SamplingTime) %>% #(h)
  mutate(Concentration = Concentration/0.76)  #correcting for the recovery (recovery was 72-80% => averaged that to 76% and then corrected the concentrations for that)

## Plotting Simulated Vs Experimental data for Zearalenone ####
## ------------------------------------------------------

# Experimental data for only Zearalenone
POData.ZEN <- DataInvivo %>%
  filter(SampleCompound == "ZEN")

# Blood
InvivoBlood <- POData.ZEN %>% filter(Tissue == "Plasma") %>%
  mutate(Concentration = Concentration*ZEN.dat$BP) %>% #(ng/ml) blood concentration: changing from plasma concentrations to blood concentrations by multiplying with BP
  mutate(Tissue = recode(Tissue, "Plasma" = "Blood")) #renaming plasma to blood


BLOOD_CONCENTRATIONS.ZEN %>% ggplot() +
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
SimulatedIntestine <- CONCENTRATIONS.ZEN %>% filter(Tissue == "Intestine") %>%
  filter(Time>=0.1)

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

