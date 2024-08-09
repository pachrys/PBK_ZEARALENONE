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
WhatAmITesting <- "ChangingUptakeToLiver"
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
    dAInL <- - AInL - kTIn*AInL + kHydro*AInL.glu # Assuming everything gets passively absorbed to the tissue; -ka*AInL 
    
    # Intestinal tissue
    dAIn <-  - ka*AIn + AInL - VmaxMet.In*CIn/(KmGluc.In+CIn) + QIn*(CBl - CvIn)   # Assuming ZEN also leaves from the intestine to the Liver via Papp
    
    # # Feces
    dAFe <- kTIn*AInL
    #       fromIntestine
    
    # Liver
    dALi <- + ka*AIn + QIn*CvIn + QHa*CBl - QLi*CvLi - VmaxMet.Li*CLi/(KmGluc.Li+CLi)
    
    # Kidney
    dAKi <- QKi*(CBl - CvKi) - CL_Re*CKi_free #CL_Re is equal to rat GER
    #       fromtoblood        toUrine
    
    dAUr <- CL_Re*CKi_free
    
    # Adipose
    dAAd <- QAd*(CBl - CvAd)
    
    # Slowly Perfused Tissue
    dASl <- QSl*(CBl - CvSl)
    
    # Rapidly Perfused Tissue
    dARa <- QRa*(CBl - CvRa)
    
    dABl <- QLi*CvLi + QKi*CvKi + QAd*CvAd + QSl*CvSl + QRa*CvRa +
      - (QIn + QHa + QKi + QAd + QSl + QRa)*CBl
    
    # ## Glucuronide---------------------- ##
    
    # Intestine
    dAIn.glu <- VmaxMet.In*CIn/(KmGluc.In+CIn) - AIn.glu + QIn*CBl.glu - QIn*CvIn.glu 
    
    # Intestinal lumen
    dAInL.glu <- AIn.glu - kHydro*AInL.glu - kTIn*AInL.glu + kTge*ABi.glu 
    
    # Feces
    dAFe.glu <- kTIn*AInL.glu
    
    # Liver
    dALi.glu <- VmaxMet.Li*CLi/(KmGluc.Li+CLi) + QIn*CvIn.glu + QHa*CBl.glu - QLi*CvLi.glu - CL_Bi*CLi.glu
    
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
      - ((QIn + QHa + QKi + QAd + QSl + QRa) * CBl.glu)
    
    
    #----Output---------------------------------------------------------------#
    
    
    return(list(c(dAInL = dAInL, 
                  dAIn = dAIn, 
                  dAFe = dAFe,
                  dALi = dALi, 
                  dAKi = dAKi, 
                  dAUr = dAUr, 
                  dAAd = dAAd, 
                  dASl = dASl, 
                  dARa = dARa, 
                  dABl = dABl,
                  dAIn.glu = dAIn.glu,
                  dAInL.glu = dAInL.glu, 
                  dAFe.glu = dAFe.glu,
                  dALi.glu = dALi.glu, 
                  dABi.glu = dABi.glu, 
                  dAKi.glu = dAKi.glu, 
                  dAUr.glu = dAUr.glu,
                  dAAd.glu = dAAd.glu, 
                  dASl.glu = dASl.glu, 
                  dARa.glu = dARa.glu, 
                  dABl.glu = dABl.glu),
                list(c(CInL = CInL, 
                       CIn = CIn, 
                       CLi = CLi ,
                       CKi = CKi, 
                       CAd = CAd, 
                       CSl = CSl, 
                       CRa = CRa, 
                       CBl = CBl, 
                       CIn.glu = CIn.glu,
                       CInL.glu = CInL.glu, 
                       CLi.glu = CLi.glu, 
                       CBi.glu = CBi.glu, 
                       CKi.glu = CKi.glu, 
                       CAd.glu = CAd.glu, 
                       CSl.glu = CSl.glu, 
                       CRa.glu = CRa.glu, 
                       CBl.glu = CBl.glu 
                ))
    ))
  })
}



# -----------------------------------------------------------------------------
# RUN THE MODEL ####
# -----------------------------------------------------------------------------

# Initialize an empty list to store results
results_list <- list()

# Loop over each ka value
for (ka_val in ka_values) {
  
  # Update the ka parameter in the PARAMETERS list
  PARAMETERS["ka"] <- ka_val



PBK <- ode(y = INITIALstate,
           times = Tout,  #used to be times
           func = PBK_ZENandGLUC,
           parms = PARAMETERS,
           atol = 1e-50,
           #events = list(data = events), #for repeated exposure
           method = "lsoda")
head(PBK)


# Convert output to data frame
PBK_df <- as.data.frame(PBK)

# Add ka value to the data frame for reference
PBK_df$ka <- ka_val

# Store the result in the list
results_list[[as.character(ka_val)]] <- PBK_df
}

# Combine all results into a single data frame
final_results <- do.call(rbind, results_list)

# Save the combined results to a CSV file
write.csv(final_results, file = here(TodaysOuput, "OUTPUT_ka_range.csv"), row.names=FALSE)


# -----------------------------------------------------------------------------
# RESULTS ####
# -----------------------------------------------------------------------------
# plot(PBK)
# ggsave(here(PlotFolderIs, "Overall.png"), dpi = 300)

# Output in amounts (umol)
OUTPUT.A <- as.data.frame(PBK) %>% select("time", "AInL" : "ABl.glu")
write.csv(OUTPUT.A, file = here(TodaysOuput, "OUTPUT.A.csv"), row.names=FALSE)


## Mass Balance / Error ####
## ------------------------------------------------------

CALCULATED <- c( #!Change accordingly depending on the compartments in the function#
  "AInL", "AIn", "AFe",
  "ALi", "AKi", "AUr", "AAd", "ASl", "ARa", "ABl",
  "AIn.glu", "AInL.glu", "AFe.glu",
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
#OUTPUT.C <- as.data.frame(PBK) %>% select(!c("AInL" : "ABl.glu"))
OUTPUT.C <- as.data.frame(final_results) %>% select(!c("AInL" : "ABl.glu"))

write.csv(OUTPUT.C, file = here(TodaysOuput, "OUTPUT.C.csv"), row.names=FALSE)


# To create df of simulated data for plotting
CONCENTRATIONS <- OUTPUT.C %>% #!change accordingly depending on the compartments in the model
  rename(#!change accordingly depending on the compartments in the model
    "Time" = "time",
    "Intestine_Lumen" = "CInL",
    "Intestine" = "CIn",
    "Liver" = "CLi",
    "Kidney" = "CKi",
    "Adipose" = "CAd",
    "Slowly_perfused" = "CSl",
    "Rapidly_perfused" = "CRa",
    "Blood" = "CBl",
    "Intestine_glucuronide" = "CIn.glu",
    "Intestine_Lumen_glucuronide" = "CInL.glu",
    "Liver_glucuronide" = "CLi.glu",
    "Bile_glucuronide" = "CBi.glu",
    "Kidney_glucuronide" = "CKi.glu",
    "Adipose_glucuronide" = "CAd.glu",
    "Slowly_perfused_glucuronide" = "CSl.glu",
    "Rapidly_perfused_glucuronide" = "CRa.glu",
    "Blood_glucuronide" = "CBl.glu"
  ) %>% pivot_longer(names_to = "Tissue",
                     values_to = "Concentration",
                     Intestine_Lumen:Blood_glucuronide)

# Simulated data for Zearalenone
CONCENTRATIONS.ZEN <- CONCENTRATIONS %>%
  filter(!str_ends(Tissue, "_glucuronide")) %>%
  mutate(Concentration = Concentration * INPUT$MW) %>%  #converting umol/L to ng/ml
  pivot_wider(names_from = Tissue, values_from = Concentration) %>%
  mutate(Intestine = Intestine_Lumen + Intestine ) %>%
  select(!"Intestine_Lumen") %>%
  pivot_longer(names_to = "Tissue", values_to = "Concentration", Intestine:Blood)

# Simulated data for Zearalenone glucuronide
CONCENTRATIONS.ZENGluc <- CONCENTRATIONS %>%
  filter(str_ends(Tissue, "_glucuronide")) %>%
  mutate(Concentration = Concentration * INPUT$MW.glu) %>% #converting umol/L to ng/ml
  mutate(Tissue = str_replace(Tissue, "_glucuronide$", "")) %>%
  pivot_wider(names_from = Tissue, values_from = Concentration) %>%
  mutate(Intestine = Intestine_Lumen + Intestine) %>%
  select(!"Intestine_Lumen") %>%
  pivot_longer(names_to = "Tissue", values_to = "Concentration", Intestine:Blood)



# To plot Zearalenone blood concentration
CONCENTRATIONS.ZEN %>%
  filter(Tissue == "Blood") %>%
  ggplot(aes(Time, Concentration, colour = ka)) +
  facet_wrap(~ka) +
  geom_path(linewidth = 0.5) +
  labs(title = "Zearalenone blood concentration",
       x = "Time (h)",
       y = "Concentration (ng/ml)") +
  CP_theme()
ggsave(here(PlotFolderIs, "ZEN.blood.simulateddependingonkTIn.png"), dpi = 300)


# To plot Zearalenone tissue concentrations
CONCENTRATIONS.ZEN %>%
  filter(Tissue != "Intestine") %>%
  ggplot(aes(Time, Concentration, colour = Tissue)) +
  facet_wrap(~ka) +
  geom_path(linewidth = 0.5) +
  labs(title = "Zearalenone tissue concentrations",
       x = "Time (h)",
       y = "Concentration (ng/ml)") +
  CP_theme()
ggsave(here(PlotFolderIs, "ZEN.simulatedtissue.png"), dpi = 300)



# To plot Zearalenone and Zearalenone glucuronide blood concentrations
BLOOD_CONCENTRATIONS.ZEN <- CONCENTRATIONS.ZEN %>% filter(Tissue == "Blood") %>% mutate(Compound = "Zearalenone")
BLOOD_CONCENTRATIONS.ZENGluc <- CONCENTRATIONS.ZENGluc %>% filter(Tissue == "Blood") %>% mutate(Compound = "Z.Glucuronide")
BLOOD_CONCENTRATIONS <- rbind(BLOOD_CONCENTRATIONS.ZEN, BLOOD_CONCENTRATIONS.ZENGluc)

BLOOD_CONCENTRATIONS %>%
  ggplot(aes(Time, Concentration)) +
  geom_path(linewidth = 0.5) +
  facet_wrap(~Compound) +
  labs(title = "Zearalenone blood concentrations",
       x = "Time (h)",
       y = "Concentration (ng/ml)") +
  CP_theme()
ggsave(here(PlotFolderIs, "ZEN.blood.simulated.png"), dpi = 300)


# To plot Zearalenone glucuronide tissue concentrations
CONCENTRATIONS.ZENGluc %>%
  filter(Tissue != "Bile") %>%
  ggplot(aes(Time, Concentration, colour = Tissue)) +
  geom_path(linewidth = 0.5) +
  facet_wrap(~ka) +
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
  facet_wrap(~ka) +
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
  facet_wrap(~ka) +
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

# Blood
InvivoBlood <- POData.ZENGluc %>% filter(Tissue == "Plasma") %>%
  #mutate(Concentration = Concentration*ZEN.dat$BP) %>% # assuming it's one as I don't have an experimental value
  mutate(Tissue = recode(Tissue, "Plasma" = "Blood")) #renaming plasma to blood


BLOOD_CONCENTRATIONS.ZENGluc %>% ggplot() +
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
