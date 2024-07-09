# Libraries-------------------------------------------------------------------#

library(deSolve)
library(tidyverse)
library(ggplot2)
library(here)
library(showtext)
font_add(family = "Garamond", regular = "GARA.TTF")
showtext_auto()

# Working directory-----------------------------------------------------------#
WhatAmITesting <- "StepByStep12"
TodaysOuput <- here("Output/Data", Sys.Date(), WhatAmITesting)
dir.create(TodaysOuput, WhatAmITesting,recursive = TRUE)
PlotFolderIs <- here("Output/Plots", Sys.Date(), WhatAmITesting)
dir.create(PlotFolderIs, WhatAmITesting,recursive = TRUE)


#*****************************************************************************#
#---Run script----------------------------------------------------------------#
#*****************************************************************************#

#----Bring PBK function-------------------------------------------------------#
source(here("Script/ZENPBK_Function.R"))


#----Extracting input---------------------------------------------------------#
# Script for CALCULATED input
INPUT <- read.csv(here("InputData/CALCULATED.INPUT.csv"), header=1) %>%  #The function for calculating input can be called here: source(here("Script/ZENPBK_Input"))                                         #call input script
  pivot_wider(names_from = "Parameter", 
              values_from = "Value")

#----Set Parameters for the simulation----------------------------------------#

# Tout
STARTtime <- 0 # start of simulation (h)
STOPtime <- 24 # end of simmulation (h) (can be from 24h to years)
Dt <- 1/100 # steps of simulation (now it per minute) (could be 1/100, 1/200 or 1 for simulations over years)
times <- seq(STARTtime,STOPtime,by=Dt)

# ExPOsure
PO_CONC = 1 #mg/kg day, Videmann
PO_DOSE = PO_CONC*INPUT$BW/INPUT$MW*1000 #umol


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
             PInb.glu = INPUT$PIn.glu, 
             PLib.glu = INPUT$PLi.glu,
             PKib.glu = INPUT$PKi.glu, 
             PAdb.glu = INPUT$PAd.glu,
             PSlb.glu = INPUT$PSl.glu,
             PRab.glu = INPUT$PRa.glu,
             fub.glu = INPUT$fub.glu,
             QIn = INPUT$QIn,
             QPo = INPUT$QPo,
             QLi = INPUT$QLi,
             QHa = INPUT$QLi - INPUT$QPo, #hepatic artery (Qliver is sum of both Qhepatic artery and QPOrtal vein)
             QKi = INPUT$QKi,
             QAd = INPUT$QAd, 
             QSl = INPUT$QSl,  
             QRa = INPUT$QRa,
             ka = INPUT$ka.is,
             ka.is = INPUT$ka.is,
             ka.sl = INPUT$ka.sl, 
             ka.fa = INPUT$ka.fa, 
             ka.fit = INPUT$ka.fit, 
             kex = INPUT$kex, 
             kfil = INPUT$kfil, 
             kbil = INPUT$kbil, 
             ka.sl.glu = INPUT$ka.sl.glu, 
             ka.fa.glu = INPUT$ka.fa.glu, 
             ka.fit.glu = INPUT$ka.fit.glu, 
             kex.glu = INPUT$kex.glu, 
             kfil.glu = INPUT$kfil.glu, 
             kbil.glu = INPUT$kbil.glu, 
             VmaxAbs = INPUT$VmaxAbs, 
             KmAbs = INPUT$KmAbs, 
             VmaxSec = INPUT$VmaxSec, 
             KmSec = INPUT$KmSec, 
             kTIn = INPUT$kTIn,
             kTge = INPUT$kTge, #gastric emptying rate, assumed to mimic bile secretion rate
             CL_Re = INPUT$CL_Re,
             CL_Bi = INPUT$CL_Bi, #biliary clearance
             VmaxMet.In = INPUT$VmaxMet.In, #intestinal glucuronidation
             KmGluc.In = INPUT$KmGluc.In*1000,  #affinity constant zearalenone to UGT
             VmaxMet.Li = INPUT$VmaxMet.Li, #liver glucuronidation
             KmGluc.Li = INPUT$KmGluc.Li,
             VmaxExc.In.glu = INPUT$VmaxExc.In.glu, #active transPOrt from intestinal tissue to lumen
             KmGluc.Exc = INPUT$KmGluc.Exc, #affinity constant zearalenone glucuronide to efflux transPOrter
             VmaxAbs.In.glu = INPUT$VmaxAbs.In.glu, #active transPOrt from lumen to tissue
             KmGluc.Abs = INPUT$KmGluc.Abs, #affinity constant zearalenone glucuronide to influx transPOrer
             kHydro = 1000 # /hhydrolysis constant; assumption as in Veronique
             ) 

#Initial state
INITIALstate=c(
  AInL = PO_DOSE, 
  # AInL.f = 0,
  # AInL.s = 0,
  AIn = 0, 
  AFe = 0,
  # ALi = 0,
  # #ABi = 0,
  # AKi = 0, 
  # AUr = 0, 
  # AAd = 0, 
  # ASl = 0, 
  # ARa = 0, 
  # ABl = 0,
  AIn.glu = 0,
  AInL.glu = 0
  # AFe.glu = 0, 
  # ABl.glu = 0, 
  # ALi.glu = 0,
  # ABi.glu = 0,
  # AKi.glu = 0, 
  # AUr.glu = 0, 
  # AAd.glu = 0, 
  # ASl.glu = 0, 
  # ARa.glu = 0 
  )


#----Simulation---------------------------------------------------------------#
PBK <- ode(y = INITIALstate, 
           times = times,  
           func = PBK_ZENandGLUC,
           parms = PARAMETERS, 
           atol = 1e-50,
           #events = list(data = dosing),
           method = "lsoda")
head(PBK)

#----Results------------------------------------------------------------------#
plot(PBK)
ggsave(here(PlotFolderIs, "Overall.png"), dpi = 300)


OUTPUT.A <- as.data.frame(PBK) #(umol)
write.csv(OUTPUT.A, file = here(TodaysOuput, "OUTPUT.A.csv"), row.names=FALSE)


# Mass Balance
CALCULATED <- c( #!Change accordingly depending on the compartments in the function#
  "AInL",
  # "AInL.f",
  # "AInL.s",
  "AIn",
  "AFe",
  # "ALi",
  # "AKi",
  # "AUr",
  # "AAd",
  # "ASl",
  # "ARa",
  # "ABl",
  "AIn.glu",
  "AInL.glu"
  # "AFe.glu",
  # "ABl.glu", 
  # "ALi.glu",
  # "ABi.glu",
  # "AKi.glu",
  # "AUr.glu",
  # "AAd.glu",
  # "ASl.glu",
  # "ARa.glu"
  )

# Plotting function
source(here("Script/Functions/PlottingFunction.R"))                               #call plotting function


#-----Mass Balance Calculations----------------------------------------------#
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




#**STEP 5**#
#****************************************************************************#

#------------------------Plotting Simulated Output---------------------------#


# From amounts to concentrations
OUTPUT.C <- OUTPUT.A %>% 
  mutate(
    AInL = AInL/INPUT$VInL, 
    AIn = AIn/INPUT$VIn, 
    ALi = ALi/INPUT$VLi,
    AKi = AKi/INPUT$VKi,
    AAd = AAd/INPUT$VAd, 
    ASl = ASl/INPUT$VSl,
    ARa = ARa/INPUT$VRa,
    ABl = ABl/INPUT$VBl,
    AIn.glu = AIn.glu/INPUT$VIn,
    AInL.glu = AInL.glu/INPUT$VInL,
    ALi.glu = ALi.glu/INPUT$VLi,
    ABi.glu = ABi.glu/INPUT$VBi,
    AKi.glu = AKi.glu/INPUT$VKi,
    AAd.glu = AAd.glu/INPUT$VAd,
    ASl.glu = ASl.glu/INPUT$VSl,
    ARa.glu = ARa.glu/INPUT$VRa,
    ABl.glu = ABl.glu/INPUT$VBl
  ) %>% 
  select(!c("AFe", "AUr", "AFe.glu", "AUr.glu"))
  
names(OUTPUT.C) <- sub("^A", "C", names(OUTPUT.C))

write.csv(OUTPUT.C, file = here(TodaysOuput, "OUTPUT.C.csv"), row.names=FALSE)


# To create df of simulated data for plotting
CONCENTRATIONS <- OUTPUT.C %>% #!change accordingly depending on the compartments in the model
  rename(#!change accordingly depending on the compartments in the model
    "Time" = "time",
    "Intestinal_lumen" = "CInL",
    "Intestine" = "CIn",
    "Liver" = "CLi",
    "Kidney" = "CKi",
    "Adipose" = "CAd",
    "Slowly_perfused" = "CSl",
    "Rapidly_perfused" = "CRa",
    "Blood" = "CBl",
    "Intestine_glucuronide" = "CIn.glu",
    "Intestinal_lumen_glucuronide" = "CInL.glu",
    "Blood_glucuronide" = "CBl.glu",
    "Liver_glucuronide" =  "CLi.glu",
    "Bile_glucuronide" = "CBi.glu",
    "Kidney_glucuronide" = "CKi.glu",
    "Adipose_glucuronide" = "CAd.glu",
    "Slowly_perfused_glucuronide" = "CSl.glu",
    "Rapidly_perfused_glucuronide" = "CRa.glu"
  ) %>% 
  pivot_longer(names_to = "Tissue",
               values_to = "Concentration",
               Intestinal_lumen:Rapidly_perfused_glucuronide,
  ) %>%
  mutate(Concentration = Concentration*INPUT$MW)  #converting umol/L to ng/ml



# To plot Zearalenone blood concentration
CONCENTRATIONS %>%
  filter(str_starts(Tissue, "Blood"),
         !str_ends(Tissue, "_glucuronide")) %>%
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

#---End-----------------------------------------------------------------------#
#*****************************************************************************#


