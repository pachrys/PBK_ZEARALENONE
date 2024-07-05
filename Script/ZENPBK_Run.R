# Libraries-------------------------------------------------------------------#

library(deSolve)
library(tidyverse)
library(ggplot2)
library(here)
library(showtext)
font_add(family = "Garamond", regular = "GARA.TTF")
showtext_auto()

# Working directory-----------------------------------------------------------#
WhatAmITesting <- "RunWithDeSolve"
TodaysOuput <- here("Output/Data", Sys.Date(), WhatAmITesting)      #folder for saving results
dir.create(TodaysOuput, recursive = TRUE)


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

PO_TIME = 0.01 #h time it takes to administer the dose to the intestinal lumen; I checked that the amount at t=0.01 in InL is 0.781408707582861

# dosing = data.frame(var = c("AInL"),
#                    time = c(seq(0, 24, by = 24)), 
#                    value = c(PO_DOSE),
#                    method = c("add"))

# Model parameters 
PARAMETERS=c(PO_DOSE = PO_DOSE, 
             PO_TIME = PO_TIME,
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
             QHa = INPUT$QLi - INPUT$QPo, #hepatic artery (Qliver is sum of both Qhepatic artery and QPOrtal vein)
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
             VmaxExc.In.glu = INPUT$VmaxExc.In.glu, #active transPOrt from intestinal tissue to lumen
             KmGluc.Exc = INPUT$KmGluc.Exc, #affinity constant zearalenone glucuronide to efflux transPOrter
             VmaxAbs.In.glu = INPUT$VmaxAbs.In.glu, #active transPOrt from lumen to tissue
             KmGluc.Abs = INPUT$KmGluc.Abs #affinity constant zearalenone glucuronide to influx transPOrer
             ) 


#Initial state
INITIALstate=c(
  AInL = 0, 
  AIn = 0, 
  AFe = 0, 
  ALi = 0,
  #ABi = 0,
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
  ABl.glu = 0)


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

OUTPUT.A <- as.data.frame(PBK) #(umol)
write.csv(OUTPUT.A, file = here(TodaysOuput, "OUTPUT.A2.csv"), row.names=FALSE)


# Mass Balance
CALCULATED <- c( #!Change accordingly depending on the compartments in the function#
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
CALCULATED <- rowSums(OUTPUT.A[ ,CALCULATED])
TOTAL <- 0.7814087


MASSBALANCE <- TOTAL - CALCULATED + 1
ERROR <- (TOTAL - CALCULATED) / (TOTAL + 1*10^(-30)) * 100                     #mass error (% of mass lost)

MASSBALANCE_Error <- cbind("time" = OUTPUT.A[, "time"],                       #data frame
                           MASSBALANCE, ERROR) %>%
  as.data.frame()

# Mass Balance and Error plot
MASSBALANCE_Error %>%
  filter(time >= 0.01) %>% #as PO administration time is 0,01h
  ggplot(aes(time, MASSBALANCE)) +
  geom_line(color = "grey") +
  geom_line(aes(time, ERROR), color = "grey6") +
  labs(x = "Time (h)", y = "Value",
       caption = "MASSBALANCE in grey, Error in black") +
  CP_theme()
ggsave(here(TodaysOuput, "MB_ERROR.png"), dpi = 300)             #saving the plot




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
  select(!c("In", "AFe", "AUr", "AFe.glu", "AUr.glu"))
  
names(OUTPUT.C) <- sub("^A", "C", names(OUTPUT.C))

# To create df of simulated data for plotting
CONCENTRATIONS <- OUTPUT.C %>% #!change accordingly depending on the compartments in the model
  rename(#!change accordingly depending on the compartments in the model
    "Time" = "time",
    "Intestinal_lumen" = "CInL",
    "Intestinal_tissue" = "CIn",
    "Liver" = "CLi",
    "Kidney" = "CKi",
    "Adipose" = "CAd",
    "Slowly_perfused" = "CSl",
    "Rapidly_perfused" = "CRa",
    "Blood" = "CBl",
    "Intestinal_tissue_glucuronide" = "CIn.glu",
    "Intestinal_lumen_glucuronide" = "CInL.glu",
    "Liver_glucuronide" =  "CLi.glu",
    "Bile_glucuronide" = "CBi.glu",
    "Kidney_glucuronide" = "CKi.glu",
    "Adipose_glucuronide" = "CAd.glu",
    "Slowly_perfused_glucuronide" = "CSl.glu",
    "Rapidly_perfused_glucuronide" = "CRa.glu",
    "Blood_glucuronide" = "CBl.glu" 
  ) %>% 
  pivot_longer(names_to = "Tissue",
               values_to = "Concentration",
               Intestinal_lumen:Blood_glucuronide,
  ) %>%
  mutate(Concentration = Concentration*INPUT$MW)  #converting umol/L to ng/ml

# To plot only simulated data
CONCENTRATIONS %>%
  filter(Tissue %in% c("Liver", "Blood", "Kidney", "Intestine",
                       "Adipose", "PIorlyPerfused", "RapidlyPerfused")) %>%
  ggplot(aes(Time, Concentration, color = Tissue)) +
  geom_line() +
  labs(title = "Simulated zearalenone tissue concentrations",
       x = "Time (h)",
       y = "Concentration (ng/ml)",
       caption = paste("Administration type:", "Dose:1 mg/kg bw")) +
  CP_theme() #+
  #scale_colour_manual(values=Chrysa_colors)
ggsave(here(PlotFolderIs, "Simulated.png"), dpi = 300)

# To plot blood concentrations
CONCENTRATIONS %>%
  filter(Tissue == "Blood") %>%
  ggplot(aes(Time, Concentration)) +
  geom_path(linewidth = 0.5) +
  labs(title = "Zearalenone blood concentrations",
       x = "Time (h)",
       y = "Concentration (ng/ml)") +
  CP_theme() 
ggsave(here(PlotFolderIs, "Simulated_blood.png"), dpi = 300)


# To plot glucuronide
CONCENTRATIONS %>%
  filter(Tissue %in% c("Intestine", "Liver", "Liver_glu")) %>% #"Intestine_glu", )) %>%
  ggplot(aes(Time, Concentration, color = Tissue)) +
  geom_path(linewidth = 0.5) +
  labs(title = "Zearalenone and zearalenone-glucuronide intestine",
       x = "Time (h)",
       y = "Concentration (ng/ml)",
       caption = paste(AdminType)) +
  CP_theme() +
  scale_colour_manual(values=Chrysa_colors)
ggsave(here(PlotFolderIs, "Simulated_GlucTissueConc.png"), dpi = 300)



#---End-----------------------------------------------------------------------#
#*****************************************************************************#


