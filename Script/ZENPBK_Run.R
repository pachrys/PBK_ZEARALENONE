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
# Script for calculated input
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
  In = 0,
  #D = 0,
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

OUTPUT <- as.data.frame(PBK) #(umol)
write.csv(OUTPUT, file = here(TodaysOuput, "OUTPUT.amounts21.csv"), row.names=FALSE)


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
Total <- 7.814339e-01


Massbalance <- Total - Calculated + 1
ErrPOl <- (Total - Calculated) / (Total + 1*10^(-30)) * 100                     #mass error (% of mass lost)

MassBalance_Error <- cbind("time" = OUTPUT[, "time"],                       #data frame
                           Massbalance, ErrPOl) %>%
  as.data.frame()

# Mass Balance and Error plot
MassBalance_Error %>%
  filter(time >= 0.01) %>% #as PO administration time is 0,01h
  ggplot(aes(time, Massbalance)) +
  geom_line(color = "grey") +
  geom_line(aes(time, ErrPOl), color = "grey6") +
  labs(x = "Time (h)", y = "Value",
       caption = "Massbalance in grey, Error in black") +
  CP_theme()
ggsave(here(TodaysOuput, "MassbalanceError.png"), dpi = 300)             #saving the plot




#**STEP 5**#
#****************************************************************************#

#------------------------Plotting Simulated Output---------------------------#
# 
# 
# # From amounts to concentrations
# OUTPUT %>% 
#   select(c("AInL":"ABl.glu")) %>% 
#   mutate(
#     AInL =   
#   )
# 
# # To create df of simulated data for plotting
# PlotConc <- OUTPUT %>% #!change accordingly depending on the compartments in the model
#   rename("Time" = "time",
#          "Intestine" = "CIn",
#          "Liver" = "CLi",
#          "Liver_glu" = "CLi.glu",
#          "Kidney" = "CKi",
#          "AdiPOse" = "CAd",
#          "POorlyPerfused" = "CSl",
#          "RapidlyPerfused" = "CRa",
#          "Blood" = "CBl",
#          "LiverVenous" = "CvLi",
#          "KidneyVenous" = "CvKi",
#          "AdiPOseVenous" = "CvAd",
#          "POorlyPerfusedVenous" = "CvSl",
#          "RapidlyPerfusedVenous" = "CvRa") %>%
#   select(Time, Intestine:RapidlyPerfusedVenous) %>% #!change accordingly depending on the compartments in the model
#   pivot_longer(names_to = "Tissue",
#                values_to = "Concentration",
#                Intestine:RapidlyPerfusedVenous,
#   ) %>%
#   mutate(Concentration = Concentration*INPUT$MW)  #converting umol/L to ng/ml
# 
# # To plot only simulated data
# PlotConc %>%
#   filter(Tissue %in% c("Liver", "Blood", "Kidney", "Intestine",
#                        "AdiPOse", "POorlyPerfused", "RapidlyPerfused")) %>%
#   ggplot(aes(Time, Concentration, color = Tissue)) +
#   geom_line() +
#   labs(title = "Simulated zearalenone tissue concentrations",
#        x = "Time (h)",
#        y = "Concentration (ng/ml)",
#        caption = paste("Administration type:", AdminType, "Dose:", PO_DOSE.mgkg, "mg/kg bw")) +
#   CP_theme() +
#   scale_colour_manual(values=Chrysa_colors)
# ggsave(here(PlotFolderIs, "Simulated.png"), dpi = 300)
# 
# # To plot blood concentrations
# PlotConc %>%
#   filter(Tissue == "Blood") %>%
#   ggplot(aes(Time, Concentration)) +
#   geom_path(linewidth = 0.5) +
#   labs(title = "Zearalenone blood concentrations",
#        x = "Time (h)",
#        y = "Concentration (ng/ml)",
#        caption = paste(AdminType)) +
#   CP_theme() +
#   scale_colour_manual(values=Chrysa_colors)
# ggsave(here(PlotFolderIs, "Simulated_blood.png"), dpi = 300)
# 
# 
# # To plot glucuronide
# PlotConc %>%
#   filter(Tissue %in% c("Intestine", "Liver", "Liver_glu")) %>% #"Intestine_glu", )) %>%
#   ggplot(aes(Time, Concentration, color = Tissue)) +
#   geom_path(linewidth = 0.5) +
#   labs(title = "Zearalenone and zearalenone-glucuronide intestine",
#        x = "Time (h)",
#        y = "Concentration (ng/ml)",
#        caption = paste(AdminType)) +
#   CP_theme() +
#   scale_colour_manual(values=Chrysa_colors)
# ggsave(here(PlotFolderIs, "Simulated_GlucTissueConc.png"), dpi = 300)
# 
# 
# 
# #---End-----------------------------------------------------------------------#
# #*****************************************************************************#
# 
# 
