# Libraries-------------------------------------------------------------------#

library(deSolve)
library(tidyverse)
library(ggplot2)
library(here)
library(showtext)
font_add(family = "Garamond", regular = "GARA.TTF")
showtext_auto()

# Working directory-----------------------------------------------------------#
WhatAmITesting <- "SimpleJustLiverMetabolism"
TodaysOuput <- here("Output/Data", Sys.Date(), WhatAmITesting)
dir.create(TodaysOuput, WhatAmITesting,recursive = TRUE)
PlotFolderIs <- here("Output/Plots", Sys.Date(), WhatAmITesting)
dir.create(PlotFolderIs, WhatAmITesting,recursive = TRUE)


#*****************************************************************************#
#---Run script----------------------------------------------------------------#
#*****************************************************************************#

#----Extracting input---------------------------------------------------------#
# Script for CALCULATED input
INPUT <- read.csv(here("InputData/CALCULATED.INPUT.csv"), header=1) %>%  #The function for calculating input can be called here: source(here("Script/ZENPBK_Input"))                                         #call input script
  pivot_wider(names_from = "Parameter", 
              values_from = "Value")
ZEN.dat <- read.csv(here("InputData/Zen.input.csv"), header=1)     #from {Agahi, 2020 #456} and chemicalize, BP from Shin, Hong et al.2009


#----Set Parameters for the simulation----------------------------------------#

# Tout
STARTtime <- 0 # start of simulation (h)
STOPtime <- 24 # end of simmulation (h) (can be from 24h to years)
Dt <- 1 / 100 # steps of simulation (now it per minute) (could be 1/100, 1/200 or 1 for simulations over years)
times <- seq(STARTtime,STOPtime,by=Dt)

# Exposure
PO_CONC = 1 #mg/kg day, Videmann
PO_DOSE = PO_CONC*INPUT$BW / INPUT$MW * 1000 #umol


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
             # kTIn = INPUT$kTIn,#total intestinal transit time
             kTsi = 1 / ZEN.dat$kTsi, 
             kTco = 1 / ZEN.dat$kTco,
             kTge = INPUT$kTge, #gastric emptying rate, "Assumed to mimic bile secretion rate
             CL_Re = INPUT$CL_Re,
             # CL_Bi = INPUT$CL_Bi, #biliary clearance
             # VmaxMet.In = INPUT$VmaxMet.In, #intestinal glucuronidation
             # KmGluc.In = INPUT$KmGluc.In,  #affinity constant zearalenone to UGT
             VmaxMet.Li = INPUT$VmaxMet.Li, #liver glucuronidation
             KmGluc.Li = INPUT$KmGluc.Li
             # VmaxExc.In.glu = INPUT$VmaxExc.In.glu, #active transPOrt from intestinal tissue to lumen
             # KmGluc.Exc = INPUT$KmGluc.Exc, #affinity constant zearalenone glucuronide to efflux transPOrter
             # VmaxAbs.In.glu = INPUT$VmaxAbs.In.glu, #active transPOrt from lumen to tissue
             # KmGluc.Abs = INPUT$KmGluc.Abs, #affinity constant zearalenone glucuronide to influx transPOrer
             # kHydro = 1000 # /hhydrolysis constant; assumption as in Veronique
) 

PARAMETERS

# Initial state
INITIALstate = c(
  ASt = PO_DOSE, AIn = 0, ACo = 0, AFe = 0, ALi = 0, ALi.glu = 0, AKi = 0, AUr = 0, 
  AAd = 0, ASl = 0, ARa = 0, ABl = 0
)

INITIALstate

#----PBK function-------------------------------------------------------------#
# source(here("Script/ZENPBK_Function.R"))

PBK_ZENandGLUC <- function(t,state,parameters){
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
    
    
    #----Differential equations-----------------------------------------------#
    # Rate of change in amounts in each tissue (umol/h)
    
    # Zearalenone-------------------- ##
    dASt <- -(kTge*ASt)                                                      #(umol/h) rate of change of amount reaching the intestine (correcting for gastric emptying time)
    
    dAIn <- (kTge*ASt) - (ka*AIn) - (kTsi*AIn)                               #(umol/h) rate of change in amount of ZEN in the intestine: the dose + the amount of zen formed after microbial hydrolysis of zen-gluc - fecal excretion
    dACo <- (kTsi*AIn) - (kTco*ACo)                                          #(umol/h) rate of change in amount of ZEN in colon
    dAFe <- kTco*ACo                                                         #(umol/h) rate of change in amount of ZEN eliminated in the feces
    
    #Differential equations perfused tissues
    dALi <- (ka*AIn) + QLi*(CBl - CvLi) +
      - ((VmaxMet.Li * CLi_free) / (KmGluc.Li + CLi_free))                            #(umol/h) rate of change in amount of ZEN in the liver: what is absorbed from the intestine + what is absorbed and leaves from/to the blood - rate of metabolism (umol/h) is metabolised to glucuronides
    dALi.glu <- ((VmaxMet.Li * CLi_free) / (KmGluc.Li + CLi_free))                                                           #(umol/h) rate of change in amount in ZEN glucuronides the liver: what is formed after glucuronidation minus what is secreted to the bile 
    
    dAKi <- QKi*(CBl - CvKi) - (CL_Re*CKi_free)                               #(umol/h) rate of change in amount in the kidney
    dAUr <- CL_Re*CKi_free                                                                #(umol/h) rate of change in amount ZEN eliminated in urine
    
    dAAd <- QAd*(CBl - CvAd)                                                 #(umol/h) rate of change in amount in the adipose tissue
    dASl <- QSl*(CBl - CvSl)                                                 #(umol/h) rate of change in amount in the slowly perfused
    dARa <- QRa*(CBl - CvRa)                                                 #(umol/h) rate of change in amount in the rapidly perfused
    
    dABl <- (QLi*CvLi) + (QKi*CvKi) +
      (QAd*CvAd) + (QSl*CvSl) + (QRa*CvRa) +
      -(QLi + QKi + QAd + QSl + QRa)*(CBl)                                       #(umol/h) rate of change in amount in the blood
    
    
    #----Output---------------------------------------------------------------#
    
    return(list(c(dASt, dAIn, dACo, dAFe, dALi, dALi.glu, dAKi, dAUr, 
                  dAAd, dASl, dARa, dABl),   
    CIn, #22, 
    CLi, #23
    CKi, #24
    CAd, #25
    CSl, #26
    CRa, #27
    CBl #28
    ))
  })
}

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
  "ASt", "AIn", "ACo", "AFe", "ALi", "ALi.glu", "AKi", "AUr.glu", 
  "AAd", "ASl", "ARa", "ABl"
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
    AIn = AIn/INPUT$VIn, 
    ALi = ALi/INPUT$VLi,
    AKi = AKi/INPUT$VKi,
    AAd = AAd/INPUT$VAd, 
    ASl = ASl/INPUT$VSl,
    ARa = ARa/INPUT$VRa,
    ABl = ABl/INPUT$VBl
  ) %>% 
  select(!c("AFe", "AUr"))
  
names(OUTPUT.C) <- sub("^A", "C", names(OUTPUT.C))

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
    "Blood" = "CBl"
  ) %>% 
  pivot_longer(names_to = "Tissue",
               values_to = "Concentration",
               Intestine:Blood,
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

#---End-----------------------------------------------------------------------#
#*****************************************************************************#


