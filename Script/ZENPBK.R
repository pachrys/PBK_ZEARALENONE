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

source(here("Script/Functions/PlottingFunction.R"))                               #call plotting function


# -----------------------------------------------------------------------------
# PBK MODEL FOR ZEARALENONE ####
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
STOPtime <- 24 # end of simmulation (h) (can be from 24h to years)
Dt <- 1 / 100 # steps of simulation (now it per minute) (could be 1/100, 1/200 or 1 for simulations over years)
times <- seq(STARTtime,STOPtime,by=Dt)

# Exposure
PO_CONC = 1 #mg/kg day, Videmann
PO_DOSE = PO_CONC*INPUT$BW / INPUT$MW * 1000 #umol


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


## Initials ####
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
    ),  
    CInL, #"22" 
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
    ))
  })
}



# -----------------------------------------------------------------------------
# RUN THE MODEL ####
# -----------------------------------------------------------------------------

PBK <- ode(y = INITIALstate, 
           times = times,  
           func = PBK_ZENandGLUC,
           parms = PARAMETERS, 
           atol = 1e-50,
           #events = list(data = dosing),
           method = "lsoda")
head(PBK)



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
OUTPUT.C <- as.data.frame(PBK) %>% select(!c("AInL" : "ABl.glu")) %>% 
  rename(
    "Time" = "time",
    "Intestinal_Lumen" = "22", 
    "Intestine" = "23", 
    "Liver" = "24", 
    "Kidney" = "25", 
    "Adipose" = "26", 
    "Slowly_perfused" = "27", 
    "Rapidly_perfused" = "28", 
    "Blood" = "29", 
    "Intestine_glucuronide" = "30", 
    "Intestinal_lumen_glucuronide" = "31", 
    "Liver_glucuronide" = "32",
    "Bile_glucuronide" = "33",
    "Kidney_glucuronide" = "34", 
    "Adipose_glucuronide" = "35", 
    "Slowly_perfused_glucuronide" = "36", 
    "Rapidly_perfused_glucuronide" = "37", 
    "Blood_glucuronide" = "38" 
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
               Intestinal_Lumen:Blood_glucuronide,
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
DataInvivo.ZEN <- read.csv(here("InputData/InVivoData.csv"), header=1) %>% 
  filter(Author == "Videmann", #I chose Videmann as it had the best recoveries, also Videmann separates Zearalenone from it's metabolites
         Gender != "Pregnant.F", #pregnancy not included yet
         TestCompound == "ZEN", #other options are ZEN_glu and a_ZEL
         AdminType == "PO", #"IV"), #Videmann doesn't have IV, but subcutaneous!!
         Tissue %in% c("Plasma", "Intestine", "Liver", "Adipose", "Kidney", "Urine", "Colon"), 
         SampleCompound == "ZEN") %>%  #only for ZEN 
  select(c("Gender", "SamplingTime", "Tissue", "Concentration"))%>% 
  rename(Time = SamplingTime) %>% #(h)
  mutate(Concentration = Concentration*ZEN.dat$BP) %>% #(ng/ml) blood concentration: changing from plasma concentrations to blood concentrations by multiplying with BP
  mutate(Concentration = Concentration/0.76) %>%  #correcting for the recovery (recovery was 72-80% => averaged that to 76% and then corrected the concentrations for that)
  mutate(Tissue = recode(Tissue, "Plasma" = "Blood")) #renaming plasma to blood 

# To create a file for either PO or IV admin
POData <- DataInvivo.ZEN 



## Plotting Simulated Vs Experimental data ####
## ------------------------------------------------------

# Simulated data for Zearalenone 
CONCENTRATIONS.ZEN <- CONCENTRATIONS %>% 
  filter(!str_ends(Tissue, "_glucuronide"))


# Creating dataframes with the same tissues for both in vivo and simulated
Invivo <- POData %>% 
  filter(Tissue %in% c("Blood",
                       "Intestine", "Liver", "Kidney")) 

Simulated <- CONCENTRATIONS.ZEN %>%                                                         
  filter(Tissue %in% c("Blood",                                                  
                       "Intestine", "Liver", "Kidney")) 

# Simulated VS Experimental all tissue AUC 
Simulated %>% # choose compartment accordingly
  ggplot() +
  geom_point(data = Invivo, aes(Time, Concentration, color = Gender)) + 
  geom_line(data = Invivo, aes(Time, Concentration, color = Gender), 
            linewidth = 0.1) + 
  facet_wrap(~Tissue) +
  geom_path(aes(Time, Concentration, color = "black")) +
  facet_wrap(~Tissue) +
  #scale_colour_manual(values = Chrysa_colors) +
  labs(title = "Zearalenone tissue concentrations",
       x = "Time (h)",
       y = "Concentration (ng/ml)",
       color = "Gender") +
  CP_theme()
ggsave(here(PlotFolderIs, "AllTissueSimulatedInvivo.png"), dpi = 300)

# Separate organs
# Blood
InvivoBlood <- POData %>% filter(Tissue =="Blood") 

SimulatedBlood <- CONCENTRATIONS.ZEN %>% filter(Tissue =="Blood") 

# Simulated VS Experimental all tissue AUC 
SimulatedBlood %>% # choose compartment accordingly
  ggplot() +
  geom_point(data = InvivoBlood, aes(Time, Concentration, color = Gender)) + 
  geom_line(data = InvivoBlood, aes(Time, Concentration, color = Gender), 
            linewidth = 0.1) + 
  geom_path(aes(Time, Concentration)) +
  #scale_colour_manual(values = Chrysa_colors) +
  labs(title = "Zearalenone blood concentrations",
       x = "Time (h)",
       y = "Concentration (ng/ml)",
       color = "Gender") +
  CP_theme()
ggsave(here(PlotFolderIs, "BloodSimulatedInvivo.png"), dpi = 300)

# Intestine
InvivoIntestine <- POData %>% filter(Tissue =="Intestine") 

SimulatedIntestine <- CONCENTRATIONS.ZEN %>% filter(Tissue =="Intestine") 

# Simulated VS Experimental all tissue AUC 
SimulatedIntestine %>% # choose compartment accordingly
  ggplot() +
  geom_point(data = InvivoIntestine, aes(Time, Concentration, color = Gender)) + 
  geom_line(data = InvivoIntestine, aes(Time, Concentration, color = Gender), 
            linewidth = 0.1) + 
  geom_path(aes(Time, Concentration, color = "black")) +
  #scale_colour_manual(values = Chrysa_colors) +
  labs(title = "Zearalenone Intestine concentrations",
       x = "Time (h)",
       y = "Concentration (ng/ml)",
       color = "Gender") +
  CP_theme()
ggsave(here(PlotFolderIs, "IntestineSimulatedInvivo.png"), dpi = 300)

# Liver
InvivoLiver <- POData %>% filter(Tissue =="Liver") 

SimulatedLiver <- CONCENTRATIONS.ZEN %>% filter(Tissue =="Liver") 

# Simulated VS Experimental all tissue AUC 
SimulatedLiver %>% # choose compartment accordingly
  ggplot() +
  geom_point(data = InvivoLiver, aes(Time, Concentration, color = Gender)) + 
  geom_line(data = InvivoLiver, aes(Time, Concentration, color = Gender), 
            linewidth = 0.1) + 
  geom_path(aes(Time, Concentration, color = "black")) +
  #scale_colour_manual(values = Chrysa_colors) +
  labs(title = "Zearalenone Liver concentrations",
       x = "Time (h)",
       y = "Concentration (ng/ml)",
       color = "Gender") +
  CP_theme()
ggsave(here(PlotFolderIs, "LiverSimulatedInvivo.png"), dpi = 300)

# Kidney
InvivoKidney <- POData %>% filter(Tissue =="Kidney") 

SimulatedKidney <- CONCENTRATIONS.ZEN %>% filter(Tissue =="Kidney") 

# Simulated VS Experimental all tissue AUC 
SimulatedKidney %>% # choose compartment accordingly
  ggplot() +
  geom_point(data = InvivoKidney, aes(Time, Concentration, color = Gender)) + 
  geom_line(data = InvivoKidney, aes(Time, Concentration, color = Gender), 
            linewidth = 0.1) + 
  geom_path(aes(Time, Concentration, color = "black")) +
  #scale_colour_manual(values = Chrysa_colors) +
  labs(title = "Zearalenone Kidney concentrations",
       x = "Time (h)",
       y = "Concentration (ng/ml)",
       color = "Gender") +
  CP_theme()
ggsave(here(PlotFolderIs, "KidneySimulatedInvivo.png"), dpi = 300)
