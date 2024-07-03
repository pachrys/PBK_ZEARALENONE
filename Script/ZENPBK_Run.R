# Libraries-------------------------------------------------------------------#

# library(deSolve)
# library(tidyverse)
# library(ggplot2)
# library(here)
# library(showtext)
# font_add(family = "Garamond", regular = "GARA.TTF")
# showtext_auto() 

# Working directory-----------------------------------------------------------#
HOME <- "C:/Users/pacho003/OneDrive - Wageningen University & Research/C Channel/R/PFOA"
setwd(HOME)

newday <- file.path('C:/Users/pacho003/OneDrive - Wageningen University & Research/C Channel/R/PFOA/Results', Sys.Date())
dir.create(newday)

#*****************************************************************************#
#---Run script----------------------------------------------------------------#
#*****************************************************************************#

#----Extracting input---------------------------------------------------------#
# Script for calculated input
INPUT <- read.csv(here("InputData/CALCULATED.INPUT.csv"), header=1) %>%  #The function for calculating input can be called here: source(here("Script/ZENPBK_Input"))                                         #call input script
  pivot_wider(names_from = "Parameter", values_from = "Value")

#----Set Parameters for the simulation----------------------------------------#
STARTtime <- 0 # start of simulation (h)
STOPtime <- 24 # end of simmulation (h) (can be from 24h to years)
Dt <- 1/60 # steps of simulation (now it per minute) (could be 1/100, 1/200 or 1 for simulations over years)
times <- seq(STARTtime,STOPtime,by=Dt)

#Exposure
PO_CONC=1 #mg/kg day, Videmann
PO_DOSE=PO_CONC*INPUT$BW/INPUT$MW*1000 #umol


#Model parameters 
PARAMETERS=c(VInL = INPUT$VInL,
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
             Km_glu = INPUT$,  #affinity constant zearalenone to UGT
             VmaxMet.Li, #liver glucuronidation
             VmaxExcIn.glu, #active transport from intestinal tissue to lumen
             KmExc.glu, #affinity constant zearalenone glucuronide to efflux transporter
             VmaxAbsIn.glu, #active transport from lumen to tissue
             KmAbs.glu #affinity constant zearalenone glucuronide to influx transporer
             ) 



#Initial state
INITIALstate=c(
  dAInL = 0.0, 
  dAIn = 0.0, 
  dAFe = 0.0, 
  dALi = 0.0, 
  dAKi = 0.0, 
  dAUr = 0.0, 
  dAAd = 0.0, 
  dASl = 0.0, 
  dARa = 0.0, 
  dABl = 0.0,
  dAIn.glu = 0.0, 
  dAInL.glu = 0.0, 
  dAFe.glu = 0.0, 
  ALi.glu = 0.0, 
  dABi.glu = 0.0, 
  dAKi.glu = 0.0, 
  dAUr.glu = 0.0, 
  dAAd.glu = 0.0, 
  dASl.glu = 0.0, 
  dARa.glu = 0.0, 
  dABl.glu = 0.0)

#----Simulation---------------------------------------------------------------#
PBK <- ode(y=INITIALstate,times=times, func=PBK_ZENandGLUC,parms=PARAMETERS, method="lsoda")
head(PBK)

#----Results------------------------------------------------------------------#
plot(PBK)

OUTPUTamounts <- as.data.frame(PBK) #(umol)






