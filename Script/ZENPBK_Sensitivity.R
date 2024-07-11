
# Libraries-------------------------------------------------------------------#

library(deSolve)
library(tidyverse)
library(ggplot2)
library(here)
library(showtext)
font_add(family = "Garamond", regular = "GARA.TTF")
showtext_auto()
library(data.table)
library(openxlsx)
library(writexl)
library(sensitivity)
library(pksensi)
library(PKNCA)

# Working directory-----------------------------------------------------------#
WhatAmITesting <- "AddLiverMet"
TodaysOuput <- here("Output/Data", Sys.Date(), WhatAmITesting)
dir.create(TodaysOuput, WhatAmITesting,recursive = TRUE)
PlotFolderIs <- here("Output/Plots", Sys.Date(), WhatAmITesting)
dir.create(PlotFolderIs, WhatAmITesting,recursive = TRUE)


#*****************************************************************************#
#---Sensitivity analysis script-----------------------------------------------#
#*****************************************************************************#


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

# Exposure
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
             kTIn = INPUT$kTIn,
             kTge = INPUT$kTge, #gastric emptying rate, assumed to mimic bile secretion rate
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
             kHydro = 1000 # /hhydrolysis constant; assumption as in Veronique
) 

# Initial state
INITIALstate = c(
  AInL = PO_DOSE, 
  # AInL.f = 0,
  # AInL.s = 0,
  AIn = 0, 
  AFe = 0,
  ALi = 0,
  # ABi = 0,
  AKi = 0,
  AUr = 0,
  AAd = 0,
  ASl = 0,
  ARa = 0,
  ABl = 0,
  AIn.glu = 0,
  AInL.glu = 0,
  AFe.glu = 0,
  ABl.glu = 0,
  ALi.glu = 0,
  ABi.glu = 0,
  AKi.glu = 0,
  AUr.glu = 0,
  AAd.glu = 0,
  ASl.glu = 0,
  ARa.glu = 0
)

#----PBK function-------------------------------------------------------------#
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
    dAIn <-  + ka*AInL - VmaxMet.In*CIn/(KmGluc.In+CIn) - QPo*CvIn + QIn*CBl   #+ VmaxAbs*CInL/(KmAbs+CInL) - VmaxSec*CIn/(KmSec+CIn) 
    #        fromLumenActive      toglucuronide                     fromBlood toLiver             toLumen                   fromLum                         
    
    # # Feces
    dAFe <- kTIn*AInL
    #       fromIntestine
    
    # Liver
    dALi <- QPo*CvIn + QHa*CBl - QLi*CvLi - VmaxMet.Li*CLi/(KmGluc.Li+CLi)
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
    
    # Intestinal lumen
    dAInL.glu <- AIn.glu - kHydro*AInL.glu - kTIn*AInL.glu + kTge*CBi.glu #+ VmaxExc.In.glu*CIn.glu/(KmGluc.Exc+CIn.glu) - VmaxAbs.In.glu*CInL.glu/(KmGluc.Abs+CInL.glu)  #(umol/h)
    #            fromInTis toZen             toFeces         fromBile
    # Intestine
    dAIn.glu <- VmaxMet.In*CIn/(KmGluc.In+CIn) - AIn.glu + QIn*CBl.glu - QPo*CvIn.glu # -VmaxExc.In.glu*CIn.glu/(KmGluc.Exc+CIn.glu) + VmaxAbs.In.glu*CInL.glu/(KmGluc.Abs+CInL.glu)  #(umol/h)
    #           fromIntZen                       toInL     fromblood     toportal                                        fromInL                                         
    
    # Feces
    dAFe.glu <- kTIn*AInL.glu
    #           fromInLumen
    
    # Blood
    dABl.glu <- QLi*CvLi.glu + QKi*CvKi.glu + QAd*CvAd.glu + QSl*CvSl.glu + QRa*CvRa.glu +
      - ((QIn + QHa + QKi + QAd + QSl + QRa) * CBl.glu)
    
    # Liver
    dALi.glu <- VmaxMet.Li*CLi/(KmGluc.Li+CLi) + QPo*CvIn.glu - QLi*CvLi.glu + QHa*CBl.glu - CL_Bi*CLi.glu
    #           fromZen                          fromPortal     toBlood        fromHepaticA  toBile
    
    # Bile
    dABi.glu <- CL_Bi*CLi.glu - kTge*CBi.glu
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
    
    
    
    #----Output---------------------------------------------------------------#
    
    return(list(c(dAInL, dAIn, dAFe, 
                  dALi, dAKi, dAUr, dAAd, dASl, dARa, dABl,
                  dAIn.glu, dAInL.glu, 
                  dAFe.glu,
                  dABl.glu, dALi.glu, dABi.glu, dAKi.glu, dAUr.glu, dAAd.glu,
                  dASl.glu, dARa.glu
    ),  CInL, #22, 
    CIn, #23, 
    CLi, #24
    CKi, #25
    CAd, #26
    CSl, #27
    CRa, #28
    CBl, #29
    CInL.glu, #30
    CBi.glu, #31
    CIn.glu, #32
    CLi.glu, #33
    CKi.glu, #34
    CAd.glu, #35
    CSl.glu, #36
    CRa.glu, #37
    CBl.glu #38
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


OUTPUT.A <- as.data.frame(PBK) #(umol)
write.csv(OUTPUT.A, file = here(TodaysOuput, "OUTPUT.A.csv"), row.names=FALSE)

outputs = c("CBl") #This will be changed to measure the uncertainty/sensitivity in other organs ?!

# Mass Balance
CALCULATED <- c( #!Change accordingly depending on the compartments in the function#
  "AInL",
  # "AInL.f",
  # "AInL.s",
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
  "ABl.glu",
  "ALi.glu",
  "ABi.glu",
  "AKi.glu",
  "AUr.glu",
  "AAd.glu",
  "ASl.glu",
  "ARa.glu"
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
ggsave(here(TodaysOuput, "MB_ERROR.png"), dpi = 300)             


#*****************************************************************************#
#----Sensitivity analysis using pksensi---------------------------------------#


## Define the distribution of the parameters that you will analyse in the sensitivity test 
## q: vector that contains all the parameters to be tested and defines their distribution 
## For now the distribution is assumed to be unified for all of them, but I propably need to change this
q <- c("qunif", # VInL
       "qunif", # VBi
       "qunif", # VIn
       "qunif", # VLi
       "qunif", # VKi
       "qunif", # VAd
       "qunif", # VSl
       "qunif", # VRa
       "qunif", # VBl
       "qunif", # PInb
       "qunif", # PLib
       "qunif", # PKib
       "qunif", # PAdb
       "qunif", # PSlb
       "qunif", # PRab
       "qunif", # fub
       "qunif", # PInb.glu
       "qunif", # PLib.glu
       "qunif", # PKib.glu
       "qunif", # PAdb.glu
       "qunif", # PSlb.glu
       "qunif", # PRab.glu
       "qunif", # fub.glu 
       "qunif", # QIn
       "qunif", # QPo
       "qunif", # QLi
       "qunif", # QHa
       "qunif", # QKi
       "qunif", # QAd
       "qunif", # QSl
       "qunif", # QRa
       "qunif", # ka
       "qunif", # kTIn
       "qunif", # kTge
       "qunif", # CL_Re
       "qunif", # CL_Bi
       "qunif", # VmaxMet.In
       "qunif", # KmGluc.In
       "qunif", # VmaxMet.Li
       "qunif", # KmGluc.Li
       "qunif"  # kHydro
       )


## Set parameter distribution ##
# we use 10% change in all parameters

LL <- 0.9 # 10% lower limit
UL <- 1.1 # 10% upper limit

## q.arg: lists all the values of each parameter to be tested
q.arg <- list(list(min = PARAMETERS["VInL"]*LL, max= PARAMETERS["VInL"]*UL),
              list(min = PARAMETERS["VBi"]*LL, max= PARAMETERS["VBi"]*UL),
              list(min = PARAMETERS["VIn"]*LL, max =  PARAMETERS["VIn"]*UL),
              list(min = PARAMETERS["VLi"]*LL, max =  PARAMETERS["VLi"]*UL),
              list(min = PARAMETERS["VKi"]*LL, max =  PARAMETERS["VKi"]*UL),
              list(min = PARAMETERS["VAd"]*LL, max = PARAMETERS["VAd"]*UL),
              list(min = PARAMETERS["VSl"]*LL, max = PARAMETERS["VSl"]*UL),
              list(min = PARAMETERS["VRa"]*LL, max = PARAMETERS["VRa"]*UL),
              list(min = PARAMETERS["VBl"]*LL, max = PARAMETERS["VBl"]*UL),
              list(min = PARAMETERS["PInb"]*LL, max = PARAMETERS["PInb"]*UL),
              list(min = PARAMETERS["PLib"]*LL, max = PARAMETERS["PLib"]*UL),
              list(min = PARAMETERS["PKib"]*LL, max = PARAMETERS["PKib"]*UL),
              list(min = PARAMETERS["PAdb"]*LL, max = PARAMETERS["PAdb"]*UL),
              list(min = PARAMETERS["PSlb"]*LL, max = PARAMETERS["PSlb"]*UL),
              list(min = PARAMETERS["PRab"]*LL, max = PARAMETERS["PRab"]*UL),
              list(min = PARAMETERS["fub"]*LL, max = PARAMETERS["fub"]*UL),
              list(min = PARAMETERS["PInb.glu"]*LL, max = PARAMETERS["PInb.glu"]*UL),
              list(min = PARAMETERS["PLib.glu"]*LL, max = PARAMETERS["PLib.glu"]*UL),
              list(min = PARAMETERS["PKib.glu"]*LL, max = PARAMETERS["PKib.glu"]*UL),
              list(min = PARAMETERS["PAdb.glu"]*LL, max = PARAMETERS["PAdb.glu"]*UL),
              list(min = PARAMETERS["PSlb.glu"]*LL, max = PARAMETERS["PSlb.glu"]*UL),
              list(min = PARAMETERS["PRab.glu"]*LL, max = PARAMETERS["PRab.glu"]*UL),
              list(min = PARAMETERS["fub.glu"]*LL, max = PARAMETERS["fub.glu"]*UL), 
              list(min = PARAMETERS["QIn"]*LL, max = PARAMETERS["QIn"]*UL),
              list(min = PARAMETERS["QPo"]*LL, max = PARAMETERS["QPo"]*UL),
              list(min = PARAMETERS["QLi"]*LL, max = PARAMETERS["QLi"]*UL),
              list(min = PARAMETERS["QHa"]*LL, max = PARAMETERS["QHa"]*UL),
              list(min = PARAMETERS["QKi"]*LL, max = PARAMETERS["QKi"]*UL),
              list(min = PARAMETERS["QAd"]*LL, max = PARAMETERS["QAd"]*UL),
              list(min = PARAMETERS["QSl"]*LL, max = PARAMETERS["QSl"]*UL),
              list(min = PARAMETERS["QRa"]*LL, max = PARAMETERS["QRa"]*UL),
              list(min = PARAMETERS["ka"]*LL, max = PARAMETERS["ka"]*UL),
              list(min = PARAMETERS["kTIn"]*LL, max = PARAMETERS["kTIn"]*UL),
              list(min = PARAMETERS["kTge"]*LL, max = PARAMETERS["kTge"]*UL),
              list(min = PARAMETERS["CL_Re"]*LL, max = PARAMETERS["CL_Re"]*UL),
              list(min = PARAMETERS["CL_Bi"]*LL, max = PARAMETERS["CL_Bi"]*UL),
              list(min = PARAMETERS["VmaxMet.In"]*LL, max = PARAMETERS["VmaxMet.In"]*UL),
              list(min = PARAMETERS["KmGluc.In"]*LL, max = PARAMETERS["KmGluc.In"]*UL),
              list(min = PARAMETERS["VmaxMet.Li"]*LL, max = PARAMETERS["VmaxMet.Li"]*UL),
              list(min = PARAMETERS["KmGluc.Li"]*LL, max = PARAMETERS["KmGluc.Li"]*UL),
              list(min = PARAMETERS["kHydro"]*LL, max = PARAMETERS["kHydro"]*UL)
)



## Create parameter matrix ##  
params <- c("VInL",
            "VBi",
            "VIn",
            "VLi",
            "VKi",
            "VAd",
            "VSl",
            "VRa",
            "VBl",
            "PInb",
            "PLib",
            "PKib",
            "PAdb",
            "PSlb",
            "PRab",
            "fub",
            "PInb.glu",
            "PLib.glu",
            "PKib.glu",
            "PAdb.glu",
            "PSlb.glu",
            "PRab.glu",
            "fub.glu", 
            "QIn",
            "QPo",
            "QLi",
            "QHa",
            "QKi",
            "QAd",
            "QSl",
            "QRa",
            "ka",
            "kTIn",
            "kTge",
            "CL_Re",
            "CL_Bi",
            "VmaxMet.In",
            "KmGluc.In",
            "VmaxMet.Li",
            "KmGluc.Li",
            "kHydro")

length(params)==length(q)  # Check if lengths of params and q match

set.seed(1234) # For reproducibility ?! not sure what this is
x <- rfast99(params = params, n = 200, q = q, q.arg = q.arg, rep = 10, conf = 0.95) # Generate the parameter matrix
dim(x$a) # The array of c(model evaluation, replication, parameters)


## Conduct simulation ##
# out <- solve_fun(x, time=times, func = PBPKmodPFOA, initState = yini, outnames = outputs)

out <- solve_fun(x,  
           times = times,  
           func = PBK_ZENandGLUC,
           initState = INITIALstate, 
           outnames = outputs
           )
# head(PBK)


saveRDS(object=out, file="out_scaled.rds") # Save the output
#out <- readRDS("out_scaled.rds")

## Output of the Uncertainty analysis ##
pdf("out_scaled.pdf") # Open a PDF device
pksim(out) # Output the results of the uncertainty analysis
dev.off() # Close the PDF device

## Output from the sensitivity analysis ##
pdf("out_scaled.pdf")
plot(out)
dev.off()

ResultsSI <- as.data.frame(print(out["tSI"])) # Extract and convert sensitivity indices to a data frame
ResultsSI$Times <- rownames(ResultsSI)  # Get times as row names
ResultsSI$Times <- as.numeric(as.character(ResultsSI$Times)) # Convert times to numeric

write.xlsx(ResultsSI,
           file =file.path(newday,"ResultsSI.xlsx"),
           colNames = TRUE, borders = "rows"
)

write.xlsx(ResultsSI,
           file = "Results/2023-09-24/ResultsSI.xlsx",
           colNames = TRUE
)

check(out) #Perform checks on the output?!
pdf("heat_check_CI_scaled.pdf") #Open a PDF device
heat_check(out, index = "CI") #Plot all heat checks
dev.off() #Close the PDF device

pdf("heat_check_all.pdf") #Open a PDF device
heat_check(out, show.all = TRUE) #Plot all heat checks
dev.off() #Close the PDF device


