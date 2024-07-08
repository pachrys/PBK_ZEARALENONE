# This scipt is made to easily calculate the parameters needed for the PBK.
# This can be used to add/delete tissues, recalculate partition coefficients
# organ volumes etc. 
# Absoprtion, metabolism and transport parameters are also calculated here.

library(here)
library(tidyverse)

#*****************************************************************************#
#---Script to calculate input-------------------------------------------------#
#*****************************************************************************#

#----Physiology---------------------------------------------------------------#

BW <- 0.250 #(kg)
QC <- 80 * 60 / 1000 #(L/h)

Rat.input <- read.csv(here("InputData/Rat.input.csv"), header=1) #tissue composition mainly as R&R composition, unified by Utsey and optimised as per Musther
Rat.dat <- Rat.input %>% filter(!tissue %in% c("RBCs", "Plasma")) #df for all tissues except for adipose and RBCs

# Flow rates
## Extracting fractional flow rate (constant)
names <- substr(Rat.dat$tissue, 1, 2)
F_QC.names <- paste0("FQ", names)
F_QC.dat <- cbind(F_QC.names, select(Rat.dat, F_QC))

## Calculating tissue blood flow rates (L/h)
QC.dat <- F_QC.dat %>%
  mutate(F_QC = F_QC * QC) %>% #QC = QCconstant*QC
  rename(QTissue = F_QC.names) %>%
  rename(QValue = F_QC) %>%
  mutate(QTissue = sub("^F", "", QTissue)) %>%
  pivot_wider(names_from = QTissue, values_from = QValue) %>%
  as.data.frame()

# Tissue volumes
## Extracting fractional volumes (constant)
F_Vo.names <- paste0("FV", names)
F_Vo.dat <- cbind(F_Vo.names, select(Rat.dat, F_Vo))

## Calculating tissue volumes (L) (assumption tissue density of 1 kg/L)            fractional volumes come from Musther (insignificant differences from Brown, Punt, Ruark)
Vo.dat <- F_Vo.dat %>%
  mutate(F_Vo = F_Vo * BW) %>%
  rename(VTissue = F_Vo.names) %>%
  rename(VValue = F_Vo) %>%
  mutate(VTissue = sub("F", "", VTissue)) %>%
  pivot_wider(names_from = VTissue, values_from = VValue) %>%
  as.data.frame()



#----Physicochemical parameters-----------------------------------------------#

ZEN.dat <- read.csv(here("InputData/Zen.input.csv"), header=1)     #from {Agahi, 2020 #456} and chemicalize, BP from Shin, Hong et al.2009


# Tissue partitioning (based on Schmitt)

## Calculating fraction unbound in plasma #Reference: Lobell M, Sivarajah V. In silico prediction of aqueous solubility, human plasma protein binding and volume of distribution of compounds from calculated pKa and AlogP98 values. Mol Divers. 2003;7(1):69-87. doi: 10.1023/b:modi.0000006562.93049.36. PMID: 14768905.
fup <- 1 / (1 + 10^(0.4485 * ZEN.dat$logP - 0.4782))

## Calculating ionisation status
## Assumption that Zearalenone is a monoprotic acid; in reality it's a diprotic acid with pKa2=10.38. Insignificant difference in tissue:plasma partition coefficients was observed when diprotic acid equations were tested. Check.As.Diprotic.RnR.R file.
## Zearalenone is not ionized at physiological pH. RnR's equation for neutrals was overestimating the partition in all tissues by 7%. Check As.Neutral.RnR.R
K_n_pl <- 10^ZEN.dat$logP #neutral phospholipids:water partition coefficient (logP here is logMA of Shmitt)
K_protein <- ((0.81 + 0.11 * K_n_pl) / 24.92) * 5 #protein:water partition; Schmitt, Walter (2008)
pH <- Rat.dat$pH #tissue pH
alpha <- 1e-3 #ratio between ditribution coefficient at given pH (D) and that in neutral form (D0)

W <- 10^(Rat.dat$pH - ZEN.dat$pKa) #Calculates ZEN ionisation status depending on tissue pH
K_n_l <- K_n_pl * (((1 - alpha) / (1 + W)) + alpha) #Calculates neutral lipids:neutral phospholipids partition coefficient
K_a_pl <- K_n_pl * ((1 / (1 + W)) + 0.05 * (1 - (1 / (1 + W)))) #Calculates acid phospholipids:neutral lipids partition coefficient based on ionisation status

## Calculating Tissue:Plasma partition coefficients
kp <- (Rat.dat$f_water +
         (K_n_l * Rat.dat$f_n_l) +
         (K_n_pl * Rat.dat$f_n_pl) +
         (K_a_pl * Rat.dat$f_a_pl) +
         (K_protein * Rat.dat$f_proteins)) * fup

PC.names <- paste0("P", names)
PC.dat <- data.frame(PC.names, kp = as.numeric(kp))


## Calculate Tissue:Blood partition coefficient 
PCB.dat <- PC.dat %>%
  mutate(kp = kp * ZEN.dat$BP) %>%
  rename(PValue = kp) %>%
  rename(PTissue = PC.names) %>%
  pivot_wider(names_from = PTissue, values_from = PValue) %>%
  as.data.frame()


# Calculate free fraction in blood
fub = fup * ZEN.dat$BP 



#----Toxicokinetic parameters------------------------------------------------#

# Gastric emptying rate (1/h)
kTge <- 1 / ZEN.dat$kTge

# Intestinal transit times (1/h)
kTIn <- 1 / (ZEN.dat$kTsi + ZEN.dat$kTco) #total intestinal transit time including small and large intestine

# Absorption constant (1/h)
## TO CHANGE!!! Scaling in vitro to in vivo absorption parameters 
Peff.human <- (10^(0.4926 * log10(ZEN.dat$Papp) - 0.1454)) #human Peff(cm^-4/s)
Peff.rat <- Peff.human / 11.04 #rat Peff(cm^-4/s)
Peff <- Peff.rat / 10000 * 3600 #conversion to cm/h
ka = Peff * 2 / ZEN.dat$R #/h  

# Clearance (L/h)
CL_Bi <- ZEN.dat$CLBi * 60 / 1000 #L/h, biliary clearance; scaling from ml/min to L/h; CHANGED FROM 0.4 TO 6.3; BASED ON BAOJIAN WU DOI 10.1002/jps.22827
CL_Re <- ZEN.dat$CLRe * 60 / 1000 #L/h, renal clearance; scaling to L/h


# Metabolism 

## Intestine
VmaxGluc.SI <- 8.617 #nmol/min/mg; S9 incubations Diana ;other input 5nmol/min/mg
KmGluc.In <- 6.198 #umol/L; S9 incubations Diana
VmaxMet.In <- VmaxGluc.SI * ZEN.dat$SF.In * Vo.dat$VGu * 1000 * 60 / 1000 #umol/h, scaling from nmol/min/mg protein to nmol/min for whole liver; equation from {Grandoni, 2019 #457}; in vivo = in vitro * scaling to microsomal protein per gram of rat liver * scaling to total weight (assuming density=1) of rat intestine; SF = 15.5mg microsomal protein/g gut ({Musther, 2017 #458})


## Liver
VmaxGluc.LI <- 7.029 #nmol/min/mg; S9 incubations Diana ;other input 5nmol/min/mg
KmGluc.Li <- 6.753 #umol/L; S9 incubations Diana
VmaxMet.Li <- ZEN.dat$VmaxMet.Li * ZEN.dat$SF.Li * Vo.dat$VLi * 1000 * 60 / 1000 #umol/h, scaling from nmol/min/mg protein to nmol/min for whole liver; equation from {Grandoni, 2019 #457}; in vivo = in vitro * scaling to microsomal protein per gram of rat liver * scaling to total weight (assuming density=1) of rat liver; SF = 46mg microsomal protein/g liver ({Musther, 2017 #458})


# Active Transport
# Scaling in vitro to in vivo transport parameters; according to {Videmann, 2008 #384}; input values are expressed in pmol/h and uM; cell concentration for transport experiments is not mentioned!!

## Absorption (intestinal lumen to intestinal tissue)
VmaxGlucAbs.In <- 40.1 #pmol/h Videmann
KmGluc.Abs <- 23.6 #umol/L
VmaxAbs.In.glu <- VmaxGlucAbs.In / 1000000 #umol/h, scaling from pmol/h


## Excretion (intestinal tissue to intestinal lumen)
VmaxGluc.Exc.In <- 33.3 #pmol/h
KmGluc.Exc <- 40 #umol/L
VmaxExc.In.glu <- VmaxGluc.Exc.In / 1000000 #umol/h, scaling from pmol/h




#----Extracting all data-----------------------------------------------------#

F_Vo.dat <- F_Vo.dat %>%
  pivot_wider(names_from = F_Vo.names, values_from = F_Vo)

MW = ZEN.dat$MW
logP = ZEN.dat$logP
BW = BW
QC = QC
QIn = QC.dat$QGu #Intestinal blood flow
QPo = QC.dat$QGu + QC.dat$QSp #Portal blood flow = Σ small intestine, colon, pancreas, stomach, spleen
QLi = QC.dat$QLi #Liver blood flow (this includes portal! to be changed later! to add Q hepatic artery only)
QKi = QC.dat$QKi #Kidney blood flow
QAd = QC.dat$QAd #Adipose blood flow
QSl = QC.dat$QBo + QC.dat$QBr + QC.dat$QMu + QC.dat$QSk #Slowly perfused blood flow
QRa = QC - (QIn + QLi + QKi + QAd + QSl) #Rapidly perfused blood flow
MassBalanceQC = QC - (QIn + QLi + QKi + QAd + QSl + QRa)  #check blood flow mass balance; should be 0
VIn = Vo.dat$VGu #Intestinal volume
VLi = Vo.dat$VLi #Liver volume
VKi = Vo.dat$VKi #Kidney volume
VAd = Vo.dat$VAd #Adipose volume
VSl = Vo.dat$VBo + Vo.dat$VBr + Vo.dat$VMu + Vo.dat$VSk #Slowly perfused volume
VRa = Vo.dat$VHe + Vo.dat$VLu + Vo.dat$VSp #Rapidly perfused volume
VBl = BW*0.06 # 0.06 is the fractional volume of blood based on Musther
VBi = BW*0.00023495 #bile volume total volume was 46.97 +- 4.21ul; male Fisher rats 180-220g => fractional volume of bile is 0.00023495 https://ncbi.nlm.nih.gov/pmc/articles/PMC1891969/
VInL = 2/1000 #total volume of lumen = 2ml
PIn = PCB.dat$PGu #Intestine:Blood partition coefficient
PLi = PCB.dat$PLi #Liver:Blood partition coefficient
PKi = PCB.dat$PKi #Kidney:Blood partition coefficient
PAd = PCB.dat$PAd #Adipose:Blood partition coefficient
PSl = (F_Vo.dat$FVBo*PCB.dat$PBo + F_Vo.dat$FVBr*PCB.dat$PBr +
         F_Vo.dat$FVMu*PCB.dat$PMu + F_Vo.dat$FVSk*PCB.dat$PSk) /
  (F_Vo.dat$FVBo + F_Vo.dat$FVBr + F_Vo.dat$FVMu + F_Vo.dat$FVSk) #slowly perfused = partition coefficient of each of these organs * fractional volume of the respective organ / sum of the fractional volume of all these organs
PRa = ZEN.dat$BP*(F_Vo.dat$FVHe*PCB.dat$PHe + F_Vo.dat$FVLu*PCB.dat$PLu +
                    F_Vo.dat$FVSp*PCB.dat$PSp) /
  (F_Vo.dat$FVHe + F_Vo.dat$FVLu + F_Vo.dat$FVSp)
fup = fup #free fraction plasma
fub = fub #free fraction blood
kTIn = kTIn #intestinal transit time (1/h)
kTge = kTge #gastric emptying rate (1/h)
ka = ka #intestinal absorption constant (1/h)
CL_Bi = CL_Bi #biliary clearance (L/h)
CL_Re = CL_Re #renal clearance (L/h)
KmGluc.In = KmGluc.In #affinity constant intestinal UGT (umol/L)
VmaxMet.In = VmaxMet.In #max metab. velocity intestinal glucuronidation (umol/h)
KmGluc.Li = KmGluc.Li #affinity constant liver UGT (umol/L)
VmaxMet.Li = VmaxMet.Li #max metab. velocity liver glucuronidation (umol/h)
KmGluc.Abs = KmGluc.Abs #affinity constant influx transporter (umol/L)
VmaxAbs.In.glu = VmaxAbs.In.glu #max transport velocity to intestinal tissue (umol/h)
KmGluc.Exc = KmGluc.Exc #affinity constant efflux transporter (umol/L)
VmaxExc.In.glu = VmaxExc.In.glu #max transport velocity from intestinal tissue (umol/h)


# MW = ZEN.dat$MW 
# logP = ZEN.dat$logP
# BW = BW
# QC = QC
# QIn = QC.dat$QGu #Intestinal blood flow
# QPo = QC.dat$QGu + QC.dat$QSp #Portal blood flow = Σ small intestine, colon, pancreas, stomach, spleen
# QLi = QC.dat$QLi #Liver blood flow (this includes portal! to be changed later! to add Q hepatic artery only)                                           
# QKi = QC.dat$QKi #Kidney blood flow                                                 
# QAd = QC.dat$QAd #Adipose blood flow                                                 
# QSl = QC.dat$QBo + QC.dat$QBr + QC.dat$QMu + QC.dat$QSk #Slowly perfused blood flow   
# QRa = QC - (QIn + QLi + QKi + QAd + QSl) #Rapidly perfused blood flow                                    
# MassBalanceQC = QC - (QIn + QLi + QKi + QAd + QSl + QRa)  #check blood flow mass balance; should be 0
# VIn = Vo.dat$VGu #Intestinal volume                                           
# VLi = Vo.dat$VLi #Liver volume                                              
# VKi = Vo.dat$VKi #Kidney volume                                              
# VAd = Vo.dat$VAd #Adipose volume                                            
# VSl = Vo.dat$VBo + Vo.dat$VBr + Vo.dat$VMu + Vo.dat$VSk #Slowly perfused volume 
# VRa = Vo.dat$VHe + Vo.dat$VLu + Vo.dat$VSp #Rapidly perfused volume
# VBl = BW*0.06 # 0.06 is the fractional volume of blood based on Musther 
# VBi = BW*0.00023495 #bile volume total volume was 46.97 +- 4.21ul; male Fisher rats 180-220g => fractional volume of bile is 0.00023495 https://ncbi.nlm.nih.gov/pmc/articles/PMC1891969/
# VInL = 2/1000 #total volume of lumen = 2ml 
# PIn = 6.51 #Intestine:Blood partition coefficient
# PLi = 4.54 #Liver:Blood partition coefficient
# PKi = 4.18 #Kidney:Blood partition coefficient                                         
# PAd = 24.06 #Adipose:Blood partition coefficient
# PSl = (F_Vo.dat$FVBo*5.42 + F_Vo.dat$FVBr*10.60 + 
#          F_Vo.dat$FVMu*2.66 + F_Vo.dat$FVSk*5.72) /
#   (F_Vo.dat$FVBo + F_Vo.dat$FVBr + F_Vo.dat$FVMu + F_Vo.dat$FVSk) #slowly perfused = partition coefficient of each of these organs * fractional volume of the respective organ / sum of the fractional volume of all these organs
# PRa = ZEN.dat$BP*(F_Vo.dat$FVHe*3.57 + F_Vo.dat$FVLu*5.17 + 
#                     F_Vo.dat$FVSp*2.49) /
#   (F_Vo.dat$FVHe + F_Vo.dat$FVLu + F_Vo.dat$FVSp)
# fup = fup #free fraction plasma
# fub = fub #free fraction blood
# kTIn = kTIn #intestinal transit time (1/h)
# kTge = kTge #gastric emptying rate (1/h)
# ka = ka #intestinal absorption constant (1/h)
# CL_Bi = CL_Bi #biliary clearance (L/h)
# CL_Re = CL_Re #renal clearance (L/h)
# KmGluc.In = KmGluc.In #affinity constant intestinal UGT (umol/L)
# VmaxMet.In = VmaxMet.In #max metab. velocity intestinal glucuronidation (umol/h)
# KmGluc.Li = KmGluc.Li #affinity constant liver UGT (umol/L)
# VmaxMet.Li = VmaxMet.Li #max metab. velocity liver glucuronidation (umol/h)
# KmGluc.Abs = KmGluc.Abs #affinity constant influx transporter (umol/L)
# VmaxAbs.In.glu = VmaxAbs.In.glu #max transport velocity to intestinal tissue (umol/h)
# KmGluc.Exc = KmGluc.Exc #affinity constant efflux transporter (umol/L)
# VmaxExc.In.glu = VmaxE-xc.In.glu #max transport velocity from intestinal tissue (umol/h)




#----Saving-------------------------------------------------------------------#

INPUT.parameters <- cbind(
  MW, 
  logP,
  BW,
  QC,
  QIn,
  QPo,
  QLi,   
  QKi, 
  QAd, 
  QSl, 
  QRa, 
  VIn, 
  VLi, 
  VKi, 
  VAd, 
  VSl, 
  VRa, 
  VBl, 
  VBi, 
  VInL, 
  PIn, 
  PLi, 
  PKi,
  PAd,
  PSl, 
  PRa,
  fup, 
  fub, 
  kTIn,
  kTge,
  ka, 
  CL_Bi, 
  CL_Re, 
  KmGluc.In, 
  VmaxMet.In,
  KmGluc.Li, 
  VmaxMet.Li,
  KmGluc.Abs,
  VmaxAbs.In.glu, 
  KmGluc.Exc, 
  VmaxExc.In.glu
) %>% as.data.frame() %>% 
  pivot_longer(names_to = "Parameter", values_to = "Value", MW:VmaxExc.In.glu)

write.csv(INPUT.parameters, file = here("InputData/", "MY.CALCULATED.INPUT.csv"), row.names=FALSE)
