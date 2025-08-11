## Admin Info ####
# ------------------------------------------------------ #
# Script originally written by: Chrysanthi Pachoulide
# User: Chrysanthi Pachoulide
# Last updated: August 2025

## Packages needed ####
# ------------------------------------------------------ #
library(here)
library(tidyverse)

# Input
logP_ZEN <- 3.580 # Molgpka
pKa1_ZEN <- 7.5 # Molgpka
pka2_ZEN <- 7.7 # Molgpka
BP_ZEN <- 0.67 # Shin et al. 2014 DOI: 10.1080/15287390903212741
fup_ZEN <- 0.156 # Mukherjee et al. 2014 doi:10.1371/journal.pone.0113632

MW_aZEL <- 320.385 # PubMed
logP_aZEL <- 3.371 # Molgpka
pKa1_aZEL <- 7.5 # Molgpka
pKa2_aZEL <- 7.7 # Molgpka

# Zearalenone 14-glucuronide
MW_ZENGLU <- 494.5 # PubMed
logP_ZENGLU <- 1.145 # Molgpka
pKa_ZENGLU <- 3 # Molgpka

## Tissue:Blood partitioning ####
## ------------------------------------------------------

# Tissue partitioning (Rodgers and Rowland as coded by kmutsey: https://github.com/metrumresearchgroup/PBPK_PC/blob/master/script/CalcKp_P%26T.R
# Assumption that Zearalenone and a-Zearalenol are diprotic acids and Zearalenone glucuronide a monoprotic acid

dat_all <- Rat.input %>% filter(!tissue %in% c("RBCs", "Adipose", "Plasma"))  #df for all tissues except for adipose, RBCs, and plasma
dat_ad <- Rat.input %>% filter(tissue == "Adipose")  #df for adipose
dat_rbc <- Rat.input %>% filter(tissue == "RBCs") #df for RBCs
dat_plas <- Rat.input %>% filter(tissue == "Plasma") #df for aplasma

pH_IW <- 7       #pH of intracellular tissue water
pH_P <- 7.4      #pH of plasma
pH_RBC <- 7.22    #pH of blood cells

## Calculate Kp for ZEN ####

logP_ZEN <- 3.580 # Molgpka
pKa1_ZEN <- 7.5 # Molgpka
pKa2_ZEN <- 7.7 # Molgpka
BP_ZEN <- 0.67 # Shin et al. 2014 DOI: 10.1080/15287390903212741
fup_ZEN <- 0.156 # Mukherjee et al. 2014 doi:10.1371/journal.pone.0113632

P <- 10^(logP_ZEN)   # octonal:water partition coeff
logP_OW <- 1.115*logP_ZEN - 1.35 #oil:water partition coeff
P_OW <- 10^(logP_OW) 
# Ka <- 10^(-pKa1_ZEN)
HCT <- 0.45 #hematocrit

Kpu_bc <- (HCT - 1 + BP_ZEN)/(HCT*fup_ZEN)

X <- 10^(pH_IW-pKa1_ZEN)+10^(2*pH_IW-pKa1_ZEN-pKa2_ZEN)       

Y <- 10^(pH_P-pKa1_ZEN)+10^(2*pH_P-pKa1_ZEN-pKa2_ZEN)     

Z <- 1 

Ka_PR <- (1/fup_ZEN - 1 - (P*dat_plas$f_n_l + (0.3*P + 0.7)*dat_plas$f_n_pl)/(1+Y))
Ka_AP <- (Kpu_bc - (1 + Z)/(1 + Y)*dat_rbc$f_iw - (P*dat_rbc$f_n_l + (0.3*P + 0.7)*dat_rbc$f_n_pl)/(1 + Y)) * (1 + Y)/dat_rbc$f_a_pl/Z

Kp_all <- (dat_all$f_ew + ((1 + X)/(1 + Y))*dat_all$f_iw + ((P*dat_all$f_n_l + (0.3*P + 0.7)*dat_all$f_n_pl))/(1 + Y) + (Ka_PR*dat_all$AR*X)/(1 + Y))*fup_ZEN  #non lipid
Kp_ad <- (dat_ad$f_ew + ((1 + X)/(1 + Y))*dat_ad$f_iw + ((P_OW*dat_ad$f_n_l + (0.3*P_OW + 0.7)*dat_ad$f_n_pl))/(1 + Y) + (Ka_PR*dat_ad$AR*X)/(1 + Y))*fup_ZEN  #lipid

nms_all <- dat_all$tissue %>% substr(1,2) %>% tolower()
nms_all <- paste("Kp", nms_all, sep="")
nms <- c("Kpad",nms_all)
Kp_ZEN <- as.list(c(Kp_ad,Kp_all)) #maybe make as dataframe
names(Kp_ZEN) <- nms

fut_ZEN <- 1/(1+((1-fup_ZEN)/fup_ZEN)*0.5) #fraction unbound in the tissues


pH <- 7.4
logD_ZEN <- logP_ZEN-log10(1+10^(pH-pKa1_ZEN)+10^(pH-pKa1_ZEN+pH-pKa2_ZEN))

frac.unionized_ZEN <- 1/(1+10^(pH-pKa1_ZEN)+10^(pH-pKa1_ZEN+pH-pKa2_ZEN))
f.frac.ionized_ZEN <- 1 - 1/(1+10^(pH-pKa1_ZEN)+10^(pH-pKa1_ZEN+pH-pKa2_ZEN))
conc.enzymes.Clint <- 165 # Ans Punt
Fuinc_ZEN <- (f.frac.ionized_ZEN+frac.unionized_ZEN)*(1/(1+conc.enzymes.Clint*10^(0.072*logD_ZEN^2+0.067*logD_ZEN-1.126)))

## Calculate Kp for aZEL ####

MW_aZEL <- 320.385 # PubMed
logP_aZEL <- 3.371 # Molgpka
pKa1_aZEL <- 7.5 # Molgpka
pKa2_aZEL <- 7.7 # Molgpka
BP_aZEL <- 0.67 # assuming the same as ZEN
fup_aZEL <- 0.156 # assuming the same as ZEN

P <- 10^(logP_aZEL)   # octonal:water partition coeff
logP_OW <- 1.115*logP_aZEL - 1.35 #oil:water partition coeff
P_OW <- 10^(logP_OW) 
# Ka <- 10^(-pKa1_aZEL)
HCT <- 0.45 #hematocrit

Kpu_bc <- (HCT - 1 + BP_aZEL)/(HCT*fup_aZEL)

X <- 10^(pH_IW-pKa1_aZEL)+10^(2*pH_IW-pKa1_aZEL-pKa2_aZEL)       

Y <- 10^(pH_P-pKa1_aZEL)+10^(2*pH_P-pKa1_aZEL-pKa2_aZEL)     

Z <- 1 

Ka_PR <- (1/fup_aZEL - 1 - (P*dat_plas$f_n_l + (0.3*P + 0.7)*dat_plas$f_n_pl)/(1+Y))
Ka_AP <- (Kpu_bc - (1 + Z)/(1 + Y)*dat_rbc$f_iw - (P*dat_rbc$f_n_l + (0.3*P + 0.7)*dat_rbc$f_n_pl)/(1 + Y)) * (1 + Y)/dat_rbc$f_a_pl/Z

Kp_all <- (dat_all$f_ew + ((1 + X)/(1 + Y))*dat_all$f_iw + ((P*dat_all$f_n_l + (0.3*P + 0.7)*dat_all$f_n_pl))/(1 + Y) + (Ka_PR*dat_all$AR*X)/(1 + Y))*fup_aZEL  #non lipid
Kp_ad <- (dat_ad$f_ew + ((1 + X)/(1 + Y))*dat_ad$f_iw + ((P_OW*dat_ad$f_n_l + (0.3*P_OW + 0.7)*dat_ad$f_n_pl))/(1 + Y) + (Ka_PR*dat_ad$AR*X)/(1 + Y))*fup_aZEL #lipid

nms_all <- dat_all$tissue %>% substr(1,2) %>% tolower()
nms_all <- paste("Kp", nms_all, sep="")
nms <- c("Kpad",nms_all)
Kp_aZEL <- as.list(c(Kp_ad,Kp_all))
names(Kp_aZEL) <- nms

fut_aZEL <- 1/(1+((1-fup_aZEL)/fup_aZEL)*0.5) #fraction unbound in the tissues

logD_aZEL <- logP_ZEN-log10(1+10^(pH-pKa1_aZEL)+10^(pH-pKa1_ZEN+pH-pKa2_aZEL))
frac.unionized_aZEL <- 1/(1+10^(pH-pKa1_aZEL)+10^(pH-pKa1_aZEL+pH-pKa2_aZEL))
f.frac.ionized_aZEL <- 1 - 1/(1+10^(pH-pKa1_aZEL)+10^(pH-pKa1_aZEL+pH-pKa2_aZEL))
conc.enzymes.Clint <- 165 # Ans Punt
Fuinc_aZEL <- (f.frac.ionized_aZEL+frac.unionized_aZEL)*(1/(1+conc.enzymes.Clint*10^(0.072*logD_aZEL^2+0.067*logD_aZEL-1.126)))



## Calculate Kp for ZENGLUC ####

MW_ZENGLU <- 494.5 # PubMed
logP_ZENGLU <- 1.145 # Molgpka
pKa_ZENGLU <- 3 # Molgpka
BP_ZENGLU <- 0.67 # assuming the same as ZEN
fup_ZENGLU <- 1 / (1 + 10^(0.4485 * logP_ZENGLU - 0.4782)) #Lobell M, Sivarajah V. et al. 2003 doi: 10.1023/b:modi.0000006562.93049.36.

P <- 10^(logP_ZENGLU)   # octonal:water partition coeff
logP_OW <- 1.115*logP_ZENGLU - 1.35 #oil:water partition coeff
P_OW <- 10^(logP_OW) 
Ka <- 10^(-pKa_ZENGLU)
HCT <- 0.45 #hematocrit


#Calculate Kp values
Kpu_bc <- (HCT - 1 + BP_ZENGLU)/(HCT*fup_ZENGLU)

X <- 10^(pH_IW-pKa_ZENGLU)
Y <- 10^(pH_P-pKa_ZENGLU)
Z <- 1

Ka_PR <- (1/fup_ZENGLU - 1 - (P*dat_plas$f_n_l + (0.3*P + 0.7)*dat_plas$f_n_pl)/(1+Y))
Ka_AP <- (Kpu_bc - (1 + Z)/(1 + Y)*dat_rbc$f_iw - (P*dat_rbc$f_n_l + (0.3*P + 0.7)*dat_rbc$f_n_pl)/(1 + Y)) * (1 + Y)/dat_rbc$f_a_pl/Z

# Multiply by fup to get Kp rather than Kpu
Kp_all <- (dat_all$f_ew + ((1 + X)/(1 + Y))*dat_all$f_iw + ((P*dat_all$f_n_l + (0.3*P + 0.7)*dat_all$f_n_pl))/(1 + Y) + (Ka_PR*dat_all$AR*X)/(1 + Y))*fup_ZENGLU  #non lipid
Kp_ad <- (dat_ad$f_ew + ((1 + X)/(1 + Y))*dat_ad$f_iw + ((P_OW*dat_ad$f_n_l + (0.3*P_OW + 0.7)*dat_ad$f_n_pl))/(1 + Y) + (Ka_PR*dat_ad$AR*X)/(1 + Y))*fup_ZENGLU #lipid

nms_all <- dat_all$tissue %>% substr(1,2) %>% tolower()
nms_all <- paste("Kp", nms_all, sep="")
nms <- c("Kpad",nms_all)
Kp_ZENGLU <- as.list(c(Kp_ad,Kp_all))
names(Kp_ZENGLU) <- nms

fut_ZENGLU <- 1/(1+((1-fup_ZENGLU)/fup_ZENGLU)*0.5) #fraction unbound in the tissues

logD_ZENGLU <- logP_ZENGLU - log10(1+10^(pH-pKa_ZENGLU))

frac.unionized_ZENGLU <- 1/(1+10^(pH-pKa_ZENGLU))
f.frac.ionized_ZENGLU <- 1-1/(1+10^(pH-pKa_ZENGLU))
conc.enzymes.Clint <- 165 # Ans Punt
Fuinc_ZENGLU <- 1/(1+conc.enzymes.Clint*10^(0.072*logD_ZENGLU^2+0.067*logD_ZENGLU-1.126)) 

# -----------------------------------------------------------------------------
# Toxicokinetic ####
# -----------------------------------------------------------------------------

# Gastric emptying rate (1/h)
kTge <- 1 / ZEN.dat$kTge

# Intestinal transit times (1/h)
kTIn <- 1 / (ZEN.dat$kTsi + ZEN.dat$kTco) #total intestinal transit time including small and large intestine

# Absorption constant (1/h)
## In vitro 
Peff.human <- (10^(0.4926 * log10(ZEN.dat$Papp) - 0.1454)) #human Peff(cm^-4/s)
Peff.rat <- Peff.human / 11.04 #rat Peff(cm^-4/s)
Peff <- Peff.rat / 10000 * 3600 #conversion to cm/h
ka = Peff * 2 / ZEN.dat$R #/h 

## Insilico  #https://pubmed.ncbi.nlm.nih.gov/17727800/
Papp.cm.power.minus.six.per.s <- 10^(-4.36 - 0.01*ZEN.dat$TPSA)*1000000 #conversion Caco-2 Papp to human Peff
Peff.human.is <- (10^(0.4926*log10(Papp.cm.power.minus.six.per.s)- 0.1454)) ##human Peff(cm^-4/s)
Peff.rat.is = Peff.human.is/11.04 #rat Peff(cm^-4/s)
Peff.is <- Peff.rat.is/10000*3600 #conversion to cm/h
ka.is <- Peff.is*2/ZEN.dat$R #/h 

## Based on Videmann Caco2 experiments
VmaxAbs <- 8056 / 1000000 #umol/h, scaling from pmol/h
KmAbs <- 579 #pmol/h
VmaxSec <- 5708 / 1000000 #umol/h, scaling from pmol/h
KmSec <- 220 #uM

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


# -----------------------------------------------------------------------------
# PHYSIOLOGY ####
# -----------------------------------------------------------------------------

## VOLUMES AND FLOW RATES ####
## ------------------------------------------------------

BW <- 0.250 #(kg)
QC <- 80 * 60 / 1000 #(L/h)

Rat.input <- read.csv(here("Input", "Rat.input.csv"), header=1) #tissue composition mainly as R&R composition, unified by Utsey and optimised as per Musther
Rat.dat <- Rat.input %>% filter(!tissue %in% c("RBCs", "Plasma")) #df for all tissues except for adipose and RBCs

# Flow rates
## Extracting fractional flow rate (constant)
names <- substr(Rat.dat$tissue, 1, 2)
F_QC.names <- paste0("FQ", names)
F_QC.dat <- cbind(F_QC.names, select(Rat.dat, F_QC))

## Calculating tissue blood flow rates (L/h)
QC.dat <- F_QC.dat %>%
  # mutate(F_QC = F_QC * QC) %>% #QC = QCconstant*QC
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
  # mutate(F_Vo = F_Vo * BW) %>%
  rename(VTissue = F_Vo.names) %>%
  rename(VValue = F_Vo) %>%
  mutate(VTissue = sub("F", "", VTissue)) %>%
  pivot_wider(names_from = VTissue, values_from = VValue) %>%
  as.data.frame()


# -----------------------------------------------------------------------------
# EXTRACT AND SAVE ALL DATA ####
# -----------------------------------------------------------------------------

F_Vo.dat <- F_Vo.dat %>%
  pivot_wider(names_from = F_Vo.names, values_from = F_Vo)

MW = ZEN.dat$MW
logP = ZEN.dat$logP
MW.glu = ZEN.dat$MW.glu
logP.glu = ZEN.dat$logP.glu
BW = BW
QC = QC
QIn = QC.dat$QGu #Intestinal blood flow
#QPo = QC.dat$QGu + QC.dat$QSp #Portal blood flow = Σ small intestine, colon, pancreas, stomach, spleen
QLi = QC.dat$QLi #Liver blood flow (this includes portal! to be changed later! to add Q hepatic artery only)
QHa = QC.dat$QLi - QC.dat$QGu + QC.dat$QSp #Hepatic artery blood flow = total liver blood flow - portal vein blood flow (Σ small intestine, colon, pancreas, stomach, spleen)
QKi = QC.dat$QKi #Kidney blood flow
QAd = QC.dat$QAd #Adipose blood flow
QSl = QC.dat$QBo + QC.dat$QBr + QC.dat$QMu + QC.dat$QSk #Slowly perfused blood flow
QRa = QC - (QLi + QKi + QAd + QSl) #Rapidly perfused blood flow
MassBalanceQC = QC - (QLi + QKi + QAd + QSl + QRa)  #check blood flow mass balance; should be 0
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
PRa = (F_Vo.dat$FVHe*PCB.dat$PHe + F_Vo.dat$FVLu*PCB.dat$PLu + F_Vo.dat$FVSp*PCB.dat$PSp) /
  (F_Vo.dat$FVHe + F_Vo.dat$FVLu + F_Vo.dat$FVSp)
fup = fup #free fraction plasma
fub = fub #free fraction blood
PIn.glu = PCB.dat.glu$PGu #Intestine:Blood partition coefficient
PLi.glu = PCB.dat.glu$PLi #Liver:Blood partition coefficient
PKi.glu = PCB.dat.glu$PKi #Kidney:Blood partition coefficient
PAd.glu = PCB.dat.glu$PAd #Adipose:Blood partition coefficient
PSl.glu = (F_Vo.dat$FVBo*PCB.dat.glu$PBo + F_Vo.dat$FVBr*PCB.dat.glu$PBr +
             F_Vo.dat$FVMu*PCB.dat.glu$PMu + F_Vo.dat$FVSk*PCB.dat.glu$PSk) /
  (F_Vo.dat$FVBo + F_Vo.dat$FVBr + F_Vo.dat$FVMu + F_Vo.dat$FVSk) #slowly perfused = partition coefficient of each of these organs * fractional volume of the respective organ / sum of the fractional volume of all these organs
PRa.glu = (F_Vo.dat$FVHe*PCB.dat.glu$PHe + F_Vo.dat$FVLu*PCB.dat.glu$PLu +
             F_Vo.dat$FVSp*PCB.dat.glu$PSp) /
  (F_Vo.dat$FVHe + F_Vo.dat$FVLu + F_Vo.dat$FVSp)
fup.glu = fup.glu #free fraction plasma
fub.glu = fup.glu #free fraction blood
kTIn = kTIn #intestinal transit time (1/h)
kTge = kTge #gastric emptying rate (1/h)
ka = ka #intestinal absorption constant (1/h)
ka.is = ka.is #intestinal absorption constant, in silico method (1/h)
ka.sl = 0.05 * 60 #slow luminal absorption rate per Mukherjee (1/h)
ka.fa = 20 * 60 #fast luminal absorption rate per Mukherjee (1/h)
ka.fit = 0.1545 * 60 #rat fitted per Mukherjee (1/h)
kex = 0.0115 * 60 #rat fecal excretion rate per Mukherjee (1/h)
kfil = 0.0016 * 60 #rat urinary filtration rate per Mukherjee (1/h)
kbil = 0.0076 * 60 #rat biliary recirculation rate per Mukherjee (1/h)
ka.is = ka.is #intestinal absorption constant, in silico method (1/h)
ka.sl.glu = 0.05 * 60 #slow luminal absorption rate per Mukherjee (1/h)
ka.fa.glu = 20 * 60 #fast luminal absorption rate per Mukherjee (1/h)
ka.fit.glu = 0.0833 * 60 #rat fitted per Mukherjee (1/h)
kex.glu = 0.053 * 60 #rat fecal excretion rate per Mukherjee (1/h)
kfil.glu = 0.0013 * 60 #rat urinary filtration rate per Mukherjee (1/h)
kbil.glu = 0.2596 * 60 #rat biliary recirculation rate per Mukherjee (1/h)
VmaxAbs <- VmaxAbs #umol/h
KmAbs <- KmAbs #pmol/h
VmaxSec <- VmaxSec #umol/h
KmSec <- KmSec #uM
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


#----Saving-------------------------------------------------------------------#

INPUT.parameters <- cbind(
  MW, 
  logP,
  MW.glu,
  logP.glu,
  BW,
  QC,
  QIn,
  QLi,
  QHa,
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
  PIn.glu, 
  PLi.glu, 
  PKi.glu,
  PAd.glu,
  PSl.glu, 
  PRa.glu,
  fup.glu, 
  fub.glu, 
  kTIn,
  kTge,
  ka, 
  ka.is,
  ka.sl,
  ka.fa,
  ka.fit,
  kex,
  kfil,
  kbil,
  ka.sl.glu,
  ka.fa.glu,
  ka.fit.glu,
  kex.glu,
  kfil.glu,
  kbil.glu,
  VmaxAbs,
  KmAbs,
  VmaxSec,
  KmSec,
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

write.csv(INPUT.parameters, file = here("Input", "CALCULATED.INPUT.csv"), row.names=FALSE)