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
Rat.dat1 <- Rat.input %>% filter(!tissue %in% c("Plasma","Adipose","RBCs")) #df for all tissues except for adipose, RBCs and plasma
Rat.dat2 <- Rat.input %>% filter(!tissue %in% c("Plasma","RBCs")) #df for all tissues expluding Plasma and RBCs
fup <- 0.156 #Mukherjee


## ZEARALENONE ##

# Tissue partitioning (Poulin and Theil http://jpharmsci.org/article/S0022-3549(16)30889-9/fulltext), as coded by kmutsey: https://github.com/metrumresearchgroup/PBPK_PC/blob/master/script/CalcKp_P%26T.R
# Assumption that Zearalenone is a monoprotic acid; in reality it's a diprotic acid with pKa2=10.38. Insignificant difference in tissue:plasma partition coefficients was observed when diprotic acid equations were tested. Check.As.Diprotic.RnR.R file.
# Zearalenone is not ionized at physiological pH. RnR's equation for neutrals was overestimating the partition in all tissues by 7%. Check As.Neutral.RnR.R
# Poulin and Theil is the best QSPR for gave the best predictions for neutral compounds 

n <- length(Rat.input$tissue)

#Plasma composition
Vwp <- Rat.input$f_water[Rat.input$tissue == "Plasma"] #fraction of water in plasma
Vnlp <- Rat.input$f_n_l[Rat.input$tissue == "Plasma"] #fraction of natural lipids in plasma
Vphp <- Rat.input$f_pl[Rat.input$tissue == "Plasma"] #fraction of phospholipids in plasma

# Partition coefficients of all tissues except adipose
Vwt <- Rat.dat2$f_water[Rat.dat2$tissue != "Adipose"] #fraction of water in all tissues except adipose
Vnlt <- Rat.dat2$f_n_l[Rat.dat2$tissue != "Adipose"] #fraction of natural lipids all tissues except adipose
Vpht <- Rat.dat2$f_pl[Rat.dat2$tissue != "Adipose"] #fraction of phospholipids in all tissues except adipose


# Partition coefficient of adipose
Vwad <- Rat.dat2$f_water[Rat.dat2$tissue == "Adipose"] #fraction of water in adipose
Vnlad <- Rat.dat2$f_n_l[Rat.dat2$tissue == "Adipose"] #fraction of natural lipids in adipose
Vphad <- Rat.dat2$f_pl[Rat.dat2$tissue == "Adipose"] #fraction of phospholipids in adipose

pH <- Rat.dat2$pH[Rat.dat2$tissue == "Adipose"] #pH in adipose
logD <- 1.115*ZEN.dat$logP-1.35 #logD is the olive oil:buffer(water) partition coefficient of nonionized species

logD_star <- logD-log10(1+10^(pH-ZEN.dat$pKa)) #as monoprotic acid
D_star <- 10^logD_star   
Kpad <- ((D_star*(Vnlad+0.3*Vphad)+(1*(Vwad+0.7*Vphad)))/(D_star*(Vnlp+0.3*Vphp)+(1*(Vwp+0.7*Vphp)))) * fup

P <- 10^ZEN.dat$logP
fut <- 1/(1+((1-fup)/fup)*0.5) #fraction unbound in the tissues
Kpt <- ((P*(Vnlt+0.3*Vpht) + #partition in the natural lipids and phospholipids of the tissue
           (1*(Vwt+0.7*Vpht))) / #partition in the water and the phospholipids of the tissue
          (P*(Vnlp+0.3*Vphp) + #partition in the natural lipids and phospholipids of plasma
             (1*(Vwp+0.7*Vphp)))) * #partition in the water and phospholipids of plasma
  (fup/fut) #related to the fraction unbound in the plasma compared to the tissue

Kp <- c(Kpad, Kpt)

PC.names <- paste0("P", names)
PC.dat <- data.frame(PC.names, kp = as.numeric(Kp))

## Calculate Tissue:Blood partition coefficient 
PCB.dat <- PC.dat %>%
  mutate(kp = kp * ZEN.dat$BP) %>% #Tissue:Blood Partition coefficients = BloodPlasma ratio (BP) * Tissue:Plasma Partition coefficient (PTi)
  rename(PValue = kp) %>%
  rename(PTissue = PC.names) %>%
  pivot_wider(names_from = PTissue, values_from = PValue) %>%
  as.data.frame()

# Calculate free fraction in blood
fub = fup * ZEN.dat$BP 



## ZEARALENONE GLUCURONIDE ##

# Calculating fraction unbound in plasma #Reference: Lobell M, Sivarajah V. In silico prediction of aqueous solubility, human plasma protein binding and volume of distribution of compounds from calculated pKa and AlogP98 values. Mol Divers. 2003;7(1):69-87. doi: 10.1023/b:modi.0000006562.93049.36. PMID: 14768905.
fup.glu <- 1 / (1 + 10^(0.4485 * ZEN.dat$logP.glu - 0.4782))

## Calculating ionisation status
## Zearalenone glucuronide is an acid and hydrophilic compound ; therefore the best QSPR for it is Schmitt. One could also consider modeling Zearalenone glucuronide kinetics are only driven by active transport and not diffusion.
K_n_pl <- 10^ZEN.dat$logP.glu #neutral phospholipids:water partition coefficient (logP here is logMA of Shmitt)
K_protein <- ((0.81 + 0.11 * K_n_pl) / 24.92) * 5 #protein:water partition; Schmitt, Walter (2008)
pH <- Rat.dat$pH #tissue pH
alpha <- 1e-3 #ratio between ditribution coefficient at given pH (D) and that in neutral form (D0)

W <- 10^(Rat.dat$pH - ZEN.dat$pKa.glu) #Calculates ZEN ionisation status depending on tissue pH
K_n_l <- K_n_pl * (((1 - alpha) / (1 + W)) + alpha) #Calculates neutral lipids:neutral phospholipids partition coefficient
K_a_pl <- K_n_pl * ((1 / (1 + W)) + 0.05 * (1 - (1 / (1 + W)))) #Calculates acid phospholipids:neutral lipids partition coefficient based on ionisation status

## Calculating Tissue:Plasma partition coefficients
Kp.glu <- (Rat.dat$f_water +
         (K_n_l * Rat.dat$f_n_l) +
         (K_n_pl * Rat.dat$f_n_pl) +
         (K_a_pl * Rat.dat$f_a_pl) +
         (K_protein * Rat.dat$f_proteins)) * fup.glu

PC.names.glu <- paste0("P", names, ".glu")
PC.dat.glu <- data.frame(PC.names.glu, kp = as.numeric(Kp.glu))


## Calculate Tissue:Blood partition coefficient 
PCB.dat.glu <- PC.dat.glu %>%
  #mutate(kp = kp * ZEN.dat$BP) %>% #assumption that BP is equal to 1
  rename(PValue = kp) %>%
  rename(PTissue = PC.names.glu) %>%
  pivot_wider(names_from = PTissue, values_from = PValue) %>%
  as.data.frame()


# Calculate free fraction in blood
#fub = fup * ZEN.dat$BP #assuming is equal as assuming BP equals 1 



#----Toxicokinetic parameters------------------------------------------------#

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




#----Extracting all data-----------------------------------------------------#

F_Vo.dat <- F_Vo.dat %>%
  pivot_wider(names_from = F_Vo.names, values_from = F_Vo)

MW = ZEN.dat$MW
logP = ZEN.dat$logP
MW.glu = ZEN.dat$MW.glu
logP.glu = ZEN.dat$logP.glu
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

write.csv(INPUT.parameters, file = here("InputData/", "CALCULATED.INPUT.csv"), row.names=FALSE)
