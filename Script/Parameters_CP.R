# --------------------------------------------------------------------------- #
# Script for calculating parameters  
# exact copy of the calculations of the PBK Model of Mendez-Catala et al. 2021 DOI: 10.1002/mnfr.202100443
# By: Chrysanthi Pachoulide
# Date: 23-06-2025
# --------------------------------------------------------------------------- #

rm(list = ls()) 

# Packages
library(here)
library(tidyverse)
library(showtext)
font_add(family = "Garamond",  regular = "GARA.TTF")
showtext_auto()

# Set output storage directory
OUTPUT = here("Output",  format(Sys.Date(),  "%Y-%m-%d"),  format(Sys.time(),  "%H-%M-%S"))
dir.create(OUTPUT,  recursive = TRUE)


# Parameters for Rat####

Rat.input <- read.csv(here("Input", "Rat.input.csv"), header=1) # tissue composition mainly as R&R composition, unified by Utsey and optimised as per Musther
Rat.dat <- Rat.input %>% filter(!tissue %in% c("RBCs", "Plasma")) # df for all tissues except for adipose and RBCs
names <- substr(Rat.dat$tissue, 1, 2)
F_Vo.names <- paste0("FV", names)
F_Vo.dat <- cbind(F_Vo.names, select(Rat.dat, F_Vo)) %>%
  pivot_wider(names_from = F_Vo.names, values_from = F_Vo)


  ## Physiological constants ####
  # Fractional Tissue volumes (reference: Punt et al. (2020))
  BW = 0.250  ##{Kg} #  body weight rat

  FV_ST = 0.012  # fraction volume of small intestine
  
  FV_SI = 0.026 # fractional volume of small intestinal tissue
  FV_Duo = 0.060
  FV_Jej = 0.769
  FV_Ile = 0.171
  
  FV_SIL = 0.014 # fractional volume of small intestinal lumen
  FV_duol = 0.183
  FV_jejl = 0.334 # calculated from Musther
  FV_ilel = 0.483 # calculated from Musther
  
  FV_COL = 0.012 # fractional volume of large intestinal lumen
  
  FV_L = 0.036  # fractional of liver tissue
  FV_R = 0.02  # 0.082-VSIc-VLc # fraction of rapidly perfused tissue
  FV_S = 0.667  # 0.737-VFc # fraction of slowly perfused tissue
  FV_F = 0.070  # fraction of fat tissue
  FV_B = 0.074  # fraction of blood
  
  # Surface area 
  SASI = 0.134*100 # surface area (*100 dm2 -> cm2)
  SALI = 1.57*100  # surface area of large intestine (*100 dm2 -> cm2) reference: Vdoviaková et al. (2016)
  
  R <- 0.18 #radius small intestine
  Rduo <- 0.275/2 #diameter / 2 = radius, cm (Musther et al)
  Rjej <- 0.4/2   #diameter / 2 = radius, cm (Musther et al)
  Rile <- 0.4/2   #diameter / 2 = radius, cm (Musther et al)

  Lduo <- 7    #length cm (Musther et al)
  Ljej <- 90   #cm (Musther et al)
  Lile <- 20   #diameter / 2 = radius, cm (Musther et al)
  
  Amp <- 120 # microvilli amplification factor, ref: https://doi.org/10.1016/S0278-6915(00)00136-8; villi and microvilli together amplify the small intestinal surface area by 60–120 times https://doi.org/10.3109/00365521.2014.898326
  Ampduo <- 1*Amp
  Ampjej <- 1*Amp #0.7
  Ampile <- 1*Amp #0.5
  
  SAduo = 2*pi*Rduo*Lduo*Ampduo #cm2
  SAjej = 2*pi*Rjej*Ljej*Ampjej #cm2
  SAile = 2*pi*Rile*Lile*Ampile #cm2
  
  # Blood flow rates (reference: Punt et al. (2020))
  QC = 4.8  ##{L/h} # cardiac output
  QSIc = 0.101  # fraction of blood flow to small intestine
  QDuoc = 0.079
  QJejc = 0.680
  QIlec = 0.241
  QLc = 0.141  # 0.242 - QSIc # fraction of blood flow to liver 101 
  QRc = 0.51  # 0.76 - QSIc - QLc # fraction of blood flow to rapidly perfused tissue
  QSc = 0.17  # 0.24 - QFc # fraction of blood flow to slowly perfused tissue
  QFc = 0.07  # fraction of blood flow to fat
  
  # Intestinal fluid flow rates; Wistar (reference: Musther et al SI; Punt et al. (2020))
  Tge <- 0.25 # stomach emptying time (h)
  Tsi <- 1.63 # small intestine transit time (h)
  Tsiduo <- 0.06*Tsi # duodenum transit time (h)
  Tsijej <- 0.192*Tsi # jejunum transit time (h)
  Tsiile <- 0.086*Tsi # ileum transit time (h)
  
  gaTR = 1/Tge # (/h)
  duoTR <- 1/Tsiduo # (/h)
  jejTR <- 1/Tsijej # (/h)
  ileTR <- 1/Tsiile # (/h)
  lilTR <- 1/19.7 #1/(0.8*19.7) # large intestinal transit time (/h)
  
  # Glomerular filtration rate: Punt et al. (2020)
  GFRc = 0.312 #L/h/kg 0.08 #(/h)
  
  ## Partition coefficients / fraction unbound / fraction unionised ####
  # Tissue partitioning (RodgaTRs and Rowland as coded by kmutsey: https://github.com/metrumresearchgroup/PBPK_PC/blob/master/script/CalcKp_P%26T.R
  # Assumption that Zearalenone and a-Zearalenol are diprotic acids and Zearalenone glucuronide a monoprotic acid
  #microsomes or S9 (Kilford et al., 2008 based on Halifax and Houston 2006)
  Rat.input <- read.csv(here("Input", "Rat.input.csv"))
  dat_all <- Rat.input %>% filter(!tissue %in% c("RBCs", "Adipose", "Plasma"))  #df for all tissues except for adipose, RBCs, and plasma
  dat_ad <- Rat.input %>% filter(tissue == "Adipose")  #df for adipose
  dat_rbc <- Rat.input %>% filter(tissue == "RBCs") #df for RBCs
  dat_plas <- Rat.input %>% filter(tissue == "Plasma") #df for aplasma
  
  pH_IW <- 7       #pH of intracellular tissue water
  pH_P <- 7.4      #pH of plasma
  pH_RBC <- 7.22    #pH of blood cells
  
  ### Calculate Kp for ZEN ####
  
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
  
  PSI_ZEN = Kp_ZEN$Kpgu #Intestine:Blood partition coefficient
  PL_ZEN = Kp_ZEN$Kpli #Liver:Blood partition coefficient
  PF_ZEN = Kp_ZEN$Kpad #Adipose:Blood partition coefficient
  PS_ZEN = (F_Vo.dat$FVBo*Kp_ZEN$Kpbo + F_Vo.dat$FVBr*Kp_ZEN$Kpbr +
              F_Vo.dat$FVMu*Kp_ZEN$Kpmu + F_Vo.dat$FVSk*Kp_ZEN$Kpsk) /
    (F_Vo.dat$FVBo + F_Vo.dat$FVBr + F_Vo.dat$FVMu + F_Vo.dat$FVSk) #slowly perfused = partition coefficient of each of these organs * fractional volume of the respective organ / sum of the fractional volume of all these organs
  PR_ZEN = (F_Vo.dat$FVHe*Kp_ZEN$Kphe + F_Vo.dat$FVKi*Kp_ZEN$Kpki + F_Vo.dat$FVLu*Kp_ZEN$Kplu + F_Vo.dat$FVSp*Kp_ZEN$Kpsp) /
    (F_Vo.dat$FVHe + F_Vo.dat$FVKi + F_Vo.dat$FVLu + F_Vo.dat$FVSp)
  
  fut_ZEN <- 1/(1+((1-fup_ZEN)/fup_ZEN)*0.5) #fraction unbound in the tissues
  
  pH <- 7.4
  logD_ZEN <- logP_ZEN-log10(1+10^(pH-pKa1_ZEN)+10^(pH-pKa1_ZEN+pH-pKa2_ZEN))
  
  frac.unionized_ZEN <- 1/(1+10^(pH-pKa1_ZEN)+10^(pH-pKa1_ZEN+pH-pKa2_ZEN))
  f.frac.ionized_ZEN <- 1 - 1/(1+10^(pH-pKa1_ZEN)+10^(pH-pKa1_ZEN+pH-pKa2_ZEN))
  
  ### Calculate Kp for aZEL ####
  
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
  
  PSI_aZEL = Kp_aZEL$Kpgu #Intestine:Blood partition coefficient
  PL_aZEL = Kp_aZEL$Kpli #Liver:Blood partition coefficient
  PF_aZEL = Kp_aZEL$Kpad #Adipose:Blood partition coefficient
  PS_aZEL = (F_Vo.dat$FVBo*Kp_aZEL$Kpbo + F_Vo.dat$FVBr*Kp_aZEL$Kpbr +
              F_Vo.dat$FVMu*Kp_aZEL$Kpmu + F_Vo.dat$FVSk*Kp_aZEL$Kpsk) /
    (F_Vo.dat$FVBo + F_Vo.dat$FVBr + F_Vo.dat$FVMu + F_Vo.dat$FVSk) #slowly perfused = partition coefficient of each of these organs * fractional volume of the respective organ / sum of the fractional volume of all these organs
  PR_aZEL = (F_Vo.dat$FVHe*Kp_aZEL$Kphe + F_Vo.dat$FVKi*Kp_aZEL$Kpki + F_Vo.dat$FVLu*Kp_aZEL$Kplu + F_Vo.dat$FVSp*Kp_aZEL$Kpsp) /
    (F_Vo.dat$FVHe + F_Vo.dat$FVKi + F_Vo.dat$FVLu + F_Vo.dat$FVSp)
  
  fut_aZEL <- 1/(1+((1-fup_aZEL)/fup_aZEL)*0.5) #fraction unbound in the tissues
  
  logD_aZEL <- logP_ZEN-log10(1+10^(pH-pKa1_aZEL)+10^(pH-pKa1_ZEN+pH-pKa2_aZEL))
  frac.unionized_aZEL <- 1/(1+10^(pH-pKa1_aZEL)+10^(pH-pKa1_aZEL+pH-pKa2_aZEL))
  f.frac.ionized_aZEL <- 1 - 1/(1+10^(pH-pKa1_aZEL)+10^(pH-pKa1_aZEL+pH-pKa2_aZEL))

  
  ### Calculate Kp for ZENGLUC ####
  
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
  
  PSI_ZENGLU = Kp_ZENGLU$Kpgu #Intestine:Blood partition coefficient
  PL_ZENGLU = Kp_ZENGLU$Kpli #Liver:Blood partition coefficient
  PF_ZENGLU = Kp_ZENGLU$Kpad #Adipose:Blood partition coefficient
  PS_ZENGLU = (F_Vo.dat$FVBo*Kp_ZENGLU$Kpbo + F_Vo.dat$FVBr*Kp_ZENGLU$Kpbr +
              F_Vo.dat$FVMu*Kp_ZENGLU$Kpmu + F_Vo.dat$FVSk*Kp_ZENGLU$Kpsk) /
    (F_Vo.dat$FVBo + F_Vo.dat$FVBr + F_Vo.dat$FVMu + F_Vo.dat$FVSk) #slowly perfused = partition coefficient of each of these organs * fractional volume of the respective organ / sum of the fractional volume of all these organs
  PR_ZENGLU = (F_Vo.dat$FVHe*Kp_ZENGLU$Kphe + F_Vo.dat$FVKi*Kp_ZENGLU$Kpki + F_Vo.dat$FVLu*Kp_ZENGLU$Kplu + F_Vo.dat$FVSp*Kp_ZENGLU$Kpsp) /
    (F_Vo.dat$FVHe + F_Vo.dat$FVKi + F_Vo.dat$FVLu + F_Vo.dat$FVSp)
  
  fut_ZENGLU <- 1/(1+((1-fup_ZENGLU)/fup_ZENGLU)*0.5) #fraction unbound in the tissues
  
  logD_ZENGLU <- logP_ZENGLU - log10(1+10^(pH-pKa_ZENGLU))
  
  frac.unionized_ZENGLU <- 1/(1+10^(pH-pKa_ZENGLU))
  f.frac.ionized_ZENGLU <- 1-1/(1+10^(pH-pKa_ZENGLU))

  # # Calculated based on R&R
  # # ZEN 
  # PSI_ZEN = 4.11 #2.64  # intestine/blood partition coefficient
  # PL_ZEN = 1.74  #2.64  # liver/blood partition coefficient
  # PR_ZEN = 1.60  #2.64  # rapidly perfused tissue/blood partition coefficient
  # PS_ZEN = 2.19 #0.76  # slowly perfused tissue/blood partition coefficient
  # PF_ZEN = 12.60 #106.93  # fat/blood partition coefficient
  # 
  # # ZEN-gluc
  # PSI_ZENGLU = 0.97
  # PL_ZENGLU = 1.66
  # PR_ZENGLU = 0.95
  # PS_ZENGLU = 1.30
  # PF_ZENGLU = 0.23
  # 
  # # α-ZEL in sub-model; need to change them
  # PIaZEL = 2.38  # intestine/blood partition coefficient
  # PLaZEL = 2.38  # liver/blood partition coefficient
  # PRaZEL = 2.38  # rapidly perfused tissue/blood partition coefficient
  # PSaZEL = 0.71  # slowly perfused tissue/blood partition coefficient
  # PFaZEL = 92.79  # fat/blood partition coefficient
  # 
  # # fraction unbound ####
  # fup_ZEN = 0.156 #0.270
  # fup_aZEL = 0.156
  # fup_ZENGLU = 0.20
  
  ### absorption/transfer rates ####
  PappCaco2_ZEN = 10.4 * 1e-6 * 3600 # (cm/s->h) reference: Pfeiffer et al. (2011) DOI 10.1002/mnfr.201000381
  PappCaco2_aZEL = 5.4 * 1e-6 * 3600 # (cm/s->h) reference: Pfeiffer et al. (2011) DOI 10.1002/mnfr.201000381

  Peff.cm.per.hr <- PappCaco2_ZEN #cm/h
  
  ## Kinetic parameters ####
  
  ### metabolism in small intestine tissue ####
  S9SI =  57 #37.1  # small intestinal S9 protein yield (mg S9 protein/gram intestine) reference: Peters et al. (2016); sum of microsomal 3.9 mg/g and cytosolic 18mg/g https://sci-hub.box/https://doi.org/10.1007/s13318-020-00666-w
  
  conc.HSD.S9I = 1.25 # mg protein/mL Malekinejad et al. (2006) doi:10.1016/j.tvjl.2005.03.004
  VmaxSI_aZELc = 32*1e-6 * 60      # assuming the same as liver; pmol->umol/min->h/mg protein reference: Malekinejad et al. (2006) doi:10.1016/j.tvjl.2005.03.004 (post mitochondrial fraction)
  KmSI_aZEL = 592                   # assuming the same as liver; uM; reference: Malekinejad et al. (2006) doi:10.1016/j.tvjl.2005.03.004
  VmaxSI_bZELc = 72*1e-6 * 60     # assuming the same as liver; pmol->umol/min->h/mg protein; reference: Malekinejad et al. (2006) doi:10.1016/j.tvjl.2005.03.004 (post mitochondrial fraction)
  KmSI_bZEL = 21                   # assuming the same as liver; uM; reference: Malekinejad et al. (2006) doi:10.1016/j.tvjl.2005.03.004
  
  conc.UGT.S9I = 0.04 # mg protein/ml Diana DOI: 10.1002/mnfr.202100443
  VmaxSI_ZENGLUc =  8.617*1e-3 * 60  # nmol->umol/min->h/mg S9; DOI: 10.1002/mnfr.202100443
  KmSI_ZENGLU = 6.198  # uM; DOI: 10.1002/mnfr.202100443
  
  VmaxSI_aZELGLUc = 6.963*1e-3 * 60 # nmol->umol/min->h/mg s9; DOI: 10.1002/mnfr.202100443
  KmSI_aZELGLU = 7.426              # uM; DOI: 10.1002/mnfr.202100443
  
  fuincS9I.HSD_ZEN <- (f.frac.ionized_ZEN+frac.unionized_ZEN)*(1/(1+conc.HSD.S9I*10^(0.072*logD_ZEN^2+0.067*logD_ZEN-1.126)))
  fuincS9I.UGT_ZEN <- (f.frac.ionized_ZEN+frac.unionized_ZEN)*(1/(1+conc.UGT.S9I*10^(0.072*logD_ZEN^2+0.067*logD_ZEN-1.126)))
  fuincS9I_aZEL <- (f.frac.ionized_aZEL+frac.unionized_aZEL)*(1/(1+conc.UGT.S9I*10^(0.072*logD_aZEL^2+0.067*logD_aZEL-1.126)))
  
  
  ### metabolism in large intestine lumen (microbiota compartment) ####
  VMB = 0.0164                    # fraction of feces of BW reference: Hoskins and Zamcheck (1968)
  VmaxLI_aZELc = 0.23*1e-6 * 60   # pmol->umol/min->h/mg feces reference: Mendez-Catala et al. (2020) https://doi.org/10.1016/j.toxrep.2020.07.010
  KmLI_aZEL = 66                  # uM; reference: Mendez-Catala et al. (2020)
  VmaxLI_bZELc = 0.14*1e-6 * 60   # pmol->umol/min->h/mg feces reference: Mendez-Catala et al. (2020)
  KmLI_bZEL = 80                  # uM; reference: Mendez-Catala et al. (2020)
  
  ### metabolism in liver ####
  S9L = 143 #125                      # liver S9 protein yield, (mg S9 protein/gram liver) references: sum of microsomal 38mg/g liver and cytosolic 87mg/g liver content, references https://doi.org/10.2174/138920007779315053 and 10.1016/j.taap.2011.03.020, https://sci-hub.box/https://doi.org/10.1007/s13318-020-00666-w
  conc.HSD.S9L = 1.25 # mg protein/mL Malekinejad et al. (2006) doi:10.1016/j.tvjl.2005.03.004
  
  VmaxL_aZELc = 32*1e-6 * 60      # pmol->umol/min->h/mg protein reference: Malekinejad et al. (2006) doi:10.1016/j.tvjl.2005.03.004 (post mitochondrial fraction)
  KmL_aZEL = 592                  # uM; reference: Malekinejad et al. (2006) doi:10.1016/j.tvjl.2005.03.004
  VmaxL_bZELc = 72*1e-6 * 60      # pmol->umol/min->h/mg protein; reference: Malekinejad et al. (2006) doi:10.1016/j.tvjl.2005.03.004 (post mitochondrial fraction)
  KmL_bZEL = 21                   # uM; reference: Malekinejad et al. (2006) doi:10.1016/j.tvjl.2005.03.004
  
  conc.UGT.S9L = 0.2 # mg protein/ml Diana DOI: 10.1002/mnfr.202100443
  VmaxL_ZENGLUc = 7.029*1e-3 * 60  # umol/h/mg S9; nmol->umol/min->h/mg S9; DOI: 10.1002/mnfr.202100443
  KmL_ZENGLU = 6.753               # uM; DOI: 10.1002/mnfr.202100443
  
  VmaxL_aZELGLUc = 6.963*1e-3 * 60 # umol/h/mg S9; nmol->umol/min->h/mg s9; DOI: 10.1002/mnfr.202100443
  KmL_aZELGLU = 7.426              # uM; DOI: 10.1002/mnfr.202100443

  fuincS9L.HSD_ZEN <- (f.frac.ionized_ZEN+frac.unionized_ZEN)*(1/(1+conc.HSD.S9L*10^(0.072*logD_ZEN^2+0.067*logD_ZEN-1.126)))
  fuincS9L.UGT_ZEN <- (f.frac.ionized_ZEN+frac.unionized_ZEN)*(1/(1+conc.UGT.S9L*10^(0.072*logD_ZEN^2+0.067*logD_ZEN-1.126)))
  fuincS9L_aZEL <- (f.frac.ionized_aZEL+frac.unionized_aZEL)*(1/(1+conc.UGT.S9L*10^(0.072*logD_aZEL^2+0.067*logD_aZEL-1.126)))
  
  ### active uptake and secretion in small intestine tissue ####
  VmaxOut_aZELc = 208*1e-6       # pmol->umol/cm2/h doi:10.1016/j.fct.2008.07.011
  KmOut_aZEL = 46.9              # uM doi:10.1016/j.fct.2008.07.011
  
  VmaxIn_aZELc = 102.8*1e-6      # pmol->umol/cm2/h doi:10.1016/j.fct.2008.07.011
  KmIn_aZEL = 14.3             # uM doi:10.1016/j.fct.2008.07.011
  
  VmaxOut_bZELc = 260.5*1e-6       # pmol->umol/cm2/h doi:10.1016/j.fct.2008.07.011
  KmOut_bZEL = 42.9              # uM doi:10.1016/j.fct.2008.07.011
  
  VmaxIn_bZELc = 51.3*1e-6      # pmol->umol/cm2/h doi:10.1016/j.fct.2008.07.011
  KmIn_bZEL = 10.6             # uM doi:10.1016/j.fct.2008.07.011
  
  VmaxOut_ZENGLUc = 40.1*1e-6       # pmol->umol/cm2/h doi:10.1016/j.fct.2008.07.011
  KmOut_ZENGLU = 23.6              # uM doi:10.1016/j.fct.2008.07.011
  
  VmaxIn_ZENGLUc = 33.3*1e-6      # pmol->umol/cm2/h doi:10.1016/j.fct.2008.07.011
  KmIn_ZENGLU = 40.0             # uM doi:10.1016/j.fct.2008.07.011
  
  VmaxOut_aZELGLUc = 53.5*1e-6       # pmol->umol/cm2/h doi:10.1016/j.fct.2008.07.011
  KmOut_aZELGLU = 10.8              # uM doi:10.1016/j.fct.2008.07.011
  
  VmaxIn_aZELGLUc = 28.7*1e-6      # pmol->umol/cm2/h doi:10.1016/j.fct.2008.07.011
  KmIn_aZELGLU = 7.6             # uM doi:10.1016/j.fct.2008.07.011
  
  ### Biliary excretion ####
  CL_Bile.fasted <- 0.00075
  CL_Bile.fed <- 0.1 #1
  
  Vmax_BSEPc <- 7                            # umol/min/mg BSEP, Average active transport of bile acids, assuming that the maximum velocity of PFOA transport by BSEP corresponds to that of bile acids ; de Bruijn et al. 2024 https://doi.org/10.14573/altex.2302011 
  Km_BSEPc <- 16.4                           # umol/L, uM, Average affinity constant of bile acids to BSEP, following the above assumption ; de Bruijn et al. 2024 https://doi.org/10.14573/altex.2302011
  
  # IVIVE 
  mgBSEP.HC <- 0.839 * 140000 * 1e-9         # mg BSEP/1e6 hepatocytes (0.839 pmole BSEP/1e6 hepatocytes (amound of BSEP per hepatocyte) * 140000 g/mole (MW BSEP) -> pg * 1e-9 -> mg) de Bruijn et al. 2024 https://doi.org/10.14573/altex.2302011]
  HC.GL <- 99                                 # 1e6 hepatocytes in the liver ; de Bruijn et al. 2024 https://doi.org/10.14573/altex.2302011
  SF_BSEP <- mgBSEP.HC * HC.GL * 1e3          # Scaling factor for BSEP mediated hepatic efflux for GCA and GCDC ; de Bruijn et al. 2024 https://doi.org/10.14573/altex.2302011
  
  f.in_out.MRP2 <- 0.42 # inside-out (active) fraction of transporter in in citro system vesicles
  f.in_out.MRP3 <- 0.58 # inside-out (active) fraction of transporter in in citro system vesicles
  pmol_mg.MRP2 <- 133.9 # protein abundance in vesicles, in pmol transporter/mg protein
  pmol_mg.MRP3 <- 239.4 # protein abundance in vesicles, in pmol transporter/mg protein
  pmolMP_mgP.LI.MRP2 <- 1.5 # pmol MRP2/mg protein in liver
  pmolMP_mgP.LI.MRP3 <- 0.5 # pmol MRP3/mg protein in liver
  pmolMP_mgP.IN.MRP2 <- 0.9 # pmol MRP2/mg protein in intestine
  pmolMP_mgP.IN.MRP3 <- 1.1 # pmol MRP3/mg protein in intestine
  pmolMP_mgP.KI.MRP2 <- 1.4 # pmol MRP2/mg protein in kidney
  pmolMP_mgP.KI.MRP3 <- 0.05 # pmol MRP3/mg protien in kidney
  mgP_mgLI <- 0.18 # mg protein per mg liver
  mgP_mgIN <- 0.133 # mg protein per mg intestine
  mgP_mgKI <- 0.17 # mg protein per mg kidney
  
  Vmax.MRP2 <- 4767e-6 # pmol -> umol/min/mg MRP2, 2935.5
  Km.MRP2 <- 109.6 # umol/L, 76.925
  Vmax.MRP3 <- 1520e-6 # pmol -> umol/min/mg , 839.675
  Km.MRP3 <- 6.0 # umol/L, 3.25
  
  CL.LI.MRP2.c <- (Vmax.MRP2 * 60 * 1e6 /Km.MRP2)/ f.in_out.MRP2 / pmol_mg.MRP2 / pmolMP_mgP.LI.MRP2 / mgP_mgLI # L/min -> h/mg -> kg liver
  CL.LI.MRP3.c <- (Vmax.MRP3 * 60 * 1e6 /Km.MRP3)/ f.in_out.MRP3 / pmol_mg.MRP3 / pmolMP_mgP.LI.MRP3 / mgP_mgLI # L/min -> h/mg -> kg liver
  
  CL.IN.MRP2.c <- (Vmax.MRP2 * 60 * 1e6 /Km.MRP2)/ f.in_out.MRP2 / pmol_mg.MRP2 / pmolMP_mgP.IN.MRP2 / mgP_mgIN # L/min -> h/mg -> kg intestine
  CL.IN.MRP3.c <- (Vmax.MRP3 * 60 * 1e6 /Km.MRP3)/ f.in_out.MRP3 / pmol_mg.MRP3 / pmolMP_mgP.IN.MRP3 / mgP_mgIN # L/min -> h/mg -> kg intestine
  
  CL.KI.MRP2.c <- (Vmax.MRP2 * 60 * 1e6 /Km.MRP2)/ f.in_out.MRP2 / pmol_mg.MRP2 / pmolMP_mgP.KI.MRP2 / mgP_mgKI # L/min -> h/mg -> kg kidney
  CL.KI.MRP3.c <- (Vmax.MRP3 * 60 * 1e6 /Km.MRP3)/ f.in_out.MRP3 / pmol_mg.MRP3 / pmolMP_mgP.KI.MRP3 / mgP_mgKI # L/min -> h/mg -> kg kidney
  
  ### Hydrolysis by GUS ####
  CL_Hydro <- 10 #L/h
  
  Parm.c.rat <- as.list(
    c(BW = BW,
      FV_ST = FV_ST,
      FV_SIL = FV_SIL,
      FV_duol = FV_duol, 
      FV_jejl = FV_jejl,
      FV_ilel = FV_ilel, 
      FV_COL = FV_COL,
      FV_SI = FV_SI,
      FV_Duo = FV_Duo,
      FV_Jej = FV_Jej,
      FV_Ile = FV_Ile,
      FV_L = FV_L,
      FV_R = FV_R,
      FV_S = FV_S,
      FV_F = FV_F,
      FV_B = FV_B,
      QC = QC,
      QSIc = QSIc,
      QDuoc = QDuoc, 
      QJejc = QJejc, 
      QIlec = QIlec, 
      QLc = QLc,
      QRc = QRc,
      QSc = QSc,
      QFc = QFc,
      GFRc = GFRc,
      SAduo = SAduo,
      SAjej = SAjej,
      SAile = SAile,
      PappCaco2_ZEN = PappCaco2_ZEN,
      S9SI = S9SI,
      VmaxSI_aZELc = VmaxSI_aZELc,
      VmaxSI_bZELc = VmaxSI_bZELc,
      VmaxSI_ZENGLUc = VmaxSI_ZENGLUc,
      VmaxSI_aZELGLUc = VmaxSI_aZELGLUc,
      VMB = VMB,
      VmaxLI_aZELc = VmaxLI_aZELc,
      VmaxLI_bZELc = VmaxLI_bZELc,
      S9L = S9L,
      VmaxL_aZELc = VmaxL_aZELc,
      VmaxL_bZELc = VmaxL_bZELc,
      VmaxL_ZENGLUc = VmaxL_ZENGLUc,
      VmaxL_aZELGLUc = VmaxL_aZELGLUc,
      CL_Bile.fed = CL_Bile.fed,
      CL_Bile.fasted = CL_Bile.fasted,
      Vmax_BSEPc = Vmax_BSEPc,
      SF_BSEP = SF_BSEP,
      Km_BSEPc = Km_BSEPc,
      CL_Hydro = CL_Hydro,
      BP_ZEN = BP_ZEN,
      BP_aZEL = BP_aZEL,
      BP_ZENGLU = BP_ZENGLU,
      gaTR = gaTR,
      duoTR = duoTR,
      jejTR = jejTR,
      ileTR = ileTR,
      lilTR = lilTR,
      KmLI_aZEL = KmLI_aZEL,
      KmLI_bZEL = KmLI_bZEL,
      KmSI_ZENGLU = KmSI_ZENGLU,
      fut_ZEN = fut_ZEN,
      KmSI_aZEL = KmSI_aZEL,
      KmSI_bZEL = KmSI_bZEL,
      fut_aZEL = fut_aZEL,
      KmSI_aZELGLU = KmSI_aZELGLU,
      KmL_aZEL = KmL_aZEL,
      KmL_bZEL = KmL_bZEL,
      KmL_ZENGLU = KmL_ZENGLU,
      fut_ZENGLU = fut_ZENGLU,
      KmL_aZELGLU = KmL_aZELGLU,
      fup_ZEN = fup_ZEN,
      fup_aZEL = fup_aZEL,
      fup_ZENGLU = fup_ZENGLU,
      PSI_ZEN = PSI_ZEN,
      PL_ZEN = PL_ZEN,
      PF_ZEN = PF_ZEN,
      PS_ZEN = PS_ZEN,
      PR_ZEN = PR_ZEN,
      PSI_aZEL = PSI_aZEL,
      PL_aZEL = PL_aZEL,
      PF_aZEL = PF_aZEL,
      PS_aZEL = PS_aZEL,
      PR_aZEL = PR_aZEL,
      PSI_ZENGLU = PSI_ZENGLU,
      PL_ZENGLU = PL_ZENGLU,
      PF_ZENGLU = PF_ZENGLU,
      PS_ZENGLU = PS_ZENGLU,
      PR_ZENGLU = PR_ZENGLU,
      CL.LI.MRP2.c = CL.LI.MRP2.c,
      CL.LI.MRP3.c = CL.LI.MRP3.c,
      CL.IN.MRP2.c = CL.IN.MRP2.c,
      CL.IN.MRP3.c = CL.IN.MRP3.c,
      CL.KI.MRP2.c = CL.KI.MRP2.c,
      CL.KI.MRP3.c = CL.KI.MRP3.c) 
  ) %>% as.data.frame()
  
  Parm.rat <- Parm.c.rat %>% 
    pivot_longer(names_to = "Parameter", 
                 values_to = "Value", 
                 cols = everything())
  
  write.csv(Parm.rat, here("Input","tbdParm.rat.csv"), row.names = FALSE)
  

# Parameters for human ####

  ## Physiological parameters ####
  # Tissue volumes, reference: Brown et al. (1997)
  BW <- 70 #{Kg} body weight rat
  VSIc <- 0.014 #0.009 # fraction of small intestine
  VLc <- 0.026 # fraction of liver tissue
  VRc <- 0.041 # 0.076-VSIc-VLc # fraction of rapidly perfused tissue
  VSc <- 0.596 # 0.81-VFc # fraction of slowly perfused tissue
  VFc <- 0.214 # fraction of fat tissue
  VBc <- 0.074 # fraction of blood
  VLIc <- 0.168 # fraction of large intestine
  
  
  # Blood flow rates, reference: Brown et al. (1997)
  QC <- 347.9 #{L/h} # cardiac output:15*BW^0.74
  QSIc <- 0.181 # fraction of blood flow to small intestine
  QLc <- 0.046 # 0.227 - QSIc # fraction of blood flow to liver
  QRc <- 0.473 # 0.7 - QSIc – QLc # fraction of blood flow to rapidly perfused tissue
  QSc <- 0.248 # 0.3 - QFc # fraction of blood flow to slowly perfused tissue
  QFc <- 0.052 # fraction of blood flow to fat
  
  ## Physicochemical parameters ####
  ### partition coefficients ####
  # calculated based on QPPR of DeJongh et al. (1997)
  
  # ZEN in main model
  PIZEN <- 16.95 #6.56 # intestine/blood partition coefficient
  PLZEN <- 12.24 #6.56 # liver/blood partition coefficient
  PRZEN <- 30.54 #6.56 # rapidly perfused tissue/blood partition coefficient
  PSZEN <- 21.70 #4.25 # slowly perfused tissue/blood partition coefficient
  PFZEN <- 218.06 #134.86 # fat/blood partition coefficient
  
  # α-ZEL in sub-model
  PIaZEL <- 6.11 # intestine/blood partition coefficient
  PLaZEL <- 6.11 # liver/blood partition coefficient
  PRaZEL <- 6.11 # rapidly perfused tissue/blood partition coefficient
  PSaZEL <- 3.99 # slowly perfused tissue/blood partition coefficient
  PFaZEL <- 131.62 # fat/blood partition coefficient
  
  ### absorption/transfer rates ####
  # Stomach emptying rate
  ksto <- 2.8 # stomach emptying rate (/h) reference: Reilly et al. (1990)
  
  # intestinal absorption and transfer rates of ZEN
  #Papp, Caco-2 ZEN <- 10.3 #( x10-5 cm/sec) fitted apparent in vitro, from Caco-2 reference: Pfeiffer et al.(2011)
  PappCaco2ZEN <- 4.9829 # Log Papp, Caco-2
  # Log (Papp,in vivo) <- 0.6836*Log(PappCaco2ZEN)-0.5579 reference: Sun et al. (2002)
  
  Vin <- 1.29 # volume for each compartment of intestines (L)
  SAin <- 10.3 # surface area (dm2)
  kin <- 2.19 # transfer rate to next compartment within the intestines (/hr)
  
  kin7 <- 0.464 # transfer rate of ZEN from small intestine to large intestine (/h) reference: Kimura and Higaki (2002)
  SAinb <- 47.12 # surface area of large intestine (dm2) reference: Vdoviaková et al. (2016)
  
  # α-ZEL
  # Papp, Caco-2 α-ZEL <- 5.4 (x10-6 cm/sec), apparent in vitro, from Caco-2 reference: Pfeiffer et al. (2011)
  PappCaco2_aZEL<--5.2676 # Log Papp, Caco-2
  #  Log (Papp,in vivo) <- 0.6836*Log(PappCaco2_aZEL)-0.5579 reference: Sun et al. (2002)
 
  # excretion (/h)
  KurZ <- 0.096 # urinary excretion of ZEN reference: Mukherjee et al. (2014)
  KfeZ <- 0.024 # fecal excretion of ZEN reference: Mukherjee et al. (2014)
  KuraZ <- 0.015 # urinary excretion of α-ZEL reference: Mukherjee et al. (2014)
  KfeaZ <- 0.024 # fecal excretion of α-ZEL reference: Mukherjee et al. (2014)
  KurG <- 59.4 # urinary excretion of α-ZEL reference: Mukherjee et al. (2014)
  
  ##  Kinetic parameters ####
  
  ### metabolism of small intestine tissue ####
  # scaling factors
  S9SI<- 37.10 # small intestinal S9 protein yield (mg S9 protein/gram intestine) S. A. Peters et al. DOI:10.1007/s40262-015-0351-6
  # metabolites ZENGLU, unscaled maximum rates of metabolism, (nmol/min/mg S9 protein)
  VmaxSIZENGLUc<- 0.491 # in vitro, S9 incubations from this study
  
  # metabolites ZENGLU, affinity constants, (μmol/L)
  KmSIZENGLU <- 1.174 # in vitro, S9 incubations from this study
  
  ### metabolism of large intestine lumen (microbiota compartment) ####
  # scaling factors
  VMB <- 0.0018 # fraction of faeces of BW, reference: Rose et al. (2015)
  #  metabolites α-ZEL and ß-ZEL, unscaled maximum rates of metabolism, (pmol/h/g faeces)
  VmaxLIaZELc<-0.054 # in vitro, anaerobic rat fecal incubations
  VmaxLIbZELc<- 0.0108 # in vitro, anaerobic rat fecal incubations
 
  #  metabolites α-ZEL and ß-ZEL, affinity constants, (μmol/L)
  KmLIaZEL <- 135 # in vitro, anaerobic rat fecal incubations
  KmLIbZEL <- 163 # in vitro, anaerobic rat fecal incubations
  
  ### metabolism of liver ####
  # scaling factors
  VLS9 <- 112.7  #  liver S9 protein yield, (mg S9 protein/gram liver) ref.  Moerenhout et al 2024. DOI:10.1021/acs.est.4c06534
  # Part 1: Zearalenone phase I metabolism α-ZEL and ß-ZEL formation
  # metabolites α-ZEL and ß-ZEL, unscaled maximum rates of metabolism, (pmol/min/mg S9 protein)
  VmaxLaZELc <- 358.7 # reference: Malekinejad et al. (2006)
  VmaxLbZELc <- 209.3 # reference: Malekinejad et al. (2006)
  
  #  metabolites α-ZEL and ß-ZEL, affinity constants, (μmol/L)
  KmLaZEL <- 9 # reference: Mendez-Catala et al. (2020)
  KmLbZEL <- 23 # reference: Mendez-Catala et al. (2020)
  
  #  Part 2: Zearalenone phase II metabolism-glucuronidation
  #  metabolites ZENGLU, unscaled maximum rates of metabolism, (nmol/min/mg S9 protein)
  VmaxLZENGLUc <- 2.97 # in vitro, S9 incubations from this study

  #  metabolites ZENGLU, affinity constants, (μmol/L)
  KmLZENGLU <- 2.043 # in vitro, S9 incubations from this study
  #  Part 3: α-ZEL phase II metabolism: glucuronidation
  #  metabolites α-ZELGLU, unscaled maximum rates of metabolism, (nmol/min/mg S9 protein)
  VmaxLaZELGLUc<- 2.983 # in vitro, S9 incubations from this study

  #  metabolites α-ZELGLU, affinity constants, (μmol/L)
  KmLaZELGLU <-2.415 # in vitro, S9 incubations from this study
  
  
  Parm.c.human <- as.list(
    c(
      BW = BW, 
      VSIc = VSIc,
      VLc = VLc,
      VRc = VRc,
      VSc = VSc,
      VFc = VFc,
      VBc = VBc,
      QC = QC,
      QSIc = QSIc,
      QLc = QLc,
      QRc = QRc,
      QSc = QSc,
      QFc = QFc,
      PIZEN = PIZEN,
      PLZEN = PLZEN,
      PRZEN = PRZEN,
      PSZEN = PSZEN,
      PFZEN = PFZEN,
      PIaZEL = PIaZEL,
      PLaZEL = PLaZEL,
      PRaZEL = PRaZEL,
      PSaZEL = PSaZEL,
      PFaZEL = PFaZEL,
      ksto = ksto,
      PappCaco2ZEN = PappCaco2ZEN,
      Vin = Vin,
      SAin = SAin,
      kin = kin,
      kin7 = kin7,
      SAinb = SAinb,
      PappCaco2_aZEL = PappCaco2_aZEL,
      KurZ = KurZ,
      KfeZ = KfeZ,
      KuraZ = KuraZ,
      KfeaZ = KfeaZ,
      S9SI = S9SI,
      VmaxSIZENGLUc = VmaxSIZENGLUc,
      KmSIZENGLU = KmSIZENGLU,
      VMB = VMB,
      VmaxLIaZELc = VmaxLIaZELc,
      VmaxLIbZELc = VmaxLIbZELc,
      KmLIaZEL = KmLIaZEL,
      KmLIbZEL = KmLIbZEL,
      VLS9 = VLS9,
      VmaxLaZELc = VmaxLaZELc,
      VmaxLbZELc = VmaxLbZELc,
      KmLaZEL = KmLaZEL,
      KmLbZEL = KmLbZEL,
      VmaxLZENGLUc = VmaxLZENGLUc,
      KmLZENGLU = KmLZENGLU,
      VmaxLaZELGLUc = VmaxLaZELGLUc,
      KmLaZELGLU = KmLaZELGLU
    )
  ) %>% as.data.frame()
  
Parm.human <- Parm.c.human %>% 
  pivot_longer(names_to = "Parameter", 
               values_to = "Value", 
               cols = everything())

write.csv(Parm.human, here("Input","Parm.human.csv"), row.names = FALSE)
