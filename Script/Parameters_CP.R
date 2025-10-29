# --------------------------------------------------------------------------- #
# Script for calculating parameters  
# exact copy of the calculations of the PBK Model of Mendez-Catala et al. 2021 Mendez-Catala et al. (2021) DOI 10.1002/mnfr.202100443
# By: Chrysanthi Pachoulide
# Date: 23-06-2025
# --------------------------------------------------------------------------- #

rm(list = ls()) 

# Packages
library(here)
library(tidyverse)
library(showtext)
library(readr)
font_add(family = "Garamond",  regular = "GARA.TTF")
showtext_auto()

# # Set output storage directory
# OUTPUT = here("Output",  format(Sys.Date(),  "%Y-%m-%d"),  format(Sys.time(),  "%H-%M-%S"))
# dir.create(OUTPUT,  recursive = TRUE)


# Parameters for Rat####
  Rat.input <- read_delim(here("Input", "TissueComposition.csv"),
                          delim = ";",
                        escape_double = FALSE, trim_ws = TRUE) %>%
  filter(Species == "Rat") # tissue composition mainly as R&R composition, unified by Utsey and optimised as per Musther et al. (2017) DOI: 10.1016/j.xphs.2017.05.003

  Rat.dat <- Rat.input %>% filter(!tissue %in% c("RBCs", "Plasma")) # df for all tissues except for adipose and RBCs
  names <- substr(Rat.dat$tissue, 1, 2)
  F_Vo.names <- paste0("FV", names)
  F_Vo.dat <- cbind(F_Vo.names, select(Rat.dat, F_Vo)) %>%
  pivot_wider(names_from = F_Vo.names, values_from = F_Vo)


  ## Physiological constants ####
  # Fractional Tissue volumes (reference: Punt et al. (2020) https://dx.doi.org/10.1021/acs.chemrestox.0c00307)
  BW = 0.250  ##{Kg} #  body weight rat

  FV_ST = 0.012  # fraction volume of small intestine
  
  # FV_SI = 0.026 # fractional volume of small intestinal tissue
  FV_SI = pull(Rat.input[Rat.input$tissue == "Gut", "F_Vo"]) # fractional volume of small intestinal tissue
  FV_Duo = 0.060
  FV_Jej = 0.769
  FV_I_IC = pull(Rat.input[Rat.input$tissue == "Gut", "f_iw"])
  FV_I_EC = pull(Rat.input[Rat.input$tissue == "Gut", "f_ew"])
  FV_Ile = 0.171
  FV_CO = 0.002 # calculated from Musther et al. (2017) DOI: 10.1016/j.xphs.2017.05.003 (total volume was 0.238ml)
  
  FV_SIL = 0.014 # fractional volume of small intestinal lumen (Punt et al. (2020) https://dx.doi.org/10.1021/acs.chemrestox.0c00307)
  FV_duol = 0.183
  FV_jejl = 0.334 # calculated from Musther et al. (2017) DOI: 10.1016/j.xphs.2017.05.003
  FV_ilel = 0.483 # calculated from Musther et al. (2017) DOI: 10.1016/j.xphs.2017.05.003
  FV_COL = 0.012 # fractional volume of large intestinal lumen (Punt et al. (2020) https://dx.doi.org/10.1021/acs.chemrestox.0c00307), agrees with the 3g of feces excreted per day in rats from Cavigelli et al. (2005) DOI https:// doi.org/10.1677/joe.1.05935
  
  # FV_L = 0.036  # fractional of liver tissue
  # FV_L_IC = pull(Rat.input[Rat.input$tissue == "Liver", "f_iw"])
  # FV_L_EC = pull(Rat.input[Rat.input$tissue == "Liver", "f_ew"])
  # FV_R = 0.02  # 0.082-VSIc-VLc # fraction of rapidly perfused tissue
  # FV_S = 0.667  # 0.737-VFc # fraction of slowly perfused tissue
  # FV_F = 0.070  # fraction of fat tissue
  # FV_B = 0.074  # fraction of blood
  
  FV_L = pull(Rat.input[Rat.input$tissue == "Liver", "F_Vo"])  # fractional of liver tissue
  FV_L_IC = pull(Rat.input[Rat.input$tissue == "Liver", "f_iw"])
  FV_L_EC = pull(Rat.input[Rat.input$tissue == "Liver", "f_ew"])
  FV_R = pull(Rat.input[Rat.input$tissue == "Heart", "F_Vo"]) +
    pull(Rat.input[Rat.input$tissue == "Kidney", "F_Vo"]) +
    pull(Rat.input[Rat.input$tissue == "Spleen", "F_Vo"]) +
    pull(Rat.input[Rat.input$tissue == "Lung", "F_Vo"]) # fraction of rapidly perfused tissue
  FV_F = pull(Rat.input[Rat.input$tissue == "Adipose", "F_Vo"])  # fraction of fat tissue
  FV_B = pull(Rat.input[Rat.input$tissue == "RBCs", "F_Vo"]) + pull(Rat.input[Rat.input$tissue == "Plasma", "F_Vo"])  # fraction of blood
  # FV_S = 0.96 - (FV_SI + FV_L + FV_R + FV_F + FV_B) # fraction of slowly perfused tissue (Bone/Skin/Adipose/Muscle/(brain and testes) (slowly equilibrating tissues with large time constants))

  R <- 0.18 #radius small intestine
  Rduo <- 0.275/2 #diameter / 2 = radius, cm (Musther et al. (2017) DOI: 10.1016/j.xphs.2017.05.003 et al)
  Rjej <- 0.4/2   #diameter / 2 = radius, cm (Musther et al. (2017) DOI: 10.1016/j.xphs.2017.05.003 et al)
  Rile <- 0.4/2   #diameter / 2 = radius, cm (Musther et al. (2017) DOI: 10.1016/j.xphs.2017.05.003 et al)
  Rcol <- 0.849/2 #diameter / 2 = radius, cm (Bartoli R et al (2016) DOI: DOI: 10.4253/wjge.v5.i5.226)
  
  Lduo <- 7    # length cm (Musther et al. (2017) DOI: 10.1016/j.xphs.2017.05.003 et al)
  Ljej <- 90   # cm (Musther et al. (2017) DOI: 10.1016/j.xphs.2017.05.003 et al)
  Lile <- 20   # cm (Musther et al. (2017) DOI: 10.1016/j.xphs.2017.05.003 et al)
  Lco <- 8.93 # cm Bartoli et al (2013) DOI: 10.4253/wjge.v5.i5.226
  
  Amp <- 120 # microvilli amplification factor, ref: Desesso et al. (2001) https://doi.org/10.1016/S0278-6915(00)00136-8; villi and microvilli together amplify the small intestinal surface area by 60–120 times https://doi.org/10.3109/00365521.2014.898326
  Ampduo <- 1*Amp
  Ampjej <- 1*Amp #0.7
  Ampile <- 1*Amp #0.5
  
  SAduo = 2*pi*Rduo*Lduo*Ampduo #cm2
  SAjej = 2*pi*Rjej*Ljej*Ampjej #cm2
  SAile = 2*pi*Rile*Lile*Ampile #cm2
  SAco = 2*pi*Rcol*Lco
  
  # Blood flow rates (reference: Punt et al. (2020) https://dx.doi.org/10.1021/acs.chemrestox.0c00307)
  QC = 4.8  ##{L/h} # cardiac output
  # QSIc = 0.101  # fraction of blood flow to small intestine Musther et al. (2017) DOI: 10.1016/j.xphs.2017.05.003 et al. SI
  # QDuoc = 0.079 # Musther et al. (2017) DOI: 10.1016/j.xphs.2017.05.003
  # QJejc = 0.680 # Musther et al. (2017) DOI: 10.1016/j.xphs.2017.05.003
  # QIlec = 0.241 # Musther et al. (2017) DOI: 10.1016/j.xphs.2017.05.003
  # # QCOc = 0.022 # fraction of blood flow to large intestine Musther et al. (2017) DOI: 10.1016/j.xphs.2017.05.003 et al. SI
  # QLc = 0.141  # 0.242 - QSIc # fraction of blood flow to liver 101 
  # QRc = 0.51  # 0.76 - QSIc - QLc # fraction of blood flow to rapidly perfused tissue
  # QSc = 0.17  # 0.24 - QFc # fraction of blood flow to slowly perfused tissue
  # QFc = 0.07  # fraction of blood flow to fat
  
  QSIc = pull(Rat.input[Rat.input$tissue == "Gut", "F_QC"])  # fraction of blood flow to small intestine
  QDuoc = 0.079 # Musther et al. (2017) DOI: 10.1016/j.xphs.2017.05.003
  QJejc = 0.680 # Musther et al. (2017) DOI: 10.1016/j.xphs.2017.05.003
  QIlec = 0.241 # Musther et al. (2017) DOI: 10.1016/j.xphs.2017.05.003
  QLc = pull(Rat.input[Rat.input$tissue == "Liver", "F_QC"])  
  QRc = pull(Rat.input[Rat.input$tissue == "Heart", "F_QC"]) +
    pull(Rat.input[Rat.input$tissue == "Kidney", "F_QC"]) +
    pull(Rat.input[Rat.input$tissue == "Spleen", "F_QC"]) #
    # pull(Rat.input[Rat.input$tissue == "Lung", "F_QC"]) # fraction of rapidly perfused tissue # fraction of blood flow to rapidly perfused tissue
  QFc = pull(Rat.input[Rat.input$tissue == "Adipose", "F_QC"])  # fraction of blood flow to fat
  # QSc = 1 - (QSIc + QLc + QRc + QFc) # fraction of blood flow to slowly perfused tissue

  # Intestinal fluid flow rates; Wistar (reference: Musther et al. (2017) DOI: 10.1016/j.xphs.2017.05.003 et al SI; Punt et al. (2020) https://dx.doi.org/10.1021/acs.chemrestox.0c00307)
  Tge <- 0.25 # stomach emptying time (h)
  Tsi <- 1.63 # small intestine transit time (h)
  Tsiduo <- 0.06*Tsi # duodenum transit time (h)
  Tsijej <- 0.768*Tsi # jejunum transit time (h) (used to be 0.192)
  Tsiile <- 0.172*Tsi # ileum transit time (h) (used to be 0.086)
  
  gaTR = 1/Tge # (/h)
  duoTR <- 1/Tsiduo # (/h)
  jejTR <- 1/Tsijej # (/h)
  ileTR <- 1/Tsiile # (/h)
  lilTR <- 1/19.7 #1/(0.8*19.7) # large intestinal transit time (/h)
  
  # Glomerular filtration rate: Punt et al. (2020) https://dx.doi.org/10.1021/acs.chemrestox.0c00307
  GFR = 0.0786 #L/h Grandoni et al . (2019) DOI:  10.5599/admet.638 
  
  ## Partition coefficients / fraction unbound / fraction unionised ####
  # Tissue partitioning (Rodgers and Rowland as coded by Utsey et al. (2020) DOI: 10.1124/dmd.120.090498: https://github.com/metrumresearchgroup/PBPK_PC/blob/master/script/CalcKp_P%26T.R
  # Assumption that Zearalenone and a-Zearalenol are diprotic acids and Zearalenone glucuronide a monoprotic acid
  #microsomes or S9 (Kilford et al., 2008 based on Halifax and Houston 2006)
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
  
  P <- 10^(logP_ZEN)   # octanol:water partition coeff
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
  Kp_ZEN <- as.list(c(Kp_ad,Kp_all)) 
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
  
  P <- 10^(logP_aZEL)   # octanol:water partition coeff
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
  
  P <- 10^(logP_ZENGLU)   # octanol:water partition coeff
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
  
  ## Kinetic parameters ####
  
  ### Absorption/transfer rates ####
  PappCaco2_ZEN = 10.4 * 1e-6 * 3600 # (cm/s->h) reference: Pfeiffer et al. (2011) DOI 10.1002/mnfr.201000381
  PappCaco2_aZEL = 5.4 * 1e-6 * 3600 # (cm/s->h) reference: Pfeiffer et al. (2011) DOI 10.1002/mnfr.201000381
  
  Peff.cm.per.hr <- PappCaco2_ZEN #cm/h
  
  ### Metabolism in small intestine tissue ####
  S9SI =  21.9 # small intestinal S9 protein yield (mg S9 protein/gram intestine) reference: Doerksen et al. (2020) doi: 10.1007/s13318-020-00666-w sum of microsomal 3.9 mg/g and cytosolic 18mg/g 
  
  # conc.HSD.S9I = 1.25 # mg protein/mL Malekinejad et al. (2006) doi:10.1016/j.tvjl.2005.03.004
  VmaxSI_aZELc = 32*1e-6 * 60      # assuming the same as liver; pmol-> umol/min->h/mg protein reference: Malekinejad et al. (2006) doi:10.1016/j.tvjl.2005.03.004 (post mitochondrial fraction)
  KmSI_aZEL = 592                   # assuming the same as liver; uM; reference: Malekinejad et al. (2006) doi:10.1016/j.tvjl.2005.03.004
  VmaxSI_bZELc = 72*1e-6 * 60     # assuming the same as liver; pmol->umol/min->h/mg protein; reference: Malekinejad et al. (2006) doi:10.1016/j.tvjl.2005.03.004 (post mitochondrial fraction)
  KmSI_bZEL = 21                   # assuming the same as liver; uM; reference: Malekinejad et al. (2006) doi:10.1016/j.tvjl.2005.03.004
  
  # conc.UGT.S9I = 0.04 # mg protein/ml Mendez-Catala et al. (2021) DOI 10.1002/mnfr.202100443
  VmaxSI_ZENGLUc =  8.617*1e-3 * 60  # nmol->umol/min->h/mg S9; Mendez-Catala et al. (2021) DOI 10.1002/mnfr.202100443
  KmSI_ZENGLU = 6.198  # uM; Mendez-Catala et al. (2021) DOI 10.1002/mnfr.202100443
  
  VmaxSI_aZELGLUc = 6.963*1e-3 * 60 # nmol->umol/min->h/mg s9; Mendez-Catala et al. (2021) DOI 10.1002/mnfr.202100443
  KmSI_aZELGLU = 7.426              # uM; Mendez-Catala et al. (2021) DOI 10.1002/mnfr.202100443
  
  # fuincS9I.HSD_ZEN <- (f.frac.ionized_ZEN+frac.unionized_ZEN)*(1/(1+conc.HSD.S9I*10^(0.072*logD_ZEN^2+0.067*logD_ZEN-1.126)))
  # fuincS9I.UGT_ZEN <- (f.frac.ionized_ZEN+frac.unionized_ZEN)*(1/(1+conc.UGT.S9I*10^(0.072*logD_ZEN^2+0.067*logD_ZEN-1.126)))
  # fuincS9I_aZEL <- (f.frac.ionized_aZEL+frac.unionized_aZEL)*(1/(1+conc.UGT.S9I*10^(0.072*logD_aZEL^2+0.067*logD_aZEL-1.126)))
  
  ### Metabolism in large intestine lumen (microbiota compartment) ####
  VMB = 0.0164                    # fraction of feces of BW reference: Hoskins and Zamcheck (1968)
  VmaxLI_aZELc = 0.23*1e-6 * 60   # pmol->umol/min->h/mg feces reference: Mendez-Catala et al. (2020) https://doi.org/10.1016/j.toxrep.2020.07.010
  KmLI_aZEL = 66                  # uM; reference: Mendez-Catala et al. (2020)
  VmaxLI_bZELc = 0.14*1e-6 * 60   # pmol->umol/min->h/mg feces reference: Mendez-Catala et al. (2020)
  KmLI_bZEL = 80                  # uM; reference: Mendez-Catala et al. (2020)
  
  ### Metabolism in liver ####
  S9L = 112.7            # liver S9 protein yield, (mg S9 protein/gram liver) reference: Moerenhout et al. (2024) DOI 10.1021/acs.est.4c06534 sum of microsomal 32mg/g liver and cytosolic 80.7mg/g liver content
  # conc.HSD.S9L = 1.25 # mg protein/mL Malekinejad et al. (2006) doi:10.1016/j.tvjl.2005.03.004
  
  VmaxL_aZELc = 32*1e-6 * 60      # pmol->umol/min->h/mg protein reference: Malekinejad et al. (2006) doi:10.1016/j.tvjl.2005.03.004 (post mitochondrial fraction)
  KmL_aZEL = 592                  # uM; reference: Malekinejad et al. (2006) doi:10.1016/j.tvjl.2005.03.004
  VmaxL_bZELc = 72*1e-6 * 60      # pmol->umol/min->h/mg protein; reference: Malekinejad et al. (2006) doi:10.1016/j.tvjl.2005.03.004 (post mitochondrial fraction)
  KmL_bZEL = 21                   # uM; reference: Malekinejad et al. (2006) doi:10.1016/j.tvjl.2005.03.004
  
  # conc.UGT.S9L = 0.2 # mg protein/ml Mendez-Catala et al. (2021) DOI 10.1002/mnfr.202100443
  VmaxL_ZENGLUc = 7.029*1e-3 * 60  # umol/h/mg S9; nmol->umol/min->h/mg S9; Mendez-Catala et al. (2021) DOI 10.1002/mnfr.202100443
  KmL_ZENGLU = 6.753               # uM; Mendez-Catala et al. (2021) DOI 10.1002/mnfr.202100443
  
  VmaxL_aZELGLUc = 6.963*1e-3 * 60 # umol/h/mg S9; nmol->umol/min->h/mg s9; Mendez-Catala et al. (2021) DOI 10.1002/mnfr.202100443
  KmL_aZELGLU = 7.426              # uM; Mendez-Catala et al. (2021) DOI 10.1002/mnfr.202100443

  # fuincS9L.HSD_ZEN <- (f.frac.ionized_ZEN+frac.unionized_ZEN)*(1/(1+conc.HSD.S9L*10^(0.072*logD_ZEN^2+0.067*logD_ZEN-1.126)))
  # fuincS9L.UGT_ZEN <- (f.frac.ionized_ZEN+frac.unionized_ZEN)*(1/(1+conc.UGT.S9L*10^(0.072*logD_ZEN^2+0.067*logD_ZEN-1.126)))
  # fuincS9L_aZEL <- (f.frac.ionized_aZEL+frac.unionized_aZEL)*(1/(1+conc.UGT.S9L*10^(0.072*logD_aZEL^2+0.067*logD_aZEL-1.126)))
  
  ### Biliary excretion ####
  Bile.fasted <- 0.00075 # in L/h; In the rat, the bile acid-dependent and-independent bile flow has been estimated to be approximately 50 and 70 μl/kg/min, respectively. in rats there is a constant flow of bile to the duodenum
  Bile.fed <- 1 # in L/h; increase of bile excretion during feeding
  
  Vmax.MRP2 <- 4767e-6  # pmol -> umol/min/mg MRP2, 2935.5
  Km.MRP2 <- 109.6      # umol/L, 76.925
  Vmax.MRP3 <- 1520e-6  # pmol -> umol/min/mg , 839.675
  Km.MRP3 <- 6          # umol/L, 3.25
  
  f.in_out.MRP2 <- 0.42 # inside-out (active) fraction of transporter in in vitro system vesicles
  f.in_out.MRP3 <- 0.58 # inside-out (active) fraction of transporter in in vitro system vesicles
  
  pmol_mg.MRP2 <- 133.9 # protein abundance in vesicles, in pmol transporter/mg protein
  pmol_mg.MRP3 <- 239.4 # protein abundance in vesicles, in pmol transporter/mg protein
  
  MWMRP2 <- 190e3 # kDa -> g/mol https://www.sciencedirect.com/science/article/abs/pii/S0308814619305898
  #assuming equal distribution of MRP2 in the small intestine
  pmolMP_mgP.IN.MRP2 <- 4.76e3/MWMRP2 # pmol transporter/mg protein; ng->pg, pg -> pmol Gavins et al. 2023 https://doi.org/10.1016/j.biopha.2023.114329 
  # pmolMP_mgP.IN.NRP2.fasted <- 4.76e3/MWMRP2 # ng->pg, og -> pmol
  # pmolMP_mgP.IN.NRP2.fed <- 2*4.76e3/MWMRP2 # ng->pg, og -> pmol
  # pmolMP_mgP.IN.MRP2 <- 0.9 # pmol MRP2/mg protein in intestine
  
  pmolMP_mgP.LI.MRP2 <- 2.08 # pmol MRP2/mg protein in liver Fallon 2016 DOI 10.1007/s11095-016-1966-5
  # pmolMP_mgP.LI.MRP2 <- 1.5 # pmol MRP2/mg protein in liver #24.1
  
  #assuming that the distribution of MRP2 in the jejunum equals that of MRP3 in colon
  MWMRP3 <- 169.4e3 # kDa -> g/mol https://pmc.ncbi.nlm.nih.gov/articles/PMC22016/
  pmolMP_mgP.IN.MRP3 <- 4.76e3/MWMRP3 # ng->pg, og -> pmol Gavins et al. 2023 https://doi.org/10.1016/j.biopha.2023.114329 
  # pmolMP_mgP.IN.MRP3 <- 1.1 # pmol MRP3/mg protein in intestine
  
  pmolMP_mgP.LI.MRP3 <- 0.32 # pmol MRP3/mg protein in liver Fallon 2016 DOI 10.1007/s11095-016-1966-5
  # pmolMP_mgP.LI.MRP3 <- 0.5 # pmol MRP3/mg protein in liver

  pmolMP_mgP.KI.MRP2 <- 0.27 # pmol MRP2/mg protein in kidney Fallon 2016 DOI 10.1007/s11095-016-1966-5
  # pmolMP_mgP.KI.MRP2 <- 1.4 # pmol MRP2/mg protein in kidney
  pmolMP_mgP.KI.MRP3 <- 1.40 # pmol MRP3/mg protien in kidney # Fallon 2016 DOI 10.1007/s11095-016-1966-5
  # pmolMP_mgP.KI.MRP3 <- 0.05 # pmol MRP3/mg protien in kidney
  mgP_mgLI <- 0.18 # mg protein per mg liver
  mgP_mgIN <- 0.133 # mg protein per mg intestine
  mgP_mgKI <- 0.17 # mg protein per mg kidney
  
  REF.IN.MRP2 <- pmolMP_mgP.IN.MRP2/pmol_mg.MRP2
  REF.LI.MRP2 <- pmolMP_mgP.LI.MRP2/pmol_mg.MRP2
  REF.KI.MRP2 <- pmolMP_mgP.KI.MRP2/pmol_mg.MRP2
  
  REF.IN.MRP3 <- pmolMP_mgP.IN.MRP3/pmol_mg.MRP3
  REF.LI.MRP3 <- pmolMP_mgP.LI.MRP3/pmol_mg.MRP3
  REF.KI.MRP3 <- pmolMP_mgP.KI.MRP3/pmol_mg.MRP3
  
  SF.IN.MRP2 <- mgP_mgIN * REF.IN.MRP2       # (mg protein / mg small intestine) #pmolMP_mgP.IN.MRP2*1e-6*MWMRP2*
  SF.LI.MRP2 <- mgP_mgLI * REF.LI.MRP2       # (mg protein / mg liver) #pmolMP_mgP.LI.MRP2*1e-6*MWMRP2*
  SF.KI.MRP2 <- mgP_mgLI * REF.KI.MRP2       # (mg protein / mg kidney) #pmolMP_mgP.KI.MRP2*1e-6*MWMRP2*
  
  SF.IN.MRP3 <- mgP_mgLI * REF.IN.MRP3       # (mg protein / mg liver) #pmolMP_mgP.IN.MRP3*1e-6*MWMRP3
  SF.LI.MRP3 <- mgP_mgLI * REF.LI.MRP3       # (mg protein / mg liver) #pmolMP_mgP.LI.MRP3*1e-6*MWMRP3*
  SF.KI.MRP3 <- mgP_mgLI * REF.KI.MRP3       # (mg protein / mg liver) #pmolMP_mgP.KI.MRP3*1e-6*MWMRP3*
  
  Vmax.MRP2.LI.c <- Vmax.MRP2 / f.in_out.MRP2 * SF.LI.MRP2  # umol/min/mg MRP2 * mg MRP2 / mg liver -> umol/min/mg liver 
  Vmax.MRP2.IN.c <- Vmax.MRP2 / f.in_out.MRP2 * SF.IN.MRP2  # umol/min/mg MRP2 * mg MRP2 / mg liver -> umol/min/mg liver
  Vmax.MRP3.LI.c <- Vmax.MRP3 / f.in_out.MRP3 * SF.LI.MRP3  # umol/min/mg MRP2 * mg MRP2 / mg liver -> umol/min/mg liver
  Vmax.MRP3.IN.c <- Vmax.MRP3 / f.in_out.MRP3 * SF.IN.MRP3  # umol/min/mg MRP2 * mg MRP2 / mg liver -> umol/min/mg liver
  
  CL.LI.MRP2.c <- (Vmax.MRP2 / f.in_out.MRP2 * SF.LI.MRP2 * 60)/Km.MRP2
  CL.IN.MRP2.c <- (Vmax.MRP2 / f.in_out.MRP2 * SF.IN.MRP2 * 60)/Km.MRP2
  CL.KI.MRP2.c <- (Vmax.MRP2 / f.in_out.MRP2 * SF.KI.MRP2 * 60)/Km.MRP2

  CL.LI.MRP3.c <- (Vmax.MRP3 / f.in_out.MRP3 * SF.LI.MRP3 * 60)/Km.MRP3
  CL.IN.MRP3.c <- (Vmax.MRP3 / f.in_out.MRP3 * SF.IN.MRP3 * 60)/Km.MRP3
  CL.KI.MRP3.c <- (Vmax.MRP3 / f.in_out.MRP3 * SF.KI.MRP3 * 60)/Km.MRP3

  # CL.LI.MRP2.c <- (Vmax.MRP2 * 60 /Km.MRP2)/ f.in_out.MRP2 / pmol_mg.MRP2 / pmolMP_mgP.LI.MRP2 / mgP_mgLI # L/min -> h/mg -> kg liver
  # CL.LI.MRP3.c <- (Vmax.MRP3 * 60 /Km.MRP3)/ f.in_out.MRP3 / pmol_mg.MRP3 / pmolMP_mgP.LI.MRP3 / mgP_mgLI # L/min -> h/mg -> kg liver
  #
  # CL.IN.MRP2.c <- (Vmax.MRP2 * 60 * 1e6 /Km.MRP2)/ f.in_out.MRP2 / pmol_mg.MRP2 / pmolMP_mgP.IN.MRP2 / mgP_mgIN # L/min -> h/mg -> kg intestine
  # CL.IN.MRP3.c <- (Vmax.MRP3 * 60 * 1e6 /Km.MRP3)/ f.in_out.MRP3 / pmol_mg.MRP3 / pmolMP_mgP.IN.MRP3 / mgP_mgIN # L/min -> h/mg -> kg intestine
  #
  # CL.KI.MRP2.c <- (Vmax.MRP2 * 60 * 1e6 /Km.MRP2)/ f.in_out.MRP2 / pmol_mg.MRP2 / pmolMP_mgP.KI.MRP2 / mgP_mgKI # L/min -> h/mg -> kg kidney
  # CL.KI.MRP3.c <- (Vmax.MRP3 * 60 * 1e6 /Km.MRP3)/ f.in_out.MRP3 / pmol_mg.MRP3 / pmolMP_mgP.KI.MRP3 / mgP_mgKI # L/min -> h/mg -> kg kidney
  # 
  ### Hydrolysis by GUS ####
  CL_Hydro <- 1 #L/h
  
  ## Final params
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
      FV_I_IC = FV_I_IC,
      FV_I_EC = FV_I_EC,
      FV_Ile = FV_Ile,
      FV_CO = FV_CO,
      FV_L = FV_L,
      FV_L_IC = FV_L_IC,
      FV_L_EC = FV_L_EC,
      FV_R = FV_R,
      # FV_S = FV_S,
      FV_F = FV_F,
      FV_B = FV_B,
      VBi = NA,
      QC = QC,
      QSIc = QSIc,
      QDuoc = QDuoc, 
      QJejc = QJejc, 
      QIlec = QIlec,
      # QCOc = QCOc,
      QLc = QLc,
      QRc = QRc,
      # QSc = QSc,
      QFc = QFc,
      GFR = GFR,
      SAduo = SAduo,
      SAjej = SAjej,
      SAile = SAile,
      SAco = SAco,
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
      Tgb = NA,
      Bile.fed = Bile.fed,
      Bile.fasted = Bile.fasted,
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
      CL.KI.MRP3.c = CL.KI.MRP3.c,
      Vmax.MRP2.LI.c = Vmax.MRP2.LI.c, 
      Vmax.MRP2.IN.c = Vmax.MRP2.IN.c,
      Vmax.MRP3.LI.c = Vmax.MRP3.LI.c, 
      Km.MRP2 = Km.MRP2, 
      Km.MRP3 = Km.MRP3) 
  ) %>% as.data.frame()
  
  Parm.rat <- Parm.c.rat %>% 
    pivot_longer(names_to = "Parameter", 
                 values_to = "Value", 
                 cols = everything())
  
  write.csv(Parm.rat, here("Input","Parm.rat.csv"), row.names = FALSE)
  

# Parameters for human ####

  
  Human.input <- read_delim(here("Input", "TissueComposition.csv"), 
                          delim = ";",
                          escape_double = FALSE, trim_ws = TRUE) %>%
    filter(Species == "Human") # tissue composition mainly as R&R composition, unified by Utsey and optimised as per Musther et al. (2017) DOI: 10.1016/j.xphs.2017.05.003
  
  
Human.dat <- Human.input %>% filter(!tissue %in% c("RBCs", "Plasma")) # df for all tissues except for adipose and RBCs
  names <- substr(Human.dat$tissue, 1, 2)
  F_Vo.names <- paste0("FV", names)
  F_Vo.dat <- cbind(F_Vo.names, select(Human.dat, F_Vo)) %>%
    pivot_wider(names_from = F_Vo.names, values_from = F_Vo)


  ## Physiological constants ####
  # Fractional Tissue volumes (reference: Punt et al. (2020) https://dx.doi.org/10.1021/acs.chemrestox.0c00307)
  BW = 70  ##{Kg} #  body weight human

  FV_ST = 0.000714  # fraction volume of stomach lumen

  FV_SI = pull(Human.input[Human.input$tissue == "Gut", "F_Vo"]) # fractional volume of small intestinal tissue
  FV_Duo = 0.037 #(Desesso et al 2001)
  FV_Jej = 0.382 #(Desesso et al 2001)
  FV_I_IC = pull(Human.input[Human.input$tissue == "Gut", "f_iw"])
  FV_I_EC = pull(Human.input[Human.input$tissue == "Gut", "f_ew"])
  FV_Ile = 0.581 #(Desesso et al 2001)

  R <- 1.5 # radius small intestine
  Rduo <- 1.75 # radius, cm (Willman et al; Helander and Fandriks 2014 10.3109/00365521.2014.898326)
  Rjej <- 1.5  # radius, cm (Willman et al; Helander and Fandriks 2014 10.3109/00365521.2014.898326)
  Rile <- 1.5    # radius, cm (Willman et al; Helander and Fandriks 2014 10.3109/00365521.2014.898326)

  Lduo <- 20    # length cm (Willman et al)
  Ljej <- 104   # cm (Willman et al)
  Lile <- 156   # cm (Willman et al)


  Amp <- 120 # microvilli amplification factor, ref: Desesso et al. (2001) https://doi.org/10.1016/S0278-6915(00)00136-8; villi and microvilli together amplify the small intestinal surface area by 60–120 times https://doi.org/10.3109/00365521.2014.898326
  Ampduo <- 1*Amp
  Ampjej <- 1*Amp #0.7
  Ampile <- 1*Amp #0.5

  SAduo = 2*pi*Rduo*Lduo*Ampduo #cm2
  SAjej = 2*pi*Rjej*Ljej*Ampjej #cm2
  SAile = 2*pi*Rile*Lile*Ampile #cm2

  Vduo = pi*Rduo^2*Lduo
  Vjej = pi*Rjej^2*Ljej
  Vile = pi*Rile^2*Lile
  TotalV = 2030.254

  FV_SIL = 0.008586 # fractional volume of small intestinal lumen (reference Punt et al. (2020) https://dx.doi.org/10.1021/acs.chemrestox.0c00307)
  FV_duol = 0.095 # calculated (Vduo/TotalV)
  FV_jejl =  0.362 # calculated (Vjej/TotalV)
  FV_ilel = 0.543 # calculated (Vile/TotalV)

  FV_COL = 0.014 # fractional volume of large intestinal lumen

  FV_L = pull(Human.input[Human.input$tissue == "Liver", "F_Vo"])  # fractional of liver tissue
  FV_L_IC = pull(Human.input[Human.input$tissue == "Liver", "f_iw"])
  FV_L_EC = pull(Human.input[Human.input$tissue == "Liver", "f_ew"])
  FV_R = pull(Human.input[Human.input$tissue == "Heart", "F_Vo"]) +
    pull(Human.input[Human.input$tissue == "Kidney", "F_Vo"]) +
    pull(Human.input[Human.input$tissue == "Spleen", "F_Vo"]) +
    pull(Human.input[Human.input$tissue == "Lung", "F_Vo"]) # fraction of rapidly perfused tissue
  FV_F = pull(Human.input[Human.input$tissue == "Adipose", "F_Vo"])  # fraction of fat tissue
  FV_B = pull(Human.input[Human.input$tissue == "RBCs", "F_Vo"]) + pull(Human.input[Human.input$tissue == "Plasma", "F_Vo"])  # fraction of blood
  FV_S = 0.96 - (FV_SI + FV_L + FV_R + FV_F + FV_B) # fraction of slowly perfused tissue (Bone/Skin/Adipose/Muscle/(brain and testes) (slowly equilibrating tissues with large time constants))

  VBi <- 100/1000 # L, volume of the bile from the ICRP, vbile is 50-100ml

  # Blood flow humans (reference: Punt et al. (2020) https://dx.doi.org/10.1021/acs.chemrestox.0c00307)
  QC = 389.99  ##{L/h} # cardiac output
  QSIc = pull(Human.input[Human.input$tissue == "Gut", "F_QC"])  # fraction of blood flow to small intestine
  QDuoc = 0.079 # Assuming the same as Rat?
  QJejc = 0.680 # Assuming the same as Rat?
  QIlec = 0.241 # Assuming the same as Rat?
  QLc = pull(Human.input[Human.input$tissue == "Liver", "F_QC"])  # 0.242 - QSIc # fraction of blood flow to liver 101
  QRc = pull(Human.input[Human.input$tissue == "Heart", "F_QC"]) +
    pull(Human.input[Human.input$tissue == "Kidney", "F_QC"]) +
    pull(Human.input[Human.input$tissue == "Spleen", "F_QC"]) # fraction of rapidly perfused tissue # fraction of blood flow to rapidly perfused tissue
  QFc = pull(Human.input[Human.input$tissue == "Adipose", "F_QC"])  # fraction of blood flow to fat
  QSc = 1 - (QSIc + QLc + QRc + QFc) # fraction of blood flow to slowly perfused tissue


  # Intestinal fluid flow humans; (reference: Grandoni et al . (2019) DOI:  10.5599/admet.638 ; Punt et al. (2020) https://dx.doi.org/10.1021/acs.chemrestox.0c00307)
  # Intestinal segment transit time is calculated by the total transit time divided by the intestinal segment length as a proportion of total small intestinal length 
  Tge <- 0.5 # stomach emptying time (h)
  Tsi <- 199.2/60 # small intestine transit time (h)
  Tsiduo <- 0.037*Tsi # duodenum transit time (h) calculated
  Tsijej <- 0.382*Tsi # jejunum transit time (h) calculated
  Tsiile <- 0.581*Tsi # ileum transit time (h) calculated

  gaTR = 1/Tge # (/h)
  duoTR <- 1/Tsiduo # (/h)
  jejTR <- 1/Tsijej # (/h)
  ileTR <- 1/Tsiile # (/h)
  lilTR <- 1/11 # large intestinal transit time (/h)

  # Glomerular filtration human: Punt et al. (2020) https://dx.doi.org/10.1021/acs.chemrestox.0c00307
  GFR = 7.5 # L/h Grandoni et al . (2019) DOI:  10.5599/admet.638

  ## Partition coefficients / fraction unbound / fraction unionised ####
  # Tissue partitioning (Rodgers and Rowland as coded by Utsey et al. (2020) DOI: 10.1124/dmd.120.090498: https://github.com/metrumresearchgroup/PBPK_PC/blob/master/script/CalcKp_P%26T.R
  # Assumption that Zearalenone and a-Zearalenol are diprotic acids and Zearalenone glucuronide a monoprotic acid
  # microsomes or S9 (Kilford et al., 2008 based on Halifax and Houston 2006)
  dat_all <- Human.input %>% filter(!tissue %in% c("RBCs", "Adipose", "Plasma"))  #df for all tissues except for adipose, RBCs, and plasma
  dat_ad <- Human.input %>% filter(tissue == "Adipose")  #df for adipose
  dat_rbc <- Human.input %>% filter(tissue == "RBCs") #df for RBCs
  dat_plas <- Human.input %>% filter(tissue == "Plasma") #df for aplasma

  pH_IW <- 7       #pH of intracellular tissue water
  pH_P <- 7.4      #pH of plasma
  pH_RBC <- 7.22    #pH of blood cells

  ### Calculate Kp for ZEN ####

  logP_ZEN <- 3.580 # Molgpka
  pKa1_ZEN <- 7.5 # Molgpka
  pKa2_ZEN <- 7.7 # Molgpka
  BP_ZEN <- 0.67 # Shin et al. 2014 DOI: 10.1080/15287390903212741
  fup_ZEN <- 0.156 # Mukherjee et al. 2014 doi:10.1371/journal.pone.0113632

  P <- 10^(logP_ZEN)   # octanol:water partition coeff
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
  Kp_ZEN <- as.list(c(Kp_ad,Kp_all)) 
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

  P <- 10^(logP_aZEL)   # octanol:water partition coeff
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

  P <- 10^(logP_ZENGLU)   # octanol:water partition coeff
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

  ## Kinetic parameters ####

  ### Absorption/transfer humans ####
  PappCaco2_ZEN = 10.4 * 1e-6 * 3600 # (cm/s->h) reference: Pfeiffer et al. (2011) DOI 10.1002/mnfr.201000381
  PappCaco2_aZEL = 5.4 * 1e-6 * 3600 # (cm/s->h) reference: Pfeiffer et al. (2011) DOI 10.1002/mnfr.201000381
  
  Peff.cm.per.hr <- PappCaco2_ZEN # cm/h
  
  
  ### Metabolism in small intestine tissue ####
  S9SI =  35.2 # 35.20 mg/g intestine, reference S. A. Peters et al. DOI:10.1007/s40262-015-0351-6
  
  # conc.HSD.S9I = 1.25 # mg protein/mL Malekinejad et al. (2006) doi:10.1016/j.tvjl.2005.03.004
  VmaxSI_aZELc = 358.7*1e-6 * 60      # assuming the same as liver; pmol->umol/min->h/mg protein reference: Malekinejad et al. (2006) doi:10.1016/j.tvjl.2005.03.004 (post mitochondrial fraction)
  KmSI_aZEL = 9                   # assuming the same as liver; uM; reference: Malekinejad et al. (2006) doi:10.1016/j.tvjl.2005.03.004
  VmaxSI_bZELc = 209.3 *1e-6 * 60     # assuming the same as liver; pmol->umol/min->h/mg protein; reference: Malekinejad et al. (2006) doi:10.1016/j.tvjl.2005.03.004 (post mitochondrial fraction)
  KmSI_bZEL = 23                   # assuming the same as liver; uM; reference: Malekinejad et al. (2006) doi:10.1016/j.tvjl.2005.03.004

  # conc.UGT.S9I = 0.04 # mg protein/ml Mendez-Catala et al. (2021) DOI 10.1002/mnfr.202100443
  VmaxSI_ZENGLUc = 0.49*1e-3 * 60  # nmol->umol/min->h/mg S9; Mendez-Catala et al. (2021) DOI 10.1002/mnfr.202100443
  KmSI_ZENGLU = 1.17   # uM; Mendez-Catala et al. (2021) DOI 10.1002/mnfr.202100443

  VmaxSI_aZELGLUc = 6.963*1e-3 * 60 # nmol->umol/min->h/mg s9; Mendez-Catala et al. (2021) DOI 10.1002/mnfr.202100443
  KmSI_aZELGLU = 7.426              # uM; Mendez-Catala et al. (2021) DOI 10.1002/mnfr.202100443

  # fuincS9I.HSD_ZEN <- (f.frac.ionized_ZEN+frac.unionized_ZEN)*(1/(1+conc.HSD.S9I*10^(0.072*logD_ZEN^2+0.067*logD_ZEN-1.126)))
  # fuincS9I.UGT_ZEN <- (f.frac.ionized_ZEN+frac.unionized_ZEN)*(1/(1+conc.UGT.S9I*10^(0.072*logD_ZEN^2+0.067*logD_ZEN-1.126)))
  # fuincS9I_aZEL <- (f.frac.ionized_aZEL+frac.unionized_aZEL)*(1/(1+conc.UGT.S9I*10^(0.072*logD_aZEL^2+0.067*logD_aZEL-1.126)))


  ### Metabolism in large intestine lumen (microbiota compartment) ####
  VMB = 0.0018                    # fraction of feces of BW reference: Tose et al (2015)
  VmaxLI_aZELc = 0.054*1e-6 * 60   # pmol->umol/min->h/mg feces reference: Mendez-Catala et al. (2020) https://doi.org/10.1016/j.toxrep.2020.07.010
  KmLI_aZEL = 135                  # uM; reference: Mendez-Catala et al. (2020)
  VmaxLI_bZELc = 0.0108*1e-6 * 60   # pmol->umol/min->h/mg feces reference: Mendez-Catala et al. (2020)
  KmLI_bZEL = 163                  # uM; reference: Mendez-Catala et al. (2020)

  ### Metabolism in liver ####
  S9L = 122                        # liver S9 protein yield, (mg S9 protein/gram liver) references: Moerenhout et al 2024: DOI:10.1021/acs.est.4c06534
  # conc.HSD.S9L = 1.25 # mg protein/mL Malekinejad et al. (2006) doi:10.1016/j.tvjl.2005.03.004

  VmaxL_aZELc = 358.7*1e-6 * 60      # pmol->umol/min->h/mg protein reference: Malekinejad et al. (2006) doi:10.1016/j.tvjl.2005.03.004 (post mitochondrial fraction)
  KmL_aZEL = 9                      # uM; reference: Malekinejad et al. (2006) doi:10.1016/j.tvjl.2005.03.004
  VmaxL_bZELc = 209.3*1e-6 * 60      # pmol->umol/min->h/mg protein; reference: Malekinejad et al. (2006) doi:10.1016/j.tvjl.2005.03.004 (post mitochondrial fraction)
  KmL_bZEL = 23                   # uM; reference: Malekinejad et al. (2006) doi:10.1016/j.tvjl.2005.03.004

  # conc.UGT.S9L = 0.2 # mg protein/ml Mendez-Catala et al. (2021) DOI 10.1002/mnfr.202100443
  VmaxL_ZENGLUc = 2.97*1e-3 * 60  # umol/h/mg S9; nmol->umol/min->h/mg S9; Mendez-Catala et al. (2021) DOI 10.1002/mnfr.202100443
  KmL_ZENGLU = 2.043               # uM; Mendez-Catala et al. (2021) DOI 10.1002/mnfr.202100443

  VmaxL_aZELGLUc = 2.983*1e-3 * 60 # umol/h/mg S9; nmol->umol/min->h/mg s9; Mendez-Catala et al. (2021) DOI 10.1002/mnfr.202100443
  KmL_aZELGLU = 2.415              # uM; Mendez-Catala et al. (2021) DOI 10.1002/mnfr.202100443

  # fuincS9L.HSD_ZEN <- (f.frac.ionized_ZEN+frac.unionized_ZEN)*(1/(1+conc.HSD.S9L*10^(0.072*logD_ZEN^2+0.067*logD_ZEN-1.126)))
  # fuincS9L.UGT_ZEN <- (f.frac.ionized_ZEN+frac.unionized_ZEN)*(1/(1+conc.UGT.S9L*10^(0.072*logD_ZEN^2+0.067*logD_ZEN-1.126)))
  # fuincS9L_aZEL <- (f.frac.ionized_aZEL+frac.unionized_aZEL)*(1/(1+conc.UGT.S9L*10^(0.072*logD_aZEL^2+0.067*logD_aZEL-1.126)))

  ### Biliary excretion ####
  # specifying the duration of gallbladder emptying (assuming 3 meals per day at 8, 14 and 20h)
  
  Tgb = 25/60 # h time for the bile to travel through he gallbladder
  Bile.fasted <- 0 # L/h
  Bile.fed <- 1    # L/h
  
  # Vmax_BSEPc <- 7                            # umol/min/mg BSEP, Average active transport of bile acids, assuming that the maximum velocity of PFOA transport by BSEP corresponds to that of bile acids ; de Bruijn et al. 2024 https://doi.org/10.14573/altex.2302011
  # Km_BSEPc <- 16.4                           # umol/L, uM, Average affinity constant of bile acids to BSEP, following the above assumption ; de Bruijn et al. 2024 https://doi.org/10.14573/altex.2302011

  # IVIVE
  # mgBSEP.HC <- 0.839 * 140000 * 1e-9         # mg BSEP/1e6 hepatocytes (0.839 pmole BSEP/1e6 hepatocytes (amound of BSEP per hepatocyte) * 140000 g/mole (MW BSEP) -> pg * 1e-9 -> mg) de Bruijn et al. 2024 https://doi.org/10.14573/altex.2302011]
  # HC.GL <- 99                                 # 1e6 hepatocytes in the liver ; de Bruijn et al. 2024 https://doi.org/10.14573/altex.2302011
  # SF_BSEP <- mgBSEP.HC * HC.GL * 1e3          # Scaling factor for BSEP mediated hepatic efflux for GCA and GCDC ; de Bruijn et al. 2024 https://doi.org/10.14573/altex.2302011

  Vmax.MRP2 <- 4767e-6 # pmol -> umol/min/mg MRP2, 2935.5
  Km.MRP2 <- 109.6 # umol/L, 76.925
  Vmax.MRP3 <- 1520e-6 # pmol -> umol/min/mg , 839.675
  Km.MRP3 <- 6.0 # umol/L, 3.25
  
  
  f.in_out.MRP2 <- 0.42 # inside-out (active) fraction of transporter in in citro system vesicles
  f.in_out.MRP3 <- 0.58 # inside-out (active) fraction of transporter in in citro system vesicles
  
  pmol_mg.MRP2 <- 133.9 # protein abundance in vesicles, in pmol transporter/mg protein
  pmol_mg.MRP3 <- 239.4 # protein abundance in vesicles, in pmol transporter/mg protein
  
  pmolMP_mgP.IN.MRP2 <- 0.9 # pmol MRP2/mg protein in intestine
  pmolMP_mgP.LI.MRP2 <- 1.5 # pmol MRP2/mg protein in liver
 
  pmolMP_mgP.IN.MRP3 <- 1.1 # pmol MRP3/mg protein in intestine
  pmolMP_mgP.LI.MRP3 <- 0.5 # pmol MRP3/mg protein in liver
  
  pmolMP_mgP.KI.MRP2 <- 1.4 # pmol MRP2/mg protein in kidney
  pmolMP_mgP.KI.MRP3 <- 0.05 # pmol MRP3/mg protien in kidney
  
  mgP_mgLI <- 0.18 # mg protein per mg liver
  mgP_mgIN <- 0.133 # mg protein per mg intestine
  mgP_mgKI <- 0.17 # mg protein per mg kidney

  CL.LI.MRP2.c <- (Vmax.MRP2 / f.in_out.MRP2 * SF.LI.MRP2 * 60)/Km.MRP2
  CL.IN.MRP2.c <- (Vmax.MRP2 / f.in_out.MRP2 * SF.IN.MRP2 * 60)/Km.MRP2
  CL.KI.MRP2.c <- (Vmax.MRP2 / f.in_out.MRP2 * SF.KI.MRP2 * 60)/Km.MRP2
  
  CL.LI.MRP3.c <- (Vmax.MRP3 / f.in_out.MRP3 * SF.LI.MRP3 * 60)/Km.MRP3
  CL.IN.MRP3.c <- (Vmax.MRP3 / f.in_out.MRP3 * SF.IN.MRP3 * 60)/Km.MRP3
  CL.KI.MRP3.c <- (Vmax.MRP3 / f.in_out.MRP3 * SF.KI.MRP3 * 60)/Km.MRP3
  
  # 
  # CL.LI.MRP2.c <- (Vmax.MRP2 * 60 * 1e6 /Km.MRP2)/ f.in_out.MRP2 / pmol_mg.MRP2 / pmolMP_mgP.LI.MRP2 / mgP_mgLI # L/min -> h/mg -> kg liver
  # CL.LI.MRP3.c <- (Vmax.MRP3 * 60 * 1e6 /Km.MRP3)/ f.in_out.MRP3 / pmol_mg.MRP3 / pmolMP_mgP.LI.MRP3 / mgP_mgLI # L/min -> h/mg -> kg liver
  # 
  # CL.IN.MRP2.c <- (Vmax.MRP2 * 60 * 1e6 /Km.MRP2)/ f.in_out.MRP2 / pmol_mg.MRP2 / pmolMP_mgP.IN.MRP2 / mgP_mgIN # L/min -> h/mg -> kg intestine
  # CL.IN.MRP3.c <- (Vmax.MRP3 * 60 * 1e6 /Km.MRP3)/ f.in_out.MRP3 / pmol_mg.MRP3 / pmolMP_mgP.IN.MRP3 / mgP_mgIN # L/min -> h/mg -> kg intestine
  # 
  # CL.KI.MRP2.c <- (Vmax.MRP2 * 60 * 1e6 /Km.MRP2)/ f.in_out.MRP2 / pmol_mg.MRP2 / pmolMP_mgP.KI.MRP2 / mgP_mgKI # L/min -> h/mg -> kg kidney
  # CL.KI.MRP3.c <- (Vmax.MRP3 * 60 * 1e6 /Km.MRP3)/ f.in_out.MRP3 / pmol_mg.MRP3 / pmolMP_mgP.KI.MRP3 / mgP_mgKI # L/min -> h/mg -> kg kidney

  ### Hydrolysis by GUS ####
  CL_Hydro <- 1 #L/h

  ## Final params ####
  
  Parm.c.human <- as.list(
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
      FV_I_IC = FV_I_IC,
      FV_I_EC = FV_I_EC,
      FV_Ile = FV_Ile,
      FV_CO = FV_CO,
      FV_L = FV_L,
      FV_L_IC = FV_L_IC,
      FV_L_EC = FV_L_EC,
      FV_R = FV_R,
      # FV_S = FV_S,
      FV_F = FV_F,
      FV_B = FV_B,
      VBi = VBi, 
      QC = QC,
      QSIc = QSIc,
      QDuoc = QDuoc, 
      QJejc = QJejc, 
      QIlec = QIlec,
      # QCOc = QCOc,
      QLc = QLc,
      QRc = QRc,
      # QSc = QSc,
      QFc = QFc,
      GFR = GFR,
      SAduo = SAduo,
      SAjej = SAjej,
      SAile = SAile,
      # SAco = SAco,
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
      Tgb = Tgb,
      Bile.fed = Bile.fed,
      Bile.fasted = Bile.fasted,
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
  
  Parm.human <- Parm.c.human %>%
    pivot_longer(names_to = "Parameter",
                 values_to = "Value",
                 cols = everything())

  # write.csv(Parm.human, here("Input","Parm.human.csv"), row.names = FALSE)
  
  Parms <- tibble(
    Parameter = Parm.human$Parameter,
    Rat = Parm.rat$Value, 
    Human = Parm.human$Value
  )
  
  write.csv(Parms, here("Input","Parms.csv"), row.names = FALSE)
  