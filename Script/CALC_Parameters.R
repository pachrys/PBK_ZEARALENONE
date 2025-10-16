# --------------------------------------------------------------------------- #
# SCRIPT FOR CALCULATING BASE PARAMETERS
# By: Chrysanthi Pachoulide
# Date: 14-04-2025
# --------------------------------------------------------------------------- #

# !! max BW in Physio_params is 76 for Female

# COMMON PARAMETERS ####
BASE_PARAMS <- function(expAGE = NULL, expBW = NULL, sex = NULL) { 
  
  
  Human.input <- read_delim(here("Input", "TissueComposition.csv"), 
                            delim = ";",
                            escape_double = FALSE, trim_ws = TRUE) %>%
    filter(Species == "Human") # tissue composition mainly as R&R composition, unified by Utsey and optimised as per Musther
  
  
  Human.dat <- Human.input %>% filter(!tissue %in% c("RBCs", "Plasma")) # df for all tissues except for adipose and RBCs
  names <- substr(Human.dat$tissue, 1, 2)
  F_Vo.names <- paste0("FV", names)
  F_Vo.dat <- cbind(F_Vo.names, select(Human.dat, F_Vo)) %>%
    pivot_wider(names_from = F_Vo.names, values_from = F_Vo)
  
  ## Function conditions ####
  if (is.na(expAGE) && is.na(expBW)) {
    expBW <- 70  # Default body weight
  } 
  
  if (!is.na(expAGE) && expAGE > 80) {
    warning("Age above 80 which is the max age in the lifestage model. Physiological parameters assumed as those of an 80 year old")
    expAGE <- 80
  } 
  
  if (is.na(sex)) {
    sex <- "M"  # Default sex is male...
  }
  # 25ml plasma per cycle (29.2 days per cycle until; 12.5 cycles per year ); Ruark (12.89 +- 9.11 cycles per year)
  
  # Menstruation, reference Ruark et al 2016 http://dx.doi.org/10.1016/j.envint.2016.11.030
  if (sex == "F" && expAGE >= 12 && expAGE < 45){  
    CL_menses <-  25e-3/29.2 # L/day
  } else if (sex == "F" && expAGE >= 45 && expAGE <= 51){
    CL_menses <- 27.1e-3/35.1 # L/day
  } else {
    CL_menses <- 0 # L/day
  }
  
  # # Serum albumin (g_HSA /L_serum, or mg/ml); Weaving et al 2016; DOI: 10.1177/0004563215593561
  # if (expAGE <= 5) {
  #   SAlb <- 43  
  # } else if (expAGE <= 15 && sex == "M") {
  #   SAlb <- 44
  # } else if (expAGE <= 25 && sex == "M") {
  #   SAlb <- 46
  # } else if (expAGE <= 35 && sex == "M") {  
  #   SAlb <- 45
  # } else if (expAGE <= 50 && sex == "M") {
  #   SAlb <- 44
  # } else if (expAGE <= 60 && sex == "M") {
  #   SAlb <- 43
  # } else if (expAGE <= 75 && sex == "M") {
  #   SAlb <- 42
  # } else if (expAGE <= 80 && sex == "M") {
  #   SAlb <- 41
  # } else if (expAGE <= 90 && sex == "M") {
  #   SAlb <- 40
  # } else if (expAGE <= 10 && sex == "F") {
  #   SAlb <- 45
  # } else if (expAGE <= 20 && sex == "F") {
  #   SAlb <- 44
  # } else if (expAGE <= 40 && sex == "F") {
  #   SAlb <- 43
  # } else if (expAGE <= 75 && sex == "F") {
  #   SAlb <- 42
  # } else if (expAGE <= 85 && sex == "F") {
  #   SAlb <- 41
  # } else if (expAGE <= 90 && sex == "F") {
  #   SAlb <- 40
  # } else if (expAGE <= 94) { 
  #   SAlb <- 39
  # } else { 
  #   SAlb <- 38
  # } 
  # 
  ## Input
  Physio_params <- Physio.c
  
  suffix <- ifelse(sex == "F", "_F", "_M")
  
  if (!is.na(expAGE)) { # filter based on age
    Physio_params <- Physio_params %>% 
      mutate(age_diff = abs(age - expAGE)) %>% 
      filter(age_diff == min(age_diff)) %>% 
      slice(1)
  } else if (!is.na(expBW)) { # filter based on BW if age is not provided
    bw_col <- paste0("BW", suffix)
    Physio_params <- Physio_params %>%
      mutate(bw_diff = abs(!!sym(bw_col) - expBW)) %>% 
      filter(bw_diff == min(bw_diff)) %>% 
      slice(1)
  }
  
  Physio_params <- Physio_params %>% select(ends_with(suffix))
  
  
  ## Physiological Parameters ####
  Physio_params <- Physio_params %>% 
    mutate(BloodFlowSum = rowSums(select(., starts_with("Q_")))) %>% # 0.9935, total blood flow as the sum of the fractional blood flows of all organs on which we have data
    mutate(VolumesSum = rowSums(select(., starts_with("V_")))) # 0.96, total volume as the sum of the fractional organ volumes of all organs on which we have data
  
  BW <- if (!is.na(expBW)) {
    expBW # BW is the actual BW at the time of measurement
  } else {
    Physio_params[[paste0("BW", suffix)]]}
  # BW <- Physio_params[[paste0("BW", suffix)]]   # Body weight calculated only from the lifestage equations
  
  QC <- Physio_params[[paste0("CardOut", suffix)]]              # L/d, This is corrected for hematocrit already so it's plasma
  Hct <- Physio_params[[paste0("Hct", suffix)]]                 # Hematocrit
  
  BH <- Physio_params[[paste0("BH", suffix)]]                   # Height (cm)
  BSA <- exp(-3.75 + 0.42*log(BH)+0.52*log(BW))*1e4             # Body Surface area (cm2)
  
  GFR <- Physio_params[[paste0("GFR", suffix)]]                 # Age/creatinine dependent glomerular filtration rate (L/day)
  
  SAlb <- Physio_params[[paste0("SAlb", suffix)]]               # Serum albumin concentration (g_HSA /L_serum, or mg/ml)
  ### Organ volumes -------------------------
  
  FV_ST = 0.000714 #to change
  
  # Intestine
  FV_SI <- Physio_params[[paste0("V_gutFraction", suffix)]] 
  
  FV_Duo = 0.037 #(Desesso et al 2001)
  FV_Jej = 0.382
  FV_I_IC = 0.475
  FV_I_EC = 0.282
  FV_Ile = 0.581 #(Desesso et al 2001)
  FV_CO = 0.002 # calculated from Musther (yotal volume was 0.238ml)
  
  R <- 1.5 # radius small intestine
  Rduo <- 1.75 # radius, cm (Willman et al; Helander and Fandriks 2014 10.3109/00365521.2014.898326)
  Rjej <- 1.5  # radius, cm (Willman et al; Helander and Fandriks 2014 10.3109/00365521.2014.898326)
  Rile <- 1.5    # radius, cm (Willman et al; Helander and Fandriks 2014 10.3109/00365521.2014.898326)
  
  Lduo <- 20    # length cm (Willman et al)
  Ljej <- 104   # cm (Willman et al)
  Lile <- 156   # cm (Willman et al)
  
  
  Amp <- 120 # microvilli amplification factor, ref: https://doi.org/10.1016/S0278-6915(00)00136-8; villi and microvilli together amplify the small intestinal surface area by 60–120 times https://doi.org/10.3109/00365521.2014.898326
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
  
  FV_SIL = 0.008586 # fractional volume of small intestinal lumen (from Punt et al. 2020)
  FV_duol = 0.095 # calculated (Vduo/TotalV)
  FV_jejl =  0.362 # calculated (Vjej/TotalV)
  FV_ilel = 0.543 # calculated (Vile/TotalV)
  
  FV_COL = 0.014 # fractional volume of large intestinal lumen
  
  
  # Liver
  FV_L <- Physio_params[[paste0("V_liverFraction", suffix)]]  # cm2, (adult of 70kg) Willmann2004 doi: 10.1021/jm030999b        
  
  # # Liver is divided in intracellular and extracellular compartments. 
  # # Extracellular compartment combines both the vascular and interstitial space
  # # Fractional volume of intracellular space was taken from Utsey et al. 2020 https://doi.org/10.1124/dmd.120.090498, https://github.com/metrumresearchgroup/PBPK_PC/blob/master/data/unified_tissue_comp.csv
  FV_L_IC = 0.573                   # Fractional volume of intracellular space in the liver
  FV_L_EC = 1 - 0.573               # Fractional volume of extracellular space in the liver
  
  # Adipose
  FV_F <- Physio_params[[paste0("V_adiposeFraction", suffix)]]
  MAc <- Physio_params[[paste0("AdiposeMass", suffix)]]
  
  FV_R <- Physio_params[[paste0("V_heartFraction", suffix)]]+
    Physio_params[[paste0("V_kidneyFraction", suffix)]]+
    Physio_params[[paste0("V_spleenFraction", suffix)]]+
    Physio_params[[paste0("V_lungFraction", suffix)]] # fraction of rapidly perfused tissue
  
  # Plasma
  FV_B <- Physio_params[[paste0("V_bloodFraction", suffix)]]
  
  VBi <- 100/1000 # L, volume of the bile from the ICRP, vbile is 50-100ml
  
  
  # Total body volume
  VTotc <- Physio_params$VolumesSum                  
  
  ### Organ blood flows -------------------------
  QTotc <- Physio_params$BloodFlowSum
  QSIc <- Physio_params[[paste0("Q_gutFraction", suffix)]]
  QDuoc = 0.079 # Assuming the same as Rat?
  QJejc = 0.680 # Assuming the same as Rat?
  QIlec = 0.241 # Assuming the same as Rat?
  QLc <- Physio_params[[paste0("Q_liverFraction", suffix)]]
  QFc <- Physio_params[[paste0("Q_adiposeFraction", suffix)]]
  QRc = Physio_params[[paste0("Q_heartFraction", suffix)]] +
    Physio_params[[paste0("Q_kidneyFraction", suffix)]] +
    Physio_params[[paste0("Q_spleenFraction", suffix)]] # fractio
  
  
  ### Other physiological flows and constants -------------------------
  # Intestinal fluid flow humans; (reference: Grandoni ; Punt et al. (2020))
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
  lilTR <- 1/11 # large intestinal transit time (/h)siological pHs in different matrices -------------------------
  
  ## Chemical Specific ####
  
  ### Partition coefficients -------------------------
  
  # Tissue partitioning (Rodgers and Rowland as coded by kmutsey: https://github.com/metrumresearchgroup/PBPK_PC/blob/master/script/CalcKp_P%26T.R
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
  
  # Multiply by fup to get Kp humanher than Kpu
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
  
  
  ### Absorption/transfer humans ####
  PappCaco2_ZEN = 10.4 * 1e-6 * 3600 # (cm/s->h) reference: Pfeiffer et al. (2011) DOI 10.1002/mnfr.201000381
  PappCaco2_aZEL = 5.4 * 1e-6 * 3600 # (cm/s->h) reference: Pfeiffer et al. (2011) DOI 10.1002/mnfr.201000381
  
  Peff.cm.per.hr <- PappCaco2_ZEN #cm/h
  
  
  ### Metabolism in small intestine tissue ####
  S9SI =  35.2 # 35.20 mg/g intestine, From S. A. Peters et al. DOI:10.1007/s40262-015-0351-6
  
  conc.HSD.S9I = 1.25 # mg protein/mL Malekinejad et al. (2006) doi:10.1016/j.tvjl.2005.03.004
  VmaxSI_aZELc = 358.7*1e-6 * 60      # assuming the same as liver; pmol->umol/min->h/mg protein reference: Malekinejad et al. (2006) doi:10.1016/j.tvjl.2005.03.004 (post mitochondrial fraction)
  KmSI_aZEL = 9                   # assuming the same as liver; uM; reference: Malekinejad et al. (2006) doi:10.1016/j.tvjl.2005.03.004
  VmaxSI_bZELc = 209.3 *1e-6 * 60     # assuming the same as liver; pmol->umol/min->h/mg protein; reference: Malekinejad et al. (2006) doi:10.1016/j.tvjl.2005.03.004 (post mitochondrial fraction)
  KmSI_bZEL = 23                   # assuming the same as liver; uM; reference: Malekinejad et al. (2006) doi:10.1016/j.tvjl.2005.03.004
  
  conc.UGT.S9I = 0.04 # mg protein/ml Diana DOI: 10.1002/mnfr.202100443
  VmaxSI_ZENGLUc = 0.49*1e-3 * 60  # nmol->umol/min->h/mg S9; DOI: 10.1002/mnfr.202100443
  KmSI_ZENGLU = 1.17   # uM; DOI: 10.1002/mnfr.202100443
  
  VmaxSI_aZELGLUc = 6.963*1e-3 * 60 # nmol->umol/min->h/mg s9; DOI: 10.1002/mnfr.202100443
  KmSI_aZELGLU = 7.426              # uM; DOI: 10.1002/mnfr.202100443
  
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
  conc.HSD.S9L = 1.25 # mg protein/mL Malekinejad et al. (2006) doi:10.1016/j.tvjl.2005.03.004
  
  VmaxL_aZELc = 358.7*1e-6 * 60      # pmol->umol/min->h/mg protein reference: Malekinejad et al. (2006) doi:10.1016/j.tvjl.2005.03.004 (post mitochondrial fraction)
  KmL_aZEL = 9                      # uM; reference: Malekinejad et al. (2006) doi:10.1016/j.tvjl.2005.03.004
  VmaxL_bZELc = 209.3*1e-6 * 60      # pmol->umol/min->h/mg protein; reference: Malekinejad et al. (2006) doi:10.1016/j.tvjl.2005.03.004 (post mitochondrial fraction)
  KmL_bZEL = 23                   # uM; reference: Malekinejad et al. (2006) doi:10.1016/j.tvjl.2005.03.004
  
  conc.UGT.S9L = 0.2 # mg protein/ml Diana DOI: 10.1002/mnfr.202100443
  VmaxL_ZENGLUc = 2.97*1e-3 * 60  # umol/h/mg S9; nmol->umol/min->h/mg S9; DOI: 10.1002/mnfr.202100443
  KmL_ZENGLU = 2.043               # uM; DOI: 10.1002/mnfr.202100443
  
  VmaxL_aZELGLUc = 2.983*1e-3 * 60 # umol/h/mg S9; nmol->umol/min->h/mg s9; DOI: 10.1002/mnfr.202100443
  KmL_aZELGLU = 2.415              # uM; DOI: 10.1002/mnfr.202100443
  
  fuincS9L.HSD_ZEN <- (f.frac.ionized_ZEN+frac.unionized_ZEN)*(1/(1+conc.HSD.S9L*10^(0.072*logD_ZEN^2+0.067*logD_ZEN-1.126)))
  fuincS9L.UGT_ZEN <- (f.frac.ionized_ZEN+frac.unionized_ZEN)*(1/(1+conc.UGT.S9L*10^(0.072*logD_ZEN^2+0.067*logD_ZEN-1.126)))
  fuincS9L_aZEL <- (f.frac.ionized_aZEL+frac.unionized_aZEL)*(1/(1+conc.UGT.S9L*10^(0.072*logD_aZEL^2+0.067*logD_aZEL-1.126)))
  
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
  
  REF.IN.MRP2 <- pmolMP_mgP.IN.MRP2/pmol_mg.MRP2
  REF.LI.MRP2 <- pmolMP_mgP.LI.MRP2/pmol_mg.MRP2
  REF.KI.MRP2 <- pmolMP_mgP.KI.MRP2/pmol_mg.MRP2
  
  REF.IN.MRP3 <- pmolMP_mgP.IN.MRP3/pmol_mg.MRP3
  REF.LI.MRP3 <- pmolMP_mgP.LI.MRP3/pmol_mg.MRP3
  REF.KI.MRP3 <- pmolMP_mgP.KI.MRP3/pmol_mg.MRP3
  
  SF.IN.MRP2 <- mgP_mgLI * REF.IN.MRP2 * 1e6       # (mg liver * 1e6 -> kg liver)
  SF.LI.MRP2 <- mgP_mgLI * REF.LI.MRP2 * 1e6       # (mg liver * 1e6 -> kg liver)
  SF.KI.MRP2 <- mgP_mgLI * REF.KI.MRP2 * 1e6       # (mg liver * 1e6 -> kg liver)
  
  SF.IN.MRP3 <- mgP_mgLI * REF.IN.MRP3 * 1e6       # (mg liver * 1e6 -> kg liver)
  SF.LI.MRP3 <- mgP_mgLI * REF.LI.MRP3 * 1e6       # (mg liver * 1e6 -> kg liver)
  SF.KI.MRP3 <- mgP_mgLI * REF.KI.MRP3 * 1e6       # (mg liver * 1e6 -> kg liver)
  
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
  
  
  
  ### Final parameter constants -------------------------
  
  parm.c <- list(
    BW = BW,
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
    CL.KI.MRP3.c = CL.KI.MRP3.c
  )
  
  return(parm.c)
}

# # PARAMETERS FOR DERMAL EXPOSURE ####
# DERMAL_PARAMS <- function(expAGE = NULL, expBW = NULL, sex = "M", base.parm.c) { 
#   
#   ## Function conditions ####
#   if (is.na(expAGE) && is.na(expBW)) {
#     expBW <- 70  # Default body weight
#   } 
#   
#   if (!is.na(expAGE) && expAGE > 80) {
#     warning("Age above 80 which is the max age in the lifestage model. Physiological parameters assumed as those of an 80 year old")
#     expAGE <- 80
#   } 
#   
#   if (is.na(sex)) {
#     sex <- "M"  # Default sex to male
#   }
#   
#   ## Input
#   Physio_params <- Physio.c
#   
#   suffix <- ifelse(sex == "F", "_F", "_M")
#   
#   
#   if (!is.na(expAGE)) { # filters based on Age
#     Physio_params <- Physio_params %>% 
#       mutate(age_diff = abs(age - expAGE)) %>%  
#       filter(age_diff == min(age_diff)) %>% 
#       slice(1)
#   } else if (!is.na(expBW)) { # filters based on BW if Age is not provided
#     bw_col <- paste0("BW", suffix)
#     Physio_params <- Physio_params %>%
#       mutate(bw_diff = abs(!!sym(bw_col) - expBW)) %>% 
#       filter(bw_diff == min(bw_diff)) %>% 
#       slice(1)
#   }
#   
#   Physio_params <- Physio_params %>% select(ends_with(suffix))
#   
#   
#   BW <- if (!is.na(expBW)) {
#     expBW
#   } else {
#     Physio_params[[paste0("BW", suffix)]]}
#   # BW <- Physio_params[[paste0("BW", suffix)]]   # Body weight calculated only from the lifestage equations
#   
#   BH <- Physio_params[[paste0("BH", suffix)]]                   # Height (cm)
#   BSA <- exp(-3.75 + 0.42*log(BH)+0.52*log(BW))*1e4             # Body Surface area (cm2), reference: Gastellu et al. 2024, 10.1016/j.envres.2024.120393 (supplementary file Physio_equations_detailed.xlsx, eq. from Pendse et al. 2020)
#   # Trine calculated the surface area based on the BodyWeight with this formula: SA_SkB = 9.1*(BW*1000)^0.666  # cm2, total body area of the skin (Husoy)
#   # Another option: BSA = 0.02350*BH^0.4226*BW^0.51456, ref. https://www.rivm.nl/bibliotheek/rapporten/090013003.pdf
# 
#   # Physiological Parameters ####
#   Physio_params <- Physio_params %>% 
#     mutate(BloodFlowSum = rowSums(select(., starts_with("Q_")))) %>% # 0.9935, total blood flow as the sum of the fractional blood flows of all organs on which we have data
#     mutate(VolumesSum = rowSums(select(., starts_with("V_")))) # 0.96, total volume as the sum of the fractional organ volumes of all organs on which we have data
#   
#   
#   # Skin
#   VSkc <- Physio_params[[paste0("V_skinFraction", suffix)]]  # fractional skin volume
#   SA_SkB <- BSA                       # cm2, total body surface area
#   # SA_SkB <- 15670                   # cm2, total body surface area except head, SCCS 2021 table 4 (https://health.ec.europa.eu/document/download/89af1a70-a2b1-44da-a868-e7d80a8e736c_en?filename=sccs_o_250.pdf)
# 
#   H_SkB <- 83.1                     # cm, average thickness of the skin barrier, 83.7 +- 16.6 J.Sandby-Moller et al. 2003, Table II DOI: 10.1080/00015550310015419
#   VSkB <- SA_SkB*H_SkB * 1e-3       # L, volume of the skin barrier 
#   
#   QSkc <- Physio_params[[paste0("Q_skinFraction", suffix)]]  # fractional skin blood flow
#   
#   # Partition coefficients -------------------------
#   
#   # Code for calculating the partition coefficients is a slidely modified version from Utsey et al. 2020 https://github.com/metrumresearchgroup/PBPK_PC/blob/master/script/CalcKp_Schmitt.R
#   
#   logML <- 3.52  # 3.52 ± 0.08 from Ebert A., Allendorf F. et al (2020), https://dx.doi.org/10.1021/acs.est.0c00175  (Liposomes composed of POPC (1-palmitoyl-2-oleoyl-glycero-3-phosphocholine))
#   logSP <- 1.61  # 1.61 ± 0.15 from Allendorf, F., Goss, K.-U. and Ulrich, N. (2021), https://doi.org/10.1002/etc.4954  (Structural proteins from chicken breast fillet (actin & myosin 60-95%), Recovery 95%)
#   logALB <- 4.33 # 4.33 ± 0.05 from Allendorf, F., Goss, K.-U. and Ulrich, N. (2021), https://doi.org/10.1002/etc.4954 (BSA (fatty acid free) Molar ratio compound to BSA < 0.1 72-96h, Recovery 94%))
#   logSL <- -1.37 # -1.37 ± 0.01 from Allendorf, F., Goss, K.-U. and Ulrich, N. (2021), https://doi.org/10.1002/etc.4954 (Olive oil with a high fraction of unsaturated fatty acids, Recovery 95%)
#   logFABP <- 4.3 # calculated by Allendorf, F., Goss, K.-U. and Ulrich, N. (2021), https://doi.org/10.1002/etc.4954
#   
#   k_ML <- 10^logML       # neutral phospholipids:water partition coefficient
#   k_SP <- 10^logSP       # structural protein:water partition coefficient
#   k_ALB <- 10^logALB     # albumin:water partition coefficient
#   k_SL <- 10^logSL       # structural lipids:water partition coefficient
#   k_FABP <- 10^logFABP   # fatty acid-binding protein:water partition coefficient
#   
#   Physio.data <- Tissue.c %>% 
#     filter(Species == "Human") %>%
#     filter(!Tissue %in% c("Comment", "NamingInUtsey")) %>%
#     select(- c(Comment, Species)) %>% 
#     mutate(across(-Tissue, as.numeric))
#   
#   Physio.Tissues <- Physio.data %>% filter(!Tissue %in% c("Blood", "Plasma")) 
#   
#   Kp <- (Physio.Tissues$f_W +
#            (k_ML*Physio.Tissues$f_ML) +
#            (k_SL*Physio.Tissues$f_SL) +
#            (k_SP*Physio.Tissues$f_SP) +
#            (k_ALB*Physio.Tissues$f_ALB) +
#            (k_FABP*Physio.Tissues$f_FABP)
#   ) #* fup Removing fup from here, Kps will be recalculated directly in the model for correct sensitivity analysis
#   
#   Names <- Physio.Tissues$Tissue %>% substr(1,2)
#   Names <- paste("Kp",Names ,sep="")
#   Kp.df <- as.data.frame(list(Tissue = Names)) %>%
#     mutate(Value = Kp) %>% 
#     pivot_wider(names_from = Tissue, values_from = Value)
#   
#   
#   # Calculate Plasma/Rest of the body partition coefficient
#   # KpRe = partition coefficient of each of the lumped organs * fractional volume of the respective organ / sum of the fractional volume of all these organs
#   KpRe <- (Kp.df$KpLu * Physio_params[[paste0("V_lungFraction", suffix)]] +
#              Kp.df$KpSt * Physio_params[[paste0("V_stomachFraction", suffix)]] +
#              Kp.df$KpBr * Physio_params[[paste0("V_brainFraction", suffix)]] +
#              Kp.df$KpHe * Physio_params[[paste0("V_heartFraction", suffix)]] + 
#              Kp.df$KpMu * Physio_params[[paste0("V_muscleFraction", suffix)]] +
#              Kp.df$KpSp * Physio_params[[paste0("V_spleenFraction", suffix)]] +
#              Kp.df$KpGo * Physio_params[[paste0("V_reproFraction", suffix)]] +
#              Kp.df$KpBo * Physio_params[[paste0("V_boneFraction", suffix)]]) / (
#                  Physio_params[[paste0("V_lungFraction", suffix)]] +
#                  Physio_params[[paste0("V_stomachFraction", suffix)]] +
#                  Physio_params[[paste0("V_brainFraction", suffix)]] +
#                  Physio_params[[paste0("V_heartFraction", suffix)]] +
#                  Physio_params[[paste0("V_muscleFraction", suffix)]] +
#                  Physio_params[[paste0("V_spleenFraction", suffix)]] +
#                  Physio_params[[paste0("V_reproFraction", suffix)]] +
#                  Physio_params[[paste0("V_boneFraction", suffix)]]) 
#   Kp.df$KpRe <- KpRe
#   
#   # Choose between calculated partition coefficients or initial ones (rat)
#   # In Kp.df, PF, PI, PK, PL, PSK and PR are the initial PCs from Kudo 2007 (rat). These were recalculated to total tissue (/fup) to enable differentiating the fup for the sensitivity analysis
#   # Correcting for fraction unbound (as it was not incorporated in the input calculating file)
#   PSkc <- Kp.df$KpSk #PSk  # Skin
#   PRc <- Kp.df$KpRe  #PR   # Rest
#   
#   base.parm.c["PRc"] <- PRc
#   
#   
#   # Uptake from skin -------------------------
#   Papp_SkB = 3.82*1e-3 * 60*60                 # cm/s, Ragnarsdottir et al. 2024 https://doi.org/10.1016/j.envint.2024.108772 (calculations 3.82*1e-3 cm/h -> *60*60 cm/s)
#   tlag = 1/6.21 * 24                           # /d, lag time for dermal absorption, Ragnarsdottir et al. 2024 https://doi.org/10.1016/j.envint.2024.108772 (calculations 6.21h *24d)
#   
#   # Final parameter constants -------------------------
#   
#   parm.c <- c(
#     base.parm.c, 
#     VSkc = VSkc, 
#     SA_SkB = SA_SkB,
#     VSkB = VSkB,
#     QSkc = QSkc,
#     PSkc = PSkc,
#     Papp_SkB = Papp_SkB,
#     tlag = tlag)
#   
#   
#   # Function out -------------------------
#   return(parm.c)
# }
# 
# # PARAMETERS FOR INHALATION EXPOSURE ####
# INHALATION_PARAMS <- function(expAGE = NULL, expBW = NULL, sex = "M", base.parm.c) { 
#   
#   
#   ## Function conditions ####
#   if (is.na(expAGE) && is.na(expBW)) {
#     expBW <- 70  # Default body weight
#   } 
#   
#   if (!is.na(expAGE) && expAGE > 80) {
#     warning("Age above 80 which is the max age in the lifestage model. Physiological parameters assumed as those of an 80 year old")
#     expAGE <- 80
#   } 
#   
#   if (is.na(sex)) {
#     sex <- "M"  # Default sex to male
#   }
#   
#   ## Input
#   Physio_params <- Physio.c
#   
#   suffix <- ifelse(sex == "F", "_F", "_M")
#   
#   
#   if (!is.na(expAGE)) { # filters based on Age
#     Physio_params <- Physio_params %>% 
#       mutate(age_diff = abs(age - expAGE)) %>%  
#       filter(age_diff == min(age_diff)) %>% 
#       slice(1)
#   } else if (!is.na(expBW)) { # filters based on BW if Age is not provided
#     bw_col <- paste0("BW", suffix)
#     Physio_params <- Physio_params %>%
#       mutate(bw_diff = abs(!!sym(bw_col) - expBW)) %>% 
#       filter(bw_diff == min(bw_diff)) %>% 
#       slice(1)
#   }
#   
#   Physio_params <- Physio_params %>% select(ends_with(suffix))
#   
# 
#   # Physiological Parameters ####
#   Physio_params <- Physio_params %>% 
#     mutate(BloodFlowSum = rowSums(select(., starts_with("Q_")))) %>% # 0.9935, total blood flow as the sum of the fractional blood flows of all organs on which we have data
#     mutate(VolumesSum = rowSums(select(., starts_with("V_")))) # 0.96, total volume as the sum of the fractional organ volumes of all organs on which we have data
#   
#   # Lungs
#   VLuc <- Physio_params[[paste0("V_lungFraction", suffix)]]
#   
#   # Plasma
#   FV_B <- Physio_params[[paste0("V_plasmaFraction", suffix)]]
#   VAPc <- 0.39*FV_B                  # Arterial plasma
#   VFV_B <- 0.61*FV_B                  # Venous plasma  
#   
#   # Partition coefficients -------------------------
#   
#   # Code for calculating the partition coefficients is a slidely modified version from Utsey et al. 2020 https://github.com/metrumresearchgroup/PBPK_PC/blob/master/script/CalcKp_Schmitt.R
#   
#   logML <- 3.52  # 3.52 ± 0.08 from Ebert A., Allendorf F. et al (2020), https://dx.doi.org/10.1021/acs.est.0c00175  (Liposomes composed of POPC (1-palmitoyl-2-oleoyl-glycero-3-phosphocholine))
#   logSP <- 1.61  # 1.61 ± 0.15 from Allendorf, F., Goss, K.-U. and Ulrich, N. (2021), https://doi.org/10.1002/etc.4954  (Structural proteins from chicken breast fillet (actin & myosin 60-95%), Recovery 95%)
#   logALB <- 4.33 # 4.33 ± 0.05 from Allendorf, F., Goss, K.-U. and Ulrich, N. (2021), https://doi.org/10.1002/etc.4954 (BSA (fatty acid free) Molar ratio compound to BSA < 0.1 72-96h, Recovery 94%))
#   logSL <- -1.37 # -1.37 ± 0.01 from Allendorf, F., Goss, K.-U. and Ulrich, N. (2021), https://doi.org/10.1002/etc.4954 (Olive oil with a high fraction of unsaturated fatty acids, Recovery 95%)
#   logFABP <- 4.3 # calculated by Allendorf, F., Goss, K.-U. and Ulrich, N. (2021), https://doi.org/10.1002/etc.4954
#   
#   k_ML <- 10^logML       # neutral phospholipids:water partition coefficient
#   k_SP <- 10^logSP       # structural protein:water partition coefficient
#   k_ALB <- 10^logALB     # albumin:water partition coefficient
#   k_SL <- 10^logSL       # structural lipids:water partition coefficient
#   k_FABP <- 10^logFABP   # fatty acid-binding protein:water partition coefficient
#   
#   Physio.data <- Tissue.c %>% 
#     filter(Species == "Human") %>%
#     filter(!Tissue %in% c("Comment", "NamingInUtsey")) %>%
#     select(- c(Comment, Species)) %>% 
#     mutate(across(-Tissue, as.numeric))
#   
#   Physio.Tissues <- Physio.data %>% filter(!Tissue %in% c("Blood", "Plasma")) 
#   
#   Kp <- (Physio.Tissues$f_W +
#            (k_ML*Physio.Tissues$f_ML) +
#            (k_SL*Physio.Tissues$f_SL) +
#            (k_SP*Physio.Tissues$f_SP) +
#            (k_ALB*Physio.Tissues$f_ALB) +
#            (k_FABP*Physio.Tissues$f_FABP)
#   ) #* fup Removing fup from here, Kps will be recalculated directly in the model for correct sensitivity analysis
#   
#   Names <- Physio.Tissues$Tissue %>% substr(1,2)
#   Names <- paste("Kp",Names ,sep="")
#   Kp.df <- as.data.frame(list(Tissue = Names)) %>%
#     mutate(Value = Kp) %>% 
#     pivot_wider(names_from = Tissue, values_from = Value)
#   
#   
#   # Calculate Plasma/Rest of the body partition coefficient
#   # KpRe = partition coefficient of each of the lumped organs * fractional volume of the respective organ / sum of the fractional volume of all these organs
#   KpRe <- (Kp.df$KpSk * Physio_params[[paste0("V_skinFraction", suffix)]] +
#              Kp.df$KpSt * Physio_params[[paste0("V_stomachFraction", suffix)]] +
#              Kp.df$KpBr * Physio_params[[paste0("V_brainFraction", suffix)]] +
#              Kp.df$KpHe * Physio_params[[paste0("V_heartFraction", suffix)]] + 
#              Kp.df$KpMu * Physio_params[[paste0("V_muscleFraction", suffix)]] +
#              Kp.df$KpSp * Physio_params[[paste0("V_spleenFraction", suffix)]] +
#              Kp.df$KpGo * Physio_params[[paste0("V_reproFraction", suffix)]] +
#              Kp.df$KpBo * Physio_params[[paste0("V_boneFraction", suffix)]]) / (
#                Physio_params[[paste0("V_skinFraction", suffix)]] +
#                  Physio_params[[paste0("V_stomachFraction", suffix)]] +
#                  Physio_params[[paste0("V_brainFraction", suffix)]] +
#                  Physio_params[[paste0("V_heartFraction", suffix)]] +
#                  Physio_params[[paste0("V_muscleFraction", suffix)]] +
#                  Physio_params[[paste0("V_spleenFraction", suffix)]] +
#                  Physio_params[[paste0("V_reproFraction", suffix)]] +
#                  Physio_params[[paste0("V_boneFraction", suffix)]])  
#   Kp.df$KpRe <- KpRe
#   
#   # Choose between calculated partition coefficients or initial ones (rat)
#   # Correcting for fraction unbound (as it was not incorporated in the input calculating file)
#   PRc <- Kp.df$KpRe   # Rest
#   PLuc <- Kp.df$KpLu  # Lungs
#   
#   base.parm.c["PRc"] <- PRc
#   
#   # Final parameter constants -------------------------
#   
#   parm.c <- c(
#     base.parm.c, 
#     VLuc = VLuc,
#     VAPc = VAPc, 
#     VFV_B = VFV_B,
#     PLuc = PLuc)
#   
#   # Function out -------------------------
#   return(parm.c)
# }
