# --------------------------------------------------------------------------- #
# Script for calculating the PBK model parameters (lifestage)
# # By: Chrysanthi Pachoulide
# Date: 29-01-2026
# --------------------------------------------------------------------------- #

# !! max BW in Physio_params is 76 for Female

# Human Physiology
PARAMS <- function(expAGE = NULL, expBW = NULL, sex = NULL) { 
  
  # Input
  Human.input <- read_delim(here("Input", "TissueComposition.csv"), 
                            delim = ";",
                            escape_double = FALSE, trim_ws = TRUE) %>%
    filter(Species == "Human") # tissue composition mainly as R&R composition, unified by Utsey and optimised as per Musther
  
  
  Human.dat <- Human.input %>% filter(!tissue %in% c("RBCs", "Plasma")) # df for all tissues except for adipose and RBCs
  names <- substr(Human.dat$tissue, 1, 2)
  F_Vo.names <- paste0("FV", names)
  F_Vo.dat <- cbind(F_Vo.names, select(Human.dat, F_Vo)) %>%
    pivot_wider(names_from = F_Vo.names, values_from = F_Vo)
  
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
  
  ## Physiological constants ----
  Physio_params <- Physio_params %>% 
    mutate(BloodFlowSum = rowSums(select(., starts_with("Q_")))) %>%            # 0.9935, total blood flow as the sum of the fractional blood flows of all organs on which we have data
    mutate(VolumesSum = rowSums(select(., starts_with("V_"))))                  # 0.96, total volume as the sum of the fractional organ volumes of all organs on which we have data
  
  # Body weight 
  BW <- if (!is.na(expBW)) {
    expBW # BW is the actual BW at the time of measurement
  } else {
    Physio_params[[paste0("BW", suffix)]]}

  ### Fractional volumes ----           
  # Reference (if not stated otherwise): Punt et al. (2020) https://dx.doi.org/10.1021/acs.chemrestox.0c00307
  
  FV_SI <- Physio_params[[paste0("V_gutFraction", suffix)]] 
  
  FV_Duo <- 0.037                                                               # Reference: Desesso et al 2001
  FV_Jej <- 0.382
  FV_Ile <- 0.581                                                               # Reference: Desesso et al 2001
  
  FV_SI_IC <- 0.475
  FV_SI_EC <- 0.282

  FV_SIL <- 0.008586                                                            # fractional volume of small intestinal lumen 
  FV_duol <- 0.095                                                              # calculated (Vduo/TotalV)
  FV_jejl <- 0.362                                                              # calculated (Vjej/TotalV)
  FV_ilel <- 0.543                                                              # calculated (Vile/TotalV)
  
  Vfe_IleL <- 0.128 * 0.1                                                       # volume of "feces" in the ileum, 0.128kg is the total volume of feces excreted per day, 0.007 is to correct for the small fraction already formed in the ileum

  FV_L <- Physio_params[[paste0("V_liverFraction", suffix)]]                          
  FV_L_IC = 0.573                                                               # fractional volume of intracellular space in the liver; Reference: Utsey et al. 2020 https://doi.org/10.1124/dmd.120.090498, https://github.com/metrumresearchgroup/PBPK_PC/blob/master/data/unified_tissue_comp.csv
  FV_L_EC = 0.161                                                               # Fractional volume of extracellular space in the liver; Reference: Utsey et al. 2020 https://doi.org/10.1124/dmd.120.090498, https://github.com/metrumresearchgroup/PBPK_PC/blob/master/data/unified_tissue_comp.csv
  FV_K <- Physio_params[[paste0("V_liverFraction", suffix)]]                          
  FV_F <- Physio_params[[paste0("V_adiposeFraction", suffix)]]
  MAc <- Physio_params[[paste0("AdiposeMass", suffix)]]
  FV_R <- Physio_params[[paste0("V_heartFraction", suffix)]]+
    Physio_params[[paste0("V_kidneyFraction", suffix)]]+
    Physio_params[[paste0("V_spleenFraction", suffix)]]+
    Physio_params[[paste0("V_lungFraction", suffix)]]                           # fraction of rapidly perfused tissue
  FV_R_gluc <- Physio_params[[paste0("V_heartFraction", suffix)]]+
    Physio_params[[paste0("V_spleenFraction", suffix)]]+
    Physio_params[[paste0("V_lungFraction", suffix)]]                           # fraction of rapidly perfused tissue
  FV_B <- Physio_params[[paste0("V_bloodFraction", suffix)]]
  
  VBile <- 100/1000                                                             # L, volume of the bile from the ICRP, vbile is 50-100ml
  
  # Total body volume
  VTotc <- Physio_params$VolumesSum                  
  
  ### To calculate the surface aeria of the intestine ----
  # Reference (if not stated otherwise): Willman et al; Helander and Fandriks 2014 10.3109/00365521.2014.898326
  
  # Radius
  R <- 1.5                                                                      # cm; radius small intestine
  Rduo <- 1.75                                                                  # cm; radius
  Rjej <- 1.5                                                                   # cm; radius 
  Rile <- 1.5                                                                   # cm; radius
  
  # Length
  Lduo <- 20                                                                    # cm; length 
  Ljej <- 104                                                                   # cm; length 
  Lile <- 156                                                                   # cm; length 
  
  # Microvilli amplification factor
  # villi and microvilli together amplify the small intestinal surface area by 60–120 times https://doi.org/10.3109/00365521.2014.898326
  # Reference: Desesso et al. (2001) https://doi.org/10.1016/S0278-6915(00)00136-8
  
  Amp <- 120                                                                    # small intestine 
  Ampduo <- 1*Amp                                                               # duodenum       
  Ampjej <- 0.7*Amp                                                             # jejunum
  Ampile <- 0.5*Amp   
  
  # Surface aeria
  SAduo = 2*pi*Rduo*Lduo*Ampduo                                                 # cm2
  SAjej = 2*pi*Rjej*Ljej*Ampjej                                                 # cm2
  SAile = 2*pi*Rile*Lile*Ampile                                                 # cm2
  
  
  ### Fractional blood flow rates ----
  QC <- Physio_params[[paste0("CardOut", suffix)]]                              # L/d, This is corrected for hematocrit already so it's plasma
  Hct <- Physio_params[[paste0("Hct", suffix)]]                                 # Hematocrit
  
  QTotc <- Physio_params$BloodFlowSum
  QSIc <- Physio_params[[paste0("Q_gutFraction", suffix)]]
  QDuoc = Lduo/(Lduo+Ljej+Lile)                                                 # fraction; calculated based on the fractional length of each section as there is a oral-to-anal gradient in blood flow along the small intestine. Reference: Granger and Kvietys 2006 doi: 10.1016/B0-12-386860-2/00778-4
  QJejc = Ljej/(Lduo+Ljej+Lile)                                                 # fraction; calculated based on the fractional length of each section as there is a oral-to-anal gradient in blood flow along the small intestine. Reference: Granger and Kvietys 2006 doi: 10.1016/B0-12-386860-2/00778-4
  QIlec = Lile/(Lduo+Ljej+Lile)                                                 # fraction; calculated based on the fractional length of each section as there is a oral-to-anal gradient in blood flow along the small intestine. Reference: Granger and Kvietys 2006 doi: 10.1016/B0-12-386860-2/00778-4
  
  QLc <- Physio_params[[paste0("Q_liverFraction", suffix)]]
  QKc <- Physio_params[[paste0("Q_kidneyFraction", suffix)]] 
  QRc <- Physio_params[[paste0("Q_heartFraction", suffix)]] +
    Physio_params[[paste0("Q_kidneyFraction", suffix)]] +
    Physio_params[[paste0("Q_spleenFraction", suffix)]] 
  QRc_gluc <- Physio_params[[paste0("Q_heartFraction", suffix)]] +
    Physio_params[[paste0("Q_spleenFraction", suffix)]] 
  QFc <- Physio_params[[paste0("Q_adiposeFraction", suffix)]]
  
  
  ### Fluid transit times (flow rates) ----
  # Reference: Grandoni et al. (2019) DOI:  10.5599/admet.638  
  #            Punt et al. (2020) https://dx.doi.org/10.1021/acs.chemrestox.0c00307)
  # Intestinal segment transit time is calculated by the total transit time divided by the intestinal segment length as a proportion of total small intestinal length 
  
  # Times
  t_GER <- 0.5 / 24                                                             # day; stomach emptying time
  t_SI <- 199.2/60 / 24                                                         # day; small intestine/residence transit time
  t_si_duo <- 0.037*t_SI / 24                                                   # day; duodenum transit/residence time
  t_si_jej <- 0.382*t_SI / 24                                                   # day; jejunum transit/residence time 
  t_si_ile <- 0.581*t_SI / 24                                                   # day; ileum transit/residence time 
  
  # Emptying/transit/residence rates
  k_GER = 1/t_GER                                                               # /day; stomach 
  k_duo <- 1/t_si_duo                                                           # /day; duodenum 
  k_jej <- 1/t_si_jej                                                           # /day; jejunum
  k_ile <- 1/t_si_ile                                                           # /day; ileum
  
  # Glomerular filtration rate
  GFR <- Physio_params[[paste0("GFR", suffix)]]                                 # Age/creatinine dependent glomerular filtration rate (L/day)

  
  ## Tissue:plasma partition coefficients and fraction unbound ----
  # Tissue partitioning is calculated using the equations developed by Rodgers and Rowland, as coded by Utsey et al. (2020) DOI: 10.1124/dmd.120.090498: https://github.com/metrumresearchgroup/PBPK_PC/blob/master/script/CalcKp_P%26T.R
  # Assumption that all substances are diprotic acids
  
  ### Partition coefficients -------------------------
  
  # Tissue partitioning is calculated using the equations developed by Rodgers and Rowland, as coded by Utsey et al. (2020) DOI: 10.1124/dmd.120.090498: https://github.com/metrumresearchgroup/PBPK_PC/blob/master/script/CalcKp_P%26T.R
  # Assumption that all substances are diprotic acids
  
  PhysChem <- read_delim("Input/PhysChem.csv", 
                         delim = ";", escape_double = FALSE, trim_ws = TRUE)
  
  calc_Kp_fus <- function(logP, pKa1, pka2, BPR, Rat.input, F_Vo.dat){
  
    dat_all <- Human.input %>% filter(
      !tissue %in% c("RBCs", "Adipose", "Plasma"))                              # df for all tissues except for adipose, RBCs, and plasma
    dat_ad <- Human.input %>% filter(tissue == "Adipose")                       # df for adipose
    dat_rbc <- Human.input %>% filter(tissue == "RBCs")                         # df for RBCs
    dat_plas <- Human.input %>% filter(tissue == "Plasma")                      # df for plasma
    
    pH_IW <- 7                                                                  # pH of intracellular tissue water
    pH_P <- 7.4                                                                 # pH of plasma
    pH_RBC <- 7.22                                                              # pH of blood cells
    
    P <- 10^(logP)                                                              # octanol:water partition coeff
    logP_OW <- 1.115*logP - 1.35                                                # oil:water partition coeff
    P_OW <- 10^(logP_OW) 
    
    HCT <- 0.45                                                                 # hematocrit
    
    # Calculate fraction unbound in plasma and tissues
    # Reference: Lobell and Sivarajah 2003
    # Equation slightly different depending on the ionisation status in plasma, therefore between glycuronides and aglycones
    
    if(pKa1 >= 7){
      fup <- 1/(1 + 10^((0.4485*logP) - 0.4782))                                # fraction unbound in plasma
    }else {
      fup <- 1/(1 + 10^((0.3649*logP) + 0.41678))                               # fraction unbound in plasma
    }
    
    fut <- 1/(1 + ((1 - fup)/fup)*0.5)                                          # fraction unbound in tissues
    
    
    # Calculate partitioning and ionisation status
    Kpu_bc <- (HCT - 1 + BPR) / (HCT*fup)                                       # partitioning to blood cells                     
    X <- 10^(pH_IW - pKa1) + 10^(2*pH_IW - pKa1 - pKa2)                         # ionisation status in cells    
    Y <- 10^(pH_P - pKa1) + 10^(2*pH_P - pKa1 - pKa2)                           # ionisation status in plasma 
    Z <- 1 
    Ka_PR <- (                                                
      1/fup - 1 - ((P*dat_plas$f_n_l + (0.3*P + 0.7)*dat_plas$f_n_pl) / (1+Y))
    )                                                                           # albumin association constant
    
    # Calculate tissue:plasma partition coefficients for non-lipid tissues
    Kp_all <- (
      dat_all$f_ew + ((1 + X)/(1 + Y))*dat_all$f_iw + 
        ((P*dat_all$f_n_l + (0.3*P + 0.7)*dat_all$f_n_pl)) / (1 + Y) + 
        (Ka_PR*dat_all$AR*X)/(1 + Y)
    )*fup  
    
    # Calculate tissue:plasma partition coefficients for the adipose tissue
    Kp_ad <- (
      dat_ad$f_ew + ((1 + X)/(1 + Y))*dat_ad$f_iw + 
        ((P_OW*dat_ad$f_n_l + (0.3*P_OW + 0.7)*dat_ad$f_n_pl)) / (1 + Y) + 
        (Ka_PR*dat_ad$AR*X)/(1 + Y)
    )*fup  
    
    nms_all <- dat_all$tissue %>% substr(1,2) %>% tolower()
    nms_all <- paste("Kp", nms_all, sep="")
    nms <- c("Kpad",nms_all)
    Kp <- as.list(c(Kp_ad,Kp_all)) 
    names(Kp) <- nms  
    
    PSI <- Kp$Kpgu                                                              # Intestine:Blood partition coefficient
    PL <- Kp$Kpli                                                               # Liver:Blood partition coefficient
    PF <- Kp$Kpad                                                               # Adipose:Blood partition coefficient
    PK <- Kp$Kpki # only for glucuronides
    PS <- (
      F_Vo.dat$FVBo*Kp$Kpbo + F_Vo.dat$FVBr*Kp$Kpbr +
        F_Vo.dat$FVMu*Kp$Kpmu + F_Vo.dat$FVSk*Kp$Kpsk
    ) /
      (F_Vo.dat$FVBo + F_Vo.dat$FVBr + 
         F_Vo.dat$FVMu + F_Vo.dat$FVSk)                                         # Slowly perfused = partition coefficient of each of these organs * fractional volume of the respective organ / sum of the fractional volume of all these organs
    PR <- (
      F_Vo.dat$FVHe*Kp$Kphe + F_Vo.dat$FVKi*Kp$Kpki + 
        F_Vo.dat$FVLu*Kp$Kplu + F_Vo.dat$FVSp*Kp$Kpsp
    ) /
      (F_Vo.dat$FVHe + F_Vo.dat$FVKi + 
         F_Vo.dat$FVLu + F_Vo.dat$FVSp)                                         # Rapidly perfused = partition coefficient of each of these organs * fractional volume of the respective organ / sum of the fractional volume of all these organs
    PR_gluc <- ( # only for glucuronides
      F_Vo.dat$FVHe*Kp$Kphe +  
        F_Vo.dat$FVLu*Kp$Kplu + F_Vo.dat$FVSp*Kp$Kpsp
    ) /
      (F_Vo.dat$FVHe + F_Vo.dat$FVKi + 
         F_Vo.dat$FVLu + F_Vo.dat$FVSp)                                         # Rapidly perfused = partition coefficient of each of these organs * fractional volume of the respective organ / sum of the fractional volume of all these organs
    
    Kp_fus_BPR <- as.list(c(PSI = PSI, 
                            PL = PL, 
                            PF = PF, 
                            PK = PK,
                            PS = PS, 
                            PR = PR,
                            PR_gluc= PR_gluc,
                            fup = fup,
                            fut = fut,
                            BPR = BPR))
    
    return(Kp_fus_BPR)
    
  }
  
  
  mylist <- list() 
  
  
  for (i in 1:nrow(PhysChem)) {
    
    logP <- PhysChem$logP[i]
    pKa1 <- PhysChem$pKa.1[i]
    pKa2 <- PhysChem$pKa.2[i] 
    BPR  <- PhysChem$BPR[i] 
    
    Kp_fus_BPR <- calc_Kp_fus(logP = logP, 
                              pKa1 = pKa1, 
                              pka2 = pKa2, 
                              BPR = BPR, 
                              Rat.input = Rat.input,
                              F_Vo.dat = F_Vo.dat) 
    
    mylist[[i]] <- Kp_fus_BPR
    
  }
  
  Kp_fus_BPR_df <- do.call("rbind",mylist) 
  rownames(Kp_fus_BPR_df) <- c("ZEN", "aZEL", "ZEN_gluc", "aZEL_gluc")
  rm(mylist)
  rm(Kp_fus_BPR)
  
  
  ## Kinetic parameters ----
  
  ### Absorption/transfer rates ---
  
  # Apparent permebility accross the intestinal membrane
  # Reference: Pfeiffer et al. (2011) DOI 10.1002/mnfr.201000381
  PappCaco2_ZEN <- 10.4 * 1e-6                                                  # cm/s 
  PappCaco2_aZEL <- 5.4 * 1e-6                                                  # cm/s
  
  
  ### Metabolism ----
  
  #### Small intestine ----
  
  S9SI <- 35.2                                                                  # small intestinal S9 protein yield (mg S9 protein/gram intestine); Reference: S. A. Peters et al. DOI:10.1007/s40262-015-0351-6
  
  # Hydrolysis
  # Assuming the same as liver
  # Reference: Malekinejad et al. (2006) doi:10.1016/j.tvjl.2005.03.004
  
  Vmax_SI_aZEL.c <- 358.7*1e-6                                                  # umol/min/mg protein; from pmol-> umol
  Km_SI_aZEL <- 9                                                               # uM
  Vmax_SI_bZEL.c <- 209.3 *1e-6                                                 # umol/min/mg protein; from pmol-> umol
  Km_SI_bZEL <- 23                                                              # uM
  
  # Glucuronidation
  # Reference: Mendez-Catala et al. (2021) DOI 10.1002/mnfr.202100443
  Vmax_SI_ZEN_gluc.c <- 0.49*1e-3                                               # umol/min/mg protein; from nmol-> umol 
  Km_SI_ZEN_gluc <- 1.17                                                        # uM
  Vmax_SI_aZEL_gluc.c <- 6.963*1e-3                                             # umol/min/mg protein; from nmol-> umol
  Km_SI_aZEL_gluc <- 7.426                                                      # uM

  # Glucuronide hydrolysis by gmGUS
  # Reference: Ebuzoeme et al. (2021) https://doi.org/  10.3390/pharmaceutics13071043 (wogonoside glucuronide)
  S9gmGUS <- 4                                                                  # mg protein/ml feces
  Vmax_Hydro <- 5.35                                                            # umol/h
  Km_Hydro <- 0.41*15                                                           # uM

  #### Liver ----
  
  S9L <- 122                                                                    # liver S9 protein yield, (mg S9 protein/gram liver) reference: Moerenhout et al. (2024) DOI 10.1021/acs.est.4c06534 sum of microsomal 32mg/g liver and cytosolic 80.7mg/g liver content

  # Hydrolysis
  # Reference: Malekinejad et al. (2006) doi:10.1016/j.tvjl.2005.03.004
  Vmax_L_aZEL.c <- 358.7*1e-6                                                   # pmol->umol/min/mg protein
  Km_L_aZEL <- 9                                                                # uM
  Vmax_L_bZEL.c <- 209.3*1e-6                                                   # pmol->umol/min/mg protein   
  Km_L_bZEL <- 2.043                                                            # uM
  
  # Glucuronidation
  # Reference: Mendez-Catala et al. (2021) DOI 10.1002/mnfr.202100443
  Vmax_L_ZEN_gluc.c <- 2.97*1e-3                                                # nmol->umol/min/mg S9
  Km_L_ZEN_gluc <- 2.043                                                        # uM
  Vmax_L_aZEL_gluc.c <- 2.983*1e-3                                              # nmol->umol/min/mg S9
  Km_L_aZEL_gluc <- 2.415                                                       # uM 0.87
  
  
  ### Excretion ----
  
  #### Biliary excretion ----
  # specifying the duration of gallbladder emptying (assuming 3 meals per day at 8, 14 and 20h)
  
  Tgb <- 25/60   / 24                                                           # day time for the bile to travel through he gallbladder
  Bile.fasted <- 0                                                              # L/day
  Bile.fed <- 1                                                                 # L/day
  
  
  #### Active transport ----
  # Reference: Yanfei Li et al. 2019 doi:10.1016/j.jsbmb.2019.03.027.
  
  ##### Constants ----
  Vmax_MRP2.c <- 3366e-6                                                        # umol/min/mg MRP2; pmol -> umol
  Km_MRP2 <- 78.9                                                               # uM
  Vmax_MRP3.c <- 986.2e-6                                                       # umol/min/mg MRP3; pmol -> umol
  Km_MRP3 <- 5.7                                                                # uM
  Vmax_BCRP.c <- 2470e-6                    
  Km_BCRP <- 676.4 

  ##### IVIVE ----
  # Inside-out (active) fraction of transporter in in vitro system vesicles
  f_in_out_MRP2 <- 0.42 
  f_in_out_MRP3 <- 0.58    
  f_in_out_BCRP <- 0.326
  
  # Protein abundance in vesicles
  pmol_mg_MRP2 <- 133.9                                                         # pmol transporter/mg total protein
  pmol_mg_MRP3 <- 239.4                                                         # pmol transporter/mg total protein  
  pmol_mg_BCRP <- 11.4                                                          # pmol transporter/mg total protein
  
  # Protein abundance in vivo (tissues)
  # Reference: Couto et al. (2020) DOI https://doi.org/10.1124/dmd.119.089656
  #            Fallon et al. (2016) DOI 10.1007/s11095-016-1966-5
  
  pmolMP_mgP_IN_MRP2 <- 3.3                                                     # pmol MRP2/mg protein in intestine 
  pmolMP_mgP_LI_MRP2 <- 0.38                                                    # pmol MRP2/mg protein in liver 
  pmolMP_mgP_KI_MRP2 <- 0.48                                                    # pmol MRP2/mg protein in kidney
  
  pmolMP_mgP_LI_MRP3 <- 0.69                                                    # pmol MRP3/mg protein in liver
  
  pmolMP_mgP_IN_BCRP <- 0.59                                                    # pmol BCRP/mg protein in intestine
  pmolMP_mgP_LI_BCRP <- 0.13                                                    # pmol BCRP/mg protein in liver
  
  # Total protein in tissues
  mgP_mgLI <- 0.18                                                              # mg protein per mg liver
  mgP_mgIN <- 0.133                                                             # mg protein per mg intestine
  mgP_mgKI <- 0.17                                                              # mg protein per mg kidney
  
  # Relative expression factor
  REF_IN_MRP2 <- pmolMP_mgP_IN_MRP2 / pmol_mg_MRP2
  REF_LI_MRP2 <- pmolMP_mgP_LI_MRP2 / pmol_mg_MRP2
  REF_KI_MRP2 <- pmolMP_mgP_KI_MRP2 / pmol_mg_MRP2
  
  REF_LI_MRP3 <- pmolMP_mgP_LI_MRP3 / pmol_mg_MRP3
  
  REF_IN_BCRP <- pmolMP_mgP_IN_BCRP / pmol_mg_BCRP
  REF_LI_BCRP <- pmolMP_mgP_LI_BCRP / pmol_mg_BCRP
  
  # Scaling factor
  SF_IN_MRP2 <- mgP_mgIN * REF_IN_MRP2 / f_in_out_MRP2                          # mg MRP2 / mg small intestine
  SF_LI_MRP2 <- mgP_mgLI * REF_LI_MRP2 / f_in_out_MRP2 *10                      # mg MRP2 / mg liver
  SF_KI_MRP2 <- mgP_mgKI * REF_KI_MRP2 / f_in_out_MRP2 *10                      # mg MRP2 / mg liver
  
  SF_LI_MRP3 <- mgP_mgLI * REF_LI_MRP3 / f_in_out_MRP3 /10                      # mg MRP3 / mg liver 
  
  SF_IN_BCRP <- mgP_mgIN * REF_IN_BCRP / f_in_out_BCRP                          # mg BCRP / mg liver 
  SF_LI_BCRP <- mgP_mgLI * REF_LI_BCRP / f_in_out_BCRP *10                      # mg BCRP / mg liver 
  
  
  
  ## Collect all parameters ----
  Kp_fus_BPR_df <- as.data.frame(Kp_fus_BPR_df) %>%
    mutate(Substance = rownames(Kp_fus_BPR_df)) %>%
    pivot_longer(!Substance, names_to = "Parameter", values_to = "Value") %>%
    unite("Parameter", Substance:Parameter) 
  
  Other_params <- tibble(
    BW = BW,  
    FV_SI = FV_SI,
    FV_Duo = FV_Duo,
    FV_Jej = FV_Jej,
    FV_Ile = FV_Ile,
    FV_SI_IC = FV_SI_IC,
    FV_SI_EC = FV_SI_EC,
    FV_SIL = FV_SIL,
    FV_duol = FV_duol,
    FV_jejl = FV_jejl,
    FV_ilel = FV_ilel,
    Vfe_IleL = Vfe_IleL,
    FV_L = FV_L,
    FV_L_IC = FV_L_IC,
    FV_L_EC = FV_L_EC,
    FV_R = FV_R,
    FV_R_gluc = FV_R_gluc, 
    FV_F = FV_F,
    FV_B = FV_B,
    FV_K = FV_K,
    VBile = VBile,
    SAduo = SAduo,
    SAjej = SAjej,
    SAile = SAile,
    QC = QC,
    QSIc = QSIc,
    QDuoc = QDuoc,
    QJejc = QJejc,
    QIlec = QIlec,
    QLc = QLc,
    QKc = QKc,
    QRc = QRc,
    QRc_gluc = QRc_gluc,
    QFc = QFc,
    k_GER = k_GER,
    k_duo = k_duo,
    k_jej = k_jej,
    k_ile = k_ile,
    GFR = GFR,
    PappCaco2_ZEN = PappCaco2_ZEN,
    PappCaco2_aZEL = PappCaco2_aZEL,
    S9SI =  S9SI,
    Vmax_SI_bZEL.c =  Vmax_SI_bZEL.c,
    Km_SI_bZEL = Km_SI_bZEL,
    Vmax_SI_aZEL.c =  Vmax_SI_aZEL.c,
    Km_SI_aZEL = Km_SI_aZEL,
    Vmax_SI_aZEL_gluc.c =  Vmax_SI_aZEL_gluc.c,
    Km_SI_aZEL_gluc = Km_SI_aZEL_gluc,
    Vmax_SI_ZEN_gluc.c =  Vmax_SI_ZEN_gluc.c,
    Km_SI_ZEN_gluc = Km_SI_ZEN_gluc,
    S9gmGUS = S9gmGUS,
    Vmax_Hydro = Vmax_Hydro, 
    Km_Hydro = Km_Hydro,
    S9L = S9L,
    Vmax_L_aZEL.c = Vmax_L_aZEL.c,
    Km_L_aZEL = Km_L_aZEL,
    Vmax_L_bZEL.c = Vmax_L_bZEL.c,
    Km_L_bZEL = Km_L_bZEL,
    Vmax_L_ZEN_gluc.c = Vmax_L_ZEN_gluc.c,
    Km_L_ZEN_gluc = Km_L_ZEN_gluc,
    Vmax_L_aZEL_gluc.c = Vmax_L_aZEL_gluc.c,
    Km_L_aZEL_gluc = Km_L_aZEL_gluc,
    Bile.fasted = Bile.fasted,
    Bile.fed = Bile.fed,
    Tgb = Tgb,
    Vmax_MRP2.c = Vmax_MRP2.c,
    Km_MRP2 = Km_MRP2,
    Vmax_MRP3.c = Vmax_MRP3.c,
    Km_MRP3 = Km_MRP3,
    Vmax_BCRP.c = Vmax_BCRP.c,
    Km_BCRP = Km_BCRP,
    SF_IN_MRP2 = SF_IN_MRP2,
    SF_LI_MRP2 = SF_LI_MRP2,
    SF_KI_MRP2 = SF_KI_MRP2,
    SF_LI_MRP3 = SF_LI_MRP3,
    SF_IN_BCRP = SF_IN_BCRP,
    SF_LI_BCRP = SF_LI_BCRP
  ) %>% 
    pivot_longer(BW:SF_LI_BCRP,names_to = "Parameter", values_to = "Value")
  
  Parameters_df <- rbind(Kp_fus_BPR_df, Other_params) %>% 
    mutate(Value = as.numeric(Value))
  
  Parameters_df <- Parameters_df %>% filter(!grepl("ZEN_PK", Parameter),
                                            !grepl("ZEN_PR_gluc", Parameter),
                                            !grepl("aZEL_PK", Parameter),
                                            !grepl("aZEL_PR_gluc", Parameter),
                                            !grepl("\\bZEN_gluc_PR\\b", Parameter),
                                            !grepl("\\baZEL_gluc_PR\\b", Parameter))
  
  write.csv(Parameters_df, here("Input", "Parameters_human.csv"), row.names = FALSE)
  
  
  # parm.c <- split(Parameters_df, seq(nrow(Parameters_df)))
  return(Parameters_df)
}

