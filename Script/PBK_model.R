# --------------------------------------------------------------------------- #
# SCRIPT FOR PBK MODEL FOR ORAL EXPOSURE
# By: Chrysanthi Pachoulide
# Date: 30-07-20266
# --------------------------------------------------------------------------- #

PBK_RUN <- function(y, parms, times, events){ 
  
  PBK.model <- function(t, state, parms){
    with(as.list(c(state, parms)), {
      
      ## Model Parameters ----
      
      # Physiological
      
      ## Volumes
      BW <- BW  
      
      VSIL <- FV_SIL*BW 
      VDuoL <- VSIL*FV_duol 
      VJejL <- VSIL*FV_jejl 
      VIleL <- VSIL*FV_ilel 
      
      VSI <- FV_SI*BW 
      VDuo <- VSI*FV_Duo
      VJej <- VSI*FV_Jej
      VJej_IC <- FV_SI_IC*VJej
      VJej_EC <- FV_SI_EC*VJej
      VIle <- VSI*FV_Ile
      
      VL <- FV_L*BW  
      VL_IC <- FV_L_IC*VL
      VL_EC <- FV_L_EC*VL
      
      VR <- FV_R*BW 
      VR_gluc <- FV_R_gluc*BW
      VF <- FV_F*BW  
      VB <- FV_B*BW   
      VS <- (0.96*BW) - (VSI + VL + VR + VF + VB)                               # fraction of slowly perfused tissue (Bone/Skin/Adipose/Muscle/(brain and testes) (slowly equilibrating tissues with large time constants))
      VS_gluc <- (0.96*BW) - (VSI + VL + VR_gluc + VF + VB)                     # fraction of slowly perfused tissue (Bone/Skin/Adipose/Muscle/(brain and testes) (slowly equilibrating tissues with large time constants))
      
      
      # For MRP3 kidney excretion use only
      VK <- FV_K*BW
      
      ## Blood flow rates 
      QC <- QC  
      
      QSI <- QSIc*QC  
      QDuo <- QSI*QDuoc
      QJej <- QSI*QJejc
      QIle <- QSI*QIlec
      
      QL <- QLc*QC  
      QK <- QKc*QC
      QR <- QRc*QC  
      QR_gluc <- QRc_gluc*QC
      QF <- QFc*QC
      QS <- QC - (QSI + QL + QR + QF) 
      QS_gluc <- QC - (QSI + QL + QR_gluc + QF) 
      
      
      # Kinetic
      
      ## Intestinal absorption 
      
      ### Passive permeability
      CLduo_ZEN <- PappCaco2_ZEN * SAduo * 3600 * 1e-3    * 24                  # L/day; from cm/s to cm3/h to 
      CLjej_ZEN <- PappCaco2_ZEN * SAjej * 3600 * 1e-3    * 24                  # L/day; from cm/s to cm3/h to 
      CLile_ZEN <- PappCaco2_ZEN * SAile * 3600 * 1e-3    * 24                  # L/day; from cm/s to cm3/h to 
      
      CLile_aZEL <- PappCaco2_aZEL * SAile * 3600 * 1e-3  * 24                  # L/day; from cm/min to cm3/h to 
      
      
      ## Metabolism
      
      ### Small intestinal tissue
      Vmax_SI_aZEL <- Vmax_SI_aZEL.c * 60 * S9SI * VJej * 1e3   * 24            # umol/day; from umol/ min->day /mg protein -> whole intestine (kg)
      Vmax_SI_bZEL <- Vmax_SI_bZEL.c * 60 * S9SI * VJej * 1e3   * 24            # umol/day; from umol/ min->day /mg protein -> whole intestine (kg)
      Vmax_SI_ZEN_gluc <- Vmax_SI_ZEN_gluc.c * 60 * 24 *
        S9SI * VJej * 1e3                                                       # umol/day; from umol/ min->day /mg protein -> whole intestine (kg)
      Vmax_SI_aZEL_gluc <- Vmax_SI_aZEL_gluc.c * 60 * 24 *
        S9SI * VJej  * 1e3                                                      # umol/day; from umol/ min->day /mg protein -> whole intestine (kg)
      
      Vmax_Hydro_IleL <- Vmax_Hydro * 60 * S9gmGUS * Vfe_IleL                   # umol/h; from umol/ min->h /mg protein -> whole intestinal micriobiome (kg)
      
      
      Vmax_Hydro_IleL <- Vmax_Hydro * 60 * S9gmGUS * Vfe_IleL                     # umol/h; from umol/ min->h /mg protein -> whole intestinal micriobiome (kg)
      
      ### Liver
      Vmax_L_aZEL <- Vmax_L_aZEL.c * 60 * S9L * VL * 1e3 * 24                   # umol/day; from umol/ min->day /mg protein -> whole liver (kg)
      Vmax_L_bZEL <- Vmax_L_bZEL.c * 60 * S9L * VL * 1e3  * 24                  # umol/day; from umol/ min->day /mg protein -> whole liver (kg)
      Vmax_L_ZEN_gluc <- Vmax_L_ZEN_gluc.c * 60 * 24 *
        S9L * VL * 1e3                                                          # umol/day; from umol/ min->day /mg protein -> whole liver (kg)
      Vmax_L_aZEL_gluc <- Vmax_L_aZEL_gluc.c * 60 * 24 *
        S9L * VL * 1e3                                                          # umol/day; from umol/ min->day /mg protein -> whole liver (kg)
      
<<<<<<< HEAD
=======
      
>>>>>>> 3121499f20b39ffe684b0f7d92dc61ba6f4f0fdc
      ## Active transport
      
      Vmax_MRP2_LI <- Vmax_MRP2.c * 60 * SF_LI_MRP2 * VL * 1e6 * 24             # umol/min/mg liver; scaling: umol/min->day/mg MRP2 * mg MRP2 / mg liver ->  liver (kg) 
      Vmax_MRP2_IN <- Vmax_MRP2.c * 60 * SF_IN_MRP2 * VSI * 1e6 * 24            # umol/min/mg intestine; scaling: umol/min->day/mg MRP2 * mg MRP2 / mg liver ->  intestine (kg)  
      
      Vmax_MRP3_LI <- Vmax_MRP3.c * 60 * SF_LI_MRP3 * VL * 1e6 * 24             # umol/min/mg liver; scaling: umol/min->day/mg MRP3 * mg MRP3 / mg liver ->  liver (kg)
      Vmax_MRP2_KI <- Vmax_MRP2.c * 60 * SF_KI_MRP2 * VK * 1e6 * 24             # umol/min/mg kidney; scaling: umol/min->day/mg MRP3 * mg MRP3 / mg liver ->  kidney (kg)
      
      Vmax_BCRP_LI <- Vmax_BCRP.c * 60 * SF_LI_BCRP * VL * 1e6 * 24             # umol/min/mg liver; scaling: umol/min->day/mg BCRP * mg BCRP / mg liver ->  liver (kg)
      Vmax_BCRP_IN <- Vmax_BCRP.c * 60 * SF_IN_BCRP * VSI * 1e6 * 24            # umol/min/mg liver; scaling: umol/min->day/mg BCRP * mg BCRP / mg liver ->  liver (kg)
      
      ## Biliary clearance
      
      CL_Bile <- ifelse(((t > (8/24) & t < ((8/24) +  Tgb)) |
                           (t > (14/24) & t < ((14/24) + Tgb)) |
                           (t > (20/24) & t < ((20/24) + Tgb))), Bile.fed, Bile.fasted)
      
      
      ## Tissue Concentrations ----
      
      ### Zearalenone ----
      CDuoL_ZEN = ADuoL_ZEN/VDuoL
      CJejL_ZEN = AJejL_ZEN/VJejL
      CIleL_ZEN = AIleL_ZEN/VIleL
      CDuo_ZEN = ADuo_ZEN/VDuo
      CJej_ZEN = AJej_ZEN/VJej
      CIle_ZEN = AIle_ZEN/VIle
      CL_ZEN = AL_ZEN/VL  
      CR_ZEN = AR_ZEN/VR   
      CS_ZEN = AS_ZEN/VS
      CF_ZEN = AF_ZEN/VF
      CB_ZEN = AB_ZEN/VB
      
      CVDuo_ZEN <- CDuo_ZEN/ZEN_PSI * ZEN_BPR
      CVJej_ZEN <- CJej_ZEN/ZEN_PSI * ZEN_BPR
      CVIle_ZEN <- CIle_ZEN/ZEN_PSI * ZEN_BPR
      CVL_ZEN <- CL_ZEN/ZEN_PL * ZEN_BPR
      CVF_ZEN <- CF_ZEN/ZEN_PF * ZEN_BPR
      CVR_ZEN <- CR_ZEN/ZEN_PR * ZEN_BPR
      CVS_ZEN <- CS_ZEN/ZEN_PS * ZEN_BPR
      
      ### Zearalenone-glucuronide ----
      CDuoL_ZEN_gluc = ADuoL_ZEN_gluc/VDuoL
      CJejL_ZEN_gluc <- AJejL_ZEN_gluc/VJejL
      CIleL_ZEN_gluc <- AIleL_ZEN_gluc/VIleL 
      CJej_IC_ZEN_gluc <- AJej_IC_ZEN_gluc/VJej_IC
      CJej_EC_ZEN_gluc <- AJej_EC_ZEN_gluc/VJej_EC
      CL_IC_ZEN_gluc <- AL_IC_ZEN_gluc/VL_IC
      CL_EC_ZEN_gluc <- AL_EC_ZEN_gluc/VL_EC
      CBile_ZEN_gluc <- ABile_ZEN_gluc/VBile
      CK_ZEN_gluc <- AK_ZEN_gluc/VK
      CR_ZEN_gluc <- AR_ZEN_gluc/VR_gluc
      CS_ZEN_gluc <- AS_ZEN_gluc/VS_gluc
      CF_ZEN_gluc <- AF_ZEN_gluc/VF
      CB_ZEN_gluc <- AB_ZEN_gluc/VB
      
      CVJej_EC_ZEN_gluc <- CJej_EC_ZEN_gluc/ZEN_gluc_PSI * ZEN_gluc_BPR
      CVL_EC_ZEN_gluc <- CL_EC_ZEN_gluc/ZEN_gluc_PL * ZEN_gluc_BPR
      CVK_ZEN_gluc <- CK_ZEN_gluc/ZEN_gluc_PK * ZEN_gluc_BPR
      CVF_ZEN_gluc <- CF_ZEN_gluc/ZEN_gluc_PF * ZEN_gluc_BPR
      CVR_ZEN_gluc <- CR_ZEN_gluc/ZEN_gluc_PR_gluc * ZEN_gluc_BPR
      CVS_ZEN_gluc <- CS_ZEN_gluc/ZEN_gluc_PS * ZEN_gluc_BPR
      
      ### α-Zearalenol ----
      CIleL_aZEL <- AIleL_aZEL/VIleL
      CJej_aZEL = AJej_aZEL/VJej
      CIle_aZEL <- AIle_aZEL/VIle
      CL_aZEL <- AL_aZEL/VL  
      CR_aZEL <- AR_aZEL/VR   
      CS_aZEL <- AS_aZEL/VS
      CF_aZEL <- AF_aZEL/VF
      CB_aZEL <- AB_aZEL/VB
      
      CVIle_aZEL <- CIle_aZEL/aZEL_PSI * aZEL_BPR
      CVJej_aZEL <- CJej_aZEL/aZEL_PSI * aZEL_BPR
      CVL_aZEL <- CL_aZEL/aZEL_PL * aZEL_BPR
      CVF_aZEL <- CF_aZEL/aZEL_PF * aZEL_BPR
      CVR_aZEL <- CR_aZEL/aZEL_PR * aZEL_BPR
      CVS_aZEL <- CS_aZEL/aZEL_PS * aZEL_BPR    
      
      ### a-Zearalenol-glucuronide ----
      
      CDuoL_aZEL_gluc <- ADuoL_aZEL_gluc/VDuoL
      CJejL_aZEL_gluc <- AJejL_aZEL_gluc/VJejL
      CIleL_aZEL_gluc <- AIleL_aZEL_gluc/VIleL 
      CJej_IC_aZEL_gluc <- AJej_IC_aZEL_gluc/VJej_IC
      CJej_EC_aZEL_gluc <- AJej_EC_aZEL_gluc/VJej_EC
      CL_IC_aZEL_gluc <- AL_IC_aZEL_gluc/VL_IC
      CL_EC_aZEL_gluc <- AL_EC_aZEL_gluc/VL_EC
      CK_aZEL_gluc <- AK_aZEL_gluc/VK
      CBile_aZEL_gluc <- ABile_aZEL_gluc/VBile
      CR_aZEL_gluc <- AR_aZEL_gluc/VR_gluc
      CS_aZEL_gluc <- AS_aZEL_gluc/VS_gluc
      CF_aZEL_gluc <- AF_aZEL_gluc/VF
      CB_aZEL_gluc <- AB_aZEL_gluc/VB
      
      CVJej_EC_aZEL_gluc <- CJej_EC_aZEL_gluc/aZEL_gluc_PSI * aZEL_gluc_BPR
      CVL_EC_aZEL_gluc <- CL_EC_aZEL_gluc/aZEL_gluc_PL * aZEL_gluc_BPR
      CVF_aZEL_gluc <- CF_aZEL_gluc/aZEL_gluc_PF * aZEL_gluc_BPR
      CVK_aZEL_gluc <- CK_aZEL_gluc/aZEL_gluc_PK * aZEL_gluc_BPR
      CVR_aZEL_gluc <- CR_aZEL_gluc/aZEL_gluc_PR_gluc * aZEL_gluc_BPR
      CVS_aZEL_gluc <- CS_aZEL_gluc/aZEL_gluc_PS * aZEL_gluc_BPR
      
      
      ## Oral exposure ----
      # if(t<expSTOP){DoseOn=1} else{DoseOn=0}
      #       Oral = exp*BW*DoseOn         
      # OralD = Oral/Tinput*(t %% tinterval<Tinput)
      # 
      
      
      ## Differential Equations ----
      
      dAOD_ZEN <- - k_GER*AOD_ZEN      #+ OralD                                  
      
      
      ### Zearalenone ----
      
      
      # Intestinal lumen 
      dADuoL_ZEN <- + k_GER*AOD_ZEN - k_duo*ADuoL_ZEN + 
        - CLduo_ZEN*CDuoL_ZEN 
      
      dAJejL_ZEN <- + k_duo*ADuoL_ZEN - k_jej*AJejL_ZEN + 
        - CLjej_ZEN*CJejL_ZEN 
      
      dAIleL_ZEN <- + k_jej*AJejL_ZEN - k_ile*AIleL_ZEN +  
        - CLile_ZEN*CIleL_ZEN  +
        + ((Vmax_Hydro_IleL*CIleL_ZEN_gluc) / 
             (Km_Hydro + (CIleL_ZEN_gluc)))
      
      # Small intestinal tissue 
      dADuo_ZEN <- CLduo_ZEN*CDuoL_ZEN + QDuo*(CB_ZEN-CVDuo_ZEN) 
      
      dAJej_ZEN <- CLjej_ZEN*CJejL_ZEN + QJej*(CB_ZEN-CVJej_ZEN) + 
        - ((Vmax_SI_ZEN_gluc*CJej_ZEN*ZEN_fut) / 
             (Km_SI_ZEN_gluc + (CJej_ZEN*ZEN_fut))) +
        - ((Vmax_SI_aZEL*CJej_ZEN*ZEN_fut) / 
             (Km_SI_aZEL + (CJej_ZEN*ZEN_fut))) +
        - ((Vmax_SI_bZEL*CJej_ZEN*ZEN_fut) / 
             (Km_SI_bZEL + (CJej_ZEN*ZEN_fut)))
      
      dAIle_ZEN <- CLile_ZEN*CIleL_ZEN + QIle*(CB_ZEN-CVIle_ZEN) 
      
      # Liver 
      dAL_ZEN <- QDuo*CVDuo_ZEN + QJej*CVJej_ZEN + QIle*CVIle_ZEN + 
        QL*CB_ZEN - (QSI+QL)*CVL_ZEN + 
        - ((Vmax_L_aZEL*CL_ZEN*ZEN_fut) / 
             (Km_L_aZEL + (CL_ZEN*ZEN_fut))) +
        - ((Vmax_L_bZEL*CL_ZEN*ZEN_fut) / 
             (Km_L_bZEL + (CL_ZEN*ZEN_fut))) +
        - ((Vmax_L_ZEN_gluc*CL_ZEN*ZEN_fut) / 
             (Km_L_ZEN_gluc + (CL_ZEN*ZEN_fut)))
      
      # Rapidly Perfused 
      dAR_ZEN <- QR*(CB_ZEN-CVR_ZEN) 
      
      # Slowly Perfused 
      dAS_ZEN <- QS*(CB_ZEN-CVS_ZEN)
      
      # Fat 
      dAF_ZEN <- QF*(CB_ZEN-CVF_ZEN) 
      
      # Blood
      dAB_ZEN <-  - (QSI + QL + QR + QS + QF)*CB_ZEN + #IV + 
        + (QSI + QL)*CVL_ZEN + 
        + QF*CVF_ZEN + QR*CVR_ZEN + QS*CVS_ZEN + 
        - GFR*CB_ZEN*ZEN_fup 
      
      # Urine
      dAUR_ZEN <- GFR*CB_ZEN*ZEN_fup 
      
      # Faeces 
      dAFE_ZEN <- k_ile*AIleL_ZEN
      
      # β-Zearalenol
      dA_bZEL <- ((Vmax_L_bZEL*CL_ZEN*ZEN_fut) / 
                    (Km_L_bZEL + (CL_ZEN*ZEN_fut)))  +
        ((Vmax_SI_bZEL*CJej_ZEN*ZEN_fut) / 
           (Km_SI_bZEL + (CJej_ZEN*ZEN_fut)))
      
      
      ### Zearalenone glucuronide ----
      
      # Intestinal lumen 
      dADuoL_ZEN_gluc <- CL_Bile*CBile_ZEN_gluc - k_duo*ADuoL_ZEN_gluc 
      
      dAJejL_ZEN_gluc <- + k_duo*ADuoL_ZEN_gluc - k_jej*AJejL_ZEN_gluc +
        + ((Vmax_MRP2_IN*CJej_IC_ZEN_gluc*ZEN_gluc_fut)/
             (Km_MRP2 + (CJej_IC_ZEN_gluc*ZEN_gluc_fut))) +
        + ((Vmax_BCRP_IN*CJej_IC_ZEN_gluc*ZEN_gluc_fut)/
             (Km_BCRP + (CJej_IC_ZEN_gluc*ZEN_gluc_fut))) 
      
      dAIleL_ZEN_gluc <- k_jej*AJejL_ZEN_gluc - k_ile*AIleL_ZEN_gluc +
        - ((Vmax_Hydro_IleL*CIleL_ZEN_gluc) / 
             (Km_Hydro + (CIleL_ZEN_gluc)))
      
      # Small Intestinal tissue
      dAJej_IC_ZEN_gluc <- ((Vmax_SI_ZEN_gluc*CJej_ZEN*ZEN_fut) / 
                              (Km_SI_ZEN_gluc + (CJej_ZEN*ZEN_fut))) +
        - ((Vmax_MRP2_IN*CJej_IC_ZEN_gluc*ZEN_gluc_fut) / 
             (Km_MRP2 + (CJej_IC_ZEN_gluc*ZEN_gluc_fut))) +
        - ((Vmax_BCRP_IN*CJej_IC_ZEN_gluc*ZEN_gluc_fut)/
             (Km_BCRP + (CJej_IC_ZEN_gluc*ZEN_gluc_fut)))
      
      dAJej_EC_ZEN_gluc <- QJej*(CB_ZEN_gluc-CVJej_EC_ZEN_gluc) 
      
      # Liver
      dAL_IC_ZEN_gluc <- ((Vmax_L_ZEN_gluc*CL_ZEN*ZEN_fut) / 
                            ((Km_L_ZEN_gluc) + (CL_ZEN*ZEN_fut))) +
        - ((Vmax_MRP2_LI*CL_IC_ZEN_gluc*ZEN_gluc_fut) / 
             (Km_MRP2 + (CL_IC_ZEN_gluc*ZEN_gluc_fut))) + 
        - ((Vmax_MRP3_LI*CL_IC_ZEN_gluc*ZEN_gluc_fut) / 
             (Km_MRP3 + (CL_IC_ZEN_gluc*ZEN_gluc_fut))) +
        - ((Vmax_BCRP_LI*CL_IC_ZEN_gluc*ZEN_gluc_fut) / 
             (Km_BCRP + (CL_IC_ZEN_gluc*ZEN_gluc_fut))) 
      
      dAL_EC_ZEN_gluc <- + QJej*CVJej_EC_ZEN_gluc + 
        + (QL+QDuo+QIle)*CB_ZEN_gluc - (QSI + QL)*CVL_EC_ZEN_gluc +
        + ((Vmax_MRP3_LI*CL_IC_ZEN_gluc*ZEN_gluc_fut) / 
             (Km_MRP3 + (CL_IC_ZEN_gluc*ZEN_gluc_fut)))
      
      # Bile
      dABile_ZEN_gluc <- ((Vmax_MRP2_LI*CL_IC_ZEN_gluc*ZEN_gluc_fut) / 
                            (Km_MRP2 + (CL_IC_ZEN_gluc*ZEN_gluc_fut))) +
        + ((Vmax_BCRP_LI*CL_IC_ZEN_gluc*ZEN_gluc_fut) / 
             (Km_BCRP + (CL_IC_ZEN_gluc*ZEN_gluc_fut))) +
        - CL_Bile*CBile_ZEN_gluc 
      
      # Rapidly Perfused 
      dAR_ZEN_gluc <- QR_gluc*(CB_ZEN_gluc-CVR_ZEN_gluc)
      
      # Slowly Perfused 
      dAS_ZEN_gluc <- QS_gluc*(CB_ZEN_gluc-CVS_ZEN_gluc)
      
      # Fat
      dAF_ZEN_gluc <- QF*(CB_ZEN_gluc-CVF_ZEN_gluc) 
      
      # Blood
      dAB_ZEN_gluc <- - (QSI + QL + QK + QR_gluc + QS_gluc + QF)*CB_ZEN_gluc +   
        + (QSI + QL)*CVL_EC_ZEN_gluc + 
        + QF*CVF_ZEN_gluc + QK*CVK_ZEN_gluc + QR_gluc*CVR_ZEN_gluc + QS_gluc*CVS_ZEN_gluc +
        - GFR*CB_ZEN_gluc*ZEN_gluc_fup 
      
      dAK_ZEN_gluc <- QK*(CB_ZEN_gluc-CVK_ZEN_gluc) +
        + GFR*CB_ZEN_gluc*ZEN_gluc_fup +
        - ((Vmax_MRP2_KI*CK_ZEN_gluc*ZEN_gluc_fut) / 
             (Km_MRP2 + (CK_ZEN_gluc*ZEN_gluc_fut)))
      
      # Urine
      dAUR_ZEN_gluc <- ((Vmax_MRP2_KI*CK_ZEN_gluc*ZEN_gluc_fut) / 
                          (Km_MRP2 + (CK_ZEN_gluc*ZEN_gluc_fut)))
      # 
      # Feces
      dAFE_ZEN_gluc <- k_ile*AIleL_ZEN_gluc
      
      ### α-Zearalenol ----
      
      # Intestinal lumen
      dAIleL_aZEL <- - k_ile*AIleL_aZEL +  
        - CLile_aZEL*CIleL_aZEL +  
        + ((Vmax_Hydro_IleL*CIleL_ZEN_gluc) / 
             (Km_Hydro + (CIleL_ZEN_gluc)))
      
      
      dAJej_aZEL <- ((Vmax_SI_aZEL*CJej_ZEN*ZEN_fut) / 
                       (Km_SI_aZEL + (CJej_ZEN*ZEN_fut))) +
        + QJej*(CB_aZEL-CVJej_aZEL) +
        - ((Vmax_SI_aZEL_gluc*CJej_aZEL*aZEL_fut) / 
             (Km_SI_aZEL_gluc + (CJej_aZEL*aZEL_fut)))
      
      dAIle_aZEL <- CLile_aZEL*CIleL_aZEL + 
        QIle*(CB_aZEL-CVIle_aZEL) 
      
      # Liver
      dAL_aZEL <- + QJej*CVJej_aZEL + QIle*CVIle_aZEL + 
        (QL+QDuo)*CB_aZEL - (QSI+QL)*CVL_aZEL + 
        + ((Vmax_L_aZEL*CL_ZEN*ZEN_fut) / 
             (Km_L_aZEL + (CL_ZEN*ZEN_fut))) +
        - ((Vmax_L_aZEL_gluc*CL_aZEL*aZEL_fut) / 
             (Km_L_aZEL_gluc + (CL_aZEL*aZEL_fut))) 
      
      # Rapidly Perfused 
      dAR_aZEL <- QR*(CB_aZEL-CVR_aZEL) 
      
      # Slowly Perfused 
      dAS_aZEL <- QS*(CB_aZEL-CVS_aZEL) 
      
      # Fat 
      dAF_aZEL <- QF*(CB_aZEL-CVF_aZEL) 
      
      # Blood
      dAB_aZEL <- - (QSI + QL + QR + QS + QF)*CB_aZEL +  
        + (QSI + QL)*CVL_aZEL + 
        + QF*CVF_aZEL + QR*CVR_aZEL + QS*CVS_aZEL +
        - GFR*CB_aZEL*aZEL_fup 
      
      # Urine
      dAUR_aZEL <- GFR*CB_aZEL*aZEL_fup
      
      # Feces
      dAFE_aZEL <- k_ile*AIleL_aZEL
      
      
      
      
      ### α-Zearalenol glucuronide ----
      
      
      # Intestinal lumen 
      dADuoL_aZEL_gluc <- CL_Bile*CBile_aZEL_gluc - k_duo*ADuoL_aZEL_gluc 
      
      dAJejL_aZEL_gluc <- + k_duo*ADuoL_aZEL_gluc - k_jej*AJejL_aZEL_gluc +
        + ((Vmax_MRP2_IN*CJej_IC_aZEL_gluc*aZEL_gluc_fut)/
             (Km_MRP2 + (CJej_IC_aZEL_gluc*aZEL_gluc_fut))) +
        + ((Vmax_BCRP_IN*CJej_IC_aZEL_gluc*aZEL_gluc_fut)/
             (Km_BCRP + (CJej_IC_aZEL_gluc*aZEL_gluc_fut)))
      
      dAIleL_aZEL_gluc <- k_jej*AJejL_aZEL_gluc - k_ile*AIleL_aZEL_gluc + 
        - ((Vmax_Hydro_IleL*CIleL_ZEN_gluc) / 
             (Km_Hydro + (CIleL_ZEN_gluc))) 
      
      # Small Intestinal tissue
      dAJej_IC_aZEL_gluc <- ((Vmax_SI_aZEL_gluc*CJej_aZEL*aZEL_fut) / 
                               (Km_SI_aZEL_gluc + (CJej_aZEL*aZEL_fut))) +
        - ((Vmax_MRP2_IN*CJej_IC_aZEL_gluc*aZEL_gluc_fut) / 
             (Km_MRP2 + (CJej_IC_aZEL_gluc*aZEL_gluc_fut))) +
        - ((Vmax_BCRP_IN*CJej_IC_aZEL_gluc*aZEL_gluc_fut)/
             (Km_BCRP + (CJej_IC_aZEL_gluc*aZEL_gluc_fut))) 
      
      dAJej_EC_aZEL_gluc <- QJej*(CB_aZEL_gluc-CVJej_EC_aZEL_gluc) 
      
      # Liver
      dAL_IC_aZEL_gluc <- ((Vmax_L_aZEL_gluc*CL_aZEL*aZEL_fut) / 
                             ((Km_L_aZEL_gluc) + (CL_aZEL*aZEL_fut))) +
        - ((Vmax_MRP2_LI*CL_IC_aZEL_gluc*aZEL_gluc_fut) / 
             (Km_MRP2 + (CL_IC_aZEL_gluc*aZEL_gluc_fut))) + 
        - ((Vmax_MRP3_LI*CL_IC_aZEL_gluc*aZEL_gluc_fut) / 
             (Km_MRP3 + (CL_IC_aZEL_gluc*aZEL_gluc_fut))) +
        - ((Vmax_BCRP_LI*CL_IC_aZEL_gluc*aZEL_gluc_fut) / 
             (Km_BCRP + (CL_IC_aZEL_gluc*aZEL_gluc_fut))) 
      
      dAL_EC_aZEL_gluc <- + QJej*CVJej_EC_aZEL_gluc +
        + (QL+QDuo+QIle)*CB_aZEL_gluc - (QSI + QL)*CVL_EC_aZEL_gluc +
        + ((Vmax_MRP3_LI*CL_IC_aZEL_gluc*aZEL_gluc_fut) / 
             (Km_MRP3 + (CL_IC_aZEL_gluc*aZEL_gluc_fut)))
      
      # Bile
      dABile_aZEL_gluc <- ((Vmax_MRP2_LI*CL_IC_aZEL_gluc*aZEL_gluc_fut) / 
                             (Km_MRP2 + (CL_IC_aZEL_gluc*aZEL_gluc_fut))) +
        + ((Vmax_BCRP_LI*CL_IC_aZEL_gluc*aZEL_gluc_fut) / 
             (Km_BCRP + (CL_IC_aZEL_gluc*aZEL_gluc_fut))) +
        - CL_Bile*CBile_aZEL_gluc 
      
      # Rapidly Perfused 
      dAR_aZEL_gluc <- QR_gluc*(CB_aZEL_gluc-CVR_aZEL_gluc)
      
      # Slowly Perfused 
      dAS_aZEL_gluc <- QS_gluc*(CB_aZEL_gluc-CVS_aZEL_gluc)
      
      # Fat
      dAF_aZEL_gluc <- QF*(CB_aZEL_gluc-CVF_aZEL_gluc) 
      
      # Blood
      dAB_aZEL_gluc <- - (QSI + QL + QK + QR_gluc + QS_gluc + QF)*CB_aZEL_gluc +   
        + (QSI + QL)*CVL_EC_aZEL_gluc + 
        + QF*CVF_aZEL_gluc + QK*CVK_aZEL_gluc + QR_gluc*CVR_aZEL_gluc + QS_gluc*CVS_aZEL_gluc +
        - GFR*CB_aZEL_gluc*aZEL_gluc_fup 
      
      dAK_aZEL_gluc <- QK*(CB_aZEL_gluc-CVK_aZEL_gluc) + 
        + GFR*CB_aZEL_gluc*aZEL_gluc_fup +
        - ((Vmax_MRP2_KI*CK_aZEL_gluc*aZEL_gluc_fut) / 
             (Km_MRP2 + (CK_aZEL_gluc*aZEL_gluc_fut)))
      
      # Urine
      dAUR_aZEL_gluc <- ((Vmax_MRP2_KI*CK_aZEL_gluc*aZEL_gluc_fut) / 
                           (Km_MRP2 + (CK_aZEL_gluc*aZEL_gluc_fut)))
      
      # Feces
      dAFE_aZEL_gluc <- k_ile*AIleL_aZEL_gluc
      
      ## Mass balance ----
      # dATotal <- OralD
      # ATotal <- 1 #ODOSEumol
      ACalculated <- AOD_ZEN + 
        
        ADuoL_ZEN + AJejL_ZEN + AIleL_ZEN + 
        ADuo_ZEN + AJej_ZEN + AIle_ZEN + 
        AL_ZEN + 
        AR_ZEN +
        AS_ZEN +
        AF_ZEN + 
        AB_ZEN + 
        AUR_ZEN + 
        AFE_ZEN + 
        
        A_bZEL +
        
        ADuoL_ZEN_gluc + AJejL_ZEN_gluc + AIleL_ZEN_gluc + 
        AJej_IC_ZEN_gluc + AJej_EC_ZEN_gluc + 
        AL_IC_ZEN_gluc + AL_EC_ZEN_gluc +
        ABile_ZEN_gluc + 
        AR_ZEN_gluc +
        AS_ZEN_gluc +
        AF_ZEN_gluc + 
        AB_ZEN_gluc + 
        AK_ZEN_gluc +
        AUR_ZEN_gluc + 
        AFE_ZEN_gluc +
        
        AIleL_aZEL + 
        AJej_aZEL + AIle_aZEL +
        AL_aZEL + 
        AR_aZEL +
        AS_aZEL +
        AF_aZEL +
        AB_aZEL +
        AUR_aZEL +
        AFE_aZEL +
        
        ADuoL_aZEL_gluc + AJejL_aZEL_gluc + AIleL_aZEL_gluc +
        AJej_IC_aZEL_gluc + AJej_EC_aZEL_gluc +
        AL_IC_aZEL_gluc + AL_EC_aZEL_gluc +
        ABile_aZEL_gluc + 
        AR_aZEL_gluc + 
        AS_aZEL_gluc +
        AF_aZEL_gluc +
        AB_aZEL_gluc +
        AK_aZEL_gluc +
        AUR_aZEL_gluc +
        AFE_aZEL_gluc 
      
      # MB <- ATotal - ACalculated + 1
      
      
      
      
      
      # End
      list(c(
        dAOD_ZEN,
        dADuoL_ZEN,dAJejL_ZEN, dAIleL_ZEN, 
        dADuo_ZEN, dAJej_ZEN, dAIle_ZEN,
        dAL_ZEN,
        dAR_ZEN,
        dAS_ZEN,
        dAF_ZEN,
        dAB_ZEN,
        dAUR_ZEN,
        dAFE_ZEN,
        
        dA_bZEL,
        
        dADuoL_ZEN_gluc, dAJejL_ZEN_gluc, dAIleL_ZEN_gluc,
        dAJej_IC_ZEN_gluc, dAJej_EC_ZEN_gluc,
        dAL_IC_ZEN_gluc, dAL_EC_ZEN_gluc ,
        dABile_ZEN_gluc,
        dAR_ZEN_gluc,
        dAS_ZEN_gluc,
        dAF_ZEN_gluc,
        dAB_ZEN_gluc,
        dAK_ZEN_gluc,
        dAUR_ZEN_gluc,
        dAFE_ZEN_gluc ,
        
        dAIleL_aZEL, 
        dAJej_aZEL, dAIle_aZEL,
        dAL_aZEL,
        dAR_aZEL,
        dAS_aZEL,
        dAF_aZEL,
        dAB_aZEL,
        dAUR_aZEL,
        dAFE_aZEL,
        
        dADuoL_aZEL_gluc, dAJejL_aZEL_gluc, dAIleL_aZEL_gluc,
        dAJej_IC_aZEL_gluc, dAJej_EC_aZEL_gluc,
        dAL_IC_aZEL_gluc, dAL_EC_aZEL_gluc,
        dABile_aZEL_gluc,
        dAR_aZEL_gluc,
        dAS_aZEL_gluc,
        dAF_aZEL_gluc,
        dAB_aZEL_gluc,
        dAK_aZEL_gluc,
        dAUR_aZEL_gluc,
        dAFE_aZEL_gluc #,
        
        # dATotal
      ),
      c(#ATotal = ATotal,
        ACalculated = ACalculated,
        #MB = MB,
        CDuoL_ZEN = CDuoL_ZEN, CJejL_ZEN = CJejL_ZEN, CIleL_ZEN = CIleL_ZEN, #CCOL_ZEN = CCOL_ZEN,
        CDuo_ZEN = CDuo_ZEN, CJej_ZEN = CJej_ZEN, CIle_ZEN = CIle_ZEN,
        CL_ZEN = CL_ZEN,
        CR_ZEN = CR_ZEN ,
        CS_ZEN = CS_ZEN,
        CF_ZEN = CF_ZEN,
        CB_ZEN = CB_ZEN,

        CDuoL_ZEN_gluc = CDuoL_ZEN_gluc, CJejL_ZEN_gluc = CJejL_ZEN_gluc, CIleL_ZEN_gluc = CIleL_ZEN_gluc, #CCOL_ZEN_gluc = CCOL_ZEN_gluc,
        CJej_ZEN_gluc = CJej_IC_ZEN_gluc, 
        CL_ZEN_gluc = CL_IC_ZEN_gluc, 
        CBile_ZEN_gluc = CBile_ZEN_gluc,
        CR_ZEN_gluc = CR_ZEN_gluc,
        CS_ZEN_gluc = CS_ZEN_gluc,
        CF_ZEN_gluc = CF_ZEN_gluc,
        CB_ZEN_gluc = CB_ZEN_gluc,
        CK_ZEN_gluc = CK_ZEN_gluc,

        CIleL_aZEL = CIleL_aZEL, 
        CJej_aZEL = CJej_aZEL, CIle_aZEL = CIle_aZEL,
        CL_aZEL = CL_aZEL,
        CR_aZEL = CR_aZEL,
        CS_aZEL = CS_aZEL,
        CF_aZEL = CF_aZEL,
        CB_aZEL = CB_aZEL,

        
        CDuoL_aZEL_gluc = CDuoL_aZEL_gluc, CJejL_aZEL_gluc = CJejL_aZEL_gluc, CIleL_aZEL_gluc = CIleL_aZEL_gluc, #CCOL_aZEL_gluc = CCOL_aZEL_gluc, #CJejL_aZEL_gluc = CJejL_aZEL_gluc
        CJej_aZEL_gluc = CJej_IC_aZEL_gluc, 
        CL_aZEL_gluc = CL_IC_aZEL_gluc, 
        CBile_aZEL_gluc = CBile_aZEL_gluc,
        CR_aZEL_gluc = CR_aZEL_gluc,
        CS_aZEL_gluc = CS_aZEL_gluc,
        CF_aZEL_gluc = CF_aZEL_gluc,
        CB_aZEL_gluc = CB_aZEL_gluc,
        CK_aZEL_gluc = CK_aZEL_gluc

      ))
      
    })
  }
  
  output_PBK <-lsoda(y = y, 
                     times = times, 
                     func = PBK.model,
                     parms = parms ,
                     events = list(data = Dosingevents)
                     )
  
  return(as.data.frame(output_PBK))
}
  
