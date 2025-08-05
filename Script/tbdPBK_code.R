## PBK model rat ####
PBK.model <- function(t, state, parms){
  with(as.list(c(state, parms)), {
    
    ## Calculated parameters ####
    
    ## Physiological constants ####
    # Fractional Tissue volumes (reference: Punt et al. (2020))
    BW <- BW  ##{Kg} #  body weight rat
    
    VST <- FV_ST*BW  # fraction volume of stomach
    VSIL <- FV_SIL*BW # fractional volume of small intestinal lumen
    VDuoL <- VSIL*0.183 # calculated from Musther
    VJejL <- VSIL*0.334 # calculated from Musther
    VIleL <- VSIL*0.483 # calculated from Musther
    
    VREL <- FV_REL*FV_LIL*BW # fractional volume of rectum lumen from total LI volume (only the lower part goes to the systemic circulation)
    VLIL <- (FV_LIL*BW)-VREL # corrected fractional volume of the rest of the colon
    
    VSI <- FV_SI*BW # fractional volume of small intestinal tissue
    VDuo <- VSI*0.060
    VJej <- VSI*0.769
    VIle <- VSI*0.171
    
    VL <- FV_L*BW  # fractional of liver tissue
    VR <- FV_R*BW  # 0.082-VSIc-VLc # fraction of rapidly perfused tissue
    VS <- FV_S*BW  # 0.737-VFc # fraction of slowly perfused tissue
    VF <- FV_F*BW  # fraction of fat tissue
    VB <- FV_B*BW  # fraction of blood
    
    # Surface area 
    SASI <- SASI # surface area (dm2)
    SARE <- FV_REL*SALI
    SALI <- SALI- SARE  # surface area of large intestine (dm2) reference: Vdoviaková et al. (2016)
    
    # Blood flow rates (reference: Punt et al. (2020))
    QC <- QC  ##{L/h} # cardiac output
    QSI <- QSIc*QC  # fraction of blood flow to small intestine
    QDuo <- QSI*0.079
    QJej <- QSI*0.680
    QIle <- QSI*0.241
    QL <- QLc*QC  # 0.242 - QSIc # fraction of blood flow to liver 101 
    QR <- QRc*QC  # 0.76 - QSIc - QLc # fraction of blood flow to rapidly perfused tissue
    QS <- QSc*QC  # 0.24 - QFc # fraction of blood flow to slowly perfused tissue
    QF <- QFc*QC  # fraction of blood flow to fat
    
    # Intestinal fluid flow rates (reference: Punt et al. (2020))
    GER <- GER # stomach emptying rate (/h)
    SITR <- SITR #/10 # small intestinal transit time (/h)
    LITR <- LITR #/10 # large intestinal transit time (/h)
    RETR <- RETR #/10
    
    # Glomerular filtration rate: Punt et al. (2020)
    GFR <- GFRc*BW #(/h)
    
    ## kinetic parameters ####
    ### absorption/transfer rates ####
    R <- 0.18 #radius small intestine
    
    Rduo <- 0.275/2 #diameter / 2 = radius
    Rjej <- 0.4/2
    Rile <- 0.4/2
    
    Tsi <- 1.63 #small intestine transit time 
    
    Tsiduo <- 0.06*Tsi
    Tsijej <- 0.192*Tsi
    Tsiile <- 0.086*Tsi
    
    kt <- 1/(Tsi/7)
    Peff.cm.per.hr <- PappCaco2_ZEN #cm/h
    ka <- Peff.cm.per.hr*2/R
    fa_duo <- 1 - (1 + ka/kt)**-7
    fa_jej <- 1 - (1 + ka/kt)**-7
    fa_ile <- 1 - (1 + ka/kt)**-7
    
    CLSI_ZEN <- PappCaco2_ZEN * SASI # (L/h) cm/h to L/h reference: Pfeiffer et al. (2011)
    CLLI_ZEN <- 0 #PappCaco2_ZEN * SALI #* 1e-3 # (L/h) cm/h to L/h reference: Pfeiffer et al. (2011)
    CLRE_ZEN <- 0 #PappCaco2_ZEN * SARE #* 1e-3 # (L/h) cm/h to L/h reference: Pfeiffer et al. (2011)
    CLLI_aZEL <- 0 #PappCaco2aZEL * SALI #* 1e-3 # (L/h) cm/h to L/h reference: Pfeiffer et al. (2011)
    
    ### metabolism in small intestine tissue ####
    VmaxSI_aZEL <- VmaxSI_aZELc * S9SI * VSI        # umol/h
    VmaxSI_bZEL <- VmaxSI_bZELc * S9SI * VSI        # umol/h
    VmaxSI_ZENGLU <- VmaxSI_ZENGLUc * S9SI * VSI  # umol/h
    VmaxSI_aZELGLU <- VmaxSI_aZELGLUc * S9SI * VSI # umol/h
    
    ### metabolism in large intestine lumen (microbiota compartment) ####
    VmaxLI_aZEL <- VmaxLI_aZELc*VMB  # umol/h 
    VmaxLI_bZEL <- VmaxLI_bZELc*VMB
    
    ### metabolism in liver ####
    VmaxL_aZEL <- VmaxL_aZELc * S9L * VL 
    VmaxL_bZEL <- VmaxL_bZELc * S9L * VL
    VmaxL_ZENGLU <- VmaxL_ZENGLUc * S9L * VL  
    VmaxL_aZELGLU <- VmaxL_aZELGLUc * S9L * VL
    
    ### conversion ZENGLU to ZEN
    CL_bile <- 1 #1e6 #2 #000
    CL_bile_ZEN <- 0 #1
    kHydro <- 1 #e-6 #1000 #7.440476e-05
    CL_SIZENGLU <- 0 #10
    
    ## Concentrations #####
    
    # Zearalenone 
    CSI_ZEN <- ASI_ZEN/VSI
    CSIL_ZEN <- ASIL_ZEN/VSI
    CLIL_ZEN <- ALIL_ZEN/VLIL 
    CREL_ZEN <- AREL_ZEN/VREL
    CRE_ZEN <- ARE_ZEN/VREL
    CL_ZEN <- AL_ZEN/VL
    CF_ZEN <- AF_ZEN/VF
    CR_ZEN <- AR_ZEN/VR
    CS_ZEN <- AS_ZEN/VS
    CB_ZEN <- AB_ZEN/VB
    
    CVSI_ZEN <- CSI_ZEN/PSI_ZEN
    CVRE_ZEN <- CRE_ZEN/PSI_ZEN
    CVL_ZEN <- CL_ZEN/PL_ZEN
    CVF_ZEN <- CF_ZEN/PF_ZEN
    CVR_ZEN <- CR_ZEN/PR_ZEN
    CVS_ZEN <- CS_ZEN/PS_ZEN
    
    
    # α-Zearalenol
    CLIL_aZEL <- ALIL_aZEL/VLIL
    CSI_aZEL <- ASI_aZEL/VSI
    CL_aZEL <- AL_aZEL/VL
    CF_aZEL <- AF_aZEL/VF
    CR_aZEL <- AR_aZEL/VR
    CS_aZEL <- AS_aZEL/VS
    CB_aZEL <- AB_aZEL/VB
    
    CVSI_aZEL <- CSI_aZEL/PIaZEL
    CVL_aZEL <- CL_aZEL/PLaZEL
    CVF_aZEL <- CF_aZEL/PFaZEL
    CVR_aZEL <- CR_aZEL/PRaZEL
    CVS_aZEL <- CS_aZEL/PSaZEL
    
    
    # β-Zearalenol
    CLIL_bZEL <- ALIL_bZEL/VLIL
    CL_bZEL <- AL_bZEL/VL
    
    
    # Zearalenone-glucuronide
    # CLIL_ZENGLU <- ALIL_ZENGLU/VLIL 
    CSI_ZENGLU <- ASI_ZENGLU/VSI 
    CSIL_ZENGLU <- ASIL_ZENGLU/VSIL
    CLIL_ZENGLU <- ALIL_ZENGLU/VLIL 
    CL_ZENGLU <- AL_ZENGLU/VL
    CF_ZENGLU <- AF_ZENGLU/VF
    CR_ZENGLU <- AR_ZENGLU/VR
    CS_ZENGLU <- AS_ZENGLU/VS
    CB_ZENGLU <- AB_ZENGLU/VB
    
    CVSI_ZENGLU <- CSI_ZENGLU/PSI_ZENGLU
    CVL_ZENGLU <- CL_ZENGLU/PL_ZENGLU
    CVF_ZENGLU <- CF_ZENGLU/PF_ZENGLU
    CVR_ZENGLU <- CR_ZENGLU/PR_ZENGLU
    CVS_ZENGLU <- CS_ZENGLU/PS_ZENGLU
    
    # a-Zearalenol glucuronide
    CL_aZELGLU <- AL_aZELGLU/VL   # liver a zel gluc
    CSI_aZELGLU <- ASI_aZELGLU/VSI
    
    
    
    ## Dynamic ####
    
    # Stomach (in) 
    dAST_ZEN <- - GER*AST_ZEN 
    
    # Small intestinal lumen 
    dADuo_ZEN <- + GER*AST_ZEN - Tsiduo*ADuo_ZEN*(1-fa_duo) - ka*ADuo_ZEN*fa_duo 
    
    dAJej_ZEN <- Tsiduo*ADuo_ZEN*(1-fa_duo) - ka*AJej_ZEN*fa_jej - Tsijej*AJej_ZEN*(1-fa_jej) 
    
    dAIle_ZEN <- Tsijej*AJej_ZEN*(1-fa_jej) - ka*AIle_ZEN*fa_ile - Tsiile*AIle_ZEN*(1-fa_ile) 
    
    dAJej_ZENGLU <- + CL_bile*CL_ZENGLU - SITR*ASIL_ZENGLU + CL_SIZENGLU*CSI_ZENGLU*fup_ZENGLU 
    
    
    # Large intestinal lumen 
    dALIL_ZEN <- Tsiile*AIle_ZEN*(1-fa_ile) +
      - VmaxLI_aZEL*CLIL_ZEN / (KmLI_aZEL + CLIL_ZEN) +
      - VmaxLI_bZEL*CLIL_ZEN / (KmLI_bZEL + CLIL_ZEN) +
      + kHydro*ALIL_ZENGLU + 
      - CLLI_ZEN*CLIL_ZEN +
      - LITR*ALIL_ZEN
    
    dALIL_ZENGLU <- + SITR*ASIL_ZENGLU - kHydro*ALIL_ZENGLU - LITR*ALIL_ZENGLU 
    
    dALIL_aZEL <- VmaxLI_aZEL*CLIL_ZEN / (KmLI_aZEL + CLIL_ZEN) +
      - (LITR+RETR)*ALIL_aZEL 
    
    dALIL_bZEL <- VmaxLI_bZEL*CLIL_ZEN / (KmLI_bZEL + CLIL_ZEN)  
    
    # Rectum lumen
    dAREL_ZEN <- LITR*ALIL_ZEN - RETR*AREL_ZEN - CLRE_ZEN*CREL_ZEN 
    dARE_ZEN <- CLRE_ZEN*CREL_ZEN + QSI*(CB_ZEN-CVRE_ZEN)
    
    # Faeces ####
    dAFE_ZEN <- RETR*AREL_ZEN 
    dAFE_aZEL <- (LITR+RETR)*ALIL_aZEL
    dAFE_ZENGLU <- LITR*ALIL_ZENGLU 
    
    
    # Small intestinal tissue ####
    dADuo_ZEN <- ka*ADuo_ZEN*fa_duo + QDuo*(CB_ZEN-CVDuo_ZEN) 
    dAJej_ZEN <- ka*AJej_ZEN*fa_jej + QJej*(CB_ZEN-CVJej_ZEN) +
      - VmaxSI_ZENGLU*CJej_ZEN*fup_ZEN / (KmSI_ZENGLU + CJej_ZEN*fup_ZEN) + 
      - VmaxSI_aZEL*CJej_ZEN*fup_ZEN / (KmSI_aZEL + CJej_ZEN*fup_ZEN) +
      - VmaxSI_bZEL*CJej_ZEN*fup_ZEN / (KmSI_bZEL + CJej_ZEN*fup_ZEN)
    
    dAIle_ZEN <- ka*AIle_ZEN*fa_ile + QIle*(CB_ZEN-CVDuo_ZEN) 
    
    dASI_ZENGLU <- VmaxSI_ZENGLU*CJej_ZEN*fup_ZEN / (KmSI_ZENGLU + CJej_ZEN*fup_ZEN) +
      + QSI*(CB_ZENGLU-CVSI_ZENGLU) - CL_SIZENGLU*CSI_ZENGLU*fup_ZENGLU
    
    dASI_aZEL <- VmaxSI_aZEL*CJej_ZEN*fup_ZEN / (KmSI_aZEL + CJej_ZEN*fup_ZEN) +
      + QSI*(CB_aZEL-CVSI_aZEL) +
      - VmaxSI_aZELGLU*CSI_aZEL*fup_aZEL / (KmSI_aZELGLU + CSI_aZEL*fup_aZEL)
    
    dASI_bZEL <- VmaxSI_bZEL*CJej_ZEN*fup_ZEN / (KmSI_bZEL + CJej_ZEN*fup_ZEN)
    
    dASI_aZELGLU <- VmaxSI_aZELGLU*CSI_aZEL*fup_aZEL / (KmSI_aZELGLU + CSI_aZEL*fup_aZEL)
    
    
    # Liver 
    dAL_ZEN <- + QDuo*CVDuo_ZEN + QJej*CVJej_ZEN + QIle*CVDuo_ZEN + QL*CB_ZEN +
      - (QSI+QL)*CVL_ZEN + 
      - VmaxL_aZEL*CL_ZEN*fup_ZEN / (KmL_aZEL + CL_ZEN*fup_ZEN) +
      - VmaxL_bZEL*CL_ZEN*fup_ZEN / (KmL_bZEL + CL_ZEN*fup_ZEN) +
      - VmaxL_ZENGLU*CL_ZEN*fup_ZEN / (KmL_ZENGLU + CL_ZEN*fup_ZEN)
    
    dAL_ZENGLU <- VmaxL_ZENGLU*CL_ZEN*fup_ZEN / (KmL_ZENGLU + CL_ZEN*fup_ZEN) +
      +  QSI*CVSI_ZENGLU + QL*CB_ZENGLU - (QSI+QL)*CVL_ZENGLU +
      - CL_bile*CL_ZENGLU 
    
    dAL_bZEL <- VmaxL_bZEL*CL_ZEN*fup_ZEN / (KmL_bZEL + CL_ZEN*fup_ZEN) 
    
    dAL_aZEL <- VmaxL_aZEL*CL_ZEN*fup_ZEN / (KmL_aZEL + CL_ZEN*fup_ZEN) +
      + QSI*CVSI_aZEL + QL*CB_aZEL - (QSI+QL)*CVL_aZEL +
      - VmaxL_aZELGLU*CL_aZEL*fup_aZEL / (KmL_aZELGLU + CL_aZEL*fup_aZEL) 
    
    dAL_aZELGLU <- VmaxL_aZELGLU*CL_aZEL*fup_aZEL / (KmL_aZELGLU + CL_aZEL*fup_aZEL) 
    
    # Blood
    dAB_ZEN <- - (QSI + QL + QR + QS + QF)*CB_ZEN +  
      + (QSI+QL)*CVL_ZEN + QF*CVF_ZEN + QR*CVR_ZEN + QS*CVS_ZEN + 
      - GFR*CB_ZEN*fup_ZEN +
      + CLRE_ZEN*CREL_ZEN 

    dAB_aZEL <- - (QSI + QL + QR + QS + QF)*CB_aZEL +
      + (QL+QSI)*CVL_aZEL + QF*CVF_aZEL + QR*CVR_aZEL + QS*CVS_aZEL +
      - GFR*CB_aZEL*fup_aZEL
    
    dAB_ZENGLU <- - (QSI + QL + QR + QS + QF)*CB_ZENGLU + 
      + (QSI+QL)*CVL_ZENGLU + QF*CVF_ZENGLU + QR*CVR_ZENGLU + QS*CVS_ZENGLU +
      - GFR*CB_ZENGLU*fup_ZENGLU 
    
    # Urine
    dAUR_ZEN <- GFR*CB_ZEN*fup_ZEN 
    dAUR_aZEL <- GFR*CB_aZEL*fup_aZEL
    dAUR_ZENGLU <- GFR*CB_ZENGLU*fup_ZENGLU
    
    # Fat 
    dAF_ZEN <- QF*(CB_ZEN-CVF_ZEN) 
    dAF_aZEL <- QF*(CB_aZEL-CVF_aZEL) 
    dAF_ZENGLU <- QF*(CB_ZENGLU-CVF_ZENGLU) 
    
    # Rapidly Perfused 
    dAR_ZEN <- QR*(CB_ZEN-CVR_ZEN)
    dAR_aZEL <- QR*(CB_aZEL-CVR_aZEL) 
    dAR_ZENGLU <- QR*(CB_ZENGLU-CVR_ZENGLU) 
    
    # Slowly Perfused 
    dAS_ZEN <- QS*(CB_ZEN-CVS_ZEN)
    dAS_aZEL <- QS*(CB_aZEL-CVS_aZEL) 
    dAS_ZENGLU <- QS*(CB_ZENGLU-CVS_ZENGLU)
    
    #mass balance calculation
    ATotal <- ODOSEumol + IVDOSEumol
    ACalculated <- AST_ZEN + ASIL_ZEN + ASIL_ZENGLU + 
      + ALIL_ZEN + ALIL_aZEL + ALIL_bZEL + ALIL_ZENGLU +
      + AREL_ZEN + ARE_ZEN +
      + AFE_ZEN + AFE_aZEL + AFE_ZENGLU +
      + ASI_ZEN + ASI_ZENGLU + ASI_bZEL + ASI_aZEL + ASI_aZELGLU +
      + AL_ZEN + AL_ZENGLU + AL_bZEL + AL_aZEL + AL_aZELGLU +
      + AB_ZEN + AB_aZEL + AB_ZENGLU + 
      + AUR_ZEN + AUR_aZEL + AUR_ZENGLU +
      + AF_ZEN + AF_aZEL + AF_ZENGLU +
      + AR_ZEN + AR_aZEL + AR_ZENGLU + 
      + AS_ZEN + AS_aZEL + AS_ZENGLU
    
    
    MB <- ATotal - ACalculated + 1
    
    # End
    # 
    list(c(
      dAST_ZEN,
      dASIL_ZEN, dASIL_ZENGLU,  
      dALIL_ZEN, dALIL_ZENGLU, dALIL_aZEL, dALIL_bZEL, 
      dAREL_ZEN, dARE_ZEN,
      dAFE_ZEN, dAFE_aZEL, dAFE_ZENGLU,
      dASI_ZEN, dASI_ZENGLU, dASI_bZEL, dASI_aZEL, dASI_aZELGLU,
      dAL_ZEN, dAL_ZENGLU, dAL_bZEL, dAL_aZEL, dAL_aZELGLU, 
      dAB_ZEN, dAB_aZEL, dAB_ZENGLU, 
      dAUR_ZEN, dAUR_aZEL, dAUR_ZENGLU,
      dAF_ZEN, dAF_aZEL, dAF_ZENGLU,
      dAR_ZEN, dAR_aZEL, dAR_ZENGLU,  
      dAS_ZEN, dAS_aZEL, dAS_ZENGLU#,  
      # dATotal
    ),
    c(ATotal = ATotal,
      ACalculated = ACalculated,
      MB = MB,
      CSI_ZEN = CSI_ZEN, CSI_ZENGLU = CSI_ZENGLU,  CSI_aZEL = CSI_aZEL, CSI_aZELGLU = CSI_aZELGLU,
      CSIL_ZEN = CSIL_ZEN, CSIL_ZENGLU = CSIL_ZENGLU,
      CLIL_ZEN = CLIL_ZEN, CLIL_ZENGLU = CLIL_ZENGLU, CLIL_aZEL = CLIL_aZEL, CLIL_bZEL = CLIL_bZEL,
      CREL_ZEN = CREL_ZEN,
      CL_ZEN = CL_ZEN, CL_ZENGLU = CL_ZENGLU, CL_bZEL = CL_bZEL, CL_aZEL = CL_aZEL, CL_aZELGLU = CL_aZELGLU,
      CF_ZEN = CF_ZEN, CF_aZEL = CF_aZEL, CF_ZENGLU = CF_ZENGLU,
      CR_ZEN = CR_ZEN, CR_aZEL = CR_aZEL, CR_ZENGLU = CR_ZENGLU, 
      CS_ZEN = CS_ZEN, CS_aZEL = CS_aZEL, CS_ZENGLU = CS_ZENGLU,
      CB_ZEN = CB_ZEN, CB_aZEL = CB_aZEL, CB_ZENGLU = CB_ZENGLU
    ))

  })
}