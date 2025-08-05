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
    VLIL <- FV_LIL*BW #(FV_LIL*BW)-VREL # corrected fractional volume of the rest of the colon
    
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
    
    # Intestinal fluid flow rates (reference: Musther et al SI; Punt et al. (2020))
    GER <- GER # stomach emptying rate (/h)
    SITR <- SITR # small intestinal transit time (/h)
    LITR <- LITR # large intestinal transit time (/h)
    RETR <- RETR #
    
    # Glomerular filtration rate: Punt et al. (2020)
    GFR <- GFRc*BW # (/h)
    
    ## kinetic parameters ####
    ### absorption/transfer rates ####
    R <- 0.18 #radius small intestine
    
    Rduo <- 0.275/2 #diameter / 2 = radius, cm (Musther et al)
    Rjej <- 0.4/2   #diameter / 2 = radius, cm (Musther et al)
    Rile <- 0.4/2   #diameter / 2 = radius, cm (Musther et al)
    # Rlil <- 0.928/2 #diameter / 2 = radius, cm (Bartoli et al. 2013: DOI: 10.4253/wjge.v5.i5.226, diameter of the proximal colon)
    
    Tsi <- 1.63 #small intestine transit time 
    
    Tsiduo <- 0.06*Tsi
    Tsijej <- 0.192*Tsi
    Tsiile <- 0.086*Tsi
    
    duoTR <- 1/Tsiduo
    jejTR <- 1/Tsijej
    ileTR <- 1/Tsiile
    # lilTR <- LITR
    
    # kt <- 1/(Tsi/7)
    Peff.cm.per.hr <- 5*PappCaco2_ZEN #cm/h
    ka_duo <- Peff.cm.per.hr/Rduo
    ka_jej <- Peff.cm.per.hr/Rjej
    ka_ile <- Peff.cm.per.hr/Rile
    # ka_lil <- Peff.cm.per.hr*2/Rlil
    
    fa_duo <- 1 - (1 + ka_duo/(1/Tsiduo))**-7
    fa_jej <- 1 - (1 + ka_jej/(1/Tsijej))**-7
    fa_ile <- 1 - (1 + ka_ile/(1/Tsiile))**-7
    # fa_lil <- 1 - (1 + ka/(1/Tsilil))**-7
    
    CLSI_ZEN <- 0 #PappCaco2_ZEN * SASI # (L/h) cm/h to L/h reference: Pfeiffer et al. (2011) # uptake ZEN SI
    CLLI_ZEN <- 0 #PappCaco2_ZEN * SALI #* 1e-3 # (L/h) cm/h to L/h reference: Pfeiffer et al. (2011) # uptake ZEN LI 
    CLRE_ZEN <- 0 #PappCaco2_ZEN * SARE #* 1e-3 # (L/h) cm/h to L/h reference: Pfeiffer et al. (2011) # uptake ZEN rectum
    CLLI_aZEL <- 0 #PappCaco2aZEL * SALI #* 1e-3 # (L/h) cm/h to L/h reference: Pfeiffer et al. (2011) # uptake aZEL LI
    CL_SIZENGLU <- 0 #10 # clearance ZENglu from intestinal tissue to intestinal lumen
    
    
    ### metabolism in small intestine tissue ####
    VmaxSI_aZEL <- VmaxSI_aZELc * S9SI * VJej        # umol/h
    VmaxSI_bZEL <- VmaxSI_bZELc * S9SI * VJej        # umol/h
    VmaxSI_ZENGLU <- 100*VmaxSI_ZENGLUc * S9SI * VJej  # umol/h
    VmaxSI_aZELGLU <- VmaxSI_aZELGLUc * S9SI * VJej # umol/h
    
    ### metabolism in large intestine lumen (microbiota compartment) ####
    VmaxLI_aZEL <- VmaxLI_aZELc*VMB  # umol/h 
    VmaxLI_bZEL <- VmaxLI_bZELc*VMB
    
    ### metabolism in liver ####
    VmaxL_aZEL <- VmaxL_aZELc * S9L * VL 
    VmaxL_bZEL <- VmaxL_bZELc * S9L * VL
    VmaxL_ZENGLU <- 10*VmaxL_ZENGLUc * S9L * VL  
    VmaxL_aZELGLU <- VmaxL_aZELGLUc * S9L * VL
    
    ### biliary clearance
    fBile <- 0.5
    Vmax_BSEPc <- 7                            # umol/min/mg BSEP, Average active transport of bile acids, assuming that the maximum velocity of PFOA transport by BSEP corresponds to that of bile acids ; de Bruijn et al. 2024 https://doi.org/10.14573/altex.2302011 
    Km_BSEPc <- 16.4                           # umol/L, uM, Average affinity constant of bile acids to BSEP, following the above assumption ; de Bruijn et al. 2024 https://doi.org/10.14573/altex.2302011
    
    # IVIVE 
    mgBSEP.HC <-  0.839 * 140000 * 1e-9         # mg BSEP/1e6 hepatocytes (0.839 pmole BSEP/1e6 hepatocytes (amound of BSEP per hepatocyte) * 140000 g/mole (MW BSEP) -> pg * 1e-9 -> mg) de Bruijn et al. 2024 https://doi.org/10.14573/altex.2302011]
    HC.GL <- 99                                 # 1e6 hepatocytes in the liver ; de Bruijn et al. 2024 https://doi.org/10.14573/altex.2302011
    SF_BSEP <- mgBSEP.HC * HC.GL * 1e3          # Scaling factor for BSEP mediated hepatic efflux for GCA and GCDC ; de Bruijn et al. 2024 https://doi.org/10.14573/altex.2302011
    
    Vmax_BSEP <- Vmax_BSEPc*60*24*SF_BSEP*VL         # umol/d
    Km_BSEP <- 10*Km_BSEPc                             # umol/L (uM)
    
    ### conversion ZENGLU to ZEN
    # CL_bile <- 7.440476e-02  # biliary clearance ZEN gluc 
    CL_Hydro <- 1 # metabolic conversion ZEN gluc to ZEN
    
    # fup_ZENGLU <- 0.15
    ## Concentrations #####
    
    # Zearalenone 
    # CSI_ZEN <- ASI_ZEN/VSI
    CDuo_ZEN <- ADuo_ZEN/VDuo
    CJej_ZEN <- AJej_ZEN/VJej
    CIle_ZEN <- AIle_ZEN/VIle
    # CSIL_ZEN <- ASIL_ZEN/VSI
    CDuoL_ZEN <- ADuoL_ZEN/VDuoL
    CJejL_ZEN <- AJejL_ZEN/VJejL
    CIleL_ZEN <- AIleL_ZEN/VIleL
    CLIL_ZEN <- ALIL_ZEN/VLIL 
    CREL_ZEN <- AREL_ZEN/VREL
    # CRE_ZEN <- ARE_ZEN/VREL
    CL_ZEN <- AL_ZEN/VL
    CF_ZEN <- AF_ZEN/VF
    CR_ZEN <- AR_ZEN/VR
    CS_ZEN <- AS_ZEN/VS
    CB_ZEN <- AB_ZEN/VB
    
    # CVSI_ZEN <- CSI_ZEN/PSI_ZEN
    CVDuo_ZEN <- CDuo_ZEN/PSI_ZEN
    CVJej_ZEN <- CJej_ZEN/PSI_ZEN
    CVIle_ZEN <- CIle_ZEN/PSI_ZEN
    # CVRE_ZEN <- CRE_ZEN/PSI_ZEN
    CVL_ZEN <- CL_ZEN/PL_ZEN
    CVF_ZEN <- CF_ZEN/PF_ZEN
    CVR_ZEN <- CR_ZEN/PR_ZEN
    CVS_ZEN <- CS_ZEN/PS_ZEN
    
    
    # α-Zearalenol
    CLIL_aZEL <- ALIL_aZEL/VLIL
    # CSI_aZEL <- ASI_aZEL/VSI
    CJej_aZEL <- AJej_aZEL/VJej
    CL_aZEL <- AL_aZEL/VL
    CF_aZEL <- AF_aZEL/VF
    CR_aZEL <- AR_aZEL/VR
    CS_aZEL <- AS_aZEL/VS
    CB_aZEL <- AB_aZEL/VB
    
    # CVSI_aZEL <- CSI_aZEL/PIaZEL
    CVJej_aZEL <- CJej_aZEL/PIaZEL
    CVL_aZEL <- CL_aZEL/PLaZEL
    CVF_aZEL <- CF_aZEL/PFaZEL
    CVR_aZEL <- CR_aZEL/PRaZEL
    CVS_aZEL <- CS_aZEL/PSaZEL
    
    
    # β-Zearalenol
    CLIL_bZEL <- ALIL_bZEL/VLIL
    CL_bZEL <- AL_bZEL/VL
    
    
    # Zearalenone-glucuronide
    # CSI_ZENGLU <- ASI_ZENGLU/VSI 
    CJej_ZENGLU <- AJej_ZENGLU/VJej
    CDuoL_ZENGLU <- ADuoL_ZENGLU/VDuoL
    CJejL_ZENGLU <- AJejL_ZENGLU/VJejL
    CIleL_ZENGLU <- AIleL_ZENGLU/VIleL 
    # CSIL_ZENGLU <- ASIL_ZENGLU/VSIL
    CLIL_ZENGLU <- ALIL_ZENGLU/VLIL 
    CL_ZENGLU <- AL_ZENGLU/VL
    CF_ZENGLU <- AF_ZENGLU/VF
    CR_ZENGLU <- AR_ZENGLU/VR
    CS_ZENGLU <- AS_ZENGLU/VS
    CB_ZENGLU <- AB_ZENGLU/VB
    
    # CVSI_ZENGLU <- CSI_ZENGLU/PSI_ZENGLU
    CVJej_ZENGLU <- CJej_ZENGLU/PSI_ZENGLU
    CVL_ZENGLU <- CL_ZENGLU/PL_ZENGLU
    CVF_ZENGLU <- CF_ZENGLU/PF_ZENGLU
    CVR_ZENGLU <- CR_ZENGLU/PR_ZENGLU
    CVS_ZENGLU <- CS_ZENGLU/PS_ZENGLU
    
    # a-Zearalenol glucuronide
    CL_aZELGLU <- AL_aZELGLU/VL   # liver a zel gluc
    # CSI_aZELGLU <- ASI_aZELGLU/VSI
    CJej_aZELGLU <- AJej_aZELGLU/VJej
    
    
    ## Dynamic ####
    
    # Stomach (in) 
    dAST_ZEN <- - GER*AST_ZEN 
    
    # Small intestinal lumen 
    dADuoL_ZEN <- + GER*AST_ZEN - duoTR*ADuoL_ZEN*(1-fa_duo) + 
      - ka_duo*ADuoL_ZEN*fa_duo
    dAJejL_ZEN <- duoTR*ADuoL_ZEN*(1-fa_duo) - jejTR*AJejL_ZEN*(1-fa_jej) +
      - ka_jej*AJejL_ZEN*fa_jej
    dAIleL_ZEN <- jejTR*AJejL_ZEN*(1-fa_jej) - ileTR*AIleL_ZEN*(1-fa_ile) + 
      - ka_ile*AIleL_ZEN*fa_ile  +
      + jejTR*(CL_Hydro*CJejL_ZENGLU) +
      + CL_Hydro*CIleL_ZENGLU
    
    dADuoL_ZENGLU <- fBile*ABile_ZENGLU + #(Vmax_BSEP/(Km_BSEP + (CL_ZENGLU*fup_ZENGLU)))*CL_ZENGLU*fup_ZENGLU +
      - duoTR*ADuoL_ZENGLU 
    dAJejL_ZENGLU <- duoTR*ADuoL_ZENGLU - jejTR*AJejL_ZENGLU + 
      - jejTR*(CL_Hydro*CJejL_ZENGLU) +
      + CL_SIZENGLU*CJej_ZENGLU*fup_ZENGLU 
    dAIleL_ZENGLU <- jejTR*AJejL_ZENGLU - ileTR*AIleL_ZENGLU +
      - CL_Hydro*CIleL_ZENGLU
    
    # Large intestinal lumen 
    dALIL_ZEN <- ileTR*AIleL_ZEN*(1-fa_ile) +
      - VmaxLI_aZEL*CLIL_ZEN / (KmLI_aZEL + CLIL_ZEN) +
      - VmaxLI_bZEL*CLIL_ZEN / (KmLI_bZEL + CLIL_ZEN) +
      + CL_Hydro*CLIL_ZENGLU + 
      # - CLLI_ZEN*CLIL_ZEN +
      - LITR*ALIL_ZEN
    
    dALIL_ZENGLU <- + ileTR*AIleL_ZENGLU - LITR*ALIL_ZENGLU +
      - CL_Hydro*CLIL_ZENGLU  
    
    dALIL_aZEL <- VmaxLI_aZEL*CLIL_ZEN / (KmLI_aZEL + CLIL_ZEN) +
      - LITR*ALIL_aZEL 
    
    dALIL_bZEL <- VmaxLI_bZEL*CLIL_ZEN / (KmLI_bZEL + CLIL_ZEN)  
    
    # Rectum lumen
    dAREL_ZEN <- 0 #LITR*ALIL_ZEN - RETR*AREL_ZEN - CLRE_ZEN*CREL_ZEN 
    
    # Faeces ####
    dAFE_ZEN <- LITR*ALIL_ZEN #RETR*AREL_ZEN 
    dAFE_aZEL <- LITR*ALIL_aZEL
    dAFE_ZENGLU <- LITR*ALIL_ZENGLU 
    
    
    # Small intestinal tissue ####
    dADuo_ZEN <- ka_duo*ADuoL_ZEN*fa_duo + QDuo*(CB_ZEN-CVDuo_ZEN) 
    dAJej_ZEN <- ka_jej*AJejL_ZEN*fa_jej + QJej*(CB_ZEN-CVJej_ZEN) +
      - VmaxSI_ZENGLU*CJej_ZEN*fup_ZEN / ((KmSI_ZENGLU) + CJej_ZEN*fup_ZEN) + 
      - VmaxSI_aZEL*CJej_ZEN*fup_ZEN / (KmSI_aZEL + CJej_ZEN*fup_ZEN) +
      - VmaxSI_bZEL*CJej_ZEN*fup_ZEN / (KmSI_bZEL + CJej_ZEN*fup_ZEN)
    
    dAIle_ZEN <- ka_ile*AIleL_ZEN*fa_ile + QIle*(CB_ZEN-CVIle_ZEN) 
    
    dAJej_ZENGLU <- VmaxSI_ZENGLU*CJej_ZEN*fup_ZEN / ((KmSI_ZENGLU) + CJej_ZEN*fup_ZEN) + 
      + QSI*(CB_ZENGLU-CVJej_ZENGLU) +
      - CL_SIZENGLU*CJej_ZENGLU*fup_ZENGLU
    
    dAJej_aZEL <- VmaxSI_aZEL*CJej_ZEN*fup_ZEN / (KmSI_aZEL + CJej_ZEN*fup_ZEN) +
      + QSI*(CB_aZEL-CVJej_aZEL) +
      - VmaxSI_aZELGLU*CJej_aZEL*fup_aZEL / (KmSI_aZELGLU + CJej_aZEL*fup_aZEL)
    
    dAJej_bZEL <- VmaxSI_bZEL*CJej_ZEN*fup_ZEN / (KmSI_bZEL + CJej_ZEN*fup_ZEN)
    
    dAJej_aZELGLU <- VmaxSI_aZELGLU*CJej_aZEL*fup_aZEL / (KmSI_aZELGLU + CJej_aZEL*fup_aZEL)
    
    
    # Liver 
    dAL_ZEN <- + QDuo*CVDuo_ZEN + QJej*CVJej_ZEN + QIle*CVIle_ZEN + QL*CB_ZEN +
      - (QSI+QL)*CVL_ZEN + 
      - VmaxL_aZEL*CL_ZEN*fup_ZEN / (KmL_aZEL + CL_ZEN*fup_ZEN) +
      - VmaxL_bZEL*CL_ZEN*fup_ZEN / (KmL_bZEL + CL_ZEN*fup_ZEN) +
      - VmaxL_ZENGLU*CL_ZEN*fup_ZEN / ((KmL_ZENGLU) + CL_ZEN*fup_ZEN)
    
    dAL_ZENGLU <- VmaxL_ZENGLU*CL_ZEN*fup_ZEN / ((KmL_ZENGLU) + CL_ZEN*fup_ZEN) +
      +  QSI*CVJej_ZENGLU + QL*CB_ZENGLU - (QSI+QL)*CVL_ZENGLU +
      - (Vmax_BSEP/(Km_BSEP + (CL_ZENGLU*fup_ZENGLU)))*CL_ZENGLU*fup_ZENGLU
    
    dAL_bZEL <- VmaxL_bZEL*CL_ZEN*fup_ZEN / (KmL_bZEL + CL_ZEN*fup_ZEN) 
    
    dAL_aZEL <- VmaxL_aZEL*CL_ZEN*fup_ZEN / (KmL_aZEL + CL_ZEN*fup_ZEN) +
      + QSI*CVJej_aZEL + QL*CB_aZEL - (QSI+QL)*CVL_aZEL +
      - VmaxL_aZELGLU*CL_aZEL*fup_aZEL / (KmL_aZELGLU + CL_aZEL*fup_aZEL) 
    
    dAL_aZELGLU <- VmaxL_aZELGLU*CL_aZEL*fup_aZEL / (KmL_aZELGLU + CL_aZEL*fup_aZEL) 
    
    # Bile
    dABile_ZENGLU <- (Vmax_BSEP/(Km_BSEP + (CL_ZENGLU*fup_ZENGLU)))*CL_ZENGLU*fup_ZENGLU +
      - fBile*ABile_ZENGLU
    
    # Blood
    dAB_ZEN <- - (QSI + QL + QR + QS + QF)*CB_ZEN +  
      + (QSI+QL)*CVL_ZEN + QF*CVF_ZEN + QR*CVR_ZEN + QS*CVS_ZEN + 
      - GFR*CB_ZEN*fup_ZEN #+
      # + CLRE_ZEN*CREL_ZEN 

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
    ACalculated <- AST_ZEN +
    ADuoL_ZEN + AJejL_ZEN + AIleL_ZEN + 
    ADuoL_ZENGLU + AJejL_ZENGLU + AIleL_ZENGLU +
    ALIL_ZEN + ALIL_ZENGLU + ALIL_aZEL + ALIL_bZEL +
    AREL_ZEN +
    AFE_ZEN + AFE_aZEL + AFE_ZENGLU + 
    ADuo_ZEN + AJej_ZEN + AIle_ZEN + 
    AJej_ZENGLU + AJej_aZEL + AJej_bZEL + AJej_aZELGLU + 
    AL_ZEN + AL_ZENGLU + AL_bZEL + AL_aZEL + AL_aZELGLU +
    ABile_ZENGLU +
    AB_ZEN + AB_aZEL + AB_ZENGLU +
    AUR_ZEN + AUR_aZEL + AUR_ZENGLU + 
    AF_ZEN + AF_aZEL + AF_ZENGLU + 
    AR_ZEN + AR_aZEL + AR_ZENGLU + 
    AS_ZEN + AS_aZEL + AS_ZENGLU
    
    
    MB <- ATotal - ACalculated + 1
    
    # End
    # 
    list(c(
      dAST_ZEN, 
      dADuoL_ZEN, dAJejL_ZEN, dAIleL_ZEN, 
      dADuoL_ZENGLU, dAJejL_ZENGLU, dAIleL_ZENGLU, 
      dALIL_ZEN, dALIL_ZENGLU, dALIL_aZEL, dALIL_bZEL, 
      dAREL_ZEN, 
      dAFE_ZEN, dAFE_aZEL, dAFE_ZENGLU,
      dADuo_ZEN, dAJej_ZEN, dAIle_ZEN,
      dAJej_ZENGLU, dAJej_aZEL, dAJej_bZEL, dAJej_aZELGLU,
      dAL_ZEN, dAL_ZENGLU, dAL_bZEL, dAL_aZEL, dAL_aZELGLU, 
      dABile_ZENGLU,
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
      CDuoL_ZEN = CDuoL_ZEN, CJejL_ZEN = CJejL_ZEN, CIleL_ZEN = CIleL_ZEN, 
      CDuoL_ZENGLU = CDuoL_ZENGLU, CJejL_ZENGLU = CJejL_ZENGLU, CIleL_ZENGLU = CIleL_ZENGLU, 
      CLIL_ZEN = CLIL_ZEN, CLIL_ZENGLU = CLIL_ZENGLU, CLIL_aZEL = CLIL_aZEL, CLIL_bZEL = CLIL_bZEL,
      CREL_ZEN = CREL_ZEN,
      CDuo_ZEN = CDuo_ZEN, CJej_ZEN = CJej_ZEN, CIle_ZEN = CIle_ZEN, 
      CJej_ZENGLU = CJej_ZENGLU,  CJej_aZEL = CJej_aZEL, CJej_aZELGLU = CJej_aZELGLU,
      CL_ZEN = CL_ZEN, CL_ZENGLU = CL_ZENGLU, CL_bZEL = CL_bZEL, CL_aZEL = CL_aZEL, CL_aZELGLU = CL_aZELGLU,
      CF_ZEN = CF_ZEN, CF_aZEL = CF_aZEL, CF_ZENGLU = CF_ZENGLU,
      CR_ZEN = CR_ZEN, CR_aZEL = CR_aZEL, CR_ZENGLU = CR_ZENGLU, 
      CS_ZEN = CS_ZEN, CS_aZEL = CS_aZEL, CS_ZENGLU = CS_ZENGLU,
      CB_ZEN = CB_ZEN, CB_aZEL = CB_aZEL, CB_ZENGLU = CB_ZENGLU
    ))

  })
}