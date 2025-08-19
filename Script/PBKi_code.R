## PBK model rat ####
PBK.model <- function(t, state, parms){
  with(as.list(c(state, parms)), {
    # MWZEN <- 318.37 
    # Si <- 0.002*10/MWZEN # mol/L, 0.002g/100ml -> *10 g/L, /MW mol/L 10.1016/s0065-2164(08)70160-6
    # 
    ## Calculated parameters ####
    
    ## Physiological constants ####
    # Fractional Tissue volumes (reference: Punt et al. (2020))
    BW <- BW  # Kg # body weight rat
    
    VST <- FV_ST*BW  # fraction volume of stomach
    VSIL <- FV_SIL*BW # fractional volume of small intestinal lumen
    VDuoL <- VSIL*0.183 # calculated from Musther
    VJejL <- VSIL*0.334 # calculated from Musther
    VIleL <- VSIL*0.483 # calculated from Musther
    VCOL <- FV_COL*BW #(FV_COL*BW)-VREL # corrected fractional volume of the rest of the colon
    
    VSI <- FV_SI*BW # fractional volume of small intestinal tissue
    VDuo <- VSI*0.060
    VJej <- VSI*0.769
    VIle <- VSI*0.171
    
    VL <- FV_L*BW  # fractional of liver tissue
    VR <- FV_R*BW  # 0.082-VSIc-VLc # fraction of rapidly perfused tissue
    VS <- FV_S*BW  # 0.737-VFc # fraction of slowly perfused tissue
    VF <- FV_F*BW  # fraction of fat tissue
    VB <- FV_B*BW  # fraction of blood
    
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
    
    # Glomerular filtration rate: Punt et al. (2020)
    GFR <- GFRc*BW # (/h)
    
    ## kinetic parameters ####
    ### absorption/transfer rates ####
    
    # ka_duo <- Peff.cm.per.hr/Rduo
    # ka_jej <- Peff.cm.per.hr/Rjej
    # ka_ile <- Peff.cm.per.hr/Rile
    
    # fa_duo <- 1 - (1 + ka_duo/(1/Tsiduo))**-7
    # fa_jej <- 1 - (1 + ka_jej/(1/Tsijej))**-7
    # fa_ile <- 1 - (1 + ka_ile/(1/Tsiile))**-7
    
    CLduo_ZEN <- PappCaco2_ZEN * SAduo * 1e-3 # cm/h to cm3/h to L/h
    CLjej_ZEN <- PappCaco2_ZEN * SAjej * 1e-3 # cm/h to cm3/h to L/h
    CLile_ZEN <- PappCaco2_ZEN * SAile * 1e-3 # cm/h to cm3/h to L/h
    
    ### metabolism in small intestine tissue ####
    VmaxSI_aZEL <- VmaxSI_aZELc * S9SI * VJej * 1e3       # umol/h, /fuincS9I.HSD_ZEN
    VmaxSI_bZEL <- VmaxSI_bZELc * S9SI * VJej * 1e3       # umol/h, /fuincS9I.HSD_ZEN
    VmaxSI_ZENGLU <- 0.1*VmaxSI_ZENGLUc * S9SI * VJej * 1e3   # umol/h, /fuincS9I.UGT_ZEN 
    VmaxSI_aZELGLU <- 0.1*VmaxSI_aZELGLUc * S9SI * VJej  * 1e3   # umol/h, /fuincS9I_aZEL
    
    ### metabolism in large intestine lumen (microbiota compartment) ####
    VmaxLI_aZEL <- VmaxLI_aZELc*VMB  # umol/h 
    VmaxLI_bZEL <- VmaxLI_bZELc*VMB
    
    ### metabolism in liver ####
    VmaxL_aZEL <- VmaxL_aZELc * S9L * VL * 1e3 #/fuincS9L.HSD_ZEN
    VmaxL_bZEL <- VmaxL_bZELc * S9L * VL * 1e3 #/fuincS9L.HSD_ZEN
    VmaxL_ZENGLU <- 0.1*VmaxL_ZENGLUc * S9L * VL * 1e3 #/fuincS9L.UGT_ZEN
    VmaxL_aZELGLU <- 0.1*VmaxL_aZELGLUc * S9L * VL * 1e3 #/fuincS9L_aZEL 
    
    ### biliary clearance
    CL_Bile <- ifelse(t > 3.5, CL_Bile.fed, CL_Bile.fasted) # 0.00075 L/h; derived from 50 and 70 μl/kg/min kg/min
    Vmax_BSEP <- Vmax_BSEPc*60*24*SF_BSEP*VL      # umol/d
    Km_BSEP <- Km_BSEPc                           # umol/L (uM)
    
    CL.LI.MRP2 = CL.LI.MRP2.c*VL # L/h
    CL.LI.MRP3 = CL.LI.MRP3.c*VL # L/h
    CL.IN.MRP2 = CL.IN.MRP2.c*VSI # L/h
    CL.IN.MRP3 = CL.IN.MRP3.c*VSI # L/h
    
    ## Concentrations #####
    
    # Zearalenone 
    CDuo_ZEN <- ADuo_ZEN/VDuo
    CJej_ZEN <- AJej_ZEN/VJej
    CIle_ZEN <- AIle_ZEN/VIle
    CDuoL_ZEN <- ADuoL_ZEN/VDuoL
    CJejL_ZEN <- AJejL_ZEN/VJejL
    CIleL_ZEN <- AIleL_ZEN/VIleL
    CCOL_ZEN <- ACOL_ZEN/VCOL 
    CL_ZEN <- AL_ZEN/VL
    CF_ZEN <- AF_ZEN/VF
    CR_ZEN <- AR_ZEN/VR
    CS_ZEN <- AS_ZEN/VS
    CB_ZEN <- AB_ZEN/VB
    
    CVDuo_ZEN <- CDuo_ZEN/PSI_ZEN * BP_ZEN
    CVJej_ZEN <- CJej_ZEN/PSI_ZEN * BP_ZEN
    CVIle_ZEN <- CIle_ZEN/PSI_ZEN * BP_ZEN
    CVL_ZEN <- CL_ZEN/PL_ZEN * BP_ZEN
    CVF_ZEN <- CF_ZEN/PF_ZEN * BP_ZEN
    CVR_ZEN <- CR_ZEN/PR_ZEN * BP_ZEN
    CVS_ZEN <- CS_ZEN/PS_ZEN * BP_ZEN
    
    # α-Zearalenol
    CCOL_aZEL <- ACOL_aZEL/VCOL
    CJej_aZEL <- AJej_aZEL/VJej
    CL_aZEL <- AL_aZEL/VL
    CF_aZEL <- AF_aZEL/VF
    CR_aZEL <- AR_aZEL/VR
    CS_aZEL <- AS_aZEL/VS
    CB_aZEL <- AB_aZEL/VB
    
    CVJej_aZEL <- CJej_aZEL/PSI_aZEL * BP_aZEL
    CVL_aZEL <- CL_aZEL/PL_aZEL * BP_aZEL
    CVF_aZEL <- CF_aZEL/PF_aZEL * BP_aZEL
    CVR_aZEL <- CR_aZEL/PR_aZEL * BP_aZEL
    CVS_aZEL <- CS_aZEL/PS_aZEL * BP_aZEL
    
    # β-Zearalenol
    CCOL_bZEL <- ACOL_bZEL/VCOL
    CL_bZEL <- AL_bZEL/VL
    
    # Zearalenone-glucuronide
    CJej_ZENGLU <- AJej_ZENGLU/VJej
    CDuoL_ZENGLU <- ADuoL_ZENGLU/VDuoL
    CJejL_ZENGLU <- AJejL_ZENGLU/VJejL
    CIleL_ZENGLU <- AIleL_ZENGLU/VIleL 
    CCOL_ZENGLU <- ACOL_ZENGLU/VCOL 
    CL_ZENGLU <- AL_ZENGLU/VL
    CF_ZENGLU <- AF_ZENGLU/VF
    CR_ZENGLU <- AR_ZENGLU/VR
    CS_ZENGLU <- AS_ZENGLU/VS
    CB_ZENGLU <- AB_ZENGLU/VB
    
    CVJej_ZENGLU <- CJej_ZENGLU/PSI_ZENGLU * BP_ZENGLU
    CVL_ZENGLU <- CL_ZENGLU/PL_ZENGLU * BP_ZENGLU
    CVF_ZENGLU <- CF_ZENGLU/PF_ZENGLU * BP_ZENGLU
    CVR_ZENGLU <- CR_ZENGLU/PR_ZENGLU * BP_ZENGLU
    CVS_ZENGLU <- CS_ZENGLU/PS_ZENGLU * BP_ZENGLU
    
    # a-Zearalenol glucuronide
    CL_aZELGLU <- AL_aZELGLU/VL   # liver a zel gluc
    CJej_aZELGLU <- AJej_aZELGLU/VJej
    
    ## Dynamic ####
    
    # Stomach (in) 
    dAST_ZEN <- - gaTR*AST_ZEN 
    
    # Small intestinal lumen 
    
    dADuoL_ZEN <- + gaTR*AST_ZEN - duoTR*ADuoL_ZEN + #*(1-fa_duo) + 
      - CLduo_ZEN*CDuoL_ZEN #- ka_duo*ADuoL_ZEN*fa_duo
    dAJejL_ZEN <- + duoTR*ADuoL_ZEN + #*(1-fa_duo) +
      - jejTR*AJejL_ZEN + #*(1-fa_jej) +
      - CLjej_ZEN*CJejL_ZEN + #+#ka_jej*AJejL_ZEN*fa_jej
      + CL_Hydro*CJejL_ZENGLU
    dAIleL_ZEN <- + jejTR*AJejL_ZEN + #*(1-fa_jej) +
      - ileTR*AIleL_ZEN + #*(1-fa_ile) + 
      - CLile_ZEN*CIleL_ZEN + #ka_ile*AIleL_ZEN*fa_ile +
      + CL_Hydro*CIleL_ZENGLU
    
    dADuoL_ZENGLU <- CL_Bile*ABile_ZENGLU + 
      - duoTR*ADuoL_ZENGLU 
    dAJejL_ZENGLU <- + duoTR*ADuoL_ZENGLU + 
      + CL.IN.MRP2*CJej_ZENGLU*fut_ZENGLU + #((VmaxSI_ZENGLU*CJej_ZEN*fut_ZEN) / ((KmSI_ZENGLU) + (CJej_ZEN*fut_ZEN))) +
      - jejTR*AJejL_ZENGLU + 
      - CL_Hydro*CJejL_ZENGLU
    
    dAIleL_ZENGLU <- jejTR*AJejL_ZENGLU - ileTR*AIleL_ZENGLU +
      - CL_Hydro*CIleL_ZENGLU
    
    dACOL_ZEN <- + ileTR*AIleL_ZEN + #*(1-fa_ile) +
      - ((VmaxLI_aZEL*CCOL_ZEN) / (KmLI_aZEL + (CCOL_ZEN))) +
      - ((VmaxLI_bZEL*CCOL_ZEN) / (KmLI_bZEL + (CCOL_ZEN))) +
      + CL_Hydro*CCOL_ZENGLU + 
      - lilTR*ACOL_ZEN 
    
    dACOL_ZENGLU <- ileTR*AIleL_ZENGLU +
      - lilTR*ACOL_ZENGLU +
      - CL_Hydro*CCOL_ZENGLU  
    
    dACOL_aZEL <- ((VmaxLI_aZEL*CCOL_ZEN) / (KmLI_aZEL + (CCOL_ZEN))) +
      - lilTR*ACOL_aZEL 
    
    dACOL_bZEL <- ((VmaxLI_bZEL*CCOL_ZEN) / (KmLI_bZEL + (CCOL_ZEN)))  
    
    # Faeces ####
    dAFE_ZEN <- lilTR*ACOL_ZEN #*(1-fa_lil) 
    dAFE_aZEL <- lilTR*ACOL_aZEL
    dAFE_ZENGLU <- lilTR*ACOL_ZENGLU 
    
    
    # Small intestinal tissue ####
    dADuo_ZEN <- CLduo_ZEN*CDuoL_ZEN + QDuo*(CB_ZEN-CVDuo_ZEN) #ka_duo*ADuoL_ZEN*fa_duo +
    dAJej_ZEN <- CLjej_ZEN*CJejL_ZEN + QJej*(CB_ZEN-CVJej_ZEN) + #ka_jej*AJejL_ZEN*fa_jej + 
      - ((VmaxSI_ZENGLU*CJej_ZEN*fut_ZEN) / ((KmSI_ZENGLU) + (CJej_ZEN*fut_ZEN))) +
      - ((VmaxSI_aZEL*CJej_ZEN*fut_ZEN) / (KmSI_aZEL + (CJej_ZEN*fut_ZEN))) +
      - ((VmaxSI_bZEL*CJej_ZEN*fut_ZEN) / (KmSI_bZEL + (CJej_ZEN*fut_ZEN)))
    
    dAIle_ZEN <- CLile_ZEN*CIleL_ZEN + QIle*(CB_ZEN-CVIle_ZEN) #ka_ile*AIleL_ZEN*fa_ile 
    
    dAJej_ZENGLU <- ((VmaxSI_ZENGLU*CJej_ZEN*fut_ZEN) / ((KmSI_ZENGLU) + (CJej_ZEN*fut_ZEN))) +
      - CL.IN.MRP2*CJej_ZENGLU*fut_ZENGLU +
      + QSI*(CB_ZENGLU-CVJej_ZENGLU) #
    
    dAJej_aZEL <- ((VmaxSI_aZEL*CJej_ZEN*fut_ZEN) / (KmSI_aZEL + (CJej_ZEN*fut_ZEN))) +
      + QSI*(CB_aZEL-CVJej_aZEL) +
      - ((VmaxSI_aZELGLU*CJej_aZEL*fut_aZEL) / (KmSI_aZELGLU + (CJej_aZEL*fut_aZEL)))
    
    dAJej_bZEL <- ((VmaxSI_bZEL*CJej_ZEN*fut_ZEN) / (KmSI_bZEL + (CJej_ZEN*fut_ZEN)))
    
    dAJej_aZELGLU <- ((VmaxSI_aZELGLU*CJej_aZEL*fut_aZEL) / (KmSI_aZELGLU + (CJej_aZEL*fut_aZEL)))
    
    
    # Liver 
    dAL_ZEN <- QDuo*CVDuo_ZEN + QJej*CVJej_ZEN + QIle*CVIle_ZEN + QL*CB_ZEN + #ka_lil*ACO_ZEN*fa_lil + 
      - (QSI+QL)*CVL_ZEN + 
      - ((VmaxL_aZEL*CL_ZEN*fut_ZEN) / (KmL_aZEL + (CL_ZEN*fut_ZEN))) +
      - ((VmaxL_bZEL*CL_ZEN*fut_ZEN) / (KmL_bZEL + (CL_ZEN*fut_ZEN))) +
      - ((VmaxL_ZENGLU*CL_ZEN*fut_ZEN) / ((KmL_ZENGLU) + (CL_ZEN*fut_ZEN)))
    
    dAL_ZENGLU <- +  QSI*CVJej_ZENGLU + QL*CB_ZENGLU - (QSI+QL)*CVL_ZENGLU + 
      + ((VmaxL_ZENGLU*CL_ZEN*fut_ZEN) / ((KmL_ZENGLU) + (CL_ZEN*fut_ZEN))) +
      - CL.LI.MRP2*CL_ZENGLU*fut_ZENGLU #(Vmax_BSEP/(Km_BSEP + (CL_ZENGLU*fut_ZENGLU)))*CL_ZENGLU*fut_ZENGLU
    
    dAL_bZEL <- ((VmaxL_bZEL*CL_ZEN*fut_ZEN) / (KmL_bZEL + (CL_ZEN*fut_ZEN))) 
    
    dAL_aZEL <- + QSI*CVJej_aZEL + QL*CB_aZEL - (QSI+QL)*CVL_aZEL + 
      + ((VmaxL_aZEL*CL_ZEN*fut_ZEN) / (KmL_aZEL + (CL_ZEN*fut_ZEN))) +
      - ((VmaxL_aZELGLU*CL_aZEL*fut_aZEL) / (KmL_aZELGLU + (CL_aZEL*fut_aZEL))) 
    
    dAL_aZELGLU <- ((VmaxL_aZELGLU*CL_aZEL*fut_aZEL) / (KmL_aZELGLU + (CL_aZEL*fut_aZEL))) 
    
    # Bile
    dABile_ZENGLU <- CL.LI.MRP2*CL_ZENGLU*fut_ZENGLU + #(Vmax_BSEP/(Km_BSEP + (CL_ZENGLU*fut_ZENGLU)))*CL_ZENGLU*fut_ZENGLU +
      - CL_Bile*ABile_ZENGLU
    
    # Blood
    dAB_ZEN <- - (QSI + QL + QR + QS + QF)*CB_ZEN +  #
      + (QSI+QL)*CVL_ZEN + QF*CVF_ZEN + QR*CVR_ZEN + QS*CVS_ZEN + 
      - GFR*CB_ZEN*fup_ZEN 
    
    dAB_aZEL <- - (QSI + QL + QR + QS + QF)*CB_aZEL +  #
     + (QL+QSI)*CVL_aZEL + QF*CVF_aZEL + QR*CVR_aZEL + QS*CVS_aZEL +
     - GFR*CB_aZEL*fup_aZEL
    
    dAB_ZENGLU <- - (QSI + QL + QR + QS + QF)*CB_ZENGLU +   # 
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
    ACOL_ZEN + ACOL_ZENGLU + ACOL_aZEL + ACOL_bZEL +
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
      dACOL_ZEN, dACOL_ZENGLU, dACOL_aZEL, dACOL_bZEL, 
      dAFE_ZEN, dAFE_aZEL, dAFE_ZENGLU,
      dADuo_ZEN, dAJej_ZEN, dAIle_ZEN,
      dAJej_ZENGLU, dAJej_aZEL, dAJej_bZEL, dAJej_aZELGLU,
      dAL_ZEN, dAL_ZENGLU, dAL_bZEL, dAL_aZEL, dAL_aZELGLU, 
      dABile_ZENGLU,
      dAB_ZEN, dAB_aZEL, dAB_ZENGLU, 
      dAUR_ZEN, dAUR_aZEL, dAUR_ZENGLU,
      dAF_ZEN, dAF_aZEL, dAF_ZENGLU,
      dAR_ZEN, dAR_aZEL, dAR_ZENGLU,  
      dAS_ZEN, dAS_aZEL, dAS_ZENGLU
    ),
    c(ATotal = ATotal,
      ACalculated = ACalculated,
      MB = MB,
      CDuoL_ZEN = CDuoL_ZEN, CJejL_ZEN = CJejL_ZEN, CIleL_ZEN = CIleL_ZEN, 
      CDuoL_ZENGLU = CDuoL_ZENGLU, CJejL_ZENGLU = CJejL_ZENGLU, CIleL_ZENGLU = CIleL_ZENGLU, 
      CCOL_ZEN = CCOL_ZEN, CCOL_ZENGLU = CCOL_ZENGLU, CCOL_aZEL = CCOL_aZEL, CCOL_bZEL = CCOL_bZEL,
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