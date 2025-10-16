# --------------------------------------------------------------------------- #
# SCRIPT FOR PBK MODEL FOR ORAL EXPOSURE
# By: Chrysanthi Pachoulide
# Date: 14-04-2025
# --------------------------------------------------------------------------- #

# ORAL ####
ORAL_PBK_RUN <- function(y, parms, times){ # Input for ode
  
  PBK.model <- function(t, state, parms){
    with(as.list(c(state, parms)), {
      ## Calculated parameters ####
      
      ## Physiological constants ####
      BW <- BW  
      
      VST <- FV_ST*BW  
      
      # Intestinal lumen
      VSIL <- FV_SIL*BW 
      VDuoL <- VSIL*FV_duol 
      VJejL <- VSIL*FV_jejl 
      VIleL <- VSIL*FV_ilel 
      VCOL <- FV_COL*BW 
      
      # Intestinal tissue
      VSI <- FV_SI*BW 
      VDuo <- VSI*FV_Duo
      VJej <- VSI*FV_Jej
      VJej_IC <- FV_I_IC*VJej
      VJej_EC <- FV_I_EC*VJej
      VIle <- VSI*FV_Ile
      
      # Liver 
      VL <- FV_L*BW  # fractional of liver tissue
      VL_IC <- FV_L_IC*VL
      VL_EC <- FV_L_EC*VL
      
      VBi <- 100/1000 # L from the ICRP, vbile is 50-100ml
      
      VTotc <- 0.96
      VTot <- VTotc * BW              # L, Total body volume (used for mass balance)
      
      # Rest
      VR <- FV_R*BW  
      # VS <- FV_S*BW  # 0.737-VFc # fraction of slowly perfused tissue
      VF <- FV_F*BW  # fraction of fat tissue
      VB <- FV_B*BW  # fraction of blood
      VS <- VTot - (VSI + VL + VR + VF + VB) # fraction of slowly perfused tissue (Bone/Skin/Adipose/Muscle/(brain and testes) (slowly equilibrating tissues with large time constants))
      
      
      # Blood flow rates (reference: Punt et al. (2020))
      QC <- QC  ##{L/h} # cardiac output
      QSI <- QSIc*QC  # fraction of blood flow to small intestine
      QDuo <- QSI*QDuoc
      QJej <- QSI*QJejc
      QIle <- QSI*QIlec
      # QCO <- 0 #QCOc*QC
      QL <- QLc*QC  # 0.242 - QSIc # fraction of blood flow to liver 101 
      QR <- QRc*QC  # 0.76 - QSIc - QLc # fraction of blood flow to rapidly perfused tissue
      # QS <- QSc*QC  # 0.24 - QFc # fraction of blood flow to slowly perfused tissue
      QF <- QFc*QC  # fraction of blood flow to fat
      QS <- QC - (QSI + QL + QR + QF) 
      
      # Glomerular filtration rate: Punt et al. (2020)
      GFR <- GFR # (L/h)
      
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
      # CLco_ZEN <- 0 #0.001*PappCaco2_ZEN * SAco * 1e-3 # cm/h to cm3/h to L/h
      
      ### metabolism in small intestine tissue ####
      # MISI <- 2.48
      VmaxSI_aZEL <- VmaxSI_aZELc * S9SI * VJej * 1e3       # umol/h, /fuincS9I.HSD_ZEN
      VmaxSI_bZEL <- VmaxSI_bZELc * S9SI * VJej * 1e3       # umol/h, /fuincS9I.HSD_ZEN
      VmaxJej_ZENGLU <- 0.3*VmaxSI_ZENGLUc * S9SI * VJej * 1e3   # umol/h, /fuincS9I.UGT_ZEN 
      # VmaxIle_ZENGLU <- 0 #VmaxSI_ZENGLUc * S9SI * VIle * 1e3   # umol/h, /fuincS9I.UGT_ZEN 
      # VmaxCO_ZENGLU <- 0 #VmaxSI_ZENGLUc * S9SI * VCO * 1e3   # umol/h, /fuincS9I.UGT_ZEN 
      
      VmaxSI_aZELGLU <- VmaxSI_aZELGLUc * S9SI * VJej  * 1e3   # umol/h, /fuincS9I_aZEL
      
      
      ### metabolism in large intestine lumen (microbiota compartment) ####
      VmaxLI_aZEL <- VmaxLI_aZELc*VMB  # umol/h 
      VmaxLI_bZEL <- VmaxLI_bZELc*VMB
      
      ### metabolism in liver ####
      # MILI <- 32
      VmaxL_aZEL <- VmaxL_aZELc * S9L * VL * 1e3 #/fuincS9L.HSD_ZEN
      VmaxL_bZEL <- VmaxL_bZELc * S9L * VL * 1e3 #/fuincS9L.HSD_ZEN
      VmaxL_ZENGLU <- 0.3*VmaxL_ZENGLUc * S9L * VL * 1e3 #/fuincS9L.UGT_ZEN
      VmaxL_aZELGLU <- VmaxL_aZELGLUc * S9L * VL * 1e3 #/fuincS9L_aZEL 
      
      ### biliary clearance
      # Bile.fed <- 1
      # CL_Bile <- ifelse(t > 3, Bile.fed, Bile.fasted) # 0.00075 L/h; derived from 50 and 70 μl/kg/min kg/min
      CL_Bile = ifelse(((t > (8/24) & t < ((8/24) +  Tgb)) |
                          (t > (14/24) & t < ((14/24) + Tgb)) |
                          (t > (20/24) & t < ((20/24) + Tgb))), Bile.fed, Bile.fasted)
      
      CL.LI.MRP2 = CL.LI.MRP2.c*VL # L/h
      CL.LI.MRP3 = CL.LI.MRP3.c*VL # L/h
      CL.IN.MRP2 = CL.IN.MRP2.c*VSI # L/h
      CL.IN.MRP3 = CL.IN.MRP3.c*VSI # L/h
      CL.KI.MRP3 = CL.KI.MRP3.c*VB # L/h
      
      # lilTR = lilTR
      
      ## Concentrations #####
      
      # Zearalenone 
      CDuo_ZEN <- ADuo_ZEN/VDuo
      CJej_ZEN <- AJej_ZEN/VJej
      CIle_ZEN <- AIle_ZEN/VIle
      CDuoL_ZEN <- ADuoL_ZEN/VDuoL
      CJejL_ZEN <- AJejL_ZEN/VJejL
      CIleL_ZEN <- AIleL_ZEN/VIleL
      # CCO_ZEN <- ACOL_ZEN/VCO
      CCOL_ZEN <- ACOL_ZEN/VCOL 
      CL_ZEN <- AL_ZEN/VL
      CF_ZEN <- AF_ZEN/VF
      CR_ZEN <- AR_ZEN/VR
      CS_ZEN <- AS_ZEN/VS
      CB_ZEN <- AB_ZEN/VB
      
      CVDuo_ZEN <- CDuo_ZEN/PSI_ZEN * BP_ZEN
      CVJej_ZEN <- CJej_ZEN/PSI_ZEN * BP_ZEN
      CVIle_ZEN <- CIle_ZEN/PSI_ZEN * BP_ZEN
      # CVCO_ZEN <- CCO_ZEN/PSI_ZEN * BP_ZEN
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
      CDuoL_ZENGLU <- ADuoL_ZENGLU/VDuoL
      
      CJej_IC_ZENGLU <- AJej_IC_ZENGLU/VJej_IC
      CJej_EC_ZENGLU <- AJej_EC_ZENGLU/VJej_EC
      CJejL_ZENGLU <- AJejL_ZENGLU/VJejL
      
      # CIle_IC_ZENGLU <- AIle_IC_ZENGLU/VIle_IC
      # CIle_EC_ZENGLU <- AIle_EC_ZENGLU/VIle_EC
      CIleL_ZENGLU <- AIleL_ZENGLU/VIleL 
      
      # CCO_IC_ZENGLU <- ACO_IC_ZENGLU/VCO_IC 
      # CCO_EC_ZENGLU <- ACO_EC_ZENGLU/VCO_EC 
      # 
      CCOL_ZENGLU <- ACOL_ZENGLU/VCOL 
      
      # CL_ZENGLU <- AL_ZENGLU/VL
      CL_IC_ZENGLU <- AL_IC_ZENGLU/VL_IC
      CL_EC_ZENGLU <- AL_EC_ZENGLU/VL_EC
      
      CBile_ZENGLU <- ABile_ZENGLU/VBi
      
      CF_ZENGLU <- AF_ZENGLU/VF
      CR_ZENGLU <- AR_ZENGLU/VR
      CS_ZENGLU <- AS_ZENGLU/VS
      CB_ZENGLU <- AB_ZENGLU/VB
      
      # CVJej_ZENGLU <- CJej_ZENGLU/PSI_ZENGLU * BP_ZENGLU
      CVJej_EC_ZENGLU <- CJej_EC_ZENGLU/PSI_ZENGLU * BP_ZENGLU
      # CVIle_EC_ZENGLU <- CIle_EC_ZENGLU/PSI_ZENGLU * BP_ZENGLU
      # CVCO_EC_ZENGLU <- CCO_EC_ZENGLU/PSI_ZENGLU * BP_ZENGLU
      
      # CVL_ZENGLU <- CL_ZENGLU/PL_ZENGLU * BP_ZENGLU
      CVL_EC_ZENGLU <- CL_EC_ZENGLU/PL_ZENGLU * BP_ZENGLU
      
      CVF_ZENGLU <- CF_ZENGLU/PF_ZENGLU * BP_ZENGLU
      CVR_ZENGLU <- CR_ZENGLU/PR_ZENGLU * BP_ZENGLU
      CVS_ZENGLU <- CS_ZENGLU/PS_ZENGLU * BP_ZENGLU
      
      # a-Zearalenol glucuronide
      CL_aZELGLU <- AL_aZELGLU/VL   # liver a zel gluc
      CJej_aZELGLU <- AJej_aZELGLU/VJej
      
      # Oral exposure ##
      if(t<expSTOP){DoseOn=1} else{DoseOn=0}
      
      
      Oral = expOral*BW*DoseOn         # ug, PFOA oral dose
      OralD = Oral/Tinput*(t %% tinterval<Tinput)
      
      ## Dynamic ####
      # dOD = OralD - OD             # umol/d, Oral dose input 
      
      # Stomach (in) 
      dAST_ZEN <- + OralD - gaTR*AST_ZEN 
      
      # Small intestinal lumen 
      
      dADuoL_ZEN <- + gaTR*AST_ZEN - duoTR*ADuoL_ZEN + #*(1-fa_duo) + 
        - CLduo_ZEN*CDuoL_ZEN #- ka_duo*ADuoL_ZEN*fa_duo
      dAJejL_ZEN <- + duoTR*ADuoL_ZEN + #*(1-fa_duo) +
        - jejTR*AJejL_ZEN + #*(1-fa_jej) +
        - CLjej_ZEN*CJejL_ZEN #+#ka_jej*AJejL_ZEN*fa_jej
      dAIleL_ZEN <- + jejTR*AJejL_ZEN + #*(1-fa_jej) +
        - ileTR*AIleL_ZEN + #*(1-fa_ile) + 
        - CLile_ZEN*CIleL_ZEN + #ka_ile*AIleL_ZEN*fa_ile +
        + CL_Hydro*CIleL_ZENGLU
      
      dADuoL_ZENGLU <- CL_Bile*CBile_ZENGLU + 
        - duoTR*ADuoL_ZENGLU 
      dAJejL_ZENGLU <- + duoTR*ADuoL_ZENGLU + 
        + CL.IN.MRP2*CJej_IC_ZENGLU*fut_ZENGLU + #CL.IN.MRP2*CJej_ZENGLU*fut_ZENGLU + #((VmaxSI_ZENGLU*CJej_ZEN*fut_ZEN) / ((KmSI_ZENGLU) + (CJej_ZEN*fut_ZEN))) +
        - jejTR*AJejL_ZENGLU 
      
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
        - ((VmaxJej_ZENGLU*CJej_ZEN*fut_ZEN) / ((KmSI_ZENGLU) + (CJej_ZEN*fut_ZEN))) +
        - ((VmaxSI_aZEL*CJej_ZEN*fut_ZEN) / (KmSI_aZEL + (CJej_ZEN*fut_ZEN))) +
        - ((VmaxSI_bZEL*CJej_ZEN*fut_ZEN) / (KmSI_bZEL + (CJej_ZEN*fut_ZEN)))
      
      dAIle_ZEN <- CLile_ZEN*CIleL_ZEN + QIle*(CB_ZEN-CVIle_ZEN) #+ #ka_ile*AIleL_ZEN*fa_ile 
      
      
      dAJej_IC_ZENGLU <- ((VmaxJej_ZENGLU*CJej_ZEN*fut_ZEN) / ((KmSI_ZENGLU) + (CJej_ZEN*fut_ZEN))) +
        - CL.IN.MRP2*CJej_IC_ZENGLU*fut_ZENGLU +
        - CL.IN.MRP3*CJej_IC_ZENGLU*fut_ZENGLU
      
      
      dAJej_EC_ZENGLU <- + QJej*(CB_ZENGLU-CVJej_EC_ZENGLU) +
        CL.IN.MRP3*CJej_IC_ZENGLU*fut_ZENGLU
      
      dAJej_aZEL <- ((VmaxSI_aZEL*CJej_ZEN*fut_ZEN) / (KmSI_aZEL + (CJej_ZEN*fut_ZEN))) +
        + QSI*(CB_aZEL-CVJej_aZEL) +
        - ((VmaxSI_aZELGLU*CJej_aZEL*fut_aZEL) / (KmSI_aZELGLU + (CJej_aZEL*fut_aZEL)))
      
      dAJej_bZEL <- ((VmaxSI_bZEL*CJej_ZEN*fut_ZEN) / (KmSI_bZEL + (CJej_ZEN*fut_ZEN)))
      
      dAJej_aZELGLU <- ((VmaxSI_aZELGLU*CJej_aZEL*fut_aZEL) / (KmSI_aZELGLU + (CJej_aZEL*fut_aZEL)))
      
      # Liver 
      dAL_ZEN <- QDuo*CVDuo_ZEN + QJej*CVJej_ZEN + QIle*CVIle_ZEN + #QCO*CVCO_ZEN +
        + QL*CB_ZEN + #ka_lil*ACO_ZEN*fa_lil + 
        - (QSI+QL)*CVL_ZEN + 
        - ((VmaxL_aZEL*CL_ZEN*fut_ZEN) / (KmL_aZEL + (CL_ZEN*fut_ZEN))) +
        - ((VmaxL_bZEL*CL_ZEN*fut_ZEN) / (KmL_bZEL + (CL_ZEN*fut_ZEN))) +
        - ((VmaxL_ZENGLU*CL_ZEN*fut_ZEN) / ((KmL_ZENGLU) + (CL_ZEN*fut_ZEN)))
      
      dAL_IC_ZENGLU <- + ((VmaxL_ZENGLU*CL_ZEN*fut_ZEN) / ((KmL_ZENGLU) + (CL_ZEN*fut_ZEN))) +
        - CL.LI.MRP2*CL_IC_ZENGLU*fut_ZENGLU + #(Vmax_BSEP/(Km_BSEP + (CL_ZENGLU*fut_ZENGLU)))*CL_ZENGLU*fut_ZENGLU
        - CL.LI.MRP3*CL_IC_ZENGLU*fut_ZENGLU
      
      dAL_EC_ZENGLU <- + QJej*CVJej_EC_ZENGLU + #QIle*CVIle_EC_ZENGLU + #QCO*CVCO_EC_ZENGLU + 
        + (QL+QDuo+QIle)*CB_ZENGLU - (QJej + QL)*CVL_EC_ZENGLU +
        + CL.LI.MRP3*CL_IC_ZENGLU*fut_ZENGLU
      
      dAL_bZEL <- ((VmaxL_bZEL*CL_ZEN*fut_ZEN) / (KmL_bZEL + (CL_ZEN*fut_ZEN))) 
      
      dAL_aZEL <- + QSI*CVJej_aZEL + QL*CB_aZEL - (QSI+QL)*CVL_aZEL + 
        + ((VmaxL_aZEL*CL_ZEN*fut_ZEN) / (KmL_aZEL + (CL_ZEN*fut_ZEN))) +
        - ((VmaxL_aZELGLU*CL_aZEL*fut_aZEL) / (KmL_aZELGLU + (CL_aZEL*fut_aZEL))) 
      
      dAL_aZELGLU <- ((VmaxL_aZELGLU*CL_aZEL*fut_aZEL) / (KmL_aZELGLU + (CL_aZEL*fut_aZEL))) 
      
      # Bile
      dABile_ZENGLU <- CL.LI.MRP2*CL_IC_ZENGLU*fut_ZENGLU  + #CL.LI.MRP2*CL_ZENGLU*fut_ZENGLU + #(Vmax_BSEP/(Km_BSEP + (CL_ZENGLU*fut_ZENGLU)))*CL_ZENGLU*fut_ZENGLU +
        - CL_Bile*CBile_ZENGLU
      
      # Blood
      dAB_ZEN <- - (QSI + QL + QR + QS + QF)*CB_ZEN +  #
        + (QSI + QL)*CVL_ZEN + QF*CVF_ZEN + QR*CVR_ZEN + QS*CVS_ZEN + 
        - GFR*CB_ZEN*fup_ZEN 
      
      dAB_aZEL <- - (QSI + QL + QR + QS + QF)*CB_aZEL +  #
        + (QL+QSI)*CVL_aZEL + QF*CVF_aZEL + QR*CVR_aZEL + QS*CVS_aZEL +
        - GFR*CB_aZEL*fup_aZEL
      
      dAB_ZENGLU <- - (QJej + QDuo + QIle + QL + QR + QS + QF)*CB_ZENGLU +   # 
        + (QJej + QL)*CVL_EC_ZENGLU + QF*CVF_ZENGLU + QR*CVR_ZENGLU + QS*CVS_ZENGLU +
        - GFR*CB_ZENGLU*fup_ZENGLU +
        - CL.KI.MRP3*CB_ZENGLU*fup_ZENGLU
      
      # Urine
      dAUR_ZEN <- GFR*CB_ZEN*fup_ZEN 
      dAUR_aZEL <- GFR*CB_aZEL*fup_aZEL
      dAUR_ZENGLU <- GFR*CB_ZENGLU*fup_ZENGLU + CL.KI.MRP3*CB_ZENGLU*fup_ZENGLU
      
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
      dAin <- OralD # to be used if repeated exposure
      Atot <- AST_ZEN +
        ADuoL_ZEN + AJejL_ZEN + AIleL_ZEN + 
        ADuoL_ZENGLU + AJejL_ZENGLU + AIleL_ZENGLU +
        ACOL_ZEN + ACOL_ZENGLU + ACOL_aZEL + ACOL_bZEL +
        AFE_ZEN + AFE_aZEL + AFE_ZENGLU + 
        ADuo_ZEN + AJej_ZEN + AIle_ZEN + 
        AJej_IC_ZENGLU + AJej_EC_ZENGLU + #AIle_IC_ZENGLU + AIle_EC_ZENGLU +
        AJej_aZEL + AJej_bZEL + AJej_aZELGLU + #AJej_ZENGLU
        AL_ZEN + AL_IC_ZENGLU + AL_EC_ZENGLU + AL_bZEL + AL_aZEL + AL_aZELGLU + #AL_ZENGLU
        ABile_ZENGLU +
        AB_ZEN + AB_aZEL + AB_ZENGLU +
        AUR_ZEN + AUR_aZEL + AUR_ZENGLU + 
        AF_ZEN + AF_aZEL + AF_ZENGLU + 
        AR_ZEN + AR_aZEL + AR_ZENGLU + 
        AS_ZEN + AS_aZEL + AS_ZENGLU
      
      
      MB <- Ain - Atot + 1
      
      # End
      list(c(
        dAST_ZEN, 
        dADuoL_ZEN, dAJejL_ZEN, dAIleL_ZEN, 
        dADuoL_ZENGLU, dAJejL_ZENGLU, dAIleL_ZENGLU,
        dACOL_ZEN, dACOL_ZENGLU, dACOL_aZEL, dACOL_bZEL, 
        dAFE_ZEN, dAFE_aZEL, dAFE_ZENGLU,
        dADuo_ZEN, dAJej_ZEN, dAIle_ZEN,
        dAJej_IC_ZENGLU, dAJej_EC_ZENGLU, #dAIle_IC_ZENGLU, dAIle_EC_ZENGLU,
        dAJej_aZEL, dAJej_bZEL, dAJej_aZELGLU,#dAJej_ZENGLU
        dAL_ZEN, dAL_IC_ZENGLU, dAL_EC_ZENGLU, dAL_bZEL, dAL_aZEL, dAL_aZELGLU, #dAL_ZENGLU
        dABile_ZENGLU,
        dAB_ZEN, dAB_aZEL, dAB_ZENGLU, 
        dAUR_ZEN, dAUR_aZEL, dAUR_ZENGLU,
        dAF_ZEN, dAF_aZEL, dAF_ZENGLU,
        dAR_ZEN, dAR_aZEL, dAR_ZENGLU,  
        dAS_ZEN, dAS_aZEL, dAS_ZENGLU, 
        dAin
      ),
      c(Atot = Atot,
        MB = MB,
        CDuoL_ZEN = CDuoL_ZEN, CJejL_ZEN = CJejL_ZEN, CIleL_ZEN = CIleL_ZEN, 
        CDuoL_ZENGLU = CDuoL_ZENGLU, CJejL_ZENGLU = CJejL_ZENGLU, CIleL_ZENGLU = CIleL_ZENGLU, 
        CCOL_ZEN = CCOL_ZEN, CCOL_ZENGLU = CCOL_ZENGLU, CCOL_aZEL = CCOL_aZEL, CCOL_bZEL = CCOL_bZEL,
        CDuo_ZEN = CDuo_ZEN, CJej_ZEN = CJej_ZEN, CIle_ZEN = CIle_ZEN, 
        CJej_ZENGLU = CJej_IC_ZENGLU, CJej_aZEL = CJej_aZEL, CJej_aZELGLU = CJej_aZELGLU, # CJej_ZENGLU = CJej_ZENGLU; (((CJej_IC_ZENGLU/VJej_IC) + (CJej_EC_ZENGLU/VJej_EC))/(VJej_IC + VJej_EC)); CIle_ZENGLU = CIle_IC_ZENGLU, 
        CL_ZEN = CL_ZEN, CL_ZENGLU = CL_IC_ZENGLU, CL_bZEL = CL_bZEL, CL_aZEL = CL_aZEL, CL_aZELGLU = CL_aZELGLU, #CL_ZENGLU = CL_ZENGLU; (((CL_IC_ZENGLU/VL_IC) + (CL_EC_ZENGLU/VL_EC))/(VL_IC + VL_EC))
        CF_ZEN = CF_ZEN, CF_aZEL = CF_aZEL, CF_ZENGLU = CF_ZENGLU,
        CR_ZEN = CR_ZEN, CR_aZEL = CR_aZEL, CR_ZENGLU = CR_ZENGLU, 
        CS_ZEN = CS_ZEN, CS_aZEL = CS_aZEL, CS_ZENGLU = CS_ZENGLU,
        CB_ZEN = CB_ZEN, CB_aZEL = CB_aZEL, CB_ZENGLU = CB_ZENGLU
      ))
      
    })
  }
  
  output_PBK <-lsoda(y = y, 
                      times = times, 
                      func = PBK.model,
                      parms = parms)
  
  return(as.data.frame(output_PBK))
}
  
# # DERMAL ####
# DERMAL_PBK_RUN <- function(y, parms, times){ # Input for ode
#   
#   
#   PBK.model <- function(t, state, parameters){
#     with(as.list(c(state, parameters)), {
#       
#       ### Physiological ----
#       
#       VSkB <- VSkB                    # L, Volume of skin barrier
#       VSk <- VSkc * BW                # L, Volume of skin
#       
#       VIL <- VILc * BW                # L, Volume of intestinal lumen
#       VI <- VIc * BW                  # L, Volume of intestine
#       SA_SI <- SA_SIc * BW
#       
#       
#       VL <- VLc * BW                  # L, Volume of liver
#       VL_ic <- VL_icc*VL              # L, Volume liver intracellular                 
#       VL_ec <- VL_ecc*VL              # L, Volume liver extracellular
#       
#       VK <- VKc*BW                    # L, Volume of kidney
#       VPT <- VPTc*VK                  # L, Volume of proximal tubule
#       VPTT <- VPT*VPTTc               # L, Volume of proximal tubule tissue
#       VPTL <- VPT*VPTLc               # L, Volume of proximal tubule lumen
#       VRKL <- (VK - VPT)*VRKLc        # L, Volume of rest of kidney lumen
#       VRKT <- VK - VPTT - VPTL - VRKL # L, Volume of rest of kidney tissue, used in the model
#       
#       VA <- VAc * BW/0.9 + MAc/0.9    # L, Volume of adipose
#       
#       VP <- VPc * (1-Hct) * BW                  # L, Volume of plasma
#       
#       VTotc <- 0.96
#       VTot <- VTotc * BW              # L, Total body volume (used for mass balance)
#       
#       VR <- VTot - (VSk + VI + VL + VK + VA + VP)             # L, Volume of the lumped rest compartment
#       
#       QTotc <- 0.988
#       QSk <- QSkc/QTotc * QC                # L/d, Skin 
#       QI <- QIc/QTotc * QC                  # L/d, Intestinal
#       QL <- QLc/QTotc * QC                  # L/d, Liver
#       QK <- QKc/QTotc * QC                  # L/d, Kidney
#       QA <- QAc/QTotc * QC                  # L/d, Adipose 
#       QR <- QC - (QA + QI + QK + QL + QSk)  # L/d, Rest
#       
#       QUr <- QUrc * BW              # L/d, Urine flow rate to the bladder 22 mL/kg BW/d [ICRP 89 page 161]
#       GFR <- GFR                    # L/d calculated based on age and creatinine levels
#       QT <- QTc*GFR                 # L/d, Proximal tubule fluid flow
#       
#       tco <- tco                    # /d, Bowel residence times in the colon
#       
#       
#       ### Physicochemical ----
#       MW <- 414.07
#       
#       # Fraction unbound in plasma
#       fup <- fup # 0.061/100 #
#       
#       PSk <- PSkc * fup  # Skin
#       PI <- PIc * fup    # Intestinal
#       PL <- PLc * fup    # Liver
#       PK <- PKc * fup    # Kidney
#       PA <- PAc * fup    # Adipose
#       PR <- PRc * fup    # Rest
#       
#       # Fraction unbound, calculated based on Poulin and Haddad, 2018 https://doi.org/10.1016/j.xphs.2018.03.012
#       # Equation was adapted to not account for fraction unionised
#       # OAT and OATP transporters transport the ionised compound, given that the ratio of fraction ionised at plasma to cellular pH is 1, this can be ignored (fraction unionised of PFOA is 0.9999923 at pH 7.4 and 0.9999963 at pH 7)
#       fu_PTL <- R_PTL*fup/(1 + ((R_PTL-1)*fup)) # Proximal tubule lumen
#       fu_Lic <- R_L_ec*fup/(1 + ((R_L_ec-1)*fup)) # Liver intracellular space
#       
#       
#       
#       ### Kinetic ----
#       
#       # Skin uptake
#       CL_SkBtSk <- Papp_SkB*SA_SkB*1e-3*60*60*24           # L/d, Skin barrier to skin (calculations: cm/s -> L/s * 1e-3 -> L/d *60*60*24)
#       
#       # Gastro-intestinal uptake
#       CL_IL <- (Papp_SI*SA_SI*1e-3)*60*60*24    # L/d, Intestinal lumen to intestinal tissue (calculations: cm/s -> L/s /1000 -> L/d *60*60*24) 
#       
#       VmaxOATP2B1 <- Vmax_OATP2B1c*MW*60*24*SF_OATP2B1*VIL # ug/d
#       Km_OATP2B1 <- Km_OATP2B1c*MW                         # ug/L (uM -> ug/L)
#       
#       # Liver uptake
#       Vmax_OATP1B1 <- Vmax_OATP1B1c*MW*60*24*SF_OATP1B1*VL_ec             # ug/d
#       Km_OATP1B1 <- Km_OATP1B1c*MW                                        # ug/L (uM -> ug/L)
#       Vmax_OATP1B3 <- Vmax_OATP1B3c*MW*60*24*SF_OATP1B3*VL_ec             # ug/d
#       Km_OATP1B3 <- Km_OATP1B3c*MW                                        # ug/L (uM -> ug/L)
#       
#       # Biliary excretion
#       Vmax_BSEP <- Vmax_BSEPc*MW*60*24*SF_BSEP*VL_ic         # ug/d
#       Km_BSEP <- Km_BSEPc*MW                                 # ug/L (uM -> ug/L)
#       
#       # Renal clearance
#       Vmax_OAT4 = Vmax_OAT4c*MW*60*24*SF_OAT*VPT           # ug/d (umol -> ug, min -> d)
#       Km_OAT4 = Km_OAT4c*MW                                # ug/L, scaled from uM, Louisse et al. 2024 doi.org/10.1016/j.tox.2024.153961
#       
#       # Clearance via menstruation
#       CL_menses = CL_menses     # L/d, plasma menstrual loss
#       
#       
#       ## Dose -------------------------
#       
#       if(t<expSTOP){DoseOn=1} else{DoseOn=0}
#       
#       ## Dermal exposure ##
#       DDermal = expDermal*BW*DoseOn         # ug, PFOA dermal dose
#       DermalD = DDermal/Tinput*(t %% tinterval<Tinput)
#       
#       ## Concentrations -------------------------
#       
#       CSkB <- ASkB/VSkB            # ug/L, Skin barrier
#       CSk <- ASk/VSk               # ug/L, Skin
#       CVSk <- CSk/PSk              # ug/L, Skin venous 
#       
#       CIL <- AIL/VIL               # ug/L, Intestinal lumen
#       CI <- AI/VI                  # ug/L, Intestine 
#       CVI <- CI/PI                 # ug/L, Intestine venous
#       
#       CL_ec <- AL_ec/VL_ec         # ug/L, Liver extracellular
#       CL_ic <- AL_ic/VL_ic         # ug/L, Liver intracellular
#       CVL_ec <- CL_ec/PL           # ug/L, Liver extracellular venous
#       
#       CPTT <- APTT/VPTT            # ug/L, Kidney, proximal tubule tissue
#       CPTL <- APTL/VPTL            # ug/L, Kidney, proximal tubule lumen
#       CRKT <- ARKT/VRKT            # ug/L, Rest of kidney tissue
#       CRKL <- ARKL/VRKL            # ug/L, Rest of kidney lumen
#       CVRKT <- CRKT/PK             # ug/L, Rest of kidney venous
#       
#       CA <- AA/VA                  # ug/L, Adipose
#       CVA <- CA/PA                 # ug/L, Adipose venous
#       
#       CR <- AR/VR                  # ug/L, Rest
#       CVR <- CR/PR                 # ug/L, Rest venous
#       
#       CP <- AP/VP                  # ug/L, Plasma
#       
#       
#       ## Differential equations -------------------------
#       
#       dDD = DermalD - DD                             # ug/d, Dermal dose input
#       
#       dASkB <- DD - CL_SkBtSk*CSkB                                # ug/d, Skin barrier
#       
#       dASk <- + CL_SkBtSk*CSkB + QSk*(CP-CVSk)                   # ug/d, Skin
#       
#       
#       dAIL <- - tco*AIL - CL_IL*CIL +
#         - (VmaxOATP2B1/(Km_OATP2B1 + CIL))*CIL +
#         + (Vmax_BSEP/(Km_BSEP + (CL_ic*fu_Lic)))*CL_ic*fu_Lic          # ug/d, Intestine lumen
#       
#       dAI <- QI*(CP - CVI) + CL_IL*CIL + 
#         + (VmaxOATP2B1/(Km_OATP2B1 + CIL))*CIL                     # ug/d, Intestinal
#       
#       
#       dAFe <-  tco*AIL                                       # ug/d, Feces
#       
#       
#       dAL_ec <- + QI*CVI + QL*CP - (QI+QL)*CVL_ec + 
#         - (Vmax_OATP1B1/(Km_OATP1B1 + (CL_ec*fup)))*CL_ec*fup +
#         - (Vmax_OATP1B3/(Km_OATP1B3 + (CL_ec*fup)))*CL_ec*fup            # ug/d, Liver extracellular space (vascular + interstitial space)
#       
#       dAL_ic <- (Vmax_OATP1B1/(Km_OATP1B1 + (CL_ec*fup)))*CL_ec*fup +
#         + (Vmax_OATP1B3/(Km_OATP1B3 + (CL_ec*fup)))*CL_ec*fup +
#         - (Vmax_BSEP/(Km_BSEP + (CL_ic*fu_Lic)))*CL_ic*fu_Lic                       # ug/d, Liver intracellular space
#       
#       dAPTT <- QK*(CP - CPTT) + 
#         + (Vmax_OAT4/(Km_OAT4+(CPTL*fu_PTL)))*CPTL*fu_PTL        # ug/d, Proximal tubule tissue 
#       
#       dAPTL <- + fup*GFR*CP - QT*CPTL +
#         - (Vmax_OAT4/(Km_OAT4+(CPTL*fu_PTL)))*CPTL*fu_PTL        # ug/d, Proximal tubule lumen    
#       
#       dARKT <- QK*(CPTT - CVRKT)                               # ug/d, Rest of kidney
#       
#       dARKL <- QT*CPTL - QUr*CRKL                              # ug/d, Rest of kidney lumen
#       
#       dAUr <- QUr*CRKL                                         # ug/d, Urine
#       
#       
#       dAA <- QA*(CP-CVA)                                      # ug/d, Adipose
#       
#       
#       dAR <- QR*(CP-CVR)                                      # ug/d, Rest
#       
#       dAP <- - (QSk + QI + QL + QA + QR + QK)*CP - fup*GFR*CP +     # ug/d, Arterial Plasma
#         + (QL+QI)*CVL_ec + QSk*CVSk + QK*CVRKT + QA*CVA + QR*CVR +   # ug/d, Venous Plasma
#       - CL_menses*CP # ug/d, Venous Plasma
#       
#       dAMp <- + CL_menses*CP
#       
#       
#       # Mass Balance
#       Atot <- DD + 
#         ASk +
#         AIL + AI + AFe + 
#         AL_ec + AL_ic +
#         APTT + APTL + ARKT + ARKL + AUr +
#         AA + 
#         AR +
#         AP +
#         AMp
#       
#       dAin <- DermalD # to be used if repeated exposure
#       MB <- Ain - Atot + 1    # to be used if repeated exposure
#       # MB <- DOral - Atot + 1
#       
#       # End
#       
#       list(c(dDD,
#              dASkB,
#              dASk, 
#              dAIL,
#              dAI, 
#              dAFe,
#              dAL_ec,
#              dAL_ic,
#              dAPTT,
#              dAPTL,
#              dARKT,
#              dARKL,
#              dAUr,
#              dAA,
#              dAR, 
#              dAP, 
#              dAMp,
#              dAin
#       ), 
#       c(CSkB = CSkB, 
#         CSk = CSk, 
#         CIL = CIL,
#         CI = CI, CVI = CVI, 
#         CL_ec = CL_ec,
#         CL_ic = CL_ic,
#         CVL_ec = CVL_ec,
#         CPTT = CPTT,
#         CPTL = CPTL,
#         CRKT = CRKT,
#         CRKL = CRKL,
#         CA = CA, CVA = CVA,
#         CR = CR, CVR = CVR,
#         CP = CP,
#         Atot = Atot, 
#         MB = MB
#       )
#       )
#     })
#   }
#   
#   output_PBK <-lsoda(y = y, 
#                       times = times, 
#                       func = PBK.model,
#                       parms = parms)
#   
#   return(as.data.frame(output_PBK))
# }
# 
# 
# 
# 
# # ORAL & DERMAL ####
# ORAL_DERMAL_PBK_RUN <- function(y, parms, times){ # Input for ode
#   
#   
#   PBK.model <- function(t, state, parameters){
#     with(as.list(c(state, parameters)), {
#       
#       ### Physiological ----
#       
#       VSkB <- VSkB                    # L, Volume of skin barrier
#       VSk <- VSkc * BW                # L, Volume of skin
#       
#       VIL <- VILc * BW                # L, Volume of intestinal lumen
#       VI <- VIc * BW                  # L, Volume of intestine
#       SA_SI <- SA_SIc * BW
#       
#       
#       VL <- VLc * BW                  # L, Volume of liver
#       VL_ic <- VL_icc*VL              # L, Volume liver intracellular                 
#       VL_ec <- VL_ecc*VL              # L, Volume liver extracellular
#       
#       VK <- VKc*BW                    # L, Volume of kidney
#       VPT <- VPTc*VK                  # L, Volume of proximal tubule
#       VPTT <- VPT*VPTTc               # L, Volume of proximal tubule tissue
#       VPTL <- VPT*VPTLc               # L, Volume of proximal tubule lumen
#       VRKL <- (VK - VPT)*VRKLc        # L, Volume of rest of kidney lumen
#       VRKT <- VK - VPTT - VPTL - VRKL # L, Volume of rest of kidney tissue, used in the model
#       
#       VA <- VAc * BW/0.9 + MAc/0.9    # L, Volume of adipose
#       
#       VP <- VPc * (1-Hct) * BW                  # L, Volume of plasma
#       
#       VTotc <- 0.96
#       VTot <- VTotc * BW              # L, Total body volume (used for mass balance)
#       
#       VR <- VTot - (VSk + VI + VL + VK + VA + VP)             # L, Volume of the lumped rest compartment
#       
#       QTotc <- 0.988
#       QSk <- QSkc/QTotc * QC                # L/d, Skin 
#       QI <- QIc/QTotc * QC                  # L/d, Intestinal
#       QL <- QLc/QTotc * QC                  # L/d, Liver
#       QK <- QKc/QTotc * QC                  # L/d, Kidney
#       QA <- QAc/QTotc * QC                  # L/d, Adipose 
#       QR <- QC - (QA + QI + QK + QL + QSk)  # L/d, Rest
#       
#       QUr <- QUrc * BW              # L/d, Urine flow rate to the bladder 22 mL/kg BW/d [ICRP 89 page 161]
#       GFR <- GFR
#       QT <- QTc*GFR                 # L/d, Proximal tubule fluid flow
#       
#       tco <- tco                    # /d, Bowel residence times in the colon
#       
#       
#       ### Physicochemical ----
#       MW <- 414.07
#       
#       # Fraction unbound in plasma
#       fup <- fup #0.061/100 #
#       
#       PSk <- PSkc * fup  # Skin
#       PI <- PIc * fup    # Intestinal
#       PL <- PLc * fup    # Liver
#       PK <- PKc * fup    # Kidney
#       PA <- PAc * fup    # Adipose
#       PR <- PRc * fup    # Rest
#       
#       # Fraction unbound, calculated based on Poulin and Haddad, 2018 https://doi.org/10.1016/j.xphs.2018.03.012
#       # Equation was adapted to not account for fraction unionised
#       # OAT and OATP transporters transport the ionised compound, given that the ratio of fraction ionised at plasma to cellular pH is 1, this can be ignored (fraction unionised of PFOA is 0.9999923 at pH 7.4 and 0.9999963 at pH 7)
#       fu_PTL <- R_PTL*fup/(1 + ((R_PTL-1)*fup)) # Proximal tubule lumen
#       fu_Lic <- R_L_ec*fup/(1 + ((R_L_ec-1)*fup)) # Liver intracellular space
#       
#       
#       
#       ### Kinetic ----
#       
#       # Skin uptake
#       CL_SkBtSk <- Papp_SkB*SA_SkB*1e-3*60*60*24           # L/d, Skin barrier to skin (calculations: cm/s -> L/s * 1e-3 -> L/d *60*60*24)
#       
#       # Gastro-intestinal uptake
#       CL_IL <- (Papp_SI*SA_SI*1e-3)*60*60*24    # L/d, Intestinal lumen to intestinal tissue (calculations: cm/s -> L/s /1000 -> L/d *60*60*24) 
#       
#       VmaxOATP2B1 <- Vmax_OATP2B1c*MW*60*24*SF_OATP2B1*VIL # ug/d
#       Km_OATP2B1 <- Km_OATP2B1c*MW                         # ug/L (uM -> ug/L)
#       
#       # Liver uptake
#       Vmax_OATP1B1 <- Vmax_OATP1B1c*MW*60*24*SF_OATP1B1*VL_ec             # ug/d
#       Km_OATP1B1 <- Km_OATP1B1c*MW                                        # ug/L (uM -> ug/L)
#       Vmax_OATP1B3 <- Vmax_OATP1B3c*MW*60*24*SF_OATP1B3*VL_ec             # ug/d
#       Km_OATP1B3 <- Km_OATP1B3c*MW                                        # ug/L (uM -> ug/L)
#       
#       # Biliary excretion
#       Vmax_BSEP <- Vmax_BSEPc*MW*60*24*SF_BSEP*VL_ic         # ug/d
#       Km_BSEP <- Km_BSEPc*MW                                 # ug/L (uM -> ug/L)
#       
#       # Renal clearance
#       Vmax_OAT4 = Vmax_OAT4c*MW*60*24*SF_OAT*VPT           # ug/d (umol -> ug, min -> d)
#       Km_OAT4 = Km_OAT4c*MW                                # ug/L, scaled from uM, Louisse et al. 2024 doi.org/10.1016/j.tox.2024.153961
#       
#       # Clearance via menstruation
#       CL_menses = CL_menses     # L/d, plasma menstrual loss
#       
#       ## Dose -------------------------
#       
#       if(t<expSTOP){DoseOn=1} else{DoseOn=0}
#       
#       ## Oral exposure ##
#       DOral = expOral*BW*DoseOn         # ug, PFOA oral dose
#       OralD = DOral/Tinput*(t %% tinterval<Tinput)
#       
#       ## Dermal exposure ##
#       DDermal = expDermal*BW*DoseOn         # ug, PFOA dermal dose
#       DermalD = DDermal/Tinput*(t %% tinterval<Tinput)
#       
#       ## Concentrations -------------------------
#       
#       CSkB <- ASkB/VSkB            # ug/L, Skin barrier
#       CSk <- ASk/VSk               # ug/L, Skin
#       CVSk <- CSk/PSk              # ug/L, Skin venous 
#       
#       CIL <- AIL/VIL               # ug/L, Intestinal lumen
#       CI <- AI/VI                  # ug/L, Intestine 
#       CVI <- CI/PI                 # ug/L, Intestine venous
#       
#       CL_ec <- AL_ec/VL_ec         # ug/L, Liver extracellular
#       CL_ic <- AL_ic/VL_ic         # ug/L, Liver intracellular
#       CVL_ec <- CL_ec/PL           # ug/L, Liver extracellular venous
#       
#       CPTT <- APTT/VPTT            # ug/L, Kidney, proximal tubule tissue
#       CPTL <- APTL/VPTL            # ug/L, Kidney, proximal tubule lumen
#       CRKT <- ARKT/VRKT            # ug/L, Rest of kidney tissue
#       CRKL <- ARKL/VRKL            # ug/L, Rest of kidney lumen
#       CVRKT <- CRKT/PK             # ug/L, Rest of kidney venous
#       
#       CA <- AA/VA                  # ug/L, Adipose
#       CVA <- CA/PA                 # ug/L, Adipose venous
#       
#       CR <- AR/VR                  # ug/L, Rest
#       CVR <- CR/PR                 # ug/L, Rest venous
#       
#       CP <- AP/VP                  # ug/L, Plasma
#       
#       
#       ## Differential equations -------------------------
#       
#       dOD = OralD - OD             # ug/d, Oral dose input 
#       dDD = DermalD - DD                             # ug/d, Dermal dose input
#       
#       
#       dASkB <- DD - CL_SkBtSk*CSkB                                # ug/d, Skin barrier
#       
#       dASk <- + CL_SkBtSk*CSkB + QSk*(CP-CVSk)                   # ug/d, Skin
#       
#       
#       dAIL <- + OD - tco*AIL - CL_IL*CIL +
#         - (VmaxOATP2B1/(Km_OATP2B1 + CIL))*CIL +
#         + (Vmax_BSEP/(Km_BSEP + (CL_ic*fu_Lic)))*CL_ic*fu_Lic          # ug/d, Intestine lumen
#       
#       dAI <- QI*(CP - CVI) + CL_IL*CIL + 
#         + (VmaxOATP2B1/(Km_OATP2B1 + CIL))*CIL                     # ug/d, Intestinal
#       
#       dAFe <-  tco*AIL                                       # ug/d, Feces
#       
#       
#       dAL_ec <- + QI*CVI + QL*CP - (QI+QL)*CVL_ec + 
#         - (Vmax_OATP1B1/(Km_OATP1B1 + (CL_ec*fup)))*CL_ec*fup +
#         - (Vmax_OATP1B3/(Km_OATP1B3 + (CL_ec*fup)))*CL_ec*fup            # ug/d, Liver extracellular space (vascular + interstitial space)
#       
#       dAL_ic <- (Vmax_OATP1B1/(Km_OATP1B1 + (CL_ec*fup)))*CL_ec*fup +
#         + (Vmax_OATP1B3/(Km_OATP1B3 + (CL_ec*fup)))*CL_ec*fup +
#         - (Vmax_BSEP/(Km_BSEP + (CL_ic*fu_Lic)))*CL_ic*fu_Lic                       # ug/d, Liver intracellular space
#       
#       
#       dAPTT <- QK*(CP - CPTT) + 
#         + (Vmax_OAT4/(Km_OAT4+(CPTL*fu_PTL)))*CPTL*fu_PTL        # ug/d, Proximal tubule tissue 
#       
#       dAPTL <- + fup*GFR*CP - QT*CPTL +
#         - (Vmax_OAT4/(Km_OAT4+(CPTL*fu_PTL)))*CPTL*fu_PTL        # ug/d, Proximal tubule lumen    
#       
#       dARKT <- QK*(CPTT - CVRKT)                               # ug/d, Rest of kidney
#       
#       dARKL <- QT*CPTL - QUr*CRKL                              # ug/d, Rest of kidney lumen
#       
#       dAUr <- QUr*CRKL                                         # ug/d, Urine
#       
#       
#       dAA <- QA*(CP-CVA)                                      # ug/d, Adipose
#       
#       
#       dAR <- QR*(CP-CVR)                                      # ug/d, Rest
#       
#       dAP <- - (QSk + QI + QL + QA + QR + QK)*CP - fup*GFR*CP +     # ug/d, Arterial Plasma
#         + (QL+QI)*CVL_ec + QSk*CVSk + QK*CVRKT + QA*CVA + QR*CVR +    # ug/d, Venous Plasma
#       - CL_menses*CP # ug/d, Venous Plasma
#       
#       dAMp <- + CL_menses*CP
#       
#       # Mass Balance
#       Atot <- OD + DD + 
#         ASk +
#         AIL + AI + AFe + 
#         AL_ec + AL_ic +
#         APTT + APTL + ARKT + ARKL + AUr +
#         AA + 
#         AR +
#         AMp +
#         AP  
#       
#       dAin <- OralD + DermalD # to be used if repeated exposure
#       MB <- Ain - Atot + 1    # to be used if repeated exposure
#       # End
#       
#       list(c(dOD,
#              dDD,
#              dASkB,
#              dASk, 
#              dAIL,
#              dAI, 
#              dAFe,
#              dAL_ec,
#              dAL_ic,
#              dAPTT,
#              dAPTL,
#              dARKT,
#              dARKL,
#              dAUr,
#              dAA,
#              dAR, 
#              dAP, 
#              dAMp,
#              dAin
#       ), 
#       c(CSkB = CSkB, 
#         CSk = CSk, 
#         CIL = CIL,
#         CI = CI, CVI = CVI, 
#         CL_ec = CL_ec,
#         CL_ic = CL_ic,
#         CVL_ec = CVL_ec,
#         CPTT = CPTT,
#         CPTL = CPTL,
#         CRKT = CRKT,
#         CRKL = CRKL,
#         CA = CA, CVA = CVA,
#         CR = CR, CVR = CVR,
#         CP = CP,
#         Atot = Atot, 
#         MB = MB
#       )
#       )
#     })
#   }
#   
#   output_PBK <-lsoda(y = y, 
#                       times = times, 
#                       func = PBK.model,
#                       parms = parms)
#   
#   return(as.data.frame(output_PBK))
# }
# 
# 
# # INHALATION ####
# INHALATION_PBK_RUN <- function(y, parms, times){ # Input for ode
#   
#   
#   PBK.model <- function(t, state, parameters){
#     with(as.list(c(state, parameters)), {
#       
#       ### Physiological ----
#       
#       VIL <- VILc * BW                # L, Volume of intestinal lumen
#       VI <- VIc * BW                  # L, Volume of intestine
#       SA_SI <- SA_SIc * BW
#       
#       VL <- VLc * BW                  # L, Volume of liver
#       VL_ic <- VL_icc*VL              # L, Volume liver intracellular                 
#       VL_ec <- VL_ecc*VL              # L, Volume liver extracellular
#       
#       VK <- VKc*BW                    # L, Volume of kidney
#       VPT <- VPTc*VK                  # L, Volume of proximal tubule
#       VPTT <- VPT*VPTTc               # L, Volume of proximal tubule tissue
#       VPTL <- VPT*VPTLc               # L, Volume of proximal tubule lumen
#       VRKL <- (VK - VPT)*VRKLc        # L, Volume of rest of kidney lumen
#       VRKT <- VK - VPTT - VPTL - VRKL # L, Volume of rest of kidney tissue, used in the model
#       
#       VA <- VAc * BW/0.9 + MAc/0.9    # L, Volume of adipose
#       
#       VLu <- VLuc * BW                # L, Volume of lungs
#       
#       Fr_art_blood = 0.0178 / (0.0178 + 0.0533) #fraction of arterial blood (corrected for plasma), ref. A. Ratier et al. 2024 doi: 10.1016/j.envint.2024.108621
#       VP <- VPc * (1-Hct) * BW                  # L, Volume of plasma
#       VAP <- VP*Fr_art_blood         # L, Volume of arterial plasma
#       VVP <- VP - VAP                 # L, Volume of venous plasma
#       
#       VTotc <- 0.96
#       VTot <- VTotc * BW              # L, Total body volume (used for mass balance)
#       
#       VR <- VTot - (VI + VL + VK + VA + VLu + VAP + VVP)             # L, Volume of the lumped rest compartment
#       
#       QTotc <- 0.988
#       QI <- QIc/QTotc * QC                  # L/d, Intestinal
#       QL <- QLc/QTotc * QC                  # L/d, Liver
#       QK <- QKc/QTotc * QC                  # L/d, Kidney
#       QA <- QAc/QTotc * QC                  # L/d, Adipose 
#       QLu <- QC                             # L/d, Lungs
#       QR <- QC - (QA + QI + QK + QL)        # L/d, Rest
#       
#       QUr <- QUrc * BW              # L/d, Urine flow rate to the bladder 22 mL/kg BW/d [ICRP 89 page 161]
#       GFR <- GFR
#       QT <- QTc*GFR                      # L/d, Proximal tubule fluid flow
#       
#       tco <- tco                    # /d, Bowel residence times in the colon
#       
#       
#       ### Physicochemical ----
#       MW <- 414.07
#       
#       # Fraction unbound in plasma
#       fup <- fup # 0.061/100 #
#       
#       PI <- PIc * fup     # Intestinal
#       PL <- PLc * fup     # Liver
#       PK <- PKc * fup     # Kidney
#       PA <- PAc * fup     # Adipose
#       PLu <- PLuc * fup   # Lungs
#       PR <- PRc * fup     # Rest
#       
#       # Fraction unbound, calculated based on Poulin and Haddad, 2018 https://doi.org/10.1016/j.xphs.2018.03.012
#       # Equation was adapted to not account for fraction unionised
#       # OAT and OATP transporters transport the ionised compound, given that the ratio of fraction ionised at plasma to cellular pH is 1, this can be ignored (fraction unionised of PFOA is 0.9999923 at pH 7.4 and 0.9999963 at pH 7)
#       fu_PTL <- R_PTL*fup/(1 + ((R_PTL-1)*fup)) # Proximal tubule lumen
#       fu_Lic <- R_L_ec*fup/(1 + ((R_L_ec-1)*fup)) # Liver intracellular space
#       
#       
#       ### Kinetic ----
#       
#       # Gastro-intestinal uptake
#       CL_IL <- (Papp_SI*SA_SI*1e-3)*60*60*24    # L/d, Intestinal lumen to intestinal tissue (calculations: cm/s -> L/s /1000 -> L/d *60*60*24) 
#       
#       VmaxOATP2B1 <- Vmax_OATP2B1c*MW*60*24*SF_OATP2B1*VIL # ug/d
#       Km_OATP2B1 <- Km_OATP2B1c*MW                         # ug/L (uM -> ug/L)
#       
#       # Liver uptake
#       Vmax_OATP1B1 <- Vmax_OATP1B1c*MW*60*24*SF_OATP1B1*VL_ec             # ug/d
#       Km_OATP1B1 <- Km_OATP1B1c*MW                                        # ug/L (uM -> ug/L)
#       Vmax_OATP1B3 <- Vmax_OATP1B3c*MW*60*24*SF_OATP1B3*VL_ec             # ug/d
#       Km_OATP1B3 <- Km_OATP1B3c*MW                                        # ug/L (uM -> ug/L)
#       
#       # Biliary excretion
#       Vmax_BSEP <- Vmax_BSEPc*MW*60*24*SF_BSEP*VL_ic         # ug/d
#       Km_BSEP <- Km_BSEPc*MW                                 # ug/L (uM -> ug/L)
#       
#       # Renal clearance
#       Vmax_OAT4 = Vmax_OAT4c*MW*60*24*SF_OAT*VPT           # ug/d (umol -> ug, min -> d)
#       Km_OAT4 = Km_OAT4c*MW                                # ug/L (uM -> ug/L)
#       
#       # Clearance via menstruation
#       CL_menses = CL_menses     # L/d, plasma menstrual loss
#       
#       
#       ## Dose -------------------------
#       
#       if(t<expSTOP){DoseOn=1} else{DoseOn=0}
#       
#       ## Inhalation exposure ##
#       DLung = expLung*BW*DoseOn         # ug, PFOA inhalation dose
#       LungD = DLung/Tinput*(t %% tinterval<Tinput)
#       
#       ## Concentrations -------------------------
#       
#       CIL <- AIL/VIL               # ug/L, Intestinal lumen
#       CI <- AI/VI                  # ug/L, Intestine 
#       CVI <- CI/PI                 # ug/L, Intestine venous
#       
#       CL_ec <- AL_ec/VL_ec         # ug/L, Liver extracellular
#       CL_ic <- AL_ic/VL_ic         # ug/L, Liver intracellular
#       CVL_ec <- CL_ec/PL           # ug/L, Liver extracellular venous
#       
#       CPTT <- APTT/VPTT            # ug/L, Kidney, proximal tubule tissue
#       CPTL <- APTL/VPTL            # ug/L, Kidney, proximal tubule lumen
#       CRKT <- ARKT/VRKT            # ug/L, Rest of kidney tissue
#       CRKL <- ARKL/VRKL            # ug/L, Rest of kidney lumen
#       CVRKT <- CRKT/PK             # ug/L, Rest of kidney venous
#       
#       CA <- AA/VA                  # ug/L, Adipose
#       CVA <- CA/PA                 # ug/L, Adipose venous
#       
#       CR <- AR/VR                  # ug/L, Rest
#       CVR <- CR/PR                 # ug/L, Rest venous
#       
#       CLu <- ALu/VLu               # ug/L, Lungs
#       CVLu <- CLu/PLu              # ug/L, Lungs venous
#       
#       CAP <- AAP/VAP               # ug/L, Arterial plasma
#       CVP <- AVP/VVP               # ug/L, Venous plasma
#       
#       
#       ## Differential equations -------------------------
#       
#       dLuD = LungD - LuD             # ug/d, Inhaled dose input 
#       
#       dALu <- LuD + QC*(CVP - CVLu)                            # ug/d, Lungs
#       
#       
#       dAIL <- - tco*AIL - CL_IL*CIL +
#         - (VmaxOATP2B1/(Km_OATP2B1 + CIL))*CIL +
#         + (Vmax_BSEP/(Km_BSEP + (CL_ic*fu_Lic)))*CL_ic*fu_Lic          # ug/d, Intestine lumen
#       
#       dAI <- QI*(CP - CVI) + CL_IL*CIL + 
#         + (VmaxOATP2B1/(Km_OATP2B1 + CIL))*CIL                     # ug/d, Intestinal
#       
#       
#       dAFe <-  tco*AIL                                       # ug/d, Feces
#       
#       
#       dAL_ec <- + QI*CVI + QL*CAP - (QI+QL)*CVL_ec + 
#         - (Vmax_OATP1B1/(Km_OATP1B1 + (CL_ec*fup)))*CL_ec*fup +
#         - (Vmax_OATP1B3/(Km_OATP1B3 + (CL_ec*fup)))*CL_ec*fup            # ug/d, Liver extracellular space (vascular + interstitial space)
#       
#       dAL_ic <- (Vmax_OATP1B1/(Km_OATP1B1 + (CL_ec*fup)))*CL_ec*fup +
#         + (Vmax_OATP1B3/(Km_OATP1B3 + (CL_ec*fup)))*CL_ec*fup +
#         - (Vmax_BSEP/(Km_BSEP + (CL_ic*fu_Lic)))*CL_ic*fu_Lic             # ug/d, Liver intracellular space
#       
#       
#       dAPTT <- QK*(CAP - CPTT) + 
#         + (Vmax_OAT4/(Km_OAT4+(CPTL*fu_PTL)))*CPTL*fu_PTL      # ug/d, Proximal tubule tissue 
#       
#       dAPTL <- + fup*GFR*CAP - QT*CPTL +
#         - (Vmax_OAT4/(Km_OAT4+(CPTL*fu_PTL)))*CPTL*fu_PTL      # ug/d, Proximal tubule lumen    
#       
#       dARKT <- QK*(CPTT - CVRKT)                               # ug/d, Rest of kidney
#       
#       dARKL <- QT*CPTL - QUr*CRKL                              # ug/d, Rest of kidney lumen
#       
#       dAUr <- QUr*CRKL                                         # ug/d, Urine
#       
#       
#       dAA <- QA*(CAP-CVA)                                      # ug/d, Adipose
#       
#       
#       dAR <- QR*(CAP-CVR)                                      # ug/d, Rest
#       
#       
#       dAAP <- - (QI + QL + QA + QR + QK)*CAP + QC*CVLu - fup*GFR*CAP    #+ QSt       # ug/d, Arterial Plasma
#       dAVP <- (QL+QI)*CVL_ec + QK*CVRKT + QA*CVA + QR*CVR - QC*CVP +      # ug/d, Venous Plasma
#         - CL_menses*CP # ug/d, Venous Plasma
#       
#       dAMp <- + CL_menses*CP
#       
#       # Mass Balance
#       Atot <- LuD +
#         ALu +
#         AIL + AI + AFe + 
#         AL_ec + AL_ic +
#         APTT + APTL + ARKT + ARKL + AUr +
#         AA + 
#         AR +
#         AAP + AVP +
#         AMp
#       
#       dAin <- LungD 
#       MB <- Ain - Atot + 1    
#       
#       
#       # End
#       
#       list(c(dLuD,
#              dALu, 
#              dAIL,
#              dAI, 
#              dAFe,
#              dAL_ec,
#              dAL_ic,
#              dAPTT,
#              dAPTL,
#              dARKT,
#              dARKL,
#              dAUr,
#              dAA,
#              dAR, 
#              dAAP,
#              dAVP,
#              dAMp,
#              dAin
#       ), 
#       c(CLu = CLu,
#         CIL = CIL,
#         CI = CI, CVI = CVI, 
#         CL_ec = CL_ec,
#         CL_ic = CL_ic,
#         CVL_ec = CVL_ec,
#         CPTT = CPTT,
#         CPTL = CPTL,
#         CRKT = CRKT,
#         CRKL = CRKL,
#         CA = CA, CVA = CVA,
#         CR = CR, CVR = CVR,
#         CP = CAP + CVP,
#         Atot = Atot, 
#         MB = MB
#       )
#       )
#     })
#   }
#   
#   output_PBK <-lsoda(y = y, 
#                       times = times, 
#                       func = PBK.model,
#                       parms = parms)
#   
#   return(as.data.frame(output_PBK))
# }
# 
# 
