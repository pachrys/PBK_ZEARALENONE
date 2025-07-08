## PBK model rat ####
PBK.model <- function(t, state, parms){
  with(as.list(c(state, parms)), {
    
    ## Calculated parameters ####
    
    ## Physiological constants ####
    # Fractional Tissue volumes (reference: Punt et al. (2020))
    BW = BW  ##{Kg} #  body weight rat
    
    VST = FV_ST*BW  # fraction volume of small intestine
    VSIL = FV_SIL*BW # fractional volume of small intestinal lumen
    VREL = FV_REL*FV_LIL*BW # fractional volume of rectum lumen from total LI volume (only the lower part goes to the systemic circulation)
    VLIL = (FV_LIL*BW)-VREL
    
    VSI = FV_SI*BW # fractional volume of small intestinal tissue
    VL = FV_L*BW  # fractional of liver tissue
    VR = FV_R*BW  # 0.082-VSIc-VLc # fraction of rapidly perfused tissue
    VS = FV_S*BW  # 0.737-VFc # fraction of slowly perfused tissue
    VF = FV_F*BW  # fraction of fat tissue
    VB = FV_B*BW  # fraction of blood
    
    # Surface area 
    SASI = SASI # surface area (dm2)
    SARE = FV_REL*SALI
    SALI = SALI- SARE  # surface area of large intestine (dm2) reference: Vdoviaková et al. (2016)
    
    # Blood flow rates (reference: Punt et al. (2020))
    QC = QC  ##{L/h} # cardiac output
    QSI = QSIc*QC  # fraction of blood flow to small intestine
    QL = QLc*QC  # 0.242 - QSIc # fraction of blood flow to liver 101 
    QR = QRc*QC  # 0.76 - QSIc - QLc # fraction of blood flow to rapidly perfused tissue
    QS = QSc*QC  # 0.24 - QFc # fraction of blood flow to slowly perfused tissue
    QF = QFc*QC  # fraction of blood flow to fat
    
    # Intestinal fluid flow rates (reference: Punt et al. (2020))
    GER = GER # stomach emptying rate (/h)
    SITR = SITR # small intestinal transit time (/h)
    LITR = LITR # large intestinal transit time (/h)
    
    # Glomerular filtration rate: Punt et al. (2020)
    GFR = GFR #(/h)
    
    ## Kinetic parameters ####
    ### absorption/transfer rates ####
    CLSI_ZEN = PappCaco2_ZEN * SASI * 1e-3 # (L/h) cm/h to L/h reference: Pfeiffer et al. (2011)
    CLLI_ZEN = PappCaco2_ZEN * SALI * 1e-3 # (L/h) cm/h to L/h reference: Pfeiffer et al. (2011)
    CLRE_ZEN = PappCaco2_ZEN * SARE * 1e-3 # (L/h) cm/h to L/h reference: Pfeiffer et al. (2011)
    
    CLLI_aZEL = PappCaco2aZEL * SALI * 1e-3 # (L/h) cm/h to L/h reference: Pfeiffer et al. (2011)
    
    
    ### metabolism in small intestine tissue ####
    VmaxSI_ZENGLU <- VmaxSI_ZENGLU * S9SI * VSI * 1000 * 60/1000 #umol/h, scaling from nmol/min/mg protein to nmol/min for whole liver; equation from {Grandoni, 2019 #457}; in vivo = in vitro * scaling to microsomal protein per gram of rat liver * scaling to total weight (assuming density=1) of rat intestine; SF = 15.5mg microsomal protein/g gut ({Musther, 2017 #458})
    
    S9SI =  57 #37.1  # small intestinal S9 protein yield (mg S9 protein/gram intestine) reference: Peters et al. (2016); sum of microsomal 3.9 mg/g and cytosolic 18mg/g https://sci-hub.box/https://doi.org/10.1007/s13318-020-00666-w
    VmaxSI_ZENGLU <- VmaxSIZENGLUc/1000*60* S9SI *SI*BW / 5
    
    VmaxSI_ZENGLUc =  8.617  # in vitro  S9 incubations from this study
    KmSI_ZENGLU = 6.198  # in vitro  S9 incubations from this study
    
    ### metabolism in large intestine lumen (microbiota compartment) ####
    VMB = 0.0164  # fraction of feces of BW reference: Hoskins and Zamcheck (1968)
    VmaxLI_aZELc = 0.0138  # reference: Mendez-Catala et al. (2020)
    KmLI_aZEL = 66  # reference: Mendez-Catala et al. (2020)
    VmaxLI_bZELc =  0.0084  # reference: Mendez-Catala et al. (2020)
    KmLI_bZEL = 80  # reference: Mendez-Catala et al. (2020)
    
    ### metabolism in liver ####
    VLS9 <- 125 #  liver S9 protein yield, (mg S9 protein/gram liver) references: sum of microsomal 38mg/g liver and cytosolic 87mg/g liver content, references https://doi.org/10.2174/138920007779315053 and 10.1016/j.taap.2011.03.020, https://sci-hub.box/https://doi.org/10.1007/s13318-020-00666-w
    VmaxL_aZELc = 32  # reference: Malekinejad et al. (2006)
    KmL_aZEL = 592  # reference: Malekinejad et al. (2006)
    VmaxL_bZELc = 72  # reference: Malekinejad et al. (2006)
    KmL_bZEL = 21  # reference: Malekinejad et al. (2006)
    
    VmaxL_ZENGLUc = 7.029  # in vitro  S9 incubations from this study
    KmL_ZENGLU = 6.753  # in vitro  S9 incubations from this study
    
    VmaxL_aZELGLUc =  6.963  # in vitro  S9 incubations from this study
    KmL_aZELGLU  = 7.426 # in vitro  S9 incubations from this study
    
    ### active transport in small intestine tissue ####
    VmaxUptake_ZENGLUC = 40.1
    KmUptake_ZENGLUC = 23.6
    
 
    
    # metabolites α-ZEL and ß-ZEL scaled maximum rates of metabolism (μmol/h)
    VmaxLI_aZEL <-  VmaxLIaZELc/1000000*1000 *VMB*BW
    VmaxLI_bZEL <-  VmaxLIbZELc/1000000*1000 *VMB*BW
    
    # liver metabolism
    L <- VLc*1000  # (gram/kg BW) # liver
    # metabolites α-ZEL and ß-ZEL scaled maximum rates of metabolism (μmol/h)
    VmaxL_aZEL <- VmaxLaZELc/1000000*60* VLS9 *L*BW
    VmaxL_bZEL <-  VmaxLbZELc/1000000*60* VLS9 *L*BW
    
    # metabolites ZENGLU scaled maximum rates of metabolism (μmol/h)
    VmaxL_ZENGLU <- VmaxLZENGLUc /1000*60* VLS9 *L*BW / 5
    
    # metabolites α-ZELGLU scaled maximum rates of metabolism (μmol/h)
    VmaxL_aZELGLU <- VmaxLaZELGLUc/1000*60* VLS9 *L*BW
    
    # my additions
    fup <- 0.156
    
    
    
    ## Concentrations #####
    
    
    # Zearalenone 
    CSI_ZEN <- ASI_ZEN/VSI
    CSIL_ZEN <- ASIL_ZEN/VSI
    CLI_ZEN <- ALI_ZEN/(VMB*BW) 
    CL_ZEN <- AL_ZEN/VL
    CF_ZEN <- AF_ZEN/VF
    CR_ZEN <- AR_ZEN/VR
    CS_ZEN <- AS_ZEN/VS
    CB_ZEN <- AB_ZEN/VB
    
    CVSI_ZEN <- CSI_ZEN/PIZEN
    CVL_ZEN <- CL_ZEN/PLZEN
    CVF_ZEN <- CF_ZEN/PFZEN
    CVR_ZEN <- CR_ZEN/PRZEN
    CVS_ZEN <- CS_ZEN/PSZEN
    
    
    # α-Zearalenol
    CLI_aZEL <- ALI_aZEL/VLI #(VMB*BW)
    CL_aZEL <- AL_aZEL/VL
    CF_aZEL <- AF_aZEL/VF
    CR_aZEL <- AR_aZEL/VR
    CS_aZEL <- AS_aZEL/VS
    CB_aZEL <- AB_aZEL/VB
    
    CVL_aZEL <- CL_aZEL/PLaZEL
    CVF_aZEL <- CF_aZEL/PFaZEL
    CVR_aZEL <- CR_aZEL/PRaZEL
    CVS_aZEL <- CS_aZEL/PSaZEL
    
    
    # β-Zearalenol
    CLI_bZEL <- ALI_bZEL/VL
    CL_bZEL <- AL_bZEL/VL
    
    
    # Zearalenone-glucuronide
    CLI_ZENGLU <- ALI_ZENGLU/VLI #(VMB*BW)
    CSI_ZENGLU <- ASI_ZENGLU/VSI        
    CL_ZENGLU <- AL_ZENGLU/VL
    CF_ZENGLU <- AF_ZENGLU/VF
    CR_ZENGLU <- AR_ZENGLU/VR
    CS_ZENGLU <- AS_ZENGLU/VS
    CB_ZENGLU <- AB_ZENGLU/VB
    
    CVL_ZENGLU <- CL_ZENGLU/PLZENGLU
    CVF_ZENGLU <- CF_ZENGLU/PFZENGLU
    CVR_ZENGLU <- CR_ZENGLU/PRZENGLU
    CVS_ZENGLU <- CS_ZENGLU/PSZENGLU
    
    # a-Zearalenol glucuronide
    CL_aZELGLU <- AL_aZELGLU/(VMB*BW)   # liver a zel gluc
    CLI_bZEL <- ALI_bZEL/VL             # liver b zel gluc
    CB_ZENGLU <- AB_ZENGLU/VB
    
    
    
    ## Dynamic ####
    
    # Stomach (in) 
    dAST_ZEN <- - GER*AST_ZEN # amount in stomach #ODOSEumol
    
    # Fat 
    dAF_ZEN <- QF*(CB_ZEN-CVF_ZEN) # amount of ZEN in fat tissue, (μmol)
    dAF_aZEL <- QF*(CB_aZEL-CVF_aZEL) # amount of α-ZEL in fat tissue, (μmol)
    dAF_ZENGLUC <- QF*(CB_ZENGLUC-CVF_ZENGLUC) # amount of ZEN-gluc in fat tissue, (μmol)
    
    # Rapidly Perfused 
    dAR_ZEN <- QR*(CB_ZEN-CVR_ZEN) # amount of ZEN in rapidly perfused tissue, (μmol)
    dAR_aZEL <- QR*(CB_aZEL-CVR_aZEL) # amount of α-ZEL in rapidly perfused tissue, (μmol)
    dAR_ZENGLUC <- QR*(CB_ZENGLUC-CVR_ZENGLUC) # amount of α-ZEL in rapidly perfused tissue, (μmol)
    
    # Slowly Perfused 
    dAS_ZEN <- QS*(CB_ZEN-CVS_ZEN) # amount of ZEN in slowly perfused tissue, (μmol)
    dAS_aZEL <- QS*(CB_aZEL-CVS_aZEL) # amount of α-ZEL in slowly perfused tissue, (μmol)
    dAS_ZENGLUC <- QS*(CB_ZENGLUC-CVS_ZENGLUC) # amount of α-ZEL in slowly perfused tissue, (μmol)
    
    # Small intestinal lumen 
    dASI_ZEN <- GER*AST_ZEN - kin*ASI_ZEN - kabin*CSI_ZEN # small intestine lumen compartment (μmol)
    dASI_ZENGLUC <- + 0.5*AL_ZENGLU - kin*ASI_ZENGLUC 
    
    # Large intestinal lumen 
    dALI_ZEN <- kin*ASIL_ZEN + 
      - VmaxLI_aZEL*CLI_ZEN / (KmLIaZEL + CLI_ZEN) +
      - VmaxLI_bZEL*CLI_ZEN / (KmLIbZEL + CLI_ZEN) +
      - kco*ALI_ZEN - KfeZ*ALI_ZEN +
      kin*ASIL_ZENGLUC  # amount of a-ZEL in large intestine lumen, (μmol)
    
    dALI_aZEL <- VmaxLI_aZEL*CLI_ZEN / (KmLI_aZEL + CLI_ZEN) +
      - kabin8*(ALI_aZEL/VMB*BW) - KfeaZ*ALI_aZEL # amount of α-ZEL formed due to gut microbial activity, (μmol)
    
    dALI_bZEL <- VmaxLI_bZEL*CLI_ZEN / (KmLI_bZEL + CLI_ZEN) # amount of b-ZEL formed due to gut microbial activity, (μmol)
    
    # Faeces ####
    dAFE_ZEN <- KfeZ*ALI_ZEN
    dAFE_aZEL <- KfeaZ*ALI_aZEL
    
    
    # Small intestinal tissue ####
    dASIL_ZEN <- kabin*CSIL_ZEN + QSI*(CB_ZEN-CVSI_ZEN) +
      - VmaxSI_ZENGLU*CVSI_ZEN*fup / (KmSI_ZENGLU + CVSI_ZEN*fup)
     # + (ASI_ZENGLU) # amount of ZEN in small intestinal tissue, (μmol)
    
    dASIL_ZENGLU <- VmaxSI_ZENGLU*CVSI_ZEN*fup / (KmSI_ZENGLU + CVSI_ZEN*fup) +
      - VmaxUptake_ZENGLU*CVSI_ZEN*fup/(KmUptake_ZENGLU + CVSI_ZEN*fup) +
      + QSI*(CB_ZENGLU-CVSI_ZENGLU)
    #+ ALZENGLU - (ASIZENGLU) # amount of ZEN metabolized to metabolite ZENGLU, (μmol)
    
    
    # Liver 
    dAL_ZEN <- QSI*CVSI_ZEN + QL*CB_ZEN - (QSI+QL)*CVL_ZEN +
      - VmaxL_aZEL*CVL_ZEN*fup / (KmLaZEL + CVL_ZEN*fup) +
      - VmaxL_bZEL*CVL_ZEN*fup / (KmLbZEL + CVL_ZEN*fup) +
      - VmaxL_ZENGLU*CVL_ZEN*fup / (KmLZENGLU + CVL_ZEN*fup) #kco*ALI_ZEN  # amount of ZEN in liver, (μmol)
    
    dAL_ZENGLU <- VmaxL_ZENGLU*CVL_ZEN*fup / (KmLZENGLU + CVL_ZEN*fup) +
      - 0.5*AL_ZENGLU +
      - (QSI+QL)*CVL_ZENGLU + QL*CB_ZENGLU + QSI*CVSI_ZENGLU +
      + VmaxUptake_ZENGLU*CVSI_ZEN*fup/(KmUptake_ZENGLU + CVSI_ZEN*fup) # amount of ZEN metabolized to metabolite ZENGLU in liver, (μmol)
    
    dAL_bZEL <- VmaxL_bZEL*CVL_ZEN*fup / (KmLbZEL + CVL_ZEN*fup) # amount of ZEN metabolized to metabolite ß-ZEL in liver, (μmol)
    
    dAL_aZEL <- VmaxL_aZEL*CVL_ZEN*fup / (KmLaZEL + CVL_ZEN*fup) +
      + kabin8*(ALI_aZEL/VMB*BW) +
      - VmaxL_aZELGLU*CVL_aZEL*fup / (KmLaZELGLU + CVL_aZEL*fup) +
      - (QL+QSI)*CVL_aZEL + (QL+QSI)*CB_aZEL # amount of ZEN metabolized to metabolite α-ZEL in liver, (μmol)
    
    dAL_aZELGLU <- VmaxL_aZELGLU*CVL_aZEL*fup / (KmLaZELGLU + CVL_aZEL*fup) # amount of α-ZEL glucuronide in liver
    
    # Blood
    dAB_ZEN <- - (QSI + QL + QR + QS + QF)*CB_ZEN + 
      + (QSI+QL)*CVL_ZEN + QF*CVF_ZEN + QR*CVR_ZEN + QS*CVS_ZEN +
      - GFR*AB_ZEN*fup +
      + kco*ALI_ZEN  # amount of ZEN in blood, (μmol) IVDOSEumol
    
    dAB_aZEL <- - (QSI + QL + QR + QS + QF)*CB_aZEL +
      + (QL+QSI)*CVL_aZEL + QF*CVF_aZEL + QR*CVR_aZEL + QS*CVS_aZEL +
      - GFR*AB_aZEL*fup # amount of α-ZELin blood, (μmol)
    
    dAB_ZENGLU <- - (QSI + QL + QR + QS + QF)*CB_ZENGLU + 
      + (QSI+QL)*CVL_ZENGLUC + QF*CVF_ZENGLUC + QR*CVR_ZENGLUC + QS*CVS_ZENGLUC +
      - GFR*AB_ZENGLUC*fup 
    
    # Urine
    dAUR_ZEN <- GFR*AB_ZEN*fup # urinary excretion
    dAUR_aZEL <- GFR*AB_aZEL*fup
    dAUR_ZENGLU <- GFR*CB_ZENGLU*fup
    
    
    #mass balance calculation
    dATotal <- AST_ZEN #ODOSEumol + IVDOSEumol
    ACalculated <-  + AF_ZEN + AF_aZEL + AF_ZENGLUC +
      + AR_ZEN + AR_aZEL + AR_ZENGLUC + 
      + AS_ZEN + AS_aZEL + AS_ZENGLUC + 
      + ASIL_ZEN + ASIL_ZENGLUC + 
      + ALI_ZEN + ALI_aZEL + ALI_bZEL + 
      + AFE_ZEN + AFE_aZEL + 
      + ASI_ZEN + ASI_ZENGLU +
      + AL_ZEN + AL_ZENGLU + AL_bZEL + AL_aZEL + AL_aZELGLU +
      + AB_ZEN + AB_aZEL + AB_ZENGLU + 
      + AUR_ZEN + AUR_aZEL + AUR_ZENGLU
    
    
    MB <- ATotal - ACalculated + 1
    
    # End
    # 
    list(c(
      AF_ZEN, AF_aZEL, AF_ZENGLUC,
      AR_ZEN, AR_aZEL, AR_ZENGLUC,  
      AS_ZEN, AS_aZEL, AS_ZENGLUC,  
      ASIL_ZEN, ASIL_ZENGLUC,  
      ALI_ZEN, ALI_aZEL,ALI_bZEL, 
      AFE_ZEN, AFE_aZEL, 
      ASI_ZEN, ASI_ZENGLU,
      AL_ZEN, AL_ZENGLU, AL_bZEL, AL_aZEL, AL_aZELGLU, 
      AB_ZEN, AB_aZEL, AB_ZENGLU, 
      AUR_ZEN, AUR_aZEL, AUR_ZENGLU,
      dATotal
    ),
    c(ACalculated = ACalculated,
      MB = MB,
      CSI_ZEN = CSI_ZEN,
      CSIL_ZEN = CSIL_ZEN,
      CL_ZEN = CL_ZEN,
      CF_ZEN = CF_ZEN,
      CR_ZEN = CR_ZEN,
      CS_ZEN = CS_ZEN,
      CB_ZEN = CB_ZEN,
      CVSI_ZEN = CVSI_ZEN,
      CVL_ZEN = CVL_ZEN,
      CVF_ZEN = CVF_ZEN,
      CVR_ZEN = CVR_ZEN,
      CVS_ZEN = CVS_ZEN,
      CLI_aZEL = CLI_aZEL,
      CL_aZEL = CL_aZEL,
      CF_aZEL = CF_aZEL,
      CR_aZEL = CR_aZEL,
      CS_aZEL = CS_aZEL,
      CB_aZEL = CB_aZEL,
      CVL_aZEL = CVL_aZEL,
      CVF_aZEL = CVF_aZEL,
      CVR_aZEL = CVR_aZEL,
      CVS_aZEL = CVS_aZEL,
      CLI_bZEL = CLI_bZEL,
      CL_bZEL = CL_bZEL,
      CLI_ZENGLU = CLI_ZENGLU,
      CSIL_ZENGLU = CSIL_ZENGLU,  
      CL_ZENGLU = CL_ZENGLU,
      CF_ZENGLU = CF_ZENGLU,
      CR_ZENGLU = CR_ZENGLU,
      CS_ZENGLU = CS_ZENGLU,
      CB_ZENGLU = CB_ZENGLU,
      CVL_ZENGLU = CVL_ZENGLU,
      CVF_ZENGLU = CVF_ZENGLU,
      CVR_ZENGLU = CVR_ZENGLU,
      CVS_ZENGLU = CVS_ZENGLU,
      CL_aZELGLU = CL_aZELGLU,
      CLI_bZEL = CLI_bZEL,
      CB_ZENGLU = CB_ZENGLU
    ))

  })
}