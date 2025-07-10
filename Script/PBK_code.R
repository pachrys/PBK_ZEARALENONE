## PBK model rat ####
PBK.model <- function(t, state, parms){
  with(as.list(c(state, parms)), {
    
    ## Calculated parameters ####
    
    ### Physiological parameters ####
    VSI <- VSIc*BW # (L or Kg) # volume of small intestine tissue (calculated)
    VL <- VLc*BW # (L or Kg) # volume of liver tissue (calculated)
    VR <- VRc*BW # (L or Kg) # volume of rapidly perfused tissue (calculated)
    VS <- VSc*BW # (L or Kg) # volume of slowly perfused tissue (calculated)
    VF <- VFc*BW # (L or Kg) # volume of fat tissue (calculated)
    VB <- VBc*BW # (L or Kg) # volume of blood (calculated)
    
    QSI <- QSIc*QC # (L/h) # blood flow to small intestine tissue (calculated)
    QL <- QLc*QC # (L/h) # blood flow to liver tissue (calculated)
    QR <- QRc*QC # (L/h) # blood flow to rapidly perfused tissue (calculated)
    QS <- QSc*QC # (L/h) # blood flow to slowly perfused tissue (calculated)
    QF <- QFc*QC # (L/h) # blood flow to fat tissue (calculated)
    
    
    # Intestinal sub compartments
    Vin1 <- Vin # volume of intestine compartment 1 (L)
    SASIL1 <- SAin # surface area of intestine compartment 1 (dm2)
    kin1 <- kin # transfer rate to intestine compartment 2 (/hr)
    
    Vin2 <- Vin # volume of intestine compartment 2 (L)
    SASIL2 <- SAin # surface area of intestine compartment 2 (dm2)
    kin2 <- kin # transfer rate to intestine compartment 3 (/hr)
    
    Vin3 <- Vin # volume of intestine compartment 3 (L)
    SASIL3 <-  SAin # surface area of intestine compartment 3 (dm2)
    kin3 <- kin # transfer rate to intestine compartment 4 (/hr)
    
    Vin4 <- Vin # volume of intestine compartment 4 (L)
    SASIL4 <- SAin # surface area of intestine compartment 4 (dm2)
    kin4 <- kin # transfer rate to intestine compartment 5 (/hr)
    
    Vin5 <- Vin # volume of intestine compartment 5 (L)
    SASIL5 <- SAin # surface area of intestine compartment 5 (dm2)
    kin5 <- kin # transfer rate to intestine compartment 6 (/hr)
    
    Vin6 <- Vin # volume of intestine compartment 6 (L)
    SASIL6 <- SAin # surface area of intestine compartment 6 (dm2)
    kin6 <- kin # transfer rate to intestine compartment 7 (hr)
    
    Vin7 <- Vin # volume of intestine compartment 7 (L)
    SASIL7 <- SAin # surface area of intestine compartment 7 (dm2)
    
    Vin8 <- Vin # volume of intestine compartment 7 (L)
    SASIL8 <- SAinb # surface area of intestine compartment 7 (dm2)
    
    
    ### Kinetic #####
    
    # intestinal absorption
    PappZEN <- 10^(0.6836*PappCaco2ZEN-0.5579)*3600/10 #(dm/h)
    kabin1 <- PappZEN*SASIL1 # absorption rate constant of intestine compartment 1 (L/hr)
    kabin2 <- PappZEN*SASIL2 # absorption rate constant of intestine compartment 2 (L/hr)
    kabin3 <- PappZEN*SASIL3 # absorption rate constant of intestine compartment 3 (L/hr)
    kabin4 <- PappZEN*SASIL4 # absorption rate constant of intestine compartment 4 (L/hr)
    kabin5 <- PappZEN*SASIL5 # absorption rate constant of intestine compartment 5 (L/hr)
    kabin6 <- PappZEN*SASIL6 # absorption rate constant of intestine compartment 6 (L/hr)
    kabin7 <- PappZEN*SASIL7 # absorption rate constant of intestine compartment 7 (L/hr)
    
    Kb <- PappZEN*SASIL8 # absorption rate constant from the large intestine (L/hr)
    
    PappaZEL <- 10^(0.6836*PappCaco2aZEL-0.5579)*3600/10 #(dm/h)
    kabin8 <- PappaZEL*SASIL8 # absorption rate constant of α-ZEL from LI to liver (L/hr) excretion (/h)
    
    # intestinal metabolism
    SI <- VSIc*1000  # small intestine in body weight (gram/kg BW)
    # metabolites ZENGLU scaled maximum rates of metabolism (μmol/h)
    VmaxSI_ZENGLU <- VmaxSIZENGLUc/1000*60* S9SI *SI*BW
    
    # metabolites α-ZEL and ß-ZEL scaled maximum rates of metabolism (μmol/h)
    VmaxLIL_aZEL <-  VmaxLIaZELc/1000000*1000 *VMB*BW
    VmaxLIL_bZEL <-  VmaxLIbZELc/1000000*1000 *VMB*BW
    
    # liver metabolism
    L <- VLc*1000  # (gram/kg BW) # liver
    # metabolites α-ZEL and ß-ZEL scaled maximum rates of metabolism (μmol/h)
    VmaxLaZEL <- VmaxLaZELc/1000000*60* VLS9 *L*BW
    VmaxLbZEL <-  VmaxLbZELc/1000000*60* VLS9 *L*BW
    
    # metabolites ZENGLU scaled maximum rates of metabolism (μmol/h)
    VmaxL_ZENGLU <- VmaxLZENGLUc /1000*60* VLS9 *L*BW
    
    # metabolites α-ZELGLU scaled maximum rates of metabolism (μmol/h)
    VmaxL_aZELGLU <- VmaxLaZELGLUc/1000*60* VLS9 *L*BW
    
    fup <- 0.156
    
    ## Concentrations #####
    
    CSIL_ZEN1 <- ASIL_ZEN1/Vin1
    CSIL_ZEN2 <- ASIL_ZEN2/Vin2
    CSIL_ZEN3 <- ASIL_ZEN3/Vin3
    CSIL_ZEN4 <- ASIL_ZEN4/Vin4
    CSIL_ZEN5 <- ASIL_ZEN5/Vin5
    CSIL_ZEN6 <- ASIL_ZEN6/Vin6
    CSIL_ZEN7 <- ASIL_ZEN7/Vin7
    
    CSI_ZEN <- ASI_ZEN/VSI
    CLIL_ZEN <- ALIL_ZEN/(VMB*BW) 
    CL_ZEN <- AL_ZEN/VL
    CF <- AF/VF
    CR <- AR/VR
    CS <- AS/VS
    CB <- AB/VB
    
    CVSI_ZEN <- CSI_ZEN/PIZEN
    CVL_ZEN <- CL_ZEN/PLZEN
    CVF_ZEN <- CF_ZEN/PFZEN
    CVR <- CR/PRZEN
    CVS <- CS/PSZEN
    
    
    # alpha Zearalenol
    
    CLIL_aZEL <- ALIL_aZEL/(VMB*BW)
    CLaZEL <- ALaZEL/VL
    CFaZEL <- AFaZEL/VF
    CRaZEL <- ARaZEL/VR
    CSaZEL <- ASaZEL/VS
    CBaZEL <- ABaZEL/VB
    
    CVLaZEL <- CLaZEL/PLaZEL
    CVFaZEL <- CFaZEL/PFaZEL
    CVRaZEL <- CRaZEL/PRaZEL
    CVSaZEL <- CSaZEL/PSaZEL
    
    # beta Zearalenol
    CLIL_bZEL <- ALIL_bZEL/VL
    CLbZEL <- ALbZEL/VL
    
    # Zearalenone glucuronide
    CSI_ZENGLU <- ASI_ZENGLU/VSI # small intestine zen gluc
    CL_ZENGLU <- AL_ZENGLU/VL    # liver zen gluc
    
    # Zearalenol glucuronide
    CL_aZELGLU <- AL_aZELGLU/(VMB*BW)  # liver a zel gluc
    CLIL_bZEL <- ALIL_bZEL/VL      # liver b zel gluc
    
    # Newly added
    CVL_ZENGLU <- CL_ZENGLU/PLaZEL    # liver zen gluc
    CB_ZENGLU <- AB_ZENGLU/VB
    
    
    ## Dynamic ####
    
    # Stomach:
    dAST_ZEN <- - ksto*AST_ZEN # amount in stomach ODOSEumol
    
    # small intestine lumen compartment
    # intestines, divided in 7 compartments
    dASIL_ZEN1 <- ksto*AST_ZEN - kin1*ASIL_ZEN1 - kabin1*CSIL_ZEN1 # Amount ZEN in intestine compartment 1 (μmol)
    dASIL_ZEN2 <- kin1*ASIL_ZEN1 - kin2*ASIL_ZEN2 - kabin2*CSIL_ZEN2 # Amount ZEN in intestine compartment 2 (μmol)
    dASIL_ZEN3 <- kin2*ASIL_ZEN2 - kin3*ASIL_ZEN3 - kabin3*CSIL_ZEN3 # Amount ZEN in intestine compartment 3 (μmol)
    dASIL_ZEN4 <- kin3*ASIL_ZEN3 - kin4*ASIL_ZEN4 - kabin4*CSIL_ZEN4 # Amount ZEN in intestine compartment 4 (μmol)
    dASIL_ZEN5 <- kin4*ASIL_ZEN4 - kin5*ASIL_ZEN5 - kabin5*CSIL_ZEN5 # Amount ZEN in intestine compartment 5 (μmol)
    dASIL_ZEN6 <- kin5*ASIL_ZEN5 - kin6*ASIL_ZEN6 - kabin6*CSIL_ZEN6 # Amount ZEN in intestine compartment 6 (μmol)
    dASIL_ZEN7 <- kin6*ASIL_ZEN6 - kin7*ASIL_ZEN7 - kabin7*CSIL_ZEN7 # Amount ZEN in intestine compartment 7 (μmol)
    
    # small intestine tissue compartment
    dASI_ZEN <- kabin1*CSIL_ZEN1 + kabin2*CSIL_ZEN2 + kabin3*CSIL_ZEN3 +
      + kabin4*CSIL_ZEN4 + kabin5*CSIL_ZEN5 + kabin6*CSIL_ZEN6 + kabin7*CSIL_ZEN7 +
      + QSI*(CB-CVSI_ZEN) +
      - VmaxSI_ZENGLU*CVSI_ZEN*fup / (KmSIZENGLU + CVSI_ZEN*fup)
      # + (ASI_ZENGLU) # amount of ZEN in small intestinal tissue, (μmol)
    
    
    dASI_ZENGLU <- VmaxSI_ZENGLU*CVSI_ZEN*fup / (KmSIZENGLU + CVSI_ZEN*fup) #+ AL_ZENGLU - (ASI_ZENGLU) # amount of ZEN metabolized to metabolite ZENGLU, (μmol)
    
    # large intestine lumen compartment: microbial activity
    dALIL_ZEN <- kin7*ASIL_ZEN7 + 
      - VmaxLIL_aZEL*CLIL_ZEN / (KmLIaZEL + CLIL_ZEN) +
      - VmaxLIL_bZEL*CLIL_ZEN / (KmLIbZEL + CLIL_ZEN) +
      - Kb*CLIL_ZEN - KfeZ*ALIL_ZEN #+
      + 0.3*AL_ZENGLU # amount of ZEN in large intestine lumen, (μmol)
    
    dALIL_aZEL <- VmaxLIL_aZEL*CLIL_ZEN / (KmLIaZEL + CLIL_ZEN) +
      - kabin8*(ALIL_aZEL/VMB*BW) - KfeaZ*ALIL_aZEL # amount of α-ZEL formed due to gut microbial activity, (μmol)
    
    dALIL_bZEL <- VmaxLIL_bZEL*CLIL_ZEN / (KmLIbZEL + CLIL_ZEN) # amount of ß-ZEL formed due to gut microbial activity, (μmol)
    
    dAFE_ZEN <- KfeZ*ALIL_ZEN
    dAFE_aZEL <- KfeaZ*ALIL_aZEL
    
    #  liver compartment
    dAL_ZEN <- QSI*CVSI_ZEN + QL*CB - (QSI+QL)*CVL_ZEN +
      - VmaxLaZEL*CVL_ZEN*fup / (KmLaZEL + CVL_ZEN*fup) +
      - VmaxLbZEL*CVL_ZEN*fup / (KmLbZEL + CVL_ZEN*fup) +
      - VmaxL_ZENGLU*CVL_ZEN*fup / (KmLZENGLU + CVL_ZEN*fup) #+
      # + Kb*CLIL_ZEN # amount of ZEN in liver, (μmol)
    
    dALaZEL <- VmaxLaZEL*CVL_ZEN*fup / (KmLaZEL + CVL_ZEN*fup) +
      + kabin8*(ALIL_aZEL/VMB*BW) +
      - VmaxL_aZELGLU*CVLaZEL*fup / (KmLaZELGLU + CVLaZEL*fup) +
      - (QL+QSI)*CVLaZEL + (QL+QSI)*CBaZEL # amount of ZEN metabolized to metabolite α-ZEL in liver, (μmol)
    
    dALbZEL <- VmaxLbZEL*CVL_ZEN*fup / (KmLbZEL + CVL_ZEN*fup) # amount of ZEN metabolized to metabolite ß-ZEL in liver, (μmol)
   
    dAL_ZENGLU <- VmaxL_ZENGLU*CVL_ZEN*fup / (KmLZENGLU + CVL_ZEN*fup) +
      - 0.3*AL_ZENGLU +
      - (QSI+QL)*CVL_ZENGLU + (QSI+QL)*CB_ZENGLU # amount of ZEN metabolized to metabolite ZENGLU in liver, (μmol)
    
    dAB_ZENGLU <- (QSI+QL)*CVL_ZENGLU - (QSI+QL)*CB_ZENGLU - (0.1*CB_ZENGLU*fup)
    dAUR_ZENGLU <- (0.1*CB_ZENGLU*fup)
    
    dAL_aZELGLU <- VmaxL_aZELGLU*CVLaZEL*fup / (KmLaZELGLU + CVLaZEL*fup) # amount of α-ZEL glucuronide in liver
    
    #  fat compartment
    dAF_ZEN <- QF*(CB-CVF_ZEN) # amount of ZEN in fat tissue, (μmol)
    
    # rapidly perfused tissue
    dAR <- QR*(CB-CVR) # amount of ZEN in rapidly perfused tissue, (μmol)
    
    # slowly perfused tissue
    dAS <- QS*(CB-CVS) # amount of ZEN in slowly perfused tissue, (μmol)
    
    #  blood compartment
    dAB <- - QSI*CB - QL*CB + (QSI+QL)*CVL_ZEN - QF*CB + QF*CVF_ZEN +
      - QR*CB + QR*CVR - QS*CB + QS*CVS - KurZ*AB*fup + Kb*CLIL_ZEN# amount of ZEN in blood, (μmol) IVDOSEumol
    
    # urinary excretion
    dAZur <- KurZ*AB*fup
    
    # fat compartment
    dAFaZEL <- QF*(CBaZEL-CVFaZEL) # amount of α-ZEL in fat tissue, (μmol)
    
    # rapidly perfused tissue
    dARaZEL <- QR*(CBaZEL-CVRaZEL) # amount of α-ZEL in rapidly perfused tissue, (μmol)
    
    #slowly perfused tissue
    dASaZEL <- QS*(CBaZEL-CVSaZEL) # amount of α-ZEL in slowly perfused tissue, (μmol)
    
    #blood compartment
    dABaZEL <- + (QL+QSI)*CVLaZEL - (QL+QSI)*CBaZEL +
      - QF*(CBaZEL-CVFaZEL) - QR*(CBaZEL-CVRaZEL) - QS*(CBaZEL-CVSaZEL) +
      - KurZ*ABaZEL*fup # amount of α-ZELin blood, (μmol)
    
    # urinary excretion
    dAaZur <- KurZ*ABaZEL*fup
    
    #mass balance calculation
    dATotal <- AST_ZEN #ODOSEumol + IVDOSEumol
    ACalculated <- AST_ZEN+ ASIL_ZEN1 + ASIL_ZEN2 + ASIL_ZEN3 + ASIL_ZEN4 + ASIL_ZEN5 + ASIL_ZEN6 + ASIL_ZEN7 +
      + ASI_ZEN + ASI_ZENGLU +
      + ALIL_ZEN + ALIL_aZEL + ALIL_bZEL + AFE_ZEN +
      + AL_ZEN + ALaZEL + ALbZEL + 
      + AL_ZENGLU + AB_ZENGLU + AUR_ZENGLU +
      + AL_aZELGLU +
      + AF_ZEN + AR + AS +
      + AB + AZur +
      + AFaZEL + ARaZEL + ASaZEL +
      + ABaZEL +
      + AaZur 
    MB <- ATotal - ACalculated + 1
    
    # End
    # 
    list(c(
      dAST_ZEN,
      dASIL_ZEN1,
      dASIL_ZEN2,
      dASIL_ZEN3,
      dASIL_ZEN4,
      dASIL_ZEN5,
      dASIL_ZEN6,
      dASIL_ZEN7,
      dASI_ZEN,
      dASI_ZENGLU,
      dALIL_ZEN,
      dALIL_aZEL,
      dALIL_bZEL,
      dAFE_ZEN,
      dAFE_aZEL,
      dAL_ZEN,
      dALaZEL,
      dALbZEL,
      dAL_ZENGLU,
      dAB_ZENGLU, #new
      dAUR_ZENGLU, #new
      dAL_aZELGLU,
      dAF_ZEN,
      dAR,
      dAS,
      dAB,
      dAZur,
      dAFaZEL,
      dARaZEL,
      dASaZEL,
      dABaZEL,
      dAaZur,
      dATotal
    ),
    c(ACalculated = ACalculated,
      MB = MB,
      CSIL_ZEN1 = CSIL_ZEN1,
      CSIL_ZEN2 = CSIL_ZEN2,
      CSIL_ZEN3 = CSIL_ZEN3,
      CSIL_ZEN4 = CSIL_ZEN4,
      CSIL_ZEN5 = CSIL_ZEN5,
      CSIL_ZEN6 = CSIL_ZEN6,
      CSIL_ZEN7 = CSIL_ZEN7,
      CSI_ZEN = CSI_ZEN, CVSI_ZEN = CVSI_ZEN,
      CSI_ZENGLU = CSI_ZENGLU,
      CLIL_ZEN = CLIL_ZEN, 
      CLIL_aZEL = CLIL_aZEL, 
      CLIL_bZEL = CLIL_bZEL,
      CL_ZEN = CL_ZEN, CVL_ZEN = CVL_ZEN,
      CLaZEL = CLaZEL, CVLaZEL = CVLaZEL,
      CLbZEL = CLbZEL,
      CL_ZENGLU = CL_ZENGLU,
      CB_ZENGLU = CB_ZENGLU,
      CL_aZELGLU = CL_aZELGLU,
      CF = CF, CVF = CVF,
      CR = CR, CVR = CVR,
      CS = CS, CVS = CVS,
      CB = CB,
      CFaZEL = CFaZEL, CVFaZEL = CVFaZEL,
      CRaZEL = CRaZEL, CVRaZEL = CVRaZEL,
      CSaZEL = CSaZEL, CVSaZEL = CVSaZEL,
      CBaZEL = CBaZEL
    ))

  })
}