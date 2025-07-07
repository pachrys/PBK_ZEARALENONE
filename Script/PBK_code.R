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
    
    
    # Intestinal lumen 
    # Small intestine
    Vin <- 0.01424*BW # volume of intestine compartment 1 (L)
    SAin <- SAin*7 # surface area of intestine compartment 1 (dm2)
    kin <- kin # transfer rate to intestine compartment 2 (/hr)
    
    Vin8 <- Vin # volume of intestine compartment 7 (L)
    SAin8 <- SAinb # surface area of intestine compartment 7 (dm2)
    
    # Large intestine
    Vco <- 0.012*BW
    
    ### Kinetic #####
    
    # intestinal absorption
    PappZEN <- 10^(0.6836*PappCaco2ZEN-0.5579)*3600/10 #(dm/h)
    kabin <- PappZEN*SAin # absorption rate constant (L/hr)
    
    Kb <- PappZEN*SAinb # absorption rate constant from the large intestine (L/hr)
    
    PappaZEL <- 10^(0.6836*PappCaco2aZEL-0.5579)*3600/10 #(dm/h)
    kabin8 <- PappaZEL*SAin8 # absorption rate constant of α-ZEL from LI to liver (L/hr) excretion (/h)
    
    # intestinal metabolism
    SI <- VSIc*1000  # small intestine in body weight (gram/kg BW)
    # metabolites ZENGLU scaled maximum rates of metabolism (μmol/h)
    VmaxSIZENGLU <- VmaxSIZENGLUc/1000*60* S9SI *SI*BW
    
    # metabolites α-ZEL and ß-ZEL scaled maximum rates of metabolism (μmol/h)
    VmaxLIaZEL <-  VmaxLIaZELc/1000000*1000 *VMB*BW
    VmaxLIbZEL <-  VmaxLIbZELc/1000000*1000 *VMB*BW
    
    # liver metabolism
    L <- VLc*1000  # (gram/kg BW) # liver
    # metabolites α-ZEL and ß-ZEL scaled maximum rates of metabolism (μmol/h)
    VmaxLaZEL <- VmaxLaZELc/1000000*60* VLS9 *L*BW
    VmaxLbZEL <-  VmaxLbZELc/1000000*60* VLS9 *L*BW
    
    # metabolites ZENGLU scaled maximum rates of metabolism (μmol/h)
    VmaxLZENGLU <- VmaxLZENGLUc /1000*60* VLS9 *L*BW
    
    # metabolites α-ZELGLU scaled maximum rates of metabolism (μmol/h)
    VmaxLaZELGLU <- VmaxLaZELGLUc/1000*60* VLS9 *L*BW
    
    # my additions
    fup <- 0.156
    
    ## Concentrations #####
    
    Cin <- Ain/Vin
   
    CSIZEN <- ASIZEN/VSI
    CLIZEN <- ALIZEN/(VMB*BW) 
    CLZEN <- ALZEN/VL
    CF <- AF/VF
    CR <- AR/VR
    CS <- AS/VS
    CB <- AB/VB
    
    CVSIZEN <- CSIZEN/PIZEN
    CVLZEN <- CLZEN/PLZEN
    CVF <- CF/PFZEN
    CVR <- CR/PRZEN
    CVS <- CS/PSZEN
    
    
    # alpha Zearalenol
    
    CLIaZEL <- ALIaZEL/(VMB*BW)
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
    CLIbZEL <- ALIbZEL/VL
    CLbZEL <- ALbZEL/VL
    
    # Zearalenone glucuronide
    CSIZENGLU <- ASIZENGLU/VSI # small intestine zen gluc
    CLZENGLU <- ALZENGLU/VL    # liver zen gluc
    
    # Zearalenol glucuronide
    CLaZELGLU <- ALaZELGLU/(VMB*BW)  # liver a zel gluc
    CLIbZEL <- ALIbZEL/VL      # liver b zel gluc
    
    # Newly added
    CVLZENGLU <- CLZENGLU/PLaZEL    # liver zen gluc
    CBZENGLU <- ABZENGLU/VB
    
    
    ## Dynamic ####
    
    # Stomach:
    dAst <- - ksto*Ast # amount in stomach ODOSEumol
    
    # small intestine lumen compartment
    # intestines, divided in 7 compartments
    dAin <- ksto*Ast - kin*Ain - kabin*Cin # Amount ZEN in intestine compartment 1 (μmol)
    
    # small intestine tissue compartment
    dASIZEN <- kabin*Cin + 
      + QSI*(CB-CVSIZEN) +
      - VmaxSIZENGLU*CVSIZEN*fup / (KmSIZENGLU + CVSIZEN*fup)
     # + (ASIZENGLU) # amount of ZEN in small intestinal tissue, (μmol)
    
    
    dASIZENGLU <- VmaxSIZENGLU*CVSIZEN*fup / (KmSIZENGLU + CVSIZEN*fup) #+ ALZENGLU - (ASIZENGLU) # amount of ZEN metabolized to metabolite ZENGLU, (μmol)
    
    # large intestine lumen compartment: microbial activity
    dALIZEN <- kin*Ain + 
      - VmaxLIaZEL*CLIZEN / (KmLIaZEL + CLIZEN) +
      - VmaxLIbZEL*CLIZEN / (KmLIbZEL + CLIZEN) +
      - Kb*CLIZEN - KfeZ*ALIZEN +
      + ALZENGLU # amount of ZEN in large intestine lumen, (μmol)
    
    dALIaZEL <- VmaxLIaZEL*CLIZEN / (KmLIaZEL + CLIZEN) +
      - kabin8*(ALIaZEL/VMB*BW) - KfeaZ*ALIaZEL # amount of α-ZEL formed due to gut microbial activity, (μmol)
    
    dALIbZEL <- VmaxLIbZEL*CLIZEN / (KmLIbZEL + CLIZEN) # amount of ß-ZEL formed due to gut microbial activity, (μmol)
    
    dAfeZ <- KfeZ*ALIZEN
    dAfeaZ <- KfeaZ*ALIaZEL
    
    #  liver compartment
    dALZEN <- QSI*CVSIZEN + QL*CB - (QSI+QL)*CVLZEN +
      - VmaxLaZEL*CVLZEN*fup / (KmLaZEL + CVLZEN*fup) +
      - VmaxLbZEL*CVLZEN*fup / (KmLbZEL + CVLZEN*fup) +
      - VmaxLZENGLU*CVLZEN*fup / (KmLZENGLU + CVLZEN*fup) #+
      # + Kb*CLIZEN # amount of ZEN in liver, (μmol)
    
    dALaZEL <- VmaxLaZEL*CVLZEN*fup / (KmLaZEL + CVLZEN*fup) +
      + kabin8*(ALIaZEL/VMB*BW) +
      - VmaxLaZELGLU*CVLaZEL*fup / (KmLaZELGLU + CVLaZEL*fup) +
      - (QL+QSI)*CVLaZEL + (QL+QSI)*CBaZEL # amount of ZEN metabolized to metabolite α-ZEL in liver, (μmol)
    
    dALbZEL <- VmaxLbZEL*CVLZEN*fup / (KmLbZEL + CVLZEN*fup) # amount of ZEN metabolized to metabolite ß-ZEL in liver, (μmol)
   
    dALZENGLU <- VmaxLZENGLU*CVLZEN*fup / (KmLZENGLU + CVLZEN*fup) +
      - ALZENGLU +
      - (QSI+QL)*CVLZENGLU + (QSI+QL)*CBZENGLU # amount of ZEN metabolized to metabolite ZENGLU in liver, (μmol)
    
    dABZENGLU <- (QSI+QL)*CVLZENGLU - (QSI+QL)*CBZENGLU - (0.1*CBZENGLU*fup)
    dAZENGLUur <- (0.1*CBZENGLU*fup)
    
    dALaZELGLU <- VmaxLaZELGLU*CVLaZEL*fup / (KmLaZELGLU + CVLaZEL*fup) # amount of α-ZEL glucuronide in liver
    
    #  fat compartment
    dAF <- QF*(CB-CVF) # amount of ZEN in fat tissue, (μmol)
    
    # rapidly perfused tissue
    dAR <- QR*(CB-CVR) # amount of ZEN in rapidly perfused tissue, (μmol)
    
    # slowly perfused tissue
    dAS <- QS*(CB-CVS) # amount of ZEN in slowly perfused tissue, (μmol)
    
    #  blood compartment
    dAB <- - QSI*CB - QL*CB + (QSI+QL)*CVLZEN - QF*CB + QF*CVF +
      - QR*CB + QR*CVR - QS*CB + QS*CVS - KurZ*AB*fup + Kb*CLIZEN# amount of ZEN in blood, (μmol) IVDOSEumol
    
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
    dATotal <- Ast #ODOSEumol + IVDOSEumol
    ACalculated <- Ast + Ain +
      + ASIZEN + ASIZENGLU +
      + ALIZEN + ALIaZEL + ALIbZEL + AfeZ +
      + ALZEN + ALaZEL + ALbZEL + 
      + ALZENGLU + ABZENGLU + AZENGLUur +
      + ALaZELGLU +
      + AF + AR + AS +
      + AB + AZur +
      + AFaZEL + ARaZEL + ASaZEL +
      + ABaZEL +
      + AaZur 
    MB <- ATotal - ACalculated + 1
    
    # End
    # 
    list(c(
      dAst,
      dAin,
      dASIZEN,
      dASIZENGLU,
      dALIZEN,
      dALIaZEL,
      dALIbZEL,
      dAfeZ,
      dAfeaZ,
      dALZEN,
      dALaZEL,
      dALbZEL,
      dALZENGLU,
      dABZENGLU, #new
      dAZENGLUur, #new
      dALaZELGLU,
      dAF,
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
      Cin = Cin,
      CSIZEN = CSIZEN, CVSIZEN = CVSIZEN,
      CSIZENGLU = CSIZENGLU,
      CLIZEN = CLIZEN, 
      CLIaZEL = CLIaZEL, 
      CLIbZEL = CLIbZEL,
      CLZEN = CLZEN, CVLZEN = CVLZEN,
      CLaZEL = CLaZEL, CVLaZEL = CVLaZEL,
      CLbZEL = CLbZEL,
      CLZENGLU = CLZENGLU,
      CBZENGLU = CBZENGLU,
      CLaZELGLU = CLaZELGLU,
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