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
    SAin1 <- SAin # surface area of intestine compartment 1 (dm2)
    kin1 <- kin # transfer rate to intestine compartment 2 (/hr)
    
    Vin2 <- Vin # volume of intestine compartment 2 (L)
    SAin2 <- SAin # surface area of intestine compartment 2 (dm2)
    kin2 <- kin # transfer rate to intestine compartment 3 (/hr)
    
    Vin3 <- Vin # volume of intestine compartment 3 (L)
    SAin3 <-  SAin # surface area of intestine compartment 3 (dm2)
    kin3 <- kin # transfer rate to intestine compartment 4 (/hr)
    
    Vin4 <- Vin # volume of intestine compartment 4 (L)
    SAin4 <- SAin # surface area of intestine compartment 4 (dm2)
    kin4 <- kin # transfer rate to intestine compartment 5 (/hr)
    
    Vin5 <- Vin # volume of intestine compartment 5 (L)
    SAin5 <- SAin # surface area of intestine compartment 5 (dm2)
    kin5 <- kin # transfer rate to intestine compartment 6 (/hr)
    
    Vin6 <- Vin # volume of intestine compartment 6 (L)
    SAin6 <- SAin # surface area of intestine compartment 6 (dm2)
    kin6 <- kin # transfer rate to intestine compartment 7 (hr)
    
    Vin7 <- Vin # volume of intestine compartment 7 (L)
    SAin7 <- SAin # surface area of intestine compartment 7 (dm2)
    
    Vin8 <- Vin # volume of intestine compartment 7 (L)
    SAin8 <- SAinb # surface area of intestine compartment 7 (dm2)
    
    
    ### Kinetic #####
    
    # intestinal absorption
    PappZEN <- 10^(0.6836*PappCaco2ZEN-0.5579)*3600/10 #(dm/h)
    kabin1 <- PappZEN*SAin # absorption rate constant (L/hr)
    kabin1 <- PappZEN*SAin1 # absorption rate constant of intestine compartment 1 (L/hr)
    kabin2 <- PappZEN*SAin2 # absorption rate constant of intestine compartment 2 (L/hr)
    kabin3 <- PappZEN*SAin3 # absorption rate constant of intestine compartment 3 (L/hr)
    kabin4 <- PappZEN*SAin4 # absorption rate constant of intestine compartment 4 (L/hr)
    kabin5 <- PappZEN*SAin5 # absorption rate constant of intestine compartment 5 (L/hr)
    kabin6 <- PappZEN*SAin6 # absorption rate constant of intestine compartment 6 (L/hr)
    kabin7 <- PappZEN*SAin7 # absorption rate constant of intestine compartment 7 (L/hr)
    
    Kb <- PappZEN*SAinb # absorption rate constant from the large intestine (L/hr)
    
    PappaZEL <- 10^(0.6836*PappCaco2aZEL-0.5579)*3600/10 #(dm/h)
    kabin8 <- PappaZEL*SAin8 # absorption rate constant of α-ZEL from LI to liver (L/hr) excretion (/h)
    
    # intestinal metabolism
    SI <- VSIc*1000  # small intestine in body weight (gram/kg BW)
    # metabolites ZENGLU scaled maximum rates of metabolism (μmol/h)
    VmaxSIZENGLU <- VmaxSIZENGLUc/1000*60*S9SI*SI*BW
    
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
    
    
    ## Concentrations #####
    
    Cin1 <- Ain1/Vin1
    Cin2 <- Ain2/Vin2
    Cin3 <- Ain3/Vin3
    Cin4 <- Ain4/Vin4
    Cin5 <- Ain5/Vin5
    Cin6 <- Ain6/Vin6
    Cin7 <- Ain7/Vin7
    
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
    
    ## Dynamic ####
    
    # Stomach:
    dAst <- - ksto*Ast # amount in stomach ODOSEumol
    
    # small intestine lumen compartment
    # intestines, divided in 7 compartments
    dAin1 <- ksto*Ast - kin1*Ain1 - kabin1*Cin1 # Amount ZEN in intestine compartment 1 (μmol)
    dAin2 <- kin1*Ain1 - kin2*Ain2 - kabin2*Cin2 # Amount ZEN in intestine compartment 2 (μmol)
    dAin3 <- kin2*Ain2 - kin3*Ain3 - kabin3*Cin3 # Amount ZEN in intestine compartment 3 (μmol)
    dAin4 <- kin3*Ain3 - kin4*Ain4 - kabin4*Cin4 # Amount ZEN in intestine compartment 4 (μmol)
    dAin5 <- kin4*Ain4 - kin5*Ain5 - kabin5*Cin5 # Amount ZEN in intestine compartment 5 (μmol)
    dAin6 <- kin5*Ain5 - kin6*Ain6 - kabin6*Cin6 # Amount ZEN in intestine compartment 6 (μmol)
    dAin7 <- kin6*Ain6 - kin7*Ain7 - kabin7*Cin7 # Amount ZEN in intestine compartment 7 (μmol)
    
    # small intestine tissue compartment
    dASIZEN <- kabin1*Cin1 + kabin2*Cin2 + kabin3*Cin3 +
      + kabin4*Cin4 + kabin5*Cin5 + kabin6*Cin6 + kabin7*Cin7 +
      + QSI*(CB-CVSIZEN) +
      - VmaxSIZENGLU*CVSIZEN / (KmSIZENGLU+ CVSIZEN) # amount of ZEN in small intestinal tissue, (μmol)
    
    dASIZENGLU <- VmaxSIZENGLU*CVSIZEN / (KmSIZENGLU+ CVSIZEN) # amount of ZEN metabolized to metabolite ZENGLU, (μmol)
    
    # large intestine lumen compartment: microbial activity
    dALIZEN <- kin7*Ain7 + 
      - VmaxLIaZEL*CLIZEN / (KmLIaZEL + CLIZEN) +
      - VmaxLIbZEL*CLIZEN / (KmLIbZEL + CLIZEN) +
      - Kb*CLIZEN - KfeZ*ALIZEN # amount of ZEN in large intestine lumen, (μmol)
    
    dALIaZEL <- VmaxLIaZEL*CLIZEN / (KmLIaZEL + CLIZEN) +
      - kabin8*(ALIaZEL/VMB*BW) - KfeaZ*ALIaZEL # amount of α-ZEL formed due to gut microbial activity, (μmol)
    
    dALIbZEL <- VmaxLIbZEL*CLIZEN / (KmLIbZEL + CLIZEN) # amount of ß-ZEL formed due to gut microbial activity, (μmol)
    
    dAfeZ <- KfeZ*ALIZEN
    dAfeaZ <- KfeaZ*ALIaZEL
    
    #  liver compartment
    dALZEN <- QSI*CVSIZEN + QL*CB - (QSI+QL)*CVLZEN +
      - VmaxLaZEL*CVLZEN / (KmLaZEL + CVLZEN) +
      - VmaxLbZEL*CVLZEN / (KmLbZEL + CVLZEN) +
      - VmaxLZENGLU*CVLZEN / (KmLZENGLU + CVLZEN) +
      + Kb*CLIZEN # amount of ZEN in liver, (μmol)
    
    dALaZEL <- VmaxLaZEL*CVLZEN / (KmLaZEL + CVLZEN) +
      + kabin8*(ALIaZEL/VMB*BW) +
      - VmaxLaZELGLU*CVLaZEL / (KmLaZELGLU + CVLaZEL) +
      - (QL+QSI)*CVLaZEL + (QL+QSI)*CBaZEL # amount of ZEN metabolized to metabolite α-ZEL in liver, (μmol)
    
    dALbZEL <- VmaxLbZEL*CVLZEN / (KmLbZEL + CVLZEN) # amount of ZEN metabolized to metabolite ß-ZEL in liver, (μmol)
    dALZENGLU <- VmaxLZENGLU*CVLZEN / (KmLZENGLU + CVLZEN) # amount of ZEN metabolized to metabolite ZENGLU in liver, (μmol)
    
    dALaZELGLU <- VmaxLaZELGLU*CVLaZEL / (KmLaZELGLU + CVLaZEL) # amount of α-ZEL glucuronide in liver
    
    #  fat compartment
    dAF <- QF*(CB-CVF) # amount of ZEN in fat tissue, (μmol)
    
    # rapidly perfused tissue
    dAR <- QR*(CB-CVR) # amount of ZEN in rapidly perfused tissue, (μmol)
    
    # slowly perfused tissue
    dAS <- QS*(CB-CVS) # amount of ZEN in slowly perfused tissue, (μmol)
    
    #  blood compartment
    dAB <- - QSI*CB - QL*CB + (QSI+QL)*CVLZEN - QF*CB + QF*CVF +
      - QR*CB + QR*CVR - QS*CB + QS*CVS - KurZ*AB # amount of ZEN in blood, (μmol) IVDOSEumol
    
    # urinary excretion
    dAZur <- KurZ*AB
    
    # fat compartment
    dAFaZEL <- QF*(CBaZEL-CVFaZEL) # amount of α-ZEL in fat tissue, (μmol)
    
    # rapidly perfused tissue
    dARaZEL <- QR*(CBaZEL-CVRaZEL) # amount of α-ZEL in rapidly perfused tissue, (μmol)
    
    #slowly perfused tissue
    dASaZEL <- QS*(CBaZEL-CVSaZEL) # amount of α-ZEL in slowly perfused tissue, (μmol)
    
    #blood compartment
    dABaZEL <- + (QL+QSI)*CVLaZEL - (QL+QSI)*CBaZEL +
      - QF*(CBaZEL-CVFaZEL) - QR*(CBaZEL-CVRaZEL) - QS*(CBaZEL-CVSaZEL) +
      - KurZ*ABaZEL # amount of α-ZELin blood, (μmol)
    
    # urinary excretion
    dAaZur <- KurZ*ABaZEL
    
    #mass balance calculation
    dATotal <- ODOSEumol + IVDOSEumol
    ACalculated <- Ast+ Ain1 + Ain2 + Ain3 + Ain4 + Ain5 + Ain6 + Ain7 +
      + ASIZEN + ASIZENGLU +
      + ALIZEN + ALIaZEL + ALIbZEL + AfeZ +
      + ALZEN + ALaZEL + ALbZEL + ALZENGLU + ALaZELGLU +
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
      dAin1,
      dAin2,
      dAin3,
      dAin4,
      dAin5,
      dAin6,
      dAin7,
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
      Cin1 = Cin1,
      Cin2 = Cin2,
      Cin3 = Cin3,
      Cin4 = Cin4,
      Cin5 = Cin5,
      Cin6 = Cin6,
      Cin7 = Cin7,
      CSIZEN = CSIZEN, CVSIZEN = CVSIZEN,
      CSIZENGLU = CSIZENGLU,
      CLIZEN = CLIZEN, 
      CLIaZEL = CLIaZEL, 
      CLIbZEL = CLIbZEL,
      CLZEN = CLZEN, CVLZEN = CVLZEN,
      CLaZEL = CLaZEL, CVLaZEL = CVLaZEL,
      CLbZEL = CLbZEL,
      CLZENGLU = CLZENGLU,
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