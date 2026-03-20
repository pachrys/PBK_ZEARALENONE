# --------------------------------------------------------------------------- #
# Script for the Initial PBK Model 
# Based on the SI of Mendez-Catala et al. 2021 DOI: 10.1002/mnfr.202100443
# By: Chrysanthi Pachoulide
# Date: 23-06-2025
# --------------------------------------------------------------------------- #

rm(list=ls()) 

# Packages
library(here)
library(tidyverse)
library(deSolve)
library(showtext)
library(plotly)
library(readr)
font_add(family = "Garamond", regular = "GARA.TTF")
showtext_auto()

# For plotting

CP_theme <- theme_minimal() +
  theme(
    plot.background = element_rect(fill = 'white'), #element_blank(), #
    plot.margin = margin(0.3, 0.3, 0.3, 0.3, "cm"), 
    plot.title = element_text(hjust = 0, size = 48, face = "bold"),
    plot.title.position = "plot", 
    panel.grid = element_line(linewidth = 0.1),  
    axis.title = element_text(size = 44, face = "bold"),
    axis.text = element_text(size = 44),
    axis.line = element_line(linewidth = 0.5),
    legend.text = element_text(size = 44),
    legend.title = element_text(size = 44), 
    legend.key.spacing = unit(0.1, "cm"),
    legend.key.size = unit(0.5, "lines"),
    legend.spacing.y = unit(0.2, "cm"), 
    legend.spacing.x = unit(0.2, "cm"),
    legend.margin = margin(0,0,0,0),
    legend.position = "right"
  )


# Set output storage directory

OUTPUT <- here("Output", format(Sys.Date(), "%Y-%m-%d")) 
dir.create(OUTPUT, recursive = TRUE)


# -------------------------------------------------------------------------#
# -------------------------------------------------------------------------#





# Calculate PBK model parameters ####
# -------------------------------------------------------------------------#

## Physiological constants ####
# Tissue volumes (reference: Brown et al. (1997))
BW = 0.250  ##{Kg} #  body weight rat
VSIc = 0.014  # fraction of small intestine
VLc = 0.034  # fraction of liver tissue
VRc = 0.034  # 0.082-VSIc-VLc # fraction of rapidly perfused tissue
VSc = 0.667  # 0.737-VFc # fraction of slowly perfused tissue
VFc = 0.070  # fraction of fat tissue
VBc = 0.074  # fraction of blood

# Blood flow rates (reference: Brown et al. (1997))
QC = 5.38  ##{L/h} # cardiac output:15*BW^0.74
QSIc = 0.151  # fraction of blood flow to small intestine
QLc = 0.099  # 0.25 - QSIc # fraction of blood flow to liver
QRc = 0.51  # 0.76 - QSIc - QLc # fraction of blood flow to rapidly perfused tissue
QSc = 0.17  # 0.24 - QFc # fraction of blood flow to slowly perfused tissue
QFc = 0.07  # fraction of blood flow to fat

## Physicochemical parameters ####

### partition coefficients ####
# calculated based on QPPR of DeJongh et al. (1997)
# ZEN in main model
PIZEN = 2.64  # intestine/blood partition coefficient
PLZEN = 2.64  # liver/blood partition coefficient
PRZEN = 2.64  # rapidly perfused tissue/blood partition coefficient
PSZEN = 0.76  # slowly perfused tissue/blood partition coefficient
PFZEN = 106.93  # fat/blood partition coefficient
# α-ZEL in sub-model
PIaZEL = 2.38  # intestine/blood partition coefficient
PLaZEL = 2.38  # liver/blood partition coefficient
PRaZEL = 2.38  # rapidly perfused tissue/blood partition coefficient
PSaZEL = 0.71  # slowly perfused tissue/blood partition coefficient
PFaZEL = 92.79  # fat/blood partition coefficient

### absorption/transfer rates ####
# Stomach emptying rate
ksto = 2.8  # stomach emptying rate (/h) reference: Reilly et al. (1990)

# intestinal absorption and transfer rates of ZEN
PappCaco2ZEN = 1.01 #3.9838  # fitted Log Papp  Caco-2 (cm/h) reference: Pfeiffer et al. (2011)
# Log (Papp in vivo) = 0.6836*Log(PappCaco2ZEN)-0.5579 # reference: Sun et al. (2002)

Vin = 0.0012  # volume for each compartment of intestines (L)
SAin = 0.134  # surface area (dm2)
kin = 4.17  # transfer rate to next compartment within the intestines (/hr)
kin7 = 0.464  #(/h) # transfer rate of ZEN from small intestine to large intestine reference:  Kimura and Higaki (2002)
Vinb = 0.0338  # volume of large intestine (L) reference: Vdoviaková et al. (2016)
SAinb = 1.57  # surface area of large intestine (dm2) reference: Vdoviaková et al. (2016)

# α-ZEL
PappCaco2aZEL = 5.2676  #Log Papp  Caco-2 reference: Pfeiffer et al. (2011) (cm/h)
# Log (Papp in vivo) = 0.6836*Log(PappCaco2aZEL)-0.5579 reference: Sun et al. (2002)

# urinary excretion
KurZ = 0.096  # urinary excretion of ZEN reference: Mukherjee et al. (2014)
KfeZ = 0.69  # fecal excretion of ZEN reference: Mukherjee et al. (2014)
KuraZ = 0.012  # urinary excretion of α-ZEL reference: Mukherjee et al. (2014)
KfeaZ = 0.69  # fecal excretion of α-ZEL reference: Mukherjee et al. (2014)


## Kinetic parameters ####

### metabolism in small intestine tissue ####
# scaling factors
S9SI =  37.1  # small intestinal S9 protein yield (mg S9 protein/gram intestine) reference: Peters et al. (2016)
# metabolites ZENGLU  unscaled maximum rates of metabolism  (nmol/min/mg S9 protein)
VmaxSIZENGLUc =  8.617  # in vitro  S9 incubations from this study

# metabolites ZENGLU  affinity constants  (μmol/L)
KmSIZENGLU = 6.198  # in vitro  S9 incubations from this study

### metabolism in large intestine lumen (microbiota compartment) ####
# scaling factors
VMB = 0.0164  # fraction of feces of BW reference: Hoskins and Zamcheck (1968)
# metabolites α-ZEL and ß-ZEL  unscaled maximum rates of metabolism  (pmol/h/g faeces)
VmaxLIaZELc = 0.0138  # reference: Mendez-Catala et al. (2020)
VmaxLIbZELc =  0.0084  # reference: Mendez-Catala et al. (2020)

# metabolites α-ZEL and ß-ZEL  affinity constants  (μmol/L)
KmLIaZEL = 66  # reference: Mendez-Catala et al. (2020)
KmLIbZEL = 80  # reference: Mendez-Catala et al. (2020)

### metabolism in liver ####
# scaling factors
VLS9 = 143  # liver S9 protein yield (mg S9 protein/gram liver) reference: Punt et al. (2008)
# Part 1: Zearalenone phase I metabolism α-ZEL and β-ZEL formation
# metabolites α-ZEL and ß-ZEL  unscaled maximum rates of metabolism  (pmol/min/mg S9 protein)
VmaxLaZELc = 32  # reference: Malekinejad et al. (2006)
VmaxLbZELc = 72  # reference: Malekinejad et al. (2006)

# metabolites α-ZEL and ß-ZEL  affinity constants  (μmol/L)
KmLaZEL = 592  # reference: Malekinejad et al. (2006)
KmLbZEL = 21  # reference: Malekinejad et al. (2006)
# Part 2: Zearalenone phase II metabolism-glucuronidation
# metabolites ZENGLU  unscaled maximum rates of metabolism  (nmol/min/mg S9 protein)
VmaxLZENGLUc = 7.029  # in vitro  S9 incubations from this study

# metabolites ZENGLU  affinity constants  (μmol/L)
KmLZENGLU = 6.753  # in vitro  S9 incubations from this study
# Part 3: α-ZEL phase II metabolism: glucuronidation
# metabolites α-ZELGLU  unscaled maximum rates of metabolism  (nmol/min/mg S9 protein)
VmaxLaZELGLUc =  6.963  # in vitro  S9 incubations from this study

# metabolites α-ZELGLU  affinity constants  (μmol/L)
KmLaZELGLU  = 7.426 # in vitro  S9 incubations from this study



Parm.c.rat <- as.list(
  c(
    BW = BW,
    VSIc = VSIc, 
    VLc = VLc, 
    VRc = VRc, 
    VSc = VSc, 
    VFc = VFc, 
    VBc = VBc, 
    QC = QC, 
    QSIc = QSIc, 
    QLc = QLc, 
    QRc = QRc, 
    QSc = QSc, 
    QFc = QFc, 
    PIZEN = PIZEN, 
    PLZEN = PLZEN, 
    PRZEN = PRZEN, 
    PSZEN = PSZEN,
    PFZEN = PFZEN, 
    PIaZEL = PIaZEL, 
    PLaZEL = PLaZEL, 
    PRaZEL = PRaZEL, 
    PSaZEL = PSaZEL, 
    PFaZEL = PFaZEL, 
    ksto = ksto, 
    PappCaco2ZEN = PappCaco2ZEN, 
    Vin = Vin, 
    SAin = SAin, 
    kin = kin, 
    kin7 = kin7, 
    SAinb = SAinb, 
    PappCaco2aZEL = PappCaco2aZEL, 
    KurZ = KurZ, 
    KfeZ = KfeZ, 
    KuraZ = KuraZ, 
    KfeaZ = KfeaZ, 
    S9SI =  S9SI, 
    VmaxSIZENGLUc =  VmaxSIZENGLUc, 
    KmSIZENGLU = KmSIZENGLU, 
    VMB = VMB, 
    VmaxLIaZELc = VmaxLIaZELc,   
    VmaxLIbZELc = VmaxLIbZELc, 
    KmLIaZEL = KmLIaZEL, 
    KmLIbZEL = KmLIbZEL, 
    VLS9 = VLS9, 
    VmaxLaZELc = VmaxLaZELc, 
    VmaxLbZELc = VmaxLbZELc, 
    KmLaZEL = KmLaZEL, 
    KmLbZEL = KmLbZEL, 
    VmaxLZENGLUc = VmaxLZENGLUc, 
    KmLZENGLU = KmLZENGLU, 
    VmaxLaZELGLUc =  VmaxLaZELGLUc, 
    KmLaZELGLU  = KmLaZELGLU
  ) 
) %>% as.data.frame()

Parm.rat <- Parm.c.rat %>% 
  pivot_longer(names_to = "Parameter", 
               values_to = "Value", 
               cols = everything())

write.csv(Parm.rat, here("Input","Parm.rat.csv"), row.names = FALSE)

# PBK Model ----

R.PBK.model <- function(t, state, parms){
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

# Initialise model ----


## Parameters ----

Parameters <- as.list(Parameters_df$Value)
names(Parameters) <- Parameters_df$Parameter
MW_ZEN <- pull(PhysChem[PhysChem$Substance == "Zearalenone", 
                        "MolecularWeight"])

## Exposure ---- 

Select_AdminType <- "PO" # IV # PO
Select_Author <- "Videmann" # Shin # Q.Lu # Videmann 

BW <- if(Select_Author == "Shin"){0.26
}else if(Select_Author == "Q.Lu"){0.2
}else if(Select_Author == "Videmann"){0.26}

# Oral dose
ODOSEmg <- if(Select_Author == "Shin"){8
}else if(Select_Author == "Q.Lu"){6.65
}else if(Select_Author == "Videmann"){1}                                        # mg/kg bw IV dose
ODOSEumol <- ((ODOSEmg*1000)/MW_ZEN)*BW                                         # μmol

# IV dose
IVDOSEmg <- if(Select_AdminType != "IV"){0
}else if(Select_Author == "Shin"){1 #4 #2 # 4 #8#1.13#8
}else if(Select_Author == "Q.Lu"){6.65}

IVDOSEumol <-((IVDOSEmg*1000)/MW_ZEN)*BW                                        # μmol

names(ODOSEumol) <- "ODOSEumol"
names(IVDOSEumol) <- "IVDOSEumol"
parm.c <- c(Parameters, c(ODOSEumol, IVDOSEumol))

timeSeq <- seq(0, 24, by=1/100)                                                 # per minute


## Initials ----

A_init <- c(
  Ast = ODOSEumol,
  Ain1 = 0,
  Ain2 = 0,
  Ain3 = 0,
  Ain4 = 0,
  Ain5 = 0,
  Ain6 = 0,
  Ain7 = 0,
  ASIZEN = 0,
  ASIZENGLU = 0,
  ALIZEN = 0,
  ALIaZEL = 0,
  ALIbZEL = 0,
  AfeZ = 0,
  AfeaZ = 0,
  ALZEN = 0,
  ALaZEL = 0,
  ALbZEL = 0,
  ALZENGLU = 0,
  ABZENGLU = 0,
  AZENGLUur = 0,
  ALaZELGLU = 0,
  AF = 0,
  AR = 0,
  AS = 0,
  AB = IVDOSEumol,
  AZur = 0,
  AFaZEL = 0,
  ARaZEL = 0,
  ASaZEL = 0,
  ABaZEL = 0,
  AaZur = 0,
  ATotal = 0
)

## Run the model ----
PBK.out <- lsoda(y = A_init,
                 times = timeSeq, 
                 func =  R.PBK.model,
                 parms = parm.c
)


PBK_OUTPUT <- as.data.frame(PBK.out)



# Results ####

# Mass Balance
MB.df <- PBK_OUTPUT %>% select(time, ATotal, ACalculated, MB)                   
MB.df$MB <- round(MB.df$MB, 10)                                                 # rounding significance points
MB.df$ERROR <- (MB.df$ATotal - MB.df$ACalculated) / MB.df$ACalculated * 100
MB.df$ERROR <- round(MB.df$ERROR, 10)                                           # rounding significance points
MB.plot <- ggplot(data = MB.df)+
  geom_line(aes(x = time, y = ERROR, color = "ERROR")) +
  geom_line(aes(x = time, y = MB, color = "MB")) +
  scale_color_manual(values = c("ERROR" = "blue", "MB" = "black"),
                     name = NULL) +
  # ylim(0,1) +
  CP_theme +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12)
  )+
  ylab("MB / ERROR")
MB.plot
ggsave(
  filename = here(OUTPUT,"MB.Rat.png"),
  dpi = 600,
  width = 11,
  height = 7,
  units = "cm"
)
