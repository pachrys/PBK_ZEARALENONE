# --------------------------------------------------------------------------- #
# Script for calculating parameters  
# exact copy of the calculations of the PBK Model of Mendez-Catala et al. 2021 DOI: 10.1002/mnfr.202100443
# By: Chrysanthi Pachoulide
# Date: 23-06-2025
# --------------------------------------------------------------------------- #

# --------------------------------------------------------------------------- #
# Script for the PBK Model of DIANA et al. 
# By: Chrysanthi Pachoulide
# Date: 23-06-2025
# --------------------------------------------------------------------------- #

rm(list = ls()) 

# Packages
library(here)
library(tidyverse)
library(deSolve)
library(PKNCA)
library(pracma)
library(showtext)
font_add(family = "Garamond",  regular = "GARA.TTF")
showtext_auto()

# Set output storage directory
OUTPUT = here("Output",  format(Sys.Date(),  "%Y-%m-%d"),  format(Sys.time(),  "%H-%M-%S"))
dir.create(OUTPUT,  recursive = TRUE)


# Parameters for Rat####

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
  

# Parameters for human ####

  ## Physiological parameters ####
  # Tissue volumes, reference: Brown et al. (1997)
  BW <- 70 #{Kg} body weight rat
  VSIc <- 0.009 # fraction of small intestine
  VLc <- 0.026 # fraction of liver tissue
  VRc <- 0.041 # 0.076-VSIc-VLc # fraction of rapidly perfused tissue
  VSc <- 0.596 # 0.81-VFc # fraction of slowly perfused tissue
  VFc <- 0.214 # fraction of fat tissue
  VBc <- 0.074 # fraction of blood
  
  # Blood flow rates, reference: Brown et al. (1997)
  QC <- 347.9 #{L/h} # cardiac output:15*BW^0.74
  QSIc <- 0.181 # fraction of blood flow to small intestine
  QLc <- 0.046 # 0.227 - QSIc # fraction of blood flow to liver
  QRc <- 0.473 # 0.7 - QSIc – QLc # fraction of blood flow to rapidly perfused tissue
  QSc <- 0.248 # 0.3 - QFc # fraction of blood flow to slowly perfused tissue
  QFc <- 0.052 # fraction of blood flow to fat
  
  ## Physicochemical parameters ####
  ### partition coefficients ####
  # calculated based on QPPR of DeJongh et al. (1997)
  
  # ZEN in main model
  PIZEN <- 6.56 # intestine/blood partition coefficient
  PLZEN <- 6.56 # liver/blood partition coefficient
  PRZEN <- 6.56 # rapidly perfused tissue/blood partition coefficient
  PSZEN <- 4.25 # slowly perfused tissue/blood partition coefficient
  PFZEN <- 134.86 # fat/blood partition coefficient
  
  # α-ZEL in sub-model
  PIaZEL <- 6.11 # intestine/blood partition coefficient
  PLaZEL <- 6.11 # liver/blood partition coefficient
  PRaZEL <- 6.11 # rapidly perfused tissue/blood partition coefficient
  PSaZEL <- 3.99 # slowly perfused tissue/blood partition coefficient
  PFaZEL <- 131.62 # fat/blood partition coefficient
  
  ### absorption/transfer rates ####
  # Stomach emptying rate
  ksto <- 2.8 # stomach emptying rate (/h) reference: Reilly et al. (1990)
  
  # intestinal absorption and transfer rates of ZEN
  #Papp, Caco-2 ZEN <- 10.3 #( x10-5 cm/sec) fitted apparent in vitro, from Caco-2 reference: Pfeiffer et al.(2011)
  PappCaco2ZEN <- 4.9829 # Log Papp, Caco-2
  # Log (Papp,in vivo) <- 0.6836*Log(PappCaco2ZEN)-0.5579 reference: Sun et al. (2002)
  
  Vin <- 1.29 # volume for each compartment of intestines (L)
  SAin <- 10.3 # surface area (dm2)
  kin <- 2.19 # transfer rate to next compartment within the intestines (/hr)
  
  kin7 <- 0.464 # transfer rate of ZEN from small intestine to large intestine (/h) reference: Kimura and Higaki (2002)
  SAinb <- 47.12 # surface area of large intestine (dm2) reference: Vdoviaková et al. (2016)
  
  # α-ZEL
  # Papp, Caco-2 α-ZEL <- 5.4 (x10-6 cm/sec), apparent in vitro, from Caco-2 reference: Pfeiffer et al. (2011)
  PappCaco2aZEL<--5.2676 # Log Papp, Caco-2
  #  Log (Papp,in vivo) <- 0.6836*Log(PappCaco2aZEL)-0.5579 reference: Sun et al. (2002)
 
  # excretion (/h)
  KurZ <- 0.096 # urinary excretion of ZEN reference: Mukherjee et al. (2014)
  KfeZ <- 0.024 # fecal excretion of ZEN reference: Mukherjee et al. (2014)
  KuraZ <- 0.015 # urinary excretion of α-ZEL reference: Mukherjee et al. (2014)
  KfeaZ <- 0.024 # fecal excretion of α-ZEL reference: Mukherjee et al. (2014)
  KurG <- 59.4 # urinary excretion of α-ZEL reference: Mukherjee et al. (2014)
  
  ##  Kinetic parameters ####
  
  ### metabolism of small intestine tissue ####
  # scaling factors
  S9SI<- 35.2 # small intestinal S9 protein yield (mg S9 protein/gram intestine) reference: (Peters et al., 2016)
  # metabolites ZENGLU, unscaled maximum rates of metabolism, (nmol/min/mg S9 protein)
  VmaxSIZENGLUc<- 0.491 # in vitro, S9 incubations from this study
 
  # metabolites ZENGLU, affinity constants, (μmol/L)
  KmSIZENGLU <- 1.174 # in vitro, S9 incubations from this study
  
  ### metabolism of large intestine lumen (microbiota compartment) ####
  # scaling factors
  VMB <- 0.0018 # fraction of faeces of BW, reference: Rose et al. (2015)
  #  metabolites α-ZEL and ß-ZEL, unscaled maximum rates of metabolism, (pmol/h/g faeces)
  VmaxLIaZELc<-0.054 # in vitro, anaerobic rat fecal incubations
  VmaxLIbZELc<- 0.0108 # in vitro, anaerobic rat fecal incubations
 
  #  metabolites α-ZEL and ß-ZEL, affinity constants, (μmol/L)
  KmLIaZEL <- 135 # in vitro, anaerobic rat fecal incubations
  KmLIbZEL <- 163 # in vitro, anaerobic rat fecal incubations
  
  ### metabolism of liver ####
  # scaling factors
  VLS9 <- 120.7 #  liver S9 protein yield, (mg S9 protein/gram liver) reference: Cubitt et al. (2011)
  
  # Part 1: Zearalenone phase I metabolism α-ZEL and ß-ZEL formation
  # metabolites α-ZEL and ß-ZEL, unscaled maximum rates of metabolism, (pmol/min/mg S9 protein)
  VmaxLaZELc <- 358.7 # reference: Malekinejad et al. (2006)
  VmaxLbZELc <- 209.3 # reference: Malekinejad et al. (2006)
  
  #  metabolites α-ZEL and ß-ZEL, affinity constants, (μmol/L)
  KmLaZEL <- 9 # reference: Mendez-Catala et al. (2020)
  KmLbZEL <- 23 # reference: Mendez-Catala et al. (2020)
  
  #  Part 2: Zearalenone phase II metabolism-glucuronidation
  #  metabolites ZENGLU, unscaled maximum rates of metabolism, (nmol/min/mg S9 protein)
  VmaxLZENGLUc <- 2.97 # in vitro, S9 incubations from this study

  #  metabolites ZENGLU, affinity constants, (μmol/L)
  KmLZENGLU <- 2.043 # in vitro, S9 incubations from this study
  #  Part 3: α-ZEL phase II metabolism: glucuronidation
  #  metabolites α-ZELGLU, unscaled maximum rates of metabolism, (nmol/min/mg S9 protein)
  VmaxLaZELGLUc<- 2.983 # in vitro, S9 incubations from this study

  #  metabolites α-ZELGLU, affinity constants, (μmol/L)
  KmLaZELGLU <-2.415 # in vitro, S9 incubations from this study
  
  
  Parm.c.human <- as.list(
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
      S9SI = S9SI,
      VmaxSIZENGLUc = VmaxSIZENGLUc,
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
      VmaxLaZELGLUc = VmaxLaZELGLUc,
      KmLaZELGLU = KmLaZELGLU
    )
  ) %>% as.data.frame()
  
Parm.human <- Parm.c.human %>% 
  pivot_longer(names_to = "Parameter", 
               values_to = "Value", 
               cols = everything())

write.csv(Parm.human, here("Input","Parm.human.csv"), row.names = FALSE)
