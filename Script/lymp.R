# --------------------------------------------------------------------------- #
# Script for the PBK Model of DIANA et al. 
# By: Chrysanthi Pachoulide
# Date: 23-06-2025
# --------------------------------------------------------------------------- #

rm(list=ls()) 

# Packages
library(here)
library(tidyverse)
library(deSolve)
library(PKNCA)
library(pracma)
library(showtext)
library(plotly)
font_add(family = "Garamond", regular = "GARA.TTF")
showtext_auto()

CP_theme <- theme_minimal() +
  theme(
    # plot.background = element_rect(fill = 'white', color = NA),
    # panel.background = element_rect(fill = 'white', color = NA),
    axis.title = element_text(size = 44, face = "bold"),
    axis.text = element_text(size = 40),
    plot.title = element_text(hjust = 0, size = 48, face = "bold"),
    plot.margin = margin(0.3, 0.3, 0.3, 0.3, "cm"), 
    axis.line = element_line(size = 0.1),
    panel.grid = element_line(size = 0.3),
    plot.caption = element_text(hjust = 0, size = 40, face= "italic"),
    plot.title.position = "plot", 
    plot.caption.position = "plot",
    legend.position = "none"
  ) 

# Set output storage directory

OUTPUT <- here("Output", format(Sys.Date(), "%Y-%m-%d")) #, 
# format(Sys.time(), "%H-%M-%S"))
dir.create(OUTPUT, recursive = TRUE)


# Input ####

# Parameters

Parameters <- read.csv(here("Input", "Parms.csv")) %>% 
  pivot_longer(Rat:Human, 
               names_to = "Species", values_to = "Value") %>% 
  filter(Species == "Rat") %>% 
  select(c(Parameter, Value))

MWZEN <- 318.37 # molecular weight ZEN
MWZENGLUC <- 494.493 # molecular weight ZENGLUC
MWaZEL <- 320.4 # molecular weight aZEL

# Exposure ####

Select_AdminType <- "PO"
Select_Author <- "Videmann"

# Oral dose
ODOSEmg <- if(Select_Author == "Shin"){8
}else if(Select_Author == "Q.Lu"){6.65
}else if(Select_Author == "Videmann"){1}
# ODOSEmg <-  8 # mg/kg bw oral dose, variable Shin is 8, Videmann = 1, Q.Lu = 6.65
ODOSEumol <- ((ODOSEmg*1000)/MWZEN)*Parameters[1, "Value"][[1]] # μmol

# IV dose
IVDOSEmg <- 0 # mg/kg bw IV dose, variable Shin is 8, Q.Lu = 6.65
IVDOSEumol <-((IVDOSEmg*1000)/MWZEN)*Parameters[1, "Value"][[1]] # μmol

Dosing <- data.frame(Parameter = c("ODOSEumol", "IVDOSEumol"),
                     Value = c(ODOSEumol, IVDOSEumol))

Parameters <- rbind(Parameters, Dosing)
parm.c <- as.list(setNames(Parameters$Value, Parameters$Parameter)) # %>% unlist()

timeSeq <- seq(0, 24, by=1/60)

# Initials ####

A_init <- if(Select_AdminType == "PO"){
  
  c(AST_ZEN = ODOSEumol, 
    
    ASIL_ZEN = 0, 
    ASIL_ZENGLU = 0, 
    
    ACOL_ZEN = 0, ACOL_ZENGLU = 0, ACOL_aZEL = 0, ACOL_bZEL = 0, 
    AFE_ZEN = 0, AFE_aZEL = 0, AFE_ZENGLU = 0, 
    
    ASI_ZEN = 0,
    
    ASI_IC_ZENGLU = 0, ASI_EC_ZENGLU = 0, #AIle_IC_ZENGLU = 0, AIle_EC_ZENGLU = 0,
    ASI_aZEL = 0, ASI_bZEL = 0, ASI_aZELGLU = 0, #AJej_ZENGLU = 0, 
    
    AL_ZEN = 0, AL_IC_ZENGLU = 0, AL_EC_ZENGLU = 0, AL_bZEL = 0, AL_aZEL = 0, AL_aZELGLU = 0, #AL_ZENGLU = 0
    ABile_ZENGLU = 0,
    
    AB_ZEN = 0, AB_aZEL = 0, AB_ZENGLU = 0, 
    
    AUR_ZEN = 0, AUR_aZEL = 0, AUR_ZENGLU = 0, 
    AF_ZEN = 0, AF_aZEL = 0, AF_ZENGLU = 0, 
    AR_ZEN = 0, AR_aZEL = 0, AR_ZENGLU = 0, 
    AS_ZEN = 0, AS_aZEL = 0, AS_ZENGLU = 0
    
  )
  
}else{
  
  c(
    
    AST_ZEN = 0, 
    
    ADuoL_ZEN = 0, AJejL_ZEN = 0, AIleL_ZEN = 0, 
    ADuoL_ZENGLU = 0, AJejL_ZENGLU = 0, AIleL_ZENGLU = 0,
    
    ACOL_ZEN = 0, ACOL_ZENGLU = 0, ACOL_aZEL = 0, ACOL_bZEL = 0, 
    AFE_ZEN = 0, AFE_aZEL = 0, AFE_ZENGLU = 0,
    
    ADuo_ZEN = 0, AJej_ZEN = 0, AIle_ZEN = 0, 
    
    AJej_IC_ZENGLU = 0, AJej_EC_ZENGLU = 0, #AIle_IC_ZENGLU = 0, AIle_EC_ZENGLU = 0,
    AJej_aZEL = 0, AJej_bZEL = 0, AJej_aZELGLU = 0, #AJej_ZENGLU = 0,
    
    AL_ZEN = 0, AL_IC_ZENGLU = 0, AL_EC_ZENGLU = 0, AL_bZEL = 0, AL_aZEL = 0, AL_aZELGLU = 0, #AL_ZENGLU = 0
    ABile_ZENGLU = 0,
    
    AB_ZEN = IVDOSEumol, AB_aZEL = 0, AB_ZENGLU = 0, 
    
    AUR_ZEN = 0, AUR_aZEL = 0, AUR_ZENGLU = 0, 
    
    AF_ZEN = 0, AF_aZEL = 0, AF_ZENGLU = 0, 
    AR_ZEN = 0, AR_aZEL = 0, AR_ZENGLU = 0,
    AS_ZEN = 0, AS_aZEL = 0, AS_ZENGLU = 0
    
  )
}


# PBK model ####

R.PBK.model <- function(t, state, parms){
  with(as.list(c(state, parms)), {
    ## Calculated parameters 
    
    ## Physiological constants
    BW <- BW  
    
    VST <- FV_ST*BW  
    
    # Intestinal lumen
    VSIL <- FV_SIL*BW 
    VCOL <- FV_COL*BW 
    
    # Intestinal tissue
    VSI <- FV_SI*BW 
    VSI_IC <- FV_I_IC*VSI
    VSI_EC <- FV_I_EC*VSI
    VIle <- VSI*FV_Ile
    
    # Liver 
    VL <- FV_L*BW  # fraction of liver tissue
    VL_IC <- FV_L_IC*VL
    VL_EC <- FV_L_EC*VL
    
    # Rest
    VR <- FV_R*BW  
    # VS <- FV_S*BW  # 0.737-VFc # fraction of slowly perfused tissue
    VF <- FV_F*BW  # fraction of fat tissue
    VB <- FV_B*BW  # fraction of blood
    VS <- (0.96*BW) - (VSI + VL + VR + VF + VB) # fraction of slowly perfused tissue (Bone/Skin/Adipose/Muscle/(brain and testes) (slowly equilibrating tissues with large time constants))
    
    # Blood flow rates (reference: Punt et al. (2020))
    QC <- QC  ##{L/h} # cardiac output
    QSI <- QSIc*QC  # fraction of blood flow to small intestine
   
    # QCO <- 0 #QCOc*QC
    QL <- QLc*QC  # 0.242 - QSIc # fraction of blood flow to liver 101 
    QR <- QRc*QC  # 0.76 - QSIc - QLc # fraction of blood flow to rapidly perfused tissue
    # QS <- QSc*QC  # 0.24 - QFc # fraction of blood flow to slowly perfused tissue
    QF <- QFc*QC  # fraction of blood flow to fat
    QS <- QC - (QSI + QL + QR + QF) 
    
    # Glomerular filtration rate: Punt et al. (2020)
    GFR <- GFR # (L/h)
    
    Tsi <- 1/1.63 # small intestine transit time (/h)
    
    ## Kinetic parameters 
    ### Absorption/transfer rates
    R <- 0.18 #radius small intestine
    ka <- PappCaco2_ZEN/R
    
    fa <- 1 - (1 + ka/(1/1.63))**-7
   
    ### Metabolism in small intestine tissue
    # MISI <- 2.48
    VmaxSI_aZEL <- VmaxSI_aZELc * S9SI * VSI * 1e3       # umol/h, /fuincS9I.HSD_ZEN
    VmaxSI_bZEL <- VmaxSI_bZELc * S9SI * VSI * 1e3       # umol/h, /fuincS9I.HSD_ZEN
    VmaxSI_ZENGLU <- VmaxSI_ZENGLUc * S9SI * VSI * 1e3   # umol/h, /fuincS9I.UGT_ZEN 
    # VmaxIle_ZENGLU <- 0 #VmaxSI_ZENGLUc * S9SI * VIle * 1e3   # umol/h, /fuincS9I.UGT_ZEN 
    # VmaxCO_ZENGLU <- 0 #VmaxSI_ZENGLUc * S9SI * VCO * 1e3   # umol/h, /fuincS9I.UGT_ZEN 
    
    VmaxSI_aZELGLU <- VmaxSI_aZELGLUc * S9SI * VSI  * 1e3   # umol/h, /fuincS9I_aZEL
    
    
    ### metabolism in large intestine lumen (microbiota compartment)
    VmaxLI_aZEL <- 0 #VmaxLI_aZELc*VMB  # umol/h 
    VmaxLI_bZEL <- 0 #VmaxLI_bZELc*VMB
    
    ### Metabolism in liver
    # MILI <- 32
    VmaxL_aZEL <- VmaxL_aZELc * S9L * VL * 1e3 #/fuincS9L.HSD_ZEN
    VmaxL_bZEL <- VmaxL_bZELc * S9L * VL * 1e3 #/fuincS9L.HSD_ZEN
    VmaxL_ZENGLU <- VmaxL_ZENGLUc * S9L * VL * 1e3 #/fuincS9L.UGT_ZEN
    VmaxL_aZELGLU <- VmaxL_aZELGLUc * S9L * VL * 1e3 #/fuincS9L_aZEL 
    
    ### Biliary clearance
    Bile.fed <- 12.7e-6*60 # L/h
    Bile.fasted <- 6.35e-6*60 # L/h
    # CL_Bile <- ifelse(t > 3, Bile.fed, Bile.fasted) # 0.00075 L/h; derived from 50 and 70 μl/kg/min -> L/h
    CL_Bile <- ifelse(t > 3 && t < 20, Bile.fed, Bile.fasted) # 0.00075 L/h; derived from 50 and 70 μl/kg/min -> L/h
    # CL_Hydro <- 1000
    # CL_Hydro <- 1 #1e6
    ka_Lymph <- ifelse(t > 3, 0.33, 0) # 0.00075 L/h; derived from 50 and 70 μl/kg/min -> L/h
    
    
    CL.LI.MRP2 = CL.LI.MRP2.c*VL # L/h
    CL.LI.MRP3 = CL.LI.MRP3.c*VL # L/h
    CL.IN.MRP2 = CL.IN.MRP2.c*VSI # L/h
    # CL.IN.MRP3 = 0 #CL.IN.MRP3.c*VCO # L/h
    CL.KI.MRP3 = CL.KI.MRP3.c*VB # L/h
    
    # lilTR = lilTR
    
    ## Concentrations 
    
    # Zearalenone 
    CSI_ZEN <- ASI_ZEN/VSI
    CSIL_ZEN <- ASIL_ZEN/VSIL
    CCOL_ZEN <- ACOL_ZEN/VCOL 
    CL_ZEN <- AL_ZEN/VL
    CF_ZEN <- AF_ZEN/VF
    CR_ZEN <- AR_ZEN/VR
    CS_ZEN <- AS_ZEN/VS
    CB_ZEN <- AB_ZEN/VB
    
    CVSI_ZEN <- CSI_ZEN/PSI_ZEN * BP_ZEN
    CVL_ZEN <- CL_ZEN/PL_ZEN * BP_ZEN
    CVF_ZEN <- CF_ZEN/PF_ZEN * BP_ZEN
    CVR_ZEN <- CR_ZEN/PR_ZEN * BP_ZEN
    CVS_ZEN <- CS_ZEN/PS_ZEN * BP_ZEN
    
    # α-Zearalenol
    CCOL_aZEL <- ACOL_aZEL/VCOL
    CSI_aZEL <- ASI_aZEL/VSI
    CL_aZEL <- AL_aZEL/VL
    CF_aZEL <- AF_aZEL/VF
    CR_aZEL <- AR_aZEL/VR
    CS_aZEL <- AS_aZEL/VS
    CB_aZEL <- AB_aZEL/VB
    
    CVSI_aZEL <- CSI_aZEL/PSI_aZEL * BP_aZEL
    CVL_aZEL <- CL_aZEL/PL_aZEL * BP_aZEL
    CVF_aZEL <- CF_aZEL/PF_aZEL * BP_aZEL
    CVR_aZEL <- CR_aZEL/PR_aZEL * BP_aZEL
    CVS_aZEL <- CS_aZEL/PS_aZEL * BP_aZEL
    
    # β-Zearalenol
    CCOL_bZEL <- ACOL_bZEL/VCOL
    CL_bZEL <- AL_bZEL/VL
    
    # Zearalenone-glucuronide
    CSIL_ZENGLU <- ASIL_ZENGLU/VSIL
    
    CSI_IC_ZENGLU <- ASI_IC_ZENGLU/VSI_IC
    CSI_EC_ZENGLU <- ASI_EC_ZENGLU/VSI_EC
    CSIL_ZENGLU <- ASIL_ZENGLU/VSIL
    
    # CIle_IC_ZENGLU <- AIle_IC_ZENGLU/VIle_IC
    # CIle_EC_ZENGLU <- AIle_EC_ZENGLU/VIle_EC
    CSIL_ZENGLU <- ASIL_ZENGLU/VSIL 
    
    # CCO_IC_ZENGLU <- ACO_IC_ZENGLU/VCO_IC 
    # CCO_EC_ZENGLU <- ACO_EC_ZENGLU/VCO_EC 
    # 
    CCOL_ZENGLU <- ACOL_ZENGLU/VCOL 
    
    # CL_ZENGLU <- AL_ZENGLU/VL
    CL_IC_ZENGLU <- AL_IC_ZENGLU/VL_IC
    CL_EC_ZENGLU <- AL_EC_ZENGLU/VL_EC
    
    CF_ZENGLU <- AF_ZENGLU/VF
    CR_ZENGLU <- AR_ZENGLU/VR
    CS_ZENGLU <- AS_ZENGLU/VS
    CB_ZENGLU <- AB_ZENGLU/VB
    
    VBile <- 16000e-6 # L from ul
    CBile_ZENGLU <- ABile_ZENGLU/VBile
    
    # CVJej_ZENGLU <- CJej_ZENGLU/PSI_ZENGLU * BP_ZENGLU
    CVSI_EC_ZENGLU <- CSI_EC_ZENGLU/PSI_ZENGLU * BP_ZENGLU
    # CVIle_EC_ZENGLU <- CIle_EC_ZENGLU/PSI_ZENGLU * BP_ZENGLU
    # CVCO_EC_ZENGLU <- CCO_EC_ZENGLU/PSI_ZENGLU * BP_ZENGLU
    
    # CVL_ZENGLU <- CL_ZENGLU/PL_ZENGLU * BP_ZENGLU
    CVL_EC_ZENGLU <- CL_EC_ZENGLU/PL_ZENGLU * BP_ZENGLU
    
    CVF_ZENGLU <- CF_ZENGLU/PF_ZENGLU * BP_ZENGLU
    CVR_ZENGLU <- CR_ZENGLU/PR_ZENGLU * BP_ZENGLU
    CVS_ZENGLU <- CS_ZENGLU/PS_ZENGLU * BP_ZENGLU
    
    # a-Zearalenol glucuronide
    CL_aZELGLU <- AL_aZELGLU/VL   # liver a zel gluc
    CSI_aZELGLU <- ASI_aZELGLU/VSI
    
    ## Dynamic 
    
    # Stomach (in) 
    dAST_ZEN <- - gaTR*AST_ZEN 
    
    # Small intestinal lumen 
    dASIL_ZEN <- + gaTR*AST_ZEN +
      - fa*ASIL_ZEN +
      - (1-fa)*ka_Lymph*ASIL_ZEN +
      - Tsi*ASIL_ZEN 
    
    dASIL_ZENGLU <- CL.LI.MRP2*CL_IC_ZENGLU*fut_ZENGLU*CL_Bile + 
      - Tsi*ASIL_ZENGLU 
    
    dACOL_ZEN <- + Tsi*ASIL_ZEN + #*(1-fa_ile) +
      + CL_Hydro*CCOL_ZENGLU + 
      - lilTR*ACOL_ZEN
    
    dACOL_ZENGLU <- Tsi*ASIL_ZENGLU  +
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
    dASI_ZEN <- fa*ASIL_ZEN * QSI*(CB_ZEN-CVSI_ZEN) +
      - ((VmaxSI_ZENGLU*CSI_ZEN*fut_ZEN) / ((KmSI_ZENGLU) + (CSI_ZEN*fut_ZEN))) +
      - ((VmaxSI_aZEL*CSI_ZEN*fut_ZEN) / (KmSI_aZEL + (CSI_ZEN*fut_ZEN))) +
      - ((VmaxSI_bZEL*CSI_ZEN*fut_ZEN) / (KmSI_bZEL + (CSI_ZEN*fut_ZEN)))
    
    dASI_IC_ZENGLU <- ((VmaxSI_ZENGLU*CSI_ZEN*fut_ZEN) / 
                         ((KmSI_ZENGLU) + (CSI_ZEN*fut_ZEN))) +
      - CL.IN.MRP2*CSI_IC_ZENGLU*fut_ZENGLU 
    
    dASI_EC_ZENGLU <- + QSI*(CB_ZENGLU-CVSI_EC_ZENGLU) #CL.IN.MRP3*CJej_IC_ZENGLU*fut_ZENGLU
    
    dASI_aZEL <- ((VmaxSI_aZEL*CSI_ZEN*fut_ZEN) / (KmSI_aZEL + (CSI_ZEN*fut_ZEN))) +
      + QSI*(CB_aZEL-CVSI_aZEL) +
      - ((VmaxSI_aZELGLU*CSI_aZEL*fut_aZEL) / (KmSI_aZELGLU + (CSI_aZEL*fut_aZEL)))
    
    dASI_bZEL <- ((VmaxSI_bZEL*CSI_ZEN*fut_ZEN) / (KmSI_bZEL + (CSI_ZEN*fut_ZEN)))
    
    dASI_aZELGLU <- ((VmaxSI_aZELGLU*CSI_aZEL*fut_aZEL) / (KmSI_aZELGLU + (CSI_aZEL*fut_aZEL)))
    
    # Liver 
    dAL_ZEN <- QSI*CVSI_ZEN + QL*CB_ZEN + #ka_lil*ACO_ZEN*fa_lil + 
      - (QSI+QL)*CVL_ZEN + 
      - ((VmaxL_aZEL*CL_ZEN*fut_ZEN) / (KmL_aZEL + (CL_ZEN*fut_ZEN))) +
      - ((VmaxL_bZEL*CL_ZEN*fut_ZEN) / (KmL_bZEL + (CL_ZEN*fut_ZEN))) +
      - ((VmaxL_ZENGLU*CL_ZEN*fut_ZEN) / ((KmL_ZENGLU) + (CL_ZEN*fut_ZEN)))
    
    dAL_IC_ZENGLU <- + ((VmaxL_ZENGLU*CL_ZEN*fut_ZEN) / ((KmL_ZENGLU) + (CL_ZEN*fut_ZEN))) +
      - CL.LI.MRP2*CL_IC_ZENGLU*fut_ZENGLU*CL_Bile + #(Vmax_BSEP/(Km_BSEP + (CL_ZENGLU*fut_ZENGLU)))*CL_ZENGLU*fut_ZENGLU
      - CL.LI.MRP3*CL_IC_ZENGLU*fut_ZENGLU
    
    dAL_EC_ZENGLU <- + QSI*CVSI_EC_ZENGLU + #QIle*CVIle_EC_ZENGLU + #QCO*CVCO_EC_ZENGLU + 
      + (QL+ QSI)*CB_ZENGLU - (QSI + QL)*CVL_EC_ZENGLU +
      + CL.LI.MRP3*CL_IC_ZENGLU*fut_ZENGLU
    
    dAL_bZEL <- ((VmaxL_bZEL*CL_ZEN*fut_ZEN) / (KmL_bZEL + (CL_ZEN*fut_ZEN))) 
    
    dAL_aZEL <- + QSI*CSI_aZEL + QL*CB_aZEL - (QSI+QL)*CVL_aZEL + 
      + ((VmaxL_aZEL*CL_ZEN*fut_ZEN) / (KmL_aZEL + (CL_ZEN*fut_ZEN))) +
      - ((VmaxL_aZELGLU*CL_aZEL*fut_aZEL) / (KmL_aZELGLU + (CL_aZEL*fut_aZEL))) 
    
    dAL_aZELGLU <- ((VmaxL_aZELGLU*CL_aZEL*fut_aZEL) / (KmL_aZELGLU + (CL_aZEL*fut_aZEL))) 
    
    # Bile
    dABile_ZENGLU <- 0 #CL.LI.MRP2*CL_IC_ZENGLU*fut_ZENGLU + #CL.LI.MRP2*CL_ZENGLU*fut_ZENGLU + #(Vmax_BSEP/(Km_BSEP + (CL_ZENGLU*fut_ZENGLU)))*CL_ZENGLU*fut_ZENGLU +
    # - CL_Bile*CBile_ZENGLU #ABile_ZENGLU #CBile_ZENGLU
    
    # Blood
    dAB_ZEN <- - (QSI + QL + QR + QS + QF)*CB_ZEN +  #
      + (QSI + QL)*CVL_ZEN + QF*CVF_ZEN + QR*CVR_ZEN + QS*CVS_ZEN + 
      - GFR*CB_ZEN*fup_ZEN +
      + ka_Lymph*ASIL_ZEN  
    
    dAB_aZEL <- - (QSI + QL + QR + QS + QF)*CB_aZEL +  #
      + (QL+QSI)*CVL_aZEL + QF*CVF_aZEL + QR*CVR_aZEL + QS*CVS_aZEL +
      - GFR*CB_aZEL*fup_aZEL
    
    dAB_ZENGLU <- - (QSI + QL + QR + QS + QF)*CB_ZENGLU +   # 
      + (QSI + QL)*CVL_EC_ZENGLU + QF*CVF_ZENGLU + QR*CVR_ZENGLU + QS*CVS_ZENGLU +
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
    ATotal <- ODOSEumol + IVDOSEumol
    ACalculated <- AST_ZEN +
      ASIL_ZEN + ASIL_ZENGLU +
      ACOL_ZEN + ACOL_ZENGLU + ACOL_aZEL + ACOL_bZEL +
      AFE_ZEN + AFE_aZEL + AFE_ZENGLU + 
      ASI_ZEN + 
      ASI_IC_ZENGLU + ASI_EC_ZENGLU + #AIle_IC_ZENGLU + AIle_EC_ZENGLU +
      ASI_aZEL + ASI_bZEL + ASI_aZELGLU + #AJej_ZENGLU
      AL_ZEN + AL_IC_ZENGLU + AL_EC_ZENGLU + AL_bZEL + AL_aZEL + AL_aZELGLU + #AL_ZENGLU
      ABile_ZENGLU +
      AB_ZEN + AB_aZEL + AB_ZENGLU +
      AUR_ZEN + AUR_aZEL + AUR_ZENGLU + 
      AF_ZEN + AF_aZEL + AF_ZENGLU + 
      AR_ZEN + AR_aZEL + AR_ZENGLU + 
      AS_ZEN + AS_aZEL + AS_ZENGLU
    
    
    MB <- ATotal - ACalculated + 1
    
    # End
    list(c(
      dAST_ZEN,
      dASIL_ZEN, 
      dASIL_ZENGLU, 
      dACOL_ZEN, dACOL_ZENGLU, dACOL_aZEL, dACOL_bZEL, 
      dAFE_ZEN, dAFE_aZEL, dAFE_ZENGLU,
      dASI_ZEN, 
      dASI_IC_ZENGLU, dASI_EC_ZENGLU, #dAIle_IC_ZENGLU, dAIle_EC_ZENGLU,
      dASI_aZEL, dASI_bZEL, dASI_aZELGLU,#dAJej_ZENGLU
      dAL_ZEN, dAL_IC_ZENGLU, dAL_EC_ZENGLU, dAL_bZEL, dAL_aZEL, dAL_aZELGLU, #dAL_ZENGLU
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
      CSIL_ZEN = CSIL_ZEN,
      CSIL_ZENGLU = CSIL_ZENGLU, 
      CCOL_ZEN = CCOL_ZEN, CCOL_ZENGLU = CCOL_ZENGLU, CCOL_aZEL = CCOL_aZEL, CCOL_bZEL = CCOL_bZEL,
      CSI_ZEN = CSI_ZEN,
      CSI_IC_ZENGLU = CSI_IC_ZENGLU, CSI_EC_ZENGLU = CSI_EC_ZENGLU, 
      CSI_aZEL = CSI_aZEL, CSI_aZELGLU = CSI_aZELGLU, # CJej_ZENGLU = CJej_ZENGLU; (((CJej_IC_ZENGLU/VJej_IC) + (CJej_EC_ZENGLU/VJej_EC))/(VJej_IC + VJej_EC)); CIle_ZENGLU = CIle_IC_ZENGLU, 
      CL_ZEN = CL_ZEN, CL_IC_ZENGLU = CL_IC_ZENGLU, CL_EC_ZENGLU = CL_EC_ZENGLU,
      CL_bZEL = CL_bZEL, CL_aZEL = CL_aZEL, CL_aZELGLU = CL_aZELGLU, #CL_ZENGLU = CL_ZENGLU; (((CL_IC_ZENGLU/VL_IC) + (CL_EC_ZENGLU/VL_EC))/(VL_IC + VL_EC))
      CF_ZEN = CF_ZEN, CF_aZEL = CF_aZEL, CF_ZENGLU = CF_ZENGLU,
      CR_ZEN = CR_ZEN, CR_aZEL = CR_aZEL, CR_ZENGLU = CR_ZENGLU, 
      CS_ZEN = CS_ZEN, CS_aZEL = CS_aZEL, CS_ZENGLU = CS_ZENGLU,
      CB_ZEN = CB_ZEN, CB_aZEL = CB_aZEL, CB_ZENGLU = CB_ZENGLU
    ))
    
  })
}

PBK.out <- lsoda(y = A_init,
                 times = timeSeq,#per minute
                 func =  R.PBK.model,
                 parms = parm.c
)


PBK_OUTPUT <- as.data.frame(PBK.out)

# save(PBK_OUTPUT, file = here(OUTPUT, "RawData_PBK.RAT.RData"))


# Post run analysis ####
MB.df <- PBK_OUTPUT %>% select(time, ATotal, ACalculated, MB)  # change days to years if needed
MB.df$MB <- round(MB.df$MB, 10) # rounding significance points
MB.df$ERROR <- (MB.df$ATotal - MB.df$ACalculated) / MB.df$ACalculated * 100
MB.df$ERROR <- round(MB.df$ERROR, 10) # rounding significance points
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
# ggsave(
#   filename = here(OUTPUT,"MB.Rat.png"),
#   dpi = 600,
#   width = 11,
#   height = 7,
#   units = "cm"
# )

BW <- parm.c$BW  ##{Kg} #  body weight rat
VST <-parm.c$FV_ST*BW  # fraction volume of small intestine
VSIL <- parm.c$FV_SIL*BW # fractional volume of small intestinal lumen
# VCOL <- (parm.c$FV_LIL*BW)-VREL
VSI <- parm.c$FV_SI*BW # fractional volume of small intestinal tissue

VSI_IC <- parm.c$FV_I_IC*VSI
VSI_EC <- parm.c$FV_I_EC*VSI

C.Pre <- PBK_OUTPUT %>%
  transmute(Time = time,
            Intestine_ZEN = ((CSI_ZEN/VSI) + (CSIL_ZEN/VSIL))/(VSI + VSIL),
            Liver_ZEN = CL_ZEN, 
            Fat_ZEN = CF_ZEN, 
            Rapidly_ZEN = CR_ZEN, Slowly_ZEN = CS_ZEN,
            Blood_ZEN = CB_ZEN, 
            
            Intestine_aZEL = 0,  
            Liver_aZEL = CL_aZEL,
            Fat_aZEL = CF_aZEL, Rapidly_aZEL = CR_aZEL, Slowly_aZEL = CS_aZEL,
            Blood_aZEL = CB_aZEL,
            
            Intestine_ZENGLU = ((CSI_IC_ZENGLU/VSI_IC) + (CSI_EC_ZENGLU/VSI_EC) +
                                (CSIL_ZENGLU/VSIL))/ (VSI_IC + VSI_EC + VSIL),
            LIntestineL_ZENGLU = CCOL_ZENGLU, 
            LiverIC_ZENGLU = CL_IC_ZENGLU, LiverEC_ZENGLU  = CL_EC_ZENGLU, 
            Fat_ZENGLU = CF_ZENGLU, Rapidly_ZENGLU = CR_ZENGLU, Slowly_ZENGLU = CS_ZENGLU,
            Blood_ZENGLU = CB_ZENGLU
  )

C.Pre_long <- C.Pre %>%
  pivot_longer(names_to = c("Tissue", "Compound"),
               IntestinalT_ZEN:Blood_ZENGLU,
               values_to = "Concentration",
               names_sep = "_") %>%
  mutate(Concentration = ifelse(Compound == "ZEN", Concentration*MWZEN, # correct units
                                ifelse(Compound == "ZENGLU",Concentration*MWZENGLUC,
                                       Concentration*MWaZEL)))



C.Pre_wide <- C.Pre_long %>%
  pivot_wider(names_from = "Tissue", values_from = "Concentration") %>%
  select(c("Compound","Time", 
           "IntestinalT", "IntestinalL", 
           "LIntestineL", "Liver", "Blood", "Fat", "Rapidly", "Slowly")) %>%
  # mutate("IntestineL" = (((DuodenumL*VDuoL) + (JejunumL*VJejL) + (IleumL*VIleL) + (LIntestineL*VCOL))/(VDuoL+VJejL+VIleL+VCOL))) %>%
  # mutate("Intestine" = (((SIntestineT*VSI) + (IntestineL*(VSIL+VCOL+VREL)))/(VSI + VSIL+VCOL+VREL))) %>%
  mutate("Intestine" = ((IntestinalT*VSI) + (IntestinalL*VSIL))/(VSI + VSIL)) %>% #L+VCOL+VRE  +VDuoL+VJejL+VIle
  # rename("Intestine" = "SIntestineT") %>%
  select(c("Compound", "Time", "Intestine", "Liver", "Blood", "Fat", "Rapidly", "Slowly")) 

# Plot Blood
CBlood.plot <- C.Pre_long %>% filter(Tissue == "Blood") %>%
  ggplot() + geom_path(aes(x = Time, y = Concentration)) +
  facet_wrap(~ Compound, scales = "free") +
  CP_theme +
  theme(strip.text = element_text(size = 45),
        axis.text = element_text(size = 40)) +
  ylab("Blood Concentration (ng/ml)") +
  xlab("Hours")+
  theme(legend.position = "bottom",
        legend.title = element_text(size = 35),
        legend.text = element_text(size = 35),
        legend.key.size = unit(0.5, "lines"),
        legend.spacing.y = unit(0, "cm"),
        legend.margin = margin(0, 0, 0, 0))
CBlood.plot
# ggsave(
#   filename = here(OUTPUT,"BloodC.Rat.png"),
#   dpi = 600,
#   width = 11,
#   height = 7,
#   units = "cm"
# )

# Plot all tissues
C.Tissues <- C.Pre_wide %>% 
  pivot_longer(Intestine:Slowly, 
               names_to = "Tissue", values_to = "Concentration")

CTissues.plot <- C.Tissues %>% 
  ggplot() + geom_path(aes(x = Time, y = Concentration, colour = Compound)) +
  scale_color_brewer(palette = "Dark2") +
  facet_wrap(~ Tissue, scales = "free") +
  CP_theme +
  theme(strip.text = element_text(size = 45),
        axis.text = element_text(size = 40)) +
  ylab("Concentration (ng/ml)") +
  xlab("Hours")+
  theme(legend.position = "bottom",
        legend.title = element_text(size = 35),
        legend.text = element_text(size = 35),
        legend.key.size = unit(0.5, "lines"),
        legend.spacing.y = unit(0, "cm"),
        legend.margin = margin(0, 0, 0, 0))
CTissues.plot
# ggsave(
#   filename = here(OUTPUT,"TissueC.Rat.png"),
#   dpi = 600,
#   width = 11,
#   height = 7,
#   units = "cm"
# )


# C.Tissues %>% 
#   ggplot() + geom_path(aes(x = Time, y = log10(Concentration), colour = Compound)) +
#   scale_color_brewer(palette = "Dark2") +
#   facet_wrap(~ Tissue, scales = "free") +
#   CP_theme +
#   theme(strip.text = element_text(size = 45),
#         axis.text = element_text(size = 40)) +
#   ylab("Concentration log(ng/ml)") +
#   xlab("Hours")+
#   theme(legend.position = "bottom",
#         legend.title = element_text(size = 35),
#         legend.text = element_text(size = 35),
#         legend.key.size = unit(0.5, "lines"),
#         legend.spacing.y = unit(0, "cm"),
#         legend.margin = margin(0, 0, 0, 0))
# 
# ggsave(
#   filename = here(OUTPUT,"TissueC.log.Rat.png"),
#   dpi = 600,
#   width = 11,
#   height = 7,
#   units = "cm"
# )

# write_csv(C.Tissues, here(OUTPUT, "C.Tissues.Rat.csv"))

Results.Rat <- list(
  Metadata.Run.Rat = list(Select_AdminType = Select_AdminType,
                          Select_Author = Select_Author,
                          ODOSEmg = ODOSEmg, 
                          ODOSEumol = ODOSEumol,
                          IVDOSEmg = IVDOSEmg,
                          IVDOSEumol = IVDOSEumol,
                          timeSeq = timeSeq,
                          MWZEN = 318.37, # molecular weight ZEN
                          MWZENGLUC = 494.493, # molecular weight ZENGLUC
                          MWaZEL = 320.4, # molecular weight aZEL
                          BP = 0.67 # assuming the same for all ZEN, aZEL and ZENGLUC
  ), 
  Raw.Data.Rat = PBK_OUTPUT,
  C.Tissues.Rat = C.Tissues,
  Plots = list(MB.plot, CBlood.plot, CTissues.plot)
)

save(Results.Rat, file = here(OUTPUT, "Results.Rat.RData"))
