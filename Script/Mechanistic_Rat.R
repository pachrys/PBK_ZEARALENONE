# --------------------------------------------------------------------------- #
# Script for the Refined PBK Model 
# By: Chrysanthi Pachoulide
# Date: 01-05-2025
# --------------------------------------------------------------------------- #

rm(list=ls()) 

# Packages
library(here)
library(tidyverse)
library(deSolve)
library(showtext)
library(plotly)
library(readr)
library(PKNCA)
library(pracma)
font_add(family = "Garamond", regular = "GARA.TTF")
showtext_auto()

# For plotting
test = "IV_Q.Lu"

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

timeSeq <- seq(0, 24, by=1/100)                                                  # per minute



# Calculate PBK model parameters ####
# -------------------------------------------------------------------------#

# Rat physiology
Rat.input <- read_delim(here("Input", "TissueComposition.csv"),
                        delim = ";",
                        escape_double = FALSE, trim_ws = TRUE) %>%
  filter(Species == "Rat")                                                      # tissue composition mainly as R&R composition, unified by Utsey and optimised as per Musther et al. (2017) DOI: 10.1016/j.xphs.2017.05.003

Rat.dat <- Rat.input %>% filter(!tissue %in% c("RBCs", "Plasma"))          
names <- substr(Rat.dat$tissue, 1, 2)
F_Vo.names <- paste0("FV", names)
F_Vo.dat <- cbind(F_Vo.names, select(Rat.dat, F_Vo)) %>%
  pivot_wider(names_from = F_Vo.names, values_from = F_Vo)
# Bile_Lymph.df <- read_delim(here("Input", "Bile_LymphN.csv"),
#                             delim = ",",
#                             escape_double = FALSE, trim_ws = TRUE)

## Physiological constants ----
# Reference (if not stated otherwise): Punt et al. (2020) https://dx.doi.org/10.1021/acs.chemrestox.0c00307

# Body weight 
BW <- 0.250                                                                     # Kg

### Fractional volumes ----           

FV_SI <- pull(Rat.input[Rat.input$tissue == "Gut", "F_Vo"])                     # small intestinal tissue
FV_Duo <- 0.060                                                                 # fraction of small intestine that is duodenum
FV_Jej <- 0.769                                                                 # fraction of small intestine that is jejunum
FV_Ile <- 0.171                                                                 # fraction of small intestine that is ileum

FV_SI_IC <- pull(Rat.input[Rat.input$tissue == "Gut", "f_iw"])                  # fraction of small intestine that is intracellular water
FV_SI_EC <- pull(Rat.input[Rat.input$tissue == "Gut", "f_ew"])                  # fraction of small intestine that is extracellular water

FV_SIL <- 0.014                                                                 # small intestine lumen 
FV_duol <- 0.183                                                                # duodenum lumen; calculated from Musther et al. (2017) DOI: 10.1016/j.xphs.2017.05.003
FV_jejl <- 0.334                                                                # jejunum lumen; calculated from Musther et al. (2017) DOI: 10.1016/j.xphs.2017.05.003
FV_ilel <- 0.483                                                                # ileum lumen; calculated from Musther et al. (2017) DOI: 10.1016/j.xphs.2017.05.003

Vfe_IleL <- 4.1e-3 * 0.1                                                        # volume of "feces" in the ileum, 4.1g is the total volume of feces excreted per day, 0.007 is to correct for the small fraction already formed in the ileum

FV_L <- pull(Rat.input[Rat.input$tissue == "Liver", "F_Vo"])                    # liver tissue
FV_L_IC <- pull(Rat.input[Rat.input$tissue == "Liver", "f_iw"])                 # fraction of liver that is intracellular water
FV_L_EC <- pull(Rat.input[Rat.input$tissue == "Liver", "f_ew"])                 # fraction of liver that is extracellular water
FV_R <- pull(Rat.input[Rat.input$tissue == "Heart", "F_Vo"]) +
  pull(Rat.input[Rat.input$tissue == "Kidney", "F_Vo"]) +
  pull(Rat.input[Rat.input$tissue == "Spleen", "F_Vo"]) +
  pull(Rat.input[Rat.input$tissue == "Lung", "F_Vo"])                           # rapidly perfused tissues
FV_R_gluc <- pull(Rat.input[Rat.input$tissue == "Heart", "F_Vo"]) +
  pull(Rat.input[Rat.input$tissue == "Spleen", "F_Vo"]) +
  pull(Rat.input[Rat.input$tissue == "Lung", "F_Vo"])                           # rapidly perfused tissues
FV_F <- pull(Rat.input[Rat.input$tissue == "Adipose", "F_Vo"])                  # fat tissue
FV_B <- pull(Rat.input[Rat.input$tissue == "RBCs", "F_Vo"]) + 
  pull(Rat.input[Rat.input$tissue == "Plasma", "F_Vo"])                         # blood

FV_K <- pull(Rat.input[Rat.input$tissue == "Kidney", "F_Vo"])                   # kidney tissue

VBile <- 16e-3                                                                  # L from ml; Reference: Vonk, Doorn, and Strubbe (1978) DOI: 10.1042/cs0550253 

### To calculate the surface aeria of the intestine ----
# Reference (if not stated otherwise): Musther et al. (2017) DOI: 10.1016/j.xphs.2017.05.003 

# Radius
R <- 0.18                                                                       # cm; radius small intestine
Rduo <- 0.275/2                                                                 # cm; radius duodenum; (diameter / 2 = radius)
Rjej <- 0.4/2                                                                   # cm; radius jejunum; (diameter / 2 = radius)
Rile <- 0.4/2                                                                   # cm; radius ileum; (diameter / 2 = radius)

# Length
Lduo <- 7                                                                       # cm duodenum
Ljej <- 90                                                                      # cm jejunum
Lile <- 20                                                                      # cm ileum

# Villi amplification factor
# Reference: Desesso et al. (2001) https://doi.org/10.1016/S0278-6915(00)00136-8
# Reference: Stenling and Helander 1981 10.1007/BF00233821, rat villi amplification factor at a non-fasted state

Amp <- 10                                                                       # small intestine 
Ampduo <- 1*Amp                                                                 # duodenum       
Ampjej <- 0.8*Amp                                                               # jejunum
Ampile <- 0.5*Amp                                                               # ileum

# Surface aeria
SAduo = 2*pi*Rduo*Lduo*Ampduo                                                   # cm2; duodenum
SAjej = 2*pi*Rjej*Ljej*Ampjej                                                   # cm2; jejunum
SAile = 2*pi*Rile*Lile*Ampile                                                   # cm2; ileum


# Cardiac output
QC = 4.8                                                                        # L/h

### Fractional blood flow rates ----
# Reference: Punt et al. (2020) https://dx.doi.org/10.1021/acs.chemrestox.0c00307)

QSIc <- pull(Rat.input[Rat.input$tissue == "Gut", "F_QC"])                      # small intestine
QDuoc <- 0.079                                                                  # fraction of small intestinal flow rate that goes to the duodenum; Musther et al. (2017) DOI: 10.1016/j.xphs.2017.05.003
QJejc <- 0.680                                                                  # fraction of small intestinal flow rate that goes to the jejunum; Musther et al. (2017) DOI: 10.1016/j.xphs.2017.05.003
QIlec <- 0.241                                                                  # fraction of small intestinal flow rate that goes to the ileum; Musther et al. (2017) DOI: 10.1016/j.xphs.2017.05.003

QLc <- pull(Rat.input[Rat.input$tissue == "Liver", "F_QC"])                     # liver  
QKc <- pull(Rat.input[Rat.input$tissue == "Kidney", "F_QC"])                    # kidney
QRc <- pull(Rat.input[Rat.input$tissue == "Heart", "F_QC"]) +
  pull(Rat.input[Rat.input$tissue == "Kidney", "F_QC"]) +
  pull(Rat.input[Rat.input$tissue == "Spleen", "F_QC"])                         # rapidly perfused tissues
QRc_gluc <- pull(Rat.input[Rat.input$tissue == "Heart", "F_QC"]) +
  pull(Rat.input[Rat.input$tissue == "Spleen", "F_QC"])                         # rapidly perfused tissues
QFc <- pull(Rat.input[Rat.input$tissue == "Adipose", "F_QC"])                   # fat


### Fluid transit times (flow rates) ----
# Reference: Musther et al. (2017) DOI: 10.1016/j.xphs.2017.05.003 SI and
#            Punt et al. (2020) https://dx.doi.org/10.1021/acs.chemrestox.0c00307)

# Times
t_GER <- 0.25                                                                   # h; stomach emptying time
t_SI <- 1.63                                                                    # h; small intestine/residence transit time
t_si_duo <- 0.06*t_SI                                                           # h; duodenum transit/residence time
t_si_jej <- 0.769*t_SI                                                          # h; jejunum transit/residence time 
t_si_ile <- 0.171*t_SI                                                          # h; ileum transit/residence time 

# Emptying/transit/residence rates
k_GER = 1/t_GER                                                                 # /h; stomach 
k_duo <- 1/t_si_duo                                                             # /h; duodenum 
k_jej <- 1/t_si_jej                                                             # /h; jejunum
k_ile <- 1/t_si_ile                                                             # /h; ileum

# Glomerular filtration rate
GFR <- 0.0786                                                                   # L/h; Reference: Grandoni et al. (2019) DOI:  10.5599/admet.638 



## Tissue:plasma partition coefficients and fraction unbound ----

# Tissue partitioning is calculated using the equations developed by Rodgers and Rowland, as coded by Utsey et al. (2020) DOI: 10.1124/dmd.120.090498: https://github.com/metrumresearchgroup/PBPK_PC/blob/master/script/CalcKp_P%26T.R
# Assumption that all substances are diprotic acids

PhysChem <- read_delim("Input/PhysChem.csv", 
                       delim = ";", escape_double = FALSE, trim_ws = TRUE)

calc_Kp_fus <- function(logP, pKa1, pka2, BPR, Rat.input, F_Vo.dat){
  
  dat_all <- Rat.input %>% filter(!tissue %in% c("RBCs", "Adipose", "Plasma"))  # df for all tissues except for adipose, RBCs, and plasma
  dat_ad <- Rat.input %>% filter(tissue == "Adipose")                           # df for adipose
  dat_rbc <- Rat.input %>% filter(tissue == "RBCs")                             # df for RBCs
  dat_plas <- Rat.input %>% filter(tissue == "Plasma")                          # df for plasma
  
  pH_IW <- 7                                                                    # pH of intracellular tissue water
  pH_P <- 7.4                                                                   # pH of plasma
  pH_RBC <- 7.22                                                                # pH of blood cells
  
  P <- 10^(logP)                                                                # octanol:water partition coeff
  logP_OW <- 1.115*logP - 1.35                                                  # oil:water partition coeff
  P_OW <- 10^(logP_OW) 
  
  HCT <- 0.45                                                                   # hematocrit
  
  # Calculate fraction unbound in plasma and tissues
  # Reference: Lobell and Sivarajah 2003
  # Equation slightly different depending on the ionisation status in plasma, therefore between glycuronides and aglycones
  
  if(pKa1 >= 7){
    fup <- 1/(1 + 10^((0.4485*logP) - 0.4782))                                  # fraction unbound in plasma
  }else {
    fup <- 1/(1 + 10^((0.3649*logP) + 0.41678))                                 # fraction unbound in plasma
  }

  fut <- 1/(1 + ((1 - fup)/fup)*0.5)                                            # fraction unbound in tissues
   
  
  # Calculate partitioning and ionisation status
  Kpu_bc <- (HCT - 1 + BPR) / (HCT*fup)                                         # partitioning to blood cells                     
  X <- 10^(pH_IW - pKa1) + 10^(2*pH_IW - pKa1 - pKa2)                           # ionisation status in cells    
  Y <- 10^(pH_P - pKa1) + 10^(2*pH_P - pKa1 - pKa2)                             # ionisation status in plasma 
  Z <- 1 
  Ka_PR <- (                                                
    1/fup - 1 - ((P*dat_plas$f_n_l + (0.3*P + 0.7)*dat_plas$f_n_pl) / (1+Y))
  )                                                                             # albumin association constant
  
  # Calculate tissue:plasma partition coefficients for non-lipid tissues
  Kp_all <- (
    dat_all$f_ew + ((1 + X)/(1 + Y))*dat_all$f_iw + 
      ((P*dat_all$f_n_l + (0.3*P + 0.7)*dat_all$f_n_pl)) / (1 + Y) + 
      (Ka_PR*dat_all$AR*X)/(1 + Y)
  )*fup  
  
  # Calculate tissue:plasma partition coefficients for the adipose tissue
  Kp_ad <- (
    dat_ad$f_ew + ((1 + X)/(1 + Y))*dat_ad$f_iw + 
      ((P_OW*dat_ad$f_n_l + (0.3*P_OW + 0.7)*dat_ad$f_n_pl)) / (1 + Y) + 
      (Ka_PR*dat_ad$AR*X)/(1 + Y)
  )*fup  
  
  nms_all <- dat_all$tissue %>% substr(1,2) %>% tolower()
  nms_all <- paste("Kp", nms_all, sep="")
  nms <- c("Kpad",nms_all)
  Kp <- as.list(c(Kp_ad,Kp_all)) 
  names(Kp) <- nms  
  
  PSI <- Kp$Kpgu                                                                # Intestine:Blood partition coefficient
  PL <- Kp$Kpli                                                                 # Liver:Blood partition coefficient
  PF <- Kp$Kpad                                                                 # Adipose:Blood partition coefficient
  PK <- Kp$Kpki # only for glucuronides
  PS <- (
    F_Vo.dat$FVBo*Kp$Kpbo + F_Vo.dat$FVBr*Kp$Kpbr +
      F_Vo.dat$FVMu*Kp$Kpmu + F_Vo.dat$FVSk*Kp$Kpsk
  ) /
    (F_Vo.dat$FVBo + F_Vo.dat$FVBr + 
       F_Vo.dat$FVMu + F_Vo.dat$FVSk)                                           # Slowly perfused = partition coefficient of each of these organs * fractional volume of the respective organ / sum of the fractional volume of all these organs
  PR <- (
    F_Vo.dat$FVHe*Kp$Kphe + F_Vo.dat$FVKi*Kp$Kpki + 
      F_Vo.dat$FVLu*Kp$Kplu + F_Vo.dat$FVSp*Kp$Kpsp
  ) /
    (F_Vo.dat$FVHe + F_Vo.dat$FVKi + 
       F_Vo.dat$FVLu + F_Vo.dat$FVSp)                                           # Rapidly perfused = partition coefficient of each of these organs * fractional volume of the respective organ / sum of the fractional volume of all these organs
  PR_gluc <- ( # only for glucuronides
    F_Vo.dat$FVHe*Kp$Kphe +  
      F_Vo.dat$FVLu*Kp$Kplu + F_Vo.dat$FVSp*Kp$Kpsp
  ) /
    (F_Vo.dat$FVHe + F_Vo.dat$FVKi + 
       F_Vo.dat$FVLu + F_Vo.dat$FVSp)                                           # Rapidly perfused = partition coefficient of each of these organs * fractional volume of the respective organ / sum of the fractional volume of all these organs
  
  Kp_fus_BPR <- as.list(c(PSI = PSI, 
                          PL = PL, 
                          PF = PF, 
                          PK = PK,
                          PS = PS, 
                          PR = PR,
                          PR_gluc= PR_gluc,
                          fup = fup,
                          fut = fut,
                          BPR = BPR))
  
  return(Kp_fus_BPR)
  
}


mylist <- list() 


for (i in 1:nrow(PhysChem)) {
  
  logP <- PhysChem$logP[i]
  pKa1 <- PhysChem$pKa.1[i]
  pKa2 <- PhysChem$pKa.2[i] 
  BPR  <- PhysChem$BPR[i] 

  Kp_fus_BPR <- calc_Kp_fus(logP = logP, 
                            pKa1 = pKa1, 
                            pka2 = pKa2, 
                            BPR = BPR, 
                            Rat.input = Rat.input,
                            F_Vo.dat = F_Vo.dat) 
  
  mylist[[i]] <- Kp_fus_BPR
  
}


Kp_fus_BPR_df <- do.call("rbind",mylist) 
rownames(Kp_fus_BPR_df) <- c("ZEN", "aZEL", "ZEN_gluc", "aZEL_gluc")
rm(mylist)
rm(Kp_fus_BPR)


## Kinetic parameters ----

### Absorption/transfer rates ---

# Apparent permeability across the intestinal membrane
# Reference: Pfeiffer et al. (2011) DOI 10.1002/mnfr.201000381
PappCaco2_ZEN <- 10.4 * 1e-6                                                    # cm/s 
PappCaco2_aZEL <- 5.4 * 1e-6                                                    # cm/s

# Absorption to (& release from) the lymph
ka_Lymph <- 0.55                                                                # /h

### Metabolism ----

#### Small intestine ----

S9SI <- 21.9                                                                    # small intestinal S9 protein yield (mg S9 protein/gram intestine); Reference: Doerksen et al. (2020) doi: 10.1007/s13318-020-00666-w sum of microsomal 3.9 mg/g and cytosolic 18mg/g 

# Hydrolysis
# Assuming the same as liver
# Reference: Malekinejad et al. (2006) doi:10.1016/j.tvjl.2005.03.004
Vmax_SI_aZEL.c <- 32*1e-6                                                       # umol/min/mg protein; from pmol-> umol
Km_SI_aZEL <- 592                                                               # uM 
Vmax_SI_bZEL.c <- 72*1e-6                                                       # umol/min/mg protein; from pmol-> umol
Km_SI_bZEL <- 21                                                                # uM 

# Glucuronidation
# Reference: Mendez-Catala et al. (2021) DOI 10.1002/mnfr.202100443
Alamethicin_correction <- 1/3                                                   # correction for the overestimation of UGD metabolism in S9, Miners et al. (2021) https://doi.org/10.1016/j.pharmthera.2020.107689
Vmax_SI_ZEN_gluc.c <-  8.617*1e-3 * Alamethicin_correction                      # umol/min/mg protein; from nmol-> umol 
Km_SI_ZEN_gluc <- 6.198                                                         # uM 
Vmax_SI_aZEL_gluc.c <- 6.963*1e-3 * Alamethicin_correction                      # umol/min/mg protein; from nmol-> umol
Km_SI_aZEL_gluc <- 7.426                                                        # uM

# Glucuronide hydrolysis by gmGUS
# Reference: Ebuzoeme et al. (2021) https://doi.org/  10.3390/pharmaceutics13071043 (wogonoside glucuronide)
S9gmGUS <- 4                                                                    # mg protein/ml feces
Vmax_Hydro <- 4.42                                                              # umol/h
Km_Hydro <- 3.08                                                                # uM

#### Liver ----

S9L <- 87                                                                       # liver S9 protein yield, (mg S9 protein/gram liver) reference: Mendez-catala et al. (2020) https://doi.org/10.1016/j.toxrep.2020.07.010

# Hydrolysis
# Reference: Malekinejad et al. (2006) doi:10.1016/j.tvjl.2005.03.004
Vmax_L_aZEL.c <- 32*1e-6                                                        # pmol->umol/min/mg protein
Km_L_aZEL <- 592                                                                # uM
Vmax_L_bZEL.c <- 72*1e-6                                                        # pmol->umol/min/mg protein   
Km_L_bZEL <- 21                                                                 # uM

# Glucuronidation
# Reference: Mendez-Catala et al. (2021) DOI 10.1002/mnfr.202100443
Alamethicin_correction <- 1/3                                                   # correction for the overestimation of UGD metabolism in S9, Miners et al. (2021) https://doi.org/10.1016/j.pharmthera.2020.107689
Vmax_L_ZEN_gluc.c <- 7.029*1e-3  /3                                               # nmol->umol/min/mg S9
Km_L_ZEN_gluc <- 6.753                                                          # uM
Vmax_L_aZEL_gluc.c <- 6.963*1e-3  /3                                            # nmol->umol/min/mg S9
Km_L_aZEL_gluc <- 7.426                                                         # uM 0.87


### Excretion ----

#### Biliary excretion ----

Bile.fasted <- 6.35e-6 * 60                                                     # L/h
# BileIPf <- approxfun(Bile_Lymph.df$Time,                                        # L/h
#                      Bile_Lymph.df$BileContentsExcretionR*1e-6*60,
#                      rule = 2)
# plot(BileIPf(timeSeq))

#### Active transport ----
# Reference: Yanfei Li et al. 2019 doi:10.1016/j.jsbmb.2019.03.027.

##### Constants ----
Vmax_MRP2.c <- 3366e-6                                                          # umol/min/mg MRP2; pmol -> umol
Km_MRP2 <- 78.9                                                                 # uM
Vmax_MRP3.c <- 986.2e-6                                                         # umol/min/mg MRP3; pmol -> umol
Km_MRP3 <- 5.7                                                                  # uM
Vmax_BCRP.c <- 2470e-6                    
Km_BCRP <- 676.4 

##### IVIVE ----
# Inside-out (active) fraction of transporter in in vitro system vesicles
f_in_out_MRP2 <- 0.42 
f_in_out_MRP3 <- 0.58    
f_in_out_BCRP <- 0.326

# Protein abundance in vesicles
pmol_mg_MRP2 <- 133.9                                                           # pmol transporter/mg total protein
pmol_mg_MRP3 <- 239.4                                                           # pmol transporter/mg total protein  
pmol_mg_BCRP <- 11.4                                                            # pmol transporter/mg total protein

# Protein abundance in vivo (tissues)
MW_MRP2 <- 190e3                                                                # kDa -> g/mol; Reference: https://www.sciencedirect.com/science/article/abs/pii/S0308814619305898
pmolMP_mgP_IN_MRP2 <- 4.76e3/MW_MRP2                                            # pmol MRP2/mg protein in intestine; ng->pg, pg -> pmol; Reference: Gavins et al. 2023 https://doi.org/10.1016/j.biopha.2023.114329 
pmolMP_mgP_LI_MRP2 <- 2.08                                                      # pmol MRP2/mg protein in liver; Reference: Fallon 2016 DOI 10.1007/s11095-016-1966-5
pmolMP_mgP_KI_MRP2 <- 1.4                                                       # pmol MRP2/mg protein in kidney; Reference: Fallon 2016 DOI 10.1007/s11095-016-1966-5

MW_MRP3 <- 169.4e3                                                              # kDa -> g/mol; Reference: https://pmc.ncbi.nlm.nih.gov/articles/PMC22016/
pmolMP_mgP_LI_MRP3 <- 0.32                                                      # pmol MRP3/mg protein in liver; Reference: Fallon 2016 DOI 10.1007/s11095-016-1966-5

MW_BCRP <- 75e3                                                                 # kDa -> g/mol; Reference: Rosenberg et al. 2010 doi: 10.1016/j.str.2010.01.017
pmolMP_mgP_IN_BCRP <- 12.51e3/MW_BCRP                                           # pmol BCRP/mg protein in intestine; ng->pg, pg -> pmol; Reference: Gavins et al. 2023 https://doi.org/10.1016/j.biopha.2023.114329 
pmolMP_mgP_LI_BCRP <- 0.17                                                      # pmol BCRP/mg protein in liver; Reference: Fallon 2016 DOI 10.1007/s11095-016-1966-5


# Total protein in tissues
mgP_mgIN <- 0.133                                                               # mg protein/mg intestine
mgP_mgLI <- 0.21                                                                # mg protein/mg liver
mgP_mgKI <- 0.18                                                                # mg protein/mg kidney

# Relative expression factor
REF_IN_MRP2 <- pmolMP_mgP_IN_MRP2 / pmol_mg_MRP2
REF_LI_MRP2 <- pmolMP_mgP_LI_MRP2 / pmol_mg_MRP2
REF_KI_MRP2 <- pmolMP_mgP_KI_MRP2 / pmol_mg_MRP2

REF_LI_MRP3 <- pmolMP_mgP_LI_MRP3 / pmol_mg_MRP3

REF_IN_BCRP <- pmolMP_mgP_IN_BCRP / pmol_mg_BCRP
REF_LI_BCRP <- pmolMP_mgP_LI_BCRP / pmol_mg_BCRP

# Scaling factor
SF_IN_MRP2 <- mgP_mgIN * REF_IN_MRP2 / f_in_out_MRP2                            # mg MRP2 / mg small intestine
SF_LI_MRP2 <- mgP_mgLI * REF_LI_MRP2 / f_in_out_MRP2 *10                        # mg MRP2 / mg liver
SF_KI_MRP2 <- mgP_mgKI * REF_KI_MRP2 / f_in_out_MRP2 *10                        # mg MRP2 / mg liver

SF_LI_MRP3 <- mgP_mgLI * REF_LI_MRP3 / f_in_out_MRP3                            # mg MRP3 / mg liver 

SF_IN_BCRP <- mgP_mgIN * REF_IN_BCRP / f_in_out_BCRP                            # mg BCRP / mg liver 
SF_LI_BCRP <- mgP_mgLI * REF_LI_BCRP / f_in_out_BCRP *10                        # mg BCRP / mg liver 

## Collect all parameters ----
Kp_fus_BPR_df <- as.data.frame(Kp_fus_BPR_df) %>%
  mutate(Substance = rownames(Kp_fus_BPR_df)) %>%
  pivot_longer(!Substance, names_to = "Parameter", values_to = "Value") %>%
  unite("Parameter", Substance:Parameter) 

Other_params <- tibble(
  BW = BW,  
  FV_SI = FV_SI,
  FV_Duo = FV_Duo,
  FV_Jej = FV_Jej,
  FV_Ile = FV_Ile,
  FV_SI_IC = FV_SI_IC,
  FV_SI_EC = FV_SI_EC,
  FV_SIL = FV_SIL,
  FV_duol = FV_duol,
  FV_jejl = FV_jejl,
  FV_ilel = FV_ilel,
  Vfe_IleL = Vfe_IleL,
  FV_L = FV_L,
  FV_L_IC = FV_L_IC,
  FV_L_EC = FV_L_EC,
  FV_R = FV_R,
  FV_R_gluc = FV_R_gluc, 
  FV_F = FV_F,
  FV_B = FV_B,
  FV_K = FV_K,
  VBile = VBile,
  SAduo = SAduo,
  SAjej = SAjej,
  SAile = SAile,
  QC = QC,
  QSIc = QSIc,
  QDuoc = QDuoc,
  QJejc = QJejc,
  QIlec = QIlec,
  QLc = QLc,
  QKc = QKc,
  QRc = QRc,
  QRc_gluc = QRc_gluc,
  QFc = QFc,
  k_GER = k_GER,
  k_duo = k_duo,
  k_jej = k_jej,
  k_ile = k_ile,
  GFR = GFR,
  PappCaco2_ZEN = PappCaco2_ZEN,
  PappCaco2_aZEL = PappCaco2_aZEL,
  ka_Lymph = ka_Lymph,
  S9SI =  S9SI,
  Vmax_SI_aZEL.c = Vmax_SI_aZEL.c,
  Km_SI_aZEL = Km_SI_aZEL,
  Vmax_SI_bZEL.c = Vmax_SI_bZEL.c,
  Km_SI_bZEL = Km_SI_bZEL,
  Vmax_SI_ZEN_gluc.c =  Vmax_SI_ZEN_gluc.c,
  Km_SI_ZEN_gluc = Km_SI_ZEN_gluc,
  Vmax_SI_aZEL_gluc.c = Vmax_SI_aZEL_gluc.c,
  Km_SI_aZEL_gluc = Km_SI_aZEL_gluc,
  S9gmGUS = S9gmGUS,
  Vmax_Hydro = Vmax_Hydro,
  Km_Hydro = Km_Hydro,
  S9L = S9L,
  Vmax_L_aZEL.c = Vmax_L_aZEL.c,
  Km_L_aZEL = Km_L_aZEL,
  Vmax_L_bZEL.c = Vmax_L_bZEL.c,
  Km_L_bZEL = Km_L_bZEL,
  Vmax_L_ZEN_gluc.c = Vmax_L_ZEN_gluc.c,
  Km_L_ZEN_gluc = Km_L_ZEN_gluc,
  Vmax_L_aZEL_gluc.c = Vmax_L_aZEL_gluc.c,
  Km_L_aZEL_gluc = Km_L_aZEL_gluc,
  Bile.fasted = Bile.fasted,
  Vmax_MRP2.c = Vmax_MRP2.c,
  Km_MRP2 = Km_MRP2,
  Vmax_MRP3.c = Vmax_MRP3.c,
  Km_MRP3 = Km_MRP3,
  Vmax_BCRP.c = Vmax_BCRP.c,
  Km_BCRP = Km_BCRP,
  SF_IN_MRP2 = SF_IN_MRP2,
  SF_LI_MRP2 = SF_LI_MRP2,
  SF_LI_MRP3 = SF_LI_MRP3,
  SF_KI_MRP2 = SF_KI_MRP2,
  SF_IN_BCRP = SF_IN_BCRP,
  SF_LI_BCRP = SF_LI_BCRP
) %>% 
  pivot_longer(BW:SF_LI_BCRP,names_to = "Parameter", values_to = "Value")

Parameters_df <- rbind(Kp_fus_BPR_df, Other_params) %>% 
  mutate(Value = as.numeric(Value))

Parameters_df <- Parameters_df %>% filter(!grepl("ZEN_PK", Parameter),
                                          !grepl("ZEN_PR_gluc", Parameter),
                                          !grepl("aZEL_PK", Parameter),
                                          !grepl("aZEL_PR_gluc", Parameter),
                                          !grepl("\\bZEN_gluc_PR\\b", Parameter),
                                          !grepl("\\baZEL_gluc_PR\\b", Parameter))

write.csv(Parameters_df, here("Input", "Parameters_rat.csv"), row.names = FALSE)

# PBK Model ----

R.PBK.model <- function(t, state, parms){
  with(as.list(c(state, parms)), {
    
    ## Model Parameters ----
    
    # Physiological
    
    ## Volumes
    BW <- BW  
    
    VSIL <- FV_SIL*BW 
    VDuoL <- VSIL*FV_duol 
    VJejL <- VSIL*FV_jejl 
    VIleL <- VSIL*FV_ilel 

    VSI <- FV_SI*BW 
    VDuo <- VSI*FV_Duo
    VJej <- VSI*FV_Jej
    VJej_IC <- FV_SI_IC*VJej
    VJej_EC <- FV_SI_EC*VJej
    VIle <- VSI*FV_Ile
    
    VL <- FV_L*BW  
    VL_IC <- FV_L_IC*VL
    VL_EC <- FV_L_EC*VL
    
    VR <- FV_R*BW 
    VR_gluc <- FV_R_gluc*BW
    VF <- FV_F*BW  
    VB <- FV_B*BW   
    VS <- (0.96*BW) - (VSI + VL + VR + VF + VB)                                 # fraction of slowly perfused tissue (Bone/Skin/Adipose/Muscle/(brain and testes) (slowly equilibrating tissues with large time constants))
    VS_gluc <- (0.96*BW) - (VSI + VL + VR_gluc + VF + VB)                       # fraction of slowly perfused tissue (Bone/Skin/Adipose/Muscle/(brain and testes) (slowly equilibrating tissues with large time constants))
    
    
    # For MRP3 kidney excretion use only
    VK <- FV_K*BW

    ## Blood flow rates 
    QC <- QC  
    
    QSI <- QSIc*QC  
    QDuo <- QSI*QDuoc
    QJej <- QSI*QJejc
    QIle <- QSI*QIlec
    
    QL <- QLc*QC  
    QK <- QKc*QC
    QR <- QRc*QC  
    QR_gluc <- QRc_gluc*QC
    QF <- QFc*QC
    QS <- QC - (QSI + QL + QR + QF) 
    QS_gluc <- QC - (QSI + QL + QR_gluc + QF) 
    
    
    # Kinetic
    
    ## Intestinal absorption 
    
    ### Passive permeability
    CLduo_ZEN <- PappCaco2_ZEN * SAduo * 3600 * 1e-3                            # L/h; from cm/s to cm3/h to 
    CLjej_ZEN <- PappCaco2_ZEN * SAjej * 3600 * 1e-3                            # L/h; from cm/s to cm3/h to 
    CLile_ZEN <- PappCaco2_ZEN * SAile * 3600 * 1e-3                            # L/h; from cm/s to cm3/h to 
    
    CLile_aZEL <- PappCaco2_aZEL * SAile * 3600 * 1e-3                          # L/h; from cm/min to cm3/h to 
    
    
    ### Lymphatic absorption
    LymphIPf    <- splinefun(c(0,
                               1.5,
                               2.4,
                               3.5,
                               4.5,
                               5.5,
                               7.5,
                               9.5,
                               15,
                               20,
                               24), ka_Lymph*c(0.0084,
                                               0.0503,
                                               0.3659,
                                               1,
                                               0.5810,
                                               0.3715,
                                               0.1704,
                                               0.1369,
                                               0.0100,
                                               0.0090,
                                               0.0084),
                             method = "natural")
    ka_Lymph <- LymphIPf(t)
    
    ## Metabolism
    
    ### Small intestinal tissue
    Vmax_SI_aZEL <- Vmax_SI_aZEL.c * 60 * S9SI * VJej * 1e3                     # umol/h; from umol/ min->h /mg protein -> whole intestine (kg)
    Vmax_SI_bZEL <- Vmax_SI_bZEL.c * 60 * S9SI * VJej * 1e3                     # umol/h; from umol/ min->h /mg protein -> whole intestine (kg)
    Vmax_SI_ZEN_gluc <- Vmax_SI_ZEN_gluc.c * 60 * 
      S9SI * VJej * 1e3                                                         # umol/h; from umol/ min->h /mg protein -> whole intestine (kg)
    Vmax_SI_aZEL_gluc <- Vmax_SI_aZEL_gluc.c * 60 * 
      S9SI * VJej  * 1e3                                                        # umol/h; from umol/ min->h /mg protein -> whole intestine (kg)
    
    Vmax_Hydro_IleL <- Vmax_Hydro * 60 * S9gmGUS * Vfe_IleL                     # umol/h; from umol/ min->h /mg protein -> whole intestinal micriobiome (kg)

    
    ### Liver
    Vmax_L_aZEL <- Vmax_L_aZEL.c * 60 * S9L * VL * 1e3                          # umol/h; from umol/ min->h /mg protein -> whole liver (kg)
    Vmax_L_bZEL <- Vmax_L_bZEL.c * 60 * S9L * VL * 1e3                          # umol/h; from umol/ min->h /mg protein -> whole liver (kg)
    Vmax_L_ZEN_gluc <- Vmax_L_ZEN_gluc.c * 60 * 
      S9L * VL * 1e3                                                            # umol/h; from umol/ min->h /mg protein -> whole liver (kg)
    Vmax_L_aZEL_gluc <- Vmax_L_aZEL_gluc.c * 60 *
      S9L * VL * 1e3                                                            # umol/h; from umol/ min->h /mg protein -> whole liver (kg)
    
    
    ## Active transport
    
    Vmax_MRP2_LI <- Vmax_MRP2.c * 60 * SF_LI_MRP2 * VL * 1e6                    # umol/min/mg liver; scaling: umol/min->h/mg MRP2 * mg MRP2 / mg liver -> whole liver (kg) 
    Vmax_MRP2_IN <- Vmax_MRP2.c * 60 * SF_IN_MRP2 * VSI * 1e6                   # umol/min/mg intestine; scaling: umol/min->h/mg MRP2 * mg MRP2 / mg liver -> whole intestine (kg)  
    
    Vmax_MRP3_LI <- Vmax_MRP3.c * 60 * SF_LI_MRP3 * VL * 1e6                    # umol/min/mg liver; scaling: umol/min->h/mg MRP3 * mg MRP3 / mg liver -> whole liver (kg)
    Vmax_MRP2_KI <- Vmax_MRP2.c * 60 * SF_KI_MRP2 * VK * 1e6                    # umol/min/mg kidney; scaling: umol/min->h/mg MRP3 * mg MRP3 / mg liver -> whole VK_IC (kg)
    
    Vmax_BCRP_LI <- Vmax_BCRP.c * 60 * SF_LI_BCRP * VL * 1e6                    # umol/min/mg liver; scaling: umol/min->h/mg BCRP * mg BCRP / mg liver -> whole liver (kg)
    Vmax_BCRP_IN <- Vmax_BCRP.c * 60 * SF_IN_BCRP * VSI * 1e6                   # umol/min/mg liver; scaling: umol/min->h/mg BCRP * mg BCRP / mg liver -> whole liver (kg)
    
    ## Biliary clearance
    
    # CL_Bile <- ifelse(t > 3 && t < 6, Bile.fasted, 1) #
    BileIPf    <- splinefun(c(0,
                              2,
                              4,
                              5,
                              7,
                              10,
                              13,
                              16,
                              20,
                              24), Bile.fasted*c(1.2368,
                                                 1.0394,
                                                 0.8743,
                                                 0.8793,
                                                 0.9297,
                                                 1.2685,
                                                 1.7844,
                                                 2.0877,
                                                 1.7761,
                                                 1.5435),
                            method = "natural")
    CL_Bile <- BileIPf(Times)
    
    ## Tissue Concentrations ----
    
    ### Zearalenone ----
    CDuoL_ZEN = ADuoL_ZEN/VDuoL
    CJejL_ZEN = AJejL_ZEN/VJejL
    CIleL_ZEN = AIleL_ZEN/VIleL
    CDuo_ZEN = ADuo_ZEN/VDuo
    CJej_ZEN = AJej_ZEN/VJej
    CIle_ZEN = AIle_ZEN/VIle
    CL_ZEN = AL_ZEN/VL  
    CR_ZEN = AR_ZEN/VR   
    CS_ZEN = AS_ZEN/VS
    CF_ZEN = AF_ZEN/VF
    CB_ZEN = AB_ZEN/VB
    
    CVDuo_ZEN <- CDuo_ZEN/ZEN_PSI * ZEN_BPR
    CVJej_ZEN <- CJej_ZEN/ZEN_PSI * ZEN_BPR
    CVIle_ZEN <- CIle_ZEN/ZEN_PSI * ZEN_BPR
    CVL_ZEN <- CL_ZEN/ZEN_PL * ZEN_BPR
    CVF_ZEN <- CF_ZEN/ZEN_PF * ZEN_BPR
    CVR_ZEN <- CR_ZEN/ZEN_PR * ZEN_BPR
    CVS_ZEN <- CS_ZEN/ZEN_PS * ZEN_BPR
    
    ### Zearalenone-glucuronide ----
    CDuoL_ZEN_gluc = ADuoL_ZEN_gluc/VDuoL
    CJejL_ZEN_gluc <- AJejL_ZEN_gluc/VJejL
    CIleL_ZEN_gluc <- AIleL_ZEN_gluc/VIleL 
    CJej_IC_ZEN_gluc <- AJej_IC_ZEN_gluc/VJej_IC
    CJej_EC_ZEN_gluc <- AJej_EC_ZEN_gluc/VJej_EC
    CL_IC_ZEN_gluc <- AL_IC_ZEN_gluc/VL_IC
    CL_EC_ZEN_gluc <- AL_EC_ZEN_gluc/VL_EC
    CBile_ZEN_gluc <- ABile_ZEN_gluc/VBile
    CK_ZEN_gluc <- AK_ZEN_gluc/VK
    CR_ZEN_gluc <- AR_ZEN_gluc/VR_gluc
    CS_ZEN_gluc <- AS_ZEN_gluc/VS_gluc
    CF_ZEN_gluc <- AF_ZEN_gluc/VF
    CB_ZEN_gluc <- AB_ZEN_gluc/VB
    
    CVJej_EC_ZEN_gluc <- CJej_EC_ZEN_gluc/ZEN_gluc_PSI * ZEN_gluc_BPR
    CVL_EC_ZEN_gluc <- CL_EC_ZEN_gluc/ZEN_gluc_PL * ZEN_gluc_BPR
    CVK_ZEN_gluc <- CK_ZEN_gluc/ZEN_gluc_PK * ZEN_gluc_BPR
    CVF_ZEN_gluc <- CF_ZEN_gluc/ZEN_gluc_PF * ZEN_gluc_BPR
    CVR_ZEN_gluc <- CR_ZEN_gluc/ZEN_gluc_PR_gluc * ZEN_gluc_BPR
    CVS_ZEN_gluc <- CS_ZEN_gluc/ZEN_gluc_PS * ZEN_gluc_BPR
    
    ### α-Zearalenol ----
    CIleL_aZEL = AIleL_aZEL/VIleL
    CJej_aZEL = AJej_aZEL/VJej
    CIle_aZEL = AIle_aZEL/VIle
    CL_aZEL = AL_aZEL/VL  
    CR_aZEL = AR_aZEL/VR   
    CS_aZEL = AS_aZEL/VS
    CF_aZEL = AF_aZEL/VF
    CB_aZEL = AB_aZEL/VB
    
    CVIle_aZEL <- CIle_aZEL/aZEL_PSI * aZEL_BPR
    CVJej_aZEL <- CJej_aZEL/aZEL_PSI * aZEL_BPR
    CVL_aZEL <- CL_aZEL/aZEL_PL * aZEL_BPR
    CVF_aZEL <- CF_aZEL/aZEL_PF * aZEL_BPR
    CVR_aZEL <- CR_aZEL/aZEL_PR * aZEL_BPR
    CVS_aZEL <- CS_aZEL/aZEL_PS * aZEL_BPR    
    
    ### a-Zearalenol-glucuronide ----
    
    CDuoL_aZEL_gluc = ADuoL_aZEL_gluc/VDuoL
    CJejL_aZEL_gluc <- AJejL_aZEL_gluc/VJejL
    CIleL_aZEL_gluc <- AIleL_aZEL_gluc/VIleL 
    CJej_IC_aZEL_gluc <- AJej_IC_aZEL_gluc/VJej_IC
    CJej_EC_aZEL_gluc <- AJej_EC_aZEL_gluc/VJej_EC
    CL_IC_aZEL_gluc <- AL_IC_aZEL_gluc/VL_IC
    CL_EC_aZEL_gluc <- AL_EC_aZEL_gluc/VL_EC
    CK_aZEL_gluc <- AK_aZEL_gluc/VK
    CBile_aZEL_gluc <- ABile_aZEL_gluc/VBile
    CR_aZEL_gluc <- AR_aZEL_gluc/VR_gluc
    CS_aZEL_gluc <- AS_aZEL_gluc/VS_gluc
    CF_aZEL_gluc <- AF_aZEL_gluc/VF
    CB_aZEL_gluc <- AB_aZEL_gluc/VB
    
    CVJej_EC_aZEL_gluc <- CJej_EC_aZEL_gluc/aZEL_gluc_PSI * aZEL_gluc_BPR
    CVL_EC_aZEL_gluc <- CL_EC_aZEL_gluc/aZEL_gluc_PL * aZEL_gluc_BPR
    CVF_aZEL_gluc <- CF_aZEL_gluc/aZEL_gluc_PF * aZEL_gluc_BPR
    CVK_aZEL_gluc <- CK_aZEL_gluc/aZEL_gluc_PK * aZEL_gluc_BPR
    CVR_aZEL_gluc <- CR_aZEL_gluc/aZEL_gluc_PR_gluc * aZEL_gluc_BPR
    CVS_aZEL_gluc <- CS_aZEL_gluc/aZEL_gluc_PS * aZEL_gluc_BPR
   
    
    # ## IV Perf ----
    # 
    # if(t<=6){DoseOn=1} else{DoseOn=0}
    # IV = IVDOSEumol*DoseOn                                                    
    
    ## Differential Equations ----
    
    dAOD_ZEN <- - k_GER*AOD_ZEN                                        
    
    
    ### Zearalenone ----
    
    
    # Intestinal lumen 
    dADuoL_ZEN <- + k_GER*AOD_ZEN - k_duo*ADuoL_ZEN + 
      - CLduo_ZEN*CDuoL_ZEN - ka_Lymph*ADuoL_ZEN
    
    dAJejL_ZEN <- + k_duo*ADuoL_ZEN - k_jej*AJejL_ZEN + 
      - CLjej_ZEN*CJejL_ZEN - ka_Lymph*AJejL_ZEN
    
    dAIleL_ZEN <- + k_jej*AJejL_ZEN - k_ile*AIleL_ZEN +  
      - CLile_ZEN*CIleL_ZEN - ka_Lymph*AIleL_ZEN +  
      + ((Vmax_Hydro_IleL*CIleL_ZEN_gluc) / 
           (Km_Hydro + (CIleL_ZEN_gluc)))

    # Small intestinal tissue 
    dADuo_ZEN <- CLduo_ZEN*CDuoL_ZEN + QDuo*(CB_ZEN-CVDuo_ZEN) 
    
    dAJej_ZEN <- CLjej_ZEN*CJejL_ZEN + QJej*(CB_ZEN-CVJej_ZEN) + 
      - ((Vmax_SI_ZEN_gluc*CJej_ZEN*ZEN_fut) / 
           (Km_SI_ZEN_gluc + (CJej_ZEN*ZEN_fut))) +
      - ((Vmax_SI_aZEL*CJej_ZEN*ZEN_fut) / 
           (Km_SI_aZEL + (CJej_ZEN*ZEN_fut))) +
      - ((Vmax_SI_bZEL*CJej_ZEN*ZEN_fut) / 
           (Km_SI_bZEL + (CJej_ZEN*ZEN_fut)))
    
    dAIle_ZEN <- CLile_ZEN*CIleL_ZEN + QIle*(CB_ZEN-CVIle_ZEN) 
    
    # Liver 
    dAL_ZEN <- QDuo*CVDuo_ZEN + QJej*CVJej_ZEN + QIle*CVIle_ZEN + 
      QL*CB_ZEN - (QSI+QL)*CVL_ZEN + 
      - ((Vmax_L_aZEL*CL_ZEN*ZEN_fut) / 
           (Km_L_aZEL + (CL_ZEN*ZEN_fut))) +
      - ((Vmax_L_bZEL*CL_ZEN*ZEN_fut) / 
           (Km_L_bZEL + (CL_ZEN*ZEN_fut))) +
      - ((Vmax_L_ZEN_gluc*CL_ZEN*ZEN_fut) / 
           (Km_L_ZEN_gluc + (CL_ZEN*ZEN_fut)))
    
    # Rapidly Perfused 
    dAR_ZEN <- QR*(CB_ZEN-CVR_ZEN) 
    
    # Slowly Perfused 
    dAS_ZEN <- QS*(CB_ZEN-CVS_ZEN)
    
    # Fat 
    dAF_ZEN <- QF*(CB_ZEN-CVF_ZEN) 
    
    # Blood
    dAB_ZEN <-  - (QSI + QL + QR + QS + QF)*CB_ZEN + #IV + 
      + (QSI + QL)*CVL_ZEN + 
      + QF*CVF_ZEN + QR*CVR_ZEN + QS*CVS_ZEN + 
      - GFR*CB_ZEN*ZEN_fup +
      + ka_Lymph*(ADuoL_ZEN + AJejL_ZEN + AIleL_ZEN)
    
    # Urine
    dAUR_ZEN <- GFR*CB_ZEN*ZEN_fup 
    
    # Faeces 
    dAFE_ZEN <- k_ile*AIleL_ZEN 
    
    # β-Zearalenol
    dA_bZEL <- ((Vmax_SI_bZEL*CJej_ZEN*ZEN_fut) / 
                  (Km_SI_bZEL + (CJej_ZEN*ZEN_fut))) + 
      ((Vmax_L_bZEL*CL_ZEN*ZEN_fut) / 
         (Km_L_bZEL + (CL_ZEN*ZEN_fut))) 
   
    
    ### Zearalenone glucuronide ----
    
    # Intestinal lumen 
    dADuoL_ZEN_gluc <- CL_Bile*CBile_ZEN_gluc - k_duo*ADuoL_ZEN_gluc 
    
    dAJejL_ZEN_gluc <- + k_duo*ADuoL_ZEN_gluc - k_jej*AJejL_ZEN_gluc +
      + ((Vmax_MRP2_IN*CJej_IC_ZEN_gluc*ZEN_gluc_fut)/
           (Km_MRP2 + (CJej_IC_ZEN_gluc*ZEN_gluc_fut))) +
      + ((Vmax_BCRP_IN*CJej_IC_ZEN_gluc*ZEN_gluc_fut)/
           (Km_BCRP + (CJej_IC_ZEN_gluc*ZEN_gluc_fut))) 
    
    dAIleL_ZEN_gluc <- k_jej*AJejL_ZEN_gluc - k_ile*AIleL_ZEN_gluc +
      - ((Vmax_Hydro_IleL*CIleL_ZEN_gluc) / 
           (Km_Hydro + (CIleL_ZEN_gluc)))
    
    # Small Intestinal tissue
    dAJej_IC_ZEN_gluc <- ((Vmax_SI_ZEN_gluc*CJej_ZEN*ZEN_fut) / 
                            (Km_SI_ZEN_gluc + (CJej_ZEN*ZEN_fut))) +
      - ((Vmax_MRP2_IN*CJej_IC_ZEN_gluc*ZEN_gluc_fut) / 
           (Km_MRP2 + (CJej_IC_ZEN_gluc*ZEN_gluc_fut))) +
      - ((Vmax_BCRP_IN*CJej_IC_ZEN_gluc*ZEN_gluc_fut)/
           (Km_BCRP + (CJej_IC_ZEN_gluc*ZEN_gluc_fut)))
    
    dAJej_EC_ZEN_gluc <- QJej*(CB_ZEN_gluc-CVJej_EC_ZEN_gluc) 
    
    # Liver
    dAL_IC_ZEN_gluc <- ((Vmax_L_ZEN_gluc*CL_ZEN*ZEN_fut) / 
                          ((Km_L_ZEN_gluc) + (CL_ZEN*ZEN_fut))) +
      - ((Vmax_MRP2_LI*CL_IC_ZEN_gluc*ZEN_gluc_fut) / 
           (Km_MRP2 + (CL_IC_ZEN_gluc*ZEN_gluc_fut))) + 
      - ((Vmax_MRP3_LI*CL_IC_ZEN_gluc*ZEN_gluc_fut) / 
           (Km_MRP3 + (CL_IC_ZEN_gluc*ZEN_gluc_fut))) +
      - ((Vmax_BCRP_LI*CL_IC_ZEN_gluc*ZEN_gluc_fut) / 
           (Km_BCRP + (CL_IC_ZEN_gluc*ZEN_gluc_fut))) 
    
    dAL_EC_ZEN_gluc <- + QJej*CVJej_EC_ZEN_gluc + 
      + (QL+QDuo+QIle)*CB_ZEN_gluc - (QSI + QL)*CVL_EC_ZEN_gluc +
      + ((Vmax_MRP3_LI*CL_IC_ZEN_gluc*ZEN_gluc_fut) / 
           (Km_MRP3 + (CL_IC_ZEN_gluc*ZEN_gluc_fut)))
    
    # Bile
    dABile_ZEN_gluc <- ((Vmax_MRP2_LI*CL_IC_ZEN_gluc*ZEN_gluc_fut) / 
                          (Km_MRP2 + (CL_IC_ZEN_gluc*ZEN_gluc_fut))) +
      + ((Vmax_BCRP_LI*CL_IC_ZEN_gluc*ZEN_gluc_fut) / 
         (Km_BCRP + (CL_IC_ZEN_gluc*ZEN_gluc_fut))) +
      - CL_Bile*CBile_ZEN_gluc 
    
    # Rapidly Perfused 
    dAR_ZEN_gluc <- QR_gluc*(CB_ZEN_gluc-CVR_ZEN_gluc)
    
    # Slowly Perfused 
    dAS_ZEN_gluc <- QS_gluc*(CB_ZEN_gluc-CVS_ZEN_gluc)
    
    # Fat
    dAF_ZEN_gluc <- QF*(CB_ZEN_gluc-CVF_ZEN_gluc) 
    
    # Blood
    dAB_ZEN_gluc <- - (QSI + QL + QK + QR_gluc + QS_gluc + QF)*CB_ZEN_gluc +   
      + (QSI + QL)*CVL_EC_ZEN_gluc + 
      + QF*CVF_ZEN_gluc + QK*CVK_ZEN_gluc + QR_gluc*CVR_ZEN_gluc + QS_gluc*CVS_ZEN_gluc +
      - GFR*CB_ZEN_gluc*ZEN_gluc_fup 
    
    dAK_ZEN_gluc <- QK*(CB_ZEN_gluc-CVK_ZEN_gluc) +
      + GFR*CB_ZEN_gluc*ZEN_gluc_fup +
      - ((Vmax_MRP2_KI*CK_ZEN_gluc*ZEN_gluc_fut) / 
         (Km_MRP2 + (CK_ZEN_gluc*ZEN_gluc_fut)))
    
    # Urine
    dAUR_ZEN_gluc <- ((Vmax_MRP2_KI*CK_ZEN_gluc*ZEN_gluc_fut) / 
         (Km_MRP2 + (CK_ZEN_gluc*ZEN_gluc_fut)))
    
    # Feces
    dAFE_ZEN_gluc <- k_ile*AIleL_ZEN_gluc 
    
    ### α-Zearalenol ----
    
    # Intestinal lumen
    dAIleL_aZEL <- - k_ile*AIleL_aZEL +  
      - CLile_aZEL*CIleL_aZEL - ka_Lymph*AIleL_aZEL +  
      + ((Vmax_Hydro_IleL*CIleL_ZEN_gluc) / 
           (Km_Hydro + (CIleL_ZEN_gluc)))
    
    # Small intestinal tissue 
    dAJej_aZEL <- ((Vmax_SI_aZEL*CJej_ZEN*ZEN_fut) / 
                     (Km_SI_aZEL + (CJej_ZEN*ZEN_fut))) +
      + QJej*(CB_aZEL-CVJej_aZEL) +
      - ((Vmax_SI_aZEL_gluc*CJej_aZEL*aZEL_fut) / 
           (Km_SI_aZEL_gluc + (CJej_aZEL*aZEL_fut)))
    
    dAIle_aZEL <- CLile_aZEL*CIleL_aZEL + 
      QIle*(CB_aZEL-CVIle_aZEL) 
    
    # Liver
    dAL_aZEL <- + QJej*CVJej_aZEL + QIle*CVIle_aZEL + 
      (QL+QDuo)*CB_aZEL - (QSI+QL)*CVL_aZEL + 
      + ((Vmax_L_aZEL*CL_ZEN*ZEN_fut) / 
           (Km_L_aZEL + (CL_ZEN*ZEN_fut))) +
      - ((Vmax_L_aZEL_gluc*CL_aZEL*aZEL_fut) / 
           (Km_L_aZEL_gluc + (CL_aZEL*aZEL_fut))) 
    
    # Rapidly Perfused 
    dAR_aZEL <- QR*(CB_aZEL-CVR_aZEL) 
    
    # Slowly Perfused 
    dAS_aZEL <- QS*(CB_aZEL-CVS_aZEL) 
    
    # Fat 
    dAF_aZEL <- QF*(CB_aZEL-CVF_aZEL) 
    
    # Blood
    dAB_aZEL <- - (QSI + QL + QR + QS + QF)*CB_aZEL +  
      + (QSI + QL)*CVL_aZEL + 
      + QF*CVF_aZEL + QR*CVR_aZEL + QS*CVS_aZEL +
      - GFR*CB_aZEL*aZEL_fup +
      + ka_Lymph*AIleL_aZEL
    
    # Urine
    dAUR_aZEL <- GFR*CB_aZEL*aZEL_fup
    
    # Feces
    dAFE_aZEL <- k_ile*AIleL_aZEL
    
    
    ### α-Zearalenol glucuronide ----
    
    
    # Intestinal lumen 
    dADuoL_aZEL_gluc <- CL_Bile*CBile_aZEL_gluc - k_duo*ADuoL_aZEL_gluc 
    
    dAJejL_aZEL_gluc <- + k_duo*ADuoL_aZEL_gluc - k_jej*AJejL_aZEL_gluc +
      + ((Vmax_MRP2_IN*CJej_IC_aZEL_gluc*aZEL_gluc_fut)/
           (Km_MRP2 + (CJej_IC_aZEL_gluc*aZEL_gluc_fut))) +
      + ((Vmax_BCRP_IN*CJej_IC_aZEL_gluc*aZEL_gluc_fut)/
           (Km_BCRP + (CJej_IC_aZEL_gluc*aZEL_gluc_fut))) 
    
    dAIleL_aZEL_gluc <- k_jej*AJejL_aZEL_gluc - k_ile*AIleL_aZEL_gluc +
      - ((Vmax_Hydro_IleL*CIleL_ZEN_gluc) / 
           (Km_Hydro + (CIleL_ZEN_gluc)))
    
    # Small Intestinal tissue
    dAJej_IC_aZEL_gluc <- ((Vmax_SI_aZEL_gluc*CJej_aZEL*aZEL_fut) / 
                             (Km_SI_aZEL_gluc + (CJej_aZEL*aZEL_fut))) +
      - ((Vmax_MRP2_IN*CJej_IC_aZEL_gluc*aZEL_gluc_fut) / 
           (Km_MRP2 + (CJej_IC_aZEL_gluc*aZEL_gluc_fut))) +
      - ((Vmax_BCRP_IN*CJej_IC_aZEL_gluc*aZEL_gluc_fut)/
           (Km_BCRP + (CJej_IC_aZEL_gluc*aZEL_gluc_fut))) 
    
    dAJej_EC_aZEL_gluc <- QJej*(CB_aZEL_gluc-CVJej_EC_aZEL_gluc) 
    
    # Liver
    dAL_IC_aZEL_gluc <- ((Vmax_L_aZEL_gluc*CL_aZEL*aZEL_fut) / 
                           ((Km_L_aZEL_gluc) + (CL_aZEL*aZEL_fut))) +
      - ((Vmax_MRP2_LI*CL_IC_aZEL_gluc*aZEL_gluc_fut) / 
           (Km_MRP2 + (CL_IC_aZEL_gluc*aZEL_gluc_fut))) + 
      - ((Vmax_MRP3_LI*CL_IC_aZEL_gluc*aZEL_gluc_fut) / 
           (Km_MRP3 + (CL_IC_aZEL_gluc*aZEL_gluc_fut))) +
      - ((Vmax_BCRP_LI*CL_IC_aZEL_gluc*aZEL_gluc_fut) / 
           (Km_BCRP + (CL_IC_aZEL_gluc*aZEL_gluc_fut))) 
    
    dAL_EC_aZEL_gluc <- + QJej*CVJej_EC_aZEL_gluc + 
      + (QL+QDuo+QIle)*CB_aZEL_gluc - (QSI + QL)*CVL_EC_aZEL_gluc +
      + ((Vmax_MRP3_LI*CL_IC_aZEL_gluc*aZEL_gluc_fut) / 
           (Km_MRP3 + (CL_IC_aZEL_gluc*aZEL_gluc_fut)))
    
    # Bile
    dABile_aZEL_gluc <- ((Vmax_MRP2_LI*CL_IC_aZEL_gluc*aZEL_gluc_fut) / 
                           (Km_MRP2 + (CL_IC_aZEL_gluc*aZEL_gluc_fut))) +
      + ((Vmax_BCRP_LI*CL_IC_aZEL_gluc*aZEL_gluc_fut) / 
           (Km_BCRP + (CL_IC_aZEL_gluc*aZEL_gluc_fut))) +
      - CL_Bile*CBile_aZEL_gluc 
    
    # Rapidly Perfused 
    dAR_aZEL_gluc <- QR_gluc*(CB_aZEL_gluc-CVR_aZEL_gluc)
    
    # Slowly Perfused 
    dAS_aZEL_gluc <- QS_gluc*(CB_aZEL_gluc-CVS_aZEL_gluc)
    
    # Fat
    dAF_aZEL_gluc <- QF*(CB_aZEL_gluc-CVF_aZEL_gluc) 
    
    # Blood
    dAB_aZEL_gluc <- - (QSI + QL + QK + QR_gluc + QS_gluc + QF)*CB_aZEL_gluc +   
      + (QSI + QL)*CVL_EC_aZEL_gluc + 
      + QF*CVF_aZEL_gluc + QK*CVK_aZEL_gluc + QR_gluc*CVR_aZEL_gluc + QS_gluc*CVS_aZEL_gluc +
      - GFR*CB_aZEL_gluc*aZEL_gluc_fup 
    
    dAK_aZEL_gluc <- QK*(CB_aZEL_gluc-CVK_aZEL_gluc) + 
      + GFR*CB_aZEL_gluc*aZEL_gluc_fup +
      - ((Vmax_MRP2_KI*CK_aZEL_gluc*aZEL_gluc_fut) / 
           (Km_MRP2 + (CK_aZEL_gluc*aZEL_gluc_fut)))
    
    # Urine
    dAUR_aZEL_gluc <- ((Vmax_MRP2_KI*CK_aZEL_gluc*aZEL_gluc_fut) / 
                         (Km_MRP2 + (CK_aZEL_gluc*aZEL_gluc_fut)))
    
    # Feces
    dAFE_aZEL_gluc <- k_ile*AIleL_aZEL_gluc
    
    
    ## Mass balance ----
    ATotal <- ODOSEumol + IVDOSEumol
    ACalculated <- AOD_ZEN + 
      
      ADuoL_ZEN + AJejL_ZEN + AIleL_ZEN + 
      ADuo_ZEN + AJej_ZEN + AIle_ZEN + 
      AL_ZEN + 
      AR_ZEN +
      AS_ZEN +
      AF_ZEN + 
      AB_ZEN + 
      AUR_ZEN + 
      AFE_ZEN + 
      
      A_bZEL +
      
      ADuoL_ZEN_gluc + AJejL_ZEN_gluc + AIleL_ZEN_gluc + 
      AJej_IC_ZEN_gluc + AJej_EC_ZEN_gluc + 
      AL_IC_ZEN_gluc + AL_EC_ZEN_gluc +
      ABile_ZEN_gluc + 
      AR_ZEN_gluc +
      AS_ZEN_gluc +
      AF_ZEN_gluc + 
      AB_ZEN_gluc + 
      AK_ZEN_gluc +
      AUR_ZEN_gluc + 
      AFE_ZEN_gluc +
      
      AIleL_aZEL +
      AJej_aZEL + AIle_aZEL +
      AL_aZEL + 
      AR_aZEL +
      AS_aZEL +
      AF_aZEL +
      AB_aZEL +
      AUR_aZEL +
      AFE_aZEL +
      
      ADuoL_aZEL_gluc + AJejL_aZEL_gluc + AIleL_aZEL_gluc + 
      AJej_IC_aZEL_gluc + AJej_EC_aZEL_gluc +
      AL_IC_aZEL_gluc + AL_EC_aZEL_gluc +
      ABile_aZEL_gluc + 
      AR_aZEL_gluc + 
      AS_aZEL_gluc +
      AF_aZEL_gluc +
      AB_aZEL_gluc +
      AK_aZEL_gluc +
      AUR_aZEL_gluc +
      AFE_aZEL_gluc 
    
    MB <- ATotal - ACalculated + 1
    
    
    
    
    
    # End
    list(c(
      dAOD_ZEN,
      dADuoL_ZEN,dAJejL_ZEN, dAIleL_ZEN,
      dADuo_ZEN, dAJej_ZEN, dAIle_ZEN,
      dAL_ZEN,
      dAR_ZEN,
      dAS_ZEN,
      dAF_ZEN,
      dAB_ZEN,
      dAUR_ZEN,
      dAFE_ZEN,
      
      dA_bZEL,
      
      dADuoL_ZEN_gluc, dAJejL_ZEN_gluc, dAIleL_ZEN_gluc, 
      dAJej_IC_ZEN_gluc, dAJej_EC_ZEN_gluc,
      dAL_IC_ZEN_gluc, dAL_EC_ZEN_gluc ,
      dABile_ZEN_gluc,
      dAR_ZEN_gluc,
      dAS_ZEN_gluc,
      dAF_ZEN_gluc,
      dAB_ZEN_gluc,
      dAK_ZEN_gluc,
      dAUR_ZEN_gluc,
      dAFE_ZEN_gluc ,
      
      dAIleL_aZEL, 
      dAJej_aZEL, dAIle_aZEL,
      dAL_aZEL,
      dAR_aZEL,
      dAS_aZEL,
      dAF_aZEL,
      dAB_aZEL,
      dAUR_aZEL,
      dAFE_aZEL,
      
      dADuoL_aZEL_gluc, dAJejL_aZEL_gluc, dAIleL_aZEL_gluc, 
      dAJej_IC_aZEL_gluc, dAJej_EC_aZEL_gluc,
      dAL_IC_aZEL_gluc, dAL_EC_aZEL_gluc,
      dABile_aZEL_gluc,
      dAR_aZEL_gluc,
      dAS_aZEL_gluc,
      dAF_aZEL_gluc,
      dAB_aZEL_gluc,
      dAK_aZEL_gluc,
      dAUR_aZEL_gluc,
      dAFE_aZEL_gluc 
    ),
    c(ATotal = ATotal,
      ACalculated = ACalculated,
      MB = MB,
      CDuoL_ZEN = CDuoL_ZEN, CJejL_ZEN = CJejL_ZEN, CIleL_ZEN = CIleL_ZEN, 
      CDuo_ZEN = CDuo_ZEN, CJej_ZEN = CJej_ZEN, CIle_ZEN = CIle_ZEN,
      CL_ZEN = CL_ZEN,
      CR_ZEN = CR_ZEN ,
      CS_ZEN = CS_ZEN,
      CF_ZEN = CF_ZEN,
      CB_ZEN = CB_ZEN,
      
      CDuoL_ZEN_gluc = CDuoL_ZEN_gluc, CJejL_ZEN_gluc = CJejL_ZEN_gluc, CIleL_ZEN_gluc = CIleL_ZEN_gluc, 
      CJej_ZEN_gluc = CJej_IC_ZEN_gluc, 
      CL_ZEN_gluc = CL_IC_ZEN_gluc, 
      CBile_ZEN_gluc = CBile_ZEN_gluc,
      CR_ZEN_gluc = CR_ZEN_gluc,
      CS_ZEN_gluc = CS_ZEN_gluc,
      CF_ZEN_gluc = CF_ZEN_gluc,
      CB_ZEN_gluc = CB_ZEN_gluc,
      CK_ZEN_gluc = CK_ZEN_gluc,
      
      CIleL_aZEL = CIleL_aZEL, 
      CJej_aZEL = CJej_aZEL, CIle_aZEL = CIle_aZEL,
      CL_aZEL = CL_aZEL,
      CR_aZEL = CR_aZEL,
      CS_aZEL = CS_aZEL,
      CF_aZEL = CF_aZEL,
      CB_aZEL = CB_aZEL,
      
      CDuoL_aZEL_gluc = CDuoL_aZEL_gluc, CJejL_aZEL_gluc = CJejL_aZEL_gluc, CIleL_aZEL_gluc = CIleL_aZEL_gluc, 
      CJej_aZEL_gluc = CJej_IC_aZEL_gluc, 
      CL_aZEL_gluc = CL_IC_aZEL_gluc, 
      CBile_aZEL_gluc = CBile_aZEL_gluc,
      CR_aZEL_gluc = CR_aZEL_gluc,
      CS_aZEL_gluc = CS_aZEL_gluc,
      CF_aZEL_gluc = CF_aZEL_gluc,
      CB_aZEL_gluc = CB_aZEL_gluc,
      CK_aZEL_gluc = CK_aZEL_gluc
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

Select_AdminType <- "IV" # IV # PO
Select_Author <- "Q.Lu" # Shin # Q.Lu # Videmann # Malli

BW <- if(Select_Author == "Shin"){0.26
}else if(Select_Author == "Q.Lu"){0.2
}else if(Select_Author == "Malli"){0.2
}else if(Select_Author == "Videmann"){0.26}

# Oral dose
ODOSEmg <- if(Select_Author == "Shin"){8
}else if(Select_Author == "Q.Lu"){6.65
}else if(Select_Author == "Malli"){3
}else if(Select_Author == "Videmann"){1}                                        # mg/kg bw IV dose
ODOSEumol <- ((ODOSEmg*1000)/MW_ZEN)*BW                                         # μmol

# IV dose
IVDOSEmg <- if(Select_AdminType != "IV"){0
}else if(Select_Author == "Shin"){2.25 #1 #4 #2 # 4 #8 #1.13 #2.25
}else if(Select_Author == "Q.Lu"){6.65}

IVDOSEumol <-((IVDOSEmg*1000)/MW_ZEN)*BW                                        # μmol

names(ODOSEumol) <- "ODOSEumol"
names(IVDOSEumol) <- "IVDOSEumol"
parm.c <- c(Parameters, c(ODOSEumol, IVDOSEumol))

timeSeq <- timeSeq                                                 # per minute


## Initials ----

A_init <- c(
  AOD_ZEN = ifelse(Select_AdminType != "IV", ODOSEumol, 0), 
  
  ADuoL_ZEN = 0, AJejL_ZEN = 0, AIleL_ZEN = 0, 
  ADuo_ZEN = 0, AJej_ZEN = 0, AIle_ZEN = 0, 
  AL_ZEN = 0, 
  AR_ZEN = 0,
  AS_ZEN = 0,
  AF_ZEN = 0, 
  AB_ZEN = ifelse(Select_AdminType != "PO", IVDOSEumol, 0), 
  AUR_ZEN = 0, 
  AFE_ZEN = 0, 
  
  A_bZEL = 0,
  
  ADuoL_ZEN_gluc = 0, AJejL_ZEN_gluc = 0, AIleL_ZEN_gluc = 0, 
  AJej_IC_ZEN_gluc = 0, AJej_EC_ZEN_gluc = 0, 
  AL_IC_ZEN_gluc = 0, AL_EC_ZEN_gluc = 0,
  ABile_ZEN_gluc = 0, 
  AR_ZEN_gluc = 0,
  AS_ZEN_gluc = 0,
  AF_ZEN_gluc = 0, 
  AB_ZEN_gluc = 0, 
  AK_ZEN_gluc = 0,
  AUR_ZEN_gluc = 0, 
  AFE_ZEN_gluc = 0,
  
  AIleL_aZEL = 0, 
  AJej_aZEL = 0, AIle_aZEL = 0,
  AL_aZEL = 0, 
  AR_aZEL = 0,
  AS_aZEL = 0,
  AF_aZEL = 0,
  AB_aZEL = 0,
  AUR_aZEL = 0,
  AFE_aZEL = 0,
  
  ADuoL_aZEL_gluc = 0, AJejL_aZEL_gluc = 0, AIleL_aZEL_gluc = 0, 
  AJej_IC_aZEL_gluc = 0, AJej_EC_aZEL_gluc = 0,
  AL_IC_aZEL_gluc = 0, AL_EC_aZEL_gluc = 0,
  ABile_aZEL_gluc = 0, 
  AR_aZEL_gluc = 0, 
  AS_aZEL_gluc = 0,
  AF_aZEL_gluc = 0,
  AB_aZEL_gluc = 0,
  AK_aZEL_gluc = 0,
  AUR_aZEL_gluc = 0,
  AFE_aZEL_gluc = 0
  
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
# ggsave(
#   filename = here(OUTPUT,"MB.Rat.png"),
#   dpi = 600,
#   width = 11,
#   height = 7,
#   units = "cm"
# )

MW_ZEN_gluc <- pull(PhysChem[PhysChem$Substance == "Zearalenone-14O-glucuronide", 
                        "MolecularWeight"])
MW_aZEL <- pull(PhysChem[PhysChem$Substance == "a-Zearalenol", 
                        "MolecularWeight"])
MW_aZEL_gluc <- pull(PhysChem[PhysChem$Substance == "a-Zearalenol-glucuronide", 
                         "MolecularWeight"])

replacements <- c(
  "CDuoL" = "Duodenum Lumen",
  "CDuo"  = "Duodenum",
  "CIleL" = "Ileum Lumen",
  "CIle" = "Ileum",
  "CJejL" = "Jejunum Lumen", 
  "CJej" = "Jejunum", 
  "CL" = "Liver", 
  "CF" = "Fat", 
  "CR" = "Rapidly perfused", 
  "CS" = "Slowly perfused",
  "CBile" = "Bile",
  "CB" = "Blood", 
  "CK" = "Kidney",
  "_gluc" = "-gluc",
  "time" = "Time"
)


Concentrations <- PBK_OUTPUT %>% 
  select(c("time", starts_with("C"))) %>% 
  set_names(~(.) %>%
              str_replace_all(replacements)) %>% 
  pivot_longer( !Time , 
               names_to = c("Tissue", "Substance"),
               names_sep = c("_"),
               values_to = "Concentration") %>% 
  mutate(Concentration = ifelse(Substance == "ZEN", Concentration*MW_ZEN,       # correcting units
                                ifelse(Substance == "ZEN_gluc",Concentration*MW_ZEN_gluc,
                                       ifelse(Substance == "aZEL", Concentration*MW_aZEL,
                                              Concentration*MW_aZEL_gluc)
                                       )))

C.Pre_long <- Concentrations

VSI <- FV_SI*BW 
VDuo <- VSI*FV_Duo
VJej <- VSI*FV_Jej
VIle <- VSI*FV_Ile

VSIL <- FV_SIL*BW 
VDuoL <- VSIL*FV_duol 
VJejL <- VSIL*FV_jejl 
VIleL <- VSIL*FV_ilel 



Clean.replacements <- c(
  "CIn" = "Intestine",
  "CL" = "Liver", 
  "CF" = "Fat", 
  "CR" = "Rapidly perfused", 
  "CS" = "Slowly perfused",
  "CBile" = "Bile",
  "CB" = "Blood", 
  "CK" = "Kidney",
  "_gluc" = "-gluc",
  "time" = "Time"
)


Clean.C <- PBK_OUTPUT %>% 
  select(c("time", starts_with("C"))) %>% 
  mutate("CInT_ZEN" = ((CDuo_ZEN*VDuo) + (CJej_ZEN*VJej) + (CIle_ZEN*VIle))/(
    VDuo + VJej + VIle),
    "CInL_ZEN" = ((CDuoL_ZEN*VDuoL) + (CJejL_ZEN*VJejL) + (CIleL_ZEN*VIleL))/(
      VDuoL + VJejL + VIleL),
    "CInT_ZEN_gluc" = ((CJej_ZEN_gluc*VJej))/(
      VDuo + VJej + VIle),
    "CInL_ZEN_gluc" = ((CDuoL_ZEN_gluc*VDuoL) + (CJejL_ZEN_gluc*VJejL) + (CIleL_ZEN_gluc*VIleL))/(
      VDuoL + VJejL + VIleL),
    "CInT_aZEL" = ((CJej_aZEL*VJej) + (CIle_aZEL*VIle))/(
      VDuo + VJej + VIle),
    "CInL_aZEL" = ((CIleL_aZEL*VIleL))/(
      VDuoL + VJejL + VIleL),
    "CInT_aZEL_gluc" = ((CJej_aZEL_gluc*VJej))/(
      VDuo + VJej + VIle),
    "CInL_aZEL_gluc" = ((CDuoL_aZEL_gluc*VDuoL) + (CJejL_aZEL_gluc*VJejL) + (CIleL_aZEL_gluc*VIleL))/(
      VDuoL + VJejL + VIleL)
    ) %>% 
  mutate(
    "CIn_ZEN" = CInT_ZEN, 
    "CIn_ZEN_gluc" = CInT_ZEN_gluc, 
    "CIn_aZEL" = CInT_aZEL, 
    "CIn_aZEL_gluc" = CInT_aZEL_gluc, 
  ) %>% 
  select(!c(starts_with("CDuo"), 
            starts_with("CIle"), 
            starts_with("CJej"), 
            starts_with("CInT"),
            starts_with("CInL"))) %>% 
  set_names(~(.) %>%
              str_replace_all(Clean.replacements)) %>% 
  pivot_longer( !Time , 
                names_to = c("Tissue", "Substance"),
                names_sep = c("_"),
                values_to = "Concentration") %>% 
  mutate(Concentration = ifelse(Substance == "ZEN", Concentration*MW_ZEN,       # correcting units
                                ifelse(Substance == "ZEN_gluc",Concentration*MW_ZEN_gluc,
                                       ifelse(Substance == "aZEL", Concentration*MW_aZEL,
                                              Concentration*MW_aZEL_gluc)
                                )))


# Plot Blood
CBlood.plot <- C.Pre_long %>% filter(Tissue == "Blood") %>%
  ggplot() + geom_path(aes(x = Time, y = Concentration)) +
  facet_wrap(~ Substance, scales = "free") +
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
ggsave(
  filename = here(OUTPUT,paste0(test, "BloodC.Rat.png")),
  dpi = 600,
  width = 11,
  height = 7,
  units = "cm"
)

# Plot all tissues

CTissues.plot <- Clean.C %>% 
  ggplot() + 
  geom_path(aes(x = Time, y = Concentration, colour = Substance)) +
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
#   filename = here(OUTPUT, paste0(test, "TissueC.Rat.png")),
#   dpi = 600,
#   width = 11,
#   height = 7,
#   units = "cm"
# )

# write_csv(Clean.C, here(OUTPUT, "C.Tissues.Rat.csv"))

# ZEN
AUC_ZEN <- trapz(PBK_OUTPUT[ , "time"], PBK_OUTPUT[ , "CB_ZEN"])  # umol*h/L
print(AUC_ZEN)

time <- PBK_OUTPUT[ , "time"] # h
conc_ZEN <- PBK_OUTPUT[ , "CB_ZEN"] # uM 
Cmax_ZEN <- max(conc_ZEN)
Tmax_ZEN <- time[which.max(conc_ZEN)]
tlast_ZEN <- max(time[conc_ZEN > 0])
half_life_ZEN <- pk.calc.half.life(
  conc_ZEN,
  time,
  Tmax_ZEN,
  tlast_ZEN
)
HalfLife_ZEN <- half_life_ZEN$half.life  # half-life in h
print(HalfLife_ZEN)

# ZEN gluc
AUC_ZEN_gluc <- trapz(PBK_OUTPUT[ , "time"], PBK_OUTPUT[ , "CB_ZEN_gluc"])  # umol*h/L
print(AUC_ZEN_gluc)

time <- PBK_OUTPUT[ , "time"] # h
conc_ZEN_gluc <- PBK_OUTPUT[ , "CB_ZEN_gluc"] # uM 
Cmax_ZEN_gluc <- max(conc_ZEN_gluc)
Tmax_ZEN_gluc <- time[which.max(conc_ZEN_gluc)]
tlast_ZEN_gluc <- max(time[conc_ZEN_gluc > 0])
half_life_ZEN_gluc <- pk.calc.half.life(
  conc_ZEN_gluc,
  time,
  Tmax_ZEN_gluc,
  tlast_ZEN_gluc
)
HalfLife_ZEN_gluc <- half_life_ZEN_gluc$half.life  # half-life in h
print(HalfLife_ZEN_gluc)

# CB_aZEL
AUC_aZEL <- trapz(PBK_OUTPUT[ , "time"], PBK_OUTPUT[ , "CB_aZEL"])  # umol*h/L
print(AUC_aZEL)

time <- PBK_OUTPUT[ , "time"] # h
conc_aZEL <- PBK_OUTPUT[ , "CB_aZEL"] # uM 
Cmax_aZEL <- max(conc_aZEL)
Tmax_aZEL <- time[which.max(conc_aZEL)]
tlast_aZEL <- max(time[conc_aZEL > 0])
half_life_aZEL <- pk.calc.half.life(
  conc_aZEL,
  time,
  Tmax_aZEL,
  tlast_aZEL
)
HalfLife_aZEL <- half_life_aZEL$half.life  # half-life in h
print(HalfLife_aZEL)

# aZEL_gluc
AUC_aZEL_gluc <- trapz(PBK_OUTPUT[ , "time"], PBK_OUTPUT[ , "CB_aZEL_gluc"])  # umol*h/L
print(AUC_aZEL_gluc)

time <- PBK_OUTPUT[ , "time"] # h
conc_aZEL_gluc <- PBK_OUTPUT[ , "CB_aZEL_gluc"] # uM 
Cmax_aZEL_gluc <- max(conc_aZEL_gluc)
Tmax_aZEL_gluc <- time[which.max(conc_aZEL_gluc)]
tlast_aZEL_gluc <- max(time[conc_aZEL_gluc > 0])
half_life_aZEL_gluc <- pk.calc.half.life(
  conc_aZEL_gluc,
  time,
  Tmax_aZEL_gluc,
  tlast_aZEL_gluc
)
HalfLife_aZEL_gluc <- half_life_aZEL_gluc$half.life  # half-life in h
print(HalfLife_aZEL_gluc)

TK.data <- tibble(
  Mycotoxin = c("ZEN", "ZEN_gluc", "aZEL", "aZEL_gluc"),
  Cmax = c(Cmax_ZEN, Cmax_ZEN_gluc, Cmax_aZEL, Cmax_aZEL_gluc), 
  Tmax = c(Tmax_ZEN, Tmax_ZEN_gluc, Tmax_aZEL, Tmax_aZEL_gluc),
  AUC = c(AUC_ZEN, AUC_ZEN_gluc, AUC_aZEL, AUC_aZEL_gluc),
  Half_life = c(HalfLife_ZEN, HalfLife_ZEN_gluc, HalfLife_aZEL, HalfLife_aZEL_gluc), 
  Test = test
)

# Results

Results.Rat <- list(
  Metadata.Run.Rat = list(Select_AdminType = Select_AdminType,
                          Select_Author = Select_Author,
                          ODOSEmg = ODOSEmg, 
                          ODOSEumol = ODOSEumol,
                          IVDOSEmg = IVDOSEmg,
                          IVDOSEumol = IVDOSEumol,
                          timeSeq = timeSeq,
                          ZEN_BPR = pull(Kp_fus_BPR_df[Kp_fus_BPR_df$Parameter == "ZEN_BPR", "Value"]),
                          aZEL_BPR = pull(Kp_fus_BPR_df[Kp_fus_BPR_df$Parameter == "aZEL_BPR", "Value"]),
                          ZEN_gluc_BPR = pull(Kp_fus_BPR_df[Kp_fus_BPR_df$Parameter == "ZEN_gluc_BPR", "Value"]),
                          aZEL_gluc_BPR =  pull(Kp_fus_BPR_df[Kp_fus_BPR_df$Parameter == "aZEL_gluc_BPR", "Value"])), 
  Raw.Data.Rat = PBK_OUTPUT,
  C.Tissues.Rat = Clean.C,
  C.All.Rat = C.Pre_long,
  Plots = list(CBlood.plot, CTissues.plot), #MB.plot, 
  TK.data = TK.data
)

save(Results.Rat, file = here(OUTPUT, paste0(test, "Results.Rat.RData")))
