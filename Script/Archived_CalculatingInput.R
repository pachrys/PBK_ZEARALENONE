
# ---Packages------------------------------------------------------------------#

# library(rxode2)
# library(tidyverse)
# library(ggplot2)
# library(showtext)
# font_add(family = "Garamond", regular = "GARA.TTF")
# showtext_auto()
# library("here")


#*****************************************************************************#
#---Input calculations--------------------------------------------------------#
#*****************************************************************************#

#----Extracting input---------------------------------------------------------#

# Rat physiological input
Rat.datall <- read.csv(here("InputData/TisCompRat_CP.csv"), header=1)            #tissue composition mainly as R&R composition, unified by Utsey and optimised as per Musther
Rat.dat <- Rat.datall %>% filter(!tissue %in% c("RBCs", "Plasma"))               #df for all tissues except for adipose and RBCs

# Extracting Fractional flow rates F_QC
F_QC.names <- Rat.dat$tissue %>% substr(1,2)
F_QC.names <- paste("FQ", F_QC.names, sep="")
F_QC.dat <- cbind(F_QC.names, select(Rat.dat, F_QC))
F_QC.dat <- F_QC.dat %>%
  pivot_wider(names_from = F_QC.names, values_from = F_QC )

# Exrtacting Fractional volumes F_Vo
F_Vo.names <- Rat.dat$tissue %>% substr(1,2)
F_Vo.names <- paste("FV", F_Vo.names, sep="")
F_Vo.dat <- cbind(F_Vo.names, select(Rat.dat, F_Vo))
F_Vo.dat <- F_Vo.dat %>%
  pivot_wider(names_from = F_Vo.names, values_from = F_Vo )

# ZEN physicochemical input                                                       #ZEN as monoprotic acid. It has 2 pKas 8.54 and 10.38. At pH 7.4 ZEN is >90 persent not ionized.
ZEN.dat <- read.csv(here("InputData/Zen_input.csv"), header=1)     #from {Agahi, 2020 #456} and chemicalize, BP from Shin, Hong et al.2009


#----Calculating-------------------------------------------------------------#

# Calculating fraction unbound in plasma                                          #Reference: Lobell M, Sivarajah V. In silico prediction of aqueous solubility, human plasma protein binding and volume of distribution of compounds from calculated pKa and AlogP98 values. Mol Divers. 2003;7(1):69-87. doi: 10.1023/b:modi.0000006562.93049.36. PMID: 14768905.
fup <- 1/(1 + 10^(0.4485*ZEN.dat$logP - 0.4782))

# Calculating ionization status                                                   #Based on Schmitt for rat https://www.sciencedirect.com/science/article/pii/S0887233307002573?via%3Dihub, as firstly coded in: Reference: https://github.com/metrumresearchgroup/PBPK_PC/blob/master/script/CalcKp_R%26R.R
#                                                                                 #Assumption that Zearalenone is a monoprotic acid; in reality it's a diprotic acid with pKa2=10.38. Insignificant difference in tissue:plasma partition coefficients was observed when diprotic acid equations were tested. Check.As.Diprotic.RnR.R file.
#                                                                                 #Zearalenone is not ionized at physiological pH. RnR's equation for neutrals was overestimating the partition in all tissues by 7%. Check As.Neutral.RnR.R
K_n_pl <- 10^ZEN.dat$logP                                                         #Calculates neutral phospholipids:water partition coefficient (logP here is logMA of Shmitt)
K_protein <- ((0.81 + 0.11 * K_n_pl)/24.92)*5                                     #protein:water partition; Schmitt, Walter (2008)
pH <- Rat.datall$pH                                                               #tissue pH; list of 11 values (per tissue)
alpha <- 1e-3                                                                     #ratio between ditribution coefficient at given pH (D) and that in neutral form (D0)

W <- 10^(Rat.dat$pH -ZEN.dat$pKa)                                                 #Calculates ZEN ionisation status depending on tissue pH
K_n_l <- K_n_pl*(((1-alpha)/(1+W))+alpha)                                         #Calculates neutral lipids:neutral phospholipids partition coefficient
K_a_pl <- K_n_pl*((1/(1+W))+0.05*(1-(1/(1+W))))                                   #Calculates acid phospholipids:neutral lipids partition coefficient based on ionisation status

# Calculating Tissue:Plasma partition coefficients
kp <- (Rat.dat$f_water +
         (K_n_l*Rat.dat$f_n_l) +
         (K_n_pl*Rat.dat$f_n_pl) +
         (K_a_pl*Rat.dat$f_a_pl) +
         (K_protein*Rat.dat$f_proteins))*fup

name <- Rat.dat$tissue %>% substr(1,2)                                           #creates a new vector names by taking all the values from the tissue variable in Rat.datall df and then substracts the first two characters from each value
name <- paste("P", name, sep="")                                                 #Pastes "P" to the names
PartiCoef <- as.list(c(kp, fup))
names(PartiCoef) <- c(name, "fup")

Calc.init.Input <- cbind(ZEN.dat,data.frame(PartiCoef),F_QC.dat, F_Vo.dat) 
#col.names <- c("Ionisation", "Tissue", "pH")
#Ionisation <- cbind(W, pH,name)

#----Calculating Distribution Parameters--------------------------------------#
# Tissue:Blood partition coefficients                                             #Tissue:Blood Partition coefficients = BloodPlasma ratio (BP) * Tissue:Plasma Partition coefficient (PTi)
PLib =ZEN.dat$BP*Calc.init.Input$PLi                                             #liver
PKib =ZEN.dat$BP*Calc.init.Input$PKi                                             #kidney
PAdb =ZEN.dat$BP*Calc.init.Input$PAd                                             #adipose
PInb =ZEN.dat$BP*Calc.init.Input$PGu                                             #intestine
PSlb =ZEN.dat$BP*(((Calc.init.Input$FVBo*Calc.init.Input$PBo) + 
                     (Calc.init.Input$FVBr*Calc.init.Input$PBr) +
              (Calc.init.Input$FVMu*Calc.init.Input$PMu) + 
              (Calc.init.Input$FVSk*Calc.init.Input$PSk)) /
             (Calc.init.Input$FVBo + 
                Calc.init.Input$FVBr + 
                Calc.init.Input$FVMu + 
                Calc.init.Input$FVSk))                                          #slowly perfused = partition coefficient of each of these organs * fractional volume of the respective organ / sum of the fractional volume of all these organs
PRab =ZEN.dat$BP*(((Calc.init.Input$FVHe*Calc.init.Input$PHe) + 
              (Calc.init.Input$FVLu*Calc.init.Input$PLu) + 
              (Calc.init.Input$FVSp*Calc.init.Input$PSp) + 
              (Calc.init.Input$FVGu*Calc.init.Input$PKi)) /
             (Calc.init.Input$FVHe + 
                Calc.init.Input$FVLu + 
                Calc.init.Input$FVSp + 
                Calc.init.Input$FVKi))                                          #rapidly perfused; calculated similarly to rapidly perfused

# Free fraction in blood
fub = Calc.init.Input$fu*ZEN.dat$BP                                              #equation from {Grandoni, 2019 #457}



#---------------Collect all parameters together-------------------------#
CalcInput <- cbind(
  Calc.init.Input,
  PLib=PLib,
  PInb=PInb,
  PKib=PKib,
  PAdb=PAdb,
  PSlb=PSlb,
  PRab=PRab,
  fub=fub
)

print(CalcInput)

#---End-----------------------------------------------------------------------#
#*****************************************************************************#
