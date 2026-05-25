rm(list=ls()) # to clear out the global environment

# Load packages
library(here)
library(ggplot2)
library(deSolve)
library(tidyverse)

source(here("Script", "GSA","PBK_model_human.R"))
Parameters <- read.csv(here("Input", "Parameters_human.csv"))
P <- as.list(Parameters$Value)
names(P) <- Parameters$Parameter
RPF <- 63

NOEL <- 10/1000 # mg/kg bw
TDI <- 0.25/1000 # mg/kg bw

targetCmaxP <- c(0.490,
                 0.359,
                 0.058,
                 1.580,
                 3.810,
                 1.600,
                 1.639,
                 0.310)/1000 #uM

  process_dose <- function(testODOSEmg.BW) {
    
    A_init <- c(
      AOD_ZEN = ((testODOSEmg.BW*1000)/318.15)* 70.63106831,  
      
      ADuoL_ZEN = 0, AJejL_ZEN = 0, AIleL_ZEN = 0, 
      ADuo_ZEN = 0, AJej_ZEN = 0, AIle_ZEN = 0, 
      AL_ZEN = 0, 
      AR_ZEN = 0,
      AS_ZEN = 0,
      AF_ZEN = 0, 
      AB_ZEN = 0, 
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
    
    PBK.out <- lsoda(y = A_init,
                     times = seq(0, 1,by=1/100),
                     func =  H.PBK.model, 
                     parms = P)
    
    CmaxP_ZEN = (max(PBK.out[,"CB_ZEN"]) / 0.67)
    CmaxP_ZENeq_aZEL = max(PBK.out[,"CB_ZEN"])/0.67 + ((max(PBK.out[,"CB_aZEL"])/0.67)*RPF)
    CmaxP_ZENeq_aZEL_gluc  = (max(PBK.out[,"CB_ZEN"])/0.67) + 
      (max(PBK.out[,"CB_ZEN_gluc"])/0.55) + 
      ((max(PBK.out[,"CB_aZEL"]/0.67)*RPF)) + 
      ((max(PBK.out[,"CB_aZEL_gluc"])/0.55)*RPF)
    
    max_vals <- c(
      CmaxP_ZEN = CmaxP_ZEN,
      CmaxP_ZENeq_aZEL = CmaxP_ZENeq_aZEL,
      CmaxP_ZENeq_aZEL_gluc  = CmaxP_ZENeq_aZEL_gluc
    )
    
    # Return differences matrix
    CloseOrNot <- as.data.frame(sapply(max_vals, function(m) targetCmaxP/m))
    rownames(CloseOrNot) <- as.character(targetCmaxP)
    
    return(list(max_vals, CloseOrNot))
  }
  
  # 1. save predictions with TDI and OD estimations to meet target CmaxP
  
  CloseToInVitro <- list()
  ZENeqCmaxP <- list()
  testODmg.BW <- TDI #testODmg.BW#TDI #c("TDI" = TDI, "NOEL" = NOEL)
  for (i in seq_along(testODmg.BW)) {
    
    print(i)
    
    testODOSEmg.BW <- testODmg.BW[[i]]
    print(testODOSEmg.BW)
    out <- process_dose(testODOSEmg.BW)
    
    print(str(out))  # See what structure out has
    print(length(out))  # Check how many elements
    
    CloseToInVitro[[i]] <- out[[2]]
    ZENeqCmaxP[[i]] <- out[[1]]
    
    
  }
  
  
  # TDI_Cmax <- ZENeqCmaxP[[1]]
  # TDI_CloseToInVitro <- CloseToInVitro[[1]] 
  # TargetOD <- c(TDI_CloseToInVitro$CmaxP_ZEN*TDI,
  #               TDI_CloseToInVitro$CmaxP_ZENeq_aZEL*TDI,
  #               TDI_CloseToInVitro$CmaxP_ZENeq_aZEL_gluc*TDI)
  # write.csv(TargetOD, "TargetOD.csv") 

  # 2. predict with estimated OD
 
  CloseToInVitro <- list()
  ZENeqCmaxP <- list()
  testODmg.BW <- c(0.226064355,
                   0.165546622,
                   0.026543386,
                   0.724695818,
                   1.710424482,
                   0.733750284,
                   0.751395909,
                   0.142892616,
                   0.020631542,
                   0.015019469,
                   0.002386707,
                   0.069122382,
                   0.17266403,
                   0.07003313,
                   0.071810608,
                   0.012936907,
                   0.000102914,
                   7.54E-05,
                   1.22E-05,
                   0.000331943,
                   0.000800895,
                   0.000336147,
                   0.000344344,
                   6.59E-05
  ) 
  
  for (i in seq_along(testODmg.BW)) {
    
    # print(i)
    
    testODOSEmg.BW <- testODmg.BW[[i]]
    print(testODOSEmg.BW)
    out <- process_dose(testODOSEmg.BW)
    
    # print(str(out))  # See what structure out has
    # print(length(out))  # Check how many elements
    
    CloseToInVitro[[i]] <- out[[2]]
    ZENeqCmaxP[[i]] <- out[[1]]
    
    
  }
  
  
  
  
  Method <- c(rep("CmaxP_ZEN",8), rep("CmaxP_ZENeq_aZEL", 8), rep("CmaxP_ZENeq_aZEL_gluc", 8))
  testODE <- TargetOD 
  targetCmax <- rep(targetCmaxP,3)
    # c(rep(targetCmaxP[[1]],3), rep(targetCmaxP[[10]],3),
    #               rep(targetCmaxP[[3]],3), rep(targetCmaxP[[4]],3),
    #               rep(targetCmaxP[[5]],3), rep(targetCmaxP[[6]],3),
    #               rep(targetCmaxP[[7]],3), rep(targetCmaxP[[8]],3))
  predCmax <- c(as.data.frame(ZENeqCmaxP[[1]])["CmaxP_ZEN",1], 
                as.data.frame(ZENeqCmaxP[[2]])["CmaxP_ZEN",1], 
                as.data.frame(ZENeqCmaxP[[3]])["CmaxP_ZEN",1], 
                as.data.frame(ZENeqCmaxP[[4]])["CmaxP_ZEN",1], 
                as.data.frame(ZENeqCmaxP[[5]])["CmaxP_ZEN",1], 
                as.data.frame(ZENeqCmaxP[[6]])["CmaxP_ZEN",1], 
                as.data.frame(ZENeqCmaxP[[7]])["CmaxP_ZEN",1], 
                as.data.frame(ZENeqCmaxP[[8]])["CmaxP_ZEN",1], 
                as.data.frame(ZENeqCmaxP[[1]])["CmaxP_ZENeq_aZEL",1], 
                as.data.frame(ZENeqCmaxP[[2]])["CmaxP_ZENeq_aZEL",1], 
                as.data.frame(ZENeqCmaxP[[3]])["CmaxP_ZENeq_aZEL",1], 
                as.data.frame(ZENeqCmaxP[[4]])["CmaxP_ZENeq_aZEL",1], 
                as.data.frame(ZENeqCmaxP[[5]])["CmaxP_ZENeq_aZEL",1], 
                as.data.frame(ZENeqCmaxP[[6]])["CmaxP_ZENeq_aZEL",1], 
                as.data.frame(ZENeqCmaxP[[7]])["CmaxP_ZENeq_aZEL",1], 
                as.data.frame(ZENeqCmaxP[[8]])["CmaxP_ZENeq_aZEL",1], 
                as.data.frame(ZENeqCmaxP[[1]])["CmaxP_ZENeq_aZEL_gluc",1], 
                as.data.frame(ZENeqCmaxP[[2]])["CmaxP_ZENeq_aZEL_gluc",1], 
                as.data.frame(ZENeqCmaxP[[3]])["CmaxP_ZENeq_aZEL_gluc",1], 
                as.data.frame(ZENeqCmaxP[[4]])["CmaxP_ZENeq_aZEL_gluc",1], 
                as.data.frame(ZENeqCmaxP[[5]])["CmaxP_ZENeq_aZEL_gluc",1], 
                as.data.frame(ZENeqCmaxP[[6]])["CmaxP_ZENeq_aZEL_gluc",1], 
                as.data.frame(ZENeqCmaxP[[7]])["CmaxP_ZENeq_aZEL_gluc",1], 
                as.data.frame(ZENeqCmaxP[[8]])["CmaxP_ZENeq_aZEL_gluc",1]
                )
  Closeness <- c(as.data.frame(CloseToInVitro[[1]])[1, "CmaxP_ZEN"],
                 as.data.frame(CloseToInVitro[[2]])[2, "CmaxP_ZEN"],
                 as.data.frame(CloseToInVitro[[3]])[3, "CmaxP_ZEN"],
                 as.data.frame(CloseToInVitro[[4]])[4, "CmaxP_ZEN"],
                 as.data.frame(CloseToInVitro[[5]])[5, "CmaxP_ZEN"],
                 as.data.frame(CloseToInVitro[[6]])[6, "CmaxP_ZEN"],
                 as.data.frame(CloseToInVitro[[7]])[7, "CmaxP_ZEN"],
                 as.data.frame(CloseToInVitro[[8]])[8, "CmaxP_ZEN"],
                 
                 as.data.frame(CloseToInVitro[[9]])[1, "CmaxP_ZENeq_aZEL"],
                 as.data.frame(CloseToInVitro[[10]])[2, "CmaxP_ZENeq_aZEL"],
                 as.data.frame(CloseToInVitro[[11]])[3, "CmaxP_ZENeq_aZEL"],
                 as.data.frame(CloseToInVitro[[12]])[4, "CmaxP_ZENeq_aZEL"],
                 as.data.frame(CloseToInVitro[[13]])[5, "CmaxP_ZENeq_aZEL"],
                 as.data.frame(CloseToInVitro[[14]])[6, "CmaxP_ZENeq_aZEL"],
                 as.data.frame(CloseToInVitro[[15]])[7, "CmaxP_ZENeq_aZEL"],
                 as.data.frame(CloseToInVitro[[16]])[8, "CmaxP_ZENeq_aZEL"],
                 
                 as.data.frame(CloseToInVitro[[17]])[1, "CmaxP_ZENeq_aZEL_gluc"],
                 as.data.frame(CloseToInVitro[[18]])[2, "CmaxP_ZENeq_aZEL_gluc"],
                 as.data.frame(CloseToInVitro[[19]])[3, "CmaxP_ZENeq_aZEL_gluc"],
                 as.data.frame(CloseToInVitro[[20]])[4, "CmaxP_ZENeq_aZEL_gluc"],
                 as.data.frame(CloseToInVitro[[21]])[5, "CmaxP_ZENeq_aZEL_gluc"],
                 as.data.frame(CloseToInVitro[[22]])[6, "CmaxP_ZENeq_aZEL_gluc"],
                 as.data.frame(CloseToInVitro[[23]])[7, "CmaxP_ZENeq_aZEL_gluc"],
                 as.data.frame(CloseToInVitro[[24]])[8, "CmaxP_ZENeq_aZEL_gluc"]
  )
  
  Collect <- tibble(
    "Method" = Method,
    "testODE" = round(testODE,5),
    "targetCmax" = round(targetCmax,5),
    "predCmax" = round(predCmax,5),
    "Closeness" = round(Closeness,5)
  )
  
  # 3. Add predictions with the TDI
  CloseToInVitro <- list()
  ZENeqCmaxP <- list()
  testODmg.BW <- c(TDI, NOEL) 
  
  for (i in seq_along(testODmg.BW)) {
    
    # print(i)
    
    testODOSEmg.BW <- testODmg.BW[[i]]
    print(testODOSEmg.BW)
    out <- process_dose(testODOSEmg.BW)
    
    # print(str(out))  # See what structure out has
    # print(length(out))  # Check how many elements
    
    CloseToInVitro[[i]] <- out[[2]]
    ZENeqCmaxP[[i]] <- out[[1]]
    
    
  }
  
  Collectn <- tibble(
    "Method" = c("CmaxP_ZEN", "CmaxP_ZENeq_aZEL", "CmaxP_ZENeq_aZEL_gluc",
                 "CmaxP_ZEN", "CmaxP_ZENeq_aZEL", "CmaxP_ZENeq_aZEL_gluc"),
    "testODE" = c(rep(TDI,3), rep(NOEL,3)),
    "targetCmax" = NA,
    "predCmax" = c(as.data.frame(ZENeqCmaxP[[1]])["CmaxP_ZEN",1],
                   as.data.frame(ZENeqCmaxP[[1]])["CmaxP_ZENeq_aZEL",1],
                   as.data.frame(ZENeqCmaxP[[1]])["CmaxP_ZENeq_aZEL_gluc",1],
                   as.data.frame(ZENeqCmaxP[[2]])["CmaxP_ZEN",1],
                   as.data.frame(ZENeqCmaxP[[2]])["CmaxP_ZENeq_aZEL",1],
                   as.data.frame(ZENeqCmaxP[[2]])["CmaxP_ZENeq_aZEL_gluc",1]),
    "Closeness" = NA
  )
  
  All <- rbind(Collect, Collectn)
  
  write.csv(All, "N_Results_ReverseDosimetry.csv")
  save(ZENeqCmaxP_OED, file = "ZENeqCmaxP_OED.RData")
  save(ZENeqCmaxP, file = "ZENeqCmaxP_TDINOEL.RData")
  save(CloseToInVitro_OED, file = "CloseToInVitro_OED.RData")
  save(CloseToInVitro, file = "CloseToInVitro_TDINOEL.RData")
  