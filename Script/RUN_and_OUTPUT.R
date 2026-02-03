# --------------------------------------------------------------------------- #
# Script for generating PBK model results
# By: Chrysanthi Pachoulide
# Date: 30-01-2026
# --------------------------------------------------------------------------- #

# Without physiological changes with lifestage (depending on age) ----
RUNandOUT <- function(exp, 
                      Tinput, tinterval, expSTOP, 
                      # exp1STOP, exp2STOP, exp3STOP,
                      expAGE = NULL, expBW = NULL, sex = "M", 
                      Tstart, Tstop, Dt) {
  
  
  # Load R files 
  source(here("Script", "CALC_Parameters.R"), local = TRUE)
  source(here("Script", "PBK_model.R"), local = TRUE)  
  
  ## Add state ####
  A_init <- c(
    
    AOD_ZEN = 0, #exp*expBW, #0, 
    
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
    AFE_aZEL_gluc = 0#,
    
    # ATotal = 0
    
  )
  
  
  # Calculate parameters
  parms = PARAMS(expAGE, expBW, sex) 
  # Dosing <- c(Paramter = "ODOSEumol",
  #             Value = exp)
  # parms <- rbind(parms, Dosing)
  write.csv(parms, file = here(OUTPUT, "ModelParameters.csv"), row.names = FALSE)
  
  parm.c <- as.numeric(parms$Value)
  names(parm.c) <- parms$Parameter
  
  ODOSEumol <-  exp*parm.c[["BW"]]
  Dosingevents <- data.frame(var = "AOD_ZEN",
                             time =  seq(Tstart,Tstop,by=Dt)) 
  
  Dosingevents$value <- ifelse(
    Dosingevents$time == expSTOP, 
    # Dosingevents$time %in% c(expSTOP, exp1STOP, exp2STOP, exp3STOP), 
    ODOSEumol, 0) 
  Dosingevents$method <- "add"
    # mutate(value = ifelse(time <= expSTOP,
    #                       seq(0,ODOSEumol, length.out = sum(time <= expSTOP)),
    #                       0)) %>%  
    # mutate(method = "add")
    # 
  ## Run the PBK model ####
  PBK_OUTPUT <- PBK_RUN(y = A_init, 
                        parms = parm.c, 
                        times = seq(Tstart,Tstop,by=Dt),
                        events = Dosingevents
                        )
  
  ## Plot Mass Balance/Error ####
  MB.df <- PBK_OUTPUT %>% select(time, ACalculated)
  Dosingevents <- Dosingevents %>% filter(row_number() <= n()-1)
  MB.df$ATotal <- c(0,cumsum(Dosingevents$value))
  MB.df$MB <- MB.df$ATotal - MB.df$ACalculated + 1                                               # rounding significance points
  MB.df$ERROR <- (MB.df$ATotal - MB.df$ACalculated) / MB.df$ACalculated * 100
  MB.df$ERROR <- round(MB.df$ERROR, 10)                                           # rounding significance points
  MB.plot <- ggplot(data = MB.df)+
    geom_line(aes(x = time, y = ERROR, color = "ERROR")) +
    geom_line(aes(x = time, y = MB, color = "MB")) +
    scale_color_manual(values = c("ERROR" = "blue", "MB" = "black"),
                       name = NULL) +
    CP_theme +
    ylab("MB / ERROR")
  MB.plot
  ggsave(
    filename = here(OUTPUT,"MB.png"),
    dpi = 600,
    width = 11,
    height = 7,
    units = "cm"
  )
  
  ## Plot organ concentrations ####
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
    "CCOL" = "Colon",
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
    mutate(Concentration = ifelse(Substance == "ZEN", Concentration*MW_ZEN,     # correcting units
                                  ifelse(Substance == "ZEN_gluc",Concentration*MW_ZEN_gluc,
                                         ifelse(Substance == "aZEL", Concentration*MW_aZEL,
                                                Concentration*MW_aZEL_gluc)
                                  )))
  
  C.Pre_long <- Concentrations
  
  BW <- parm.c[["BW"]]                                         # kg
  
  VSI <- parm.c[["FV_SI"]]*BW 
  VDuo <- VSI*parm.c[["FV_Duo"]]
  VJej <- VSI*parm.c[["FV_Jej"]]
  VIle <- VSI*parm.c[["FV_Ile"]]
  
  VSIL <- parm.c[["FV_SIL"]]*BW # fractional volume of small intestinal lumen
  VDuoL <- VSIL*parm.c[["FV_duol"]] 
  VJejL <- VSIL*parm.c[["FV_jejl"]] 
  VIleL <- VSIL*parm.c[["FV_ilel"]] 
  
  
  Clean.replacements <- c(
    "CIn" = "Intestine",
    "CL" = "Liver", 
    "CF" = "Fat", 
    "CR" = "Rapidly perfused", 
    "CS" = "Slowly perfused",
    "CBile" = "Bile",
    "CB" = "Blood", 
    "CK" = "Kidney",
    "CUrine" = "Urine",
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
      "CInT_aZEL" = ((CIle_aZEL*VIle))/( 
        VDuo + VJej + VIle),
      "CInL_aZEL" = ((CIleL_aZEL*VIleL))/(
        VDuoL + VJejL + VIleL),
      "CInL_aZEL_gluc" = ((CDuoL_aZEL_gluc*VDuoL) + (CIleL_aZEL_gluc*VIleL))/( #(CJejL_aZEL_gluc*VJejL)
        VDuoL + VJejL + VIleL)
    ) %>% 
    mutate(
      "CIn_ZEN" = CInT_ZEN, 
      "CIn_ZEN_gluc" = CInT_ZEN_gluc, 
      "CIn_aZEL" = CInT_aZEL, 
      "CIn_aZEL_gluc" = 0, 
    ) %>% 
    select(!c(starts_with("CDuo"), 
              starts_with("CIle"), 
              starts_with("CJej"), 
              starts_with("CCOL"),
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
    xlab("Days")+
    theme(legend.position = "bottom",
          legend.title = element_text(size = 35),
          legend.text = element_text(size = 35),
          legend.key.size = unit(0.5, "lines"),
          legend.spacing.y = unit(0, "cm"),
          legend.margin = margin(0, 0, 0, 0))
  CBlood.plot
  ggsave(
    filename = here(OUTPUT,"BloodC.png"),
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
    xlab("Days")+
    theme(legend.position = "bottom",
          legend.title = element_text(size = 35),
          legend.text = element_text(size = 35),
          legend.key.size = unit(0.5, "lines"),
          legend.spacing.y = unit(0, "cm"),
          legend.margin = margin(0, 0, 0, 0))
  CTissues.plot
  ggsave(
    filename = here(OUTPUT,"TissueC.png"),
    dpi = 600,
    width = 11,
    height = 7,
    units = "cm"
  )
  
  write_csv(Clean.C, here(OUTPUT, "C.Tissues.csv"))
  
  
  return(list(
    CALC_Parameters = parm.c,
    PBK_OUTPUT = PBK_OUTPUT,
    C.Tissues = Clean.C,
    C.All = C.Pre_long,
    Dosing = Dosingevents,
    Meta_data = list(exp = paste(round(exp, digits = 10),"ug/kg/day", sep = "_"),
                         expAGE = paste(expAGE, "years", sep = "_"),
                         Tstart = paste(round(Tstart/365, digits = 10), "years", sep = "_"),
                         expSTOP = paste(round(expSTOP/365, digits = 10), "years", sep = "_"),
                         Tstop = paste(round(Tstop/365, digits = 10), "years", sep = "_"),
                         expBW = paste(round(expBW, digits = 10), "kg", sep = "_"),
                         sex = paste(sex)),
    Plots = list(MB.plot, CBlood.plot, CTissues.plot) # 
  ))
  
  
}


# With physiological changes with lifestage (depending on age) ----
RUNandOUT_lifestage <- function(exp, 
                      Tinput, tinterval, expSTOP, 
                      expAGE = NULL, expBW = NULL, sex = "M", 
                      Tstart, Tstop, Dt) {
  
  
  source(here("Script", "CALC_Parameters.R"), local = TRUE)
  source(here("Script", "PBK_model.R"), local = TRUE)  
  
  
  ## Add state ####
  A_init <- c(
    
    AOD_ZEN = 0, #exp*expBW, #0, 
    
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
    AFE_aZEL_gluc = 0#,
    
    
  )
  
  
  # Initial input for lifestage changes
  fullPBK_OUTPUT <- NULL
  newTstart <- Tstart
  newA_init <- A_init
  newAGE <- expAGE
  newBW <- expBW
  
  while (newTstart < Tstop) {
    
    yearlyTstop <- min(newTstart + 365, Tstop) 
    
    
    ## Calculare parameters, per newAGE ####
    parms = PARAMS(newAGE, newBW, sex) 
    write.csv(parms, file = here(OUTPUT, "ModelParameters.csv"), row.names = FALSE)
    
    parm.c <- as.numeric(parms$Value)
    names(parm.c) <- parms$Parameter
    
    ODOSEumol <-  exp*parm.c[["BW"]]
    Dosingevents <- data.frame(var = "AOD_ZEN",
                               time =  seq(newTstart,yearlyTstop,by=Dt)) %>% 
      mutate(value = ifelse(time <= expSTOP, ODOSEumol, 0),
             # value = if( time <= expSTOP){ODOSEumol}else{0},
             method = "add")
    
    
    
    ## Run the PBK model, per current age ####
    yearlyPBK_OUTPUT <- PBK_RUN(y = newA_init, 
                                parms = parm.c, 
                                times = seq(newTstart, yearlyTstop, by=Dt),
                                events = Dosingevents)
    
    if(is.null(fullPBK_OUTPUT)){
      fullDosingevents <- Dosingevents
      fullPBK_OUTPUT <- yearlyPBK_OUTPUT 
    } else {
      fullDosingevents <- rbind(fullDosingevents, Dosingevents[-1,])
      fullPBK_OUTPUT <- rbind(fullPBK_OUTPUT, yearlyPBK_OUTPUT[-1,]) 
    }
    
    
    # Update input for next iteration (# Last state values)
    newA_init <- round(unlist(yearlyPBK_OUTPUT[nrow(yearlyPBK_OUTPUT), 2:52]), 10)  
    newTstart <- yearlyTstop
    newAGE <- expAGE + (newTstart/365)
    newBW <- NA
    
  }

  Dosingevents <- fullDosingevents
  PBK_OUTPUT <- fullPBK_OUTPUT
  
  ## Plot Mass Balance/Error ####
  MB.df <- PBK_OUTPUT %>% select(time, ACalculated)
  Dosingevents <- Dosingevents %>% filter(row_number() <= n()-1)
  MB.df$ATotal <- c(0,cumsum(Dosingevents$value))
  MB.df$MB <- MB.df$ATotal - MB.df$ACalculated + 1                                               # rounding significance points
  MB.df$ERROR <- (MB.df$ATotal - MB.df$ACalculated) / MB.df$ACalculated * 100
  MB.df$ERROR <- round(MB.df$ERROR, 10)                                           # rounding significance points
  MB.plot <- ggplot(data = MB.df)+
    geom_line(aes(x = time, y = ERROR, color = "ERROR")) +
    geom_line(aes(x = time, y = MB, color = "MB")) +
    scale_color_manual(values = c("ERROR" = "blue", "MB" = "black"),
                       name = NULL) +
    CP_theme +
    ylab("MB / ERROR")
  MB.plot
  ggsave(
    filename = here(OUTPUT,"MB.png"),
    dpi = 600,
    width = 11,
    height = 7,
    units = "cm"
  )
  
  ## Plot organ concentrations ####
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
    "CCOL" = "Colon",
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
    mutate(Concentration = ifelse(Substance == "ZEN", Concentration*MW_ZEN,     # correcting units
                                  ifelse(Substance == "ZEN_gluc",Concentration*MW_ZEN_gluc,
                                         ifelse(Substance == "aZEL", Concentration*MW_aZEL,
                                                Concentration*MW_aZEL_gluc)
                                  )))
  
  C.Pre_long <- Concentrations
  
  BW <- parm.c[["BW"]]                                         # kg
  
  VSI <- parm.c[["FV_SI"]]*BW 
  VDuo <- VSI*parm.c[["FV_Duo"]]
  VJej <- VSI*parm.c[["FV_Jej"]]
  VIle <- VSI*parm.c[["FV_Ile"]]
  
  VSIL <- parm.c[["FV_SIL"]]*BW # fractional volume of small intestinal lumen
  VDuoL <- VSIL*parm.c[["FV_duol"]] 
  VJejL <- VSIL*parm.c[["FV_jejl"]] 
  VIleL <- VSIL*parm.c[["FV_ilel"]] 
  
  
  Clean.replacements <- c(
    "CIn" = "Intestine",
    "CL" = "Liver", 
    "CF" = "Fat", 
    "CR" = "Rapidly perfused", 
    "CS" = "Slowly perfused",
    "CBile" = "Bile",
    "CB" = "Blood", 
    "CK" = "Kidney",
    "CUrine" = "Urine",
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
      "CInT_aZEL" = ((CIle_aZEL*VIle))/( 
        VDuo + VJej + VIle),
      "CInL_aZEL" = ((CIleL_aZEL*VIleL))/(
        VDuoL + VJejL + VIleL),
      "CInL_aZEL_gluc" = ((CDuoL_aZEL_gluc*VDuoL) + (CIleL_aZEL_gluc*VIleL))/( #(CJejL_aZEL_gluc*VJejL)
        VDuoL + VJejL + VIleL)
    ) %>% 
    mutate(
      "CIn_ZEN" = CInT_ZEN, 
      "CIn_ZEN_gluc" = CInT_ZEN_gluc, 
      "CIn_aZEL" = CInT_aZEL, 
      "CIn_aZEL_gluc" = 0, 
    ) %>% 
    select(!c(starts_with("CDuo"), 
              starts_with("CIle"), 
              starts_with("CJej"), 
              starts_with("CCOL"),
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
    xlab("Days")+
    theme(legend.position = "bottom",
          legend.title = element_text(size = 35),
          legend.text = element_text(size = 35),
          legend.key.size = unit(0.5, "lines"),
          legend.spacing.y = unit(0, "cm"),
          legend.margin = margin(0, 0, 0, 0))
  CBlood.plot
  ggsave(
    filename = here(OUTPUT,"BloodC.png"),
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
    xlab("Days")+
    theme(legend.position = "bottom",
          legend.title = element_text(size = 35),
          legend.text = element_text(size = 35),
          legend.key.size = unit(0.5, "lines"),
          legend.spacing.y = unit(0, "cm"),
          legend.margin = margin(0, 0, 0, 0))
  CTissues.plot
  ggsave(
    filename = here(OUTPUT,"TissueC.png"),
    dpi = 600,
    width = 11,
    height = 7,
    units = "cm"
  )
  
  write_csv(Clean.C, here(OUTPUT, "C.Tissues.csv"))
  
  
  return(list(
    CALC_Parameters = parm.c,
    PBK_OUTPUT = PBK_OUTPUT,
    C.Tissues = Clean.C,
    C.All = C.Pre_long,
    Meta_data = list(exp = paste(round(exp, digits = 10),"ug/kg/day", sep = "_"),
                     expAGE = paste(expAGE, "years", sep = "_"),
                     Tstart = paste(round(Tstart/365, digits = 10), "years", sep = "_"),
                     expSTOP = paste(round(expSTOP/365, digits = 10), "years", sep = "_"),
                     Tstop = paste(round(Tstop/365, digits = 10), "years", sep = "_"),
                     expBW = paste(round(expBW, digits = 10), "kg", sep = "_"),
                     sex = paste(sex)),
    Plots = list(MB.plot, CBlood.plot, CTissues.plot)
  ))
  
}

# Post-run analysis #
# 
# POST.Run <- function(exp,
#                      Tinput, tinterval, expSTOP,
#                      expAGE = NULL, expBW = NULL, sex = "M",
#                      Tstart, Tstop, Dt,
#                      RawData){
# 
#   ## Input ####
#   PBK_OUTPUT <- RawData$PBK_OUTPUT
# 
#   ## Plot Mass Balance/Error ####
#   MB.df <- PBK_OUTPUT %>% select(time, ATotal, ACalculated, MB)                   
#   MB.df$MB <- round(MB.df$MB, 10)                                                 # rounding significance points
#   MB.df$ERROR <- (MB.df$ATotal - MB.df$ACalculated) / MB.df$ACalculated * 100
#   MB.df$ERROR <- round(MB.df$ERROR, 10)                                           # rounding significance points
#   MB.plot <- ggplot(data = MB.df)+
#     geom_line(aes(x = time, y = ERROR, color = "ERROR")) +
#     geom_line(aes(x = time, y = MB, color = "MB")) +
#     scale_color_manual(values = c("ERROR" = "blue", "MB" = "black"),
#                        name = NULL) +
#     CP_theme +
#     ylab("MB / ERROR")
#   MB.plot
#   ggsave(
#     filename = here(OUTPUT,"MB.png"),
#     dpi = 600,
#     width = 11,
#     height = 7,
#     units = "cm"
#   )
# 
#   ## Plot organ concentrations ####
#   MW_ZEN_gluc <- pull(PhysChem[PhysChem$Substance == "Zearalenone-14O-glucuronide", 
#                                "MolecularWeight"])
#   MW_aZEL <- pull(PhysChem[PhysChem$Substance == "a-Zearalenol", 
#                            "MolecularWeight"])
#   MW_aZEL_gluc <- pull(PhysChem[PhysChem$Substance == "a-Zearalenol-glucuronide", 
#                                 "MolecularWeight"])
#   
#   replacements <- c(
#     "CDuoL" = "Duodenum Lumen",
#     "CDuo"  = "Duodenum",
#     "CIleL" = "Ileum Lumen",
#     "CIle" = "Ileum",
#     "CJejL" = "Jejunum Lumen", 
#     "CJej" = "Jejunum", 
#     "CCOL" = "Colon",
#     "CL" = "Liver", 
#     "CF" = "Fat", 
#     "CR" = "Rapidly perfused", 
#     "CS" = "Slowly perfused",
#     "CBile" = "Bile",
#     "CB" = "Blood", 
#     "CK" = "Kidney",
#     "_gluc" = "-gluc",
#     "time" = "Time"
#   )
#   
#   
#   
#   Concentrations <- PBK_OUTPUT %>% 
#     select(c("time", starts_with("C"))) %>% 
#     set_names(~(.) %>%
#                 str_replace_all(replacements)) %>% 
#     pivot_longer( !Time , 
#                   names_to = c("Tissue", "Substance"),
#                   names_sep = c("_"),
#                   values_to = "Concentration") %>% 
#     mutate(Concentration = ifelse(Substance == "ZEN", Concentration*MW_ZEN,     # correcting units
#                                   ifelse(Substance == "ZEN_gluc",Concentration*MW_ZEN_gluc,
#                                          ifelse(Substance == "aZEL", Concentration*MW_aZEL,
#                                                 Concentration*MW_aZEL_gluc)
#                                   )))
#   
#   C.Pre_long <- Concentrations
#   
#   BW <- RawData$CALC_Parameters[["BW"]]                                         # kg
#   
#   VSI <- RawData$CALC_Parameters[["FV_SI"]]*BW 
#   VDuo <- VSI*RawData$CALC_Parameters[["FV_Duo"]]
#   VJej <- VSI*RawData$CALC_Parameters[["FV_Jej"]]
#   VIle <- VSI*RawData$CALC_Parameters[["FV_Ile"]]
#  
#   VSIL <- RawData$CALC_Parameters[["FV_SIL"]]*BW # fractional volume of small intestinal lumen
#   VDuoL <- VSIL*RawData$CALC_Parameters[["FV_duol"]] 
#   VJejL <- VSIL*RawData$CALC_Parameters[["FV_jejl"]] 
#   VIleL <- VSIL*RawData$CALC_Parameters[["FV_ilel"]] 
#   
#   
#   Clean.replacements <- c(
#     "CIn" = "Intestine",
#     "CL" = "Liver", 
#     "CF" = "Fat", 
#     "CR" = "Rapidly perfused", 
#     "CS" = "Slowly perfused",
#     "CBile" = "Bile",
#     "CB" = "Blood", 
#     "CK" = "Kidney",
#     "CUrine" = "Urine",
#     "_gluc" = "-gluc",
#     "time" = "Time"
#   )
#   
#   
#   Clean.C <- PBK_OUTPUT %>% 
#     select(c("time", starts_with("C"))) %>% 
#     mutate("CInT_ZEN" = ((CDuo_ZEN*VDuo) + (CJej_ZEN*VJej) + (CIle_ZEN*VIle))/(
#       VDuo + VJej + VIle),
#       "CInL_ZEN" = ((CDuoL_ZEN*VDuoL) + (CJejL_ZEN*VJejL) + (CIleL_ZEN*VIleL))/(
#         VDuoL + VJejL + VIleL),
#       "CInT_ZEN_gluc" = ((CJej_ZEN_gluc*VJej))/(
#         VDuo + VJej + VIle),
#       "CInL_ZEN_gluc" = ((CDuoL_ZEN_gluc*VDuoL) + (CJejL_ZEN_gluc*VJejL) + (CIleL_ZEN_gluc*VIleL))/(
#         VDuoL + VJejL + VIleL),
#       "CInT_aZEL" = ((CIle_aZEL*VIle))/( 
#         VDuo + VJej + VIle),
#       "CInL_aZEL" = ((CIleL_aZEL*VIleL))/(
#         VDuoL + VJejL + VIleL),
#       "CInL_aZEL_gluc" = ((CDuoL_aZEL_gluc*VDuoL) + (CIleL_aZEL_gluc*VIleL))/( #(CJejL_aZEL_gluc*VJejL)
#         VDuoL + VJejL + VIleL)
#     ) %>% 
#     mutate(
#       "CIn_ZEN" = CInT_ZEN, 
#       "CIn_ZEN_gluc" = CInT_ZEN_gluc, 
#       "CIn_aZEL" = CInT_aZEL, 
#       "CIn_aZEL_gluc" = 0, 
#     ) %>% 
#     select(!c(starts_with("CDuo"), 
#               starts_with("CIle"), 
#               starts_with("CJej"), 
#               starts_with("CCOL"),
#               starts_with("CInT"),
#               starts_with("CInL"))) %>% 
#     set_names(~(.) %>%
#                 str_replace_all(Clean.replacements)) %>% 
#     pivot_longer( !Time , 
#                   names_to = c("Tissue", "Substance"),
#                   names_sep = c("_"),
#                   values_to = "Concentration") %>% 
#     mutate(Concentration = ifelse(Substance == "ZEN", Concentration*MW_ZEN,       # correcting units
#                                   ifelse(Substance == "ZEN_gluc",Concentration*MW_ZEN_gluc,
#                                          ifelse(Substance == "aZEL", Concentration*MW_aZEL,
#                                                 Concentration*MW_aZEL_gluc)
#                                   )))
#   
#   
#   # Plot Blood
#   CBlood.plot <- C.Pre_long %>% filter(Tissue == "Blood") %>%
#     ggplot() + geom_path(aes(x = Time, y = Concentration)) +
#     facet_wrap(~ Substance, scales = "free") +
#     CP_theme +
#     theme(strip.text = element_text(size = 45),
#           axis.text = element_text(size = 40)) +
#     ylab("Blood Concentration (ng/ml)") +
#     xlab("Days")+
#     theme(legend.position = "bottom",
#           legend.title = element_text(size = 35),
#           legend.text = element_text(size = 35),
#           legend.key.size = unit(0.5, "lines"),
#           legend.spacing.y = unit(0, "cm"),
#           legend.margin = margin(0, 0, 0, 0))
#   CBlood.plot
#   ggsave(
#     filename = here(OUTPUT,"BloodC.png"),
#     dpi = 600,
#     width = 11,
#     height = 7,
#     units = "cm"
#   )
#   
#   # Plot all tissues
#   
#   CTissues.plot <- Clean.C %>% 
#     ggplot() + 
#     geom_path(aes(x = Time, y = Concentration, colour = Substance)) +
#     scale_color_brewer(palette = "Dark2") +
#     facet_wrap(~ Tissue, scales = "free") +
#     CP_theme +
#     theme(strip.text = element_text(size = 45),
#           axis.text = element_text(size = 40)) +
#     ylab("Concentration (ng/ml)") +
#     xlab("Days")+
#     theme(legend.position = "bottom",
#           legend.title = element_text(size = 35),
#           legend.text = element_text(size = 35),
#           legend.key.size = unit(0.5, "lines"),
#           legend.spacing.y = unit(0, "cm"),
#           legend.margin = margin(0, 0, 0, 0))
#   CTissues.plot
#   ggsave(
#     filename = here(OUTPUT,"TissueC.png"),
#     dpi = 600,
#     width = 11,
#     height = 7,
#     units = "cm"
#   )
#   
#   write_csv(Clean.C, here(OUTPUT, "C.Tissues.csv"))
#   
# 
#   return(list(
#     ANALYSED_data = list(exp = paste(round(exp, digits = 10),"ug/kg/day", sep = "_"),
#                          expAGE = paste(expAGE, "years", sep = "_"),
#                          Tstart = paste(round(Tstart/365, digits = 10), "years", sep = "_"),
#                          expSTOP = paste(round(expSTOP/365, digits = 10), "years", sep = "_"),
#                          Tstop = paste(round(Tstop/365, digits = 10), "years", sep = "_"),
#                          expBW = paste(round(expBW, digits = 10), "kg", sep = "_"),
#                          sex = paste(sex)#,
#                          # ZEN_BPR = pull(Kp_fus_BPR_df[Kp_fus_BPR_df$Parameter == "ZEN_BPR", "Value"]),
#                          # aZEL_BPR = pull(Kp_fus_BPR_df[Kp_fus_BPR_df$Parameter == "aZEL_BPR", "Value"]),
#                          # ZEN_gluc_BPR = pull(Kp_fus_BPR_df[Kp_fus_BPR_df$Parameter == "ZEN_gluc_BPR", "Value"]),
#                          # aZEL_gluc_BPR =  pull(Kp_fus_BPR_df[Kp_fus_BPR_df$Parameter == "aZEL_gluc_BPR", "Value"])
#                          ),
#     Raw.Data = PBK_OUTPUT,
#     C.Tissues = Clean.C,
#     C.All = C.Pre_long,
#     Plots = list(MB.plot, CBlood.plot, CTissues.plot)
#   ))
# }