# --------------------------------------------------------------------------- #
# SCRIPT FOR CALCULATING MODEL RESULTS
# By: Chrysanthi Pachoulide
# Date: 14-04-2025
# --------------------------------------------------------------------------- #

# Without physiological changes with lifestage (depending on age) ####
RUNandOUT <- function(exposure_type, 
                      exp_Oral = NULL, exp_Dermal = NULL, exp, 
                      Tinput, tinterval, expSTOP, 
                      expAGE = NULL, expBW = NULL, sex = "M", 
                      Tstart, Tstop, Dt) {
  
  
  # Load R files 
  source(here("Script", "CALC_Parameters.R"), local = TRUE)
  source(here("Script", "PBK_model.R"), local = TRUE)  
  
  ## Add state ####
  Astate <- switch (exposure_type,
                    "Oral" =  c(
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
                      
                      AB_ZEN = 0, AB_aZEL = 0, AB_ZENGLU = 0, 
                      
                      AUR_ZEN = 0, AUR_aZEL = 0, AUR_ZENGLU = 0, 
                      AF_ZEN = 0, AF_aZEL = 0, AF_ZENGLU = 0, 
                      AR_ZEN = 0, AR_aZEL = 0, AR_ZENGLU = 0, 
                      AS_ZEN = 0, AS_aZEL = 0, AS_ZENGLU = 0, 
                      
                      Ain = 0
                    ) ,
                    "Dermal" =  c(
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
                      
                      AB_ZEN = 0, AB_aZEL = 0, AB_ZENGLU = 0, 
                      
                      AUR_ZEN = 0, AUR_aZEL = 0, AUR_ZENGLU = 0, 
                      AF_ZEN = 0, AF_aZEL = 0, AF_ZENGLU = 0, 
                      AR_ZEN = 0, AR_aZEL = 0, AR_ZENGLU = 0, 
                      AS_ZEN = 0, AS_aZEL = 0, AS_ZENGLU = 0,
                      
                      Ain = 0
                      
                    ),
                    "Oral_Dermal" =   c(
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
                      
                      AB_ZEN = 0, AB_aZEL = 0, AB_ZENGLU = 0, 
                      
                      AUR_ZEN = 0, AUR_aZEL = 0, AUR_ZENGLU = 0, 
                      AF_ZEN = 0, AF_aZEL = 0, AF_ZENGLU = 0, 
                      AR_ZEN = 0, AR_aZEL = 0, AR_ZENGLU = 0, 
                      AS_ZEN = 0, AS_aZEL = 0, AS_ZENGLU = 0, 
                      
                      Ain = 0
                      
                    ),
                    "Inhalation" =  c(
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
                      
                      AB_ZEN = 0, AB_aZEL = 0, AB_ZENGLU = 0, 
                      
                      AUR_ZEN = 0, AUR_aZEL = 0, AUR_ZENGLU = 0, 
                      AF_ZEN = 0, AF_aZEL = 0, AF_ZENGLU = 0, 
                      AR_ZEN = 0, AR_aZEL = 0, AR_ZENGLU = 0, 
                      AS_ZEN = 0, AS_aZEL = 0, AS_ZENGLU = 0, 
                      
                      Ain = 0
                      
                    )
  )
  
  
  # Calculate parameters
  base.parm.c = BASE_PARAMS(expAGE, expBW, sex) # common parameters
  parm.c <- switch(exposure_type,
                   "Oral" = c(base.parm.c, 
                              expSTOP = expSTOP, 
                             expOral = exp, 
                              Tinput = Tinput, 
                              tinterval = tinterval),  
                   "Dermal" = c(DERMAL_PARAMS(expAGE, expBW, sex, base.parm.c), 
                                list(expSTOP = expSTOP, 
                                   expDermal = exp, 
                                     Tinput = Tinput, 
                                     tinterval = tinterval)),
                   "Oral_Dermal" = c(DERMAL_PARAMS(expAGE, expBW, sex = "M", base.parm.c), 
                                     list(expSTOP = expSTOP, 
                                         expOral = exp_Oral,
                                        expDermal = exp_Dermal, 
                                          Tinput = Tinput, 
                                          tinterval = tinterval)),
                   "Inhalation" = c(INHALATION_PARAMS(expAGE, expBW, sex = "M", base.parm.c), 
                                    list(expSTOP = expSTOP, 
                                      expLung = exp, 
                                         Tinput = Tinput, 
                                         tinterval = tinterval))
  )
  # write.csv(parm.c, file = here(OUTPUT, "ModelParameters.csv"), row.names = FALSE)
  
  
  
  ## Run the PBK model ####
  PBK_OUTPUT <- switch(exposure_type,
                        "Oral" = ORAL_PBK_RUN(y = Astate, 
                                              parms = parm.c, 
                                              times = seq(Tstart,Tstop,by=Dt)),
                        "Dermal" = DERMAL_PBK_RUN(y = Astate, 
                                                  parms = parm.c, 
                                                  times = seq(Tstart,Tstop,by=Dt)),
                        "Oral_Dermal" = ORAL_DERMAL_PBK_RUN(y = Astate, 
                                                            parms = parm.c, 
                                                            times = seq(Tstart, Tstop, by=Dt)),
                        "Inhalation" = INHALATION_PBK_RUN(y = Astate, 
                                                          parms = parm.c, 
                                                          times = seq(Tstart,Tstop,by=Dt))
  )

  return(list(
    CALC_Parameters = parm.c,
    PBK_OUTPUT = PBK_OUTPUT
  ))
  
  
}


# With physiological changes with lifestage (depending on age) ####
RUNandOUT_lifestage <- function(exposure_type, 
                      exp_Oral = NULL, exp_Dermal = NULL, exp, 
                      Tinput, tinterval, expSTOP, 
                      expAGE = NULL, expBW = NULL, sex = "M", 
                      Tstart, Tstop, Dt) {
  
  
  source(here("Script", "CALC_Parameters.R"), local = TRUE)
  source(here("Script", "PBK_model.R"), local = TRUE)  
  
  
  ## Add state ####
  Astate <- switch (exposure_type,
                    "Oral" =   c(
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
                      
                      AB_ZEN = 0, AB_aZEL = 0, AB_ZENGLU = 0, 
                      
                      AUR_ZEN = 0, AUR_aZEL = 0, AUR_ZENGLU = 0, 
                      AF_ZEN = 0, AF_aZEL = 0, AF_ZENGLU = 0, 
                      AR_ZEN = 0, AR_aZEL = 0, AR_ZENGLU = 0, 
                      AS_ZEN = 0, AS_aZEL = 0, AS_ZENGLU = 0, 
                      
                      Ain = 0
                      
                    ),
                    "Dermal" =   c(
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
                      
                      AB_ZEN = 0, AB_aZEL = 0, AB_ZENGLU = 0, 
                      
                      AUR_ZEN = 0, AUR_aZEL = 0, AUR_ZENGLU = 0, 
                      AF_ZEN = 0, AF_aZEL = 0, AF_ZENGLU = 0, 
                      AR_ZEN = 0, AR_aZEL = 0, AR_ZENGLU = 0, 
                      AS_ZEN = 0, AS_aZEL = 0, AS_ZENGLU = 0, 
                      
                      Ain = 0
                      
                    ),
                    "Oral_Dermal" =   c(
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
                      
                      AB_ZEN = 0, AB_aZEL = 0, AB_ZENGLU = 0, 
                      
                      AUR_ZEN = 0, AUR_aZEL = 0, AUR_ZENGLU = 0, 
                      AF_ZEN = 0, AF_aZEL = 0, AF_ZENGLU = 0, 
                      AR_ZEN = 0, AR_aZEL = 0, AR_ZENGLU = 0, 
                      AS_ZEN = 0, AS_aZEL = 0, AS_ZENGLU = 0, 
                      
                      Ain = 0
                      
                    ),
                    "Inhalation" =   c(
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
                      
                      AB_ZEN = 0, AB_aZEL = 0, AB_ZENGLU = 0, 
                      
                      AUR_ZEN = 0, AUR_aZEL = 0, AUR_ZENGLU = 0, 
                      AF_ZEN = 0, AF_aZEL = 0, AF_ZENGLU = 0, 
                      AR_ZEN = 0, AR_aZEL = 0, AR_ZENGLU = 0, 
                      AS_ZEN = 0, AS_aZEL = 0, AS_ZENGLU = 0, 
                      
                      Ain = 0
                      
                    )
  )
  
  # Initial input for lifestage changes
  fullPBK_OUTPUT <- NULL
  newTstart <- Tstart
  newAstate <- Astate
  newAGE <- expAGE
  newBW <- expBW
  
  while (newTstart < Tstop) {
    
    yearlyTstop <- min(newTstart + 365, Tstop) 
    
    
    ## Calculare parameters, per newAGE ####
    
    base.parm.c = BASE_PARAMS(newAGE, newBW, sex) # common parameters
    parm.c <- switch(exposure_type,
                     "Oral" = c(base.parm.c, 
                                expSTOP = expSTOP, 
                               expOral = exp, 
                                Tinput = Tinput, 
                                tinterval = tinterval),  
                     "Dermal" = c(DERMAL_PARAMS(newAGE, newBW, sex, base.parm.c), 
                                  list(expSTOP = expSTOP, 
                                     expDermal = exp, 
                                       Tinput = Tinput, 
                                       tinterval = tinterval)),
                     "Oral_Dermal" = c(DERMAL_PARAMS(newAGE, newBW, sex, base.parm.c), 
                                       list(expSTOP = expSTOP, 
                                           expOral = exp_Oral,
                                          expDermal = exp_Dermal, 
                                            Tinput = Tinput, 
                                            tinterval = tinterval)),
                     "Inhalation" = c(INHALATION_PARAMS(newAGE, newBW, sex, base.parm.c), 
                                      list(expSTOP = expSTOP, 
                                        expLung = exp, 
                                           Tinput = Tinput, 
                                           tinterval = tinterval))
    )
    # write.csv(parm.c, file = here(OUTPUT, "ModelParameters.csv"), row.names = FALSE)
    
    
    
    ## Run the PBK model, per current age ####
    yearlyPBK_OUTPUT <- switch(exposure_type,
                               "Oral" = ORAL_PBK_RUN(y = newAstate, 
                                                     parms = parm.c, 
                                                     times = seq(newTstart, yearlyTstop, by=Dt)),
                               "Dermal" = DERMAL_PBK_RUN(y = newAstate, 
                                                         parms = parm.c, 
                                                         times = seq(newTstart, yearlyTstop, by=Dt)),
                               "Oral_Dermal" = ORAL_DERMAL_PBK_RUN(y = newAstate, 
                                                                   parms = parm.c, 
                                                                   times = seq(newTstart, yearlyTstop, by=Dt)),
                               "Inhalation" = INHALATION_PBK_RUN(y = newAstate, 
                                                                 parms = parm.c, 
                                                                 times = seq(newTstart, yearlyTstop, by=Dt))
    )
    
    if(is.null(fullPBK_OUTPUT)){
      fullPBK_OUTPUT <- yearlyPBK_OUTPUT # this is for the first iteration
    } else {
      fullPBK_OUTPUT <- rbind(fullPBK_OUTPUT, yearlyPBK_OUTPUT[-1,]) # binding all iteration results 
    }
    
    
    # Update input for next iteration
    newAstate <- switch(exposure_type, # Last state values
                        "Oral" =round(unlist(yearlyPBK_OUTPUT[nrow(yearlyPBK_OUTPUT), 2:46]), 10),  
                        "Dermal" = unlist(yearlyPBK_OUTPUT[nrow(yearlyPBK_OUTPUT), 2:19]),
                        "Oral_Dermal" = unlist(yearlyPBK_OUTPUT[nrow(yearlyPBK_OUTPUT), 2:20]),
                        "Inhalation" = unlist(yearlyPBK_OUTPUT[nrow(yearlyPBK_OUTPUT), 2:19])
    )    
    newTstart <- yearlyTstop
    newAGE <- expAGE + (newTstart/365)
    newBW <- NA
    
  }

  PBK_OUTPUT <- fullPBK_OUTPUT
  
  return(list(
    CALC_Parameters = parm.c,
    PBK_OUTPUT = PBK_OUTPUT # time is in days
  ))
  
  
}

# # Post-run analysis and plots for 1 individual ####
# 
# POST.Run <- function(exposure_type,
#                      exp_Oral = NULL, exp_Dermal = NULL, exp,
#                      Tinput, tinterval, expSTOP,
#                      expAGE = NULL, expBW = NULL, sex = "M",
#                      Tstart, Tstop, Dt,
#                      RawData){
# 
#   ## Input ####
# #   PBK_OUTPUT <- RawData$PBK_OUTPUT
# # 
# #   ## Plot Mass Balance/Error ####
# #   # MB.df <- PBK_OUTPUT %>% select(time, Ain, Atot, MB)  # change days to years if needed
# #   # MB.df$MB <- round(MB.df$MB, 10) # rounding significance points
# #   # MB.df$ERROR <- (MB.df$Ain - MB.df$Atot) / MB.df$Atot * 100
# #   # MB.df$ERROR <- round(MB.df$ERROR, 10) # rounding significance points
# #   # MB_plot <- ggplot(data = MB.df)+
# #   #   geom_line(aes(x = time, y = ERROR, color = "ERROR")) +
# #   #   geom_line(aes(x = time, y = MB, color = "MB")) +
# #   #   scale_color_manual(values = c("ERROR" = "blue", "MB" = "black"), 
# #   #                      name = NULL) +
# #   #   # ylim(0,1) +
# #   #   CP_theme +
# #   #   # theme(axis.text = element_text(size = 10),
# #   #   #       axis.title = element_text(size = 12)
# #   #   # )+
# #   #   ylab("MB / ERROR")
# #   # MB_plot
# #   
# #   ## Plot organ concentrations ####
# #   MWZEN <- 318.37
# #   MWZENGLUC <- 494.493
# #   MWaZEL <- 320.4
# #   BP <- 0.67 #assuming the same for all zearalenone, a-zearalenol and zearalenone glucuronide
# # 
# #   BW <- RawData$CALC_Parameters$BW  ##{Kg} #  body weight rat
# #   VST <-RawData$CALC_Parameters$FV_ST*BW  # fraction volume of small intestine
# #   VSIL <- RawData$CALC_Parameters$FV_SIL*BW # fractional volume of small intestinal lumen
# #   # VCOL <- (RawData$CALC_Parameters$FV_LIL*BW)-VREL
# #   VSI <- RawData$CALC_Parameters$FV_SI*BW # fractional volume of small intestinal tissue
# #   VDuoL <- VSIL*RawData$CALC_Parameters$FV_duol 
# #   VJejL <- VSIL*RawData$CALC_Parameters$FV_jejl 
# #   VIleL <- VSIL*RawData$CALC_Parameters$FV_ilel 
# #   VCOL <- RawData$CALC_Parameters$FV_COL*BW 
# #   
# #   VDuo <- VSI*RawData$CALC_Parameters$FV_Duo
# #   VJej <- VSI*RawData$CALC_Parameters$FV_Jej
# #   VJej_IC <- RawData$CALC_Parameters$FV_I_IC*VJej
# #   VJej_EC <- RawData$CALC_Parameters$FV_I_EC*VJej
# #   VIle <- VSI*RawData$CALC_Parameters$FV_Ile
# #   
# #   # Create df with organ concentrations
# #   
# #   C.Pre <- PBK_OUTPUT %>%
# #     transmute(Time = time/365,
# #               DuodenumT_ZEN = CDuo_ZEN, DuodenumL_ZEN = CDuoL_ZEN,
# #               JejunumT_ZEN = CJej_ZEN, JejunumL_ZEN = CJejL_ZEN,
# #               IleumT_ZEN = CIle_ZEN, IleumL_ZEN = CIleL_ZEN,
# #               LIntestineL_ZEN = CCOL_ZEN,
# #               Liver_ZEN = CL_ZEN, 
# #               Fat_ZEN = CF_ZEN, Rapidly_ZEN = CR_ZEN, Slowly_ZEN = CS_ZEN,
# #               Blood_ZEN = CB_ZEN, DuodenumT_aZEL = 0, DuodenumL_aZEL = 0,
# #               JejunumT_aZEL = CJej_aZEL, JejunumL_aZEL = 0,
# #               IleumT_aZEL = 0, IleumL_aZEL = 0,
# #               LIntestineL_aZEL = 0, Liver_aZEL = CL_aZEL,
# #               Fat_aZEL = CF_aZEL, Rapidly_aZEL = CR_aZEL, Slowly_aZEL = CS_aZEL,
# #               Blood_aZEL = CB_aZEL,
# #               DuodenumT_ZENGLU = 0, DuodenumL_ZENGLU = CDuoL_ZENGLU,
# #               JejunumT_ZENGLU = CJej_ZENGLU, JejunumL_ZENGLU = CJejL_ZENGLU,
# #               IleumT_ZENGLU = 0, IleumL_ZENGLU = CIleL_ZENGLU, LIntestineL_ZENGLU = CCOL_ZENGLU, JejunumT_ZENGLU = CJej_ZENGLU,
# #               Liver_ZENGLU = CL_ZENGLU, 
# #               Fat_ZENGLU = CF_ZENGLU, Rapidly_ZENGLU = CR_ZENGLU, Slowly_ZENGLU = CS_ZENGLU,
# #               Blood_ZENGLU = CB_ZENGLU)
# #   
# #   
# #   C.Pre_long.in <- C.Pre %>%
# #     pivot_longer(names_to = c("Tissue", "Compound"),
# #                  DuodenumT_ZEN:Blood_ZENGLU,
# #                  values_to = "Concentration",
# #                  names_sep = "_") %>%
# #     mutate(Concentration = ifelse(Compound == "ZEN", Concentration*MWZEN, # correct units
# #                                   ifelse(Compound == "ZENGLU",Concentration*MWZENGLUC,
# #                                          Concentration*MWaZEL)))
# #   
# #   
# #   
# #   C.Pre_wide <- C.Pre_long.in %>%
# #     pivot_wider(names_from = "Tissue", values_from = "Concentration") %>%
# #     select(c("Compound","Time", "DuodenumT", "JejunumT", "IleumT", "DuodenumL", "JejunumL", "IleumL", "LIntestineL", "Liver", "Fat", "Blood")) %>%
# #     # mutate("IntestineL" = (((DuodenumL*VDuoL) + (JejunumL*VJejL) + (IleumL*VIleL) + (LIntestineL*VCOL))/(VDuoL+VJejL+VIleL+VCOL))) %>%
# #     # mutate("Intestine" = (((SIntestineT*VSI) + (IntestineL*(VSIL+VCOL+VREL)))/(VSI + VSIL+VCOL+VREL))) %>%
# #     mutate("Intestine" = ((DuodenumT*VDuo) + (JejunumT*VJej) + (IleumT*VIle))/(# + (DuodenumL*VDuoL) + (JejunumL*VJejL) + (IleumL*VIleL)+ (LIntestineL*VCOL) + (RectumL*VREL)
# #       VDuo + VJej + VIle)) %>% #L+VCOL+VRE  +VDuoL+VJejL+VIle
# #     # rename("Intestine" = "SIntestineT") %>%
# #     select(c("Compound", "Time", "Intestine", "Liver", "Fat", "Blood")) 
# #   
# #   C.Pre_long <- C.Pre_wide %>% 
# #     pivot_longer(Intestine:Blood, values_to = "Concentration", names_to = "Tissue")
# #   # Plot organ concentrations
# # #   Plot.C.Pre <- 
# # #     
# # # # C.Pre_long %>%
# #   ggplot(C.Pre_long) +
# #     geom_path(aes(x = Time*365, y = Concentration, colour = Tissue)) +
# #     facet_wrap(~Compound, scales = "free") +
# #     CP_theme +
# #     theme(strip.text = element_text(size = 45),
# #           axis.text = element_text(size = 40)) +
# #     ylab("Tissue Concentration (ng/ml)") +
# #     xlab("Days")+
# #     theme(legend.position = "bottom",  
# #           legend.title = element_text(size = 35),
# #           legend.text = element_text(size = 35),
# #           legend.key.size = unit(0.5, "lines"),
# #           legend.spacing.y = unit(0, "cm"),
# #           legend.margin = margin(0, 0, 0, 0))
# #   # print(Plot.C.Pre)
# #   ggsave(
# #     filename = here(OUTPUT,"ConcentrationAll.14DExp.30Sim.png"),
# #     dpi = 600,
# #     width = 11,
# #     height = 7,
# #     units = "cm"
# #   )
# #   
# #   # C_Blood <- 
# #   C.Pre_long %>% filter(Tissue == "Blood") %>% 
# #     ggplot() + geom_path(aes(x = Time*365, y = Concentration)) +
# #     facet_wrap(~ Compound, scales = "free") +
# #     CP_theme +
# #     theme(strip.text = element_text(size = 45),
# #           axis.text = element_text(size = 40)) +
# #     ylab("Blood Concentration (ng/ml)") +
# #     xlab("Days")+
# #     theme(legend.position = "bottom",  
# #           legend.title = element_text(size = 35),
# #           legend.text = element_text(size = 35),
# #           legend.key.size = unit(0.5, "lines"),
# #           legend.spacing.y = unit(0, "cm"),
# #           legend.margin = margin(0, 0, 0, 0))
# # 
# #   ggsave(
# #     filename = here(OUTPUT,"BloodC.14DExp.30DSim.png"),
# #     dpi = 600,
# #     width = 11,
# #     height = 7,
# #     units = "cm"
# #   )
# # #   
# #   
# #   message("Post run analysis finished")
# #   
# #   return(list(
# #     ANALYSED_data = data.frame(exposure_type = paste(exposure_type),
# #                                exp = paste(round(exp, digits = 10),"ug/kg/day", sep = "_"),
# #                                expAGE = paste(expAGE, "years", sep = "_"),
# #                                Tstart = paste(round(Tstart/365, digits = 10), "years", sep = "_"),
# #                                expSTOP = paste(round(expSTOP/365, digits = 10), "years", sep = "_"),
# #                                Tstop = paste(round(Tstop/365, digits = 10), "years", sep = "_"),
# #                                expBW = paste(round(expBW, digits = 10), "kg", sep = "_"),
# #                                sex = paste(sex),
# #     OUT_Plots = list(Plot.C.Pre, C_Blood) # could be removed if it's too heavy for R
# #   )))
# # }
# # 
# # 
# # # Post-run analysis and plots per subject ####
# # 
# # Pers.POST.Run <- function(exposure_type, 
# #                      exp_Oral = NULL, exp_Dermal = NULL, exp, 
# #                      Tinput, tinterval, expSTOP, 
# #                      expAGE = NULL, expBW = NULL, sex = "M", 
# #                      Tstart, Tstop, Dt,
# #                      RawData){
# #   
# #   ## Input ####
# #   PBK_OUTPUT <- RawData #$PBK_OUTPUT
# # 
# #   ## Plot Mass Balance/Error ####
# #   # MB.df <- PBK_OUTPUT %>% select(time, Ain, Atot, MB)  # change days to years if needed
# #   # MB.df$MB <- round(MB.df$MB, 10) # rounding significance points
# #   # MB.df$ERROR <- (MB.df$Ain - MB.df$Atot) / MB.df$Atot * 100
# #   # MB.df$ERROR <- round(MB.df$ERROR, 10) # rounding significance points
# #   # MB_plot <- ggplot(data = MB.df)+
# #   #   geom_line(aes(x = time, y = ERROR, color = "ERROR")) +
# #   #   geom_line(aes(x = time, y = MB, color = "MB")) +
# #   #   scale_color_manual(values = c("ERROR" = "blue", "MB" = "black"), 
# #   #                      name = NULL) +
# #   #   # ylim(0,1) +
# #   #   CP_theme +
# #   #   # theme(axis.text = element_text(size = 10),
# #   #   #       axis.title = element_text(size = 12)
# #   #   # )+
# #   #   ylab("MB / ERROR")
# #   # MB_plot
# #   
# #   ## Plot organ concentrations ####
# #   MWZEN <- 318.37
# #   MWZENGLUC <- 494.493
# #   MWaZEL <- 320.4
# #   BP <- 0.67 #assuming the same for all zearalenone, a-zearalenol and zearalenone glucuronide
# #   
# #   BW <- RawData$CALC_Parameters$BW  ##{Kg} #  body weight rat
# #   VST <-RawData$CALC_Parameters$FV_ST*BW  # fraction volume of small intestine
# #   VSIL <- RawData$CALC_Parameters$FV_SIL*BW # fractional volume of small intestinal lumen
# #   # VCOL <- (RawData$CALC_Parameters$FV_LIL*BW)-VREL
# #   VSI <- RawData$CALC_Parameters$FV_SI*BW # fractional volume of small intestinal tissue
# #   VDuoL <- VSIL*RawData$CALC_Parameters$FV_duol 
# #   VJejL <- VSIL*RawData$CALC_Parameters$FV_jejl 
# #   VIleL <- VSIL*RawData$CALC_Parameters$FV_ilel 
# #   VCOL <- RawData$CALC_Parameters$FV_COL*BW 
# #   
# #   VDuo <- VSI*RawData$CALC_Parameters$FV_Duo
# #   VJej <- VSI*RawData$CALC_Parameters$FV_Jej
# #   VJej_IC <- RawData$CALC_Parameters$FV_I_IC*VJej
# #   VJej_EC <- RawData$CALC_Parameters$FV_I_EC*VJej
# #   VIle <- VSI*RawData$CALC_Parameters$FV_Ile
# #   
# #   # Create df with organ concentrations
# #   
# #   C.Pre <- PBK_OUTPUT %>%
# #     
# #     C_organs.df <- switch (exposure_type,
# #                            "Oral" = PBK_OUTPUT %>% 
# #                              transmute(Time = time,
# #                                        DuodenumT_ZEN = CDuo_ZEN, DuodenumL_ZEN = CDuoL_ZEN,
# #                                        JejunumT_ZEN = CJej_ZEN, JejunumL_ZEN = CJejL_ZEN,
# #                                        IleumT_ZEN = CIle_ZEN, IleumL_ZEN = CIleL_ZEN,
# #                                        LIntestineL_ZEN = CCOL_ZEN,
# #                                        Liver_ZEN = CL_ZEN, 
# #                                        Fat_ZEN = CF_ZEN, Rapidly_ZEN = CR_ZEN, Slowly_ZEN = CS_ZEN,
# #                                        Blood_ZEN = CB_ZEN, DuodenumT_aZEL = 0, DuodenumL_aZEL = 0,
# #                                        JejunumT_aZEL = CJej_aZEL, JejunumL_aZEL = 0,
# #                                        IleumT_aZEL = 0, IleumL_aZEL = 0,
# #                                        LIntestineL_aZEL = 0, Liver_aZEL = CL_aZEL,
# #                                        Fat_aZEL = CF_aZEL, Rapidly_aZEL = CR_aZEL, Slowly_aZEL = CS_aZEL,
# #                                        Blood_aZEL = CB_aZEL,
# #                                        DuodenumT_ZENGLU = 0, DuodenumL_ZENGLU = CDuoL_ZENGLU,
# #                                        JejunumT_ZENGLU = CJej_ZENGLU, JejunumL_ZENGLU = CJejL_ZENGLU,
# #                                        IleumT_ZENGLU = 0, IleumL_ZENGLU = CIleL_ZENGLU, LIntestineL_ZENGLU = CCOL_ZENGLU, JejunumT_ZENGLU = CJej_ZENGLU,
# #                                        Liver_ZENGLU = CL_ZENGLU, 
# #                                        Fat_ZENGLU = CF_ZENGLU, Rapidly_ZENGLU = CR_ZENGLU, Slowly_ZENGLU = CS_ZENGLU,
# #                                        Blood_ZENGLU = CB_ZENGLU), 
# #                            "Dermal" = PBK_OUTPUT %>% 
# #                              transmute(Time = time,
# #                                        DuodenumT_ZEN = CDuo_ZEN, DuodenumL_ZEN = CDuoL_ZEN,
# #                                        JejunumT_ZEN = CJej_ZEN, JejunumL_ZEN = CJejL_ZEN,
# #                                        IleumT_ZEN = CIle_ZEN, IleumL_ZEN = CIleL_ZEN,
# #                                        LIntestineL_ZEN = CCOL_ZEN,
# #                                        Liver_ZEN = CL_ZEN, 
# #                                        Fat_ZEN = CF_ZEN, Rapidly_ZEN = CR_ZEN, Slowly_ZEN = CS_ZEN,
# #                                        Blood_ZEN = CB_ZEN, DuodenumT_aZEL = 0, DuodenumL_aZEL = 0,
# #                                        JejunumT_aZEL = CJej_aZEL, JejunumL_aZEL = 0,
# #                                        IleumT_aZEL = 0, IleumL_aZEL = 0,
# #                                        LIntestineL_aZEL = 0, Liver_aZEL = CL_aZEL,
# #                                        Fat_aZEL = CF_aZEL, Rapidly_aZEL = CR_aZEL, Slowly_aZEL = CS_aZEL,
# #                                        Blood_aZEL = CB_aZEL,
# #                                        DuodenumT_ZENGLU = 0, DuodenumL_ZENGLU = CDuoL_ZENGLU,
# #                                        JejunumT_ZENGLU = CJej_ZENGLU, JejunumL_ZENGLU = CJejL_ZENGLU,
# #                                        IleumT_ZENGLU = 0, IleumL_ZENGLU = CIleL_ZENGLU, LIntestineL_ZENGLU = CCOL_ZENGLU, JejunumT_ZENGLU = CJej_ZENGLU,
# #                                        Liver_ZENGLU = CL_ZENGLU, 
# #                                        Fat_ZENGLU = CF_ZENGLU, Rapidly_ZENGLU = CR_ZENGLU, Slowly_ZENGLU = CS_ZENGLU,
# #                                        Blood_ZENGLU = CB_ZENGLU), 
# #                            "Oral_Dermal" = PBK_OUTPUT %>% 
# #                              transmute(Time = time,
# #                                        DuodenumT_ZEN = CDuo_ZEN, DuodenumL_ZEN = CDuoL_ZEN,
# #                                        JejunumT_ZEN = CJej_ZEN, JejunumL_ZEN = CJejL_ZEN,
# #                                        IleumT_ZEN = CIle_ZEN, IleumL_ZEN = CIleL_ZEN,
# #                                        LIntestineL_ZEN = CCOL_ZEN,
# #                                        Liver_ZEN = CL_ZEN, 
# #                                        Fat_ZEN = CF_ZEN, Rapidly_ZEN = CR_ZEN, Slowly_ZEN = CS_ZEN,
# #                                        Blood_ZEN = CB_ZEN, DuodenumT_aZEL = 0, DuodenumL_aZEL = 0,
# #                                        JejunumT_aZEL = CJej_aZEL, JejunumL_aZEL = 0,
# #                                        IleumT_aZEL = 0, IleumL_aZEL = 0,
# #                                        LIntestineL_aZEL = 0, Liver_aZEL = CL_aZEL,
# #                                        Fat_aZEL = CF_aZEL, Rapidly_aZEL = CR_aZEL, Slowly_aZEL = CS_aZEL,
# #                                        Blood_aZEL = CB_aZEL,
# #                                        DuodenumT_ZENGLU = 0, DuodenumL_ZENGLU = CDuoL_ZENGLU,
# #                                        JejunumT_ZENGLU = CJej_ZENGLU, JejunumL_ZENGLU = CJejL_ZENGLU,
# #                                        IleumT_ZENGLU = 0, IleumL_ZENGLU = CIleL_ZENGLU, LIntestineL_ZENGLU = CCOL_ZENGLU, JejunumT_ZENGLU = CJej_ZENGLU,
# #                                        Liver_ZENGLU = CL_ZENGLU, 
# #                                        Fat_ZENGLU = CF_ZENGLU, Rapidly_ZENGLU = CR_ZENGLU, Slowly_ZENGLU = CS_ZENGLU,
# #                                        Blood_ZENGLU = CB_ZENGLU), 
# #                            "Inhalation" = PBK_OUTPUT %>% 
# #                              transmute(Time = time,
# #                                        DuodenumT_ZEN = CDuo_ZEN, DuodenumL_ZEN = CDuoL_ZEN,
# #                                        JejunumT_ZEN = CJej_ZEN, JejunumL_ZEN = CJejL_ZEN,
# #                                        IleumT_ZEN = CIle_ZEN, IleumL_ZEN = CIleL_ZEN,
# #                                        LIntestineL_ZEN = CCOL_ZEN,
# #                                        Liver_ZEN = CL_ZEN, 
# #                                        Fat_ZEN = CF_ZEN, Rapidly_ZEN = CR_ZEN, Slowly_ZEN = CS_ZEN,
# #                                        Blood_ZEN = CB_ZEN, DuodenumT_aZEL = 0, DuodenumL_aZEL = 0,
# #                                        JejunumT_aZEL = CJej_aZEL, JejunumL_aZEL = 0,
# #                                        IleumT_aZEL = 0, IleumL_aZEL = 0,
# #                                        LIntestineL_aZEL = 0, Liver_aZEL = CL_aZEL,
# #                                        Fat_aZEL = CF_aZEL, Rapidly_aZEL = CR_aZEL, Slowly_aZEL = CS_aZEL,
# #                                        Blood_aZEL = CB_aZEL,
# #                                        DuodenumT_ZENGLU = 0, DuodenumL_ZENGLU = CDuoL_ZENGLU,
# #                                        JejunumT_ZENGLU = CJej_ZENGLU, JejunumL_ZENGLU = CJejL_ZENGLU,
# #                                        IleumT_ZENGLU = 0, IleumL_ZENGLU = CIleL_ZENGLU, LIntestineL_ZENGLU = CCOL_ZENGLU, JejunumT_ZENGLU = CJej_ZENGLU,
# #                                        Liver_ZENGLU = CL_ZENGLU, 
# #                                        Fat_ZENGLU = CF_ZENGLU, Rapidly_ZENGLU = CR_ZENGLU, Slowly_ZENGLU = CS_ZENGLU,
# #                                        Blood_ZENGLU = CB_ZENGLU)
# #     ) %>% mutate(time = time/365) #transforming time in years
# #   
# #   C.Pre_long <- C.Pre %>%
# #     pivot_longer(names_to = c("Tissue", "Compound"),
# #                  DuodenumT_ZEN:Blood_ZENGLU,
# #                  values_to = "Concentration",
# #                  names_sep = "_") %>%
# #     mutate(Concentration = ifelse(Compound == "ZEN", Concentration*MWZEN, # correct units
# #                                   ifelse(Compound == "ZENGLU",Concentration*MWZENGLUC,
# #                                          Concentration*MWaZEL)))
# #   
# #   # Plot organ concentrations
# #   Plot.C.Pre <- C.Pre_long %>%
# #     ggplot() +
# #     geom_path(aes(x = Time, y = Concentration, colour = Compound), linewidth = 0.5) +
# #     facet_wrap(~Tissue, scales = "free")+
# #     # theme_minimal() +
# #     CP_theme +
# #     theme(strip.text = element_text(size = 45),
# #           axis.text = element_text(size = 40)) +
# #     ylab("Tissue Concentration (ng/ml)") +
# #     xlab("Time (hours)")
# #   print(Plot.C.Pre)
# #   ggsave(
# #     filename = here(OUTPUT,"ConcentrationAll.png"),
# #     dpi = 600,
# #     width = 11,
# #     height = 7,
# #     units = "cm"
# #   )
# #   
# #   C_Blood <- C.Pre_long %>%
# #     filter(Tissue == "Blood") %>%
# #     ggplot()+
# #     geom_path(aes(x = Time, y = Concentration, colour = Compound), linewidth = 0.5) +
# #     theme_minimal()+
# #     # CP_theme +
# #     ylab("Blood Concentration (ng/ml)") +
# #     xlab("Time (hours)")
# #   print(C_Blood)
# #   ggsave(
# #     filename = here(OUTPUT,"C_Blood.png"),
# #     dpi = 600,
# #     width = 11,
# #     height = 7,
# #     units = "cm"
# #   )
# #   
# #   
#   message("Post run analysis finished")
# 
#   return(list(
#     ANALYSED_data = data.frame(exposure_type = paste(exposure_type),
#                                exp = paste(round(exp, digits = 10),"ug/kg/day", sep = "_"),
#                                expAGE = paste(expAGE, "years", sep = "_"),
#                                Tstart = paste(round(Tstart/365, digits = 10), "years", sep = "_"),
#                                expSTOP = paste(round(expSTOP/365, digits = 10), "years", sep = "_"),
#                                Tstop = paste(round(Tstop/365, digits = 10), "years", sep = "_"),
#                                expBW = paste(round(expBW, digits = 10), "kg", sep = "_"),
#                                sex = paste(sex),,
#     OUT_Plots = list(OUT_Plots = list(Plot.C.Pre, C_Blood) # could be removed if it's too heavy for R
#                      ) # could be removed if it's too heavy for R
#                                )))
# }
# 
# 
# # Post-run analysis and plots per population ####
# 
# Pop.POST.Run <- function(Input,
#                          RawData,
#                          Pers.POST.RUN_RESULTS){
#   
#   
#   ## Input ####
#   
#   PBK_OUT <- RawData$PBK_OUTPUT # PBK model results per subject
#   ANALYSED_data <- lapply(Pers.POST.RUN_RESULTS, function(x) x$ANALYSED_data)
# 
#   MWZEN <- 318.37
#   MWZENGLUC <- 494.493
#   MWaZEL <- 320.4
#   BP <- 0.67 #assuming the same for all zearalenone, a-zearalenol and zearalenone glucuronide
#   
#   Parameters <- RawData$CALC_Parameters
#   
#   BW <- RawData$CALC_Parameters$BW  ##{Kg} #  body weight rat
#   VST <-RawData$CALC_Parameters$FV_ST*BW  # fraction volume of small intestine
#   VSIL <- RawData$CALC_Parameters$FV_SIL*BW # fractional volume of small intestinal lumen
#   # VCOL <- (RawData$CALC_Parameters$FV_LIL*BW)-VREL
#   VSI <- RawData$CALC_Parameters$FV_SI*BW # fractional volume of small intestinal tissue
#   VDuoL <- VSIL*RawData$CALC_Parameters$FV_duol 
#   VJejL <- VSIL*RawData$CALC_Parameters$FV_jejl 
#   VIleL <- VSIL*RawData$CALC_Parameters$FV_ilel 
#   VCOL <- RawData$CALC_Parameters$FV_COL*BW 
#   
#   VDuo <- VSI*RawData$CALC_Parameters$FV_Duo
#   VJej <- VSI*RawData$CALC_Parameters$FV_Jej
#   VJej_IC <- RawData$CALC_Parameters$FV_I_IC*VJej
#   VJej_EC <- RawData$CALC_Parameters$FV_I_EC*VJej
#   VIle <- VSI*RawData$CALC_Parameters$FV_Ile
#   
#   
#   BW <- sapply(Parameters, function(x) x$BW)
#   VST <- sapply(Parameters, function(x) x$VST)*BW
#   VSIL <- sapply(Parameters, function(x) x$FV_SIL)*BW
#   VSI <- sapply(Parameters, function(x) x$FV_SI)*BW
#   VDuoL <- sapply(Parameters, function(x) x$FV_duol)
#   VJejL <- sapply(Parameters, function(x) x$FV_jejl)
#   VIleL <- sapply(Parameters, function(x) x$FV_ilel)
#   VCOL <- sapply(Parameters, function(x) x$VCOL)*BW
#   VDuo <- sapply(Parameters, function(x) x$FV_Duo)*VSI
#   VJej <- sapply(Parameters, function(x) x$FV_Jej)*VSI
#   VJej_IC <- sapply(Parameters, function(x) x$FV_I_IC)*VJej
#   VJej_EC <- sapply(Parameters, function(x) x$FV_I_EC)*VJej
#   VIle <- sapply(Parameters, function(x) x$FV_Ile)*VSI
#   
#   
#   
#   PredictedObserved <- data.frame(
#     Idcode = Input$Idcode,
#     exp = as.numeric(Input$exp),
#     exp_Oral = as.numeric(Input$exp_Oral),
#     exp_Dermal = as.numeric(Input$exp_Dermal),
#     expSTOP = as.numeric(Input$expSTOP), # this is the time of the stop of exposure in days
#     CP_observed = as.numeric(Input$CP_measured), # measured plasma PFOA concentration
#     samplingT = as.numeric(Input$samplingT), # time at which the plasma concentration was measured
#     HL_observed = as.numeric(Input$HL_observed), # half life reported in the HBM study
#     sex = Input$sex,
#     expAGE = as.numeric(Input$expAGE),
#     BW = BW,
#     GFR = GFR
#     ) 
#   PredictedObserved <- PredictedObserved %>% 
#     mutate(log_exp = log10(exp),
#            log_exp_Oral = log10(exp_Oral),
#            log_exp_Dermal = log10(exp_Dermal))
#   
#   # Add predicted concentration at the sampling time 
#   PredictedObserved <- PredictedObserved %>% 
#     mutate(
#       CP_predicted = map2_dbl( 
#         PBK_OUT,samplingT, ~ {
#           idx <- which.min(abs(.x$time - .y))
#           .x$CP[idx]
#         }
#       )
#     ) 
#   
#   # Add predicted half life
#   PredictedObserved <- PredictedObserved %>% mutate(
#     Idcode = seq_along(ANALYSED_data),
#     HalfLife = sapply(ANALYSED_data, function(x) x$HalfLife)) %>%
#     separate(col = HalfLife, into = c("HL_predicted", "unit"), sep = "_") %>%
#     mutate(HL_predicted = as.numeric(HL_predicted), 
#            HL_predicted = round(HL_predicted,1))
#   
#   ## Plots ####
#   
#   PredictedObserved_long <- PredictedObserved %>% 
#     select(exp, exp_Oral, exp_Dermal, log_exp, log_exp_Oral, log_exp_Dermal,
#            CP_observed, HL_observed) %>%
#     pivot_longer(exp:HL_observed, 
#                  names_to = "Variable", 
#                  values_to = "Value")
#   
#   # Create faceted histograms
#   Histograms_variables <- ggplot(PredictedObserved_long, aes(x = Value)) +
#     geom_histogram() +
#     facet_wrap(~ Variable, scales = "free") +
#     CP_theme +
#     theme(axis.text = element_text(size = 5),
#           axis.title = element_text(size = 7),
#           plot.margin = margin(0, 0, 0, 0, "cm") 
#     )
#   Histograms_variables
#   
#   ### Exposure estimate vs plasma concentration ####
#   
#   if(all(PredictedObserved$CP_observed >0) &&
#      !any(is.na(PredictedObserved$CP_observed))) {
#     
#     Plot_exp_vs_O_CP <- PredictedObserved %>% 
#       ggplot(aes(exp, CP_observed)) +
#       geom_point(color = "black", size = 0.5) +
#       CP_theme +
#       labs(title="Exposure vs observed plasma concentration",
#            x="\n Exposure (\u03BCg/kg bw/dayL)", 
#            y="Observed concentration (\u03BCg/L)\n") 
#     
#     Plot_exp_vs_Pr_CP <- PredictedObserved %>% 
#       ggplot(aes(exp, CP_predicted)) +
#       geom_point(color = "black", size = 0.5) +
#       CP_theme +
#       labs(title="Exposure vs predicted plasma concentration",
#            x="\n Exposure (\u03BCg/kg bw/dayL)", 
#            y="Predicted concentration (\u03BCg/L)\n") 
#     
#     Plot_exp_vs_CP <- Plot_exp_vs_O_CP | Plot_exp_vs_Pr_CP
#     
#   } else{ 
#     
#     Plot_exp_vs_P_CP <- PredictedObserved %>% 
#     ggplot(aes(exp, CP_predicted)) +
#     geom_point(color = "black", size = 0.5) +
#     CP_theme +
#     labs(title="Exposure vs predicted plasma concentration",
#          x="\n Exposure (\u03BCg/kg bw/dayL)", 
#          y="Predicted concentration (\u03BCg/L)\n") 
#     
#     Plot_exp_vs_CP <- Plot_exp_vs_P_CP
#     
#     }
#   
#   ### Exposure estimate vs half life ####
#   
#   if(all(PredictedObserved$HL_observed >0) &&
#      !any(is.na(PredictedObserved$HL_observed))) {
#     
#     Plot_exp_vs_O_HL <- PredictedObserved %>% 
#       ggplot(aes(exp, HL_observed)) +
#       geom_point(color = "black", size = 1) +
#       CP_theme +
#       labs(title="Exposure vs observed half life",
#            x="\n Exposure (\u03BCg/kg bw/dayL)", 
#            y="Predicted half life (years)") 
#     
#     Plot_exp_vs_Pr_HL <- PredictedObserved %>% 
#       ggplot(aes(exp, HL_predicted)) +
#       geom_point(color = "black", size = 1) +
#       CP_theme +
#       labs(title="Exposure vs predicted half life",
#            x="\n Exposure (\u03BCg/kg bw/dayL)", 
#            y="Predicted half life (years)") 
#     
#     Plot_exp_vs_HL <- Plot_exp_vs_O_HL | Plot_exp_vs_Pr_HL
#     
#   } else{ 
#     
#     Plot_exp_vs_P_HL <- PredictedObserved %>% 
#       ggplot(aes(exp, HL_predicted)) +
#       geom_point(color = "black", size = 1) +
#       CP_theme +
#       labs(title="Exposure vs predicted half life",
#            x="\n Exposure (\u03BCg/kg bw/dayL)", 
#            y="Predicted half life (years)") 
#     
#     Plot_exp_vs_HL <- Plot_exp_vs_P_HL
#     
#   }
#   
#   ### Exposure age vs predictions ####
#   
#   Plot_age_vs_P_CP <- PredictedObserved %>% 
#     ggplot(aes(expAGE, CP_predicted)) +
#     geom_point(color = "black", size = 0.5) +
#     CP_theme +
#     labs(title="Age at the start of exposure vs predicted plasma concentration",
#          x="\n Exposure age (years)", 
#          y="Predicted concentration (\u03BCg/L)\n") 
#   
#   Plot_age_vs_P_HL <- PredictedObserved %>% 
#     ggplot(aes(expAGE, HL_predicted)) +
#     geom_point(color = "black", size = 1) +
#     CP_theme +
#     labs(title="Age at the start of exposure vs predicted half life",
#          x="\n Exposure age (years)", 
#          y="Predicted half life (years)") 
#   
#   Plot_age <- Plot_age_vs_P_CP | Plot_age_vs_P_HL
#   
#   ### GFR  vs predictions ####
#   
#   Plot_GFR_vs_P_CP <- PredictedObserved %>% 
#     ggplot() +
#     geom_point(aes(GFR, CP_predicted), color = "black", size = 1) +
#     # geom_point(aes(GFR.ff, CP_predicted), color = "blue", size = 1) +
#     CP_theme +
#     labs(title="Glomerular filtration rate vs predicted plasma concentration",
#          x="\n GFR (L/d)", 
#          y="Predicted concentration (\u03BCg/L)\n") 
#   
#   Plot_GFR_vs_P_HL <- PredictedObserved %>% 
#     ggplot() +
#     geom_point(aes(GFR, HL_predicted), color = "black", size = 1) +
#     # geom_point(aes(GFR.ff, HL_predicted), color = "blue", size = 1) +
#     CP_theme +
#     labs(title="Glomerular filtration rate vs predicted half life",
#          x="\n GFR (L/d)", 
#          y="Predicted half life (years)") 
# 
#   Plot_GFR <- Plot_GFR_vs_P_CP | Plot_GFR_vs_P_HL
#   
#   ### Plot predicted half lives over those reported in literature ####
#   
#   HL_literature <- Lit.HalfLifes %>%
#     filter(species == "human", chemical == "pfoa", parameter == "HalfLife") %>%
#     select(c(value_average, n, sex)) %>%
#     rename(HL_observed = value_average) %>%
#     mutate(HL_observed = as.numeric(HL_observed),  # years
#            n = as.numeric(n)) %>% 
#     filter(sex %in% c("F", "M")) %>% 
#     mutate(value = case_when(sex == "F" ~ 0.75,  
#                              sex == "M" ~ 1.5),
#            Origin = "Observed") 
#   
#   HL_predicted <- data.frame(
#     HL_predicted = PredictedObserved$HL_predicted, 
#     sex = PredictedObserved$sex) %>% 
#     mutate(value = case_when(sex == "F" ~ 0.75,  
#                              sex == "M" ~ 1.5),
#            Origin = "Predicted")
#   
#   HL_violin_plot <- 
#     ggplot() +
#     geom_violin(data = filter(HL_literature, sex == "M"), 
#                 aes(x = 1.5, y = HL_observed), 
#                 fill = "grey89", color = NA, width = 0.5, trim = FALSE) +
#     geom_point(data = filter(HL_literature, sex == "M"),
#                aes(x = 1.5, y = HL_observed, size = n),
#                color = "grey70", alpha = 0.5, position = position_jitter(width = 0.05)) +
#     geom_point(data = filter(HL_predicted, sex == "M"),
#                aes(x = 1.5, y = HL_predicted),
#                shape = 18, size = 5,  alpha = 0.9, position = position_jitter(width = 0.01)) +
#     
#     geom_violin(data = filter(HL_literature, sex == "F"), 
#                 aes(x = 0.75, y = HL_observed), 
#                 fill = "grey89", color = NA, width = 0.5, trim = FALSE) +
#     geom_point(data = filter(HL_literature, sex == "F"),
#                aes(x = 0.75, y = HL_observed, size = n),
#                color = "grey70", alpha = 0.5, position = position_jitter(width = 0.05)) +
#     geom_point(data = filter(HL_predicted, sex == "F"),
#                aes(x = 0.75, y = HL_predicted),
#                shape = 18, size = 5,  alpha = 0.9, position = position_jitter(width = 0.01)) +
#     scale_x_continuous(breaks = c(0.75, 1.5),       
#                        labels = c("Female", "Male")) + 
#     scale_size_continuous(range = c(1, 5)) +
#     
#     labs(title = "Predicted half lifes plotted over the range of observed",
#          x = "", y = "Half life (years)") +
#     guides(size = guide_legend(title = "HBM sample size")) +
#     CP_theme +
#     theme(legend.position = "right")
#   HL_violin_plot
#   
#   ### Plot predicted Vs Observed CP ####
#   
#   if(all(PredictedObserved$CP_observed > 0) &&
#      !any(is.na(PredictedObserved$CP_observed > 0))){
#     ### Observed vs Predicted plasma concentrations ####
#     CP_Regression <- lm(CP_predicted~CP_observed, data=PredictedObserved) #lm(y~x) (y is the dependent variable)
#     summary(CP_Regression)
#     
#     CP_Regression_plot <- CP_Regression %>%
#       ggplot(aes(log(CP_observed), log(CP_predicted))) +
#       # geom_smooth(method='lm', color = "black", se = TRUE) +
#       geom_abline(intercept = 0, slope = 1, linetype = "solid", linewidth = 0.3, color = "grey50") +  
#       geom_abline(intercept = log(2), slope = 1, linetype = "dashed", linewidth = 0.2, color = "grey50") + #2 fold
#       geom_abline(intercept = log(0.5), slope = 1, linetype = "dashed", linewidth = 0.2, color = "grey50") + #2 fold
#       geom_abline(intercept = log(5), slope = 1, linetype = "dotted", linewidth = 0.1, color = "grey50") +  #5 fold
#       geom_abline(intercept = log(1/5), slope = 1, linetype = "dotted", linewidth = 0.1, color = "grey50") + #5 fold
#       geom_point(color = "black", size = 1) +
#       CP_theme +
#       labs(title="Observed vs predicted plasma concentrations",
#            x="\n log10 (Observed concentration) (\u03BCg/L)", 
#            y=" log10 (Predicted concentration) (\u03BCg/L)\n") 
#     CP_Regression_plot
#   } else{ CP_Regression_plot <- "No regression plot as no observed data"
#     message("Predicted Vs Observed Plasma concentrations cannot be plotted as observed plasma concentrations have not been provided")} 
#   
#   ### Plot predicted Vs Observed HL ####
#   
#   if(all(PredictedObserved$HL_observed >0) &&
#      !any(is.na(PredictedObserved$HL_observed))) {
#     ### Observed vs Predicted half lives ####
#     HL_Regression <- lm(HL_predicted~HL_observed, data=PredictedObserved)
#     summary(HL_Regression)
#     
#     HL_Regression_plot <- HL_Regression %>%
#       ggplot(aes(log10(HL_observed), log10(HL_predicted))) +
#       # geom_smooth(method='lm', color = "black", se = TRUE) +
#       geom_abline(intercept = 0, slope = 1, linetype = "solid", linewidth = 0.3, color = "grey50") +  
#       geom_abline(intercept = log10(2), slope = 1, linetype = "dashed", linewidth = 0.2, color = "grey50") +  #2 fold
#       geom_abline(intercept = log10(0.5), slope = 1, linetype = "dashed", linewidth = 0.2, color = "grey50") + #2 fold
#       geom_abline(intercept = log10(5), slope = 1, linetype = "dotted", linewidth = 0.1, color = "grey50") +  #5 fold
#       geom_abline(intercept = log10(1/5), slope = 1, linetype = "dotted", linewidth = 0.1, color = "grey50") + #5 fold
#       geom_point(color = "black", size = 1) +
#       CP_theme +
#       labs(title="Observed vs predicted half lives",
#            x = "log10 (Observed Half life) (years)", 
#            y = "log10 (Predicted Half life) (years)")
#     HL_Regression_plot
#   } else{ HL_Regression_plot <- "No regression plot as no observed data"
#   message("Predicted Vs Observed half lives cannot be plotted as observed half lives have not been provided")}
#   
#   message("Population post run analysis finished")
#   
#   PredictedObserved.df <- PredictedObserved %>% select(c(CP_observed, CP_predicted, HL_observed, HL_predicted))
#   
#   
#   return(list(Histograms_variables,
#               Plot_exp_vs_CP,
#               Plot_exp_vs_HL,
#               Plot_age,
#               Plot_GFR,
#               CP_Regression_plot,
#               HL_Regression_plot,
#               HL_violin_plot,
#               PredictedObserved.df
#               ))
#   
# 
#   }
