PBK_OUTPUT <- RawData$PBK_OUTPUT

## Plot Mass Balance/Error ####
# MB.df <- PBK_OUTPUT %>% select(time, Ain, Atot, MB)  # change days to years if needed
# MB.df$MB <- round(MB.df$MB, 10) # rounding significance points
# MB.df$ERROR <- (MB.df$Ain - MB.df$Atot) / MB.df$Atot * 100
# MB.df$ERROR <- round(MB.df$ERROR, 10) # rounding significance points
# MB_plot <- ggplot(data = MB.df)+
#   geom_line(aes(x = time, y = ERROR, color = "ERROR")) +
#   geom_line(aes(x = time, y = MB, color = "MB")) +
#   scale_color_manual(values = c("ERROR" = "blue", "MB" = "black"), 
#                      name = NULL) +
#   # ylim(0,1) +
#   CP_theme +
#   # theme(axis.text = element_text(size = 10),
#   #       axis.title = element_text(size = 12)
#   # )+
#   ylab("MB / ERROR")
# MB_plot

## Plot organ concentrations ####
MWZEN <- 318.37
MWZENGLUC <- 494.493
MWaZEL <- 320.4
BP <- 0.67 #assuming the same for all zearalenone, a-zearalenol and zearalenone glucuronide

BW <- RawData$CALC_Parameters$BW  ##{Kg} #  body weight rat
VST <-RawData$CALC_Parameters$FV_ST*BW  # fraction volume of small intestine
VSIL <- RawData$CALC_Parameters$FV_SIL*BW # fractional volume of small intestinal lumen
# VCOL <- (RawData$CALC_Parameters$FV_LIL*BW)-VREL
VSI <- RawData$CALC_Parameters$FV_SI*BW # fractional volume of small intestinal tissue
VDuoL <- VSIL*RawData$CALC_Parameters$FV_duol 
VJejL <- VSIL*RawData$CALC_Parameters$FV_jejl 
VIleL <- VSIL*RawData$CALC_Parameters$FV_ilel 
VCOL <- RawData$CALC_Parameters$FV_COL*BW 

VDuo <- VSI*RawData$CALC_Parameters$FV_Duo
VJej <- VSI*RawData$CALC_Parameters$FV_Jej
VJej_IC <- RawData$CALC_Parameters$FV_I_IC*VJej
VJej_EC <- RawData$CALC_Parameters$FV_I_EC*VJej
VIle <- VSI*RawData$CALC_Parameters$FV_Ile

# Create df with organ concentrations

C.Pre <- PBK_OUTPUT %>%
  transmute(Time = time/365,
            DuodenumT_ZEN = CDuo_ZEN, DuodenumL_ZEN = CDuoL_ZEN,
            JejunumT_ZEN = CJej_ZEN, JejunumL_ZEN = CJejL_ZEN,
            IleumT_ZEN = CIle_ZEN, IleumL_ZEN = CIleL_ZEN,
            LIntestineL_ZEN = CCOL_ZEN,
            Liver_ZEN = CL_ZEN, 
            Fat_ZEN = CF_ZEN, Rapidly_ZEN = CR_ZEN, Slowly_ZEN = CS_ZEN,
            Blood_ZEN = CB_ZEN, DuodenumT_aZEL = 0, DuodenumL_aZEL = 0,
            JejunumT_aZEL = CJej_aZEL, JejunumL_aZEL = 0,
            IleumT_aZEL = 0, IleumL_aZEL = 0,
            LIntestineL_aZEL = 0, Liver_aZEL = CL_aZEL,
            Fat_aZEL = CF_aZEL, Rapidly_aZEL = CR_aZEL, Slowly_aZEL = CS_aZEL,
            Blood_aZEL = CB_aZEL,
            DuodenumT_ZENGLU = 0, DuodenumL_ZENGLU = CDuoL_ZENGLU,
            JejunumT_ZENGLU = CJej_ZENGLU, JejunumL_ZENGLU = CJejL_ZENGLU,
            IleumT_ZENGLU = 0, IleumL_ZENGLU = CIleL_ZENGLU, LIntestineL_ZENGLU = CCOL_ZENGLU, JejunumT_ZENGLU = CJej_ZENGLU,
            Liver_ZENGLU = CL_ZENGLU, 
            Fat_ZENGLU = CF_ZENGLU, Rapidly_ZENGLU = CR_ZENGLU, Slowly_ZENGLU = CS_ZENGLU,
            Blood_ZENGLU = CB_ZENGLU)

# 
# C.Pre_long.in <- C.Pre %>%
#   pivot_longer(names_to = c("Tissue", "Compound"),
#                DuodenumT_ZEN:Blood_ZENGLU,
#                values_to = "Concentration",
#                names_sep = "_") %>%
#   mutate(Concentration = ifelse(Compound == "ZEN", Concentration*MWZEN, # correct units
#                                 ifelse(Compound == "ZENGLU",Concentration*MWZENGLUC,
#                                        Concentration*MWaZEL)))
# 
# 
# 
# C.Pre_wide <- C.Pre_long.in %>%
#   pivot_wider(names_from = "Tissue", values_from = "Concentration") %>%
#   select(c("Compound","Time", "DuodenumT", "JejunumT", "IleumT", "DuodenumL", "JejunumL", "IleumL", "LIntestineL", "Liver", "Fat", "Blood")) %>%
#   # mutate("IntestineL" = (((DuodenumL*VDuoL) + (JejunumL*VJejL) + (IleumL*VIleL) + (LIntestineL*VCOL))/(VDuoL+VJejL+VIleL+VCOL))) %>%
#   # mutate("Intestine" = (((SIntestineT*VSI) + (IntestineL*(VSIL+VCOL+VREL)))/(VSI + VSIL+VCOL+VREL))) %>%
#   mutate("Intestine" = ((DuodenumT*VDuo) + (JejunumT*VJej) + (IleumT*VIle))/(# + (DuodenumL*VDuoL) + (JejunumL*VJejL) + (IleumL*VIleL)+ (LIntestineL*VCOL) + (RectumL*VREL)
#     VDuo + VJej + VIle)) %>% #L+VCOL+VRE  +VDuoL+VJejL+VIle
#   # rename("Intestine" = "SIntestineT") %>%
#   select(c("Compound", "Time", "Intestine", "Liver", "Fat", "Blood")) 
# 
# C.Pre_long <- C.Pre_wide %>% 
#   pivot_longer(Intestine:Blood, values_to = "Concentration", names_to = "Tissue")
# # Plot organ concentrations
# #   Plot.C.Pre <- 
# #     
# # # C.Pre_long %>%
# ggplot(C.Pre_long) +
#   geom_path(aes(x = Time*365, y = Concentration, colour = Tissue)) +
#   facet_wrap(~Compound, scales = "free") +
#   CP_theme +
#   theme(strip.text = element_text(size = 45),
#         axis.text = element_text(size = 40)) +
#   ylab("Tissue Concentration (ng/ml)") +
#   xlab("Days")+
#   theme(legend.position = "bottom",  
#         legend.title = element_text(size = 35),
#         legend.text = element_text(size = 35),
#         legend.key.size = unit(0.5, "lines"),
#         legend.spacing.y = unit(0, "cm"),
#         legend.margin = margin(0, 0, 0, 0))
# # print(Plot.C.Pre)
# ggsave(
#   filename = here(OUTPUT,"ConcentrationAll.png"),
#   dpi = 600,
#   width = 11,
#   height = 7,
#   units = "cm"
# )
# 
# # C_Blood <- 
# C.Pre_long %>% filter(Tissue == "Blood") %>% 
#   ggplot() + geom_path(aes(x = Time*365, y = Concentration)) +
#   facet_wrap(~ Compound, scales = "free") +
#   CP_theme +
#   theme(strip.text = element_text(size = 45),
#         axis.text = element_text(size = 40)) +
#   ylab("Blood Concentration (ng/ml)") +
#   xlab("Days")+
#   theme(legend.position = "bottom",  
#         legend.title = element_text(size = 35),
#         legend.text = element_text(size = 35),
#         legend.key.size = unit(0.5, "lines"),
#         legend.spacing.y = unit(0, "cm"),
#         legend.margin = margin(0, 0, 0, 0))
# 
# ggsave(
#   filename = here(OUTPUT,"BloodC.png"),
#   dpi = 600,
#   width = 11,
#   height = 7,
#   units = "cm"
# )
# 
# 
# C.Pre_long %>% filter(Tissue == "Fat") %>% 
#   ggplot() + geom_path(aes(x = Time*365, y = Concentration)) +
#   facet_wrap(~ Compound, scales = "free") +
#   CP_theme +
#   theme(strip.text = element_text(size = 45),
#         axis.text = element_text(size = 40)) +
#   ylab("Fat Concentration (ng/ml)") +
#   xlab("Days")+
#   theme(legend.position = "bottom",  
#         legend.title = element_text(size = 35),
#         legend.text = element_text(size = 35),
#         legend.key.size = unit(0.5, "lines"),
#         legend.spacing.y = unit(0, "cm"),
#         legend.margin = margin(0, 0, 0, 0))
# 
# ggsave(
#   filename = here(OUTPUT,"FatC.png"),
#   dpi = 600,
#   width = 11,
#   height = 7,
#   units = "cm"
# )
# 
# Test <- C.Pre_long %>% mutate(Time = Time*365) %>% filter(Tissue == "Blood",
#                                                           Compound=="ZEN") %>% 
#   mutate(initialT = Time/365)
# 
# C.Pre_long %>% filter(Tissue == "Intestine") %>% 
#   ggplot() + geom_path(aes(x = Time*365, y = Concentration)) +
#   facet_wrap(~ Compound, scales = "free") +
#   CP_theme +
#   theme(strip.text = element_text(size = 45),
#         axis.text = element_text(size = 40)) +
#   ylab("Intestine Concentration (ng/ml)") +
#   xlab("Days")+
#   theme(legend.position = "bottom",  
#         legend.title = element_text(size = 35),
#         legend.text = element_text(size = 35),
#         legend.key.size = unit(0.5, "lines"),
#         legend.spacing.y = unit(0, "cm"),
#         legend.margin = margin(0, 0, 0, 0))
# 
# ggsave(
#   filename = here(OUTPUT,"IntestineC.png"),
#   dpi = 600,
#   width = 11,
#   height = 7,
#   units = "cm"
# )
# 

C_ZEN <- last(C.Pre$Blood_ZEN)* 1000 #umol to nmol/L* MWZEN # * 1000 #* 1000
print(C_ZEN)
C_aZEL <- last(C.Pre$Blood_aZEL) * 1000 #umol to nmol/L* MWaZEL #* 1000 #* 1000
print(C_aZEL)
C_ZENGluc <- last(C.Pre$Blood_ZENGLU) * 1000 #umol to nmol/L* MWZENGLUC
print(C_ZENGluc)

Cmax_ZEN <- max(C.Pre$Blood_ZEN) * 1000 #umol to nmol/L* MWZEN # * 1000 #* 1000
print(Cmax_ZEN)
Cmax_aZEL <- max(C.Pre$Blood_aZEL) * 1000 #umol to nmol/L* MWaZEL #* 1000 #* 1000
print(Cmax_aZEL)
Cmax_ZENGluc <- max(C.Pre$Blood_ZENGLU) * 1000 #umol to nmol/L* MWZENGLUC
print(C_ZENGluc)
