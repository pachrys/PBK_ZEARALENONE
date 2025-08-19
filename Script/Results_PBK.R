# # --------------------------------------------------------------------------- #
# # Script for producing the results of the Zearalenone PBK model 
# # By: Chrysanthi Pachoulide
# # Date: 23-06-2025
# # --------------------------------------------------------------------------- #
# rm(list=ls()) 
# 
# # Packages
# library(here)
# library(tidyverse)
# library(deSolve)
# library(PKNCA)
# library(pracma)
# library(showtext)
# library(plotly)
# font_add(family = "Garamond", regular = "GARA.TTF")
# showtext_auto()

CP_theme <- theme_minimal() +
  theme(
    plot.background = element_rect(fill = 'white', color = NA),
    panel.background = element_rect(fill = 'white', color = NA),
    axis.text = element_text(size = 56),
    axis.title = element_text(size = 62),
    plot.title = element_text(size = 65, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 56, hjust = 0.5),
    legend.position = "bottom",
    legend.text = element_text(size = 45),
    legend.title = element_text(size = 50),
    legend.margin = margin(0, 0, 0, 0)
  )
MWZEN <- 318.37
MWZENGLUC <- 494.493
MWaZEL <- 320.4
BP <- 0.67 #assuming the same for all zearalenone, a-zearalenol and zearalenone glucuronide

# load(here(OUTPUT, "RawData_PBK.RData"))
# save(PBK_OUTPUT, file = here(OUTPUT, "RawData_PBK.RData"))
# PBK_OUTPUT <- PBK.out %>% as.data.frame()

C.Pre <- PBK_OUTPUT %>%
  transmute(Time = time,
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

C.Pre_long <- C.Pre %>%
  pivot_longer(names_to = c("Tissue", "Compound"),
               DuodenumT_ZEN:Blood_ZENGLU,
               values_to = "Concentration",
               names_sep = "_") %>%
  mutate(Concentration = ifelse(Compound == "ZEN", Concentration*MWZEN, # correct units
                                ifelse(Compound == "ZENGLU",Concentration*MWZENGLUC,
                                       Concentration*MWaZEL)))

BW <- parm.c$BW  ##{Kg} #  body weight rat
VST <-parm.c$FV_ST*BW  # fraction volume of small intestine
VSIL <- parm.c$FV_SIL*BW # fractional volume of small intestinal lumen
# VCOL <- (parm.c$FV_LIL*BW)-VREL
VSI <- parm.c$FV_SI*BW # fractional volume of small intestinal tissue
VDuoL <- VSIL*0.183 # calculated from Musther
VJejL <- VSIL*0.334 # calculated from Musther
VIleL <- VSIL*0.483 # calculated from Musther
VDuo <- VSI*0.060
VJej <- VSI*0.769
VIle <- VSI*0.171

C.Pre_wide <- C.Pre_long %>%
  pivot_wider(names_from = "Tissue", values_from = "Concentration") %>%
  select(c("Compound","Time", "DuodenumT", "JejunumT", "IleumT", "DuodenumL", "JejunumL", "IleumL", "LIntestineL", "Liver", "Blood")) %>%
  # mutate("IntestineL" = (((DuodenumL*VDuoL) + (JejunumL*VJejL) + (IleumL*VIleL) + (LIntestineL*VCOL))/(VDuoL+VJejL+VIleL+VCOL))) %>%
  # mutate("Intestine" = (((SIntestineT*VSI) + (IntestineL*(VSIL+VCOL+VREL)))/(VSI + VSIL+VCOL+VREL))) %>%
  mutate("Intestine" = ((DuodenumT*VDuo) + (JejunumT*VJej) + (IleumT*VIle))/(# + (DuodenumL*VDuoL) + (JejunumL*VJejL) + (IleumL*VIleL)+ (LIntestineL*VCOL) + (RectumL*VREL)
    VDuo + VJej + VIle)) %>% #L+VCOL+VRE  +VDuoL+VJejL+VIle
  # rename("Intestine" = "SIntestineT") %>%
  select(c("Compound", "Time", "Intestine", "Liver", "Blood")) %>% 
  filter(Time >= 0.03)

# Plot against in vivo ####

## ZEN ##
DataInvivo <- read.csv(here("Input", "InVivoData.csv"), header=1)

Obs <- DataInvivo %>%
  filter(AdminType == Select_AdminType,
           Gender != "Pregnant.F", #pregnancy not included yet
           TestCompound == "ZEN", #other options are ZEN_glu and a_ZEL
           Author == Select_Author,
           SampleCompound == "ZEN") %>% 
  select(c("Author", "Gender","Tissue", "SamplingTime", "Concentration", "Cmin", "Cmax")) %>%
  rename(Time = SamplingTime,
         Sex = Gender) %>%
  na.omit() %>% 
  mutate(Time = as.numeric(Time), 
         Concentration = as.numeric(Concentration),
         Concentration = ifelse(Tissue == "Plasma",
                               Concentration * BP, Concentration), # (ng/ml) blood concentration: changing from plasma concentrations to blood concentrations by multiplying with BP
         Cmin = ifelse(Tissue == "Plasma",
                                Cmin * BP, Cmin),
         Cmax = ifelse(Tissue == "Plasma",
                       Cmax * BP, Cmax),
         Tissue = ifelse(Tissue == "Plasma", "Blood", Tissue))

if(Select_Author == "Videmann"){
  Cobs <- Obs %>%
    filter(Tissue %in% c("Blood", "Intestine", "Liver"))
  Cpre <- C.Pre_wide %>% filter(Compound == "ZEN") %>%
    pivot_longer(Intestine:Blood, names_to = "Tissue", values_to = "Concentration")
}else{
  Cobs <- Obs %>% filter(Tissue == "Blood")
  Cpre <- C.Pre_wide %>% filter(Compound == "ZEN") %>%
  select(Compound, Time, Blood) %>% mutate("Tissue" = "Blood") %>%
  rename("Concentration" = "Blood")
}


# Plots
PredvsObsBlood <- ggplot(data = filter(Cobs, Tissue == "Blood"),
         aes(x = Time, y = Concentration)) +
  geom_ribbon(aes(ymin = Cmin, ymax = Cmax, fill = Sex),
              alpha = 0.3, color = NA) +
  scale_fill_manual(values = c("#008475", "#8631A7"))+
  geom_point(aes(color = Sex),size = 1) +
  scale_color_manual(values = c("#008475", "#8631A7"))+
  geom_path(data = filter(Cpre, Tissue == "Blood"),
            aes(Time, Concentration)) +
  labs(title = "Zearalenone blood concentrations",
       x = "Time (h)",
       y = "Concentration (ng/ml)") +
  guides(color = guide_legend(override.aes = list(fill = NA))) +
  # theme_minimal() +
  CP_theme
PredvsObsBlood
ggsave(
  filename = here(OUTPUT,"PredvsObsBlood.png"),
  dpi = 1000,
  width = 11,
  height = 7,
  units = "cm"
)

if(Select_Author == "Videmann"){

  PredvsObsAll <- ggplot(data = Cobs,
                         aes(x = Time, y = Concentration)) +
    facet_wrap(~Tissue, scales = "free")+
    geom_ribbon(aes(ymin = Cmin, ymax = Cmax, fill = Sex),
                alpha = 0.3, color = NA) +
    scale_fill_manual(values = c("#008475", "#8631A7"))+
    geom_point(aes(color = Sex),size = 1) +
    scale_color_manual(values = c("#008475", "#8631A7"))+
    geom_path(data = Cpre,
              aes(Time, Concentration)) +
    facet_wrap(~Tissue, scales = "free")+
    geom_path(data = Cpre,
              aes(Time, Concentration)) +
    facet_wrap(~Tissue, scales = "free")+
    labs(title = "Zearalenone concentrations",
         x = "Time (h)",
         y = "Concentration (ng/ml)") +
    guides(color = guide_legend(override.aes = list(fill = NA))) +
    # theme_minimal() +
    CP_theme +
    theme(strip.text = element_text(size = 45))
  PredvsObsAll
  ggsave(
    filename = here(OUTPUT,"PredvsObsAll.png"),
    dpi = 1000,
    width = 11,
    height = 7,
    units = "cm"
  )

}else{print("only blood this time")}

## ZEN0GLUC ##

if(Select_Author == "Videmann"){

  Obs.glu <- DataInvivo %>%
    filter(AdminType == Select_AdminType,
           Gender != "Pregnant.F", #pregnancy not included yet
           TestCompound == "ZEN", #other options are ZEN_glu and a_ZEL
           Author == Select_Author,
           SampleCompound == "ZEN_glu")   %>% 
    select(c("Author", "Gender","Tissue", "SamplingTime", "Concentration")) %>% #, "Cmin", "Cmax"
    rename(Time = SamplingTime,
           Sex = Gender) %>%
    na.omit() %>% 
    mutate(Time = as.numeric(Time), 
           Concentration = as.numeric(Concentration),
           Concentration = ifelse(Tissue == "Plasma",
                                  Concentration * BP, Concentration), # (ng/ml) blood concentration: changing from plasma concentrations to blood concentrations by multiplying with BP
           # Cmin = ifelse(Tissue == "Plasma",
           #               Cmin * BP, Cmin),
           # Cmax = ifelse(Tissue == "Plasma",
           #               Cmax * BP, Cmax),
           Tissue = ifelse(Tissue == "Plasma", "Blood", Tissue))
  
  Cobs.glu <- Obs.glu %>%
      filter(Tissue %in% c("Blood", "Intestine", "Liver"))
  Cpre.glu <- C.Pre_wide %>% filter(Compound == "ZENGLU") %>%
      pivot_longer(Intestine:Blood, names_to = "Tissue", values_to = "Concentration")

  PredvsObsBlood.GLU <- ggplot() +
    geom_point(data = filter(Cobs.glu, Tissue == "Blood"),
               aes(Time, Concentration, color = Sex)) +
    geom_path(data = filter(Cpre.glu, Tissue == "Blood"),
              aes(Time, Concentration)) +
    labs(title = "Zearalenone glucuronide blood concentrations",
         x = "Time (h)",
         y = "Concentration (ng/ml)")+
    # theme_minimal() +
    CP_theme + 
    theme(strip.text = element_text(size = 45))
  PredvsObsBlood.GLU
  ggsave(
    filename = here(OUTPUT,"PredvsObsBloodZENGLUC.png"),
    dpi = 1000,
    width = 11,
    height = 7,
    units = "cm"
  )

  PredvsObsAll.GLU <- ggplot() +
    geom_point(data = Cobs.glu,
               aes(Time, Concentration, color = Sex)) +
    facet_wrap(~Tissue, scales = "free")+
    geom_path(data = Cpre.glu,
              aes(Time, Concentration)) +
    facet_wrap(~Tissue, scales = "free")+
    labs(title = "Zearalenone glucuronide concentrations",
         x = "Time (h)",
         y = "Concentration (ng/ml)")+
    # theme_minimal() +
    CP_theme +
    guides(color = guide_legend(
      override.aes = list(shape = c(4, 5),
                          color = c("blueviolet", "darkblue")
      ))
    ) +
    theme(strip.text = element_text(size = 45))
    PredvsObsAll.GLU
  ggsave(
    filename = here(OUTPUT,"PredvsObsAllZENGLU.png"),
    dpi = 1000,
    width = 11,
    height = 7,
    units = "cm"
  )

}else{ print("Only ZEN nothing else ;)")}

# Regression lines ####

Cobs <- Cobs %>% filter(Tissue == "Blood")
Cpre <- Cpre %>% filter(Tissue == "Blood")

Obstimepoints <- Cobs$Time
Predicted <- Cpre %>% filter(Time %in% Obstimepoints)
if(Select_Author == "Videmann"){

  Observed.M <- Cobs %>% filter(Time %in% Predicted$Time) %>% filter(Sex == "Male")
  Observed.F <- Cobs %>% filter(Time %in% Predicted$Time) %>% filter(Sex == "Female")

  C_PredictedObserved.M <- data.frame(Observed = Observed.M$Concentration,
                                            Predicted = Predicted$Concentration)
  CRegression.M <- lm(Predicted~Observed, data = C_PredictedObserved.M) #lm(y~x) (y = dependent, x = independent)
  summary(CRegression.M)

  C_PredictedObserved.F <- data.frame(Observed = Observed.F$Concentration,
                                      Predicted = Predicted$Concentration)
  CRegression.F <- lm(Predicted~Observed, data = C_PredictedObserved.F)
  summary(CRegression.F)


  CRegression_plot <- CRegression.M %>%
    ggplot(aes(x = log(Observed), y = log(Predicted))) +
    # geom_smooth(method='lm', color = "black", se = TRUE)  +
    geom_abline(intercept = 0, slope = 1, linetype = "solid", linewidth = 0.3, color = "grey50") +
    geom_abline(intercept = log(2), slope = 1, linetype = "dashed", linewidth = 0.2, color = "grey50") +
    geom_abline(intercept = log(0.5), slope = 1, linetype = "dashed", linewidth = 0.2, color = "grey50") +
    geom_abline(intercept = log(5), slope = 1, linetype = "dotted", linewidth = 0.1, color = "grey50") +
    geom_abline(intercept = log(1/5), slope = 1, linetype = "dotted", linewidth = 0.1, color = "grey50") +
    geom_point(color = "blueviolet", shape = 4) +
    geom_point(data = CRegression.F, aes(x = log(Observed), y = log(Predicted)),
               color = "darkblue", shape = 5) +
    scale_color_manual(name = "Sex",
      values = c("Male" = "blueviolet", "Female" = "darkblue"),
      labels = c("Male", "Female")
    ) +
    scale_shape_manual(name = "Sex",
      values = c("Male" = 4, "Female" = 5),
      labels = c("Male", "Female")
    ) +
    geom_abline(intercept = 0, slope = 1, linetype = "solid", linewidth = 0.3, color = "grey50") +
    geom_abline(intercept = log(2), slope = 1, linetype = "dashed", linewidth = 0.2, color = "grey50") +
    geom_abline(intercept = log(0.5), slope = 1, linetype = "dashed", linewidth = 0.2, color = "grey50") +
    geom_abline(intercept = log(5), slope = 1, linetype = "dotted", linewidth = 0.1, color = "grey50") +
    geom_abline(intercept = log(1/5), slope = 1, linetype = "dotted", linewidth = 0.1, color = "grey50") +
    CP_theme +
    labs(x = "log 10 (Observed Blood Concentration) (ng/ml)", y = "log10 (Predicted Blood Concentration) (ng/ml)")
  CRegression_plot

  CRegression_plot <- ggplot() +
    geom_abline(intercept = 0, slope = 1, linetype = "solid", linewidth = 0.3, color = "grey50") +
    geom_abline(intercept = log(2), slope = 1, linetype = "dashed", linewidth = 0.2, color = "grey50") +
    geom_abline(intercept = log(0.5), slope = 1, linetype = "dashed", linewidth = 0.2, color = "grey50") +
    geom_abline(intercept = log(5), slope = 1, linetype = "dotted", linewidth = 0.1, color = "grey50") +
    geom_abline(intercept = log(1/5), slope = 1, linetype = "dotted", linewidth = 0.1, color = "grey50") +
    geom_point(data = CRegression.M,
      aes(x = log(Observed), y = log(Predicted), color = "Male", shape = "Male"),
      size = 1) +
    geom_point(
      data = CRegression.F,
      aes(x = log(Observed), y = log(Predicted), color = "Female", shape = "Female"),
      size = 1) +
    scale_color_manual(name = "Sex", values = c("Male" = "blueviolet", "Female" = "darkblue")) +
    scale_shape_manual(name = "Sex", values = c("Male" = 4, "Female" = 5)) +
    CP_theme +
    labs(x = "log 10 (Observed Blood Concentration) (ng/ml)",
         y = "log10 (Predicted Blood Concentration) (ng/ml)") +
    guides(color = guide_legend(override.aes = list(shape = c(4, 5))))
  CRegression_plot
  ggsave(filename = here(OUTPUT, "CRegression_plot.png"),
         plot = CRegression_plot,
         dpi = 1000,
         width = 9,height = 7, units = "cm")

  }else{

  Observed <- Cobs %>% filter(Time %in% Predicted$Time)
  C_PredictedObserved <- data.frame(Observed = Observed$Concentration,
                                      Predicted = Predicted$Concentration)
  CRegression <- lm(Predicted~Observed, data = C_PredictedObserved)
  summary(CRegression)

  CRegression_plot <- CRegression %>%
    ggplot(aes(x = log(Observed), y = log(Predicted))) +
    geom_abline(intercept = 0, slope = 1, linetype = "solid", linewidth = 0.3, color = "grey50") +
    geom_abline(intercept = log(2), slope = 1, linetype = "dashed", linewidth = 0.2, color = "grey50") +
    geom_abline(intercept = log(0.5), slope = 1, linetype = "dashed", linewidth = 0.2, color = "grey50") +
    geom_abline(intercept = log(5), slope = 1, linetype = "dotted", linewidth = 0.1, color = "grey50") +
    geom_abline(intercept = log(1/5), slope = 1, linetype = "dotted", linewidth = 0.1, color = "grey50") +
    geom_point(color = "blueviolet", shape = 4) +
    CP_theme +
    labs(x = "log 10 (Observed Blood Concentration) (ng/ml)", y = "log10 (Predicted Blood Concentration) (ng/ml)")
  CRegression_plot
  ggsave(filename = here(OUTPUT, "CRegression_plot.png"),
         plot = CRegression_plot,
         dpi = 1000,
         width = 9,height = 7, units = "cm")

  }


# Plot model predictions ####

A.Pre <- PBK_OUTPUT %>%
  transmute(Time = time,
            DuodenumT_ZEN = ADuo_ZEN,
            DuodenumL_ZEN = ADuoL_ZEN,
            JejunumT_ZEN = AJej_ZEN,
            JejunumL_ZEN = AJejL_ZEN,
            IleumT_ZEN = AIle_ZEN,
            IleumL_ZEN = AIleL_ZEN,
            LIntestineL_ZEN = ACOL_ZEN,
            Feces_ZEN = AFE_ZEN,
            Liver_ZEN = AL_ZEN,
            Fat_ZEN = AF_ZEN,
            Rapidly_ZEN = AR_ZEN,
            Slowly_ZEN = AS_ZEN,
            Blood_ZEN = AB_ZEN,
            Urine_ZEN = AUR_ZEN,
            DuodenumT_aZEL = 0,
            DuodenumL_aZEL = 0,
            JejunumT_aZEL = AJej_aZEL,
            JejunumL_aZEL = 0,
            IleumT_aZEL = 0,
            IleumL_aZEL = 0,
            LIntestineL_aZEL = ACOL_aZEL,
            Feces_aZEL = AFE_aZEL,
            Liver_aZEL = AL_aZEL,
            Fat_aZEL = AF_aZEL,
            Rapidly_aZEL = AR_aZEL,
            Slowly_aZEL = AS_aZEL,
            Blood_aZEL = AB_aZEL,
            Urine_aZEL = AUR_aZEL,
            DuodenumT_ZENGLU = 0,
            DuodenumL_ZENGLU = ADuoL_ZENGLU,
            JejunumT_ZENGLU = AJej_IC_ZENGLU + AJej_EC_ZENGLU,#AJej_ZENGLU,
            JejunumL_ZENGLU = AJejL_ZENGLU,
            IleumT_ZENGLU = 0,
            IleumL_ZENGLU = AIleL_ZENGLU,
            LIntestineL_ZENGLU = ACOL_ZENGLU,
            Feces_ZENGLU = AFE_ZENGLU,
            Liver_ZENGLU = AL_IC_ZENGLU + AL_EC_ZENGLU, #AL_ZENGLU,
            Bile_ZENGLU = ABile_ZENGLU,
            Fat_ZENGLU = AF_ZENGLU,
            Rapidly_ZENGLU = AR_ZENGLU,
            Slowly_ZENGLU = AS_ZENGLU,
            Blood_ZENGLU = AB_ZENGLU,
            Urine_ZENGLU = AUR_ZENGLU)

A.Pre_long <- A.Pre %>%
  pivot_longer(names_to = c("Tissue", "Compound"),
               DuodenumT_ZEN:Urine_ZENGLU,
               values_to = "Amount",
               names_sep = "_")


Plot.A.Pre <- A.Pre_long %>%
  ggplot() +
  geom_path(aes(x = Time, y = Amount, colour = Compound), linewidth = 0.5) +
  facet_wrap(~Tissue, scales = "free")+
  # theme_minimal() +
  CP_theme +
  theme(strip.text = element_text(size = 45),
        axis.text = element_text(size = 40)) +
  ylab("Tissue Amount (umol)") +
  xlab("Time (hours)")
print(Plot.A.Pre)
ggsave(
  filename = here(OUTPUT,"AmountAll.png"),
  dpi = 1000,
  width = 11,
  height = 7,
  units = "cm"
)

Plot.C.Pre <- C.Pre_long %>%
  ggplot() +
  geom_path(aes(x = Time, y = Concentration, colour = Compound), linewidth = 0.5) +
  facet_wrap(~Tissue, scales = "free")+
  # theme_minimal() +
  CP_theme +
  theme(strip.text = element_text(size = 45),
        axis.text = element_text(size = 40)) +
  ylab("Tissue Concentration (ng/ml)") +
  xlab("Time (hours)")
print(Plot.C.Pre)
ggsave(
  filename = here(OUTPUT,"ConcentrationAll.png"),
  dpi = 1000,
  width = 11,
  height = 7,
  units = "cm"
)

C_Blood <- C.Pre_long %>%
  filter(Tissue == "Blood") %>%
  ggplot()+
  geom_path(aes(x = Time, y = Concentration, colour = Compound), linewidth = 0.5) +
  # theme_minimal()+
  CP_theme +
  ylab("Blood Concentration (ng/ml)") +
  xlab("Time (hours)")
print(C_Blood)
ggsave(
  filename = here(OUTPUT,"C_Blood.png"),
  dpi = 1000,
  width = 11,
  height = 7,
  units = "cm"
)

# Mass Balance ####
MB.df <- PBK_OUTPUT %>% select(time, ATotal, ACalculated, MB)  # change days to years if needed
ATOTAL <- ODOSEumol + IVDOSEumol
MB.df$MB <- (ATOTAL- MB.df$ACalculated) + 1
MB.df$MB <- round(MB.df$MB, 10) # rounding significance points
# MB.df$ERROR <- (ATOTAL - MB.df$ACalculated) / MB.df$ACalculated * 100
MB.df$ERROR <- (MB.df$ATotal - MB.df$ACalculated) / MB.df$ACalculated * 100
MB.df$ERROR <- round(MB.df$ERROR, 10) # rounding significance points
MB_plot <- ggplot(data = MB.df)+
  geom_line(aes(x = time, y = ERROR, color = "ERROR")) +
  geom_line(aes(x = time, y = MB, color = "MB")) +
  scale_color_manual(values = c("ERROR" = "blue", "MB" = "black"),
                     name = NULL) +
  # ylim(0,1) +
  theme_minimal() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12)
  )+
  ylab("MB / ERROR")
print(MB_plot)
ggsave(
  filename = here(OUTPUT,"MB_plot.png"),
  dpi = 1000,
  width = 11,
  height = 7,
  units = "cm"
)