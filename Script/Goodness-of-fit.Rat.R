# --------------------------------------------------------------------------- #
# Script for evaluating the goodness-of-fit the results of the Zearalenone PBK model
# By: Chrysanthi Pachoulide
# Date: 23-06-2025
# --------------------------------------------------------------------------- #

rm(list=ls())

# Packages
library(here)
library(tidyverse)
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
load(here(OUTPUT, "Results_Rat.RData"))
Select_AdminType <- Results.Rat$Metadata.Run.Rat$Select_AdminType
Select_Author <- Results.Rat$Metadata.Run.Rat$Select_Author
BP <- Results.Rat$Metadata.Run.Rat$BP

# load(here(OUTPUT, "RawData_PBK.Rat.RData"))
Predicted.data <- Results.Rat$C.Tissues.Rat
Observed.data <- read.csv(here("Input", "ObservedData.Rat.csv"), header=1)

Obs <- Observed.data %>%
  filter(AdminType == Select_AdminType,
           Sex != "Pregnant.F", #pregnancy not included
           TestCompound == "ZEN", #other options are ZEN_glu and a_ZEL
           Author == Select_Author,
           SampleCompound == "ZEN") %>% 
  select(c("Author", "Sex","Tissue", "SamplingTime", "Concentration", "Cmin", "Cmax")) %>%
  rename(Time = SamplingTime) %>%
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
  Cpre <- Predicted.data %>% filter(Compound == "ZEN", 
                                    Tissue %in% c("Blood", "Intestine", "Liver")) 
}else{
  Cobs <- Obs %>% filter(Tissue == "Blood")
  Cpre <- Predicted.data %>% filter(Compound == "ZEN", 
                                    Tissue == "Blood")%>% 
    select(Compound, Time, Concentration) #%>% rename("Blood" = "Concentratio")
}


# Plots
if(Select_Author == "Videmann"){
  
  ggplot(
    data = filter(Cobs, Tissue == "Blood"),
    aes(x = Time, y = Concentration)) +
    geom_ribbon(aes(ymin = Cmin, ymax = Cmax, fill = Sex),
                alpha = 0.3, color = NA) +
    scale_fill_manual(values = c("#008475", "#8631A7"))+
    geom_point(aes(color = Sex),size = 1) +
    scale_color_manual(values = c("#008475", "#8631A7")) +
    geom_path(data = filter(Cpre, Tissue == "Blood"),
              aes(Time, Concentration)) +
    labs(title = "Zearalenone blood concentrations",
         x = "Time (h)",
         y = "Concentration (ng/ml)") +
    guides(color = guide_legend(override.aes = list(fill = NA))) +
    CP_theme
  
}else{
  
  ggplot(
  data = Cobs, #filter(Cobs, Tissue == "Blood"),
  aes(x = Time, y = Concentration)) +
  geom_ribbon(aes(ymin = Cmin, ymax = Cmax, fill = Sex),
              alpha = 0.3, color = NA) +
  scale_fill_manual(values = c("#008475", "#8631A7"))+
  geom_point(aes(color = Sex),size = 1) +
  scale_color_manual(values = c("#008475", "#8631A7")) +
  geom_path(data = Cpre, #filter(Cpre, Tissue == "Blood"),
            aes(Time, Concentration)) +
  labs(title = "Zearalenone blood concentrations",
       x = "Time (h)",
       y = "Concentration (ng/ml)") +
  guides(color = guide_legend(override.aes = list(fill = NA))) +
  CP_theme
  
  }

ggsave(
  filename = here(OUTPUT,"PredvsObsBlood.png"),
  dpi = 600,
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
    dpi = 600,
    width = 11,
    height = 7,
    units = "cm"
  )

}else{print("only blood this time")}

## ZEN0GLUC ##

if(Select_Author == "Videmann"){

  Obs.glu <- Observed.data %>%
    filter(AdminType == Select_AdminType,
           Sex != "Pregnant.F", #pregnancy not included yet
           TestCompound == "ZEN", #other options are ZEN_glu and a_ZEL
           Author == Select_Author,
           SampleCompound == "ZEN_glu")   %>% 
    select(c("Author", "Sex","Tissue", "SamplingTime", "Concentration")) %>% #, "Cmin", "Cmax"
    rename(Time = SamplingTime) %>%
    na.omit() %>% 
    mutate(Time = as.numeric(Time), 
           Concentration = as.numeric(Concentration),
           Concentration = ifelse(Tissue == "Plasma",
                                  Concentration * BP, Concentration), # (ng/ml) blood concentration: changing from plasma concentrations to blood concentrations by multiplying with BP
           Tissue = ifelse(Tissue == "Plasma", "Blood", Tissue))
  
  Cobs.glu <- Obs.glu %>%
      filter(Tissue %in% c("Blood", "Intestine", "Liver"))
  Cpre.glu <- Predicted.data %>% filter(Compound == "ZENGLU", 
                                        Tissue %in% c("Blood", "Intestine", "Liver")) 
  
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
    dpi = 600,
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
    dpi = 600,
    width = 11,
    height = 7,
    units = "cm"
  )

}else{ print("Only ZEN nothing else ;)")}

# Regression lines ####
# 
# Cobs <- Cobs %>% filter(Tissue == "Blood")
# Cpre <- Cpre %>% filter(Tissue == "Blood")

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
    ggplot(aes(x = log10(Observed), y = log10(Predicted))) +
    # geom_smooth(method='lm', color = "black", se = TRUE)  +
    geom_abline(intercept = 0, slope = 1, linetype = "solid", linewidth = 0.3, color = "grey50") +
    geom_abline(intercept = log10(2), slope = 1, linetype = "dashed", linewidth = 0.2, color = "grey50") +
    geom_abline(intercept = log10(0.5), slope = 1, linetype = "dashed", linewidth = 0.2, color = "grey50") +
    geom_abline(intercept = log10(10), slope = 1, linetype = "dotted", linewidth = 0.2, color = "grey50") +
    geom_abline(intercept = log10(1/10), slope = 1, linetype = "dotted", linewidth = 0.2, color = "grey50") +
    geom_point(color = "blueviolet", shape = 4) +
    geom_point(data = CRegression.F, aes(x = log10(Observed), y = log10(Predicted)),
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
    geom_abline(intercept = log10(2), slope = 1, linetype = "dashed", linewidth = 0.2, color = "grey50") +
    geom_abline(intercept = log10(0.5), slope = 1, linetype = "dashed", linewidth = 0.2, color = "grey50") +
    geom_abline(intercept = log10(10), slope = 1, linetype = "dotted", linewidth = 0.2, color = "grey50") +
    geom_abline(intercept = log10(1/10), slope = 1, linetype = "dotted", linewidth = 0.2, color = "grey50") +
    CP_theme +
    labs(x = "log10 10 (Observed Blood Concentration) (ng/ml)", y = "log1010 (Predicted Blood Concentration) (ng/ml)")
  CRegression_plot

  CRegression_plot <- ggplot() +
    geom_abline(intercept = 0, slope = 1, linetype = "solid", linewidth = 0.3, color = "grey50") +
    geom_abline(intercept = log10(2), slope = 1, linetype = "dashed", linewidth = 0.2, color = "grey50") +
    geom_abline(intercept = log10(0.5), slope = 1, linetype = "dashed", linewidth = 0.2, color = "grey50") +
    geom_abline(intercept = log10(10), slope = 1, linetype = "dotted", linewidth = 0.2, color = "grey50") +
    geom_abline(intercept = log10(1/10), slope = 1, linetype = "dotted", linewidth = 0.2, color = "grey50") +
    geom_point(data = CRegression.M,
      aes(x = log10(Observed), y = log10(Predicted), color = "Male", shape = "Male"),
      size = 1) +
    geom_point(
      data = CRegression.F,
      aes(x = log10(Observed), y = log10(Predicted), color = "Female", shape = "Female"),
      size = 1) +
    scale_color_manual(name = "Sex", values = c("Male" = "blueviolet", "Female" = "darkblue")) +
    scale_shape_manual(name = "Sex", values = c("Male" = 4, "Female" = 5)) +
    CP_theme +
    labs(x = "log10 10 (Observed Blood Concentration) (ng/ml)",
         y = "log1010 (Predicted Blood Concentration) (ng/ml)") +
    guides(color = guide_legend(override.aes = list(shape = c(4, 5))))
  CRegression_plot
  ggsave(filename = here(OUTPUT, "CRegression_plot.png"),
         plot = CRegression_plot,
         dpi = 600,
         width = 9,height = 7, units = "cm")

  }else{

  Observed <- Cobs %>% filter(Time %in% Predicted$Time)
  C_PredictedObserved <- data.frame(Observed = Observed$Concentration,
                                      Predicted = Predicted$Concentration)
  CRegression <- lm(Predicted~Observed, data = C_PredictedObserved)
  summary(CRegression)

  CRegression_plot <- CRegression %>%
    ggplot(aes(x = log10(Observed), y = log10(Predicted))) +
    geom_abline(intercept = 0, slope = 1, linetype = "solid", linewidth = 0.3, color = "grey50") +
    geom_abline(intercept = log10(2), slope = 1, linetype = "dashed", linewidth = 0.2, color = "grey50") +
    geom_abline(intercept = log10(0.5), slope = 1, linetype = "dashed", linewidth = 0.2, color = "grey50") +
    geom_abline(intercept = log10(10), slope = 1, linetype = "dotted", linewidth = 0.2, color = "grey50") +
    geom_abline(intercept = log10(1/10), slope = 1, linetype = "dotted", linewidth = 0.2, color = "grey50") +
    geom_point(color = "blueviolet", shape = 4)  +
    CP_theme +
    labs(x = expression(paste(bold("Observed (log"[10] ~ "(ng/mL))"))), 
         y = expression(paste(bold("Predicted (log"[10] ~ "(ng/mL))")))
    ) #+
    # coord_cartesian(xlim = c(0.5, 1.5), ylim = c(0.5, 1.5))
  CRegression_plot
  ggsave(filename = here(OUTPUT, "CRegression_plot.png"),
         plot = CRegression_plot,
         dpi = 600,
         width = 9,height = 7, units = "cm")

  }

# 
# # Plot model predictions ####
# 
# A.Pre <- PBK_OUTPUT %>%
#   transmute(Time = time,
#             DuodenumT_ZEN = ADuo_ZEN,
#             DuodenumL_ZEN = ADuoL_ZEN,
#             JejunumT_ZEN = AJej_ZEN,
#             JejunumL_ZEN = AJejL_ZEN,
#             IleumT_ZEN = AIle_ZEN,
#             IleumL_ZEN = AIleL_ZEN,
#             LIntestineL_ZEN = ACOL_ZEN,
#             Feces_ZEN = AFE_ZEN,
#             Liver_ZEN = AL_ZEN,
#             Fat_ZEN = AF_ZEN,
#             Rapidly_ZEN = AR_ZEN,
#             Slowly_ZEN = AS_ZEN,
#             Blood_ZEN = AB_ZEN,
#             Urine_ZEN = AUR_ZEN,
#             DuodenumT_aZEL = 0,
#             DuodenumL_aZEL = 0,
#             JejunumT_aZEL = AJej_aZEL,
#             JejunumL_aZEL = 0,
#             IleumT_aZEL = 0,
#             IleumL_aZEL = 0,
#             LIntestineL_aZEL = ACOL_aZEL,
#             Feces_aZEL = AFE_aZEL,
#             Liver_aZEL = AL_aZEL,
#             Fat_aZEL = AF_aZEL,
#             Rapidly_aZEL = AR_aZEL,
#             Slowly_aZEL = AS_aZEL,
#             Blood_aZEL = AB_aZEL,
#             Urine_aZEL = AUR_aZEL,
#             DuodenumT_ZENGLU = 0,
#             DuodenumL_ZENGLU = ADuoL_ZENGLU,
#             JejunumT_ZENGLU = AJej_IC_ZENGLU + AJej_EC_ZENGLU,#AJej_ZENGLU,
#             JejunumL_ZENGLU = AJejL_ZENGLU,
#             IleumT_ZENGLU = 0,
#             IleumL_ZENGLU = AIleL_ZENGLU,
#             LIntestineL_ZENGLU = ACOL_ZENGLU,
#             Feces_ZENGLU = AFE_ZENGLU,
#             Liver_ZENGLU = AL_IC_ZENGLU + AL_EC_ZENGLU, #AL_ZENGLU,
#             Bile_ZENGLU = ABile_ZENGLU,
#             Fat_ZENGLU = AF_ZENGLU,
#             Rapidly_ZENGLU = AR_ZENGLU,
#             Slowly_ZENGLU = AS_ZENGLU,
#             Blood_ZENGLU = AB_ZENGLU,
#             Urine_ZENGLU = AUR_ZENGLU)
# 
# A.Pre_long <- A.Pre %>%
#   pivot_longer(names_to = c("Tissue", "Compound"),
#                DuodenumT_ZEN:Urine_ZENGLU,
#                values_to = "Amount",
#                names_sep = "_")
# 
# 
# Plot.A.Pre <- A.Pre_long %>%
#   ggplot() +
#   geom_path(aes(x = Time, y = Amount, colour = Compound), linewidth = 0.5) +
#   facet_wrap(~Tissue, scales = "free")+
#   # theme_minimal() +
#   CP_theme +
#   theme(strip.text = element_text(size = 45),
#         axis.text = element_text(size = 40)) +
#   ylab("Tissue Amount (umol)") +
#   xlab("Time (hours)")
# print(Plot.A.Pre)
# ggsave(
#   filename = here(OUTPUT,"AmountAll.png"),
#   dpi = 600,
#   width = 11,
#   height = 7,
#   units = "cm"
# )
# 
# Plot.C.Pre <- C.Pre_long %>%
#   ggplot() +
#   geom_path(aes(x = Time, y = Concentration, colour = Compound), linewidth = 0.5) +
#   facet_wrap(~Tissue, scales = "free")+
#   theme_minimal() +
#   # CP_theme +
#   # theme(strip.text = element_text(size = 45),
#   #       axis.text = element_text(size = 40)) +
#   ylab("Tissue Concentration (ng/ml)") +
#   xlab("Time (hours)")
# print(Plot.C.Pre)
# ggsave(
#   filename = here(OUTPUT,"ConcentrationAll.png"),
#   dpi = 600,
#   width = 11,
#   height = 7,
#   units = "cm"
# )
# 
# C_Blood <- C.Pre_long %>%
#   filter(Tissue == "Blood") %>%
#   ggplot()+
#   geom_path(aes(x = Time, y = Concentration, colour = Compound), linewidth = 0.5) +
#   theme_minimal()+
#   # CP_theme +
#   ylab("Blood Concentration (ng/ml)") +
#   xlab("Time (hours)")
# print(C_Blood)
# ggsave(
#   filename = here(OUTPUT,"C_Blood.png"),
#   dpi = 600,
#   width = 11,
#   height = 7,
#   units = "cm"
# )
# 
# 
# 
# C_ZEN <- max(C.Pre$Blood_ZEN)/MWZEN * 1000 #* 1000
# C_aZEL <- max(C.Pre$Blood_aZEL)/MWZEN * 1000 #* 1000
# 
# 
# # Mass Balance ####
# MB.df <- PBK_OUTPUT %>% select(time, ATotal, ACalculated, MB)  # change days to years if needed
# ATOTAL <- ODOSEumol + IVDOSEumol
# MB.df$MB <- (ATOTAL- MB.df$ACalculated) + 1
# MB.df$MB <- round(MB.df$MB, 10) # rounding significance points
# # MB.df$ERROR <- (ATOTAL - MB.df$ACalculated) / MB.df$ACalculated * 100
# MB.df$ERROR <- (MB.df$ATotal - MB.df$ACalculated) / MB.df$ACalculated * 100
# MB.df$ERROR <- round(MB.df$ERROR, 10) # rounding significance points
# MB_plot <- ggplot(data = MB.df)+
#   geom_line(aes(x = time, y = ERROR, color = "ERROR")) +
#   geom_line(aes(x = time, y = MB, color = "MB")) +
#   scale_color_manual(values = c("ERROR" = "blue", "MB" = "black"),
#                      name = NULL) +
#   # ylim(0,1) +
#   theme_minimal() +
#   theme(axis.text = element_text(size = 10),
#         axis.title = element_text(size = 12)
#   )+
#   ylab("MB / ERROR")
# print(MB_plot)
# ggsave(
#   filename = here(OUTPUT,"MB_plot.png"),
#   dpi = 600,
#   width = 11,
#   height = 7,
#   units = "cm"
# )