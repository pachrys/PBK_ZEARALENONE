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

# Input ---- 
load(here(OUTPUT, "Results.Rat.RData"))
Select_AdminType <- Results.Rat$Metadata.Run.Rat$Select_AdminType
Select_Author <- Results.Rat$Metadata.Run.Rat$Select_Author

ZEN_BPR = Results.Rat$Metadata.Run.Rat$ZEN_BPR[[1]]
aZEL_BPR = Results.Rat$Metadata.Run.Rat$aZEL_BPR[[1]]
ZEN_gluc_BPR = Results.Rat$Metadata.Run.Rat$ZEN_gluc_BPR[[1]]
aZEL_gluc_BPR = Results.Rat$Metadata.Run.Rat$aZEL_gluc_BPR[[1]]


# load(here(OUTPUT, "RawData_PBK.Rat.RData"))
Predicted.data <- Results.Rat$C.Tissues.Rat %>% filter(Time >= 0.04)
Observed.data <- read.csv(here("Input", "ObservedData.Rat.csv"), header=1)

Obs <- Observed.data %>%
  filter(AdminType == Select_AdminType,
           Sex != "Pregnant.F", #pregnancy not included
           TestCompound == "ZEN", #other options are ZEN_glu and a_ZEL
           Author == Select_Author,
           SampleCompound == "ZEN") %>% 
  select(c("Author", "Sex",
           "Tissue", "SamplingTime", 
           "Concentration","Concentration_cumulative", 
           "Cmin", "Cmax")) %>%
  rename(Time = SamplingTime) %>%
  # na.omit() %>% 
  mutate(Time = as.numeric(Time), 
         Concentration = as.numeric(Concentration),
         Concentration = ifelse(Tissue == "Plasma",
                               Concentration * ZEN_BPR, Concentration), # (ng/ml) blood concentration: changing from plasma concentrations to blood concentrations by multiplying with BPR
         Cmin = ifelse(Tissue == "Plasma",
                                Cmin *ZEN_BPR, Cmin),
         Cmax = ifelse(Tissue == "Plasma",
                       Cmax * ZEN_BPR, Cmax),
         Tissue = ifelse(Tissue == "Plasma", "Blood", Tissue))

if(Select_Author == "Videmann"){
  Cobs <- Obs %>%
    filter(Tissue %in% c("Blood", "Intestine", "Liver"))
  Cpre <- Predicted.data %>% filter(Substance == "ZEN", 
                                    Tissue %in% c("Blood", "Intestine", "Liver")) 
}else{
  Cobs <- Obs %>% filter(Tissue == "Blood")
  Cpre <- Predicted.data %>% filter(Substance == "ZEN", 
                                    Tissue == "Blood")%>% 
    select(Substance, Time, Concentration) #%>% rename("Blood" = "Concentratio")
}


# Plots ---- 

## ZEN ----

### Blood ----

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

### All Tissues ----

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

### Urine ----


if(Select_Author == "Videmann"){
  
  Cobs_urine <- Obs %>%
    filter(Tissue == "Urine")
  Cpre_urine <- Predicted.data %>% filter(Substance == "ZEN", 
                                    Tissue == "Urine") 
  
  ggplot(
    data = filter(Cobs_urine),
    aes(x = Time, y = Concentration_cumulative)) +
    geom_point(aes(color = Sex),size = 1) +
    scale_color_manual(values = c("#008475", "#8631A7")) +
    geom_path(data = filter(Cpre_urine, Tissue == "Urine"),
              aes(Time, Concentration)) +
    labs(title = "Zearalenone Urine concentrations",
         x = "Time (h)",
         y = "Concentration (ng/ml)") +
    guides(color = guide_legend(override.aes = list(fill = NA))) +
    CP_theme
  
}else{
  
  ggplot(
    data = Cobs, #filter(Cobs, Tissue == "Urine"),
    aes(x = Time, y = Concentration)) +
    geom_ribbon(aes(ymin = Cmin, ymax = Cmax, fill = Sex),
                alpha = 0.3, color = NA) +
    scale_fill_manual(values = c("#008475", "#8631A7"))+
    geom_point(aes(color = Sex),size = 1) +
    scale_color_manual(values = c("#008475", "#8631A7")) +
    geom_path(data = Cpre, #filter(Cpre, Tissue == "Urine"),
              aes(Time, Concentration)) +
    labs(title = "Zearalenone Urine concentrations",
         x = "Time (h)",
         y = "Concentration (ng/ml)") +
    guides(color = guide_legend(override.aes = list(fill = NA))) +
    CP_theme
  
}

ggsave(
  filename = here(OUTPUT,"PredvsObsUrine.png"),
  dpi = 600,
  width = 11,
  height = 7,
  units = "cm"
)




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
                                  Concentration * ZEN_gluc_BPR, Concentration), # (ng/ml) blood concentration: changing from plasma concentrations to blood concentrations by multiplying with BP
           Tissue = ifelse(Tissue == "Plasma", "Blood", Tissue))
  
  Cobs.glu <- Obs.glu %>%
      filter(Tissue %in% c("Blood", "Intestine", "Liver"))
  Cpre.glu <- Predicted.data %>% filter(Substance == "ZEN-gluc", 
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


## a-ZEL ##

if(Select_Author == "Videmann"){
  
  Obs.aZEL <- Observed.data %>%
    filter(AdminType == Select_AdminType,
           Sex != "Pregnant.F", #pregnancy not included yet
           TestCompound == "ZEN", #other options are ZEN_glu and a_ZEL
           Author == Select_Author,
           SampleCompound == "a_ZEL")   %>% 
    select(c("Author", "Sex","Tissue", "SamplingTime", "Concentration")) %>% #, "Cmin", "Cmax"
    rename(Time = SamplingTime) %>%
    na.omit() %>% 
    mutate(Time = as.numeric(Time), 
           Concentration = as.numeric(Concentration),
           Concentration = ifelse(Tissue == "Plasma",
                                  Concentration * aZEL_BPR, Concentration), # (ng/ml) blood concentration: changing from plasma concentrations to blood concentrations by multiplying with BP
           Tissue = ifelse(Tissue == "Plasma", "Blood", Tissue))
  
  Cobs.aZEL <- Obs.aZEL %>%
    filter(Tissue %in% c("Blood", "Intestine", "Liver"))
  Cpre.aZEL <- Predicted.data %>% filter(Substance == "aZEL", 
                                        Tissue %in% c("Blood", "Intestine", "Liver")) 
  
  PredvsObsBlood.aZEL <- ggplot() +
    geom_point(data = filter(Cobs.aZEL, Tissue == "Blood"),
               aes(Time, Concentration, color = Sex)) +
    geom_path(data = filter(Cpre.aZEL, Tissue == "Blood"),
              aes(Time, Concentration)) +
    labs(title = "a-Zearalenol blood concentrations",
         x = "Time (h)",
         y = "Concentration (ng/ml)")+
    # theme_minimal() +
    CP_theme + 
    theme(strip.text = element_text(size = 45))
  PredvsObsBlood.aZEL
  ggsave(
    filename = here(OUTPUT,"PredvsObsBloodaZEL.png"),
    dpi = 600,
    width = 11,
    height = 7,
    units = "cm"
  )
  
  PredvsObsAll.aZEL <- ggplot() +
    geom_point(data = Cobs.aZEL,
               aes(Time, Concentration, color = Sex)) +
    facet_wrap(~Tissue, scales = "free")+
    geom_path(data = Cpre.glu,
              aes(Time, Concentration)) +
    facet_wrap(~Tissue, scales = "free")+
    labs(title = "a-Zearalenol concentrations",
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
  PredvsObsAll.aZEL
  ggsave(
    filename = here(OUTPUT,"PredvsObsAllaZEL.png"),
    dpi = 600,
    width = 11,
    height = 7,
    units = "cm"
  )
  
}else{ print("Only ZEN nothing else ;)")}



## a-ZEL glucuronide##

if(Select_Author == "Videmann"){
  
  Obs.aZEL_gluc <- Observed.data %>%
    filter(AdminType == Select_AdminType,
           Sex != "Pregnant.F", #pregnancy not included yet
           TestCompound == "ZEN", #other options are ZEN_glu and a_ZEL
           Author == Select_Author,
           SampleCompound == "a_ZEL_glu")   %>% 
    select(c("Author", "Sex","Tissue", "SamplingTime", "Concentration")) %>% #, "Cmin", "Cmax"
    rename(Time = SamplingTime) %>%
    na.omit() %>% 
    mutate(Time = as.numeric(Time), 
           Concentration = as.numeric(Concentration),
           Concentration = ifelse(Tissue == "Plasma",
                                  Concentration * aZEL_gluc_BPR, Concentration), # (ng/ml) blood concentration: changing from plasma concentrations to blood concentrations by multiplying with BP
           Tissue = ifelse(Tissue == "Plasma", "Blood", Tissue))
  
  Cobs.aZEL_gluc <- Obs.aZEL_gluc %>%
    filter(Tissue %in% c("Blood", "Intestine", "Liver"))
  Cpre.aZEL_gluc <- Predicted.data %>% filter(Substance == "aZEL-gluc", 
                                         Tissue %in% c("Blood", "Intestine", "Liver")) 
  
  PredvsObsBlood.aZEL_gluc <- ggplot() +
    geom_point(data = filter(Cobs.aZEL_gluc, Tissue == "Blood"),
               aes(Time, Concentration, color = Sex)) +
    geom_path(data = filter(Cpre.aZEL_gluc, Tissue == "Blood"),
              aes(Time, Concentration)) +
    labs(title = "a-Zearalenol glucuronide blood concentrations",
         x = "Time (h)",
         y = "Concentration (ng/ml)")+
    # theme_minimal() +
    CP_theme + 
    theme(strip.text = element_text(size = 45))
  PredvsObsBlood.aZEL_gluc
  ggsave(
    filename = here(OUTPUT,"PredvsObsBloodaZEL_gluc.png"),
    dpi = 600,
    width = 11,
    height = 7,
    units = "cm"
  )
  
  PredvsObsAll.aZEL_gluc <- ggplot() +
    geom_point(data = Cobs.aZEL_gluc,
               aes(Time, Concentration, color = Sex)) +
    facet_wrap(~Tissue, scales = "free")+
    geom_path(data = Cpre.aZEL_gluc,
              aes(Time, Concentration)) +
    facet_wrap(~Tissue, scales = "free")+
    labs(title = "a-Zearalenol glucuronide concentrations",
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
  PredvsObsAll.aZEL_gluc
  ggsave(
    filename = here(OUTPUT,"PredvsObsAllaZEL_gluc.png"),
    dpi = 600,
    width = 11,
    height = 7,
    units = "cm"
  )
  
}else{ print("Only ZEN nothing else ;)")}




# Regression lines ----
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
