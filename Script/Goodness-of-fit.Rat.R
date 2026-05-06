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
library(readxl)
library(patchwork)
font_add(family = "Garamond", regular = "GARA.TTF")
showtext_auto()

CP_theme <- theme_minimal() +
  theme(
    # plot.background = element_rect(fill = 'white', color = NA),
    # panel.background = element_rect(fill = 'white', color = NA),
    axis.title = element_text(size = 24, face = "bold"),
    axis.text = element_text(size = 20),
    plot.title = element_text(hjust = 0, size = 28, face = "bold"),
    plot.margin = margin(0.3, 0.3, 0.3, 0.3, "cm"), 
    axis.line = element_line(linewidth = 0.1),
    panel.grid = element_line(linewidth = 0.3),
    plot.caption = element_text(hjust = 0, size = 20, face= "italic"),
    plot.title.position = "plot", 
    plot.caption.position = "plot",
    legend.position = "none"
  ) 

# Set output storage directory

OUTPUT <- here("Output", format(Sys.Date(), "%Y-%m-%d")) #, 
# format(Sys.time(), "%H-%M-%S"))
dir.create(OUTPUT, recursive = TRUE)

# Input ---- 
load(here(OUTPUT, "PO_Videmann_Results.Rat.RData"))
Select_AdminType <- Results.Rat$Metadata.Run.Rat$Select_AdminType
# Select_AdminType <- "IV.PERF"
Select_Author <- Results.Rat$Metadata.Run.Rat$Select_Author

ZEN_BPR = Results.Rat$Metadata.Run.Rat$ZEN_BPR[[1]]
aZEL_BPR = Results.Rat$Metadata.Run.Rat$aZEL_BPR[[1]]
ZEN_gluc_BPR = Results.Rat$Metadata.Run.Rat$ZEN_gluc_BPR[[1]]
aZEL_gluc_BPR = Results.Rat$Metadata.Run.Rat$aZEL_gluc_BPR[[1]]
IVDOSEmg = Results.Rat$Metadata.Run.Rat$IVDOSEmg

# load(here(OUTPUT, "RawData_PBK.Rat.RData"))
Predicted.data <- Results.Rat$C.Tissues.Rat %>% filter(Time >= 0.01)#001)
Observed.data <- read_excel("Input/GoodnessOfFit.Rat.xlsx", 
                            col_types = c("text", "text", "text", 
                                          "numeric", "numeric", "text", "text", 
                                          "text", "numeric", "numeric", "numeric", "numeric", 
                                          "text", "text", "text", "numeric", 
                                          "text", "numeric", "numeric", "numeric", 
                                          "numeric", "numeric", "numeric", 
                                          "numeric", "numeric", "numeric", 
                                          "numeric", "numeric", "numeric", 
                                          "text"))
  
Predicted.urine <- Results.Rat$Raw.Data.Rat %>% filter(time >= 0.005) %>% 
  select(c("time", starts_with("AUR"))) %>% 
  rename("Time" = time) %>% 
  set_names(~(.) %>%
              str_replace_all("_gluc", "-gluc"))%>% 
  pivot_longer( !Time , 
                names_to = c("Tissue", "Substance"),
                names_sep = c("_"),
                values_to = "Amount")


Obs <- Observed.data %>%
  filter(AdminType == Select_AdminType,
         TestCompound == "ZEN", #other options are ZEN_glu and a_ZEL
         Author == Select_Author,
         SampleCompound == "ZEN"#,
         # IVPerfRate == IVDOSEmg
         # IVBolusDose == IVDOSEmg
         )%>%  
  select(c("Author", "Sex",
           "Tissue", "SamplingTime",
           "Concentration",
           "Cmin", "Cmax",
           "A_cum_umol", "A_cum_min_umol", "A_cum_max_umol"
           )) %>%
  rename(Time = SamplingTime) %>%
  # na.omit() %>% 
  mutate(Time = as.numeric(Time), 
         Concentration = as.numeric(Concentration),
         Concentration = ifelse(Tissue == "Plasma",
                               Concentration * ZEN_BPR, Concentration), # (ng/ml) blood concentration: changing from plasma concentrations to blood concentrations by dividing with BPR
         Cmin = ifelse(Tissue == "Plasma",
                                Cmin * ZEN_BPR, Cmin),
         Cmax = ifelse(Tissue == "Plasma",
                       Cmax * ZEN_BPR, Cmax),
         Tissue = ifelse(Tissue == "Plasma", "Blood", Tissue))


if(Select_Author == "Videmann"){
  Cobs <- Obs %>%
    filter(Tissue %in% c("Blood", "Intestine", "Liver"))
  Cpre <- Predicted.data %>% filter(Substance == "ZEN", 
                                    Tissue %in% c("Blood", "Intestine", "Liver"))# %>% 
    # filter(Time == 6)
}else{
  Cobs <- Obs %>% filter(Tissue == "Blood")
  Cpre <- Predicted.data %>% filter(Substance == "ZEN", 
                                    Tissue == "Blood")%>% 
    select(Substance, Time, Concentration) 
}

# AUC-halflife ####
Predicted <- Results.Rat$Raw.Data.Rat %>% select(time) %>% pull()
B <- Cobs %>% filter(Tissue == "Blood") #%>% filter(Sex == "Male")
b    <- approxfun(B$Time, B$Concentration, rule = 2)
# plot(b(Predicted))
t <- b(Predicted)
df <- tibble(Time = Predicted, Concentration = t)
AUC_ZEN <- trapz(Predicted, t)  # umol*h/L
print(AUC_ZEN)

time <- Predicted # h
conc_ZEN <- t # uM
Cmax_ZEN <- max(t)
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
  filename = here(OUTPUT,"PredvsObsBlood_ZEN.png"),
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
    filename = here(OUTPUT,"PredvsObsAll_ZEN.png"),
    dpi = 600,
    width = 11,
    height = 7,
    units = "cm"
  )

}else if(Select_Author == "Shin" && Select_AdminType == "IV.PERF"){
  
  Shin_obs <- Cobs %>% select(c(Tissue, Concentration, Cmin, Cmax)) %>% 
    mutate(Datatype = "Observed")
  Shin_pred <- Cpre %>% select(c(Tissue, Concentration)) %>% mutate(
    Cmin = Concentration, 
    Cmax = Concentration,
    Datatype = "Predicted")
  
  Shin_data <- rbind(Shin_obs, Shin_pred) 
  
  PredvsObsAll_Shin_A <- Shin_data %>%
    mutate(Tissue = reorder(Tissue, +Concentration)) %>% 
    ggplot(aes(x=Tissue, y=log10(Concentration))) +
    geom_segment(data = filter(Shin_data, Datatype == "Observed"), 
                 aes(x=Tissue, xend=Tissue, y=log10(Cmin), yend=log10(Cmax)),
                 color= "#F7D574",alpha = 0.5,  size=3 )+
    geom_point(aes(x=Tissue, y=log10(Concentration), color = Datatype), size=3 ) +
    scale_colour_manual(name = "Legend",
                        values = c( Observed = "#F7D574",
                                    Predicted = "#8631A7"),
                        labels = c( Observed = "Observed",
                                    Predicted = "Predicted"))  +
    labs(
      x = "",
      y = expression(paste(("Concentration (log"[10] ~ "(ng/ml))")))) +
    coord_cartesian(ylim = c(2, 4)) +
    guides(color = guide_legend(label.position = "right"),
           shape = guide_legend(label.position = "right"))+
    ggtitle("A.")+
    CP_theme +
    theme(legend.position = "none")
  PredvsObsAll_Shin_A
  
  PredvsObsAll_Shin_B <- Shin_data %>%
    mutate(Tissue = reorder(Tissue, +Concentration)) %>% 
    ggplot(aes(x=Tissue, y=log10(Concentration))) +
    geom_segment(data = filter(Shin_data, Datatype == "Observed"), 
                 aes(x=Tissue, xend=Tissue, y=log10(Cmin), yend=log10(Cmax)),
                 color= "#F7D574",alpha = 0.5,  size=3 )+
    geom_point(aes(x=Tissue, y=log10(Concentration), color = Datatype), size=3 ) +
    scale_colour_manual(name = "Legend",
                        values = c( Observed = "#F7D574",
                                    Predicted = "#8631A7"),
                        labels = c( Observed = "Observed",
                                    Predicted = "Predicted"))  +
    labs(
      x = "",
      y = "" )+ #expression(paste(("Concentration (log"[10] ~ "(ng/ml))")))) +
    coord_cartesian(ylim = c(2, 4)) +
    guides(color = guide_legend(label.position = "right"),
           shape = guide_legend(label.position = "right"))+
    ggtitle("B.")+
    CP_theme +
    theme(legend.position = "right",
          legend.text = element_text(size = 40),
          legend.title = element_blank(),
          legend.key.size = unit(0.5, "lines"),
          legend.spacing.y = unit(0, "cm"),
          legend.margin = margin(0, 0, 0, 0))  
  PredvsObsAll_Shin_B
  
  PredvsObsAll_ZEN_Shin <- PredvsObsAll_Shin_A + PredvsObsAll_Shin_B +
    plot_layout(widths = c(0.4, 0.4)) 
  
  ggsave(
    filename = here(OUTPUT, "PredvsObsAll_ZEN_Shin.png"),
    plot = PredvsObsAll_ZEN_Shin,
    dpi = 600,  
    width = 18,  
    height = 9, 
    units = "cm",
    limitsize = TRUE
  )
  
 
  
}else{
  print("only blood this time")}


### Urine ----

if(Select_Author == "Videmann"){
  
  
  Aobs_urine <- Obs %>%
    filter(Tissue == "Urine")
  Apre_urine <- Predicted.urine %>% filter(Substance == "ZEN") 
  
  ggplot(
    data = filter(Aobs_urine),
    aes(x = Time, y = A_cum_umol)) +
    geom_ribbon(aes(ymin = A_cum_min_umol, ymax = A_cum_max_umol, fill = Sex),
                alpha = 0.3, color = NA) +
    scale_fill_manual(values = c("#008475", "#8631A7"))+
    geom_point(aes(color = Sex),size = 1) +
    scale_color_manual(values = c("#008475", "#8631A7")) +
    geom_path(data = filter(Apre_urine, Tissue == "AUR"),
              aes(Time, Amount)) +
    labs(title = "Cumulative Amount of Zearalenone in urine",
         x = "Time (h)",
         y = "Amount (umol)") +
    guides(color = guide_legend(override.aes = list(fill = NA))) +
    CP_theme
  
  ggsave(
    filename = here(OUTPUT,"PredvsObsUrine_ZEN.png"),
    dpi = 600,
    width = 11,
    height = 7,
    units = "cm")
  
}else{
  
}




## ZENGLUC ----


### Blood ----
if(Select_Author == "Videmann"){

  Obs.glu <- Observed.data %>%
    filter(AdminType == Select_AdminType,
           Sex != "Pregnant.F", #pregnancy not included yet
           TestCompound == "ZEN", #other options are ZEN_glu and a_ZEL
           Author == Select_Author,
           SampleCompound == "ZEN_glu")   %>% 
    select(c("Author", "Sex",
             "Tissue", "SamplingTime",
             "Concentration",
             "Cmin", "Cmax",
             "A_cum_umol", "A_cum_min_umol", "A_cum_max_umol")) %>% #
    rename(Time = SamplingTime) %>%
    # na.omit() %>% 
    mutate(Time = as.numeric(Time), 
           Concentration = as.numeric(Concentration),
           Concentration = ifelse(Tissue == "Plasma",
                                  Concentration * ZEN_gluc_BPR, Concentration),
           Cmin = ifelse(Tissue == "Plasma",
                                  Cmin * ZEN_gluc_BPR, Cmin),
           Cmax = ifelse(Tissue == "Plasma",
                         Cmax * ZEN_gluc_BPR, Cmax),
           # (ng/ml) blood concentration: changing from plasma concentrations to blood concentrations by multiplying with BP
           Tissue = ifelse(Tissue == "Plasma", "Blood", Tissue))
  
  Cobs.glu <- Obs.glu %>%
      filter(Tissue %in% c("Blood","Intestine", "Liver"))
  Cpre.glu <- Predicted.data %>% filter(Substance == "ZEN-gluc", 
                                        Tissue %in% c("Blood", "Intestine", "Liver")) 
  
  PredvsObsBlood.GLU <- 
    ggplot(
      data = filter(Cobs.glu, Tissue == "Blood"),
      aes(x = Time, y = Concentration)) +
    geom_ribbon(aes(ymin = Cmin, ymax = Cmax, fill = Sex),
                alpha = 0.3, color = NA) +
    scale_fill_manual(values = c("#008475", "#8631A7"))+
    geom_point(aes(color = Sex),size = 1) +
    scale_color_manual(values = c("#008475", "#8631A7")) +
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
    filename = here(OUTPUT,"PredvsObsBlood_ZEN-GLUC.png"),
    dpi = 600,
    width = 11,
    height = 7,
    units = "cm"
  )
  
  ### All tissues ----
  PredvsObsAll.GLU <- ggplot(data = Cobs.glu,
                         aes(x = Time, y = Concentration)) +
    facet_wrap(~Tissue, scales = "free")+
    geom_ribbon(aes(ymin = Cmin, ymax = Cmax, fill = Sex),
                alpha = 0.3, color = NA) +
    scale_fill_manual(values = c("#008475", "#8631A7"))+
    geom_point(aes(color = Sex),size = 1) +
    scale_color_manual(values = c("#008475", "#8631A7"))+
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
    filename = here(OUTPUT,"PredvsObsAll_ZEN-GLUC.png"),
    dpi = 600,
    width = 11,
    height = 7,
    units = "cm"
  )

}else{ print("Only ZEN nothing else ;)")}

### Urine ----

if(Select_Author == "Videmann"){
  
  
  Aobs_urine.GLU  <- Obs.glu %>%
    filter(Tissue == "Urine")
  Apre_urine.GLU <- Predicted.urine %>% filter(Substance == "ZEN-gluc") 
  
  ggplot(
    data = filter(Aobs_urine.GLU),
    aes(x = Time, y = A_cum_umol)) +
    geom_ribbon(aes(ymin = A_cum_min_umol, ymax = A_cum_max_umol, fill = Sex),
                alpha = 0.3, color = NA) +
    scale_fill_manual(values = c("#008475", "#8631A7"))+
    geom_point(aes(color = Sex),size = 1) +
    scale_color_manual(values = c("#008475", "#8631A7")) +
    geom_path(data = filter(Apre_urine.GLU, Tissue == "AUR"),
              aes(Time, Amount)) +
    labs(title = "Cumulative Amount of Zearalenone glucuronide in urine",
         x = "Time (h)",
         y = "Amount (umol)") +
    guides(color = guide_legend(override.aes = list(fill = NA))) +
    CP_theme
  
  ggsave(
    filename = here(OUTPUT,"PredvsObsUrine_ZEN-GLUC.png"),
    dpi = 600,
    width = 11,
    height = 7,
    units = "cm")
  
  
}else{
  print("No Urine")
}



## a-ZEL ----

if(Select_Author == "Videmann"){
  
  Obs.aZEL <- Observed.data %>%
    filter(AdminType == Select_AdminType,
           Sex != "Pregnant.F", #pregnancy not included yet
           TestCompound == "ZEN", #other options are ZEN_glu and a_ZEL
           Author == Select_Author,
           SampleCompound == "a_ZEL")  %>% 
    select(c("Author", "Sex",
             "Tissue", "SamplingTime",
             "Concentration",
             "Cmin", "Cmax",
             "A_cum_umol", "A_cum_min_umol", "A_cum_max_umol")) %>% #, "Cmin", "Cmax"
    rename(Time = SamplingTime) %>%
    # na.omit() %>% 
    mutate(Time = as.numeric(Time), 
           Concentration = as.numeric(Concentration),
           Concentration = ifelse(Tissue == "Plasma",
                                  Concentration * aZEL_BPR, Concentration), # (ng/ml) blood concentration: changing from plasma concentrations to blood concentrations by multiplying with BP
           Tissue = ifelse(Tissue == "Plasma", "Blood", Tissue))
  
  Cobs.aZEL <- Obs.aZEL %>%
    filter(Tissue %in% c("Intestine", "Liver"))
  Cpre.aZEL <- Predicted.data %>% filter(Substance == "aZEL", 
                                        Tissue %in% c("Intestine", "Liver")) 
  
  ### All tissues ----
  PredvsObsAll.aZEL <- ggplot(data = Cobs.aZEL,
                             aes(x = Time, y = Concentration)) +
    facet_wrap(~Tissue, scales = "free")+
    geom_ribbon(aes(ymin = Cmin, ymax = Cmax, fill = Sex),
                alpha = 0.3, color = NA) +
    scale_fill_manual(values = c("#008475", "#8631A7"))+
    geom_point(aes(color = Sex),size = 1) +
    scale_color_manual(values = c("#008475", "#8631A7"))+
    geom_path(data = Cpre.aZEL,
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
    filename = here(OUTPUT,"PredvsObsAll_aZEL.png"),
    dpi = 600,
    width = 11,
    height = 7,
    units = "cm"
  )
  
}else{ print("Only ZEN nothing else ;)")}

### Urine ----

if(Select_Author == "Videmann"){
  
  Aobs_urine.aZEL  <- Obs.aZEL %>%
    filter(Tissue == "Urine")
  Apre_urine.aZEL <- Predicted.urine %>% filter(Substance == "aZEL") 
  
  ggplot(
    data = filter(Aobs_urine.aZEL),
    aes(x = Time, y = A_cum_umol)) +
    geom_ribbon(aes(ymin = A_cum_min_umol, ymax = A_cum_max_umol, fill = Sex),
                alpha = 0.3, color = NA) +
    scale_fill_manual(values = c("#008475", "#8631A7"))+
    geom_point(aes(color = Sex),size = 1) +
    scale_color_manual(values = c("#008475", "#8631A7")) +
    geom_path(data = filter(Apre_urine.aZEL, Tissue == "AUR"),
              aes(Time, Amount)) +
    labs(title = "Cumulative Amount of a-Zearalenol in urine",
         x = "Time (h)",
         y = "Amount (umol)") +
    guides(color = guide_legend(override.aes = list(fill = NA))) +
    CP_theme
  
  ggsave(
    filename = here(OUTPUT,"PredvsObsUrine_aZEL.png"),
    dpi = 600,
    width = 11,
    height = 7,
    units = "cm")
  
  
}else{
  print("No Urine")
  
}



## a-ZEL glucuronide ----

if(Select_Author == "Videmann"){
  
  Obs.aZEL_gluc <- Observed.data %>%
    filter(AdminType == Select_AdminType,
           Sex != "Pregnant.F", #pregnancy not included yet
           TestCompound == "ZEN", #other options are ZEN_glu and a_ZEL
           Author == Select_Author,
           SampleCompound == "a_ZEL_glu")   %>% 
    select(c("Author", "Sex",
             "Tissue", "SamplingTime",
             "Concentration",
             "Cmin", "Cmax",
             "A_cum_umol", "A_cum_min_umol", "A_cum_max_umol")) %>% #, "Cmin", "Cmax"
    rename(Time = SamplingTime) %>%
    # na.omit() %>% 
    mutate(Time = as.numeric(Time), 
           Concentration = as.numeric(Concentration),
           Concentration = ifelse(Tissue == "Plasma",
                                  Concentration * aZEL_gluc_BPR, Concentration), # (ng/ml) blood concentration: changing from plasma concentrations to blood concentrations by multiplying with BP
           Tissue = ifelse(Tissue == "Plasma", "Blood", Tissue))
  
  Cobs.aZEL_gluc <- Obs.aZEL_gluc %>%
    filter(Tissue %in% c("Intestine", "Liver"))
  Cpre.aZEL_gluc <- Predicted.data %>% filter(Substance == "aZEL-gluc", 
                                         Tissue %in% c("Intestine", "Liver")) 
  
  ### All tissues ----
  PredvsObsAll.aZEL_gluc <- ggplot(data = Cobs.aZEL_gluc,
                              aes(x = Time, y = Concentration)) +
    facet_wrap(~Tissue, scales = "free")+
    geom_ribbon(aes(ymin = Cmin, ymax = Cmax, fill = Sex),
                alpha = 0.3, color = NA) +
    scale_fill_manual(values = c("#008475", "#8631A7"))+
    geom_point(aes(color = Sex),size = 1) +
    scale_color_manual(values = c("#008475", "#8631A7"))+
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
    filename = here(OUTPUT,"PredvsObsAll_aZEL-GLUC.png"),
    dpi = 600,
    width = 11,
    height = 7,
    units = "cm"
  )
  
}else{ print("Only ZEN nothing else ;)")}

### Urine ----

if(Select_Author == "Videmann"){
  
  Aobs_urine.aZEL_gluc  <- Obs.aZEL_gluc %>%
    filter(Tissue == "Urine")
  Apre_urine.aZEL_gluc <- Predicted.urine %>% filter(Substance == "aZEL-gluc") 
  
  ggplot(
    data = filter(Aobs_urine.aZEL_gluc),
    aes(x = Time, y = A_cum_umol)) +
    geom_ribbon(aes(ymin = A_cum_min_umol, ymax = A_cum_max_umol, fill = Sex),
                alpha = 0.3, color = NA) +
    scale_fill_manual(values = c("#008475", "#8631A7"))+
    geom_point(aes(color = Sex),size = 1) +
    scale_color_manual(values = c("#008475", "#8631A7")) +
    geom_path(data = filter(Apre_urine.aZEL_gluc, Tissue == "AUR"),
              aes(Time, Amount)) +
    labs(title = "Cumulative Amount of a-Zearalenol glucuronide in urine",
         x = "Time (h)",
         y = "Amount (umol)") +
    guides(color = guide_legend(override.aes = list(fill = NA))) +
    CP_theme
  
  ggsave(
    filename = here(OUTPUT,"PredvsObsUrine_aZEL-GLUC.png"),
    dpi = 600,
    width = 11,
    height = 7,
    units = "cm")
  
}else{
  print("No Urine")
  
}




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
  ggsave(filename = here(OUTPUT, "CRegression_plot_ZEN.png"),
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
  ggsave(filename = here(OUTPUT, "CRegression_plot_ZEN.png"),
         plot = CRegression_plot,
         dpi = 600,
         width = 9,height = 7, units = "cm")

  }
