# --------------------------------------------------------------------------- #
# Script for Plotting 
# By: Chrysanthi Pachoulide
# Date: 18-08-2025
# --------------------------------------------------------------------------- #

rm(list=ls()) 

# Packages
library(colorspace)
library(here)
library(tidyverse)
library(PKNCA)
library(pracma)
library(showtext)
library(plotly)
font_add(family = "Garamond", regular = "GARA.TTF")
showtext_auto()

CP_theme <- theme_minimal() +
  theme(
    plot.background = element_rect(fill = 'white', color = NA),
    panel.background = element_rect(fill = 'white', color = NA),
    axis.text = element_text(size = 56),
    axis.title = element_text(size = 62),
    plot.title = element_text(size = 65, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 56, hjust = 0.5),
    legend.position = "right",
    legend.text = element_text(size = 45),
    legend.title = element_text(size = 50),
    legend.margin = margin(0, 0, 0, 0)
  )

# Set output storage directory
INPUT <- here("Output", "2025-08-18")
OUTPUT <- here("Output", format(Sys.Date(), "%Y-%m-%d"))
dir.create(OUTPUT, recursive = TRUE)
MWZEN <- 318.37


load(here(INPUT, "PO.Videmann.RawData_PBK.RData"))
Pred_PO.Videmann <- PBK_OUTPUT %>% 
  as.data.frame() %>% 
  select(c(time, CB_ZEN)) %>% 
  mutate(Dose = 1, 
         Author = "Videmann", 
         AdminType = "PO") 
rm(PBK_OUTPUT)

load(here(INPUT, "PO.Shin.RawData_PBK.RData"))
Pred_PO.Shin <- PBK_OUTPUT %>% 
  as.data.frame() %>% 
  select(c(time, CB_ZEN)) %>% 
  mutate(Dose = 8, 
         Author = "Shin",
         AdminType = "PO") 
rm(PBK_OUTPUT)

load(here(INPUT, "PO.Q.Lu.RawData_PBK.RData"))
Pred_PO.Q.Lu <- PBK_OUTPUT %>% 
  as.data.frame() %>% 
  select(c(time, CB_ZEN)) %>% 
  mutate(Dose = 6.65, 
         Author = "Q.Lu",
         AdminType = "PO") 
rm(PBK_OUTPUT)

load(here(INPUT, "IV.Q.Lu.RawData_PBK.RData"))
Pred_IV.Q.Lu <- PBK_OUTPUT %>% 
  as.data.frame() %>% 
  select(c(time, CB_ZEN)) %>% 
  mutate(Dose = 6.65, 
         Author = "Q.Lu",
         AdminType = "IV")
rm(PBK_OUTPUT)

load(here(INPUT, "IV.Shin.RawData_PBK.RData"))
Pred_IV.Shin <- PBK_OUTPUT %>% 
  as.data.frame() %>% 
  select(c(time, CB_ZEN)) %>% 
  mutate(Dose = 8, 
         Author = "Shin",
         AdminType = "IV") 
rm(PBK_OUTPUT)


Predicted <- rbind(Pred_PO.Q.Lu, Pred_PO.Shin, Pred_PO.Videmann,
                   Pred_IV.Q.Lu, Pred_IV.Shin) %>% 
  mutate(PredObs = "Predicted", 
         Dose = as.factor(Dose),
         CB_ZEN = CB_ZEN*MWZEN) %>% 
  rename(Concentration = CB_ZEN)

DataInvivo <- read.csv(here("Input", "InVivoData.csv"), header=1)

Observed <- DataInvivo %>%
  filter(AdminType %in% c("PO", "IV"),          ,
         Gender != "Pregnant.F", #pregnancy not included yet
         TestCompound == "ZEN", #other options are ZEN_glu and a_ZEL
         Author %in% c("Videmann", "Shin", "Q.Lu"),
         SampleCompound == "ZEN", 
         Tissue == "Plasma") %>% 
  select(c("AdminType", "Author", "Gender", "SamplingTime", "Concentration", "Cmin", "Cmax")) %>%
  rename(Time = SamplingTime,
         Sex = Gender) %>%
  na.omit() %>% 
  mutate(time = as.numeric(Time), 
         Concentration = as.numeric(Concentration),
         Concentration = Concentration * 0.67, # (ng/ml) blood concentration: changing from plasma concentrations to blood concentrations by multiplying with BP
         Cmin = Cmin * 0.67,
         Cmax = Cmax * 0.67) %>% 
  mutate(Dose = ifelse(Author == "Videmann", 1, 
                       ifelse(Author == "Shin", 8, 6.65)),
         PredObs = "Observed", 
         Dose = as.factor(Dose))


original <- "#008479"
lighter <- lighten(original, amount = 0.5)  # 20% lighter
darker  <- darken(original, amount = 0.2)   # 20% darker
print(lighter)  # e.g., "#00B4A2"
print(darker)   # e.g., "#00564F"
Dose_colors <- c("1" = "#70C7BB", "6.65" = "#008479", "8" = "#056960")
Fill_colours <- c("1" = "#70C7BB", "6.65" = "#008479", "8" = "#056960")


PredvsObsBlood.PO <- ggplot(data = filter(Predicted, AdminType == "PO"),
                         aes(x = time, y = Concentration, color = Dose)) +
  geom_path()+
  geom_point(data = filter(Observed, AdminType == "PO"),
            aes(x = time, y = Concentration, color = Dose, shape = Sex),
            size = 1) +
  geom_ribbon(data = filter(Observed, AdminType == "PO"),
              aes(ymin = Cmin, ymax = Cmax, fill = Dose),
              alpha = 0.2, color = NA) +
  scale_color_manual(values = Dose_colors) +
  scale_fill_manual(values = Dose_colors) +
  labs(title = "Observed vs Predicted Zearalenone blood concentrations",
       x = "Time (h)",
       y = "Concentration (ng/ml)") +
  guides(color = guide_legend(override.aes = list(fill = NA))) +
  CP_theme
PredvsObsBlood.PO
ggsave(
  filename = here(OUTPUT,"PredvsObsBlood.PO.png"),
  dpi = 1000,
  width = 11,
  height = 7,
  units = "cm"
)

PredvsObsBlood.IV <- ggplot(data = filter(Predicted, AdminType == "IV"),
                            aes(x = time, y = Concentration, color = Dose)) +
  geom_path()+
  geom_point(data = filter(Observed, AdminType == "IV"),
             aes(x = time, y = Concentration, color = Dose, shape = Sex),
             size = 1) +
  geom_ribbon(data = filter(Observed, AdminType == "IV"),
              aes(ymin = Cmin, ymax = Cmax, fill = Dose),
              alpha = 0.2, color = NA) +
  scale_color_manual(values = Dose_colors) +
  scale_fill_manual(values = Dose_colors) +
  labs(title = "Observed vs Predicted Zearalenone blood concentrations",
       x = "Time (h)",
       y = "Concentration (ng/ml)") +
  guides(color = guide_legend(override.aes = list(fill = NA))) +
  CP_theme
PredvsObsBlood.IV
ggsave(
  filename = here(OUTPUT,"PredvsObsBlood.IV.png"),
  dpi = 1000,
  width = 11,
  height = 7,
  units = "cm"
)



# Plasma C normalised to the dose


PredictedNormalised <- rbind(Pred_PO.Q.Lu, Pred_PO.Shin, Pred_PO.Videmann,
                   Pred_IV.Q.Lu, Pred_IV.Shin) %>% 
  mutate(PredObs = "Predicted", 
         CB_ZEN = CB_ZEN*MWZEN/Dose) %>% 
  rename(Concentration = CB_ZEN) %>% 
  mutate(Dose = as.factor(Dose))


ObservedNormalised <- DataInvivo %>%
  filter(AdminType %in% c("PO", "IV"),          ,
         Gender != "Pregnant.F", #pregnancy not included yet
         TestCompound == "ZEN", #other options are ZEN_glu and a_ZEL
         Author %in% c("Videmann", "Shin", "Q.Lu"),
         SampleCompound == "ZEN", 
         Tissue == "Plasma") %>% 
  select(c("AdminType", "Author", "Gender", "SamplingTime", "Concentration", "Cmin", "Cmax")) %>%
  rename(Time = SamplingTime,
         Sex = Gender) %>%
  na.omit() %>% 
  mutate(time = as.numeric(Time), 
         Concentration = as.numeric(Concentration),
         Concentration = Concentration * 0.67, # (ng/ml) blood concentration: changing from plasma concentrations to blood concentrations by multiplying with BP
         Cmin = Cmin * 0.67,
         Cmax = Cmax * 0.67) %>% 
  mutate(Dose = ifelse(Author == "Videmann", 1, 
                       ifelse(Author == "Shin", 8, 6.65)),
         PredObs = "Observed",
         Concentration = Concentration/Dose,
         Cmin = Cmin/Dose,
         Cmax = Cmax/Dose,
         Dose = as.factor(Dose))


original <- "#008479"
lighter <- lighten(original, amount = 0.5)  # 20% lighter
darker  <- darken(original, amount = 0.2)   # 20% darker
print(lighter)  # e.g., "#00B4A2"
print(darker)   # e.g., "#00564F"
Dose_colors <- c("1" = "#70C7BB", "6.65" = "#008479", "8" = "#056960")
Fill_colours <- c("1" = "#70C7BB", "6.65" = "#008479", "8" = "#056960")


PredvsObsBloodNormalised <- ggplot(data = filter(PredictedNormalised, AdminType == "PO"),
                         aes(x = time, y = Concentration, color = Dose)) +
  geom_path()+
  geom_point(data = filter(ObservedNormalised, AdminType == "PO"),
             aes(x = time, y = Concentration, color = Dose, shape = Sex),
             size = 1) +
  geom_ribbon(data = filter(ObservedNormalised, AdminType == "PO"),
              aes(ymin = Cmin, ymax = Cmax, fill = Dose),
              alpha = 0.2, color = NA) +
  scale_color_manual(values = Dose_colors) +
  scale_fill_manual(values = Dose_colors) +
  labs(title = "Observed vs Predicted Zearalenone blood concentrations",
       x = "Time",
       y = "Concentration normalised to the dose") +
  guides(color = guide_legend(override.aes = list(fill = NA))) +
  CP_theme
PredvsObsBloodNormalised
ggsave(
  filename = here(OUTPUT,"PredvsObsBloodNormalised.png"),
  dpi = 1000,
  width = 11,
  height = 7,
  units = "cm"
)


# Effect of CL_hydro

load(here(INPUT, "CL_hydro.0.001.RawData_PBK.RData"))
CL_hydro.0.001. <- PBK_OUTPUT %>% 
  as.data.frame() %>% 
  select(c(time, CB_ZEN)) %>% 
  mutate(CL_hydro = 0.001) 
rm(PBK_OUTPUT)

load(here(INPUT, "CL_hydro.0.01.RawData_PBK.RData"))
CL_hydro.0.01. <- PBK_OUTPUT %>% 
  as.data.frame() %>% 
  select(c(time, CB_ZEN)) %>% 
  mutate(CL_hydro = 0.01) 
rm(PBK_OUTPUT)

load(here(INPUT, "CL_hydro.0.1.RawData_PBK.RData"))
CL_hydro.0.1. <- PBK_OUTPUT %>% 
  as.data.frame() %>% 
  select(c(time, CB_ZEN)) %>% 
  mutate(CL_hydro = 0.1) 
rm(PBK_OUTPUT)

load(here(INPUT, "CL_hydro.1.RawData_PBK.RData"))
CL_hydro.1. <- PBK_OUTPUT %>% 
  as.data.frame() %>% 
  select(c(time, CB_ZEN)) %>% 
  mutate(CL_hydro = 1) 
rm(PBK_OUTPUT)

load(here(INPUT, "CL_hydro.10.RawData_PBK.RData"))
CL_hydro.10. <- PBK_OUTPUT %>% 
  as.data.frame() %>% 
  select(c(time, CB_ZEN)) %>% 
  mutate(CL_hydro = 10) 
rm(PBK_OUTPUT)

CL_hydro <- rbind(CL_hydro.0.001.,CL_hydro.0.01., CL_hydro.0.1., CL_hydro.1., CL_hydro.10.) %>% 
  mutate(CL_hydro = as.factor(CL_hydro))



original <- "#008479"
lightest <- lighten(original, amount = 0.9)  # 10% lighter
lighter  <- lighten(original, amount = 0.5)   # 10% darker
light  <- lighten(original, amount = 0.1)   # 10% darker
higher <- darken(original, amount = 0.5)
highest <- darken(original, amount = 0.9)

CL_colors <- c("0.001" = lightest, "0.01" = lighter, "0.1" = light, "1" = higher, "10" = highest)

CL_hydro.plot <-ggplot(data = CL_hydro,
                         aes(x = time, y = CB_ZEN, color = CL_hydro)) +
  geom_path() +
  scale_color_manual(values = CL_colors) +
  CP_theme +
  labs(title = "Zearalenone clood concentration at different glucuronide hydrolysis clearances",
       x = "Time (h)",
       y = "Concentration (uM)") 
CL_hydro.plot
ggsave(
  filename = here(OUTPUT,"CL_hydro.plot.png"),
  dpi = 1000,
  width = 11,
  height = 7,
  units = "cm"
)

  
# Effect of CL_bile

load(here(INPUT, "CL_bile.0.01.RawData_PBK.RData"))
CL_bile.0.01. <- PBK_OUTPUT %>% 
  as.data.frame() %>% 
  select(c(time, CB_ZEN)) %>% 
  mutate(CL_bile = 0.01) 
rm(PBK_OUTPUT)

load(here(INPUT, "CL_bile.0.1.RawData_PBK.RData"))
CL_bile.0.1. <- PBK_OUTPUT %>% 
  as.data.frame() %>% 
  select(c(time, CB_ZEN)) %>% 
  mutate(CL_bile = 0.1) 
rm(PBK_OUTPUT)

load(here(INPUT, "CL_bile.1.RawData_PBK.RData"))
CL_bile.1. <- PBK_OUTPUT %>% 
  as.data.frame() %>% 
  select(c(time, CB_ZEN)) %>% 
  mutate(CL_bile = 1) 
rm(PBK_OUTPUT)

load(here(INPUT, "CL_bile.10.RawData_PBK.RData"))
CL_bile.10. <- PBK_OUTPUT %>% 
  as.data.frame() %>% 
  select(c(time, CB_ZEN)) %>% 
  mutate(CL_bile = 10) 
rm(PBK_OUTPUT)

load(here(INPUT, "CL_bile.100.RawData_PBK.RData"))
CL_bile.100. <- PBK_OUTPUT %>% 
  as.data.frame() %>% 
  select(c(time, CB_ZEN)) %>% 
  mutate(CL_bile = 100) 
rm(PBK_OUTPUT)

CL_bile <- rbind(CL_bile.0.01.,CL_bile.0.1., CL_bile.1., CL_bile.10., CL_bile.100.) %>% 
  mutate(CL_bile = as.factor(CL_bile))

original <- "#008479"
lightest <- lighten(original, amount = 0.9)  # 10% lighter
lighter  <- lighten(original, amount = 0.5)   # 10% darker
light  <- lighten(original, amount = 0.1)   # 10% darker
higher <- darken(original, amount = 0.5)
highest <- darken(original, amount = 0.9)

CL_colors <- c("0.01" = lightest, "0.1" = lighter, "1" = light, "10" = higher, "100" = highest)

CL_bile.plot <-ggplot(data = CL_bile,
                       aes(x = time, y = CB_ZEN, color = CL_bile)) +
  geom_path() +
  scale_color_manual(values = CL_colors) +
  CP_theme +
  labs(title = "Zearalenone clood concentration at different glucuronide biliary excretion clearances",
       x = "Time (h)",
       y = "Concentration (uM)") 
CL_bile.plot
ggsave(
  filename = here(OUTPUT,"CL_bile.plot.png"),
  dpi = 1000,
  width = 11,
  height = 7,
  units = "cm"
)
