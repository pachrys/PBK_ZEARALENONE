# --------------------------------------------------------------------------- #
# Script for the PBK Model of DIANA et al. 
# By: Chrysanthi Pachoulide
# Date: 23-06-2025
# --------------------------------------------------------------------------- #

rm(list=ls()) 

# Packages
library(here)
library(tidyverse)
library(deSolve)
library(PKNCA)
library(pracma)
library(showtext)
font_add(family = "Garamond", regular = "GARA.TTF")
showtext_auto()

CP_theme <- theme_minimal() +
  theme(
    axis.text = element_text(size = 56),
    axis.title = element_text(size = 62),
    plot.title = element_text(size = 65, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 56, hjust = 0.5),
    legend.position = "none",
    plot.margin = margin(0.3, 0.3, 0.3, 0.3, "cm") 
  )

# Set output storage directory
OUTPUT <- here("Output", format(Sys.Date(), "%Y-%m-%d"), format(Sys.time(), "%H-%M-%S"))
dir.create(OUTPUT, recursive = TRUE)

# Species
Species = "Human"

# Parameters
Parameters <- if(Species == "Rat"){
  read.csv(here("Input", "Parm.rat.csv"))
}else{
  read.csv(here("Input", "Parm.human.csv"))
}

# Dosing
# molecular weight
MWZEN <- 318.37 # molecular weight ZEN

# oral dose
ODOSEmg <- 1 #{mg/kg bw} # oral dose, variable
ODOSEumol <- ODOSEmg*1000/MWZEN*Parameters[1, "Value"] # (μmol) unit change to μmol

# IV dose
IVDOSEmg1 <- 0 #{mg/kg bw} # IV dose, variable
IVDOSEumol <-IVDOSEmg1*1000/MWZEN*Parameters[1, "Value"] # unit change to μmol
# time

Dosing <- data.frame(
  Parameter = c("ODOSEumol", "IVDOSEumol"),
  Value = c(ODOSEumol, IVDOSEumol)
)


Parameters <- rbind(Parameters, Dosing)
parm.c <- as.list(setNames(Parameters$Value, Parameters$Parameter))


# Model Code
source(here("Script", "PBK_code_Mendez.R"))

# Tell Initials 
A_init <- c(
  Ast = 0,
  Ain1 = 0,
  Ain2 = 0,
  Ain3 = 0,
  Ain4 = 0,
  Ain5 = 0,
  Ain6 = 0,
  Ain7 = 0,
  ASIZEN = 0,
  ASIZENGLU = 0,
  ALIZEN = 0,
  ALIaZEL = 0,
  ALIbZEL = 0,
  AfeZ = 0,
  AfeaZ = 0,
  ALZEN = 0,
  ALaZEL = 0,
  ALbZEL = 0,
  ALZENGLU = 0,
  ALaZELGLU = 0,
  AF = 0,
  AR = 0,
  AS = 0,
  AB = 0,
  AZur = 0,
  AFaZEL = 0,
  ARaZEL = 0,
  ASaZEL = 0,
  ABaZEL = 0,
  AaZur = 0,
  ATotal = 0
  )


PBK.out <- lsoda(y = A_init,
                 times = seq(0, 24, by=1/60),#per minute
                 func =  PBK.model,
                 parms = parm.c)

PBK_OUTPUT <- as.data.frame(PBK.out)

MB.df <- PBK_OUTPUT %>% select(time, ATotal, ACalculated, MB)  # change days to years if needed
MB.df$MB <- round(MB.df$MB, 10) # rounding significance points
MB.df$ERROR <- (MB.df$ATotal - MB.df$ACalculated) / MB.df$ACalculated * 100
MB.df$ERROR <- round(MB.df$ERROR, 10) # rounding significance points
MB_plot <- ggplot(data = MB.df)+
  geom_line(aes(x = time, y = ERROR, color = "ERROR")) +
  geom_line(aes(x = time, y = MB, color = "MB")) +
  scale_color_manual(values = c("ERROR" = "blue", "MB" = "black"), 
                     name = NULL) +
  # ylim(0,1) +
  theme_minimal() +
  # theme(axis.text = element_text(size = 10),
  #       axis.title = element_text(size = 12)
  # )+
  ylab("MB / ERROR")
MB_plot

C_organs.df <- PBK_OUTPUT %>% 
  transmute(time = time,
            Intestinal.lumen_ZEN = Cin1 + Cin2 + Cin3 + Cin4 + Cin5 + Cin6 + Cin7,
            Intestinal.tissue_ZEN = CSIZEN,
            Intestinal.microb_ZEN = CLIZEN,
            Liver_ZEN = CLZEN, 
            Fat_ZEN = CF, Rapidly_ZEN = CR, Slowly_ZEN = CS, Blood_ZEN = CB, 
            Intestinal.tissue_ZENGLU = CSIZENGLU,
            Intestinal.microb_aZEL = CLIaZEL,
            Intestinal.microb_bZEL = CLIbZEL,
            Liver_aZEL = CLaZEL, 
            Liver_bZEL = CLbZEL,
            Liver_ZENGLU = CLZENGLU,
            Liver_bZEL = CLbZEL,
            Liver_aZELGLU = CLaZELGLU,
            Fat_aZEL = CFaZEL, Rapidly_aZEL = CRaZEL, Slowly_aZEL = CSaZEL, 
            Blood_aZEL = CBaZEL) 

C_organs_long <- C_organs.df %>% 
  pivot_longer(names_to = c("Organ", "Compound"), 
               Intestinal.lumen_ZEN:Blood_aZEL, 
               values_to = "Concentration",
               names_sep = "_")

AllConcentrations <- C_organs_long %>%
  ggplot() +
  geom_path(aes(x = time, y = Concentration, colour = Compound), linewidth = 0.3) +
  facet_wrap(~Organ, scales = "free")+
  # scale_color_manual(name = "Model",
  #                    values = c("Predicted" = "blueviolet", "PredictedEFSA" = "darkblue"),
  #                    labels = c("Predicted" = "Predicted NG-PBK", "PredictedEFSA" = "Predicted EFSA 2020")
  # ) +
  # scale_linetype_manual(name = "Model",
  #                       values = c("Predicted" = "solid", "PredictedEFSA" = "dashed"),
  #                       labels = c("Predicted" = "Predicted NG-PBK", "PredictedEFSA" = "Predicted EFSA 2020")) +
  # guides(color = guide_legend(override.aes = list(linetype = c("solid", "dashed"))))+
  CP_theme +
  ylab("Organ Concentration (\u03BCmol/L)") +
  xlab("Time (hours)") +
  theme(legend.position = c(0.9, 0.2),
        legend.text = element_text(size = 40),
        legend.title = element_blank(),
        legend.key.size = unit(0.5, "lines"),
        legend.spacing.y = unit(0, "cm"),
        legend.margin = margin(0, 0, 0, 0),
        strip.text = element_text(size = 45),
        # axis.text.x = element_text(margin = margin(t = 0, b = -10))
  )
AllConcentrations
ggsave(
  filename = here(OUTPUT,"AllConcentrations.png"),
  dpi = 1000,
  width = 11,
  height = 7,
  units = "cm"
)

