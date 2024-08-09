rm(list=ls()) # to clear out the global environment

library(tidyverse)
library(ggplot2)
library(here)
library(showtext)
font_add(family = "Garamond", regular = "GARA.TTF")
showtext_auto()

source(here("Script/Functions/PlottingFunction.R"))                               #call plotting function

WhatAmITesting <- "PlottingExperimental"
TodaysOuput <- here("Output/Data", Sys.Date(), WhatAmITesting)
dir.create(TodaysOuput, WhatAmITesting,recursive = TRUE)
PlotFolderIs <- here("Output/Plots", Sys.Date(), WhatAmITesting)
dir.create(PlotFolderIs, WhatAmITesting,recursive = TRUE)


## Plotting Experimental data ####
## ------------------------------------------------------

# Input in vivo data
DataInvivo <- read.csv(here("InputData/InVivoData.csv"), header=1) %>% 
  filter(AdminType == "PO", 
         Gender != "Pregnant.F", #pregnancy not included yet
         TestCompound == "ZEN", #other options are ZEN_glu and a_ZEL
         AdminType == "PO", 
         Tissue == "Plasma", 
         SampleCompound == "ZEN") %>%
  select(c("Author", "Gender", "SamplingTime", "Concentration")) %>%
  rename(Time = SamplingTime) %>% 
  mutate(Concentration = Concentration*0.67) #(ng/ml) blood concentration: changing from plasma concentrations to blood concentrations by multiplying with BP

DataInvivo %>% ggplot() +
  geom_point(aes(Time, Concentration, color = Author)) +
  geom_line(aes(Time, Concentration, color = Author),
            linewidth = 0.1) +
  facet_wrap(~Gender) +
  labs(title = "Experimental Zearalenone blood concentrations",
       x = "Time (h)",
       y = "Concentration (ng/ml)",
       color = "Author") +
  CP_theme()
ggsave(here(PlotFolderIs, "BloodExperimentalAll.png"), dpi = 300)
