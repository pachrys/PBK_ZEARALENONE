# ------------------------------------------------------------------------ ###
# Script for plotting results 
# ------------------------------------------------------------------------ ###


rm(list=ls())

library(here)
library(tidyverse)
library(showtext)
library(plotly)
library(readr)
font_add(family = "Garamond", regular = "GARA.TTF")
showtext_auto()

OUTPUT <- here("Output", format(Sys.Date(), "%Y-%m-%d")) 
dir.create(OUTPUT, recursive = TRUE)

CP_theme <- theme_minimal() +
  theme(
    plot.background = element_rect(fill = 'white'), #element_blank(), #
    plot.margin = margin(0.3, 0.3, 0.3, 0.3, "cm"), 
    plot.title = element_text(hjust = 0, size = 48, face = "bold"),
    plot.title.position = "plot", 
    # plot.caption = element_text(hjust = 0, size = 40, face= "italic"),
    # plot.caption.position = "plot",
    panel.background = element_blank(), #rect(fill = 'white', color = NA),
    panel.grid = element_line(size = 0.1),  
    axis.title = element_text(size = 44, face = "bold"),
    axis.text = element_text(size = 44),
    # axis.text.y = element_blank(),
    # axis.text.x = element_text(size = 64, face = "bold"),
    axis.line = element_line(size = 0.5),
    legend.text = element_text(size = 44),
    legend.title = element_text(size = 44), 
    legend.key.spacing = unit(0.1, "cm"),
    legend.key.size = unit(0.5, "lines"),
    legend.spacing.y = unit(0.2, "cm"), 
    legend.spacing.x = unit(0.2, "cm"),
    legend.margin = margin(0,0,0,0),
    legend.position = "right"
  )



load(here(OUTPUT, "NOEL.RawData_ZEN_PBK.RData"))
NOEL.Cmax <- MODEL_OUTPUT$Cmaxs
NOEL.exp <- MODEL_OUTPUT$Meta_data$exp %>% str_split_i('_', 1) %>% as.numeric() #umol/kg bw/day
NOEL.exp <- round(NOEL.exp*318.15, 1) #ug/kg bw/day

load(here(OUTPUT, "TDI.RawData_ZEN_PBK.RData"))
TDI.Cmax <- MODEL_OUTPUT$Cmaxs
TDI.exp <- MODEL_OUTPUT$Meta_data$exp %>% str_split_i('_', 1) %>% as.numeric() #umol/kg bw/day
TDI.exp <- round(TDI.exp*318.15, 1) #ug/kg bw/day

load(here(OUTPUT, "EC50.ZEN.0.06.RawData_ZEN_PBK.RData"))
EC50.ZEN.0.06.Cmax <- MODEL_OUTPUT$Cmaxs
EC50.ZEN.0.06.exp <- MODEL_OUTPUT$Meta_data$exp %>% str_split_i('_', 1) %>% as.numeric() #umol/kg bw/day
EC50.ZEN.0.06.exp <- round(EC50.ZEN.0.06.exp*318.15, 1) #ug/kg bw/day

load(here(OUTPUT, "EC50.ZEN_aZEL0.06.RawData_ZEN_PBK.RData"))
EC50.ZEN_aZEL.0.06.Cmax <- MODEL_OUTPUT$Cmaxs
EC50.ZEN_aZEL.0.06.exp <- MODEL_OUTPUT$Meta_data$exp %>% str_split_i('_', 1) %>% as.numeric() #umol/kg bw/day
EC50.ZEN_aZEL.0.06.exp <- round(EC50.ZEN_aZEL.0.06.exp*318.15, 1) #ug/kg bw/day

load(here(OUTPUT, "EC50.ZEN_aZEL.0.31.RawData_ZEN_PBK.RData"))
EC50.ZEN_aZEL.0.31.Cmax <- MODEL_OUTPUT$Cmaxs
EC50.ZEN_aZEL.0.31.exp <- MODEL_OUTPUT$Meta_data$exp %>% str_split_i('_', 1) %>% as.numeric() #umol/kg bw/day
EC50.ZEN_aZEL.0.31.exp <- round(EC50.ZEN_aZEL.0.31.exp*318.15, 1) #ug/kg bw/day

load(here(OUTPUT, "EC50.ZEN_aZEL.0.49.RawData_ZEN_PBK.RData"))
EC50.ZEN_aZEL.0.49.Cmax <- MODEL_OUTPUT$Cmaxs
EC50.ZEN_aZEL.0.49.exp <- MODEL_OUTPUT$Meta_data$exp %>% str_split_i('_', 1) %>% as.numeric() #umol/kg bw/day
EC50.ZEN_aZEL.0.49.exp <- round(EC50.ZEN_aZEL.0.49.exp*318.15, 1) #ug/kg bw/day

load(here(OUTPUT, "EC50.ZEN_aZEL.1.58.RawData_ZEN_PBK.RData"))
EC50.ZEN_aZEL.1.58.Cmax <- MODEL_OUTPUT$Cmaxs
EC50.ZEN_aZEL.1.58.exp <- MODEL_OUTPUT$Meta_data$exp %>% str_split_i('_', 1) %>% as.numeric() #umol/kg bw/day
EC50.ZEN_aZEL.1.58.exp <- round(EC50.ZEN_aZEL.1.58.exp*318.15, 1) #ug/kg bw/day

load(here(OUTPUT, "EC50.ZEN_aZEL.1.60.RawData_ZEN_PBK.RData"))
EC50.ZEN_aZEL.1.60.Cmax <- MODEL_OUTPUT$Cmaxs
EC50.ZEN_aZEL.1.60.exp <- MODEL_OUTPUT$Meta_data$exp %>% str_split_i('_', 1) %>% as.numeric() #umol/kg bw/day
EC50.ZEN_aZEL.1.60.exp <- round(EC50.ZEN_aZEL.1.60.exp*318.15, 1) #ug/kg bw/day

load(here(OUTPUT, "EC50.ZEN_aZEL.1.64.RawData_ZEN_PBK.RData"))
EC50.ZEN_aZEL.1.64.Cmax <- MODEL_OUTPUT$Cmaxs
EC50.ZEN_aZEL.1.64.exp <- MODEL_OUTPUT$Meta_data$exp %>% str_split_i('_', 1) %>% as.numeric() #umol/kg bw/day
EC50.ZEN_aZEL.1.64.exp <- round(EC50.ZEN_aZEL.1.64.exp*318.15, 1) #ug/kg bw/day

load(here(OUTPUT, "EC50.ZEN_aZEL.3.81.RawData_ZEN_PBK.RData"))
EC50.ZEN_aZEL.3.81.Cmax <- MODEL_OUTPUT$Cmaxs
EC50.ZEN_aZEL.3.81.exp <- MODEL_OUTPUT$Meta_data$exp %>% str_split_i('_', 1) %>% as.numeric() #umol/kg bw/day
EC50.ZEN_aZEL.3.81.exp <- round(EC50.ZEN_aZEL.3.81.exp*318.15, 1) #ug/kg bw/day

IVIVE.data <- tibble(
  HBGV = c("NOEL",                               "TDI",       "EC50",                "EC50",                    "EC50",                      "EC50",                      "EC50",                     "EC50",                     "EC50",                    "EC50"),
  Exposure = c(NOEL.exp, TDI.exp, EC50.ZEN.0.06.exp, EC50.ZEN_aZEL.0.06.exp, EC50.ZEN_aZEL.0.31.exp, EC50.ZEN_aZEL.0.49.exp, EC50.ZEN_aZEL.1.58.exp, EC50.ZEN_aZEL.1.60.exp, EC50.ZEN_aZEL.1.64.exp, EC50.ZEN_aZEL.3.81.exp),
  ZEN = c(NOEL.Cmax[1], TDI.Cmax[1], EC50.ZEN.0.06.Cmax[1], EC50.ZEN_aZEL.0.06.Cmax[1], EC50.ZEN_aZEL.0.31.Cmax[1], EC50.ZEN_aZEL.0.49.Cmax[1], EC50.ZEN_aZEL.1.58.Cmax[1], EC50.ZEN_aZEL.1.60.Cmax[1], EC50.ZEN_aZEL.1.64.Cmax[1], EC50.ZEN_aZEL.3.81.Cmax[1]),
  RPF.aZEL = c(NOEL.Cmax[5], TDI.Cmax[5], EC50.ZEN.0.06.Cmax[5], EC50.ZEN_aZEL.0.06.Cmax[5], EC50.ZEN_aZEL.0.31.Cmax[5], EC50.ZEN_aZEL.0.49.Cmax[5], EC50.ZEN_aZEL.1.58.Cmax[5], EC50.ZEN_aZEL.1.60.Cmax[5], EC50.ZEN_aZEL.1.64.Cmax[5], EC50.ZEN_aZEL.3.81.Cmax[5]),
  RPF.glucuronides = c(NOEL.Cmax[6], TDI.Cmax[6], EC50.ZEN.0.06.Cmax[6], EC50.ZEN_aZEL.0.06.Cmax[6], EC50.ZEN_aZEL.0.31.Cmax[6], EC50.ZEN_aZEL.0.49.Cmax[6], EC50.ZEN_aZEL.1.58.Cmax[6], EC50.ZEN_aZEL.1.60.Cmax[6], EC50.ZEN_aZEL.1.64.Cmax[6], EC50.ZEN_aZEL.3.81.Cmax[6])
)

IVIVE.data <- IVIVE.data %>% mutate(RPF.aZEL = round(RPF.aZEL,4),
                                    ZEN = round(ZEN,4),
                                    RPF.glucuronides = round(RPF.glucuronides,4))


IVIVE.data.long <- IVIVE.data %>% 
  pivot_longer(ZEN:RPF.glucuronides,
               names_to = "Method",
               values_to = "Cmax")

# IVIVE.data.long$Exposure <- as.factor(IVIVE.data.long$Exposure)
# IVIVE.data.long$Cmax <- as.factor(IVIVE.data.long$Cmax)

IVIVE.data.long <- IVIVE.data.long %>% group_by(HBGV)


# df$Group <- ifelse(grepl("^EC50_ZEN", df$What), "EC50_ZEN", df$What)

library(scales)
data.box <- IVIVE.data.long %>% filter(Method == "RPF.aZEL",
                                       HBGV == "EC50")
point.box <- IVIVE.data.long %>% filter(Method == "RPF.aZEL",
                                        HBGV !="EC50")

EC50plot <- ggplot() +
  geom_violin(data = data.box, aes(x = log(Cmax), y = Method), 
              fill = "grey90", color = "grey50") +
  geom_point(data = point.box, aes(x = log(Cmax), y = Method, 
                                   colour = HBGV,
                                   shape = HBGV), size = 2) +
  coord_flip() +
  scale_color_manual(values = c("NOEL" = "#1B35FF",
                                "TDI" = "#008475"))+
  labs(title = "B.",
       x = "Cmax (nM)",
       y = "") +
  CP_theme
ggsave(
  filename = here("EC50plot.png"),
  plot = EC50plot,
  dpi = 600, 
  width = 5, 
  height = 9, 
  units = "cm")  


RPF.method <- IVIVE.data.long %>% filter(Exposure <= 45)

RPFplot <- ggplot(transform(RPF.method, 
                            Method = factor(Method, levels = c("ZEN", "RPF.aZEL", "RPF.glucuronides"))),
       aes(x = as.factor(Exposure), y = Cmax, fill = HBGV)) +
  geom_col(position = "dodge") +
  geom_hline(yintercept = 0.06, linewidth = 0.5, linetype = "dashed")+
  facet_wrap(~Method, scales = "free") +
  scale_fill_manual(values = c("NOEL" = "#1B35FF",
                               "TDI" = "#008475",
                               "EC50" = "#8631A7"))+
  labs(title = "B.",
       x = "Exposure (μg/kg bw/day)",
       y = "Cmax (nM)") +
  CP_theme +
  theme(strip.text.x = element_text(size = 44))
RPFplot
ggsave(
  filename = here("RPFplot.png"),
  plot = RPFplot,
  dpi = 600, 
  width = 12, 
  height = 9, 
  units = "cm")  

