rm(list=ls()) 


library(tidyverse)
library(showtext)
library(here)
font_add(family = "Garamond", regular = "GARA.TTF")
showtext_auto()

# Set output storage directory

OUTPUT <- here("Output", format(Sys.Date(), "%Y-%m-%d")) 
dir.create(OUTPUT, recursive = TRUE)

CP_theme <- theme_minimal() +
  theme(
    plot.background = element_rect(fill = 'white'), #element_blank(), #
    plot.margin = margin(0.3, 0.3, 0.3, 0.3, "cm"),
    plot.title.position = "plot", 
    panel.grid = element_line(linewidth = 0.1),  
    axis.title = element_text(size = 34, face = "bold"),
    axis.text = element_text(size = 34),
    axis.line = element_line(linewidth = 0.5),
    legend.text = element_text(size = 34),
    legend.title = element_text(size = 34), 
    legend.key.spacing = unit(0.1, "cm"),
    legend.key.size = unit(0.5, "lines"),
    legend.spacing.y = unit(0.2, "cm"), 
    legend.spacing.x = unit(0.2, "cm"),
    legend.margin = margin(0,0,0,0),
    legend.position = "right"
  )

# Initial model (Videmann) ----
load(here(OUTPUT, paste0("Videmann", "Results.Rat.RData")))
Initial_C <- Results.Rat$C.Tissues.Rat
Initial_B <- Initial_C %>% filter(Tissue == "Blood") %>% 
  mutate(Test = "Initial", 
         Method = "On",
         graph = "On-Off Method")
Initial <- Results.Rat$TK.data %>% mutate(Test = "Initial") %>% 
  mutate(Test = "Initial", 
         Method = "On",
         graph = "On-Off Method")

# Updated ----
# off
load(here(OUTPUT, paste0("Nothing", "Results.Rat.RData")))
None_C <- Results.Rat$C.Tissues.Rat

None_BT <- None_C %>% filter(Tissue == "Blood") %>% 
  mutate(Test = "None",
         Method = "Off",
         graph = "Updated")

None_T <- Results.Rat$TK.data %>% 
  mutate(Test = "None",
         Method = "Off",
         graph = "Updated")

# on
load(here(OUTPUT, paste0("All", "Results.Rat.RData")))
Original_C <- Results.Rat$C.Tissues.Rat
Original_B <- Original_C %>% filter(Tissue == "Blood") %>% 
  mutate(Test = "All", 
         Method = "On",
         graph = "Updated")
Original <- Results.Rat$TK.data %>% 
  mutate(Test = "All", 
         Method = "On",
         graph = "Updated")

# increased
load(here(OUTPUT, paste0("Increased", "Results.Rat.RData")))
Increased_C <- Results.Rat$C.Tissues.Rat
Increased_B <- Increased_C %>% filter(Tissue == "Blood") %>% 
  mutate(Test = "Increased", 
         Method = "Increased",
         graph = "Updated")
Increased <- Results.Rat$TK.data %>% 
  mutate(Test = "Increased", 
         Method = "Increased",
         graph = "Updated")

# Lymph ----

# off 

None_BL <- None_C %>% filter(Tissue == "Blood") %>% 
  mutate(Test = "None",
         Method = "Off",
         graph = "Lymph")
rm(Results.Rat)
load(here(OUTPUT, paste0("Nothing", "Results.Rat.RData")))
None_L <- Results.Rat$TK.data %>% 
  mutate(Test = "None",
         Method = "Off",
         graph = "Lymph")

# on
load(here(OUTPUT, paste0("OnlyLymph", "Results.Rat.RData")))
OnlyLymph_C <- Results.Rat$C.Tissues.Rat
OnlyLymph_B <- OnlyLymph_C %>% filter(Tissue == "Blood") %>% 
  mutate(Test = "OnlyLymph", 
         Method = "On",
         graph = "Lymph")
OnlyLymph <- Results.Rat$TK.data %>% 
  mutate(Test = "OnlyLymph", 
         Method = "On",
         graph = "Lymph")

# Increased 
load(here(OUTPUT, paste0("Lymph10", "Results.Rat.RData")))
OnlyLymphIncreased_C <- Results.Rat$C.Tissues.Rat
OnlyLymphIncreased_B <- OnlyLymphIncreased_C %>% 
  filter(Tissue == "Blood") %>% 
  mutate(Test = "OnlyLymphIncreased", 
         Method = "Increased",
         graph = "Lymph")
OnlyLymphIncreased <- Results.Rat$TK.data %>% 
  mutate(Test = "OnlyLymphIncreased", 
         Method = "Increased",
         graph = "Lymph")


# combined
load(here(OUTPUT, paste0("LymphandBile", "Results.Rat.RData")))
LymphandBile_C <- Results.Rat$C.Tissues.Rat
LymphandBile_B <- LymphandBile_C %>% filter(Tissue == "Blood") %>% 
  mutate(Test = "LymphandBile", 
         Method = "Combined_LB",
         graph = "Lymph")
LymphandBile <- Results.Rat$TK.data %>% 
  mutate(Test = "LymphandBile", 
         Method = "Combined_LB",
         graph = "Lymph")


# Bile ----

# off
None_BB <- None_C %>% filter(Tissue == "Blood") %>% 
  mutate(Test = "None",
         Method = "Off",
         graph = "Bile")

rm(Results.Rat)
load(here(OUTPUT, paste0("Nothing", "Results.Rat.RData")))
None_B <- Results.Rat$TK.data %>% 
  mutate(Test = "None",
         Method = "Off",
         graph = "Bile")

# on
load(here(OUTPUT, paste0("OnlyBile", "Results.Rat.RData")))
OnlyBile_C <- Results.Rat$C.Tissues.Rat
OnlyBile_B <- OnlyBile_C %>% filter(Tissue == "Blood") %>% 
  mutate(Test = "OnlyBile", 
         Method = "On",
         graph = "Bile")
OnlyBile <- Results.Rat$TK.data %>% 
  mutate(Test = "OnlyBile", 
         Method = "On",
         graph = "Bile")

# increased 
load(here(OUTPUT, paste0("Bile10", "Results.Rat.RData")))
OnlyBileIncreased_C <- Results.Rat$C.Tissues.Rat
OnlyBileIncreased_B <- OnlyBileIncreased_C %>% filter(Tissue == "Blood") %>% 
  mutate(Test = "OnlyBileIncreased", 
         Method = "Increased",
         graph = "Bile")
OnlyBileIncreased <- Results.Rat$TK.data %>% 
  mutate(Test = "OnlyBileIncreased", 
         Method = "Increased",
         graph = "Bile")

# combined
load(here(OUTPUT, paste0("HydroandBile", "Results.Rat.RData")))
HydroandBile_C <- Results.Rat$C.Tissues.Rat
BileandHydro_B <- HydroandBile_C %>% filter(Tissue == "Blood") %>% 
  mutate(Test = "HydroandBile", 
         Method = "Combined_HB",
         graph = "Bile")
BileandHydro <- Results.Rat$TK.data %>% 
  mutate(Test = "HydroandBile", 
         Method = "Combined_HB",
         graph = "Bile")


# combined
load(here(OUTPUT, paste0("LymphandBile", "Results.Rat.RData")))
LymphandBile_C <- Results.Rat$C.Tissues.Rat
BileandLymph_B <- LymphandBile_C %>% filter(Tissue == "Blood") %>% 
  mutate(Test = "LymphandBile", 
         Method = "Combined_LB",
         graph = "Bile")
BileandLymph <- Results.Rat$TK.data %>% 
  mutate(Test = "LymphandBile", 
         Method = "Combined_LB",
         graph = "Bile")

# Hydrolysis ----

# off
None_BH <- None_C %>% filter(Tissue == "Blood") %>% 
  mutate(Test = "None",
         Method = "Off",
         graph = "Hydro")

rm(Results.Rat)
load(here(OUTPUT, paste0("Nothing", "Results.Rat.RData")))
None_H <- Results.Rat$TK.data %>% 
  mutate(Test = "None",
         Method = "Off",
         graph = "Hydro")

# on
load(here(OUTPUT, paste0("OnlyHydro", "Results.Rat.RData")))
OnlyHydro_C <- Results.Rat$C.Tissues.Rat
OnlyHydro_B <- OnlyHydro_C %>% filter(Tissue == "Blood") %>% 
  mutate(Test = "OnlyHydro", 
         Method = "On",
         graph = "Hydro")
OnlyHydro <- Results.Rat$TK.data %>% 
  mutate(Test = "OnlyHydro", 
         Method = "On",
         graph = "Hydro")

# increased
load(here(OUTPUT, paste0("Hydro10", "Results.Rat.RData")))
OnlyHydroIncreased_C <- Results.Rat$C.Tissues.Rat
OnlyHydroIncreased_B <- OnlyHydro_C %>% filter(Tissue == "Blood") %>% 
  mutate(Test = "OnlyHydro", 
         Method = "On",
         graph = "Hydro")
OnlyHydroIncreased <- Results.Rat$TK.data %>% 
  mutate(Test = "OnlyHydroIncreased", 
         Method = "Increased",
         graph = "Hydro")

# combined
load(here(OUTPUT, paste0("HydroandBile", "Results.Rat.RData")))
HydroandBile_C <- Results.Rat$C.Tissues.Rat
HydroandBile_B <- HydroandBile_C %>% filter(Tissue == "Blood") %>% 
  mutate(Test = "HydroandBile", 
         Method = "Combined_HB",
         graph = "Hydro")
HydroandBile <- Results.Rat$TK.data %>% 
  mutate(Test = "HydroandBile", 
         Method = "Combined_HB",
         graph = "Hydro")

# Collect data and plot ----

colouring <- c(
  "Off" = "#FF8C00",
  "On" = "#008475",
  "Increased" = "#8631A7",
  "Combined_HB" = "#1B40FF",
  "Combined_LB" = "#1B90FF"
)



Collect.B <- rbind(None_BL, OnlyLymph_B, OnlyLymphIncreased_B, LymphandBile_B,
                   None_BB, OnlyBile_B, OnlyBileIncreased_B, BileandLymph_B, BileandHydro_B,
                   None_BH, OnlyHydro_B, OnlyHydroIncreased_B, HydroandBile_B,
                   None_BT, Original_B, Increased_B, 
                   Initial_B)

Collect.B$Substance <- factor(Collect.B$Substance, 
                               levels = c("ZEN", "ZEN-gluc", 
                                          "aZEL", "aZEL-gluc"))
Collect.B$graph <- factor(Collect.B$graph, 
                              levels = c("Lymph", "Bile", "Hydro", "Updated", "On-Off Method"))
Collect.B$Method <- factor(Collect.B$Method, 
                            levels = c("Off", "On", "Increased", "Combined_HB", "Combined_LB"))


Blood <- ggplot(Collect.B, aes(Time, Concentration,
                      color = Method, linetype = Method)) +
  geom_path() +
  scale_color_manual(name = "",
                     values = colouring) +
  scale_linetype_manual(name = "",
                        values = c(
                          "Off" = "dotted",
                          "On" = "solid",
                          "Increased" = "dashed",
                          "Combined_HB" = "twodash" ,
                          "Combined_LB"= "dotdash"
                        )) +
  facet_grid(Substance ~ graph, scales = "free") +
  CP_theme +
  theme(strip.text = element_text(size = 44),
        strip.text.y = element_text(angle = 0),
        legend.position = "bottom")
Blood
ggsave(
  filename = here("Blood.png"),
  dpi = 600,
  width = 16,
  height = 16,
  units = "cm"
)


Collect.TK <- rbind(None_L, OnlyLymph, LymphandBile, #OnlyLymphIncreased,
                   None_B, OnlyBile, BileandLymph, BileandHydro, #OnlyBileIncreased,
                   None_H, OnlyHydro, HydroandBile, #OnlyHydroIncreased,
                   None_T, Original, #Increased, 
                   Initial) 

# Collect.TK <- Collect.TK %>% mutate(
#   Cmax = ifelse(Mycotoxin == "ZEN", Cmax*318.15, Cmax),
#   Cmax = ifelse(Mycotoxin == "aZEL", Cmax*320.16, Cmax),
#   Cmax = ifelse(Mycotoxin == "ZEN_gluc", Cmax*494.18, Cmax),
#   Cmax = ifelse(Mycotoxin == "aZEL_gluc", Cmax*496.191, Cmax)
# ) %>% mutate(
#   AUC = ifelse(Mycotoxin == "ZEN", AUC*318.15, AUC),
#   AUC = ifelse(Mycotoxin == "aZEL", AUC*320.16, AUC),
#   AUC = ifelse(Mycotoxin == "ZEN_gluc", AUC*494.18, AUC),
#   AUC = ifelse(Mycotoxin == "aZEL_gluc", AUC*496.191, AUC)
# )

Collect.TK <- Collect.TK %>% 
  pivot_longer(Cmax:Half_life, names_to = "TK_parameter", 
               values_to = "Value") %>% na.omit() 
Collect.TK$Mycotoxin = str_replace_all(Collect.TK$Mycotoxin, "_gluc", "-gluc")
# Collect.TK$TK_parameter = str_replace_all(Collect.TK$TK_parameter, "Cmax", "nCMax")
# Collect.TK$TK_parameter = str_replace_all(Collect.TK$TK_parameter, "AUC", "nAUC")


Collect.TK$Mycotoxin <- factor(Collect.TK$Mycotoxin, 
                              levels = c("ZEN", "ZEN-gluc", 
                                         "aZEL", "aZEL-gluc"))
Collect.TK$graph <- factor(Collect.TK$graph, 
                           levels = c("Lymph", "Bile", "Hydro", "Updated", "On-Off Method"))
Collect.TK$Method <- factor(Collect.TK$Method, 
                            levels = c("Off", "On", "Increased", "Combined_HB", "Combined_LB"))

TKParams <- ggplot(Collect.TK,
                          aes(x =graph, y = Value,
                       fill = Method)) +
  geom_col(position = "dodge") +
  scale_fill_manual(name = "",
                    values = colouring) +
  facet_grid(TK_parameter ~ Mycotoxin, scales = "free") +
  CP_theme +
  labs(y = "Parameter value", 
       x = "") +
  theme(strip.text = element_text(size = 44),
        strip.text.y = element_text(angle = 0),
        legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1))
TKParams

ggsave(
  filename = here("TKParams_All.png"),
  dpi = 600,
  width = 16,
  height = 16,
  units = "cm"
)

CalcAUC <- Collect.B %>% filter(graph %in% c("Updated", "On-Off Method"),
                                Substance == "ZEN",
                                Time <= 24) 

Upd <- CalcAUC %>% filter(graph == "Updated")
Old <- CalcAUC %>% filter(graph == "On-Off Method")
AUC_Upd <- trapz(Upd$Time, Upd$Concentration)  # ng*h/ml
AUC_Old <- trapz(Old$Time, Old$Concentration)  # ng*h/ml
Cmax_Upd <- max(Upd$Concentration)
Cmax_Old <- max(Old$Concentration)

# 
# TKParams_ZEN <- ggplot(filter(Collect.TK, Mycotoxin == "ZEN"),
#                    aes(x = TK_parameter, y = Value,
#                        fill = Method)) +
#   geom_col(position = "dodge") +
#   scale_fill_manual(name = "",
#                     values = colouring) +
#   facet_grid(~ graph, scales = "free") +
#   CP_theme +
#   labs(y = "Parameter value", 
#        x = "") +
#   theme(strip.text = element_text(size = 40),
#         strip.text.y = element_text(angle = 0),
#         legend.position = "bottom",
#         axis.text.x = element_text(angle = 45, hjust = 1))
# TKParams_ZEN
# 
# ggsave(
#   filename = here("TKParams_ZEN.png"),
#   dpi = 600,
#   width = 11,
#   height = 9,
#   units = "cm"
# )

AUC_Upd <- trapz(Upd$Time, Upd$Concentration)  # ng*h/ml
AUC_Old <- trapz(Old$Time, Old$Concentration)  # ng*h/ml
Cmax_Upd <- max(Upd$Concentration)
Cmax_Old <- max(Old$Concentration)

Observed <- tibble(
  Mycotoxin = c("ZEN"),
  nCmax = c(9.93275/1, 2.131518/8, 54.6921/6.65, Cmax_Upd, Cmax_Old),
  Tmax = c(0.01, 0.08, 0.75, 0.21, 0.22),
  nAUC = c(168.571/1, 38.4794/8, 554.7289/6.65, AUC_Upd, AUC_Old),
  Half_life = c(37.71898, 58.13978, 278.3939, NA, 2.434054),
  Identity = c("Videmann", "Shin", "Q.Lu", "Updated", "On-Off Method"),
  Method = c("Observed", "Observed","Observed", "Predicted", "Predicted")
) %>% pivot_longer(nCmax:Half_life, names_to = "TK_parameter", values_to = "Value")

Observed$Identity <- factor(Observed$Identity, 
                            levels = c("Videmann", "Shin", "Q.Lu", "Updated", "On-Off Method"))

Observed_plot <- 
  ggplot(Observed, 
                    aes(x = Identity, y = Value, fill = Method)) +
  geom_col(position = "dodge") +
  # scale_fill_manual(name = "",
  #                   values = colouring) +
  facet_wrap(~ TK_parameter, scales = "free", nrow = 1) +
  CP_theme +
  labs(y = "Parameter value", 
       x = "")  +
  theme(strip.text = element_text(size = 44),
        legend.position = "",
        axis.text.x = element_text(angle = 45, hjust = 1))

Observed_plot

ggsave(
  filename = here("Observed_plot.png"),
  dpi = 600,
  width = 12,
  height = 6,
  units = "cm"
)

