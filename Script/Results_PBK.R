# # --------------------------------------------------------------------------- #
# # Script for producing the results of the Zearalenone PBK model 
# # By: Chrysanthi Pachoulide
# # Date: 23-06-2025
# # --------------------------------------------------------------------------- #
# 
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
# 
CP_theme <- theme_minimal() +
  theme(
    axis.text = element_text(size = 56),
    axis.title = element_text(size = 62),
    plot.title = element_text(size = 65, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 56, hjust = 0.5),
    legend.position = "none",
    plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm")
  )

# load(here(OUTPUT, "RawData_PBK.RData"))
# save(PBK_OUTPUT, file = here(OUTPUT, "RawData_PBK.RData"))
# PBK_OUTPUT <- PBK.out %>% as.data.frame()

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
  # theme(axis.text = element_text(size = 10),
  #       axis.title = element_text(size = 12)
  # )+
  ylab("MB / ERROR")
print(MB_plot)

C_organs.df <- PBK_OUTPUT %>%
  transmute(time = time,
            SIntestineT_ZEN = CSI_ZEN,
            SIntestineL_ZEN = CSIL_ZEN,
            LIntestineL_ZEN = CLIL_ZEN,
            RectumL_ZEN = CREL_ZEN,
            Liver_ZEN = CL_ZEN,
            Fat_ZEN = CF_ZEN,
            Rapidly_ZEN = CR_ZEN,
            Slowly_ZEN = CS_ZEN,
            Blood_ZEN = CB_ZEN,
            LIntestineL_aZEL = CLIL_aZEL,
            SIntestineT_aZEL = CSI_aZEL,
            Liver_aZEL = CL_aZEL,
            Fat_aZEL = CF_aZEL,
            Rapidly_aZEL = CR_aZEL,
            Slowly_aZEL = CS_aZEL,
            Blood_aZEL = CB_aZEL,
            LIntestineL_bZEL = CLIL_bZEL,
            Liver_bZEL = CL_bZEL,
            SIntestineT_ZENGLU = CSI_ZENGLU,
            Liver_ZENGLU = CL_ZENGLU,
            Fat_ZENGLU = CF_ZENGLU,
            Rapidly_ZENGLU = CR_ZENGLU,
            Slowly_ZENGLU = CS_ZENGLU,
            Blood_ZENGLU = CB_ZENGLU,
            Liver_aZELGLU = CL_aZELGLU,
            SIntestineT_aZELGLU = CSI_aZELGLU
            )

C_organs_long <- C_organs.df %>%
  pivot_longer(names_to = c("Organ", "Compound"),
               SIntestineT_ZEN:SIntestineT_aZELGLU,
               values_to = "Concentration",
               names_sep = "_") %>%
  mutate(Concentration = Concentration*MWZEN)


Plasma <- C_organs_long %>%
  filter(Organ == "Blood") %>%
  ggplot()+
  geom_path(aes(x = time, y = Concentration, colour = Compound), linewidth = 0.5) +
  # scale_color_manual(name = "Model",
  #                    values = c("Predicted" = "blueviolet", "PredictedEFSA" = "darkblue"),
  #                    labels = c("Predicted" = "Predicted NG-PBK", "PredictedEFSA" = "Predicted EFSA 2020")
  # ) +
  # scale_linetype_manual(name = "Model",
  #                       values = c("Predicted" = "solid", "PredictedEFSA" = "dashed"),
  #                       labels = c("Predicted" = "Predicted NG-PBK", "PredictedEFSA" = "Predicted EFSA 2020")) +
  # guides(color = guide_legend(override.aes = list(linetype = c("solid", "dashed"))))+
  theme_minimal()+
  # CP_theme +
  ylab("Blood Concentration (ng/ml)") +
  xlab("Time (hours)") #+
  # theme(legend.position = c(0.9, 0.2),
  #       legend.text = element_text(size = 40),
  #       legend.title = element_blank(),
  #       legend.key.size = unit(0.5, "lines"),
  #       legend.spacing.y = unit(0, "cm"),
  #       legend.margin = margin(0, 0, 0, 0),
  #       strip.text = element_text(size = 45),
  #       # axis.text.x = element_text(margin = margin(t = 0, b = -10))
  # )
print(Plasma)
ggsave(
  filename = here(OUTPUT,"Blood.png"),
  dpi = 1000,
  width = 11,
  height = 7,
  units = "cm"
)


AllTissues <- C_organs_long %>%
  ggplot() +
  geom_path(aes(x = time, y = Concentration, colour = Compound), linewidth = 0.5) +
  facet_wrap(~Organ, scales = "free")+
  # scale_color_manual(name = "Model",
  #                    values = c("Predicted" = "blueviolet", "PredictedEFSA" = "darkblue"),
  #                    labels = c("Predicted" = "Predicted NG-PBK", "PredictedEFSA" = "Predicted EFSA 2020")
  # ) +
  # scale_linetype_manual(name = "Model",
  #                       values = c("Predicted" = "solid", "PredictedEFSA" = "dashed"),
  #                       labels = c("Predicted" = "Predicted NG-PBK", "PredictedEFSA" = "Predicted EFSA 2020")) +
  # guides(color = guide_legend(override.aes = list(linetype = c("solid", "dashed"))))+
  theme_minimal() +
  # CP_theme +
  ylab("Organ Concentration (ng/ml)") +
  xlab("Time (hours)") #+
  # theme(legend.position = "bottom", #c(0.9, 0.2),
  #       legend.text = element_text(size = 40),
  #       legend.title = element_blank(),
  #       legend.key.size = unit(0.5, "lines"),
  #       legend.spacing.y = unit(0, "cm"),
  #       legend.margin = margin(0, 0, 0, 0),
  #       strip.text = element_text(size = 45),
        # axis.text.x = element_text(margin = margin(t = 0, b = -10))
  # )
print(AllTissues)
# ggsave(
#   filename = here(OUTPUT,"AllTissues.png"),
#   dpi = 1000,
#   width = 11,
#   height = 7,
#   units = "cm"
# )




# Plot against in vivo ####

## ZEN ##
DataInvivo <- read.csv(here("Input", "InVivoData.csv"), header=1)
ObsPO <- DataInvivo %>%
  filter(AdminType == "PO",
         Gender != "Pregnant.F", #pregnancy not included yet
         TestCompound == "ZEN", #other options are ZEN_glu and a_ZEL
         Author == "Videmann", #Shin #Q.Lu #Videmann
         # Tissue == "Plasma",
         SampleCompound == "ZEN") %>%
  select(c("Author", "Gender","Tissue", "SamplingTime", "Concentration")) %>%
  rename(Time = SamplingTime)

# Choose tissue to compare
VidemannPO <- ObsPO %>%
  mutate(Concentration = ifelse(Tissue == "Plasma",
                               Concentration * 0.67, Concentration), #(ng/ml) blood concentration: changing from plasma concentrations to blood concentrations by multiplying with BP
  Tissue = ifelse(Tissue == "Plasma", "Blood", Tissue))

# rename(Blood = "Concentration")

# VidemannPO <- VidemannPO %>%
#   mutate(Concentration = ifelse(Tissue == "Plasma",
#                                 Concentration * 0.67, Concentration),
#          Tissue = ifelse(Tissue == "Plasma", "Blood", Tissue))

BW <- parm.c$BW  ##{Kg} #  body weight rat

VST <-parm.c$FV_ST*BW  # fraction volume of small intestine
VSIL <- parm.c$FV_SIL*BW # fractional volume of small intestinal lumen
VREL <- parm.c$FV_REL*parm.c$FV_LIL*BW # fractional volume of rectum lumen from total LI volume (only the lower part goes to the systemic circulation)
VLIL <- (parm.c$FV_LIL*BW)-VREL
VSI <- parm.c$FV_SI*BW # fractional volume of small intestinal tissue


VidemannPO_wide <- VidemannPO %>%
  pivot_wider(names_from = Tissue, values_from = Concentration) %>%
  select(c("Time", "Gender", "Blood", "Intestine", "Liver"))

VidemannPO_long <- VidemannPO_wide %>%
pivot_longer(Blood:Liver, names_to = "Tissue", values_to = "Concentration")

# Observed
C_organs_wide <- C_organs_long %>%
  pivot_wider(names_from = "Organ", values_from = "Concentration") %>%
  filter(Compound == "ZEN") %>%
  select(c("time", "SIntestineT", "SIntestineL", "LIntestineL", "RectumL", "Liver", "Blood")) %>% 
  mutate("IntestineL" = (((SIntestineL*VSIL) + (LIntestineL*VLIL) + (RectumL*VREL))/(VSIL+VLIL+VREL))) %>% 
  # mutate("TIntestine" = (((SIntestineT*VSI) + (IntestineL*(VSIL+VLIL+VREL)))/(VSI + VSIL+VLIL+VREL)))%>% 
  rename("Intestine" = "SIntestineT") %>% 
  select(c("time", "Intestine", "Liver", "Blood"))

PredPO <- C_organs_wide %>%
  rename(Concentration = Blood) # %>%
  # filter(time >= 0.05) #%>% #0.5 for Videmann, 0.05 for Q.Lu

Videmann.PredPO <- C_organs_wide %>%
  pivot_longer(Intestine:Blood, names_to = "Tissue", values_to = "Concentration") # %>%
  # filter(time >= 0.5) #%>% #0.5 for Videmann, 0.05 for Q.Lu
# filter(time <= 5)

# Plots
Blood.PredObsVidemannPlotly <- ggplot() +
  geom_point(data = filter(VidemannPO, Tissue == "Blood"), aes(Time, Concentration)) + #, color = Gender)) +
  # facet_wrap(~Tissue) +
  geom_path(data = PredPO, aes(time, Concentration)) +
  # facet_wrap(~Tissue) +
  labs(title = "Zearalenone blood concentrations",
       x = "Time (h)",
       y = "Concentration (ng/ml)",
       color = "Author")+
  theme_minimal()
  # CP_theme
# PredObsVidemannPlotly
# Blood.PredObsVidemannPlotly <- ggplotly(Blood.PredObsVidemannPlotly)
print(Blood.PredObsVidemannPlotly)
# ggsave(
#   filename = here(OUTPUT,"Blood.PredObsVidemannPlotly.png"),
#   dpi = 1000,
#   width = 11,
#   height = 7,
#   units = "cm"
# )

PredObsVidemann <- ggplot() +
  geom_point(data = VidemannPO_long, aes(Time, Concentration, color = Gender)) +
  facet_wrap(~Tissue)+
  geom_line(data = Videmann.PredPO, aes(time, Concentration))+
  facet_wrap(~Tissue)+
  labs(title = "Zearalenone blood concentrations",
       x = "Time (h)",
       y = "Concentration (ng/ml)",
       color = "Author") +
  theme_minimal()
  # CP_theme +
# theme(legend.position = c(0.9, 0.9),
#       legend.text = element_text(size = 40),
#       legend.title = element_blank(),
#       legend.key.size = unit(0.5, "lines"),
#       legend.spacing.y = unit(0, "cm"),
#       legend.margin = margin(0, 0, 0, 0),
#       strip.text = element_text(size = 45),
#       # axis.text.x = element_text(margin = margin(t = 0, b = -10))
# )
print(PredObsVidemann)
# ggsave(
#   filename = here(OUTPUT,"PredObs.PO.Videmann.png"),
#   dpi = 1000,
#   width = 11,
#   height = 7,
#   units = "cm"
# )


Blood.PredObsVidemann <- ggplot() +
  geom_point(data = filter(VidemannPO, Tissue == "Blood"), aes(Time, Concentration))+ #, color = Gender)) +
  # facet_wrap(~Tissue)+
  geom_line(data = PredPO, aes(time, Concentration))+
  # facet_wrap(~Tissue)+
  labs(title = "Zearalenone blood concentrations",
       x = "Time (h)",
       y = "Concentration (ng/ml)",
       color = "Author") +
  theme_minimal()
# CP_theme +
#   theme(legend.position = c(0.9, 0.9),
#         legend.text = element_text(size = 40),
#         legend.title = element_blank(),
#         legend.key.size = unit(0.5, "lines"),
#         legend.spacing.y = unit(0, "cm"),
#         legend.margin = margin(0, 0, 0, 0),
#         strip.text = element_text(size = 45),
#         # axis.text.x = element_text(margin = margin(t = 0, b = -10))
#   )
print(Blood.PredObsVidemann)
# ggsave(
#   filename = here(OUTPUT,"PredObs.PO.VidemannBlood.png"),
#   dpi = 1000,
#   width = 11,
#   height = 7,
#   units = "cm"
# )



## ZEN0GLUC ##
ObsPO_GLUC <- DataInvivo %>%
  filter(AdminType == "PO",
         Gender != "Pregnant.F", #pregnancy not included yet
         TestCompound == "ZEN", #other options are ZEN_glu and a_ZEL
         Author == "Videmann", #Shin #Q.Lu #Videmann
         # Tissue == "Plasma",
         SampleCompound == "ZEN_glu") %>%
  select(c("Author", "Gender","Tissue", "SamplingTime", "Concentration")) %>%
  rename(Time = SamplingTime)

VidemannPO_GLUC <- ObsPO_GLUC %>%
  mutate(#Concentration = ifelse(Tissue == "Plasma",
    #                      Concentration * 0.67, Concentration), #(ng/ml) blood concentration: changing from plasma concentrations to blood concentrations by multiplying with BP
    Tissue = ifelse(Tissue == "Plasma", "Blood", Tissue))



VidemannPO_wide_GLUC <- VidemannPO_GLUC %>%
  pivot_wider(names_from = Tissue, values_from = Concentration) %>%
  select(c("Time", "Gender", "Blood", "Intestine", "Liver"))

VidemannPO_long_GLUC <- VidemannPO_wide_GLUC %>%
  pivot_longer(Blood:Liver, names_to = "Tissue", values_to = "Concentration")

# Observed
C_organs_wide_GLUC <- C_organs_long %>%
  pivot_wider(names_from = "Organ", values_from = "Concentration") %>%
  filter(Compound == "ZENGLU") %>%
  select(c("time", "SIntestineT", "SIntestineL", "LIntestineL", "RectumL", "Liver", "Blood")) %>% 
  mutate("IntestineL" = (((SIntestineL*VSIL) + (LIntestineL*VLIL) + (RectumL*VREL))/(VSIL+VLIL+VREL))) %>% 
  # mutate("TIntestine" = (((SIntestineT*VSI) + (IntestineL*(VSIL+VLIL+VREL)))/(VSI + VSIL+VLIL+VREL)))%>% 
  rename("Intestine" = "SIntestineT") %>% 
  select(c("time", "Intestine", "Liver", "Blood"))

PredPO_GLUC <- C_organs_wide_GLUC %>%
  rename(Concentration = Blood) # %>%
# filter(time >= 0.05) #%>% #0.5 for Videmann, 0.05 for Q.Lu

Videmann.PredPO_GLUC <- C_organs_wide_GLUC %>%
  pivot_longer(Intestine:Blood, names_to = "Tissue", values_to = "Concentration") # %>%
# filter(time >= 0.5) #%>% #0.5 for Videmann, 0.05 for Q.Lu
# filter(time <= 5)

# Plots
Blood.PredObsVidemannPlotly_GLUC <- ggplot() +
  geom_point(data = filter(VidemannPO_GLUC, Tissue == "Blood"), aes(Time, Concentration)) + #, color = Gender)) +
  # facet_wrap(~Tissue) +
  geom_path(data = PredPO_GLUC, aes(time, Concentration)) +
  # facet_wrap(~Tissue) +
  labs(title = "Zearalenone glucuronide blood concentrations",
       x = "Time (h)",
       y = "Concentration (ng/ml)",
       color = "Author")+
  theme_minimal()
# CP_theme
# PredObsVidemannPlotly
# Blood.PredObsVidemannPlotly <- ggplotly(Blood.PredObsVidemannPlotly)
print(Blood.PredObsVidemannPlotly_GLUC)
# ggsave(
#   filename = here(OUTPUT,"Blood.PredObsVidemannPlotly.png"),
#   dpi = 1000,
#   width = 11,
#   height = 7,
#   units = "cm"
# )

PredObsVidemann_GLUC <- ggplot() +
  geom_point(data = VidemannPO_long_GLUC, aes(Time, Concentration, color = Gender)) +
  facet_wrap(~Tissue)+
  geom_line(data = Videmann.PredPO_GLUC, aes(time, Concentration))+
  facet_wrap(~Tissue)+
  labs(title = "Zearalenone gluc blood concentrations",
       x = "Time (h)",
       y = "Concentration (ng/ml)",
       color = "Author") +
  theme_minimal()
# CP_theme +
# theme(legend.position = c(0.9, 0.9),
#       legend.text = element_text(size = 40),
#       legend.title = element_blank(),
#       legend.key.size = unit(0.5, "lines"),
#       legend.spacing.y = unit(0, "cm"),
#       legend.margin = margin(0, 0, 0, 0),
#       strip.text = element_text(size = 45),
#       # axis.text.x = element_text(margin = margin(t = 0, b = -10))
# )
print(PredObsVidemann_GLUC)
# ggsave(
#   filename = here(OUTPUT,"PredObs.PO.Videmann.png"),
#   dpi = 1000,
#   width = 11,
#   height = 7,
#   units = "cm"
# )

# ObsIV_plot <- ObsPO %>% ggplot() +
#   geom_path(aes(Time, Concentration, color = Author)) +
#   facet_wrap(~Gender) +
#   labs(title = "Experimental Zearalenone blood concentrations",
#        x = "Time (h)",
#        y = "Concentration (ng/ml)",
#        color = "Author") +
#   theme_minimal()
# CP_theme +
#   theme(legend.position = c(0.9, 0.9),
#         legend.text = element_text(size = 40),
#         legend.title = element_blank(),
#         legend.key.size = unit(0.5, "lines"),
#         legend.spacing.y = unit(0, "cm"),
#         legend.margin = margin(0, 0, 0, 0),
#         strip.text = element_text(size = 45),
#         # axis.text.x = element_text(margin = margin(t = 0, b = -10))
#   )
# ggsave(
#   filename = here("ObservedIV.png"),
#   dpi = 1000,
#   width = 11,
#   height = 7,
#   units = "cm"
# )
#
# ObsPOPlotly <- ggplotly(ObsIV_plot)
# ObsPOPlotly
#
#
#
# ObsPred <- ggplot() +
#   geom_path(ObsPO, aes(Time, Concentration, color = Author)) +
#   facet_wrap(~Gender) +
#   labs(title = "Experimental Zearalenone blood concentrations",
#        x = "Time (h)",
#        y = "Concentration (ng/ml)",
#        color = "Author") +
#   CP_theme +
#   theme(legend.position = c(0.9, 0.9),
#         legend.text = element_text(size = 40),
#         legend.title = element_blank(),
#         legend.key.size = unit(0.5, "lines"),
#         legend.spacing.y = unit(0, "cm"),
#         legend.margin = margin(0, 0, 0, 0),
#         strip.text = element_text(size = 45),
#         # axis.text.x = element_text(margin = margin(t = 0, b = -10))
#   )