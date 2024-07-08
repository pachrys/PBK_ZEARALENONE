
#*****************************************************************************#
#---Results-Data Analysis - against only Videmann data------------------------#
#*****************************************************************************#

#source(here("Script/Functions/PlottingFunction.R"))                               #call plotting function; if not already used
# WhatAmITesting <- "CorrectedGLUCPCs"
# TodaysOuput <- here("Output/Data", Sys.Date(), WhatAmITesting)
# PlotFolderIs <- here("Output/Plots", Sys.Date(), WhatAmITesting)
# dir.create(TodaysOuput, WhatAmITesting,recursive = TRUE)

#----Plotting only experimental-----------------------------------------------#

ZEN.dat <- read.csv(here("InputData/Zen.input.csv"), header=1)     #from {Agahi, 2020 #456} and chemicalize, BP from Shin, Hong et al.2009

# Input in vivo data
DataInvivo.ZEN <- read.csv(here("InputData/InVivoData.csv"), header=1) %>% 
  filter(Author == "Videmann", #I chose Videmann as it had the best recoveries, also Videmann separates Zearalenone from it's metabolites
         Gender != "Pregnant.F", #pregnancy not included yet
         TestCompound == "ZEN", #other options are ZEN_glu and a_ZEL
         AdminType == "PO", #"IV"), #Videmann doesn't have IV, but subcutaneous!!
         Tissue %in% c("Plasma", "Intestine", "Liver", "Adipose", "Kidney", "Urine", "Colon"), 
         SampleCompound == "ZEN") %>%  #only for ZEN
  select(c("Gender", "SamplingTime", "Tissue", "Concentration"))%>% 
  rename(Time = SamplingTime) %>% #(h)
  mutate(Concentration = Concentration*ZEN.dat$BP) %>% #(ng/ml) blood concentration: changing from plasma concentrations to blood concentrations by multiplying with BP
  mutate(Concentration = Concentration/0.76) %>%  #correcting for the recovery (recovery was 72-80% => averaged that to 76% and then corrected the concentrations for that)
  mutate(Tissue = recode(Tissue, "Plasma" = "Blood")) #renaming plasma to blood 

# To create a file for either PO or IV admin
POData <- DataInvivo.ZEN 

# Input experimental data 
CONCENTRATIONS.ZEN <- CONCENTRATIONS %>% 
  filter(!str_ends(Tissue, "_glucuronide"))

# # To plot only In Vivo data
# 
# POData %>% 
#   #filter(Tissue %in% c("Blood")) %>% #choose a specific tissue if needed
#   ggplot(aes(Time, Concentration, color = Tissue)) +
#     facet_wrap(~Gender) +
#     geom_point() +
#     geom_line(linewidth = 0.1) + 
#     labs(title = "Experimental zearalenone blood concentrations",
#          x = "Time (h)",
#          y = "Concentration (ng/ml)",
#          color = "Tissue") +
#     CP_theme()
# ggsave(here(PlotFolderIs, "Invivo.png"), dpi = 300)



#----Simulated vs Experimental------------------------------------------------#

# Creating dataframes with the same tissues for both in vivo and simulated
Invivo <- POData %>% 
  filter(Tissue %in% c("Blood",
                       "Intestine", "Liver", "Kidney")) 
  
Simulated <- CONCENTRATIONS.ZEN %>%                                                         
  filter(Tissue %in% c("Blood",                                                  
                       "Intestine", "Liver", "Kidney")) 

# Simulated VS Experimental all tissue AUC 
Simulated %>% # choose compartment accordingly
  ggplot() +
  geom_point(data = Invivo, aes(Time, Concentration, color = Gender)) + 
  geom_line(data = Invivo, aes(Time, Concentration, color = Gender), 
            linewidth = 0.1) + 
  facet_wrap(~Tissue) +
  geom_path(aes(Time, Concentration, color = "black")) +
  facet_wrap(~Tissue) +
  #scale_colour_manual(values = Chrysa_colors) +
    labs(title = "Zearalenone tissue concentrations",
         x = "Time (h)",
         y = "Concentration (ng/ml)",
         color = "Gender") +
  CP_theme()
ggsave(here(PlotFolderIs, "AllTissueSimulatedInvivo.png"), dpi = 300)

# Separate organs
# Blood
InvivoBlood <- POData %>% filter(Tissue =="Blood") 

SimulatedBlood <- CONCENTRATIONS.ZEN %>% filter(Tissue =="Blood") 

# Simulated VS Experimental all tissue AUC 
SimulatedBlood %>% # choose compartment accordingly
  ggplot() +
  geom_point(data = InvivoBlood, aes(Time, Concentration, color = Gender)) + 
  geom_line(data = InvivoBlood, aes(Time, Concentration, color = Gender), 
            linewidth = 0.1) + 
  geom_path(aes(Time, Concentration)) +
  #scale_colour_manual(values = Chrysa_colors) +
  labs(title = "Zearalenone blood concentrations",
       x = "Time (h)",
       y = "Concentration (ng/ml)",
       color = "Gender") +
  CP_theme()
ggsave(here(PlotFolderIs, "BloodSimulatedInvivo.png"), dpi = 300)

# Intestine
InvivoIntestine <- POData %>% filter(Tissue =="Intestine") 

SimulatedIntestine <- CONCENTRATIONS.ZEN %>% filter(Tissue =="Intestine") 

# Simulated VS Experimental all tissue AUC 
SimulatedIntestine %>% # choose compartment accordingly
  ggplot() +
  geom_point(data = InvivoIntestine, aes(Time, Concentration, color = Gender)) + 
  geom_line(data = InvivoIntestine, aes(Time, Concentration, color = Gender), 
            linewidth = 0.1) + 
  geom_path(aes(Time, Concentration, color = "black")) +
  #scale_colour_manual(values = Chrysa_colors) +
  labs(title = "Zearalenone Intestine concentrations",
       x = "Time (h)",
       y = "Concentration (ng/ml)",
       color = "Gender") +
  CP_theme()
ggsave(here(PlotFolderIs, "IntestineSimulatedInvivo.png"), dpi = 300)

# Liver
InvivoLiver <- POData %>% filter(Tissue =="Liver") 

SimulatedLiver <- CONCENTRATIONS.ZEN %>% filter(Tissue =="Liver") 

# Simulated VS Experimental all tissue AUC 
SimulatedLiver %>% # choose compartment accordingly
  ggplot() +
  geom_point(data = InvivoLiver, aes(Time, Concentration, color = Gender)) + 
  geom_line(data = InvivoLiver, aes(Time, Concentration, color = Gender), 
            linewidth = 0.1) + 
  geom_path(aes(Time, Concentration, color = "black")) +
  #scale_colour_manual(values = Chrysa_colors) +
  labs(title = "Zearalenone Liver concentrations",
       x = "Time (h)",
       y = "Concentration (ng/ml)",
       color = "Gender") +
  CP_theme()
ggsave(here(PlotFolderIs, "LiverSimulatedInvivo.png"), dpi = 300)

# Kidney
InvivoKidney <- POData %>% filter(Tissue =="Kidney") 

SimulatedKidney <- CONCENTRATIONS.ZEN %>% filter(Tissue =="Kidney") 

# Simulated VS Experimental all tissue AUC 
SimulatedKidney %>% # choose compartment accordingly
  ggplot() +
  geom_point(data = InvivoKidney, aes(Time, Concentration, color = Gender)) + 
  geom_line(data = InvivoKidney, aes(Time, Concentration, color = Gender), 
            linewidth = 0.1) + 
  geom_path(aes(Time, Concentration, color = "black")) +
  #scale_colour_manual(values = Chrysa_colors) +
  labs(title = "Zearalenone Kidney concentrations",
       x = "Time (h)",
       y = "Concentration (ng/ml)",
       color = "Gender") +
  CP_theme()
ggsave(here(PlotFolderIs, "KidneySimulatedInvivo.png"), dpi = 300)

#*****************************************************************************#
#------------------------------Data analysis----------------------------------#
# 
# Rearranging the data--------------------------------------------------------#

# Videmann
VidemannSimple.F.data <- POData %>% #female
  filter(Tissue == "Blood", Gender == "Female") %>%
  select(Time, Concentration) %>%
  mutate(Time = ifelse(Time == 0.7, 0.705, Time)) %>% #original timepoint was 0.083, however I'm simulating by 0.005 so I'm changing it to make it comparable
  mutate(Concentration = Concentration/0.7)

VidemannSimple.M.data <- POData %>% #male
  filter(Tissue == "Blood", Gender == "Male") %>%
  select(Time, Concentration) %>%
  mutate(Time = ifelse(Time == 0.7, 0.705, Time)) #original timepoint was 0.083, however I'm simulating by 0.005 so I'm changing it to make it comparable

Videmann.TimeP <- VidemannSimple.F.data$Time #time points are the same for both male and female
SimVidemann.data <- Simulated %>% 
  filter(Tissue == "Blood",
         Time %in% Videmann.TimeP) #simulated data is the same for both genders as gender differences are not included


#Calculating AUCs------------------------------------------------------------#

gender <- "female" #choose gender "female" | "male" 

if(gender == "female"){ #selecting gender specific data from the invivo dataframe
  VidemannSimple.data <- VidemannSimple.F.data
} else if(gender == "male"){
  VidemannSimple.data <- VidemannSimple.M.data
}

Videmann.AUC <- 0 #creating empty vector for AUC values from experimental data; unit: ng.h/ml
  for (i in 2:nrow(VidemannSimple.data)) {
    Videmann.AUC <- Videmann.AUC +
      ((VidemannSimple.data$Concentration[i-1] + VidemannSimple.data$Concentration[i]) / 2) *
      (VidemannSimple.data$Time[i] - VidemannSimple.data$Time[i-1])
  }

SimVidemann.AUC <- 0 #creating empty vector for AUC values from simulated data; unit: ng.h/ml
  for (i in 2:nrow(SimVidemann.data)) {
    SimVidemann.AUC <- SimVidemann.AUC +
      ((SimVidemann.data$Concentration[i-1] + SimVidemann.data$Concentration[i]) / 2) *
      (SimVidemann.data$Time[i] - SimVidemann.data$Time[i-1])
  } 

if(gender == "female"){ #renaming depending on gender
    Videmann.F.AUC <- Videmann.AUC
  } else if(gender == "male"){
    Videmann.M.AUC <- Videmann.AUC
  }

Combined.df <- bind_cols( #creating df with combined data
  VidemannSimple.F.data, SimVidemann.data$Concentration) %>%
  rename(Experimental = Concentration, ,
         Simulated = ...3)


# Simulated/Experimental linear regression------------------------------------#
Regression <- lm(Experimental~Simulated, data=Combined.df)
summary(Regression)

Regression %>%
  ggplot(aes(Simulated, Experimental)) +
  geom_point(aes(Simulated, Experimental)) +
  geom_smooth(method='lm', color = "black") +
  CP_theme() +
  theme(legend.position = "bottom")
ggsave(here(PlotFolderIs, "ExpVsSimRegression.png"), dpi = 300)


#---End-----------------------------------------------------------------------#
#*****************************************************************************#
