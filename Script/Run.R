
# library(rxode2)
# library(tidyverse)
# library(ggplot2)
# library(here)
# library(showtext)
# font_add(family = "Garamond", regular = "GARA.TTF")
# showtext_auto()

#*****************************************************************************#

TodayIs <- format(Sys.Date())
 
#Define what you want to test! choose between PO, PO_EHC, IV, IV_EHC this will define the function that will be chosen
AdminType <- "PO_EHC"
WhatAmITesting <- "EHCNotRecoveryCorrected"


#*****************************************************************************#
#---Run script----------------------------------------------------------------#
#*****************************************************************************#

#----Extracting input---------------------------------------------------------#

source(here("Script/CalculatingInput.R"))                                         #call input script

# ZEN physicochemical input                                                       #ZEN as monoprotic acid. It has 2 pKas 8.54 and 10.38. At pH 7.4 ZEN is >90 percent not ionized.
ZEN.dat <- read.csv(here("InputData/ZearalenoneSpecificInput.csv"), header=1)     #from {Agahi, 2020 #456} and chemicalize, BP from Shin, Hong et al.2009

# In Vivo study input
PO.dat <- read.csv(here("InputData/InVivoData.csv"), header=1) %>% 
  filter(Author == "Videmann", #I chose Videmann as it had the best recoveries, also Videmann separates Zearalenone from it's metabolites
           Gender != "Pregnant.F", #pregnancy not included yet
         TestCompound == "ZEN", #other options are ZEN_gluc and a_ZEL
         AdminType == "PO", #"IV"), #Videmann doesn't have IV, but subcutaneous!!
         Tissue %in% c("Plasma", "Intestine", 
                       "Liver", "Adipose", "Kidney", "Urine", "Colon")) %>% 
  select(c("BW", "PODose")) %>%                                                   #"Author", "Gender" are the same as I now only use one study
  distinct() %>% 
  print()

#-----------------------------------------------------------------------------#


#**STEP 1**#
#*****************************************************************************#
#------------------------------Model Parameters-------------------------------#


#---Animal/Study specific parameters------------------------------------------#
BW = PO.dat$BW;                                                                   #(L)
QC = 80*60/1000;                                                                  #(L/h) units change; 80ml/min; Musther

# Tissue Volumes (L)                                                              #(these are weights; assumption that every tissue has a density of 1 kg/L); fractional volumes come from Musther (insignificant differences from Brown, Punt, Ruark)
#VSt = BW*CalcInput$FVSt;                                                          #stomach
VIn = BW*CalcInput$FVGu;                                                          #intestine
VLi = BW*CalcInput$FVLi;                                                          #liver
VKi = BW*CalcInput$FVKi;                                                          #kidney
VAd = BW*CalcInput$FVAd;                                                          #adipose
VSl = BW*(CalcInput$FVBo + CalcInput$FVBr + CalcInput$FVMu + CalcInput$FVSk);     #slowly perfused (as above)
VRa = BW*(CalcInput$FVHe + CalcInput$FVLu + CalcInput$FVSp);                      #rapidly perfused (as above)
VBl = BW*0.06;                                                                    #0.06 is the fractional volume of blood based on Musther 
VBi = BW*0.00023495; #bile volume total volume was 46.97 +- 4.21ul; male Fisher rats 180-220g => fractional volume of bile is 0.00023495 https://ncbi.nlm.nih.gov/pmc/articles/PMC1891969/
VInL = 2/1000; #total volume of lumen = 2ml 

# Tissue Blood Flow Rates (L/h)
QC = 80*60/1000;                                                                  #(L/h) units change; 80ml/min; Musther
QIn = QC*CalcInput$FQGu;
QLi = QC*CalcInput$FQLi;                                                          #(L/h) liver
QKi = QC*CalcInput$FQKi;                                                          #(L/h) kidney
QAd = QC*CalcInput$FQAd;                                                          #(L/h) adipose
QSl = QC*(CalcInput$FQBo + CalcInput$FQBr +
            CalcInput$FQMu + CalcInput$FQSk);                                     #(L/h) slowly perfused (as above)
QRa = QC - (QLi + QKi + QAd + QSl);                                               #(L/h) rapidly perfused (as above! without lungs as lungs get everything!)
MassBalanceQC = QC - (QLi + QKi + QAd + QSl + QRa);                               #check blood flow mass balance; should be 0

QGu = QC*CalcInput$FQGu

# Intestinal transit parameters
# kTge = 1/ZEN.dat$kTge, #1/h, gastric emptying rate
# kTsi = 1/ZEN.dat$kTsi.7, #1/h, small intestine transit time; change to kTsi or kTsi.7 depending on the number of intestines
# kTco = 1/ZEN.dat$kTco, #1/h, colon transit time
kTIn = 1/(ZEN.dat$kTsi+ZEN.dat$kTco) #1/h, total intestinal transit time including small and large intestine

# TO CHANGE!!! Scaling in vitro to in vivo absorption parameters 
Peff.human = (10^(0.4926*log10(ZEN.dat$Papp)- 0.1454)); #human Peff(cm^-4/s)
Peff.rat = Peff.human/11.04; #rat Peff(cm^-4/s)
Peff = Peff.rat/10000*3600; #conversion to cm/h
ka = Peff*2/ZEN.dat$R; #/h  

# Clearance parameters
CL_Bi = ZEN.dat$CL_Bi*60/1000 #L/h, biliary clearance; scaling from ml/min to L/h; CHANGED FROM 0.4 TO 6.3; BASED ON BAOJIAN WU DOI 10.1002/jps.22827
#kBi = CL_Bi/VLi #/h. biliary secretion rate
CL_Re = ZEN.dat$CL_Re*60/1000 #L/h, renal clearance; scaling to L/h


# Scaling in vitro to in vivo metabolism parameters
VmaxMet.Li = ZEN.dat$VmaxMet.Li * ZEN.dat$SF.Li * VLi*1000 * 60/1000 #umol/h, scaling from nmol/min/mg protein to nmol/min for whole liver; equation from {Grandoni, 2019 #457}; in vivo = in vitro * scaling to microsomal protein per gram of rat liver * scaling to total weight (assuming density=1) of rat liver; SF = 46mg microsomal protein/g liver ({Musther, 2017 #458})
VmaxMet.In = ZEN.dat$VmaxMet.In * ZEN.dat$SF.In * VIn*1000 * 60/1000 #umol/h, scaling from nmol/min/mg protein to nmol/min for whole liver; equation from {Grandoni, 2019 #457}; in vivo = in vitro * scaling to microsomal protein per gram of rat liver * scaling to total weight (assuming density=1) of rat intestine; SF = 15.5mg microsomal protein/g gut ({Musther, 2017 #458})

#Vmax_hyd = ZEN.dat$Vmax_hyd * ZEN.dat$SF * VIn*1000 #scaling from nmol/min/mg protein to nmol/min for whole liver; equation from {Grandoni, 2019 #457}; in vivo = in vitro * scaling to microsomal protein per gram of rat liver * scaling to total weight (assuming density=1) of rat liver
#kHydro = 1/1000 #/h according to Veronique

# Scaling in vitro to in vivo transport parameters; according to {Videmann, 2008 #384}; input values are expressed in pmol/h and uM; cell concentration for transport experiments is not mentioned!!
VmaxExc.In.glu = ZEN.dat$VmaxExc.In.glu/1000000 #umol/h, scaling from pmol/h
VmaxAbs.In.glu = ZEN.dat$VmaxAbs.In.glu/1000000 #umol/h, scaling from pmol/h

params <- cbind(
  PLib = CalcInput$PLib, #liver:plasma partition coefficient
  PKib = CalcInput$PKib, #kidney:plasma partition coefficient
  PAdb = CalcInput$PAdb, #adipose:plasma partition coefficient
  PSlb = CalcInput$PSlb, #slowlyperfusedtissue:plasma partition coefficient
  PRab = CalcInput$PRab, #rapidlyperfusedtissue:plasma partition coefficient
  fub = CalcInput$fub, #fraction unbound in blood
  VLi = VLi, #L, liver volume
  VIn = VIn, #L, intestinal volume
  VKi = VKi, #L, kidney volume
  VAd = VAd, #L, adipose volume
  VSl = VSl, #L, slowlyperfusedtissue volume
  VRa = VRa, #L. rapidlyperfusedtissue volume
  VBl = VBl, #L, blood volume
  QLi = QLi, #L/h, liver flow rate
  QKi = QKi, #L/h, kidney flow rate
  QAd = QAd, #L/h, adipose flow rate
  QSl = QSl, #L/h, slowlyperfusedtissue flow rate
  QRa = QRa, #L/h, rapidlyperfusedtissue flow rate
  ka = ka, #absorption constant
  kTIn = kTIn, #1/h, total intestinal transit time including small and large intestine
  CL_Re = CL_Re, #L/h, renal clearance; scaling to L/h
  CL_Bi = CL_Bi,
  VmaxMet.Li = VmaxMet.Li, #umol/h, liver maximum rate of metabolism; scaling from nmol/min to umol/h
  VmaxMet.In = VmaxMet.In, #umol/h, intestine maximum rate of metabolism; scaling from nmol/min to umol/h
  Km_glu = ZEN.dat$Km_glu, #uM, enzyme affinity constant; proper to the enzyme so same for both liver and intestine
  KmExc.glu = ZEN.dat$KmExc.glu, #uM, transporter affinity constant; excretion from intestinal tissue to lumen
  KmAbs.glu = ZEN.dat$KmAbs.glu #uM, transporter affinity constant; absorption from lumen to intestinal tissue
)
#-----------------------------------------------------------------------------#


#**STEP 2**#
#*****************************************************************************#
#-------------------------Model initiation information------------------------#
amount.units = "umol"
time.units = "h"
nbr.doses = 1
sampling.time.from = 0
sampling.time.to = 24
sampling.frequency = 1/200

#-----Dose--------------------------------------------------------------------#
PO_DOSE.mgkg = PO.dat$PODose;                                                     #Videmann (mg/kg) dose; to change depending on the study
PO_DOSE = PO_DOSE.mgkg*BW/CalcInput$MW*1000;                                      #(umol) unit conversion from mg to umol
  
#-----Define initial state of compartments------------------------------------#
inits <- c(
    ASt = 0,
    AIn = 0,
    ACo = 0,
    AFe = 0, 
    ALi = 0,
    ALi.glu = 0, 
    AKi = 0,
    AUr = 0,
    AAd = 0,
    ASl = 0,
    ARa = 0,
    ABl = 0)


#**STEP 3**#
#*****************************************************************************#
#-----------------------------Solving the model-------------------------------#

# Call PBK functions
source(here("Script/Functions/ModelFunctions.R"))

ModelZEN <- PO_EHC.PBKfunction() #!choose the PBK function accordingly!
print(summary(ModelZEN))

# Event table
events <- eventTable(amount.units = amount.units,
                       time.units = time.units) %>%
    add.dosing(dose = PO_DOSE, dur = 0.01, nbr.doses = nbr.doses,
               dosing.to = "ASt") %>%   #!change accordingly!              
    add.sampling(seq(from = sampling.time.from, 
                     to = sampling.time.to, 
                     by = sampling.frequency))

# Saving name
Output.title <- paste("Simulated.asPO") #!Change accordingly!

Output.Model <- rxSolve(ModelZEN, 
                        params = params, 
                        events = events,
                        inits = inits) %>% 
  as.data.frame()

DataFolderTodayIs <- here("Output/Data", TodayIs, AdminType, WhatAmITesting)      #folder for saving results
dir.create(DataFolderTodayIs, recursive = TRUE)

 

#**STEP 4**#
#*****************************************************************************#
# Plotting and saving code
source(here("Script/Functions/PlottingFunction.R"))                               #call plotting function
Chrysa_colors <- c("purple3",                                                     #choosing colours for visualisation
                   "springgreen3",
                   "royalblue2",
                   "yellow3",
                   "orange",
                   "violetred4",
                   "palegreen4",
                   "lightskyblue",
                   "tomato2")
PlotFolderIs <- here("Output/Plots", TodayIs, AdminType, WhatAmITesting)          #folder for saving results
dir.create(PlotFolderIs, recursive = TRUE)                                        #create directory only if needed


#----------------------Check Mass Balance and Error---------------------------#

Calculated <- c( #!Change accordingly depending on the compartments in the function#
    "ASt",
    "AIn",
    "ACo",
    "AFe",
    "ALi",
    "ALi.glu",
    "AKi",
    "AUr",
    "AAd",
    "ASl",
    "ARa",
    "ABl")

#-----Mass Balance Calculations----------------------------------------------#
Calculated <- rowSums(Output.Model[ ,Calculated])
Total <- PO_DOSE

Massbalance <- Total - Calculated + 1
ErrorAll <- (Total - Calculated) / (Total + 1*10^(-30)) * 100                     #mass error (% of mass lost)

MassBalance_Error <- cbind("time" = Output.Model[, "time"],                       #data frame
                           Massbalance, ErrorAll) %>%
  as.data.frame()

# Mass Balance and Error plot
MassBalance_Error %>%
  filter(time >= 0.01) %>% #as PO administration time is 0,01h
  ggplot(aes(time, Massbalance)) +
  geom_line(color = "grey") +
  geom_line(aes(time, ErrorAll), color = "grey6") +
  labs(x = "Time (h)", y = "Value",
       caption = "Massbalance in grey, Error in black") +
  CP_theme()
ggsave(here(PlotFolderIs, "MassbalanceError.png"), dpi = 300)             #saving the plot


#**STEP 5**#
#****************************************************************************#
#------------------------Plotting Simulated Output---------------------------#


# To create df of simulated data for plotting
PlotConc <- Output.Model %>% #!change accordingly depending on the compartments in the model
  rename("Time" = "time",
         "Intestine" = "CIn",
         "Liver" = "CLi",
         "Liver_glu" = "CLi.glu",
         "Kidney" = "CKi",
         "Adipose" = "CAd",
         "PoorlyPerfused" = "CSl",
         "RapidlyPerfused" = "CRa",
         "Blood" = "CBl",
         "LiverVenous" = "CvLi",
         "KidneyVenous" = "CvKi",
         "AdiposeVenous" = "CvAd",
         "PoorlyPerfusedVenous" = "CvSl",
         "RapidlyPerfusedVenous" = "CvRa") %>%
  select(Time, Intestine:RapidlyPerfusedVenous) %>% #!change accordingly depending on the compartments in the model
  pivot_longer(names_to = "Tissue",
               values_to = "Concentration",
               Intestine:RapidlyPerfusedVenous,
  ) %>%
  mutate(Concentration = Concentration*ZEN.dat$MW)  #converting umol/L to ng/ml

# To plot only simulated data
PlotConc %>%
  filter(Tissue %in% c("Liver", "Blood", "Kidney", "Intestine",
                       "Adipose", "PoorlyPerfused", "RapidlyPerfused")) %>%
  ggplot(aes(Time, Concentration, color = Tissue)) +
  geom_line() +
  labs(title = "Simulated zearalenone tissue concentrations",
       x = "Time (h)",
       y = "Concentration (ng/ml)",
       caption = paste("Administration type:", AdminType, "Dose:", PO_DOSE.mgkg, "mg/kg bw")) +
  CP_theme() +
  scale_colour_manual(values=Chrysa_colors)
ggsave(here(PlotFolderIs, "Simulated.png"), dpi = 300)           

# To plot blood concentrations
PlotConc %>%
  filter(Tissue == "Blood") %>% 
  ggplot(aes(Time, Concentration)) +
  geom_path(linewidth = 0.5) +
  labs(title = "Zearalenone blood concentrations",
       x = "Time (h)",
       y = "Concentration (ng/ml)",
       caption = paste(AdminType)) +
  CP_theme() +
  scale_colour_manual(values=Chrysa_colors)
ggsave(here(PlotFolderIs, "Simulated_blood.png"), dpi = 300)


# To plot glucuronide
PlotConc %>%
  filter(Tissue %in% c("Intestine", "Liver", "Liver_glu")) %>% #"Intestine_glu", )) %>%
  ggplot(aes(Time, Concentration, color = Tissue)) +
  geom_path(linewidth = 0.5) +
  labs(title = "Zearalenone and zearalenone-glucuronide intestine",
       x = "Time (h)",
       y = "Concentration (ng/ml)",
       caption = paste(AdminType)) +
  CP_theme() +
  scale_colour_manual(values=Chrysa_colors)
ggsave(here(PlotFolderIs, "Simulated_GlucTissueConc.png"), dpi = 300)



#---End-----------------------------------------------------------------------#
#*****************************************************************************#

