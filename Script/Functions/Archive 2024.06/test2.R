# pi_solution <- 1 - (Adissolved/Atotal) #Psachoulias
# psi_suspension = 1 - (Adissolved/Atotal)

MW <- 318.37 #g/mol
BW <-0.26 #L

GIT <- data.frame( ###rat!!
  GIT_Sections = c("St", "In1", "In2", "In3", "In4", "In5", "In6", "In7", "Co"),
  Section_volume = c(3, 0.6, 0.66, 0.66, 0.41, 0.41, 0.41, 0.41, 3),#ml, as per grandoni {Grandoni, 2019 #457}
  Section_pH = c(3, 7.1, 7.3, 7.5, 7.7, 7.9, 8, 7.4, 7.6), #as per grandoni {Grandoni, 2019 #457}
  Section_solubility = c(0.30, 0.31, 0.32, 0.33, 0.35, 0.37, 0.39, 0.32, 0.34) #mg/mlsolubility as per Chemicalize
) %>% mutate(
  Section_volume = Section_volume/1000*BW/0.25, #converting ml to L and then adjusting volumes as in grandoni BW of rat is 0.25kg instead of 0.26
  Section_solubility = Section_solubility*1000/MW #umol/L
)

ZEN.dat <- read.csv(here("InputData/ZearalenoneSpecificInput.csv"), header=1)          #from {Agahi, 2020 #456} and chemicalize, BP from Shin, Hong et al.2009

# Ionization status depending on pH
Section_ionisation <- 10^(GIT$Section_pH - ZEN.dat$pKa)
GIT <- cbind(GIT, Section_ionisation)

Vnames <- paste("V", GIT$GIT_Sections, sep="")
Sinames <- paste("Si", GIT$GIT_Sections, sep="")
Ionnames <- paste("Ion", GIT$GIT_Sections, sep="")
Anames = paste("A", GIT$GIT_Section, sep = "") 


#Step1
Parameters <- data.frame(
  names = c(Vnames, Sinames, Ionnames), 
  parameters = c(GIT$Section_volume, GIT$Section_solubility, GIT$Section_ionisation)) %>% 
  pivot_wider(names_from = names, values_from = parameters)

Initials <- data.frame(
  names = Anames,
  Amount = c(0, 0, 0, 0, 0, 0, 0, 0, 0) #in umol
) %>% pivot_wider(names_from = names, values_from = Amount)


#Step2
amount.units = "umol"
time.units = "h"
nbr.doses = 1
sampling.time.from = 0
sampling.time.to = 24
sampling.frequency = 1/200

PO_DOSE.mgkg = 1;                                             #(mg/kg) dose; to change depending on the study
PO_DOSE = PO_DOSE.mgkg*BW/MW*1000;                            #(umol) unit conversion from mg to umol

Inits <- c(
  Initials, 
  AStundis = 0, 
  AIn1undis = 0, 
  AIn2undis = 0, 
  AIn3undis = 0, 
  AIn4undis = 0, 
  AIn5undis = 0, 
  AIn6undis = 0, 
  AIn7undis = 0, 
  ACoundis = 0,
  AAbs = 0, 
  AEl = 0
)

Params1 <- c(
  kTco = 0.2083333, 
  SiCo = Parameters$SiCo, 
  SiIn7 = Parameters$SiIn7, 
  SiIn6 = Parameters$SiIn6, 
  SiIn5 = Parameters$SiIn5, 
  SiIn4 = Parameters$SiIn4, 
  SiIn3 = Parameters$SiIn3, 
  SiIn2 = Parameters$SiIn2, 
  SiIn1 = Parameters$SiIn1, 
  kTsi = 4.761905, 
  ka = 0.8216399, 
  kTge = 3.846154, 
  SiSt = Parameters$SiSt, 
  VCo = Parameters$VCo, 
  VIn7 = Parameters$VIn1, 
  VIn6 = Parameters$VIn6, 
  VIn5 = Parameters$VIn5, 
  VIn4 = Parameters$VIn4, 
  VIn3 = Parameters$VIn3, 
  VIn2 = Parameters$VIn2, 
  VIn1 = Parameters$VIn1, 
  VSt = Parameters$VSt
)

#Step3

GITModel <- rxode({
  CSt = ASt/VSt;
  CIn1 = AIn1/VIn1;
  CIn2 = AIn2/VIn2;
  CIn3 = AIn3/VIn3;
  CIn4 = AIn4/VIn4;
  CIn5 = AIn5/VIn5;
  CIn6 = AIn6/VIn6;
  CIn7 = AIn7/VIn7;
  CCo = ACo/VCo;
  
  
  d/dt(AStundis) = - AStundis*SiSt - AStundis*kTge; #what gets dissolved
  d/dt(ASt) = AStundis*SiSt - ASt*kTge - ASt*ka;
  
  d/dt(AIn1undis) = AStundis*kTsi - AIn1undis*SiIn1 - AIn1undis*kTsi; #what gets dissolved
  d/dt(AIn1) = ASt*kTsi + AIn1undis*SiIn1 - AIn1*kTsi - AIn1*ka;
  
  d/dt(AIn2undis) = AIn1undis*kTsi - AIn2undis*SiIn2 - AIn2undis*kTsi; #what gets dissolved
  d/dt(AIn2) = AIn1*kTsi + AIn2undis*SiIn2 - AIn2*kTsi - AIn2*ka;
  
  d/dt(AIn3undis) = AIn2undis*kTsi - AIn3undis*SiIn3 - AIn3undis*kTsi; #what gets dissolved
  d/dt(AIn3) = AIn2*kTsi + AIn3undis*SiIn3 - AIn3*kTsi - AIn3*ka;
  
  d/dt(AIn4undis) = AIn3undis*kTsi - AIn4undis*SiIn4 - AIn4undis*kTsi; #what gets dissolved
  d/dt(AIn4) = AIn3*kTsi + AIn4undis*SiIn4 - AIn4*kTsi - AIn4*ka;
  
  d/dt(AIn5undis) = AIn4undis*kTsi - AIn5undis*SiIn5 - AIn5undis*kTsi; #what gets dissolved
  d/dt(AIn5) = AIn4*kTsi + AIn5undis*SiIn5 - AIn5*kTsi - AIn5*ka;
  
  d/dt(AIn6undis) = AIn5undis*kTsi - AIn6undis*SiIn6 - AIn6undis*kTsi; #what gets dissolved
  d/dt(AIn6) = AIn5*kTsi + AIn6undis*SiIn6 - AIn6*kTsi - AIn6*ka;
  
  d/dt(AIn7undis) = AIn6undis*kTsi - AIn7undis*SiIn7 - AIn7undis*kTsi; #what gets dissolved
  d/dt(AIn7) = AIn6*kTsi + AIn7undis*SiIn7 - AIn7*kTsi - AIn7*ka;
  
  d/dt(ACoundis) = AIn7undis*kTsi - ACoundis*SiCo - ACoundis*kTco; #what gets dissolved
  d/dt(ACo) = AIn7*kTco + ACoundis*SiCo - ACo*kTco - ACo*ka;
  
  d/dt(AAbs) = ka * (ASt+ AIn1+ AIn2+ AIn3+ AIn4+ AIn5+ AIn6+ AIn7+ ACo)
  d/dt(AEl) = ACo*kTco + ACoundis*kTco
})


events <- eventTable(amount.units = amount.units,
                               time.units = time.units) %>%
  add.dosing(dose = PO_DOSE, dur = 0.01, nbr.doses = nbr.doses,
             dosing.to = "ASt") %>%                 
  add.sampling(seq(from = sampling.time.from, 
                   to = sampling.time.to, 
                   by = sampling.frequency))

Output.Model <- rxSolve(GITModel, 
                        params = Params1, 
                        events = events,
                        inits = Inits) %>% 
  as.data.frame()


#Step4

Output.Model %>% 
  ggplot() + 
  geom_path(aes(time, ASt), color = "purple") +
  geom_path(aes(time, AIn1), color = "slateblue1") +
  geom_path(aes(time, AIn2), color = "slateblue3") +
  geom_path(aes(time, AIn3), color = "slateblue4") +
  geom_path(aes(time, AIn4), color = "steelblue4") +
  geom_path(aes(time, AIn5), color = "skyblue4") +
  geom_path(aes(time, AIn6), color = "turquoise4") +
  geom_path(aes(time,AIn7), color = "turquoise") +
  geom_path(aes(time, ACo), color = "green") +
  geom_path(aes(time, AAbs), color = "red") +
  geom_path(aes(time, AEl), color = "orange")


