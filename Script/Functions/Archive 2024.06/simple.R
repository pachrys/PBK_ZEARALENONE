# Testing simple absorption; with gut lumen and gut tissue; what is absorbed goes to the liver

MW <- 318.37 #g/mol
BW <-0.26 #L

#Step2
amount.units = "umol"
time.units = "h"
nbr.doses = 1
sampling.time.from = 0
sampling.time.to = 24
sampling.frequency = 1/100

PO_DOSE.mgkg = 1;                                             #(mg/kg) dose; to change depending on the study
PO_DOSE = PO_DOSE.mgkg*BW/MW*1000;                            #(umol) unit conversion from mg to umol

inits <- c(
  ALumen = 0,#umol
  AAbsorbed = 0,#umol
  ALumen.frombile = 0, #umol
  AGu = 0,#umol
  ALi = 0,#umol
  Ametabolized = 0,#umol
  MassBal = 0 #umol
)

params <- c(
  ka = 0.8216399,#/h
  Vmax = 170829, #umol/h
  Km = 6750, #uM
  kBi = 2.667, #/h
  Fup = 0.069,
  VGu = 0.68, #L
  VLi = 2.01, #L
  PGu = 0.931,
  PLi = 1.16, 
  QGu = 0.4848, #L/min
  QLi = 1.1616, #L/min
  CO = 0.4848 + 1.1616
)


ZEN.PBKmodel <- rxode2({
  
  CGu = AGu / VGu; #umol/L
  CLi = ALi / VLi;
  
  CvGu = CGu / PGu;
  CvLi = ALi / PLi;
  
  Cv = CvLi * (QLi + QGu) / CO;
  Cart = Cv;
  
  Metabolism = (Vmax * (CLi * Fup)) / (Km + (CLi * Fup));
  
  d/dt(ALumen) = -ka*ALumen;
  d/dt(ALumen.frombile) = + kBi*Ametabolized - ka*ALumen.frombile;
  
  d/dt(AAbsorbed) = ka*ALumen + ka*ALumen.frombile;
  d/dt(AGu) = QGu*(Cart - CvGu) + ka*ALumen + ka*ALumen.frombile;
  
  d/dt(ALi) = QLi*(Cart - CvLi) + QGu*(CvGu - CvLi) - Metabolism;
  d/dt(Ametabolized) = Metabolism - kBi*Ametabolized;
  d/dt(Ablood) = QLi *(CvLi) + QGu*(CvLi) - QGu*(Cart) - QLi*(Cart)
  
  d/dt(MassBal) = AAbsorbed - (Ametabolized + AGu + ALi);
  
})


events <- eventTable(amount.units = amount.units,
                     time.units = time.units) %>%
  add.dosing(dose = PO_DOSE, dur = 0.01, nbr.doses = nbr.doses,
             dosing.to = "ALumen") %>%                 
  add.sampling(seq(from = sampling.time.from, 
                   to = sampling.time.to, 
                   by = sampling.frequency))

Output.Model <- rxSolve(ZEN.PBKmodel , 
                        params = params, 
                        events = events,
                        inits = inits) %>% 
  as.data.frame()


# Mass Balance and Error plot

Total <- PO_DOSE
Calculated <- rowSums(Output.Model[, c("Ametabolized", "ALi", "AGu", "ALumen", "ALumen.frombile")])
Massbalance <- Total - Calculated + 1
ErrorAll <- (Total - Calculated) / (Total + 1*10^(-30)) * 100                    #Mass error (% of mass lost)

MassBalance_Error <- cbind("time" = Output.Model[, "time"], 
                           Massbalance, ErrorAll) %>%
  as.data.frame()

# Mass Balance and Error plot
source(here("Script/Functions/PlottingFunction.R"))

MassBalance_Error %>%
  filter(time >= 0.01) %>% #as PO administration time is 0,01h
  ggplot(aes(time, Massbalance)) +
  geom_line(color = "grey") +
  geom_line(aes(time, ErrorAll), color = "grey6") +
  labs(x = "Time (h)", y = "Value",
       caption = "Massbalance in grey, Error in black") +
  CP_theme()


Output.Model %>% 
  ggplot()+ 
  geom_path(aes(time, CvGu), color = "purple") +
  geom_path(aes(time, CvLi), color = "blue") +
  geom_path(aes(time, Cv), color = "red") +
  labs(x = "Time (h)",
       y = "Concentration (ng/ml)", 
       caption = "Red coulor = [venous], Blue colour = [liver], Purpule = [gut]") +
  CP_theme()

  


