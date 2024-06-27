#Testing as Aichner

# Reservoir ?? what does he mean by that?

#Stomach concentration/amount is considered to be insoluble, only the amount solubilised is being absorbed
d/dt(AStIns) <- - (kTge*AStIns) - (VmaxStsol*CSt / (KmStsol + CSt))
d/dt(ASt) <- (VmaxStsol*CStIns / (KmStsol + CStIns)) - (kTge*ASt) #-ka*ASt


d/dt(ASiIns) <- (kTge*AStIns) - (VmaxSi*CSiIns / (KmSi + CSiIns)) - (kTsi*ASiIns)#solubilisation in the small intestine
d/dt(ASi) <- (kTge*ASt) + (VmaxSi*CSi / (KmSi + CSi)) +
  - (ka*CSi) + #absorption to small intestinal tissue
  - (kTsi*ASi) #transit to large intestine

d/dt(ACoIns) <- (kTsi*ASiIns) - (VmaxLi*CCoIns / (KmLi + CCoIns)) - (kTco*ACoIns)
d/dt(ACo) <- (kTsi*ASi) + (VmaxLi*CCoIns / (KmLi + CCoIns)) + 
  - (ka*Cco) + #absorption from colon to blood
  - (kTco*ACo) #fecal excretion


solSI = if(CSiUA < SOLUA){
  (VmaxUAsol*CSiUAinsol) / (KmUAsol + CSiUAinsol)
}else {0}

# 
# d/dt(Start) <- - pulseinsol - pulsesol
# d/dt(Start) <- PULSE(AODOSEUA - VSt*SOLUA, 24, 24) - PULSE(SOLUA*VSt, 24, 24)
# #                    umol       - L*umol/L

d/dt(ASt) <- - kTge*(ASt - VSt*SOLUA) - kTge*(VSt*SOLUA) #umol
#            what leaves as insoluble - what leaves as soluble
#            total amount in the stomach - the amount that is solubilised at max; - what is solubilised at max

if(CSiUA < SOLUA){
  d/dt(ASiUAnsol) <- kTge*(ASt - VSt*SOLUA) +  #IN insoluble
    - ASiUAinsol*Ksi + #OUT insoluble
    - (VmaxUAsol*CSiUAinsol / (KmUAsol + CSiUAinsol)) #solubilisation
}else {
  d/dt(ASiUAnsol) <- kTge*(ASt - VSt*SOLUA) - ASiUAinsol*Ksi
}

d/dt(ASiUAnsol) <- kTge*(ASt - VSt*SOLUA) - ASiUAinsol*Ksi - solSI
#                  IN insoluble          OUT insoluble     - solubilisation

d/dt(ASiUA) <- kTge*(VSt*SOLUA) + solSI + 
  entRC + #entRC=ABilUA
  - ASiUA*Ksi + #what leaves as solubilised to the gut
  - KaUA*CSiUA  #what is absorbed to the gut lumen

d/dt(AGUA) <- KaUA*CSiUA - venG + artG + 
  #           absorption to the enterocyte  - IN/OUT from blood 
  - VmaxGUAGluc*CVGUA/(KmGUAGluc + CVGUA) #gut tissue glucuronidation

d/dt(ALaUA) <- + ASiUA*Ksi - ALiUA*Kli - KbUA*CLaUA #uptake from large intestinal lumen to liver 
#                            excretion   absorption from large intestine to liver

d/dt(APUA) <- artP - venP #spleen/portal vein

d/dt(UAex) <- ALiUA*Kli + ASiUAinsol*Ksi #fecal excretion

d/dt(ALiUA) <- - venL + artL + venG + venP + KbUA*CLaUA +
  - (VmaxLUAGluc*CVLUA/(KmLUAGluc + CVLUA)) #glucuronidation 
d/dt(ABilUA) <- - ABilUA + (VmaxLUAGluc*CVLUA/(KmLUAGluc + CVLUA))


d/dt(ABlUA) <- + venL - artL - artP - artG +
  venBr - artBr +
  venS - artS +
  venR - artR +
  venF - artF +
  venK - artK + 
  venM - artM

d/dt(ABRUA) <- artBr - venBr
d/dt(ASUA) <- artS - venS #slowly perfused
d/dt(ARUA) <- artR - venR #rapidly perfused?
d/dt(AFUA) <- artF - venF 
d/dt(AKUA) <- artK - venK - urineUA
d/dt(AUCUAUr) <- urineUA
d/dt(AMUA) <- artM - venM #muscle


#Dosing
ASiUAinsol <- AODOSEUA - SOLUA*VSt
