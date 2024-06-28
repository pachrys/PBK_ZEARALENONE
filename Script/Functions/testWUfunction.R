# Testing functions as writen in DOI 10.1002/jps
# Coding what they call: Whole-Body PBPK Model A
# ka = absorption rate constant
# kg = GIT transit/degradation rate constant
# kh = GUS hydrolysis rate constant
# CL_met = intestinal metabolism to glucuronide
# CL_sec = excretion intestinal tissue to lumen
# CL_in = excretion intestinal blood to tissue/ or liver blood to tissue
# CL_ef = absorption intestinal tissue to blood/ or liver tissue to blood
# Qbile
# CL_bile = excretion to bile
# CLB = Concentration in Liver Blood
# check ref 22, 25 and 17, 20

#Volumes
VR = 6.6 #ml reservoir?! 2.64% of total body weight (250g)
VIL = 2 #ml intestinal lumen
VIB = 1.62 #ml intestinal blood
VI = 3 #ml intestinal tissue
VL = 9.07 #ml liver tissue
VLB = 1.81 #ml liver blood
VHP = 20 #ml 8.8% from total bw (value for 250g) highly perfused organs
VPP = 175 #ml 69.8% of total bw (value for 250g)


#Flows
QI = 8 #ml/min intestinal blood flow
QHA = 6 #ml/min hepatic artery blood flow
Qbile = 6.3 #ml/min bile flow rate
QHP = 28 #ml/min 13.8% of cardiac output
QPP = 8 #ml/min 48.5% of cardiac output

#Blood (R)
dAR/dt = QHP*CHP/KHP + QPP*CPP/KPP + 
  (QHA + QI)*fu*CLB +
  - (CL + QHP + QPP + QHA + QI)*fu*CR

#Higly perfused organs (HP)
dAHP/dt = QHP*(fu*CR - CHP/KHP)

#Poorly perfused organs (PP)
dAPP/dt = QPP*(fu*CR - CPP/KPP)

#Intestinal blood
dAIB/dt = QI*CR   + CL_ef1*CI            - (QI+CL_in1)*Fu*CIB
#         fromblood fromintestinaltissue   toblood totissue

#Intestinal tissue
dAI/dt = CL_in1*fu*CIB + ka*AIL    - (CL_ef1 + CL_met + CL_sec1)*fu*CI
#        fromblood       fromlumen   totissue  tomet    tolumen

#Intestinal lumen
dAIL/dt = CL_sec1*fu*CI - (ka     + kg1)*AIL + Qbile*Cbile + kh*AIL.gluc
#         fromtissue     tolumen togluc?!      frombile       ?

#Liver blood
dALB/dt = QI*fu*CIB + QHA*fu*CR    + CL_el1*CL - (QI + QHA + CL_in1)*fu*CLB
#         fromIntBl   fromHepaticArt fromLtissue toBlood     totissue?

#Liver tissue
dAL/dt = CL_in1*fu*CLB - (CL_ef1 + CL_met + CL_bile1)*CL
#        fromblood       to?!      tomet    tobile

#Bile
dABi/dt = CL_bile1*CL - Qbile*CBi

## Metabolite
## #Blood (R)
dAR.gluc/dt = QHP*CHP.gluc/KHP.gluc + QPP*CPP.gluc/KPP.gluc + 
  (QHA + QI)*fu.gluc*CLB.gluc +
  - (CL + QHP + QPP + QHA + QI)*fu.gluc*CR.gluc

#Higly perfused organs (HP)
dAHP.gluc/dt = QHP*(fu.gluc*CR.gluc - CHP.gluc/KHP.gluc)

#Poorly perfu.glucsed organs (PP)
dAPP.gluc/dt = QPP*(fu.gluc*CR.gluc - CPP.gluc/KPP.gluc)

#Intestinal blood
dAIB.gluc/dt = QI*CR.gluc   + CL_ef2*CI.gluc            - (QI    +CL_in2)*fu.gluc*CIB.gluc
#              fromblood      fromintestinaltissue        toblood totissue

#Intestinal tissue
dAI.gluc/dt = CL_in2*fu.gluc*CIB.gluc + ka2*AIL.gluc  + CL_met*fu*CI  - (CL_ef2 + CL_sec2)*fu.gluc*CI.gluc
#             fromblood                 fromlumen       fromparent      totissue  tolumen

#Intestinal lumen
dAIL.gluc/dt = CL_sec2*fu.gluc*CI.gluc - (ka2     + kg2 +kh)*AIL.gluc + Qbile*Cbile.gluc
#              fromtissue               tolumen    to?  toparent       frombile      

#Liver blood
dALB.gluc/dt = QI*fu.gluc*CIB.gluc + QHA*fu.gluc*CR.gluc + CL_ef2*CL.gluc - (QI + QHA + CL_in2)*fu.gluc*CLB.gluc
#              fromIntBl             fromHepaticArt           fromLtissue   toBlood     totissue?

#Liver tissue
dAL.gluc/dt = CL_in2*fu.gluc*CLB.gluc + CL_met*CL - (CL_ef2 + CL_bile2)*CL.gluc
#             fromblood                  to?!      tomet      tobile

#Bile
dABi.gluc/dt = CL_bile2*CL.gluc - Qbile*CBi.gluc