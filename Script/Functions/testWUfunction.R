# Testing functions as writen in DOI 10.1002/jps
# Coding what they call: Whole-Body PBPK Model A

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
dAI/dt = CL_in1*fu*CIB + ka*AIL    - (CL_ef1 + CL_met + CLsec1)*fu*CI
#        fromblood       fromlumen   totissue tomet    tolumen

#Intestinal lumen
dAIL/dt = CLsec1*fu*CI - (ka     + kg1)*AIL + Qbile*Cbile + kh*CIL.met
#         fromtissue     tolumen togluc?!     frombile      fromhydrolysis?

#Liver blood
dALB/dt = QI*fu*CIB + QHA*fu*CR    + CL_el1*CL - (QI + QHA + CL_in1)*fu*CLB
#         fromIntBl   fromHepaticArt fromLtissue toBlood     totissue?

#Liver tissue
dAL/dt = CL_in1*fu*CLB - (CL_ef1 + CL_met + CLbile)*CL
#        fromblood       to?!      tomet    tobile

#Bile
dABi/dt = CLbile*CL - Qbile*CBi