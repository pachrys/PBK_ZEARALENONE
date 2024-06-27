# Testing functions as writen in DOI 10.1002/jps
# Coding what they call: Whole-Body PBPK Model A

#Blood (R)
dAR/dt = QHP*CHP/KHP + QPP*CPP/KPP + 
  (QHA + QI)*fu*CLB +
  - (CL + QHP + QPP + QHA + QI)*Fu*CR;
