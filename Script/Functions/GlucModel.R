# f = fraction unbound
# 1 = aglycone
# 2 = glucuronide
# D = concentration parent
# M = concentration metabolite
# 
# C = A/V => A = C*V
# venous = C/KP
# 

#Blood (R in their model)
dA/dt = QHP*CvHP + QPP*CvPP + (QHA + QI)*f_free*C - (CL + QHP + QPP + QHA + QI) * f_free*CR
