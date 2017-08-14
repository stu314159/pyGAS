# ML-1 analysis script using pyGAS

import pyGAS as pg
import numpy as np

# working fluid is N2
workFluid = pg.N2();

numSp = 6;

# declare variables
h = np.zeros((numSp,),dtype=np.float64)
T = np.zeros_like(h)
T_S = np.zeros_like(h)
P = np.zeros_like(h)
h_s = np.zeros_like(h)
s = np.zeros_like(h)
s_s = np.zeros_like(h)

# problem parameters
eta_comp = 0.84 # compressor efficiency
eta_turb = 0.85 # turbine efficiency
effec_recp = 0.8 # recuperator effectiveness
effec_prec = 0.91 # pre-cooler effectiveness

# state point 0: pre-cooler outlet / compressor inlet
P[0] = 846.67; #kPa
T[0] = 343.48; # K

h[0] = workFluid.enthalpy(T[0])
s[0] = workFluid.entropy(T[0],P[0])

# compression from state point 0 to 1
P[1] = 2209.1; #kPa
T_S[1] = workFluid.isentropicWork(T[0],P[0],P[1])
h_s[1] = workFluid.enthalpy(T_S[1])
h[1] = h[0] - (h[0] - h_s[1])/eta_comp
w_comp = h[0] - h[1]

# save the recuperator until the recuperator inlet on high temp
# side is also known
P[2] = P[1] # recuperator is assumed isobaric

# Rx outlet
P[3] = P[2] # reactor core is assumed isobaric (for now...)
T[3] = 922.04
h[3] = workFluid.enthalpy(T[3])
s[3] = workFluid.entropy(T[3],P[3])

# expansion of working fluid through turbine
P[4] = P[0] # all heat exchanges assumed isentropic
T_S[4] = workFluid.isentropicWork(T[3],P[3],P[4])
h_s[4] = workFluid.enthalpy(T_S[4])
h[4] = h[3] - eta_turb*(h[3] - h_s[4])

# apply recuperator effectiveness to get h[1] and h[5]
h[2] = h[1] + effec_recp*(h[4]-h[1])
h[5] = h[4] - (h[2] - h[1]) # apply energy balance for equal mass flow rates


# now that I have all of the enthalpies; calculate
# other outputs:
q_s = h[3] - h[2] # heat supplied by the reactor
q_r = h[0] - h[5] # heat rejected in the pre-cooler

w_turb = h[3] - h[4] # turbine work

w_net = w_turb + w_comp
q_net = q_s + q_r
eta_th = w_net/q_s


print "work net = %g \n"%w_net
print "net heat = %g \n"%q_net
print "Thermal efficiency = %g percent \n"%(eta_th*100)


