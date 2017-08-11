import pyGAS as pg
import numpy as np


workFluid = pg.He()

# declare state point arrays
h = np.zeros((4,),dtype=np.float64)
T = np.zeros_like(h)
T_S = np.zeros_like(h)
P = np.zeros_like(h)
h_s = np.zeros_like(h)
s = np.zeros_like(h)
s_s = np.zeros_like(h)

r_p = 4.;
eta_comp = 1.0;
eta_turb = 1.0;

T[0] = 278.;# K
T[2] = 972.;# K
P[0] = 1000.; # kPa
P[1] = r_p*P[0]
P[3] = P[0]
P[2] = P[1]

h[0] = workFluid.enthalpy(T[0])
s[0] = workFluid.entropy(T[0],P[0])
T_S[1] = workFluid.isentropicWork(T[0],P[0],P[1])
h_s[1] = workFluid.enthalpy(T_S[1])

h[1] = h[0] - (h[0]-h_s[1])/eta_comp;
w_comp = h[0] - h[1]

# isobaric heat addition
h[2] = workFluid.enthalpy(T[2])
q_s = h[2] - h[1]
s[2] = workFluid.entropy(T[2],P[2])

# expansion in the turbine
T_S[3] = workFluid.isentropicWork(T[2],P[2],P[3])
h_s[3] = workFluid.enthalpy(T_S[3])
h[3] = h[2] - eta_turb*(h[2] - h_s[3])

w_turb = h[2]-h[3]

# isobaric heat out
q_r = h[0] - h[3]

w_net = w_turb + w_comp
q_net = q_s + q_r


print "Net Work = %g; Net Heat = %g\n"%(w_net,q_net)
print "h = ", h

eta_th = w_net/q_s
print "Thermal efficiency = %g \n"%eta_th

