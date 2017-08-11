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
T[0] = 278.;# K
T[3] = 972.;# K
P[0] = 1000.; # kPa
P[1] = r_p*P[0]
P[3] = P[0]
P[2] = P[1]

h[0] = workFluid.enthalpy(T[0])
s[0] = workFluid.entropy(T[0],P[0])
T_S[1] = workFluid.isentropicWork(T[0],P[0],P[1])
h_s[1] = workFluid.enthalpy(T_S[1])

print h
