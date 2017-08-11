# -*- coding: utf-8 -*-
"""
Created on Thu Aug 10 17:36:05 2017

@author: stu
"""

import pyGAS as pg

helium = pg.He();
To = 273.15; # K
T1 = 800; #K
P = 101.3; # kPa
cp1 = helium.Cp(T1)
h0 = helium.enthalpy(To)
h1 = helium.enthalpy(T1)
s0 = helium.entropy(To,P)
s1 = helium.entropy(T1,P)

v1 = helium.v_TP(T1,P)
T_check = helium.T_vP(v1,P)

print "specific heat at %g K = %g kJ/kg-K \n"%(T1,cp1)
print "specific enthalpy at %g K = %g kJ/kg \n"%(T1,(h1))
print "specific entropy at %g K = %g kJ/kg-K \n"%(T1,(s1))
print "specific volume at %g K, %g kPa = %g m^3/kg \n"%(T1,P,v1)
print "Temperature at %g m^3/kg, %g kPa = %g K \n"%(v1,P,T_check)

CO2 = pg.CO2();
s1 = CO2.entropy(T1,P)

P2 = 10.*P
T2 = CO2.isentropicWork(T1,P,P2)
s2 = CO2.entropy(T2,P2)
h1 = CO2.enthalpy(T1)

print "specific enthalpy at %g K = %g kJ/kg \n"%(T1,(h1))
print "specific entropy at %g K = %g kJ/kg-K \n"%(T1,(s1))

print "Cp(CO2) at 250K = %g kJ/kg-K \n"%(CO2.Cp(250))
print "h(Air) at 650K = %g kJ/kg \n"%(pg.Air().Cp(650))

print "T2 = %g K \n"%T2
print "s(CO2) at %g K and %g kPa = %g kJ/kg-K \n"%(T2,P2,s2)