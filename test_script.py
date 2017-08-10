# -*- coding: utf-8 -*-
"""
Created on Thu Aug 10 17:36:05 2017

@author: stu
"""

import pyGAS as pg

helium = pg.He();
To = 298.15; # K
T1 = 700; #K
P = 101.3; # kPa
cp1 = helium.Cp(T1)
h0 = helium.enthalpy(To)
h1 = helium.enthalpy(T1)
s0 = helium.entropy(To,P)
s1 = helium.entropy(T1,P)
print "specific heat at %g K = %g kJ/kg-K \n"%(T1,cp1)
print "specific enthalpy at %g K = %g kJ/kg \n"%(T1,(h1))
print "specific entropy at %g K = %g kJ/kg-K \n"%(T1,(s1))