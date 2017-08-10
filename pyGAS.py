# -*- coding: utf-8 -*-
"""
Created on Thu Aug 10 16:55:33 2017

@author: stu
"""
import numpy as np

class pyGAS:
    """
     base class.  Sets the interface and defines most of the
     calculations.  Derived classes will principally hold the
     data as well as customize any required calculations
    """
    def __init__(self):
        self.molecular_weight = 0.
        self.R_M = 0.
        self.bLow = np.zeros((2,),dtype=np.float64)
        self.bHigh = np.zeros((2,),dtype=np.float64)
        self.aLow = np.zeros((5,),dtype=np.float64)
        self.aHigh = np.zeros((5,),dtype=np.float64)
        self.Tmin = 300; # C 
        self.Tmax = 5000; # C
        self.R_universal = 8.314; #kJ/kmol-K
        self.To = 298.15 # K
        self.Po = 101.3 # kPa
        
    def Cp(self,T):
        """
        specific heat as a function of temperature
        kJ/kg-K
        
        """
        # to do: if T> self.Tmax throw an exception or something
        if (T > self.Tmax):
            print "Warning! Requested temperature range is high out of range for this fluid"
        
        Tvec = np.array([1., T, T**2, T**3, T**4])
        outVal = 0.
        if T>1000.:
            outVal = np.sum(np.dot(self.aHigh,Tvec))
        else:
            outVal = np.sum(np.dot(self.aLow,Tvec))
        
        return outVal*self.R_M
        
    def enthalpy(self,T):
        """
        absolute enthalpy
        """
        return self._enthalpy(T) - self._enthalpy(self.To)
    def _enthalpy(self,T):
        """
        specific enthalpy as a function of temperature
        kJ/kg
        
        """
        if (T>self.Tmax):
            print "Warning! Requested temperature is high out of range for this species!"
        
        Tvec = np.array([T,(T**2)/2.,(T**3)/3.,(T**4)/4.,(T**5)/5.])
        outVal = 0.
        
        if T>1000.:
            outVal = np.sum(np.dot(self.aHigh,Tvec)) + self.bHigh[0]
        else:
            outVal = np.sum(np.dot(self.aLow,Tvec)) + self.bLow[0]
            
        return outVal*self.R_M  
        
    def entropy(self,T,P):
        """
        absolute entropy
        """
        return self._entropy(T,P) - self._entropy(self.To,self.Po)
        
    def _entropy(self,T,P):
        """
        specific entropy as a function of temperature
        kJ/kg-K
        """
        if (T>self.Tmax):
            print "Warning! Requested temperature is high out of range for this species!"
        
        Tvec = np.array([np.log(T),T,(T**2)/2.,(T**3)/3.,(T**4)/4.])
        outVal = 0.
        
        if T>1000.:
            outVal = np.sum(np.dot(self.aHigh,Tvec)) + self.bHigh[1]
        else:
            outVal = np.sum(np.dot(self.aLow,Tvec)) + self.bLow[1]
            
        outVal-=np.log(P/self.Po)
            
        return outVal*self.R_M  
        
class He(pyGAS):
    def __init__(self):
        pyGAS.__init__(self)
        self.molecular_weight = 4.0026; #kg/kmol
        self.R_M = 2.0769; #kJ/kg-K
        self.Tmin = 200
        self.Tmax = 6000;
        self.aHigh[:] = [2.50,0,0,0,0]; #mono-atomic gas
        self.bHigh[:] = [-7.45375e2,0.92872394]
        self.aLow[:] = [2.5,0,0,0,0];
        self.bLow[:] = [2.85315086e5,1.62166556]