# -*- coding: utf-8 -*-
"""
Created on Thu Aug 10 16:55:33 2017

@author: stu
"""
import numpy as np
import scipy.optimize as opt

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
        self.To_s = 293.15 # K
        self.To_h = 298.15
        self.Po = 101.3 # kPa
        
    
    def v_TP(self,T,P):
        """
         returns specific volume in m^3/kg
         T - K
         P - kPa
         R_M - kJ/kg-K
         
        """
        return self.R_M*T/P 
        
    def T_vP(self,v,P):
        """
        returns Temperature as a function of specific
        volume and pressure
        
        """
        return P*v/self.R_M
        
    def P_vT(self,v,T):
       """
       
       """
       return self.R_M*T/v
    
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
    
    def Cv(self,T):
        """
        specific heat at constant volume.  This will be found by using
        the temperature dependent value of Cp and the material specific gas
        constant
        """
        
        return (self.Cp(T) - self.R_M)
    
    def k(self,T):
        """
        ratio of specific heats as a function of temperature
        
        """
        
        return (self.Cp(T)/self.Cv(T))
        
    def enthalpy(self,T):
        """
        absolute enthalpy
        """
        return self._enthalpy(T) - self._enthalpy(self.To_h)
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
        
    
    
    def isentropicWork(self,T0,P0,P1):
        """
        given T0 and P0, find P1
        """
        s0 = self.entropy(T0,P0)
        def wrappedFun(T):
            return (np.abs(self.entropy(T,P1) - s0))
        res = opt.minimize(wrappedFun,T0/2.)
        return res.x
        
    def entropy(self,T,P):
        """
        absolute entropy
        """
        return self._entropy(T,P) - self._entropy(self.To_s,self.Po)
        
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
        self.Tmin = 200.
        self.Tmax = 6000.;
        self.aHigh[:] = [2.50,0,0,0,0]; #mono-atomic gas
        self.bHigh[:] = [-7.45375e2,0.92872394]
        self.aLow[:] = [2.5,0,0,0,0];
        self.bLow[:] = [2.85315086e5,1.62166556]
        
class CO2(pyGAS):
    def __init__(self):
        pyGAS.__init__(self)
        self.molecular_weight = 44.0098
        self.R_M = 0.1889
        self.Tmin = 200.
        self.Tmax = 6000.
        self.aHigh[:] = [4.63659493, 2.74131991e-3,
                         -9.95828531e-7,1.60373011e-10,
                         -9.16103468e-15]
        self.bHigh[:] = [-4.90249341e4,-1.93534855]
        self.aLow[:] = [2.35677352,8.98459677e-3,-7.12356269e-6,
                        2.45919022e-9,-1.43699548e-13]
        self.bLow[:] = [-4.83719697e4,9.90105222]

class Air(pyGAS):
    """
    no data right now for b (low or high) or aHigh for air
    """
    def __init__(self):
        pyGAS.__init__(self)
        self.molecular_weight = 28.97; #kg/kmol
        self.R_M = 0.2870; #kJ/kg-K
        self.Tmin = 300;
        self.Tmax = 1000;
        self.aLow[:] = [3.653, -1.337e-3, 3.294e-6, -1.913e-9, 0.2763e-12]
        