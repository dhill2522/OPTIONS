"""
Created on Mon Nov 05 03:52:36 2018
@author: Paul
"""

### Boiler-Plate ###
from threading import Thread
import matplotlib.pylab as plt
import numpy as np
import scipy as sp
from numpy import random
import time

from Func import *
from iapws97 import _PSat_T


###############################################################################
"""""""""   All PERCS Equipment Classes   """"""""" ###########################
###############################################################################

class Pipe: # Pipes 716 & 717
    def __init__(self,len_):
        self.Di = 0.5 # ft
        self.len = len_ # ft
        self.cost = 0.0 # $
    def calc_Pipe(self):
        self.cost = 50.0*self.len # $

class Support: # PERCS Tank Support Structure
    def __init__(self,R_tank,H_tank,elev_tank):
        self.R = R_tank # ft
        self.H = H_tank # ft
        self.elev = elev_tank # ft
        self.cost = 0.0 # $
    def calc_Support(self):
        profile = np.pi*self.R**2.0 # ft^2
        cost_per_sqft = 30 # $/ft^2
        elev_factor = 0.0628*self.elev + 1.0
        self.cost = elev_factor * profile * cost_per_sqft

class HX: # Fake Heat Exchanger
    def __init__(self):
        self.P = 155.132 # bar (hot leg pressure)
        self.n = 0       # (number of tubes)
        self.A = 0.0     # m^2 (heat transfer SA)
        self.cost = 0.0  # $
    def calc_HX(self):
        # Calculate Cost
        K = np.array([4.1884,-0.2503,0.1974])
        if self.A <= 1000:
            # Purchase Cost
            C_p0 = 10.0**(K[0]+K[1]*np.log10(self.A)+K[2]*np.log10(self.A)**2.0)
        else:
            C_p0 = 11665.8957777+152.0393955*self.A
        P_g = self.P-1.0 # barg
        if P_g < 5.0: 
            C = np.array([0.0,0.0,0.0])
        else:
            C = np.array([0.03881,-0.11272,0.08183])
        # Pressure Factor
        F_P = 10.0**(C[0]+C[1]*np.log10(P_g)+C[2]*np.log10(P_g)**2.0)
        F_M = 2.75  # Material Factor (Stainless Steel)
        B_1 = 1.63
        B_2 = 1.66
        self.cost = (553.9/397.0)*C_p0*(B_1 + B_2*F_M*F_P)

class Tank: # PERCS Tank
    def __init__(self,R_,H_):
        self.P = 1.703 # bar
        self.R = R_ # ft
        self.D = self.R*2/3.28084 # m
        self.H = H_ # ft
        self.A = np.pi * (self.R**2.0) * self.H # ft^3
        self.A = self.A / 35.3147 # m^3
        self.cost = 0.0
    def calc_Tank(self):
        # Calcuate Cost
        K = np.array([3.4974,0.4485,0.1074])
        if self.A <= 520:
            C_p0 = 10.0**(K[0]+K[1]*np.log10(self.A)+K[2]*np.log10(self.A)**2.0) # Purchase Cost
        else:
            C_p0 = 637.3687*self.A-9491.036783
        P_g = self.P-1.0 # barg
        F_P = ((P_g+1.0)*self.D/(2.0*(850.0-0.6*(P_g+1.0)))+0.00315)/0.0063 # Pressure Factor
        F_M = 3.12  # Material Factor (SS)
        B_1 = 2.25
        B_2 = 1.82
        self.cost = (553.9/397.0)*C_p0*(B_1 + B_2*F_M*F_P)

class Chemical: # MgCO3 within PERCS Tank
    def __init__(self,ID):
        self.mass = 0.0 # kg
        self.ID = ID
        self.chem_costs = np.array((24.0,0.0)) # $/kg
        self.cost = 0.0
    def calc_Chemical(self):
        self.cost = self.mass * self.chem_costs[self.ID]


###############################################################################
"""""""""   All PCS Equipment Classes   """"""""" #############################
###############################################################################

class Stream:
    def __init__(self,P,T,mdot,x):
        self.y = 1
        self.P = P
        self.T = T
        self.mdot = mdot
        self.x = x

class Turbine:
    def __init__(self,Pin,Tin,mdot,x_in,Pout):
        self.y = 1
        self.Pin = Pin
        self.Tin = Tin
        self.mdot = mdot
        self.x_in = x_in
        self.Pout = Pout
        self.Tout = 0.0
        self.x_out = 0.0
        self.W = 0.0
        self.eff = 0.915
        self.cost = 0.0
    def calc_Turb(self):
        # If turbine does exist
        if self.y == 1:
            # If the incoming steam is saturated or superheated
            if self.x_in >= 1.0:
                Hin = h_pT(self.Pin,self.Tin)
            # If the incoming stream is two-phase
            else:
                Hin = h_Tx(self.Tin,self.x_in)
            # Calculate the Power Generated (W)
            Sin = S_ph(self.Pin,Hin)
            S_out_id = Sin
            H_out_id = h_pS(self.Pout,S_out_id)
            DH_id = Hin - H_out_id
            DH_real = self.eff*DH_id
            self.W = self.mdot*DH_real/1000 # MW
            # Calculate the outlet properties
            H_out = Hin - DH_real
            self.Tout = T_ph(self.Pout,H_out)
            self.x_out = x_ph(self.Pout,H_out)
            # Calcuate Cost
            A = self.W * 1000 # kW
            K = np.array([2.7051,1.4398,-0.1776])
            if A <= 9800:
                # Purchase Cost
                C_p0 = 10.0**(K[0]+K[1]*np.log10(A)+K[2]*np.log10(A)**2.0)
            else:
                C_p0 = 410763.708588+0.87078286*A
            F_P = 0.0  # Pressure Factor
            F_M = 0.0  # Material Factor (Stainless Steel)
            B_1 = 1.0
            B_2 = 0.0
            self.cost = (553.9/397.0)*C_p0*(B_1 + B_2*F_M*F_P)
        # If turbine does not exist
        else:
            self.cost = 0.0
            self.W = 0.0

class Pump:
    def __init__(self,Pin,Tin,mdot,Pout):
        self.y = 1
        self.Pin = Pin
        self.Tin = Tin
        self.mdot = mdot
        self.Pout = Pout
        self.Tout = 0.0
        self.W = 0.0
        self.eff = 0.85
        self.cost = 0.0
    def calc_Pump(self):
        # If pump does exist
        if self.y == 1:
            if self.Pin > _PSat_T(self.Tin+273.15)*10.0:
                Hin = h_pT(self.Pin,self.Tin)
            else:
                Hin = h_Tx(self.Tin,0.0)
            # Calculate the work requirement (W)
            Sin = S_ph(self.Pin,Hin)
            S_out_id = Sin
            H_out_id = h_pS(self.Pout,S_out_id)
            DH_id = H_out_id - Hin
            DH_real = DH_id / self.eff
            self.W = self.mdot*DH_real/1000 # MW
            # Calculate the outlet properties
            H_out = Hin + DH_real
            self.Tout = T_ph(self.Pout,H_out)
            # Calcuate Cost
            A = self.W * 1000 # kW
            K = np.array([3.3892,0.0536,0.1538])
            if A <= 300:
                # Purchase Cost
                C_p0 = 10.0**(K[0]+K[1]*np.log10(A)+K[2]*np.log10(A)**2.0)
            else:
                C_p0 = 5371.29236+79.50315*A
            P_g = self.Pout - 1.0 # barg
            C = np.zeros(3)
            if P_g >= 10.0:
                C = np.array([-0.3935,0.3957,-0.00226])
            # Pressure Factor
            F_P = 10.0**(C[0]+C[1]*np.log10(P_g)+C[2]*np.log10(P_g)**2.0)
            F_M = 2.28  # Material Factor (Stainless Steel)
            B_1 = 1.89
            B_2 = 1.35
            self.cost = (553.9/397.0)*C_p0*(B_1 + B_2*F_M*F_P)
        # If pump does not exist
        else:
            self.cost = 0.0
            self.W = 0.0

class PHX: # Primary Heat Exchanger (Steam Generator)
    def __init__(self,Tout):
        self.y = 1
        # Hot stream vars (Default)
        self.P_hot = 160.71 # bar
        self.Tin_hot = 328.56 # deg C
        self.Tout_hot = 296.51 # deg C
        self.Q_th = 750.0e3 # Thermal power transferred, in kW
        # Cold stream vars
        self.Pin = 64. # bar
        self.Tin = 0.0
        self.Pout = self.Pin # Zero pressure drop
        self.Tout = Tout # Optimized param
        self.xout = 0.0
        # Overall vars
        self.mdot = 0.0
        self.cost = 0.0
    def calc_PHX(self):
        # Calculate Tcold_in
        Hout = h_pT(self.Pout,self.Tout)
        Hin = Hout - self.Q_th/self.mdot
        self.Tin = T_ph(self.Pin,Hin)
        # Calculate xout (quality leaving phx)
        self.xout = x_ph(self.Pout,Hout)
        # Find required heat-transfer Area, A
        # Log mean temperature difference
        DT_1 = self.Tin_hot-self.Tout
        DT_2 = self.Tout_hot-self.Tin
        self.DT_lm = (DT_2-DT_1)/np.log(DT_2/DT_1)
        # HX calculations
        self.F = 1. # for phase change
        U = 5 #kW/m^2*K, back-calculated from B&W params
        self.A = self.Q_th / (self.F*U*self.DT_lm) # m^2
        # Calculate Cost
        K = np.array([4.1884,-0.2503,0.1974])
        if self.A <= 1000:
            # Purchase Cost
            C_p0 = 10.0**(K[0]+K[1]*np.log10(self.A)+K[2]*np.log10(self.A)**2.0)
        else:
            C_p0 = 11665.8957777+152.0393955*self.A
        P_g = self.P_hot-1.0 # barg
        if P_g < 5.0: 
            C = np.array([0.0,0.0,0.0])
        else:
            C = np.array([0.03881,-0.11272,0.08183])
        # Pressure Factor
        F_P = 10.0**(C[0]+C[1]*np.log10(P_g)+C[2]*np.log10(P_g)**2.0)
        F_M = 2.75  # Material Factor (Stainless Steel)
        B_1 = 1.63
        B_2 = 1.66
        self.cost = (553.9/397.0)*C_p0*(B_1 + B_2*F_M*F_P)

class Reheater:
    """
    1 = Reheat Stream (Hot side)
    2 = Residual Steam (Cold side)
    """
    def __init__(self,ID,Pin1,Tin1,mdot1,x_in1,Pin2,Tin2,mdot2,Satd_in2):
        self.y = 1
        # Hot stream vars
        self.ID = ID
        self.Pin1 = Pin1
        self.Tin1 = Tin1
        self.mdot1 = mdot1
        self.x_in1 = x_in1
        self.Satd_in1 = False
        self.Tout1 = 0.0
        self.Pout1 = 0.0
        self.x_out1 = 0.0
        # Cold stream vars
        self.Pin2 = Pin2
        self.Tin2 = Tin2
        self.mdot2 = mdot2
        self.Satd_in2 = Satd_in2
        self.s_lim = 0 # Fake ID
        self.Tout2 = 0.0
        self.Pout2 = 0.0
        # Overall vars
        self.DT_lm = 0.0
        self.F = 0.0
        self.A = 0.0
        self.q = 0.0 # kW
        self.pinch = False
        self.cost = 0.0
    def calc_RH(self):
        self.pinch = False
        # If reheater does exist
        if self.y == 1:
            # Check for a pinch point
            if (self.Tin1 - self.Tin2) <= 10.0:
                self.pinch = True
                self.A = 0.0
                self.q = 0.0
                self.Pout1 = self.Pin1
                self.Tout1 = self.Tin1
                self.x_out1 = self.x_in1
                self.Pout2 = self.Pin2
                self.Tout2 = self.Tin2
            # Proceed with calcs if no pinch point
            else:
                # For now, ignore pressure drops in the RHs
                self.Pout1 = self.Pin1
                self.Pout2 = self.Pin2
                # Calcs if limiting stream is 1
                To1 = self.Tin2 + 10.0
                Hin1 = 0.0 # Fake value
                if self.Satd_in1 == False and self.x_in1 == 0.0:
                    Hin1 = h_pT(self.Pin1,self.Tin1)
                elif self.Satd_in1 == True:
                    Hin1 = h_Tx(self.Tin1,self.x_in1)
                elif 0.0 < self.x_in1 < 1.0:
                    Hin1 = h_Tx(self.Tin1, self.x_in1)
                elif self.Satd_in1 == False and self.x_in1 == 1.0:
                    Hin1 = h_pT(self.Pin1,self.Tin1)
                Ho1 = h_pT(self.Pout1,To1)
                q1 = self.mdot1*(Hin1-Ho1)
                # Calcs if limiting stream is 2
                To2 = self.Tin1 - 10.0
                if self.Satd_in2 == True:
                    Hin2 = h_Tx(self.Tin2,1.0)
                else:
                    Hin2 = h_pT(self.Pin2,self.Tin2)
                Ho2 = h_pT(self.Pout2,To2)
                q2 = self.mdot2*(Ho2-Hin2)
                # Determine which stream is actually limiting
                if q1 < q2:
                    self.s_lim = 1
                else:
                    self.s_lim = 2
                # If limiting stream is 1:
                if self.s_lim == 1:
                    self.Tout1 = To1
                    self.q = q1
                    # Apply q to Turbine stream, find the new Tout2
                    DH_2 = self.q / self.mdot2
                    Hout2 = Hin2 + DH_2
                    self.Tout2 = T_ph(self.Pout2,Hout2)
                # If limiting stream is 2:
                if self.s_lim == 2:
                    self.Tout2 = To2
                    self.q = q2
                    # Apply q to Reheat stream, find the new Tout1 and x_out1
                    DH_1 = self.q / self.mdot1
                    Hout1 = Hin1 - DH_1
                    self.Tout1 = T_ph(self.Pout1,Hout1)
                    self.x_out1 = x_ph(self.Pout1,Hout1)
            # If no pinch point, Calc the Cost
            if self.pinch == False:
                # Find required heat-transfer Area, A
                DT_1 = self.Tin1-self.Tout2
                DT_2 = self.Tout1-self.Tin2
                self.DT_lm = (DT_2-DT_1)/np.log(DT_2/DT_1)
                if self.x_in1 > 0.0:
                    self.F = 1.0
                else:
                    R = (self.Tin1-self.Tout1)/(self.Tout2-self.Tin2)
                    P = (self.Tout2-self.Tin2)/(self.Tin1-self.Tin2)
                    inside = (2.0-P*(R+1.0-np.sqrt(R**2.0+1.0)))/(2.0-P*(R+1.0+np.sqrt(R**2.0+1.0)))
                    self.F = np.sqrt(R**2.0+1.0)/(R-1.0)*np.log((1.0-P)/(1.0-P*R))/np.log(inside)
                """ U = 1 kW/m^2*K until changed on 2.17.17 """
                U = 3 # kW/m^2*K
                self.A = self.q / (self.F*U*self.DT_lm) # m^2
                # Calcuate Cost
                K = np.array([4.1884,-0.2503,0.1974])
                if self.A <= 1000:
                    # Purchase Cost
                    C_p0 = 10.0**(K[0]+K[1]*np.log10(self.A)+K[2]*np.log10(self.A)**2.0)
                else:
                    C_p0 = 11665.8957777+152.0393955*self.A
                P_g = self.Pin1-1.0 # barg
                if P_g < 5.0: 
                    C = np.array([0.0,0.0,0.0])
                else:
                    C = np.array([0.03881,-0.11272,0.08183])
                # Pressure Factor
                F_P = 10.0**(C[0]+C[1]*np.log10(P_g)+C[2]*np.log10(P_g)**2.0)
                F_M = 2.75  # Material Factor (Stainless Steel)
                B_1 = 1.63
                B_2 = 1.66
                self.cost = (553.9/397.0)*C_p0*(B_1 + B_2*F_M*F_P)
            # If there was a pinch point, enact a penalty in the cost
            elif self.pinch == True:
                self.cost = 15.0e9 / 5.0
        # If reheater does not exist
        else:
            self.cost = 0.0

class MS: # Moisture Separator
    """
    V = vapor outlet
    L = liquid outlet
    """
    def __init__(self,Pin,Tin,mdot,x_in):
        self.y = 1
        self.P = Pin
        self.T = Tin
        self.mdot = mdot
        self.x_in = x_in
        self.mdot_V = 0.0
        self.mdot_L = 0.0
        self.D = 0.0 # m
        self.V = 0.0 # m^3
        self.cost = 0.0
    def calc_MS(self):
        self.mdot_V = self.mdot * self.x_in
        self.mdot_L = self.mdot * (1-self.x_in)
        # If MS does exist
        if self.y == 1:
            # Find the volume required, A
            rho_steam = rhoV_P(self.P) # kg/m^3
            rho_water = rhoL_P(self.P) # kg/m^3
            res_time = 60.0 # sec
            A = res_time*(self.mdot_V/rho_steam+self.mdot_L/rho_water) # m^3
            self.V = A
            # Calcuate Cost
            K = np.array([3.4974,0.4485,0.1074])
            if A <= 520:
                # Purchase Cost
                C_p0 = 10.0**(K[0]+K[1]*np.log10(A)+K[2]*np.log10(A)**2.0)
            else:
                C_p0 = 637.3687*A-9491.036783
            P_g = self.P-1 # barg
            self.D = A/100.0 # m
            # Pressure Factor
            F_P = ((P_g+1.0)*self.D/(2.0*(850.0-0.6*(P_g+1.0)))+0.00315)/0.0063
            F_M = 3.12  # Material Factor (Stainless Steel)
            B_1 = 2.25
            B_2 = 1.82
            self.cost = C_p0*(B_1 + B_2*F_M*F_P)
        # If MS does not exist
        else:
            self.cost = 0.0

class Condenser:
    def __init__(self,Pin,Tin,mdot,x_in):
        # Initialize vars
        self.y = 1
        self.Pin = Pin
        self.Tin = Tin
        self.mdot = mdot
        self.x_in = x_in
        self.q = 0.0
        self.Pout = 0.0
        self.Tout = 0.0
        self.x_out = 0.0
        self.A = 0.0
        self.cost = 0.0
    def calc_Condenser(self):
        # Calculate Q across HX
        Hin = h_Tx(self.Tin,self.x_in)
        self.Pout = self.Pin
        self.Tout = self.Tin
        Hout = h_Tx(self.Tout,self.x_out)
        DH = Hin - Hout
        self.q = self.mdot * DH / 1000 # MW
        # Find required heat-transfer Area, A
        DT_1 = self.Tin-30.0
        DT_2 = self.Tout-25.0
        DT_lm = (DT_2-DT_1)/np.log(DT_2/DT_1)
        R = (self.Tin-self.Tout)/(30.0-25.0)
        P = (30.0-25.0)/(self.Tin-25.0)
        inside = (2.0-P*(R+1.0-np.sqrt(R**2.0+1.0)))/(2.0-P*(R+1.0+np.sqrt(R**2.0+1.0)))
        F = np.sqrt(R**2.0+1.0)/(R-1.0)*np.log((1.0-P)/(1.0-P*R))/np.log(inside)
        U = 1 # kW/m^2*K
        self.A = self.q / (F*U*DT_lm) # m^2
        # Calcuate Cost
        K = np.array([4.8306,-0.8509,0.3187])
        if self.A <= 1000:
            # Purchase Cost
            C_p0 = 10.0**(K[0]+K[1]*np.log10(self.A)+K[2]*np.log10(self.A)**2.0)
        else:
            C_p0 = 146.21105*self.A-6225.7924
        P_g = self.Pin-1.0 # barg
        F_P = 0.0
        # Pressure Factor
        if P_g < 5.0: 
            C = np.array([0.0,0.0,0.0])
            F_P = 1.0
        else:
            C = np.array([0.03881,-0.11272,0.08183])
            F_P = 10.0**(C[0]+C[1]*np.log10(P_g)+C[2]*np.log10(P_g)**2.0)
        F_M = 2.75  # Material Factor (Stainless Steel)
        B_1 = 1.63
        B_2 = 1.66
        # If condenser does exist, which it should...
        if self.y == 1:
            self.cost = (553.9/397.0)*C_p0*(B_1 + B_2*F_M*F_P)
        else:
            self.cost = 0.0

class FWH: # Feedwater Heater
    def __init__(self,Pin1,Tin1,mdot1,x_in1,Pin2,Tin2,mdot2):
        # Initialize vars
        self.y = 1
        self.Pin1 = Pin1
        self.Tin1 = Tin1
        self.mdot1 = mdot1
        self.x_in1 = x_in1
        self.Pin2 = Pin2
        self.Tin2 = Tin2
        self.mdot2 = mdot2
        self.Pout = 0.0
        self.Tout = 0.0
        self.mdot = 0.0
        self.x_out = 0.0
        self.D = 0.0 # m
        self.V = 0.0 # m^3
        self.cost = 0.0
    def calc_FWH(self):
        # If FWH does exist
        if self.y == 1:
            # Calculate the outlet Enthalpy
            self.Pout = self.Pin2
            self.mdot = self.mdot1 + self.mdot2
            #''' Added 0.005 degC to Tin1, b/c of a once-occuring issue with IAPWS97 '''
            Hin1 = h_Tx(self.Tin1,self.x_in1)
            x1 = self.mdot1 / self.mdot
            Hin2 = h_pT(self.Pin2,self.Tin2)
            x2 = self.mdot2 / self.mdot
            Hout = x1*Hin1 + x2*Hin2
            self.Tout = T_ph(self.Pout,Hout)
            self.x_out = x_ph(self.Pout,Hout)
            # Find the volume required, A
            rho_water = 10**3.0 # kg/m^3
            res_time = 60.0 # sec
            self.V = res_time*(self.mdot/rho_water) # m^3
            A = self.V
            # Calcuate Cost
            K = np.array([3.4974,0.4485,0.1074])
            if A <= 520:
                # Purchase Cost
                C_p0 = 10.0**(K[0]+K[1]*np.log10(A)+K[2]*np.log10(A)**2.0)
            else:
                C_p0 = 637.3687*A-9491.036783
            P_g = self.Pout-1 # barg
            self.D = A/100.0 # m
            # Pressure Factor
            F_P = ((P_g+1.0)*self.D/(2.0*(850.0-0.6*(P_g+1.0)))+0.00315)/0.0063
            F_M = 3.12  # Material Factor (Stainless Steel)
            B_1 = 2.25
            B_2 = 1.82
            self.cost = (553.9/397.0)*C_p0*(B_1 + B_2*F_M*F_P)
        # If FWH does not exist
        else:
            self.cost = 0.0


###############################################################################
"""""""""   All Core Loop Equipment Classes   """"""""" #######################
###############################################################################

class Pump_rcp:
    def __init__(self,W,P_out):
        self.Pout = P_out # bar
        self.W_th = W # MW
        self.eff = 0.85
        self.cost = 0.0
    def calc_Pump(self):
        # Calcuate Cost
        A = self.W_th * 1000 # kW
        K = np.array([3.3892,0.0536,0.1538])
        if A <= 300:
            # Purchase Cost
            C_p0 = 10.0**(K[0]+K[1]*np.log10(A)+K[2]*np.log10(A)**2.0)
        else:
            C_p0 = 5371.29236+79.50315*A
        P_g = self.Pout - 1.0 # barg
        C = np.zeros(3)
        if P_g >= 10.0:
            C = np.array([-0.3935,0.3957,-0.00226])
        # Pressure Factor
        F_P = 10.0**(C[0]+C[1]*np.log10(P_g)+C[2]*np.log10(P_g)**2.0)
        F_M = 2.28  # Material Factor (Stainless Steel)
        B_1 = 1.89
        B_2 = 1.35
        self.cost = (553.9/397.0)*C_p0*(B_1 + B_2*F_M*F_P)