"""
Created on Fri Sep 28 2018
@author: Paul
"""

### Boiler-Plate ###
from threading import Thread
import matplotlib.pylab as plt
import numpy as np
import scipy as sp
from numpy import random
import time

from Class1_Eq import *
from Func import *


###############################################################################
"""""""""   PCS Superstructure Option Class   """"""""" #######################
###############################################################################


class PCS_Option:
    """
    Inputs:
        x = PCS superstructure x-optimization parameters
        y = PCS superstructure y-optimization parameters
    Parameters:
        *Individual optimization parameters (explained in __init__() function)
        s = Array of Stream instances for all 37 PCS superstructure streams
        phx = PHX instance representing the Steam Generator
        t1a, t1b, t1c, t2a, t2b = Turbines representing the diff. stages
        t1, t2 = Actual turbines (for cost purposes)
        t3, t4, t5 = Turbine instances for LPTs
        ms1, ms2 = Moisture separator instances
        rh1, rh2 = Reheater heat exchanger instances
        cond = Condenser instance
        fwh1, fwh2, fwh3, fwh4 = Feedwater heater instances
        p1, p2, p3, p4, p5, p6 = Pump instances
        eff = Obj. 1 - Thermodynamic Efficiency
        cost = Obj. 2 - Total PCS equipment cost
        fmm = Maximin fitness value
    Functions:
        eff() - Calcs & Returns the thermodynamic efficiency
        inv_eff() - Calcs & Returns the inverse of the efficiency
        cost() - Calcs & Returns the total equipment capital cost
        color() - Calcs & Returns the PCS color/superstructure configuration
        calc_Option() - Calcs correct m_dot, Calls solve_Option()
        solve_Option() - Fills out PCS superstructure & converges the cycle
    """
    def __init__(self,x_in,y_in):
        # Define the optimization parameters
        self.x = x_in
        self.To_PHX = 318.5 # degC
        self.Po_t1a = self.x[1]  # bar
        self.mf_t1a = self.x[2]
        self.Po_t1b = self.x[3]  # bar
        self.mf_t1b = self.x[4]  #
        self.Po_t1c = self.x[5]  # bar
        self.Po_t2a = self.x[6]  # bar
        self.mf_t2a = self.x[7]  #
        self.Po_t2b = self.x[8]  # bar
        self.y = y_in
        self.y_ipt = self.y[0]   # IPT
        self.y_rh1 = self.y[1]   # RH 1
        self.y_rh2 = self.y[2]   # RH 2
        self.y_s14 = self.y[3]   # s[14]
        self.y_s4  = self.y[4]   # s[4]
        self.y_s5  = self.y[5]   # s[5]
        # Initialize all Streams with zeros
        self.s = np.array([0])
        for i in range(1,37):
            self.s = np.append(self.s,Stream(0.0,0.0,0.0,0.0))
        # Create the PCS equipment w/ original opt. parameters
        self.phx = PHX(self.To_PHX)
        self.t1a = Turbine(0.0,0.0,0.0,0.0,self.Po_t1a)
        self.t1b = Turbine(0.0,0.0,0.0,0.0,self.Po_t1b)
        self.t1c = Turbine(0.0,0.0,0.0,0.0,self.Po_t1c)
        self.t1 = Turbine(0.0,0.0,0.0,0.0,self.Po_t1c)
        self.ms1 = MS(self.Po_t1c,0.0,0.0,0.0)
        self.rh1 = Reheater(1,self.Po_t1a,0.0,0.0,0.0,self.Po_t1c,0.0,0.0,False)
        self.t2a = Turbine(0.0,0.0,0.0,0.0,self.Po_t2a)
        self.t2b = Turbine(0.0,0.0,0.0,0.0,self.Po_t2b)
        self.t2 = Turbine(0.0,0.0,0.0,0.0,self.Po_t2b)
        self.ms2 = MS(self.Po_t2b,0.0,0.0,0.0)
        self.rh2 = Reheater(2,0.0,0.0,0.0,0.0,self.Po_t2b,0.0,0.0,False)
        self.t3 = Turbine(0.0,0.0,0.0,0.0,0.086)
        self.t4 = Turbine(0.0,0.0,0.0,0.0,0.086)
        self.t5 = Turbine(0.0,0.0,0.0,0.0,0.086)
        self.cond = Condenser(0.086,0.0,0.0,0.0)
        self.fwh1 = FWH(0.0,0.0,0.0,0.0,0.0,0.0,0.0)
        self.fwh2 = FWH(0.0,0.0,0.0,0.0,0.0,0.0,0.0)
        self.fwh3 = FWH(0.0,0.0,0.0,0.0,0.0,0.0,0.0)
        self.fwh4 = FWH(0.0,0.0,0.0,0.0,0.0,0.0,0.0)
        self.p1 = Pump(0.0,0.0,0.0,self.phx.Pin)
        self.p2 = Pump(0.0,0.0,0.0,self.Po_t1a)
        self.p3 = Pump(0.0,0.0,0.0,0.0)
        self.p4 = Pump(0.0,0.0,0.0,0.0)
        self.p5 = Pump(0.0,0.0,0.0,0.0)
        self.p6 = Pump(0.0,0.0,0.0,0.0)
        # Initialize fitness value
        self.fmm = 0
        # Sove the Option by performing all calcs
        self.calc_Option()
    
    # Calculate the overall efficiency of Option
    def eff(self):
        W_t1 = self.t1a.W+self.t1b.W+self.t1c.W
        W_turb = W_t1 + self.t2a.W+self.t2b.W+self.t3.W+self.t4.W+self.t5.W
        W_pump = self.p1.W+self.p2.W+self.p3.W+self.p4.W+self.p5.W
        return (W_turb - W_pump) / (self.phx.Q_th/1000)
    # Calculate the inverse overall efficiency (used for Optimization purposes)
    def inv_eff(self):
        return 1.0 - self.eff()
    # Calculate the overall cost of Option
    def cost(self):
        cost_phx = self.phx.cost
        cost_ms = self.ms1.cost + self.ms2.cost
        cost_rh = self.rh1.cost + self.rh2.cost
        cost_turb = self.t1.cost+self.t2.cost+self.t3.cost+self.t4.cost+self.t5.cost
        cost_cond = self.cond.cost
        cost_pump = self.p1.cost+self.p2.cost+self.p3.cost+self.p4.cost+self.p5.cost
        cost_fwh = self.fwh1.cost+self.fwh2.cost+self.fwh3.cost+self.fwh4.cost
        total_cost = cost_phx+cost_ms+cost_rh+cost_turb+cost_cond+cost_pump+cost_fwh
        cost_multiplier_for_nucl_safety_grade = 5.0
        return (total_cost * 10.0**-9.0) * cost_multiplier_for_nucl_safety_grade
    
    # Assign the option a color based on its y-values
    def color(self):
        c = 'black'
        z = self.y
        # If y = [0,1,0,0,0,0]
        if z[0]==0 and z[1]==1 and z[2]==0 and z[3]==0 and z[4]==0 and z[5]==0:
            c = 'red'
        # If y = [0,0,0,0,0,0]
        if z[0]==0 and z[1]==0 and z[2]==0 and z[3]==0 and z[4]==0 and z[5]==0:
            c = 'firebrick'
        # If y = [1,1,1,0,1,0]
        if z[0]==1 and z[1]==1 and z[2]==1 and z[3]==0 and z[4]==1 and z[5]==0:
            c = 'darkgreen'
        # If y = [1,1,1,0,0,1]
        if z[0]==1 and z[1]==1 and z[2]==1 and z[3]==0 and z[4]==0 and z[5]==1:
            c = 'purple'
        # If y = [1,0,1,0,0,1]
        if z[0]==1 and z[1]==0 and z[2]==1 and z[3]==0 and z[4]==0 and z[5]==1:
            c = 'deeppink'
        # If y = [1,1,1,1,0,0]
        if z[0]==1 and z[1]==1 and z[2]==1 and z[3]==1 and z[4]==0 and z[5]==0:
            c = 'blue'
        # If y = [1,0,1,1,0,0]
        if z[0]==1 and z[1]==0 and z[2]==1 and z[3]==1 and z[4]==0 and z[5]==0:
            c = 'cyan'
        # If y = [1,0,0,0,0,0]
        if z[0]==1 and z[1]==0 and z[2]==0 and z[3]==0 and z[4]==0 and z[5]==0:
            c = 'orange'
        # If y = [1,1,0,0,0,0]
        if z[0]==1 and z[1]==1 and z[2]==0 and z[3]==0 and z[4]==0 and z[5]==0:
            c = 'yellow'
        # Return color for specific y-value configuration
        return c
    
    
    def calc_Option(self):
        # Establish the test mdot
        mdot_test = 600. # kg/s
        # Assign test mdot and solve for Tin
        self.phx.mdot = mdot_test
        self.solve_Option()
        Tin = self.s[36].T
        # Calculate the real mdot
        Hin = h_pT(self.phx.Pin,Tin)
        Hout = h_pT(self.phx.Pout,self.phx.Tout)
        mdot_real = self.phx.Q_th / (Hout - Hin)
        # Assign real mdot and solve Option
        self.phx.mdot = mdot_real
        self.solve_Option()
    
    
    # Recalculate all streams and equipment from PHX-out to PHX-in
    def solve_Option(self):
        """ PHX """
        self.phx.calc_PHX()
        """ Stream 1 """
        self.s[1].P = self.phx.Pout
        self.s[1].T = self.phx.Tout
        self.s[1].mdot = self.phx.mdot
        self.s[1].x = self.phx.xout
        """ Turbine HPT_a """
        self.t1a.Pin = self.s[1].P
        self.t1a.Tin = self.s[1].T
        self.t1a.mdot = self.s[1].mdot
        self.t1a.x_in = self.s[1].x
        self.t1a.calc_Turb()
        """ Stream 2 """
        if self.y_rh1 == 1:
            self.s[2].y = 1
            self.s[2].P = self.t1a.Pout
            self.s[2].T = self.t1a.Tout
            self.s[2].mdot = self.mf_t1a * self.t1a.mdot
            self.s[2].x = self.t1a.x_out
        else:
            self.s[2].y = 0
        """ Turbine HPT_b """
        self.t1b.Pin = self.t1a.Pout
        self.t1b.Tin = self.t1a.Tout
        self.t1b.x_in = self.t1a.x_out
        if self.s[2].y == 1:
            self.t1b.mdot = (1-self.mf_t1a) * self.t1a.mdot
        else:
            self.t1b.mdot = self.t1a.mdot
        self.t1b.calc_Turb()
        """ Stream 5 """
        if self.y_s5 == 1:
            self.s[5].y = 1
            self.s[5].P = self.t1b.Pout
            self.s[5].T = self.t1b.Tout
            self.s[5].mdot = self.mf_t1b * self.t1b.mdot
            self.s[5].x = self.t1b.x_out
        else:
            self.s[5].y = 0
        """ Turbine HPT_c """
        self.t1c.Pin = self.t1b.Pout
        self.t1c.Tin = self.t1b.Tout
        self.t1c.x_in = self.t1b.x_out
        if self.s[5].y == 1:
            self.t1c.mdot = (1-self.mf_t1b) * self.t1b.mdot
        else:
            self.t1c.mdot = self.t1b.mdot
        self.t1c.calc_Turb()
        """ Turbine HPT """
        self.t1.Pin = self.t1a.Pin
        self.t1.Tin = self.t1a.Tin
        self.t1.mdot = self.t1a.mdot
        self.t1.x_in = self.t1a.x_in
        self.t1.Pout = self.t1c.Pout
        self.t1.calc_Turb()
        """ Stream 6 """
        self.s[6].P = self.t1c.Pout
        self.s[6].T = self.t1c.Tout
        self.s[6].mdot = self.t1c.mdot
        self.s[6].x = self.t1c.x_out
        """ MS 1 """
        self.ms1.P = self.s[6].P
        self.ms1.T = self.s[6].T
        self.ms1.mdot = self.s[6].mdot
        self.ms1.x_in = self.s[6].x
        self.ms1.calc_MS()
        """ Stream 7 """
        if self.y_s4==0 and self.y_s5==0:
            self.s[7].y = 1
            self.s[7].P = self.ms1.P
            self.s[7].T = self.ms1.T
            self.s[7].mdot = self.ms1.mdot_L
            self.s[7].x = 0.0
        else:
            self.s[7].y = 0
        """ Stream 8 """
        if self.y_s4==1 or self.y_s5==1:
            self.s[8].y = 1
            self.s[8].P = self.ms1.P
            self.s[8].T = self.ms1.T
            self.s[8].mdot = self.ms1.mdot_L
            self.s[8].x = 0.0
        else:
            self.s[8].y = 0
        """ Stream 9 """
        if self.y_ipt==1 and self.y_rh1==0:
            self.s[9].y = 1
            self.s[9].P = self.ms1.P
            # Add to T for the sake of h_pT(), since this stream skip RH 1
            self.s[9].T = self.ms1.T + (1e-10)
            self.s[9].mdot = self.ms1.mdot_V
            self.s[9].x = 1.0
        else:
            self.s[9].y = 0
        """ Stream 10 """
        if self.y_ipt==0 and self.y_rh1==0:
            self.s[10].y = 1
            self.s[10].P = self.ms1.P
            # Add to T for the sake of h_pT(), since this stream skip RH 1
            self.s[10].T = self.ms1.T + (1e-10)
            self.s[10].mdot = self.ms1.mdot_V
            self.s[10].x = 1.0
        else:
            self.s[10].y = 0
        """ Stream 11 """
        if self.y_rh1==1:
            self.s[11].y = 1
            self.s[11].P = self.ms1.P
            self.s[11].T = self.ms1.T
            self.s[11].mdot = self.ms1.mdot_V
            self.s[11].x = 1.0
        else:
            self.s[11].y = 0
        """ RH 1 """
        if self.y_rh1 == 1:
            self.rh1.y = 1
            self.rh1.Pin1 = self.s[2].P
            self.rh1.Tin1 = self.s[2].T
            self.rh1.mdot1 = self.s[2].mdot
            self.rh1.x_in1 = self.s[2].x
            self.rh1.Satd_in1 = False
            self.rh1.Pin2 = self.s[11].P
            self.rh1.Tin2 = self.s[11].T
            self.rh1.mdot2 = self.s[11].mdot
            self.rh1.Satd_in2 = True
        else:
            self.rh1.y = 0
        self.rh1.calc_RH()
        """ Stream 3 """
        if self.y_rh1==1 and self.y_s4==0:
            self.s[3].y = 1
            self.s[3].P = self.rh1.Pout1
            self.s[3].T = self.rh1.Tout1
            self.s[3].mdot = self.rh1.mdot1
            self.s[3].x = self.rh1.x_out1
        else:
            self.s[3].y = 0
        """ Stream 4 """
        if self.y_rh1==1 and self.y_s4==1:
            self.s[4].y = 1
            self.s[4].P = self.rh1.Pout1
            self.s[4].T = self.rh1.Tout1
            self.s[4].mdot = self.rh1.mdot1
            self.s[4].x = self.rh1.x_out1
        else:
            self.s[4].y = 0
        """ Stream 12 """
        if self.y_rh1==1 and self.y_ipt==1:
            self.s[12].y = 1
            self.s[12].P = self.rh1.Pout2
            self.s[12].T = self.rh1.Tout2
            self.s[12].mdot = self.rh1.mdot2
            self.s[12].x = 1.0
        else:
            self.s[12].y = 0
        """ Stream 13 """
        if self.y_rh1==1 and self.y_ipt==0:
            self.s[13].y = 1
            self.s[13].P = self.rh1.Pout2
            self.s[13].T = self.rh1.Tout2
            self.s[13].mdot = self.rh1.mdot2
            self.s[13].x = 1.0
        else:
            self.s[13].y = 0
        """ Turbine IPT_a """
        if self.y_ipt==1:
            self.t2a.y = 1
            id_in = 0 # Fake ID
            if self.s[9].y == 1:
                id_in = 9
            elif self.s[12].y == 1:
                id_in = 12        
            self.t2a.Pin = self.s[id_in].P
            self.t2a.Tin = self.s[id_in].T
            self.t2a.mdot = self.s[id_in].mdot
            self.t2a.x_in = self.s[id_in].x
        else:
            self.t2a.y = 0
        self.t2a.calc_Turb()
        """ Stream 14 """
        if self.y_s14==1:
            self.s[14].y = 1
            self.s[14].P = self.t2a.Pout
            self.s[14].T = self.t2a.Tout
            self.s[14].mdot = self.mf_t2a * self.t2a.mdot
            self.s[14].x = self.t2a.x_out
        else:
            self.s[14].y = 0
        """ Turbine IPT_b """
        if self.y_ipt==1:
            self.t2b.y = 1
            self.t2b.Pin = self.t2a.Pout
            self.t2b.Tin = self.t2a.Tout
            self.t2b.x_in = self.t2a.x_out
            if self.y_s14 == 1:
                self.t2b.mdot = (1-self.mf_t2a) * self.t2a.mdot
            else:
                self.t2b.mdot = self.t2a.mdot
        else:
            self.t2b.y = 0
        self.t2b.calc_Turb()
        """ Turbine IPT """
        if self.y_ipt==1:
            self.t2.y = 1
            self.t2.Pin = self.t2a.Pin
            self.t2.Tin = self.t2a.Tin
            self.t2.mdot = self.t2a.mdot
            self.t2.x_in = self.t2a.x_in
            self.t2.Pout = self.t2b.Pout
        else:
            self.t2.y = 0
        self.t2.calc_Turb()
        """ Stream 17 """
        if self.y_ipt==1:
            self.s[17].y = 1
            self.s[17].P = self.t2b.Pout
            self.s[17].T = self.t2b.Tout
            self.s[17].mdot = self.t2b.mdot
            self.s[17].x = self.t2b.x_out
        else:
            self.s[17].y = 0
        """ MS 2 """
        if self.y_ipt==1:
            self.ms2.y = 1
            self.ms2.P = self.s[17].P
            self.ms2.T = self.s[17].T
            self.ms2.mdot = self.s[17].mdot
            self.ms2.x_in = self.s[17].x
        else:
            self.ms2.y = 0
        self.ms2.calc_MS()
        """ Stream 18 """
        if self.ms2.y==1:
            self.s[18].y = 1
            self.s[18].P = self.ms2.P
            self.s[18].T = self.ms2.T
            self.s[18].mdot = self.ms2.mdot_L
            self.s[18].x = 0.0
        else:
            self.s[18].y = 0
        """ Stream 19 """
        if self.y_ipt==1 and self.y_rh2==0:
            self.s[19].y = 1
            self.s[19].P = self.ms2.P
            # Add to T for the sake of h_pT(), since this stream skip RH 2
            self.s[19].T = self.ms2.T + (1e-10)
            self.s[19].mdot = self.ms2.mdot_V
            self.s[19].x = 1.0
        else:
            self.s[19].y = 0
        """ Stream 20 """
        if self.y_ipt==1 and self.y_rh2==1:
            self.s[20].y = 1
            self.s[20].P = self.ms2.P
            self.s[20].T = self.ms2.T
            self.s[20].mdot = self.ms2.mdot_V
            self.s[20].x = 1.0
        else:
            self.s[20].y = 0
        """ RH 2 """
        if self.y_rh2 == 1:
            self.rh2.y = 1
            id1 = 0 # Fake ID
            if self.y_s4 == 1:
                id1 = 4
            elif self.y_s5 == 1:
                id1 = 5
            elif self.y_s14 == 1:
                id1 = 14
            self.rh2.Pin1 = self.s[id1].P
            self.rh2.Tin1 = self.s[id1].T
            self.rh2.mdot1 = self.s[id1].mdot
            self.rh2.x_in1 = self.s[id1].x
            self.rh2.Satd_in1 = False
            self.rh2.Pin2 = self.s[20].P
            self.rh2.Tin2 = self.s[20].T
            self.rh2.mdot2 = self.s[20].mdot
            self.rh2.Satd_in2 = True
        else:
            self.rh2.y = 0
        self.rh2.calc_RH()
        """ Stream 15 """
        if self.y_rh2==1 and self.y_s14==1:
            self.s[15].y = 1
            self.s[15].P = self.rh2.Pout1
            self.s[15].T = self.rh2.Tout1
            self.s[15].mdot = self.rh2.mdot1
            self.s[15].x = self.rh2.x_out1
        else:
            self.s[15].y = 0
        """ Stream 16 """
        if self.y_rh2==1 and self.y_s14==0:
            self.s[16].y =1
            self.s[16].P = self.rh2.Pout1
            self.s[16].T = self.rh2.Tout1
            self.s[16].mdot = self.rh2.mdot1
            self.s[16].x = self.rh2.x_out1
        else:
            self.s[16].y = 0
        """ Stream 21 """
        if self.y_rh2==1:
            self.s[21].y = 1
            self.s[21].P = self.rh2.Pout2
            self.s[21].T = self.rh2.Tout2
            self.s[21].mdot = self.rh2.mdot2
            self.s[21].x = 1.0
        else:
            self.s[21].y = 0
        """ Stream 22 """
        id_in = 0 # Fake ID
        if self.s[10].y == 1:
            id_in = 10
        elif self.s[13].y == 1:
            id_in = 13
        elif self.s[19].y == 1:
            id_in = 19
        elif self.s[21].y == 1:
            id_in = 21
        self.s[22].P = self.s[id_in].P
        self.s[22].T = self.s[id_in].T
        self.s[22].mdot = self.s[id_in].mdot
        self.s[22].x = self.s[id_in].x
        """ Turbine LPT 1 """
        self.t3.Pin = self.s[22].P
        self.t3.Tin = self.s[22].T
        self.t3.mdot = self.s[22].mdot / 3.0
        self.t3.x_in = self.s[22].x
        self.t3.calc_Turb()
        """ Turbine LPT 2 """
        self.t4.Pin = self.s[22].P
        self.t4.Tin = self.s[22].T
        self.t4.mdot = self.s[22].mdot / 3.0
        self.t4.x_in = self.s[22].x
        self.t4.calc_Turb()
        """ Turbine LPT 3 """
        self.t5.Pin = self.s[22].P
        self.t5.Tin = self.s[22].T
        self.t5.mdot = self.s[22].mdot / 3.0
        self.t5.x_in = self.s[22].x
        self.t5.calc_Turb()
        """ Stream 23 """
        self.s[23].P = self.t3.Pout
        self.s[23].T = self.t3.Tout
        self.s[23].mdot = self.t3.mdot+self.t4.mdot+self.t5.mdot
        self.s[23].x = self.t3.x_out
        """ Condenser """
        self.cond.Pin = self.s[23].P
        self.cond.Tin = self.s[23].T
        self.cond.mdot = self.s[23].mdot
        self.cond.x_in = self.s[23].x
        self.cond.calc_Condenser()
        """ Stream 24 """
        self.s[24].P = self.cond.Pout
        self.s[24].T = self.cond.Tout
        self.s[24].mdot = self.cond.mdot
        self.s[24].x = self.cond.x_out
        """ Pump 5 """
        self.p5.Pin = self.s[24].P
        self.p5.Tin = self.s[24].T
        self.p5.mdot = self.s[24].mdot
        Po_p5 = 0.0 # Fake pressure
        if self.y_ipt==0:
            Po_p5 = self.Po_t1c
        elif self.y_ipt==1:
            Po_p5 = self.Po_t2b
        self.p5.Pout = Po_p5
        self.p5.calc_Pump()
        """ Stream 25 """
        if self.y_ipt==0:
            self.s[25].y = 1
            self.s[25].P = self.p5.Pout
            self.s[25].T = self.p5.Tout
            self.s[25].mdot = self.p5.mdot
            self.s[25].x = 0.0
        else:
            self.s[25].y = 0
        """ Stream 26 """
        if self.y_ipt==1:
            self.s[26].y = 1
            self.s[26].P = self.p5.Pout
            self.s[26].T = self.p5.Tout
            self.s[26].mdot = self.p5.mdot
            self.s[26].x = 0.0
        else:
            self.s[26].y = 0
        """ FWH 4 """
        if self.y_ipt==1:
            self.fwh4.y = 1
            self.fwh4.Pin1 = self.s[18].P
            self.fwh4.Tin1 = self.s[18].T
            self.fwh4.mdot1 = self.s[18].mdot
            self.fwh4.x_in1 = self.s[18].x
            self.fwh4.Pin2 = self.s[26].P
            self.fwh4.Tin2 = self.s[26].T
            self.fwh4.mdot2 = self.s[26].mdot
        else:
            self.fwh4.y = 0
        self.fwh4.calc_FWH()
        """ Stream 27 """
        if self.fwh4.y==1:
            self.s[27].y = 1
            self.s[27].P = self.fwh4.Pout
            self.s[27].T = self.fwh4.Tout
            self.s[27].mdot = self.fwh4.mdot
            self.s[27].x = self.fwh4.x_out
        else:
            self.s[27].y = 0
        """ Pump 4 """
        if self.fwh4.y==1:
            self.p4.y = 1
            self.p4.Pin = self.s[27].P
            self.p4.Tin = self.s[27].T
            self.p4.mdot = self.s[27].mdot
            Po_p4 = 0.0 # Fake pressure
            if self.s[8].y==1 or self.s[15].y==1:
                if self.s[8].y==1:
                    Po_p4 = self.s[8].P
                else:
                    Po_p4 = self.s[15].P
            else:
                if self.s[7].y==1:
                    Po_p4 = self.s[7].P
                elif self.s[16].y==1:
                    Po_p4 = self.s[16].P
            self.p4.Pout = Po_p4
        else:
            self.p4.y = 0
        self.p4.calc_Pump()
        """ Stream 28 """
        if self.p4.y==1:
            if self.s[8].y==0 and self.s[15].y==0:
                self.s[28].y = 1
                self.s[28].P = self.p4.Pout
                self.s[28].T = self.p4.Tout
                self.s[28].mdot = self.p4.mdot
                self.s[28].x = 0.0
            else:
                self.s[28].y = 0
        else:
            self.s[28].y = 0
        """ Stream 29 """
        if self.p4.y==1:
            if self.s[8].y==1 or self.s[15].y==1:
                self.s[29].y = 1
                self.s[29].P = self.p4.Pout
                self.s[29].T = self.p4.Tout
                self.s[29].mdot = self.p4.mdot
                self.s[29].x = 0.0
            else:
                self.s[29].y = 0
        else:
            self.s[29].y = 0
        """ FWH 3 """
        if self.s[8].y==1 or self.s[15].y==1:
            self.fwh3.y = 1
            id1 = 0 # Fake ID
            if self.s[8].y==1:
                id1 = 8
            else:
                id1 = 15
            self.fwh3.Pin1 = self.s[id1].P
            self.fwh3.Tin1 = self.s[id1].T
            self.fwh3.mdot1 = self.s[id1].mdot
            self.fwh3.x_in1 = self.s[id1].x
            self.fwh3.Pin2 = self.s[29].P
            self.fwh3.Tin2 = self.s[29].T
            self.fwh3.mdot2 = self.s[29].mdot
        else:
            self.fwh3.y = 0
        self.fwh3.calc_FWH()
        """ Stream 30 """
        if self.fwh3.y==1:
            self.s[30].y = 1
            self.s[30].P = self.fwh3.Pout
            self.s[30].T = self.fwh3.Tout
            self.s[30].mdot = self.fwh3.mdot
            self.s[30].x = self.fwh3.x_out
        else:
            self.s[30].y = 0
        """ Pump 3 """
        if self.fwh3.y==1:
            self.p3.y = 1
            self.p3.Pin = self.s[30].P
            self.p3.Tin = self.s[30].T
            self.p3.mdot = self.s[30].mdot
            Po_p3 = 0.0 # Fake pressure
            if self.s[7].y==1:
                Po_p3 = self.s[7].P
            elif self.s[16].y==1:
                Po_p3 = self.s[16].P
            self.p3.Pout = Po_p3
        else:
            self.p3.y = 0
        self.p3.calc_Pump()
        """ Stream 31 """
        if self.p3.y==1:
            self.s[31].y = 1
            self.s[31].P = self.p3.Pout
            self.s[31].T = self.p3.Tout
            self.s[31].mdot = self.p3.mdot
            self.s[31].x = 0.0
        else:
            self.s[31].y = 0
        """ FWH 2 """
        id1 = 0 # Fake ID
        if self.s[7].y==1:
            id1 = 7
        elif self.s[16].y==1:
            id1 = 16
        id2 = 0 # Fake ID
        if self.s[25].y==1:
            id2 = 25
        elif self.s[28].y==1:
            id2 = 28
        elif self.s[31].y==1:
            id2 = 31
        self.fwh2.Pin1 = self.s[id1].P
        self.fwh2.Tin1 = self.s[id1].T
        self.fwh2.mdot1 = self.s[id1].mdot
        self.fwh2.x_in1 = self.s[id1].x
        self.fwh2.Pin2 = self.s[id2].P
        self.fwh2.Tin2 = self.s[id2].T
        self.fwh2.mdot2 = self.s[id2].mdot
        self.fwh2.calc_FWH()
        """ Stream 32 """
        if self.s[3].y==0:
            self.s[32].y = 1
            self.s[32].P = self.fwh2.Pout
            self.s[32].T = self.fwh2.Tout
            self.s[32].mdot = self.fwh2.mdot
            self.s[32].x = self.fwh2.x_out
        else:
            self.s[32].y = 0
        """ Stream 33 """
        if self.s[3].y==1:
            self.s[33].y = 1
            self.s[33].P = self.fwh2.Pout
            self.s[33].T = self.fwh2.Tout
            self.s[33].mdot = self.fwh2.mdot
            self.s[33].x = self.fwh2.x_out
        else:
            self.s[33].y = 0
        """ Pump 2 """
        if self.s[33].y==1:
            self.p2.y = 1
            self.p2.Pin = self.s[33].P
            self.p2.Tin = self.s[33].T
            self.p2.mdot = self.s[33].mdot
            self.p2.Pout = self.Po_t1a
        else:
            self.p2.y = 0
        self.p2.calc_Pump()
        """ Stream 34 """
        if self.p2.y==1:
            self.s[34].y = 1
            self.s[34].P = self.p2.Pout
            self.s[34].T = self.p2.Tout
            self.s[34].mdot = self.p2.mdot
            self.s[34].x = 0.0
        else:
            self.s[34].y = 0
        """ FWH 1 """
        if self.s[3].y==1:
            self.fwh1.y = 1
            self.fwh1.Pin1 = self.s[3].P
            self.fwh1.Tin1 = self.s[3].T
            self.fwh1.mdot1 = self.s[3].mdot
            self.fwh1.x_in1 = self.s[3].x
            self.fwh1.Pin2 = self.s[34].P
            self.fwh1.Tin2 = self.s[34].T
            self.fwh1.mdot2 = self.s[34].mdot
        else:
            self.fwh1.y = 0
        self.fwh1.calc_FWH()
        """ Stream 35 """
        if self.fwh1.y==1:
            self.s[35].y = 1
            self.s[35].P = self.fwh1.Pout
            self.s[35].T = self.fwh1.Tout
            self.s[35].mdot = self.fwh1.mdot
            self.s[35].x = self.fwh1.x_out
        else:
            self.s[35].y = 0
        """ Pump 1 """
        id_in = 0 # Fake ID
        if self.s[32].y==1:
            id_in = 32
        elif self.s[35].y==1:
            id_in = 35
        self.p1.Pin = self.s[id_in].P
        self.p1.Tin = self.s[id_in].T
        self.p1.mdot = self.s[id_in].mdot
        self.p1.Pout = self.phx.Pin
        self.p1.calc_Pump()
        """ Stream 36 """
        self.s[36].P = self.p1.Pout
        self.s[36].T = self.p1.Tout
        self.s[36].mdot = self.p1.mdot
        self.s[36].x = 0.0


############################
""" PCS Constraints Func """
############################

def PCS_constraints(x,y):
    """
    Inputs:
        x = Array of x-optimization parameters
        y = Array of y-optimization parameters
    Output:
        x = Array of constraint-corrected x-optimization parameters
        y = Array of constraint-corrected y-optimization parameters
    """
    # Constraints for mchx.Tout
    if x[0] < 295.0: x[0] = 295.0
    if x[0] > 307.7: x[0] = 307.7
    # Constraints for Po_t1a
    if x[1] < 8.0: x[1] = 8.0
    if x[1] > 25.0: x[1] = 25.0
    # Constraints for mf_t1a
    if x[2] < 0.05: x[2] = 0.05
    if x[2] > 0.19: x[2] = 0.19
    # Constraints for Po_t1b
    if x[3] < 6.0: x[3] = 6.0
    if x[3] > 20.0: x[3] = 20.0
    # Constraints for mf_t1b
    if x[4] < 0.05: x[4] = 0.05
    if x[4] > 0.19: x[4] = 0.19
    # Constraints for Po_t1c
    if x[5] < 4.0: x[5] = 4.0
    if x[5] > 16.0: x[5] = 16.0
    # Constraints for Po_t2a
    if x[6] < 3.0: x[6] = 3.0
    if x[6] > 13.0: x[6] = 13.0
    # Constraints for mf_t2a
    if x[7] < 0.05: x[7] = 0.05
    if x[7] > 0.19: x[7] = 0.19
    # Constraints for Po_t2b
    if x[8] < 2.0: x[8] = 2.0
    if x[8] > 11.0: x[8] = 11.0
    """ Binary constraints ensure y-values match 1 of the 9 configurations """
    # If IPT does not exist:
    if y[0] == 0:
        y[2] = 0
    # If RH 2 does not exist:
    if y[2] == 0:
        y[3] = 0
        y[4] = 0
        y[5] = 0
    # If RH 2 exists:
    if y[2] == 1:
        # If RH 1 exists:
        if y[1] == 1:
            # If s[14], s[4], and s[5] are all 0:
            if (y[3]+y[4]+y[5]) == 0:
                if random.random() < 0.33:
                    y[3] = 1
                elif random.random() < 0.5:
                    y[4] = 1
                else:
                    y[5] = 1
            # If s[14], s[4], and s[5] are all 1:
            if (y[3]+y[4]+y[5]) == 3:
                if random.random() < 0.33:
                    y[4] = 0
                    y[5] = 0
                elif random.random() < 0.5:
                    y[3] = 0
                    y[5] = 0
                else:
                    y[3] = 0
                    y[4] = 0
            # If s[14], s[4], and s[5] are a permutation of 0,1,1:
            if (y[3]+y[4]+y[5]) == 2:
                if y[3] == 0:
                    if random.random() < 0.5:
                        y[4] = 0
                    else:
                        y[5] = 0
                elif y[4] == 0:
                    if random.random() < 0.5:
                        y[3] = 0
                    else:
                        y[5] = 0
                elif y[5] == 0:
                    if random.random() < 0.5:
                        y[3] = 0
                    else:
                        y[4] = 0
        # If RH 1 does not exist:
        if y[1] == 0:
            # then s[4] should not exist either
            y[4] = 0
            # If s[14] and s[5] ar both 0:
            if (y[3]+y[5]) == 0:
                if random.random() < 0.5:
                    y[3] = 1
                else:
                    y[5] = 1
            # If s[14] and s[5] are both 1:
            if (y[3]+y[5]) == 2:
                if random.random() < 0.5:
                    y[3] = 0
                else:
                    y[5] = 0
    """
    Given the y-values, make sure x-value P's are spaced out enough to
    avoid causing pinch points in the reheaters. 
    """
    # Constraints for Turbine Combo: t1a(o),t1b(x),t1c(o),t2a(x),t2b(x)
    ''' 1 - red '''
    if y[0]==0 and y[1]==1 and y[2]==0 and y[3]==0 and y[4]==0 and y[5]==0:
        x[1],x[5] = DeltaP(x[1],x[5],4.0)
        if x[3]>=x[1] or x[3]<=x[5]: x[3] = x[5] + 1.0
    # Constraints for Turbine Combo: t1a(x),t1b(x),t1c(o),t2a(x),t2b(x)
    ''' 2 - firebrick '''
    if y[0]==0 and y[1]==0 and y[2]==0 and y[3]==0 and y[4]==0 and y[5]==0:
        if x[3]<=x[5]: x[3] = x[5] + 2.0
        if x[1]<=x[3]: x[1] = x[3] + 2.0
        if x[6]>=x[5]: x[6] = x[5] - 2.0
        if x[8]>=x[6]: x[8] = x[6] - 2.0
    # Constraints for Turbine Combo: t1a(o),t1b(x),t1c(o),t2a(x),t2b(o)
    ''' 3 - darkgreen '''
    if y[0]==1 and y[1]==1 and y[2]==1 and y[3]==0 and y[4]==1 and y[5]==0:
        x_t = np.array((2.0,x[8],4.0,x[5],4.0,x[1]))
        x_t = DeltaP4(x_t)
        x[8],x[5],x[1] = x_t[1],x_t[3],x_t[5]
        if x[3]>=x[1] or x[3]<=x[5]: x[3] = x[5] + 1.0
        if x[6]>=x[5] or x[6]<=x[8]: x[6] = x[8] + 1.0
    ''' 9 - yellow '''
    if y[0]==1 and y[1]==1 and y[2]==0 and y[3]==0 and y[4]==0 and y[5]==0:
        x_t = np.array((2.0,x[8],4.0,x[5],4.0,x[1]))
        x_t = DeltaP4(x_t)
        x[8],x[5],x[1] = x_t[1],x_t[3],x_t[5]
        if x[3]>=x[1] or x[3]<=x[5]: x[3] = x[5] + 1.0
        if x[6]>=x[5] or x[6]<=x[8]: x[6] = x[8] + 1.0
    # Constraints for Turbine Combo: t1a(o),t1b(o),t1c(o),t2a(x),t2b(o)
    ''' 4 - purple '''
    if y[0]==1 and y[1]==1 and y[2]==1 and y[3]==0 and y[4]==0 and y[5]==1:
        x_t = np.array((2.0,x[8],2.0,x[5],2.0,x[3],2.0,x[1]))
        x_t = DeltaP4(x_t)
        x[8],x[5],x[3],x[1] = x_t[1],x_t[3],x_t[5],x_t[7]
        # Make sure the dummy value for x[6] is b/t Po_t1c and Po_t2b
        if x[6]>=x[5] or x[6]<=x[8]: x[6] = x[8] + 1.0
    # Constraints for Turbine Combo: t1a(x),t1b(o),t1c(o),t2a(x),t2b(o)
    ''' 5 - deeppink '''
    if y[0]==1 and y[1]==0 and y[2]==1 and y[3]==0 and y[4]==0 and y[5]==1:
        x_t = np.array((2.0,x[8],4.0,x[5],2.0,x[3]))
        x_t = DeltaP4(x_t)
        x[8],x[5],x[3] = x_t[1],x_t[3],x_t[5]
        if x[1]<=x[3]: x[1] = x[3] + 1.0
        if x[6]>=x[5] or x[6]<=x[8]: x[6] = x[8] + 1.0
    # Constraints for Turbine Combo: t1a(o),t1b(x),t1c(o),t2a(o),t2b(o)
    ''' 6 - blue '''
    if y[0]==1 and y[1]==1 and y[2]==1 and y[3]==1 and y[4]==0 and y[5]==0:
        x_t = np.array((x[1],4.0,x[5],2.0,x[6],4.0,x[8]))
        x_t = DeltaP3(x_t)
        x[1],x[5],x[6],x[8] = x_t[0],x_t[2],x_t[4],x_t[6]
        if x[3]>=x[1] or x[3]<=x[5]: x[3] = x[5] + 1.0
    # Constraints for Turbine Combo: t1a(x),t1b(x),t1c(o),t2a(o),t2b(o)
    ''' 7 - cyan '''
    if y[0]==1 and y[1]==0 and y[2]==1 and y[3]==1 and y[4]==0 and y[5]==0:
        x_t = np.array((x[5],2.0,x[6],4.0,x[8]))
        x_t = DeltaP3(x_t)
        x[5],x[6],x[8] = x_t[0],x_t[2],x_t[4]
        if x[3]<=x[5]: x[3] = x[5] + 1.0
        if x[1]<=x[3]: x[1] = x[3] + 1.0
    # Constraints for Turbine Combo: t1a(x),t1b(x),t1c(o),t2a(x),t2b(o)
    ''' 8 - orange '''
    if y[0]==1 and y[1]==0 and y[2]==0 and y[3]==0 and y[4]==0 and y[5]==0:
        x[5],x[8] = DeltaP(x[5],x[8],4.0)
        if x[3]<=x[5]: x[3] = x[5] + 2.0
        if x[1]<=x[3]: x[1] = x[3] + 2.0
        if x[6]>=x[5] or x[6]<=x[8]: x[6] = x[8] + 1.0
    # To correct any "NaN" problems with Pump6, make sure that IPT.Pout > 1.0 bar
    if x[8]<1.0: x[8] = 1.001
    # Return corrected x,y-values
    return x,y


########################
""" PCS Maximin Func """
########################

def PCS_maximin(opt_list):
    """
    Inputs:
        opt_list = Population's 'opt_list' array of PCS_Option's
    Outputs:
        sorted_list = Array of fitness-sorted PCS_Option's
    """
    # Initialize parameters
    n = np.size(opt_list)
    fmm_min = 1000
    min_index = 50
    # Rotate through each PCS_Option in opt_list
    for i in range(n):
        # Initialize array of minimum differences
        j_mins = np.empty(0)
        # Rotate through each PCS_Option in opt_list (except j = i)
        for j in range(n):
            if j == i:
                None
            else:
                # Find min[k=1,2](f_k(x_i) - f_k(x_j))
                k_diff1 = opt_list[i].inv_eff() - opt_list[j].inv_eff()
                k_diff2 = opt_list[i].cost() - opt_list[j].cost()
                k_min = min(k_diff1,k_diff2)
                j_mins = np.append(j_mins,k_min)
        # Find the max of the j_mins and assign new f_maximin to PCS_Option i
        i_max = max(j_mins)
        opt_list[i].fmm = i_max
        # Keep track of the smallest f_maximin
        if i_max < fmm_min:
            fmm_min = i_max
            min_index = i
    # Initialize the maximin-sorted list of PCS_Option's
    sorted_list = np.array(opt_list[min_index])
    # Re-order the list of PCS_Option's in ascending "maximin" order
    for count in range(n-1):
        fmm_next = 1000
        # Rotate through the PCS_Option's of opt_list
        for i in range(n):
            # Find the next smallest f_maximin
            if i != min_index:
                # If current PCS_Option's fmm-value is less than fmm_next
                if opt_list[i].fmm < fmm_next and opt_list[i].fmm >= fmm_min:
                    # If it equals the previous minimum
                    if opt_list[i].fmm == fmm_min and i > min_index:
                        index_next = i
                        fmm_next = opt_list[i].fmm
                        break
                    else:
                        if opt_list[i].fmm == fmm_min and i < min_index:
                            None
                        else:
                            index_next = i
                            fmm_next = opt_list[i].fmm
        # Add the next best PCS_Option to the sorted list
        sorted_list = np.append(sorted_list,opt_list[index_next])
        fmm_min = fmm_next
        min_index = index_next
    # Return the maximin-sorted list of PCS_Option's
    return sorted_list


######################
""" PCS Graph Func """
######################

def PCS_Graph_Data(wb,opt_list,tab,dt):
    """
    Inputs:
        wb = Excel workbook for collecting data
        opt_list = Population's 'opt_list' array of PCS_Option's
        tab = Excel tab for collecting PCS seed population data
        dt = computational runtime
    Actions:
        Creates/Prints objective function graph in Python.
        Pastes that graph in Excel.
        Pastes important PCS_Option data in Excel.
    """
    # Declare the runtime and current time
    print ("Number of Iterations = -0") # Neg-Zero b/c it's before Iter. 0
    m = time.localtime()
    if m[3]<=12:
        hr = m[3]
        if m[3]==0: hr = 12
        ap = "AM"
    if m[3]>12: 
        hr = m[3]-12
        ap = "PM"
    print ("dt =",dt/60.0,"min  /  Time =",hr,":","%02.0f"%m[4],":","%02.0f"%m[5],ap)
    # Graph in Python
    fig = plt.figure(figsize=(6,4))
    ax = fig.gca()
    plt.xlim([0.33,0.37])
    plt.ylim([0,13])
    ax.set_xticks(np.arange(0.33,0.37,0.005))
    ax.set_yticks(np.arange(0,13,2))
    for i in range(np.size(opt_list)):
        plt.scatter(opt_list[i].eff(),opt_list[i].cost(),s=10,
                    c=opt_list[i].color(),edgecolors=opt_list[i].color())
    plt.grid()
    plt.xlabel('Efficiency')
    plt.ylabel('Cost ($1x10^9)')
    plt.show()
    # Paste graph in Excel    
    wb.sheets[tab].pictures.add(fig,name='pcs_graph')
    wb.sheets[tab].pictures[0].left = 550
    wb.sheets[tab].pictures[0].top = 80
    wb.sheets[tab].pictures[0].height = 211
    wb.sheets[tab].pictures[0].width = 316
    # Paste all pertinent PCS_Option data
    wb.sheets[tab].range('B1').value = dt/60.0 #min
    col = 3
    n_a = "-"
    for k in range(np.size(opt_list)):
        wb.sheets[tab].range(2,col).value = k + 1
        wb.sheets[tab].range(3,col).value = opt_list[k].eff()
        wb.sheets[tab].range(4,col).value = opt_list[k].cost()
        wb.sheets[tab].range(5,col).value = opt_list[k].fmm
        wb.sheets[tab].range(6,col).value = opt_list[k].color()
        wb.sheets[tab].range(7,col).value = opt_list[k].y_ipt
        wb.sheets[tab].range(8,col).value = opt_list[k].y_rh1
        wb.sheets[tab].range(9,col).value = opt_list[k].y_rh2
        wb.sheets[tab].range(10,col).value = opt_list[k].y_s14
        wb.sheets[tab].range(11,col).value = opt_list[k].y_s4
        wb.sheets[tab].range(12,col).value = opt_list[k].y_s5
        wb.sheets[tab].range(13,col).value = opt_list[k].phx.Tout
        wb.sheets[tab].range(14,col).value = opt_list[k].phx.Pout
        wb.sheets[tab].range(15,col).value = opt_list[k].phx.Tin
        wb.sheets[tab].range(16,col).value = opt_list[k].phx.Pin
        wb.sheets[tab].range(17,col).value = opt_list[k].phx.mdot
        if opt_list[k].y_rh1 == 1:
            wb.sheets[tab].range(18,col).value = opt_list[k].t1a.Pout
            wb.sheets[tab].range(19,col).value = opt_list[k].mf_t1a
        else:
            wb.sheets[tab].range(18,col).value = n_a
            wb.sheets[tab].range(19,col).value = n_a
        if opt_list[k].y_s5 == 1:
            wb.sheets[tab].range(20,col).value = opt_list[k].t1b.Pout
            wb.sheets[tab].range(21,col).value = opt_list[k].mf_t1b
        else:
            wb.sheets[tab].range(20,col).value = n_a
            wb.sheets[tab].range(21,col).value = n_a
        wb.sheets[tab].range(22,col).value = opt_list[k].t1c.Pout
        if opt_list[k].y_s14 == 1:
            wb.sheets[tab].range(23,col).value = opt_list[k].t2a.Pout
            wb.sheets[tab].range(24,col).value = opt_list[k].mf_t2a
        else:
            wb.sheets[tab].range(23,col).value = n_a
            wb.sheets[tab].range(24,col).value = n_a
        if opt_list[k].y_ipt == 1:
            wb.sheets[tab].range(25,col).value = opt_list[k].t2b.Pout
        else:
            wb.sheets[tab].range(25,col).value = n_a
        wb.sheets[tab].range(26,col).value = opt_list[k].t3.Pout
        wb.sheets[tab].range(27,col).value = opt_list[k].t1.W
        if opt_list[k].y_ipt == 1:
            wb.sheets[tab].range(28,col).value = opt_list[k].t2.W
        else:
            wb.sheets[tab].range(28,col).value = n_a
        wb.sheets[tab].range(29,col).value = opt_list[k].t3.W
        wb.sheets[tab].range(30,col).value = opt_list[k].t4.W
        wb.sheets[tab].range(31,col).value = opt_list[k].t5.W
        wb.sheets[tab].range(32,col).value = opt_list[k].p1.W
        if opt_list[k].p2.y == 1:
            wb.sheets[tab].range(33,col).value = opt_list[k].p2.W
        else:
            wb.sheets[tab].range(33,col).value = n_a
        if opt_list[k].p3.y == 1:
            wb.sheets[tab].range(34,col).value = opt_list[k].p3.W
        else:
            wb.sheets[tab].range(34,col).value = n_a
        if opt_list[k].p4.y == 1:
            wb.sheets[tab].range(35,col).value = opt_list[k].p4.W
        else:
            wb.sheets[tab].range(35,col).value = n_a
        wb.sheets[tab].range(36,col).value = opt_list[k].p5.W
        # Increment the column number between options
        col = col + 1
    # No need to return anything
    return None