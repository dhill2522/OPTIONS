"""
Created on Mon Nov 05 03:52:36 2018
@author: Paul
"""

### Boiler-Plate ###
import matplotlib.pylab as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import scipy as sp
from numpy import random
import time
import csv

from Class1_Eq import *
from Func import *


""" Change this value when changed in restart .i files """
global t_final
t_final = 10000 # seconds
global ss_fail_penalty
ss_fail_penalty = 700
global cost_multiplier_for_nucl_safety_grade
cost_multiplier_for_nucl_safety_grade = 5.0

###########################################################################
"""""""""   Tri-System Option Class   """"""""" ###########################
###########################################################################

class Option:
    """
    Inputs:
        x1 = Zion core loop x-optimization parameters
        x2 = PERCS loop x-optimization parameters
        x3 = PCS superstructure x-optimization parameters
        y = PCS superstructure y-optimization parameters
    Parameters:
        *Individual optimization parameters (explained in __init__() function)
        Core Loop:
            cards = Array of RELAP5 card numbers with core loop value changes
            i_vals = Array of column numbers for core loop value changes
            vals = Array of new values for core loop value changes
            T_fuel_cent_max = Maximum fuel centerline temperature (constraint)
            T_clad_surf_max = Maximum cladding surface temperature (constraint)
            MDNBR = Minimum departure from nucleate boiling ratio (constraint)
            T_f_over_max = [Boolean] Did fuel temperature go over the max?
            T_clad_surf_max = [Boolean] Did cladding temperature go over the max?
            MDNBR_below_1 = [Boolean] Did MDNBR go below 1.0?
            peanlized = [Boolean] Did the core loop receive a penalty?
            failed = [Boolean] Did the RELAP5 core loop model fail early?
            csvFileLocation = [String] Core's PyPost results file location
            *Parameters for T, P, m_dot, H, & x_e core data from PyPost
            k_eff = Effective multiplication factor per neutron cycle in core
            rho_0 = Initial reactivity of the core
            Bc = Cycle burn-up of the fuel [EFPD = effective full-power days]
            nBc = Discharge burn-up of the fuel
            cost_RCPs = Capital cost of RCPs
            op_cost_RCPs = Operating cost of RCPs (40 yrs)
            cost_total_fuel = Cost of UO2 fuel (40 yrs)
        PERCS Loop:
            list_card = Array of RELAP5 card numbers with PERCS value changes
            list_i_change = Array of column numbers for PERCS value changes
            list_change = Array of new values for PERCS value changes
            len_diff_717 = Parameter used to calculate length of Pipe 717
            n_tubes = Number of tubes w/in PERCS tank
            m_MgCO3 = Mass of Magnesium Carbonate w/in PERCS tank
            T_over_620 = [Boolean] Did the core outlet T go above 620K?
            T_over_635 = [Boolean] Did the core outlet T go above 635K?
            csvFileLocation2 = [String] PERCS's PyPost results file location
            *Parameters for T & alpha PERCS data from PyPost
            PERCS_failed = [Boolean] Did the PERCS RELAP5 model fail early?
            PERCS_penalty = [Boolean] Did the PERCS receive a penalty?
            cost_penalty = Multaplicative cost penalty if 'PERCS_failed' = TRUE
            ss_fail = [Boolean] Redundant of Core's 'failed'
            p716, p717 = Pipes 716 & 717 (for cost purposes)
            support = Support structure for PERCS tank (for cost purposes)
            hx = Fake heat exchanger (for cost purposes)
            tank = PERCS tank (for cost purposes)
            chemical = MgCO3 in tank (for cost purposes)
        PCS Loop:
            pinch_point = [Boolean] 
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
        Objective Functions:
            W_rcp = Core Obj. 1 - Total work of RCPs
            cost_1 = Core Obj. 2 - Total core loop costs
            obj_1_1 = Normalized W_rcp
            obj_1_2 = Normalized cost_1
            fmm_1 = Maximin fitness value for core loop
            cost_2 = PERCS Obj. 1 - Total PERCS equipment cost
            dT_int = PERCS Obj. 2 - Integral of deviation of core outlet T
            alpha = PERCS Obj. 3 - Consumption of MgCO3
            obj_2_1 = Normalized cost_2
            obj_2_2 = Normalized dT_int
            obj_2_3 = Normalized alpha
            fmm_2 = Maximin fitness value for PERCS loop
            color = [String] PCS superstructure color/configuration
            eff = PCS Obj. 1 - Thermodynamic Efficiency
            cost_3 = PCS Obj. 2 - Total PCS equipment cost
            obj_3_1 = Normalized eff
            obj_3_2 = Normalized cost_3
            fmm_3 = Maximin fitness value for PCS loop
            obj_fmm_1 = Normalized fmm_1
            obj_fmm_2 = Normalized fmm_2
            obj_fmm_3 = Normalized fmm_3
            fmm_o = Overall Maximin fitness value
    Functions:
        init_ZION_calcs() - Fills arrays to make core loop RELAP5 value changes
        init_PERCS_calcs() - Fills arrays to make PERCS RELAP5 value changes
        final_ZION_calcs() - Grabs PyPost data, Performs final core loop calcs
        final_PERCS_calcs() - Grabs PyPost data, Performs final PERCS calcs
        Alpha_calcs() - Grabs alpha PyPost data, Calcs overall Alpha
        PCS_SS_calcs() - Calls solve_PCS(), Performs final PCS calcs
        solve_PCS() - Fills out PCS superstructure & converges the cycle
    """
    def __init__(self,x1_in,x2_in,x3_in,y_in):
        self.opt_ID = 0
        self.last_sec_penalty = False
        # Define the x- and y-optimization parameter arrays
        self.x1 = x1_in  # ZION x-opt parameters
        self.x2 = x2_in  # PERCS x-opt parameters
        self.x3 = x3_in  # PCS x-opt parameters
        self.y = y_in    # PCS y-opt parameters
        # Further define the ZION Core loop opt. parameters
        self.R_f = self.x1[0]    # ft (radius of fuel per pin)
        self.H_fuel = self.x1[1] # ft (height of fuel pins)
        self.Dh_00 = self.x1[2]  # ft (hydraulic D of pipes _00)
        self.Dh_12 = self.x1[3]  # ft (hydraulic D of pipes _12)
        self.Dh_14 = self.x1[4]  # ft (hydraulic D of pipes _14)
        # Further define the PERCS loop opt. parameters
        self.R_tank = self.x2[0] # ft (radius of PERCS HX tank)
        self.pitch = self.x2[1]  # ft (pitch b/t tubes in PERCS)
        self.D_h = self.x2[2]    # ft (hydraulic D of tubes)
        self.th = self.x2[3]     # ft (thickness of tubes)
        self.Len = self.x2[4]    # ft (length of tubes / height of tank)
        self.elev = self.x2[5]   # ft (height diff. b/t core outlet & PERCS inlet)
        # Further define the PCS superstructure x-opt. parameters
        self.To_PHX = self.x3[0] # degC
        self.Po_t1a = self.x3[1] # bar
        self.mf_t1a = self.x3[2]
        self.Po_t1b = self.x3[3] # bar
        self.mf_t1b = self.x3[4]
        self.Po_t1c = self.x3[5] # bar
        self.Po_t2a = self.x3[6] # bar
        self.mf_t2a = self.x3[7]
        self.Po_t2b = self.x3[8] # bar
        # Further define the PCS superstructure y-opt. parameters
        self.y_ipt = self.y[0]   # IPT
        self.y_rh1 = self.y[1]   # RH 1
        self.y_rh2 = self.y[2]   # RH 2
        self.y_s14 = self.y[3]   # s[14]
        self.y_s4  = self.y[4]   # s[4]
        self.y_s5  = self.y[5]   # s[5]
        ################################
        """ Init stuff for ZION Core """
        ################################
        # Initialize card, i_change, and change lists for ZION
        self.cards = np.empty(119,dtype='<U32')
        self.i_vals = np.zeros(119,dtype=int)
        self.vals = np.zeros(119)
        # Initiate the Booleans that tracks thermal design limit violations
        self.T_fuel_cent_max = 2100 # degC
        self.T_clad_surf_max = 348 # degC
        self.MDNBR = 0
        self.T_f_over_max = False
        self.T_c_over_max = False
        self.MDNBR_below_1 = False
        self.penalized = False
        self.failed = False
        # Parameter data grabbed from .csv files using PyPost
        self.csvFileLocation = 'None'
        self.T_106 = 0.0 # degC
        self.T_110 = 0.0 # degC
        self.P_106 = 0.0 # bar
        self.P_110 = 0.0 # bar
        self.P_335 = np.zeros(6) # MPa
        self.P_p_out = 0.0 # bar
        self.m_dot_100 = 0.0 # kg/s
        self.m_dot_335 = 0.0 # kg/s
        self.m_dot_400 = 0.0 # kg/s
        self.m_dot_600 = 0.0 # kg/s
        self.m_dot_200 = 0.0 # kg/s
        self.H_106 = 0.0 # kJ/kg
        self.H_110 = 0.0 # kJ/kg
        self.H_335_1 = 0.0 # kJ/kg
        self.H_112_5 = 0.0 # kJ/kg
        self.H_114 = 0.0 # kJ/kg
        self.H_412_5 = 0.0 # kJ/kg
        self.H_414 = 0.0 # kJ/kg
        self.H_612_5 = 0.0 # kJ/kg
        self.H_614 = 0.0 # kJ/kg
        self.H_212_5 = 0.0 # kJ/kg
        self.H_214 = 0.0 # kJ/kg
        self.T_1336_1 = np.zeros(6) # K
        self.T_1336_17 = np.zeros(6) # K
        self.x_e_335 = np.zeros(6)
        # Other parameters that should be reported in Excel
        self.k_eff = 0.0
        self.rho_0 = 0.0
        self.Bc = 0.0 # EFPD
        self.nBc = 0.0 # yr
        # Three cost parameters that make up 'cost_1'
        self.cost_RCPs = 0.0 # $
        self.op_cost_RCPs = 0.0 # $
        self.cost_total_fuel = 0.0 # $
        ############################
        """ Init stuff for PERCS """
        ############################
        # Initialize card, i_change, and change lists for PERCS
        self.list_card = np.empty(39,dtype='<U32')
        self.list_i_change = np.zeros(39,dtype=int)
        self.list_change = np.empty(39)
        # Needed to calc the elev of Pipe 717, calc'd in Init_ZION_Calcs()
        self.len_diff_717 = 0.0 # ft
        # Initialize some stuff
        self.n_tubes = 0
        self.m_MgCO3 = 0 # kg
        # Initiate the Boolean that says whether T goes over 620 K and/or 635 K
        self.T_over_620 = False
        self.T_over_635 = False
        # Initiate the arrays for t and T and the matrix for a (alpha)
        self.csvFileLocation2 = 'None'
        self.t = np.zeros(0)
        self.T_335_6 = np.zeros(0)
        self.dT_335_6 = np.zeros(0)
        self.a_array = np.zeros(100)
        self.a = np.zeros((10,10))
        # Initiate the Boolean that says if there was a penalty for failing before t_final
        self.PERCS_failed = False
        self.PERCS_penalty = 1.0
        self.cost_penalty = 1.0
        self.ss_fail = False # Redundant
        # Initialize PERCS system equipment
        self.p716 = Pipe(self.elev)
        self.p717 = Pipe(0.0)
        self.support = Support(self.R_tank,self.Len,0.0)
        self.hx = HX()
        self.tank = Tank(self.R_tank,self.Len)
        self.chemical = Chemical(0)
        ##########################
        """ Init stuff for PCS """
        ##########################
        self.pinch_point = False
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
        ##########################################################
        """ Initiate all objective function and maximin values """
        ##########################################################
        # For ZION Core
        self.W_rcp = 0.0   # 1
        self.cost_1 = 0.0  # 2
        self.obj_1_1 = 0.0 # W_rcp
        self.obj_1_2 = 0.0 # cost_1
        self.fmm_1 = 0
        # For PERCS
        self.cost_2 = 0.0  # 1
        self.dT_int = 0.0  # 2
        self.alpha = 0.0   # 3
        self.obj_2_1 = 0.0 # cost_2
        self.obj_2_2 = 0.0 # dT_int
        self.obj_2_3 = 0.0 # consumption(alpha)
        self.fmm_2 = 0
        # For Rankine PCS
        self.color = 'black'
        self.eff = 0.0
        self.inv_eff = 0.0 # 1
        self.cost_3 = 0.0  # 2
        self.obj_3_1 = 0.0 # inv_eff
        self.obj_3_2 = 0.0 # cost_3
        self.fmm_3 = 0
        # Overall fmm-value
        self.obj_fmm_1 = 0.0 # normalized fmm_1
        self.obj_fmm_2 = 0.0 # normalized fmm_2
        self.obj_fmm_3 = 0.0 # normalized fmm_3
        self.fmm_o = 0
        #######################################################
        """ Perform the initial calculations for the Option """
        #######################################################
        self.init_ZION_calcs()
        self.init_PERCS_calcs()
    
    """
    The initial calcs take place in the init_ZION_calcs(), init_PERCS_calcs()
    function below.
    The RELAP5 and PyPost files are run from the Population.calc_Options() function.
    The obj. function and constraints calcs are run from the
    Population.final_Option_calcs() function.
    """
    
    def init_ZION_calcs(self):
        ##############################################
        """ Calcs corresponding to a change in R_f """
        ##############################################
        #-----------------------------
        """ Core Area Calculations """
        #-----------------------------
        ## Constants and Ratios
        ratio_f2m = 0.48374681 # Fuel to Moderator Ratio
        th_g = 0.002 # ft
        th_c = 0.0005 # ft
        self.n_pins = 41958.0554 # ~42,000 (value derived from RELAP5 model)
        ratio_p2D = 1.35532 # Fuel Pin Pitch to Diameter Ratio
        ## Calculations
        self.R_g = np.round(self.R_f + th_g, 4) # Gap radius [ft]
        self.R_c = np.round(self.R_f + th_g + th_c, 4) # Cladding radius [ft]
        pitch = ratio_p2D * (2.0 * self.R_c) # ft
        self.p = np.round(pitch, 4) # Fuel pin pitch [ft]
        A_f = np.pi * self.R_f**2.0 # Fuel A_c [ft^2]
        A_g = np.pi * (self.R_g**2.0 - self.R_f**2.0) # Gap A_c [ft^2]
        A_c = np.pi * (self.R_c**2.0 - self.R_g**2.0) # Cladding A_c [ft^2]
        A_p = A_f + A_g + A_c # Fuel pin A_c [ft^2]
        self.A_fuel = self.n_pins * A_f # Total fuel pin A_c [ft^2]
        self.A_gap = self.n_pins * A_g # Total gap A_c [ft^2]
        self.A_clad = self.n_pins * A_c # Total cladding A_c [ft^2]
        A_pins = self.n_pins * A_p # Total fuel pin A_c [ft^2]
        self.A_H2O = self.A_fuel / ratio_f2m # Core coolant A_c [ft^2]
        self.A_total = A_pins + self.A_H2O # Total core A_c [ft^2]
        self.A_335 = np.round(self.A_H2O,5) # Rounded core A_c [ft^2]
        A_jun_diff_335 = 2.207 # Total A_c of the baffle [ft^2]
        # Junction A_c at end of core flow segment
        self.A_jun_335 = np.round(self.A_H2O - A_jun_diff_335, 5) # ft^2
        # Hydraulic diameter of core flow segment 335 [ft]
        D_hyd = 4.0 * (pitch**2.0 - np.pi*self.R_c**2.0) / (2.0*np.pi*self.R_c)
        # Rounded hydraulic diameter of core flow segment 335
        self.Dh_335 = np.round(D_hyd,5) # ft
        # A_c of branch 336 (core above baffle) [ft^2]
        A_336 = np.round(0.272*(self.A_H2O-self.A_jun_335)+self.A_jun_335, 5)
        ## Fill the lists
        self.cards[114:117] = ['13360101','13360102','13360103']
        self.cards[78:80] = ['3350101','3350201']
        self.cards[86:88] = ['3350801','3360101']
        self.i_vals[114:117] = [3,3,3]
        self.i_vals[78:80] = [2,2]
        self.i_vals[86:88] = [3,4]
        self.vals[114:117] = [self.R_f,self.R_g,self.R_c]
        self.vals[78:80] = [self.A_335,self.A_jun_335]
        self.vals[86:88] = [self.Dh_335,A_336]
        #------------------------------------
        """ Outer Area/R_eff Calculations """
        #------------------------------------
        ## Constants and Ratios
        R_in_barrel = 6.1667 # Inner radius of the barrel [ft]
        th_baffle = 0.0937 # Thickness of the barrel [ft]
        ratio_baffle_2_core = 1.2577045 # Ratio b/t core and effective baffle
        ## Calculations
        self.R_core = np.sqrt(self.A_total/np.pi) # Radius of the entire core [ft]
        # Effective inner radius of the baffle
        Reff_in_baffle = self.R_core * ratio_baffle_2_core # ft
        # Rounded effective inner radius of the baffle
        left_bc_1335 = np.round(Reff_in_baffle, 4) # ft
        # Effective outer radius of the the baffle
        Reff_out_baffle = Reff_in_baffle + th_baffle # ft
        # Rounded effective outer radius of the baffle
        right_bc_1335 = np.round(Reff_out_baffle, 4) # ft
        # A_c taken up by the baffle
        A_baffle = np.pi * (Reff_out_baffle**2.0 - Reff_in_baffle**2.0) # ft^2
        # Total A_c of core contents (calc'd from inside out)
        A_total_plus_baffle = self.A_total + A_baffle # ft^2
        # Total A_c of core (calc'd from outside in)
        A_total_in_barrel = np.pi * R_in_barrel**2.0 # ft^2
        self.A_320_bypass = 0.0
        if (A_total_in_barrel - A_total_plus_baffle) > 18.6736:
            self.A_320_bypass = 18.6736 # ft^2
        else:
            self.A_320_bypass = A_total_in_barrel - A_total_plus_baffle # ft^2
        Dh_320 = 0.9591 # Hydraulic diameter of core bypass [ft]
        ## Fill the lists
        self.cards[106:108],self.cards[70],self.cards[77] = ['13350000','13350101'],'3200101','3200801'
        self.i_vals[106:108],self.i_vals[70],self.i_vals[77] = [6,3],2,3
        self.vals[106:108],self.vals[70],self.vals[77] = [left_bc_1335,right_bc_1335],self.A_320_bypass,Dh_320
        #################################################
        """ Calcs corresponding to a change in H_fuel """
        #################################################
        #---------------------------
        """ RPV len's and elev's """
        #---------------------------
        ## Ratios and Percentages
        # Height ratio b/t core flow segment (335) and actual fuel w/in pins
        ratio_H335_2_Hfuel = 1.1145844358
        # Length fractions per node along core flow segment (335)
        L_frac_335 = np.array((0.187389,0.1632396,0.1632396,0.1632396,0.1632396,0.1596523))
        # Length fractions per node along fuel in pins
        L_frac_pin = np.array((0.1819444,0.1819444,0.1819444,0.1819444,0.1819444,0.090278))
        ## Calculations
        # Height of core flow segment (335)
        self.H_335 = self.H_fuel * ratio_H335_2_Hfuel # ft
        # Lengths per node along core flow segment (335)
        len_335 = np.round(self.H_335 * L_frac_335, 5) # ft
        # Lengths of 'len_335' for upward-oriented RELAP5 flow segments
        Lu = [len_335[0],len_335[3],len_335[5]] # ft
        # Lengths of 'len_335' for downward-oriented RELAP5 flow segments
        Ld = [len_335[5],len_335[3],len_335[0]] # ft
        # Lengths of 'len_335' for downward-flowing RELAP5 flow segments
        nLd = [-len_335[5],-len_335[3],-len_335[0]] # ft
        len_pin = np.round(self.H_fuel * L_frac_pin, 5) # Rounded length of pin [ft]
        C_pin = 2.0*np.pi * self.R_c # Circumference of fuel pin [ft]
        # Total pin surface area on node 5
        SA_1336_5R = np.round(self.n_pins * C_pin * len_pin[4], 4) # ft^2
        # Total pin surface area on node 6
        SA_1336_6R = np.round(self.n_pins * C_pin * len_pin[5], 4) # ft^2
        ## Fill the lists
        self.cards[80:86] = ['3350301','3350302','3350303','3350701','3350702','3350703']
        self.i_vals[80:86] = [2,2,2,2,2,2]
        self.vals[80:86] = Lu+Lu
        self.cards[71:77] = ['3200301','3200302','3200303','3200701','3200702','3200703']
        self.i_vals[71:77] = [2,2,2,2,2,2]
        self.vals[71:77] = Ld+nLd
        self.cards[64:70] = ['3150301','3150302','3150303','3150701','3150702','3150703']
        self.i_vals[64:70] = [2,2,2,2,2,2]
        self.vals[64:70] = Ld+nLd
        self.cards[88:94] = ['13150501','13150502','13150503','13150601','13150602','13150603']
        self.i_vals[88:94] = [6,6,6,6,6,6]
        self.vals[88:94] = Ld+Ld
        self.cards[94:100] = ['13160501','13160502','13160503','13160601','13160602','13160603']
        self.i_vals[94:100] = [6,6,6,6,6,6]
        self.vals[94:100] = Ld+Ld
        self.cards[100:106] = ['13200501','13200502','13200503','13200601','13200602','13200603']
        self.i_vals[100:106] = [6,6,6,6,6,6]
        self.vals[100:106] = Ld+Ld
        self.cards[108:114] = ['13350501','13350502','13350503','13350601','13350602','13350603']
        self.i_vals[108:114] = [6,6,6,6,6,6]
        self.vals[108:114] = Lu+Lu
        self.cards[117:119] = ['13360601','13360602']
        self.i_vals[117:119] = [6,6]
        self.vals[117:119] = [SA_1336_5R,SA_1336_6R]
        #------------------------------
        """ PERCS p717 len and elev """
        #------------------------------
        ## Calculations
        # Deviation from original height of the fuel (for PERCS pipe 717 calc)
        self.len_diff_717 = ratio_H335_2_Hfuel * (self.H_fuel - 11.99971) # ft
        ##################################################
        """ Calcs corresponding to changes in pipe D's """
        ##################################################
        ## Calculations
        A_00 = np.round(np.pi/4.0*self.Dh_00**2.0, 3) # A_c of pipes _00 [ft^2]
        A_12 = np.round(np.pi/4.0*self.Dh_12**2.0, 3) # A_c of pipes _12 [ft^2]
        A_14 = np.round(np.pi/4.0*self.Dh_14**2.0, 3) # A_c of pipes _14 [ft^2]
        ## Fill the lists
        self.cards[0:6] = ['1000101','1000801','1020101','1020101','1040101','1040801']
        self.i_vals[0:6] = [2,3,2,9,2,3]
        self.vals[0:6] = [A_00,self.Dh_00,A_00,self.Dh_00,A_00,self.Dh_00]
        self.cards[6:10] = ['1120101','1120801','1130101','1130108']
        self.i_vals[6:10] = [2,3,2,3]
        self.vals[6:10] = [A_12,self.Dh_12,A_12,A_12]
        self.cards[10:19] = ['1130109','1140101','1140801','1160101','1160101','1161101','1162101','1180101','1180801']
        self.i_vals[10:19] = [3,2,3,2,9,4,4,2,3]
        self.vals[10:19] = [A_14,A_14,self.Dh_14,A_14,self.Dh_14,A_14,A_14,A_14,self.Dh_14]
        self.cards[19:25] = ['4000101','4000801','4120101','4120801','4130101','4130108']
        self.i_vals[19:25] = [2,3,2,3,2,3]
        self.vals[19:25] = [A_00,self.Dh_00,A_12,self.Dh_12,A_12,A_12]
        self.cards[25:34] = ['4130109','4140101','4140801','4160101','4160101','4161101','4162101','4180101','4180801']
        self.i_vals[25:34] = [3,2,3,2,9,4,4,2,3]
        self.vals[25:34] = [A_14,A_14,self.Dh_14,A_14,self.Dh_14,A_14,A_14,A_14,self.Dh_14]
        self.cards[34:40] = ['6000101','6000801','6120101','6120801','6130101','6130108']
        self.i_vals[34:40] = [2,3,2,3,2,3]
        self.vals[34:40] = [A_00,self.Dh_00,A_12,self.Dh_12,A_12,A_12]
        self.cards[40:49] = ['6130109','6140101','6140801','6160101','6160101','6161101','6162101','6180101','6180801']
        self.i_vals[40:49] = [3,2,3,2,9,4,4,2,3]
        self.vals[40:49] = [A_14,A_14,self.Dh_14,A_14,self.Dh_14,A_14,A_14,A_14,self.Dh_14]
        self.cards[49:55] = ['2000101','2000801','2120101','2120801','2130101','2130108']
        self.i_vals[49:55] = [2,3,2,3,2,3]
        self.vals[49:55] = [A_00,self.Dh_00,A_12,self.Dh_12,A_12,A_12]
        self.cards[55:64] = ['2130109','2140101','2140801','2160101','2160101','2161101','2162101','2180101','2180801']
        self.i_vals[55:64] = [3,2,3,2,9,4,4,2,3]
        self.vals[55:64] = [A_14,A_14,self.Dh_14,A_14,self.Dh_14,A_14,A_14,A_14,self.Dh_14]
    
    
    def init_PERCS_calcs(self):
        # Calc the number of tubes in PERCS
        Ac_tank = np.pi * self.R_tank**2.0 # A_c of entire tank ft^2
        Ac_hex = np.sqrt(3)/2 * self.pitch**2.0 # A_c of hexagon around tube [ft^2]
        self.n_tubes = np.round(Ac_tank / Ac_hex) # Number of PERCS tubes
        self.hx.n = self.n_tubes
        # Calc the heat transfer Surface Area in PERCS
        OD_tube = self.D_h + 2.0*self.th # Outer D of tube [ft]
        SA_tube = np.pi*OD_tube*self.Len # Surface area of tube [ft^2]
        SA_tot = SA_tube * self.n_tubes # Total surface area of tubes ft^2]
        self.hx.A = SA_tot / 10.7639 # m^2
        # Perform calcs for HX and Tank
        self.hx.calc_HX()
        self.tank.calc_Tank()
        # Calc the total cross-sectional Area of all tubes
        Ac_tube = np.pi*(self.D_h/2.0)**2 # ft^2
        Ac_tubes = np.round(Ac_tube*self.n_tubes,5) # ft^2
        # Calc the length of a single node along the tubes
        len_node = np.round((self.Len / 10.0),5) # ft
        # Calc the thickness of a single MgCO3 section (there being 10 across)
        R_hex = np.sqrt(Ac_hex/np.pi) # ft
        OR_tube = OD_tube / 2.0 # ft
        th_MgCO3 = '%.5g'%((R_hex - OR_tube)/10.0) # ft
        # Calc the heat transfer length between all tubes and MgCO3 per node
        HT_len_per_node = np.round((len_node*self.n_tubes),5) # ft
        # Calc the len and elev of Pipe 717
        self.elev_717 = np.round(-(15.62469 + self.elev - self.Len + self.len_diff_717),5) # ft  * was 15.62463
        # Now replace elev_717 values and solve Pipe 717 and Support
        self.p717.len = -self.elev_717 # ft
        self.p717.calc_Pipe()
        self.support.elev = np.round(-self.elev_717,5) # ft
        self.support.calc_Support()
        
        """ Calc the masses of MgCO3 in each of the 10 sections per axial node """
        # Define rho_MgCO3, then Calc Ac_tank, Vol_MgCO3, m_MgCO3
        rho_MgCO3 = 5.29903 # kg/ft^3
        Ac_tank = np.pi*self.R_tank**2 # ft^2
        Vol_MgCO3 = (Ac_tank - Ac_tubes)*self.Len #ft^3
        self.m_MgCO3 = rho_MgCO3 * Vol_MgCO3 # kg
        # Since we calculated the total chemical mass, calc the Chemical costs
        self.chemical.mass = self.m_MgCO3
        self.chemical.calc_Chemical()
        # Create array of radii of the all radial nodes
        radii = np.empty(11)
        for f in range(11):
            radii[f] = 0.5*self.D_h + self.th + f*float(th_MgCO3) # ft
        # Create array of annular area for all radial nodes
        areas = np.empty(10)
        for h in range(10):
            areas[h] = np.pi*radii[h+1]**2.0 - np.pi*radii[h]**2.0 # ft^2
        tot_area = sum(areas) # ft^2
        # Create array of area ratios for all radial nodes
        self.ratio = np.empty(10)
        for k in range(10):
            self.ratio[k] = areas[k]/tot_area
        # Create array of MgCO3 masses per axial node for all radial nodes
        self.masses = np.empty(10)
        for n in range(10):
            self.masses[n] = np.round(self.ratio[n]*self.m_MgCO3/10,5) # kg
        
        """ The Filling of the card, i_change, and change Lists """
        # Start filling the card, i_change, and change lists
        self.list_card[0:4] = ['7160101','7160101','7170101','7170101']
        self.list_i_change[0:4] = [3,7,3,7]
        self.list_change[0:4] = [self.elev,self.elev,self.p717.len,self.elev_717]
        # Fill lists w/ changes to PERCS tube Ac, node length, & D_h
        self.list_card[4:7] = ['7020101','7020301','7020801']
        self.list_i_change[4:7] = [2,2,3]
        self.list_change[4:7] = [Ac_tubes,len_node,self.D_h]
        # Fill lists w/ the remaining tube node lengths
        card_ = 7070301
        for k in range(9):
            self.list_card[7+k] = repr(card_)
            self.list_i_change[7+k] = 2
            self.list_change[7+k] = len_node
            card_ = card_ + 10000
        # Fill lists w/ changes to tube & MgCO3 thicknesses
        self.list_card[16:18] = ['14000101','14000102']
        self.list_i_change[16:18] = [2,2]
        self.list_change[16:18] = [self.th,th_MgCO3]
        # Fill lists w/ changes to "heat transfer length b/t all tubes and MgCO3 per node"
        self.list_card[18:20] = ['14000501','14000601']
        self.list_i_change[18:20] = [6,6]
        self.list_change[18:20] = [HT_len_per_node,HT_len_per_node]
        # Fill lists w/ changes to 9 other MgCO3 thicknesses
        card_ = 14001101
        for k in range(9):
            self.list_card[20+k] = repr(card_)
            self.list_i_change[20+k] = 2
            self.list_change[20+k] = th_MgCO3
            card_ = card_ + 1000
        # Fill lists w/ changes to MgCO3 masses for all 10 sections
        self.list_card[29] = '20507070'
        self.list_i_change[29] = 4
        self.list_change[29] = self.masses[0]
        card_ = 20514020
        for k in range(1,10):
            self.list_card[29+k] = repr(card_)
            self.list_i_change[29+k] = 4
            self.list_change[29+k] = self.masses[k]
            card_ = card_ + 10
    
    
    def final_ZION_calcs(self):
        ###################################
        """ Grab all the .csv file data """
        ###################################
        #--------------------------
        """ tempf.csv file data """
        #--------------------------
        # Read steam generator T-data from .csv file into 'tempf_data' array
        csv_name = self.csvFileLocation + "\\tempf_data.csv"
        tempf_data = np.zeros(2)
        with open(csv_name) as csvfile:
            numreader = csv.reader(csvfile)
            i = 0
            for row in numreader:
                tempf_data[i] = float(row[0]) # K
                i = i + 1
        # Assign tempf_data to the correct variables
        self.T_106 = tempf_data[0] - 273.15 # degC
        self.T_110 = tempf_data[1] - 273.15 # degC
        #----------------------
        """ p.csv file data """
        #----------------------
        # Read P-data from .csv file into 'P_106', 'P_110', and 'P_335[]'
        csv_name2 = self.csvFileLocation + "\\p_data.csv"
        p_data = np.zeros(9)
        with open(csv_name2) as csvfile2:
            numreader2 = csv.reader(csvfile2)
            i = 0
            for row in numreader2:
                p_data[i] = float(row[0]) # Pa
                i = i + 1
        # Assign p_data to the correct variables
        self.P_106 = p_data[0] / 10**5.0 # bar
        self.P_110 = p_data[1] / 10**5.0 # bar
        for i in range(6):
            self.P_335[i] = p_data[i+2] / 10**6.0 # MPa
        self.P_p_out = p_data[8] / 10**5.0 # bar
        #---------------------------
        """ mflowj.csv file data """
        #---------------------------
        # Read m_dot-data from .csv file into 1 combined array
        csv_name3 = self.csvFileLocation + "\\mflowj_data.csv"
        mflowj_data = np.zeros(29)
        with open(csv_name3) as csvfile3:
            numreader3 = csv.reader(csvfile3)
            i = 0
            for row in numreader3:
                mflowj_data[i] = float(row[0]) # kg/s
                i = i + 1
        # Assign averaged mflowj_data to the correct variables
        m_dot_100_data = mflowj_data[0:7] # kg/s
        self.m_dot_100 = np.average(m_dot_100_data) # kg/s
        m_dot_400_data = mflowj_data[7:14] # kg/s
        self.m_dot_400 = np.average(m_dot_400_data) # kg/s
        m_dot_600_data = mflowj_data[14:21] # kg/s
        self.m_dot_600 = np.average(m_dot_600_data) # kg/s
        m_dot_200_data = mflowj_data[21:28] # kg/s
        self.m_dot_200 = np.average(m_dot_200_data) # kg/s
        self.m_dot_335 = mflowj_data[28] # kg/s
        #--------------------------
        """ hvmix.csv file data """
        #--------------------------
        # Read enthalpy data from .csv file into 1 combined array
        csv_name4 = self.csvFileLocation + "\\hvmix_data.csv"
        hvmix_data = np.zeros(11)
        with open(csv_name4) as csvfile4:
            numreader4 = csv.reader(csvfile4)
            i = 0
            for row in numreader4:
                hvmix_data[i] = float(row[0]) # J/kg
                i = i + 1
        # Assign hvmix_data to the correct variables
        self.H_106 = hvmix_data[0] / 10**3.0 # kJ/kg
        self.H_110 = hvmix_data[1] / 10**3.0 # kJ/kg
        self.H_335_1 = hvmix_data[2] / 10**3.0 # kJ/kg
        self.H_112_5 = hvmix_data[3] / 10**3.0 # kJ/kg
        self.H_114 = hvmix_data[4] / 10**3.0 # kJ/kg
        self.H_412_5 = hvmix_data[5] / 10**3.0 # kJ/kg
        self.H_414 = hvmix_data[6] / 10**3.0 # kJ/kg
        self.H_612_5 = hvmix_data[7] / 10**3.0 # kJ/kg
        self.H_614 = hvmix_data[8] / 10**3.0 # kJ/kg
        self.H_212_5 = hvmix_data[9] / 10**3.0 # kJ/kg
        self.H_214 = hvmix_data[10] / 10**3.0 # kJ/kg
        #---------------------------
        """ httemp.csv file data """
        #---------------------------
        # Read fuel/cladding T-data from .csv file into 2 arrays
        csv_name5 = self.csvFileLocation + "\\httemp_data.csv"
        httemp_data = np.zeros(12)
        with open(csv_name5) as csvfile5:
            numreader5 = csv.reader(csvfile5)
            i = 0
            for row in numreader5:
                httemp_data[i] = float(row[0]) # K
                i = i + 1
        # Assign httemp_data to the correct variables
        for j in range(6):
            self.T_1336_1[j] = httemp_data[j] - 273.15 # degC
        for j in range(6):
            self.T_1336_17[j] = httemp_data[j+6] - 273.15 # degC
        #--------------------------
        """ quale.csv file data """
        #--------------------------
        # Read core quality data from .csv file into array
        csv_name6 = self.csvFileLocation + "\\quale_data.csv"
        quale_data = np.zeros(6)
        with open(csv_name6) as csvfile6:
            numreader6 = csv.reader(csvfile6)
            i = 0
            for row in numreader6:
                quale_data[i] = float(row[0])
                i = i + 1
        # Assign quale_data to the correct variables
        for j in range(6):
            self.x_e_335[j] = quale_data[j]
        #############################
        """ Thermal Design Limits """
        #############################
        #------------------------------------
        """ Enforce max fuel centerline T """
        #------------------------------------
        for T in self.T_1336_1:
            if T > self.T_fuel_cent_max:
                self.T_f_over_max = True
                self.penalized = True
        #---------------------------------
        """ Enforce max clad surface T """
        #---------------------------------
        for T in self.T_1336_17:
            if T > self.T_clad_surf_max:
                self.T_c_over_max = True
                self.penalized = True
        #-------------------------------
        """ Enforce MDNBR constraint """
        #-------------------------------
        # Calcs to find local critical heat flux per node
        A335 = self.A_335 / 10.76391 # m^2
        G = self.m_dot_335 / A335 # m^2 * s/kg
        Hf = np.zeros(6)
        for i in range(6):
            Pbar = self.P_335[i] * 10.0 # bar
            hf = h_Px(Pbar,0.0) # kJ/kg
            Hf[i] = hf # kJ/kg
        Hin = self.H_335_1 # kJ/kg
        D = self.Dh_335 * 3.28084 # m
        first = ((2.022-0.06238*self.P_335)+(0.1722-0.01427*self.P_335)*np.exp((18.177-0.5987*self.P_335)*self.x_e_335))
        second = ((0.1484-1.596*self.x_e_335+0.1729*self.x_e_335*abs(self.x_e_335))*2.326*G+3271)*(1.157-0.869*self.x_e_335)
        third = (0.2664+0.8357*np.exp(-124.1*D))*(0.8258+0.0003413*(Hf-Hin))
        qdp_cr = first * second * third # local critical heat flux per node [kW/m^2]
        # Calcs to find DNBR per node
        k_f = 0.00484 # kW/(m*K)
        k_g = 0.00153 # kW/(m*K)
        k_c = 0.01246 # kW/(m*K)
        Res_f = 1 / (4*np.pi*k_f) # (m*K)/kW
        Res_g = np.log(self.R_g/self.R_f) / (2.0*np.pi*k_g) # (m*K)/kW
        Res_c = np.log(self.R_c/self.R_g) / (2.0*np.pi*k_c) # (m*K)/kW
        qp = (self.T_1336_1 - self.T_1336_17) / (Res_f + Res_g + Res_c) # kW/m
        qdp = qp / (2.0*np.pi*self.R_c/3.28084) # kW/m^2
        DNBR = qdp_cr / qdp # DNBR per node
        # Minimum Departure from Nucleate Boiling Ratio (MDNBR)
        self.MDNBR = min(DNBR)
        # Check for need of penalty for MDNBR < 1.0
        if self.MDNBR < 1.0:
            self.MDNBR_below_1 = True
            self.penalized = True
        ##########################
        """ Find Q_PHX for PCS """
        ##########################
        # Assign the PHX's hot side stats here
        self.phx.Tin_hot = self.T_106 # degC
        self.phx.Tout_hot = self.T_110 # degC
        self.phx.P_hot = self.P_106 # bar
        # Make phx.Tout_cold = T_106 - 10 degC
        self.To_PHX = self.phx.Tin_hot - 10.0 # degC
        self.phx.Tout = self.To_PHX
        # If the new Tout will cause x_out to be 0, then lower P_c to P_sat
        if x_pT(self.phx.Pout,self.phx.Tout) == 0:
            self.phx.Pout = Psat_T(self.phx.Tout) - 1.0e-5 # bar
            self.phx.Pin = self.phx.Pout # bar
        # Then calculate Q_PHX
        DH = self.H_106 - self.H_110 # kJ/kg
        Q_PHX = self.m_dot_100 * DH # kW
        self.phx.Q_th = Q_PHX # kW
        ################################
        """ Objective Function Calcs """
        ################################
        #-------------------------
        """ Calc total W_pumps """
        #-------------------------
        W_p1 = self.m_dot_100 * (self.H_114 - self.H_112_5) / 10**3.0 # MW
        W_p2 = self.m_dot_400 * (self.H_414 - self.H_412_5) / 10**3.0 # MW
        W_p3 = self.m_dot_600 * (self.H_614 - self.H_612_5) / 10**3.0 # MW
        W_p4 = self.m_dot_200 * (self.H_214 - self.H_212_5) / 10**3.0 # MW
        W_pumps = W_p1 + W_p2 + W_p3 + W_p4 # MW
        if self.penalized == False:
            self.W_rcp = W_pumps # MW
        else:
            self.W_rcp = W_pumps * 10.0 # MW
        #----------------------
        """ Calc Pump Costs """
        #----------------------
        p1 = Pump_rcp(W_p1,self.P_p_out)
        p1.calc_Pump()
        p2 = Pump_rcp(W_p2,self.P_p_out)
        p2.calc_Pump()
        p3 = Pump_rcp(W_p3,self.P_p_out)
        p3.calc_Pump()
        p4 = Pump_rcp(W_p4,self.P_p_out)
        p4.calc_Pump()
        # Calc equipment cost of RCPs
        self.cost_RCPs = p1.cost + p2.cost + p3.cost + p4.cost # $
        # Calc operating cost of RCPs
        elec_rate = self.W_rcp * 1000 * 24 / p1.eff # kWh/day
        tot_elec = elec_rate * 365.25 * 40.0 # kWh/(40 yr)
        cost_of_elec = 0.12 # $0.12/kWh
        self.op_cost_RCPs = tot_elec * cost_of_elec # $/(40 yr)
        #-----------------------------------------
        """ Calc k_eff -> rho_0 -> nBc -> Cost """
        #-----------------------------------------
        N_A = 6.022*10**23.0 # molecule/mol
        # Define the needed molecular weights
        MW_Nb93 = 93.0 # g/mol
        MW_Sn_avg = 118.8077 # g/mol
        MW_Fe_avg = 55.9098 # g/mol
        MW_O2 = 32.0 # g/mol
        MW_Zr_avg = 91.3184 # g/mol
        MW_H2O = 18.02 # g/mol
        MW_235_UO2 = 267.0 # g/mol
        MW_238_UO2 = 270.0 # g/mol
        MW_He = 4.0 # g/mol
        # Define the wt%'s of Zircaloy
        Z_wt_Nb = 0.01
        Z_wt_Sn = 0.0095
        Z_wt_Fe = 0.0011
        Z_wt_O2 = 0.00125
        Z_wt_Zr = 0.97815
        # Define the sigma's for each isotope/element
        s_a_Zr_avg = 0.1886*10**-24.0 # cm^2
        s_a_Nb93 = 1.14236*10**-24.0 # cm^2
        s_a_Sn_avg = 0.6073085*10**-24.0 # cm^2
        s_a_Fe_avg = 2.5783*10**-24.0 # cm^2
        s_a_He_avg = 1.1*10**-34.0 # cm^2
        s_a_H1 = 0.332587*10**-24.0 # cm^2
        s_a_O16 = 0.0 # cm^2
        s_a_U238 = 2.6837*10**-24.0 # cm^2
        s_a_U235 = 686.0753*10**-24.0 # cm^2
        s_f_U235 = 586.691*10**-24.0 # cm^2
        # Define the needed densities
        rho_H2O = 1.0 # g/cm^3
        rho_Zircaloy = 6.55 # g/cm^3
        rho_He = 7.86*10**-3.0 # g/cm^3
        rho_UO2 = 10.97 # g/cm^3
        # Calc the number densities for each molecule
        N_H2O = rho_H2O*N_A/MW_H2O*(self.A_H2O/self.A_total) # molecule/cm^3
        N_Zr = Z_wt_Zr*rho_Zircaloy*N_A/MW_Zr_avg*(self.A_clad/self.A_total) # molecule/cm^3
        N_Nb = Z_wt_Nb*rho_Zircaloy*N_A/MW_Nb93*(self.A_clad/self.A_total) # molecule/cm^3
        N_Sn = Z_wt_Sn*rho_Zircaloy*N_A/MW_Sn_avg*(self.A_clad/self.A_total) # molecule/cm^3
        N_Fe = Z_wt_Fe*rho_Zircaloy*N_A/MW_Fe_avg*(self.A_clad/self.A_total) # molecule/cm^3
        N_O2 = Z_wt_O2*rho_Zircaloy*N_A/MW_O2*(self.A_clad/self.A_total) # molecule/cm^3
        N_He = rho_He*N_A/MW_He*(self.A_gap/self.A_total) # molecule/cm^3
        N_235_UO2 = 0.0495*rho_UO2*N_A/MW_235_UO2*(self.A_fuel/self.A_total) # molecule/cm^3
        N_238_UO2 = 0.9505*rho_UO2*N_A/MW_238_UO2*(self.A_fuel/self.A_total) # molecule/cm^3
        # Calc the Sigma's for each isotope/element and NF vs. F
        S_a_H1 = s_a_H1 * (2.0*N_H2O) # 1/cm
        S_a_O16 = s_a_O16 * (N_H2O + 2.0*N_O2 + 2.0*N_235_UO2 + 2.0*N_238_UO2) # 1/cm
        S_a_Zr = s_a_Zr_avg * (N_Zr) # 1/cm
        S_a_Nb93 = s_a_Nb93 * (N_Nb) # 1/cm
        S_a_Sn = s_a_Sn_avg * (N_Sn) # 1/cm
        S_a_Fe = s_a_Fe_avg * (N_Fe) # 1/cm
        S_a_He = s_a_He_avg * (N_He) # 1/cm
        S_a_U238 = s_a_U238 * (N_238_UO2) # 1/cm
        S_a_U235 = s_a_U235 * (N_235_UO2) # 1/cm
        S_a_NF = S_a_H1 + S_a_O16 + S_a_Zr + S_a_Nb93 + S_a_Sn + S_a_Fe + S_a_He + S_a_U238
        S_a_F = S_a_U235
        # Calc parts of the 6-factor formula
        ep = 1.02 * 0.8
        eta_U235 = (2.42 * s_f_U235) / (s_a_U235 + s_a_U238*(0.9505/0.0495))
        f = S_a_F / (S_a_F + S_a_NF)
        tau_H2O = 27.0 # cm^2
        L_sq_H2O = 8.1 # cm^2
        R_core = np.sqrt(self.A_total/np.pi)*30.48 # cm
        H_core = self.H_335*30.48 # cm
        B_sq = (2.405/R_core)**2.0 + (np.pi/H_core)**2.0 # cm^-2
        P_f_NL = np.exp(-B_sq*tau_H2O)
        P_th_NL = 1 / (1 + L_sq_H2O*B_sq)
        # Calculate k_eff (6-factor formula)
        self.k_eff = ep * eta_U235 * f * P_f_NL * P_th_NL
        # Calculate rho_0
        self.rho_0 = (self.k_eff - 1)/self.k_eff
        # Calculate Cycle-Burnup (Bc) and Discharge-Burnup (nBc)
        n = 3 # batches
        A1 = 8.333*10**-5.0 # 1/EFPD
        self.Bc = (n * self.rho_0) / (A1 * (n*(n+1))/2.0) # EFPD
        self.nBc = n * self.Bc / 365.2422 # yr (This is how long an entire core-full of fuel would last)
        # Calculate the total cost of fuel
        price_UO2 = 1787 # $/kg
        t_life = 40 # yr
        num_nBc = t_life / self.nBc
        m_fuel = self.n_pins * (np.pi*(self.R_f*30.48)**2.0) * H_core * rho_UO2/1000. # kg
        self.cost_total_fuel = m_fuel * price_UO2 * num_nBc # $
        #-----------------------
        """ Add up all costs """
        #-----------------------
        total_costs = self.cost_RCPs + self.op_cost_RCPs + self.cost_total_fuel # $
        if self.penalized == False:
            self.cost_1 = total_costs / (10**9.0) # $1x10^9
        else:
            self.cost_1 = total_costs / (10**9.0) * 10.0 # $1x10^9
    
    
    def final_PERCS_calcs(self):
        ###################################
        """ Grab all the .csv file data """
        ###################################
        #---------------------------
        """ T_data.csv file data """
        #---------------------------
        # Open the .csv file with core outlet T-data
        csv_name7 = self.csvFileLocation2 + "\\T_data.csv"
        first = True
        T_0 = 0.0
        with open(csv_name7) as csvfile7:
            numreader7 = csv.reader(csvfile7)
            for row in numreader7:
                # Grab and assign the first T-value
                if first == True:
                    T_0 = 600 # K
                    first = False
                # Check to see if T ever goes above 620 K
                if float(row[1]) > 620.0:
                    self.T_over_620 = True
                # Check to see if T ever goes above 635 K
                if float(row[1]) > 635.0:
                    self.T_over_635 = True
                # Fill arrays with data
                self.t = np.append(self.t,float(row[0])) # sec
                self.T_335_6 = np.append(self.T_335_6,float(row[1])) # K
                self.dT_335_6 = np.append(self.dT_335_6,float(row[1])-T_0) # delta_K
        ##########################################
        """ Check for need of a penalty factor """
        ##########################################
        # Create a value for the penalty factor
        deficit = t_final - self.t[-1]
        if abs(deficit) < 1.0:
            self.PERCS_penalty = 1.0
        else:
            self.PERCS_failed = True
            self.PERCS_penalty = np.exp(3.0*deficit/t_final)
        # Check for a failure during the steady-state RELAP run (these need to penalized heavily)
        t_step = np.round(self.t[1]-self.t[0])
        if t_step == 40.0:
            self.ss_fail = True
            self.PERCS_penalty = ss_fail_penalty * self.PERCS_penalty
        ################################
        """ Objective Function Calcs """
        ################################
        #----------------------------
        """ Calc total PERCS cost """
        #----------------------------
        # Add up all the equipment costs
        pipe_costs = (self.p716.cost + self.p717.cost) * cost_multiplier_for_nucl_safety_grade
        PERCS_costs = (self.hx.cost + self.tank.cost) * cost_multiplier_for_nucl_safety_grade
        other_costs = self.support.cost + self.chemical.cost
        tot_cost = pipe_costs + PERCS_costs + other_costs
        # Figure out the cost penalty
        if self.PERCS_failed == False:
            self.cost_penalty = 1.0
        else:
            self.cost_penalty = 75.0
        # Assign the total cost to the Option
        self.cost_2 = (tot_cost * self.cost_penalty)*10**-8.0 # $1x10^8
        #------------------
        """ Calc dT_int """
        #------------------
        # Find area under the curve using Trapezoid Method
        sum_ = 0.0
        for j in range(len(self.T_335_6)-1):
            # Only add to the total if both points are (+)
            # The obj. function will take only the (+) area under the dT curve
            if self.dT_335_6[j]>0 and self.dT_335_6[j+1]>0:
                area = 0.5*(self.dT_335_6[j]+self.dT_335_6[j+1]) * (self.t[j+1]-self.t[j])
                sum_ = sum_ + area
        integral_dT = sum_ # K*s
        # Assign the total integral of dT to the Option
        self.dT_int = integral_dT * self.PERCS_penalty
        # If T_335_6 never got up to 600 K before failing early, then dT_int = 0, which escapes a penalty
        # Fix it by forcing dT_int to be large
        if self.dT_int < 10000.0:
            self.dT_int = 500000.0
    
    def Alpha_calcs(self):
        #--------------------------
        """ alpha.csv file data """
        #--------------------------
        # Read the .csv file with alpha data into a matrix
        csv_name8 = self.csvFileLocation2 + "\\Alpha_data.csv"
        a_row = 0
        a_col = 0
        with open(csv_name8) as csvfile8:
            time.sleep(5)
            numreader8 = csv.reader(csvfile8)
            for row in numreader8:
                self.a[a_row,a_col] = float(row[0])
                a_row = a_row + 1
                if a_row % 10 == 0:
                    a_row = 0
                    a_col = a_col + 1
        #-------------------------
        """ Calc overall Alpha """
        #-------------------------
        sum_2 = 0.0
        # Rotate through col's of a-data
        for m in range(10):
            sum_2 = sum_2 + self.ratio[m]*np.average(self.a[:,m])
        Alpha = sum_2
        # Assign the overall average Alpha to the Option
        penalty_adjustment = 1.0
        if self.PERCS_failed == True and self.ss_fail == True: penalty_adjustment = 100
        self.alpha = Alpha * (self.PERCS_penalty/penalty_adjustment)
    
    def PCS_SS_calcs(self):
        # Redudantly check to make sure Alpha_calcs() worked
        if self.alpha == 0.0:
            self.Alpha_calcs()
        ##########################
        """ Converge the Cycle """
        ##########################
        # Establish the test mdot
        mdot_test = 600. # kg/s
        # Assign test mdot and solve for Tin
        self.phx.mdot = mdot_test
        self.solve_PCS()
        Tin = self.s[36].T
        # Calculate the real mdot
        Hin = h_pT(self.phx.Pin,Tin)
        Hout = h_pT(self.phx.Pout,self.phx.Tout)
        mdot_real = self.phx.Q_th / (Hout - Hin)
        # Assign real mdot and solve Option
        self.phx.mdot = mdot_real
        self.solve_PCS()
        #################################################
        """ Assign Superstructure Configuration Color """
        #################################################
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
        # Assign color for specific y-value scheme
        self.color = c
        ################################
        """ Perform obj. func. calcs """
        ################################
        #-------------------
        """ Calc PCS eff """
        #-------------------
        W_t1 = self.t1a.W+self.t1b.W+self.t1c.W # MW
        W_turb = W_t1 + self.t2a.W+self.t2b.W+self.t3.W+self.t4.W+self.t5.W # MW
        W_pump = self.p1.W+self.p2.W+self.p3.W+self.p4.W+self.p5.W # MW
        self.eff = (W_turb - W_pump) / (self.phx.Q_th/1000) # frac
        self.inv_eff = 1.0 - self.eff # frac
        #--------------------
        """ Calc PCS cost """
        #--------------------
        cost_phx = self.phx.cost
        cost_ms = self.ms1.cost + self.ms2.cost
        cost_rh = self.rh1.cost + self.rh2.cost
        cost_turb = self.t1.cost+self.t2.cost+self.t3.cost+self.t4.cost+self.t5.cost
        cost_cond = self.cond.cost
        cost_pump = self.p1.cost+self.p2.cost+self.p3.cost+self.p4.cost+self.p5.cost
        cost_fwh = self.fwh1.cost+self.fwh2.cost+self.fwh3.cost+self.fwh4.cost
        total_cost = cost_phx+cost_ms+cost_rh+cost_turb+cost_cond+cost_pump+cost_fwh
        self.cost_3 = (total_cost * 10.0**-9.0) * cost_multiplier_for_nucl_safety_grade # $1x10^9
    
    # Calculate all streams and equipment from PHX-out to PHX-in
    def solve_PCS(self):
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
            # Add to T for the sake of h_pT(), since this stream skips RH 1
            self.s[9].T = self.ms1.T + (1e-10)
            self.s[9].mdot = self.ms1.mdot_V
            self.s[9].x = 1.0
        else:
            self.s[9].y = 0
        """ Stream 10 """
        if self.y_ipt==0 and self.y_rh1==0:
            self.s[10].y = 1
            self.s[10].P = self.ms1.P
            # Add to T for the sake of h_pT(), since this stream skips RH 1
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
        # If there was a pinch in RH 1
        if self.rh1.pinch == True:
            # then mark the PCS's pinch_point as True
            self.pinch_point = True
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
            # Add to T for the sake of h_pT(), since this stream skips RH 2
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
        # If there was a pinch in RH 2
        if self.rh2.pinch == True:
            # then mark the PCS's pinch_point as True
            self.pinch_point = True
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
        



