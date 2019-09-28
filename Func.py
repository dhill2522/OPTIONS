"""
Created on Mon Nov 05 03:52:36 2018
@author: Paul
"""

from Class2_Pop import *
from xlwings import *
from iapws97 import IAPWS97


###############################################################################
"""""""""  FUNCTIONS  """"""""" ###############################################
###############################################################################

"""""""""   Steam Table Functions   """""""""

def Psat_T(T_):
    T_ = T_ + 273.15 # Convert from C to K
    example = IAPWS97(T=T_,x=1.0)
    P_ = example.P * 10.0 # Convert from MPa to bar
    return P_

def h_Px(P_,x_):
    P_ = P_/10.0 # Convert from bar to MPa
    example = IAPWS97(P=P_,x=x_)
    h_ = example.h # kJ/kg
    return h_

def h_pT(P_,T_):
    P_ = P_/10.0 # Convert from bar to MPa
    T_ = T_ + 273.15 # Convert from C to K
    example = IAPWS97(P=P_,T=T_)
    h_ = example.h # kJ/kg
    return h_

def h_pS(P_,S_):
    P_ = P_/10.0 # Convert from bar to MPa
    S_ = S_ # kJ/kg*K
    example = IAPWS97(P=P_,s=S_)
    h_ = example.h # kJ/kg
    return h_

def h_Tx(T_,x_):
    T_ = T_ + 273.15 # Convert from C to K
    example = IAPWS97(T=T_,x=x_)
    h_ = example.h # kJ/kg
    return h_

def S_pT(P_,T_):
    P_ = P_/10.0 # Convert from bar to MPa
    T_ = T_ + 273.15 # Convert from C to K
    example = IAPWS97(P=P_,T=T_)
    S_ = example.s # kJ/kg*K
    return S_

def S_ph(P_,h_):
    P_ = P_/10.0 # Convert from bar to MPa
    h_ = h_ # kJ/kg
    example = IAPWS97(P=P_,h=h_)
    S_ = example.s # kJ/kg*K
    return S_

def T_ph(P_,h_):
    P_ = P_/10.0 # Convert from bar to MPa
    h_ = h_ # kJ/kg
    example = IAPWS97(P=P_,h=h_)
    T_ = example.T - 273.15 # Convert from K to C
    return T_

def x_pT(P_,T_):
    P_ = P_/10.0 # Convert from bar to MPa
    T_ = T_ + 273.15 # Convert from C to K
    example = IAPWS97(P=P_,T=T_)
    x_ = example.x
    return x_

def x_ph(P_,h_):
    P_ = P_/10.0 # Convert from bar to MPa
    h_ = h_ # kJ/kg
    example = IAPWS97(P=P_,h=h_)
    x_ = example.x
    return x_

def rhoV_P(P_):
    P_ = P_/10.0 # Convert from bar to MPa
    x_ = 1.0 # Vapor
    example = IAPWS97(P=P_,x=x_)
    rho_ = example.rho
    return rho_

def rhoL_P(P_):
    P_ = P_/10.0 # Convert from bar to MPa
    x_ = 0.0 # Liquid
    example = IAPWS97(P=P_,x=x_)
    rho_ = example.rho
    return rho_

###########################################
"""""""""   Enforce Constraints   """""""""
###########################################

def constraints(x1,x2,x3,y):
    """
    Inputs:
        x1 = Array of core loop x-optimization parameters
        x2 = Array of PERCS x-optimization parameters
        x3 = Array of PCS x-optimization parameters
        y = Array of PCS y-optimization parameters
    Output:
        x1 = Array of constraint-corrected core loop x-optimization parameters
        x2 = Array of constraint-corrected PERCS x-optimization parameters
        x3 = Array of constraint-corrected PCS x-optimization parameters
        y = Array of unchanged PCS y-optimization parameters
    """
    #################################
    """ Constraints for x1 (Core) """
    #################################
    # Constraints for R_fuel
    if x1[0] < 0.0120: x1[0] = 0.0120
    if x1[0] > 0.0159: x1[0] = 0.0159
    # Constraints for H_fuel
    if x1[1] < 10.0: x1[1] = 10.0
    if x1[1] > 14.0: x1[1] = 14.0
    # Constraints for Dh_00
    if x1[2] < 1.5: x1[2] = 1.5
    if x1[2] > 3.5: x1[2] = 3.5
    # Constraints for Dh_12
    if x1[3] < 1.5: x1[3] = 1.5
    if x1[3] > 3.5: x1[3] = 3.5
    # Constraints for Dh_14
    if x1[4] < 1.5: x1[4] = 1.5
    if x1[4] > 3.5: x1[4] = 3.5
    ##################################
    """ Constraints for x2 (PERCS) """
    ##################################
    # Constraints for R_tank
    if x2[0] < 3.28084: x2[0] = 3.28084
    if x2[0] > 16.4042: x2[0] = 16.4042
    # Constraints for D_h
    if x2[2] < 3.28084e-2: x2[2] = 3.28084e-2
    if x2[2] > 1.64042e-1: x2[2] = 1.64042e-1
    # Constraints for pitch
    if x2[1] < 1.25*x2[2]: x2[1] = 1.25*x2[2]
    if x2[1] > 0.65617: x2[1] = 0.65617
    # Constraints for th
    if x2[3] < 3.28084e-3: x2[3] = 3.28084e-3
    if x2[3] > 8.2021e-3: x2[3] = 8.2021e-3
    # Constraints for Len
    if x2[4] < 3.28084: x2[4] = 3.28084
    if x2[4] > 49.2126: x2[4] = 49.2126
    # Constraints for elev
    if x2[5] < 16.4042: x2[5] = 16.4042
    if x2[5] < (x2[4]-15.6): x2[5] = x2[4]-15.6
    if x2[5] > 131.234: x2[5] = 131.234
    ################################
    """ Constraints for x3 (PCS) """
    ################################
    #------------------------
    """ Bound Constraints """
    #------------------------
    # Constraints for mchx.Tout
    if x3[0] < 295.0: x3[0] = 295.0
    if x3[0] > 307.7: x3[0] = 307.7
    # Constraints for Po_t1a
    if x3[1] < 8.0: x3[1] = 8.0
    if x3[1] > 25.0: x3[1] = 25.0
    # Constraints for mf_t1a
    if x3[2] < 0.05: x3[2] = 0.05
    if x3[2] > 0.19: x3[2] = 0.19
    # Constraints for Po_t1b
    if x3[3] < 6.0: x3[3] = 6.0
    if x3[3] > 20.0: x3[3] = 20.0
    # Constraints for mf_t1b
    if x3[4] < 0.05: x3[4] = 0.05
    if x3[4] > 0.19: x3[4] = 0.19
    # Constraints for Po_t1c
    if x3[5] < 4.0: x3[5] = 4.0
    if x3[5] > 16.0: x3[5] = 16.0
    # Constraints for Po_t2a
    if x3[6] < 3.0: x3[6] = 3.0
    if x3[6] > 13.0: x3[6] = 13.0
    # Constraints for mf_t2a
    if x3[7] < 0.05: x3[7] = 0.05
    if x3[7] > 0.19: x3[7] = 0.19
    # Constraints for Po_t2b
    if x3[8] < 2.0: x3[8] = 2.0
    if x3[8] > 11.0: x3[8] = 11.0
    #----------------------------
    """ Coupled P Constraints """
    #----------------------------
    # Constraints for Turbine Combo: t1a(o),t1b(x),t1c(o),t2a(x),t2b(x)
    ''' 1 - red '''
    if y[0]==0 and y[1]==1 and y[2]==0 and y[3]==0 and y[4]==0 and y[5]==0:
        x3[1],x3[5] = DeltaP(x3[1],x3[5],4.0)
        if x3[3]>=x3[1] or x3[3]<=x3[5]: x3[3] = x3[5] + 1.0
    # Constraints for Turbine Combo: t1a(x),t1b(x),t1c(o),t2a(x),t2b(x)
    ''' 2 - firebrick '''
    if y[0]==0 and y[1]==0 and y[2]==0 and y[3]==0 and y[4]==0 and y[5]==0:
        if x3[3]<=x3[5]: x3[3] = x3[5] + 2.0
        if x3[1]<=x3[3]: x3[1] = x3[3] + 2.0
        if x3[6]>=x3[5]: x3[6] = x3[5] - 2.0
        if x3[8]>=x3[6]: x3[8] = x3[6] - 2.0
    # Constraints for Turbine Combo: t1a(o),t1b(x),t1c(o),t2a(x),t2b(o)
    ''' 3 - darkgreen '''
    if y[0]==1 and y[1]==1 and y[2]==1 and y[3]==0 and y[4]==1 and y[5]==0:
        x_t = np.array((2.0,x3[8],4.0,x3[5],4.0,x3[1]))
        x_t = DeltaP4(x_t)
        x3[8],x3[5],x3[1] = x_t[1],x_t[3],x_t[5]
        if x3[3]>=x3[1] or x3[3]<=x3[5]: x3[3] = x3[5] + 1.0
        if x3[6]>=x3[5] or x3[6]<=x3[8]: x3[6] = x3[8] + 1.0
    ''' 9 - yellow '''
    if y[0]==1 and y[1]==1 and y[2]==0 and y[3]==0 and y[4]==0 and y[5]==0:
        x_t = np.array((2.0,x3[8],4.0,x3[5],4.0,x3[1]))
        x_t = DeltaP4(x_t)
        x3[8],x3[5],x3[1] = x_t[1],x_t[3],x_t[5]
        if x3[3]>=x3[1] or x3[3]<=x3[5]: x3[3] = x3[5] + 1.0
        if x3[6]>=x3[5] or x3[6]<=x3[8]: x3[6] = x3[8] + 1.0
    # Constraints for Turbine Combo: t1a(o),t1b(o),t1c(o),t2a(x),t2b(o)
    ''' 4 - purple '''
    if y[0]==1 and y[1]==1 and y[2]==1 and y[3]==0 and y[4]==0 and y[5]==1:
        x_t = np.array((2.0,x3[8],2.0,x3[5],2.0,x3[3],2.0,x3[1]))
        x_t = DeltaP4(x_t)
        x3[8],x3[5],x3[3],x3[1] = x_t[1],x_t[3],x_t[5],x_t[7]
        # Make sure the dummy value for x[6] is b/t Po_t1c and Po_t2b
        if x3[6]>=x3[5] or x3[6]<=x3[8]: x3[6] = x3[8] + 1.0
    # Constraints for Turbine Combo: t1a(x),t1b(o),t1c(o),t2a(x),t2b(o)
    ''' 5 - deeppink '''
    if y[0]==1 and y[1]==0 and y[2]==1 and y[3]==0 and y[4]==0 and y[5]==1:
        x_t = np.array((2.0,x3[8],4.0,x3[5],2.0,x3[3]))
        x_t = DeltaP4(x_t)
        x3[8],x3[5],x3[3] = x_t[1],x_t[3],x_t[5]
        if x3[1]<=x3[3]: x3[1] = x3[3] + 1.0
        if x3[6]>=x3[5] or x3[6]<=x3[8]: x3[6] = x3[8] + 1.0
    # Constraints for Turbine Combo: t1a(o),t1b(x),t1c(o),t2a(o),t2b(o)
    ''' 6 - blue '''
    if y[0]==1 and y[1]==1 and y[2]==1 and y[3]==1 and y[4]==0 and y[5]==0:
        x_t = np.array((x3[1],4.0,x3[5],2.0,x3[6],4.0,x3[8]))
        x_t = DeltaP3(x_t)
        x3[1],x3[5],x3[6],x3[8] = x_t[0],x_t[2],x_t[4],x_t[6]
        if x3[3]>=x3[1] or x3[3]<=x3[5]: x3[3] = x3[5] + 1.0
    # Constraints for Turbine Combo: t1a(x),t1b(x),t1c(o),t2a(o),t2b(o)
    ''' 7 - cyan '''
    if y[0]==1 and y[1]==0 and y[2]==1 and y[3]==1 and y[4]==0 and y[5]==0:
        x_t = np.array((x3[5],2.0,x3[6],4.0,x3[8]))
        x_t = DeltaP3(x_t)
        x3[5],x3[6],x3[8] = x_t[0],x_t[2],x_t[4]
        if x3[3]<=x3[5]: x3[3] = x3[5] + 1.0
        if x3[1]<=x3[3]: x3[1] = x3[3] + 1.0
    # Constraints for Turbine Combo: t1a(x),t1b(x),t1c(o),t2a(x),t2b(o)
    ''' 8 - orange '''
    if y[0]==1 and y[1]==0 and y[2]==0 and y[3]==0 and y[4]==0 and y[5]==0:
        x3[5],x3[8] = DeltaP(x3[5],x3[8],4.0)
        if x3[3]<=x3[5]: x3[3] = x3[5] + 2.0
        if x3[1]<=x3[3]: x3[1] = x3[3] + 2.0
        if x3[6]>=x3[5] or x3[6]<=x3[8]: x3[6] = x3[8] + 1.0
    # To correct any "NaN" problems with Pump6, make sure that IPT.Pout > 1.0 bar
    if x3[8]<1.0: x3[8] = 1.001
    # Return correct x,y-values
    return x1,x2,x3,y

def DeltaP(a,b,delP):
    """
    Inputs:
        a = Upstream pressure
        b = Downstream pressure
        delP = Minimum required pressure gap between a and b
    """
    # If pressures a and b are not at least delP apart
    if (a-b) < delP:
        # Randomly increase a or decrease b
        if random.random() < 0.5:
            a = b + delP
        else:
            b = a - delP
    # Return gap-corrected pressures
    return a,b

def DeltaP3(a):
    """
    Inputs:
        a[0] = Highest P
        a[1] = Min delta_P b/t a[0] and a[2]
        a[2] = Next P
        a[3] = Min delta_P b/t a[2] and a[4]
        a[4] = Next P
        a[...] = ...
    """
    # Rotate through the P's
    for n in range(0,len(a)-2,2):
        # If a[0] and a[2] are not at least a[1] apart
        if n==0 and (a[n]-a[n+2])<a[n+1]:
            # Randomly increase a[0] or decrease a[2]
            if random.random()<0.5:
                a[n] = a[n+2] + a[n+1]
            else:
                a[n+2] = a[n] - a[n+1]
        if n != 0:
            # If a[n] and a[n+2] are not at least a[n+1] apart,
            # and there is enough gap b/t a[n-2] and a[n]
            if (a[n]-a[n+2])<a[n+1] and (a[n-2]-a[n+2])>=(a[n-1]+a[n+1]):
                # Randomly increase a[n] or decrease a[n+2]
                if random.random()<0.5:
                    a[n] = a[n+2] + a[n+1]
                else:
                    a[n+2] = a[n] - a[n+1]
            # If a[n] and a[n+2] are not at least a[n+1] apart,
            # and there is not enough gap b/t a[n-2] and a[n]
            if (a[n]-a[n+2])<a[n+1] and (a[n-2]-a[n+2])<(a[n-1]+a[n+1]):
                # Just decrease a[n+2]
                a[n+2] = a[n] - a[n+1]
    # Return gap-corrected pressures
    return a

def DeltaP4(a):
    """
    Inputs:
        a[0] = Min P for a[1]
        a[1] = Lowest P
        a[2] = Mim delta_P b/t a[1] and a[3]
        a[3] = Next P
        a[4] = Min delta_P b/t a[3] and a[5]
        a[5] = Next P
        a[...] = ...
    """
    # Rotate through the P's
    for n in range(1,len(a)-2,2):
        # For the 1st P, when 1st and 2nd P's are not at least dP apart
        if n == 1 and (a[n+2]-a[n])<a[n+1]:
            # If (2nd P minus dP) < min for 1st P
            if (a[n+2]-a[n+1])<a[0]:
                # Enforce dP on 2nd P
                a[n+2] = a[n] + a[n+1]
            # If not so
            else:
                # Make it a 50/50 chance
                if random.random()<0.5:
                    a[n] = a[n+2] - a[n+1]
                else:
                    a[n+2] = a[n] + a[n+1]
        # For all other P's, when P and P_next are not at least dP_next apart
        if n != 1 and (a[n+2]-a[n])<a[n+1]:
            # If (P_next minus dP_next) < (P_prev + dP_prev)
            if (a[n+2]-a[n+1])<(a[n-2]+a[n-1]):
                # Enforce dP_next on P_next
                a[n+2] = a[n] + a[n+1]
            # If not so
            else:
                # Make it a 50/50 chance
                if random.random()<0.5:
                    a[n] = a[n+2] - a[n+1]
                else:
                    a[n+2] = a[n] + a[n+1]
    # Return gap-corrected pressures
    return a

    
########################################
"""""""""   Maximin Function   """""""""
########################################

def maximin(opt_list):
    """
    Inputs:
        opt_list = Population's 'opt_list' array of Option's
    Outputs:
        sorted_list = Array of fitness-sorted Option's
    """
    n = len(opt_list)
    ###############################################
    """ Normalize the objective function values """
    ###############################################
    #--------------------------
    """ Core Obj. Functions """
    #--------------------------
    # Find the largest and smallest W_rcp among the list of Options
    largest_W_rcp = 0.0
    smallest_W_rcp = 99.9e25
    for i in range(n):
        if opt_list[i].W_rcp > largest_W_rcp:
            largest_W_rcp = opt_list[i].W_rcp
        if opt_list[i].W_rcp < smallest_W_rcp:
            smallest_W_rcp = opt_list[i].W_rcp
    # Assign normalized W_rcp-values to each Option
    for i in range(n):
        opt_list[i].obj_1_1 = (opt_list[i].W_rcp-smallest_W_rcp)/(largest_W_rcp-smallest_W_rcp)
    # Find the largest and smallest Cost_1 among the list of Options
    largest_cost_1 = 0.0
    smallest_cost_1 = 99.9e25
    for i in range(n):
        if opt_list[i].cost_1 > largest_cost_1:
            largest_cost_1 = opt_list[i].cost_1
        if opt_list[i].cost_1 < smallest_cost_1:
            smallest_cost_1 = opt_list[i].cost_1
    # Assign normalized Cost_1-values to each Option
    for i in range(n):
        opt_list[i].obj_1_2 = (opt_list[i].cost_1-smallest_cost_1)/(largest_cost_1-smallest_cost_1)
    #---------------------------
    """ PERCS Obj. Functions """
    #---------------------------
    # Find the largest and smallest Cost_2 among the list of Options
    largest_cost_2 = 0.0
    smallest_cost_2 = 99.9e25
    for i in range(n):
        if opt_list[i].cost_2 > largest_cost_2:
            largest_cost_2 = opt_list[i].cost_2
        if opt_list[i].cost_2 < smallest_cost_2:
            smallest_cost_2 = opt_list[i].cost_2
    # Assign normalized Cost_2-values to each Option
    for i in range(n):
        opt_list[i].obj_2_1 = (opt_list[i].cost_2-smallest_cost_2)/(largest_cost_2-smallest_cost_2)
    # Find the largest and smallest T_dev among the list of Options
    largest_dT_int = 0.0
    smallest_dT_int = 99.9e25
    for i in range(n):
        if opt_list[i].dT_int > largest_dT_int:
            largest_dT_int = opt_list[i].dT_int
        if opt_list[i].dT_int < smallest_dT_int:
            smallest_dT_int = opt_list[i].dT_int
    # Assign normalized T_dev-values to each Option
    for i in range(n):
        opt_list[i].obj_2_2 = (opt_list[i].dT_int-smallest_dT_int)/(largest_dT_int-smallest_dT_int)
    # Find the largest and smallest Alpha among the list of Options
    largest_alpha = 0.0
    smallest_alpha = 99.9e25
    for i in range(n):
        if opt_list[i].alpha > largest_alpha:
            largest_alpha = opt_list[i].alpha
        if opt_list[i].alpha < smallest_alpha:
            smallest_alpha = opt_list[i].alpha
    # Assign normalized Alpha-values to each Option
    for i in range(n):
        opt_list[i].obj_2_3 = (opt_list[i].alpha-smallest_alpha)/(largest_alpha-smallest_alpha)
    #----------------------------
    """ PCS SS Obj. Functions """
    #----------------------------
    # Find the largest and smallest inv_eff among the list of Options
    largest_inv_eff = 0.0
    smallest_inv_eff = 99.9e25
    for i in range(n):
        if opt_list[i].inv_eff > largest_inv_eff:
            largest_inv_eff = opt_list[i].inv_eff
        if opt_list[i].inv_eff < smallest_inv_eff:
            smallest_inv_eff = opt_list[i].inv_eff
    # Assign normalized inv_eff-values to each Option
    for i in range(n):
        opt_list[i].obj_3_1 = (opt_list[i].inv_eff-smallest_inv_eff)/(largest_inv_eff-smallest_inv_eff)
    # Find the largest and smallest Cost_3 among the list of Options
    largest_cost_3 = 0.0
    smallest_cost_3 = 99.9e25
    for i in range(n):
        if opt_list[i].cost_3 > largest_cost_3:
            largest_cost_3 = opt_list[i].cost_3
        if opt_list[i].cost_3 < smallest_cost_3:
            smallest_cost_3 = opt_list[i].cost_3
    # Assign normalized Cost_3-values to each Option
    for i in range(n):
        opt_list[i].obj_3_2 = (opt_list[i].cost_3-smallest_cost_3)/(largest_cost_3-smallest_cost_3)
    ######################################################
    """ Calc the individual fmm values for each Option """
    ######################################################
    #-----------------------------
    """ Core fmm value (fmm_1) """
    #-----------------------------
    # Rotate through each Option
    for i in range(n):
        j_mins = np.empty(0)
        # Rotate through each Option (except j = i)
        for j in range(n):
            if j == i:
                None
            else:
                # Find min[k=1,2](f_k(x_i) - f_k(x_j))
                k_diff1 = opt_list[i].obj_1_1 - opt_list[j].obj_1_1
                k_diff2 = opt_list[i].obj_1_2 - opt_list[j].obj_1_2
                k_min = min(k_diff1,k_diff2)
                j_mins = np.append(j_mins,k_min)
        # Find the max of the j_mins and assign new f_maximin to Option i
        i_max = max(j_mins)
        opt_list[i].fmm_1 = i_max
    #------------------------------
    """ PERCS fmm value (fmm_2) """
    #------------------------------
    # Rotate through each Option
    for i in range(n):
        j_mins = np.empty(0)
        # Rotate through each Option (except j = i)
        for j in range(n):
            if j == i:
                None
            else:
                # Find min[k=1,2](f_k(x_i) - f_k(x_j))
                k_diff1 = opt_list[i].obj_2_1 - opt_list[j].obj_2_1
                k_diff2 = opt_list[i].obj_2_2 - opt_list[j].obj_2_2
                k_diff3 = opt_list[i].obj_2_3 - opt_list[j].obj_2_3
                k_min = min(k_diff1,k_diff2,k_diff3)
                j_mins = np.append(j_mins,k_min)
        # Find the max of the j_mins and assign new f_maximin to Option i
        i_max = max(j_mins)
        opt_list[i].fmm_2 = i_max
    #----------------------------
    """ PCS fmm value (fmm_3) """
    #----------------------------
    # Rotate through each Option
    for i in range(n):
        j_mins = np.empty(0)
        # Rotate through each Option (except j = i)
        for j in range(n):
            if j == i:
                None
            else:
                # Find min[k=1,2](f_k(x_i) - f_k(x_j))
                k_diff1 = opt_list[i].obj_3_1 - opt_list[j].obj_3_1
                k_diff2 = opt_list[i].obj_3_2 - opt_list[j].obj_3_2
                k_min = min(k_diff1,k_diff2)
                j_mins = np.append(j_mins,k_min)
        # Find the max of the j_mins and assign new f_maximin to Option i
        i_max = max(j_mins)
        opt_list[i].fmm_3 = i_max
    ###################################################
    """ Normalize the fmm objective function values """
    ###################################################
    #------------------
    """ fmm_1-value """
    #------------------
    # Find the largest and smallest fmm_1 among the list of Options
    largest_fmm_1 = -99.9e25
    smallest_fmm_1 = 99.9e25
    for i in range(n):
        if opt_list[i].fmm_1 > largest_fmm_1:
            largest_fmm_1 = opt_list[i].fmm_1
        if opt_list[i].fmm_1 < smallest_fmm_1:
            smallest_fmm_1 = opt_list[i].fmm_1
    # Assign normalized fmm_1-values to each Option
    for i in range(n):
        opt_list[i].obj_fmm_1 = (opt_list[i].fmm_1-smallest_fmm_1)/(largest_fmm_1-smallest_fmm_1)
    #------------------
    """ fmm_2-value """
    #------------------
    # Find the largest and smallest fmm_2 among the list of Options
    largest_fmm_2 = -99.9e25
    smallest_fmm_2 = 99.9e25
    for i in range(n):
        if opt_list[i].fmm_2 > largest_fmm_2:
            largest_fmm_2 = opt_list[i].fmm_2
        if opt_list[i].fmm_2 < smallest_fmm_2:
            smallest_fmm_2 = opt_list[i].fmm_2
    # Assign normalized fmm_2-values to each Option
    for i in range(n):
        opt_list[i].obj_fmm_2 = (opt_list[i].fmm_2-smallest_fmm_2)/(largest_fmm_2-smallest_fmm_2)
    #------------------
    """ fmm_3-value """
    #------------------
    # Find the largest and smallest fmm_3 among the list of Options
    largest_fmm_3 = -99.9e25
    smallest_fmm_3 = 99.9e25
    for i in range(n):
        if opt_list[i].fmm_3 > largest_fmm_3:
            largest_fmm_3 = opt_list[i].fmm_3
        if opt_list[i].fmm_3 < smallest_fmm_3:
            smallest_fmm_3 = opt_list[i].fmm_3
    # Assign normalized fmm_3-values to each Option
    for i in range(n):
        opt_list[i].obj_fmm_3 = (opt_list[i].fmm_3-smallest_fmm_3)/(largest_fmm_3-smallest_fmm_3)
    ###########################################################
    """ Calc the overall fmm values (fmm_o) for each Option """
    ###########################################################
    fmm_min = 1000
    min_index = 50
    # Rotate through each Option
    for i in range(n):
        j_mins = np.empty(0)
        # Rotate through each Option (except j = i)
        for j in range(n):
            if j == i:
                None
            else:
                # Find min[k=1,2](f_k(x_i) - f_k(x_j))
                k_diff1 = opt_list[i].obj_fmm_1 - opt_list[j].obj_fmm_1
                k_diff2 = opt_list[i].obj_fmm_2 - opt_list[j].obj_fmm_2
                k_diff3 = opt_list[i].obj_fmm_3 - opt_list[j].obj_fmm_3
                k_min = min(k_diff1,k_diff2,k_diff3)
                j_mins = np.append(j_mins,k_min)
        # Find the max of the j_mins and assign new f_maximin to Option i
        i_max = max(j_mins)
        opt_list[i].fmm_o = i_max
        # Keep track of the smallest f_maximin
        if i_max < fmm_min:
            fmm_min = i_max
            min_index = i
    ###################################################################
    """ Sort the list of Options by order of ascending fmm_o-values """
    ###################################################################
    # Initialize the maximin-sorted list of Options
    sorted_list = np.array(opt_list[min_index])
    # Re-order the list of Options in ascending "maximin" order
    for count in range(n-1):
        fmm_next = 1000
        # Rotate through the Options of opt_list
        for i in range(n):
            # Find the next smallest f_maximin
            if i != min_index:
                # If current Option's fmm_value is less than fmm_next
                if opt_list[i].fmm_o < fmm_next and opt_list[i].fmm_o >= fmm_min:
                    # If it equals the previous minimum
                    if opt_list[i].fmm_o == fmm_min and i > min_index:
                        index_next = i
                        fmm_next = opt_list[i].fmm_o
                        break
                    else:
                        if opt_list[i].fmm_o == fmm_min and i < min_index:
                            None
                        else:
                            index_next = i
                            fmm_next = opt_list[i].fmm_o
        # Add the next best Option to the sorted list
        sorted_list = np.append(sorted_list,opt_list[index_next])
        fmm_min = fmm_next
        min_index = index_next
    # Return the maximin-sorted list of Options
    return sorted_list


######################################
"""""""""   Graph the Data   """""""""
######################################

def Graph_Data(wb,Pop,num_it,dt,tot_time):
    """
    Inputs:
        wb = Excel workbook for collecting data
        Pop = Population, whose 'opt_list' array of Options will be graphed
        num_it = Current iteration number
        dt = computational run-time since last time stamp
        tot_time = total computational run-time
    Actions:
        Creates/Prints 3 objective function graphs in Python.
        Pastes PERCS and PCS graphs in Excel.
        Pastes important tri-system Option data in Excel.
    """
    # Print the time lapse and current time
    print ("Number of Iterations =", num_it)
    m = time.localtime()
    if m[3]<=12:
        hr = m[3]
        if m[3]==0: hr = 12
        ap = "AM"
    if m[3]>12: 
        hr = m[3]-12
        ap = "PM"
    print ("dt =", dt/60.0, "min  /  Time =", hr,":","%02.0f"%m[4],":","%02.0f"%m[5],ap)
    ####################################
    """ Make Individual System Plots """
    ####################################
    # Establish the tab for Excel graphing
    tab = repr(num_it)
    #----------------
    """ Core Plot """
    #----------------
    # Find the max W_RCP
    W_rcp_max = 0.0
    for i in range(len(Pop.opt_list)):
        if Pop.opt_list[i].W_rcp > W_rcp_max:
            W_rcp_max = Pop.opt_list[i].W_rcp + 0.4
    # Create figure for core loop graph
    fig = plt.figure(figsize=(6,4))
    ax = fig.gca()
    plt.xlim([11.0,W_rcp_max])
    plt.ylim([0.9,1.8])
    # Rotate through Pop's 'opt_list' of Options
    for i in range(np.size(Pop.opt_list)):
        # Plot Option on core loop graph
        plt.scatter(Pop.opt_list[i].W_rcp,Pop.opt_list[i].cost_1,s=10,c='w',edgecolors='k')
    plt.grid()
    plt.xlabel('RCP Work (MW)')
    plt.ylabel('Cost ($1x10^9)')
    plt.show()
    #-----------------
    """ PERCS Plot """
    #-----------------
    # Put obj. function data into arrays
    x = np.array((Pop.opt_list[0].cost_2))
    y = np.array((Pop.opt_list[0].dT_int))
    z = np.array((Pop.opt_list[0].alpha))
    for i in range(1,np.size(Pop.opt_list)):
        x = np.append(x,Pop.opt_list[i].cost_2)
        y = np.append(y,Pop.opt_list[i].dT_int)
        z = np.append(z,Pop.opt_list[i].alpha)
    # Create the zoomed-in PERCS graph in Python
    fig = plt.figure(figsize=(6,5))
    ax = fig.gca(projection='3d')
    ax.scatter(x,y,z,c='b',marker='o')
    ax.set_xlabel('Cost ($1x10^8)',linespacing=3.2)
    ax.set_ylabel('\ndT_int',linespacing=3.1)
    ax.set_zlabel('\nAlpha',linespacing=3.1)
    ax.view_init(azim=-110,elev=45)
    plt.show()
    # Paste the zoomed-in Python graph in Excel
    wb.sheets[tab].pictures.add(fig,name='percs_graph')
    wb.sheets[tab].pictures[0].left = 550
    wb.sheets[tab].pictures[0].top = 97
    wb.sheets[tab].pictures[0].height = 211
    wb.sheets[tab].pictures[0].width = 250
    #---------------
    """ PCS Plot """
    #---------------
    # Create figure for PCS loop graph
    fig2 = plt.figure(figsize=(6,4))
    ax = fig2.gca()
    plt.xlim([0.33,0.37])
    plt.ylim([0,13])
    ax.set_xticks(np.arange(0.33,0.37,0.005))
    ax.set_yticks(np.arange(0,13,2))
    # Rotate through Pop's 'opt_list' of Options
    for i in range(np.size(Pop.opt_list)):
        # Plot Option on PCS loop graph
        plt.scatter(Pop.opt_list[i].eff,Pop.opt_list[i].cost_3,s=10,
                    c=Pop.opt_list[i].color,edgecolors=Pop.opt_list[i].color)
    plt.grid()
    plt.xlabel('Efficiency')
    plt.ylabel('Cost ($1x10^9)')
    plt.show()
    # Paste the PCS Python graph in Excel
    wb.sheets[tab].pictures.add(fig2,name='pcs_graph')
    wb.sheets[tab].pictures[1].left = 550
    wb.sheets[tab].pictures[1].top = 322
    wb.sheets[tab].pictures[1].height = 211
    wb.sheets[tab].pictures[1].width = 316
    ############################
    """ Export Data to Excel """
    ############################
    # Export the computational run-time
    wb.sheets[tab].range('B1').value = dt/60.0 # min
    wb.sheets[tab].range('B2').value = tot_time/60.0 # min
    col = 3
    n_a = "-"
    # Paste all pertinent Option data in Excel
    for k in range(np.size(Pop.opt_list)):
        wb.sheets[tab].range(2,col).value = k + 1
        wb.sheets[tab].range(3,col).value = Pop.opt_list[k].fmm_o
        # Core obj. functions
        wb.sheets[tab].range(4,col).value = Pop.opt_list[k].W_rcp
        wb.sheets[tab].range(5,col).value = Pop.opt_list[k].cost_1
        wb.sheets[tab].range(6,col).value = Pop.opt_list[k].fmm_1
        # PERCS obj. functions
        wb.sheets[tab].range(7,col).value = Pop.opt_list[k].cost_2
        wb.sheets[tab].range(8,col).value = Pop.opt_list[k].dT_int
        wb.sheets[tab].range(9,col).value = Pop.opt_list[k].alpha
        wb.sheets[tab].range(10,col).value = Pop.opt_list[k].fmm_2
        # PCS obj. functions
        wb.sheets[tab].range(11,col).value = Pop.opt_list[k].eff
        wb.sheets[tab].range(12,col).value = Pop.opt_list[k].cost_3
        wb.sheets[tab].range(13,col).value = Pop.opt_list[k].fmm_3
        # Core optimization parameters (x1)
        wb.sheets[tab].range(14,col).value = Pop.opt_list[k].R_f
        wb.sheets[tab].range(15,col).value = Pop.opt_list[k].H_fuel
        wb.sheets[tab].range(16,col).value = Pop.opt_list[k].Dh_00
        wb.sheets[tab].range(17,col).value = Pop.opt_list[k].Dh_12
        wb.sheets[tab].range(18,col).value = Pop.opt_list[k].Dh_14
        # PERCS optimization parameters (x2)
        wb.sheets[tab].range(19,col).value = Pop.opt_list[k].R_tank
        wb.sheets[tab].range(20,col).value = Pop.opt_list[k].pitch
        wb.sheets[tab].range(21,col).value = Pop.opt_list[k].D_h
        wb.sheets[tab].range(22,col).value = Pop.opt_list[k].th
        wb.sheets[tab].range(23,col).value = Pop.opt_list[k].Len
        wb.sheets[tab].range(24,col).value = Pop.opt_list[k].elev
        # PCS optimization parameters (x3)
        wb.sheets[tab].range(25,col).value = Pop.opt_list[k].phx.Tout
        wb.sheets[tab].range(26,col).value = Pop.opt_list[k].t1a.Pout
        wb.sheets[tab].range(27,col).value = Pop.opt_list[k].mf_t1a
        if Pop.opt_list[k].y_rh1 == 0:
            wb.sheets[tab].range(26,col).color = (139,137,137)
            wb.sheets[tab].range(27,col).color = (139,137,137)
        wb.sheets[tab].range(28,col).value = Pop.opt_list[k].t1b.Pout
        wb.sheets[tab].range(29,col).value = Pop.opt_list[k].mf_t1b
        if Pop.opt_list[k].y_s5 == 0:
            wb.sheets[tab].range(28,col).color = (139,137,137)
            wb.sheets[tab].range(29,col).color = (139,137,137)
        wb.sheets[tab].range(30,col).value = Pop.opt_list[k].t1c.Pout
        wb.sheets[tab].range(31,col).value = Pop.opt_list[k].t2a.Pout
        wb.sheets[tab].range(32,col).value = Pop.opt_list[k].mf_t2a
        if Pop.opt_list[k].y_s14 == 0:
            wb.sheets[tab].range(31,col).color = (139,137,137)
            wb.sheets[tab].range(32,col).color = (139,137,137)
        wb.sheets[tab].range(33,col).value = Pop.opt_list[k].t2b.Pout
        if Pop.opt_list[k].y_ipt == 0:
            wb.sheets[tab].range(33,col).color = (139,137,137)
        # PCS optimization parameters (y)
        wb.sheets[tab].range(34,col).value = Pop.opt_list[k].color
        wb.sheets[tab].range(35,col).value = Pop.opt_list[k].y_ipt
        wb.sheets[tab].range(36,col).value = Pop.opt_list[k].y_rh1
        wb.sheets[tab].range(37,col).value = Pop.opt_list[k].y_rh2
        wb.sheets[tab].range(38,col).value = Pop.opt_list[k].y_s14
        wb.sheets[tab].range(39,col).value = Pop.opt_list[k].y_s4
        wb.sheets[tab].range(40,col).value = Pop.opt_list[k].y_s5
        # Core Constraints Data
        wb.sheets[tab].range(41,col).value = Pop.opt_list[k].T_f_over_max
        wb.sheets[tab].range(42,col).value = Pop.opt_list[k].T_c_over_max
        wb.sheets[tab].range(43,col).value = Pop.opt_list[k].MDNBR_below_1
        wb.sheets[tab].range(44,col).value = Pop.opt_list[k].penalized
        wb.sheets[tab].range(45,col).value = Pop.opt_list[k].failed
        wb.sheets[tab].range(46,col).value = max(Pop.opt_list[k].T_1336_1)
        wb.sheets[tab].range(47,col).value = max(Pop.opt_list[k].T_1336_17)
        wb.sheets[tab].range(48,col).value = Pop.opt_list[k].MDNBR
        # PERCS Constraints Data
        wb.sheets[tab].range(49,col).value = Pop.opt_list[k].T_over_620
        wb.sheets[tab].range(50,col).value = Pop.opt_list[k].T_over_635
        wb.sheets[tab].range(51,col).value = Pop.opt_list[k].PERCS_failed
        # Additional Core Data
        wb.sheets[tab].range(53,col).value = Pop.opt_list[k].m_dot_100
        wb.sheets[tab].range(54,col).value = Pop.opt_list[k].m_dot_400
        wb.sheets[tab].range(55,col).value = Pop.opt_list[k].m_dot_600
        wb.sheets[tab].range(56,col).value = Pop.opt_list[k].m_dot_200
        wb.sheets[tab].range(57,col).value = Pop.opt_list[k].m_dot_335
        wb.sheets[tab].range(58,col).value = Pop.opt_list[k].k_eff
        wb.sheets[tab].range(59,col).value = Pop.opt_list[k].rho_0
        wb.sheets[tab].range(60,col).value = Pop.opt_list[k].Bc
        wb.sheets[tab].range(61,col).value = Pop.opt_list[k].nBc
        # Additional PERCS Data
        wb.sheets[tab].range(63,col).value = Pop.opt_list[k].n_tubes
        wb.sheets[tab].range(64,col).value = Pop.opt_list[k].m_MgCO3
        wb.sheets[tab].range(65,col).value = Pop.opt_list[k].p716.cost+Pop.opt_list[k].p717.cost
        wb.sheets[tab].range(66,col).value = Pop.opt_list[k].support.cost
        wb.sheets[tab].range(67,col).value = Pop.opt_list[k].hx.cost
        wb.sheets[tab].range(68,col).value = Pop.opt_list[k].tank.cost
        # Additional PCS Data
        wb.sheets[tab].range(70,col).value = Pop.opt_list[k].phx.Pout
        wb.sheets[tab].range(71,col).value = Pop.opt_list[k].phx.Tin
        wb.sheets[tab].range(72,col).value = Pop.opt_list[k].phx.mdot
        wb.sheets[tab].range(73,col).value = Pop.opt_list[k].t3.Pout
        wb.sheets[tab].range(74,col).value = Pop.opt_list[k].t1.W
        if Pop.opt_list[k].y_ipt == 1:
            wb.sheets[tab].range(75,col).value = Pop.opt_list[k].t2.W
        else:
            wb.sheets[tab].range(75,col).value = n_a
        wb.sheets[tab].range(76,col).value = Pop.opt_list[k].t3.W
        wb.sheets[tab].range(77,col).value = Pop.opt_list[k].t4.W
        wb.sheets[tab].range(78,col).value = Pop.opt_list[k].t5.W
        wb.sheets[tab].range(79,col).value = Pop.opt_list[k].p1.W
        if Pop.opt_list[k].p2.y == 1:
            wb.sheets[tab].range(80,col).value = Pop.opt_list[k].p2.W
        else:
            wb.sheets[tab].range(80,col).value = n_a
        if Pop.opt_list[k].p3.y == 1:
            wb.sheets[tab].range(81,col).value = Pop.opt_list[k].p3.W
        else:
            wb.sheets[tab].range(81,col).value = n_a
        if Pop.opt_list[k].p4.y == 1:
            wb.sheets[tab].range(82,col).value = Pop.opt_list[k].p4.W
        else:
            wb.sheets[tab].range(82,col).value = n_a
        wb.sheets[tab].range(83,col).value = Pop.opt_list[k].p5.W
        # Increment the column number between options
        col = col + 1
    # No need to return anything
    return None


#######################################################
"""""""""   Specialized Offspring Functions   """""""""
#######################################################

def autogen_k_lists(pop,num_keepers):
    """
    Inputs:
        pop = Population with current 'opt_list' array of Options
        num_keepers = Number of single-system designs to keep per system
    Outputs:
        k_list1 = Array of 'opt_list' indices of single-system core loop designs to keep
        k_list2 = Array of 'opt_list' indices of single-system PERCS designs to keep
        k_list3 = Array of 'opt_list' indices of single-system PCS designs to keep
    """
    n = len(pop.opt_list)
    # Initiate the k_lists
    k_list1 = np.zeros(0,dtype=int)
    k_list2 = np.zeros(0,dtype=int)
    k_list3 = np.zeros(0,dtype=int)
    #----------------
    """ Core Loop """
    #----------------
    # Initiate variables for building k_list1
    fmm_min = 1000
    min_index = 100
    # Rotate through each Option
    for i in range(n):
        # Keep track of the smallest f_maximin
        if pop.opt_list[i].fmm_1 < fmm_min:
            fmm_min = pop.opt_list[i].fmm_1
            min_index = i
    # Add the min_index to k_list1
    k_list1 = np.append(k_list1,min_index)
    # Add the correct number of Options in ascending "maximin" order
    for count in range(num_keepers-1):
        fmm_next = 1000
        # Rotate through the Options
        for i in range(n):
            # Find the next smallest f_maximin
            if i != min_index:
                if pop.opt_list[i].fmm_1 < fmm_next and pop.opt_list[i].fmm_1 >= fmm_min:
                    if pop.opt_list[i].fmm_1 == fmm_min and i > min_index:
                        index_next = i
                        fmm_next = pop.opt_list[i].fmm_1
                        break
                    else:
                        if pop.opt_list[i].fmm_1 == fmm_min and i < min_index:
                            None
                        else:
                            index_next = i
                            fmm_next = pop.opt_list[i].fmm_1
        # Add the next best Option to k_list1
        k_list1 = np.append(k_list1,index_next)
        fmm_min = fmm_next
        min_index = index_next
    #-----------------
    """ PERCS Loop """
    #-----------------
    # Initiate variables for building k_list2
    fmm_min = 1000
    min_index = 100
    # Rotate through each Option
    for i in range(n):
        # Keep track of the smallest f_maximin
        if pop.opt_list[i].fmm_2 < fmm_min:
            fmm_min = pop.opt_list[i].fmm_2
            min_index = i
    # Add the min_index to k_list2
    k_list2 = np.append(k_list2,min_index)
    # Add the correct number of Options in ascending "maximin" order
    for count in range(num_keepers-1):
        fmm_next = 1000
        # Rotate through the Options
        for i in range(n):
            # Find the next smallest f_maximin
            if i != min_index:
                if pop.opt_list[i].fmm_2 < fmm_next and pop.opt_list[i].fmm_2 >= fmm_min:
                    if pop.opt_list[i].fmm_2 == fmm_min and i > min_index:
                        index_next = i
                        fmm_next = pop.opt_list[i].fmm_2
                        break
                    else:
                        if pop.opt_list[i].fmm_2 == fmm_min and i < min_index:
                            None
                        else:
                            index_next = i
                            fmm_next = pop.opt_list[i].fmm_2
        # Add the next best Option to k_list2
        k_list2 = np.append(k_list2,index_next)
        fmm_min = fmm_next
        min_index = index_next
    #---------------
    """ PCS Loop """
    #---------------
    # Initiate variables for building k_list3
    fmm_min = 1000
    min_index = 100
    # Rotate through each Option
    for i in range(n):
        # Keep track of the smallest f_maximin
        if pop.opt_list[i].fmm_3 < fmm_min:
            fmm_min = pop.opt_list[i].fmm_3
            min_index = i
    # Add the min_index to k_list3
    k_list3 = np.append(k_list3,min_index)
    # Add the correct number of Options in ascending "maximin" order
    for count in range(num_keepers-1):
        fmm_next = 1000
        # Rotate through the Options
        for i in range(n):
            # Find the next smallest f_maximin
            if i != min_index:
                if pop.opt_list[i].fmm_3 < fmm_next and pop.opt_list[i].fmm_3 >= fmm_min:
                    if pop.opt_list[i].fmm_3 == fmm_min and i > min_index:
                        index_next = i
                        fmm_next = pop.opt_list[i].fmm_3
                        break
                    else:
                        if pop.opt_list[i].fmm_3 == fmm_min and i < min_index:
                            None
                        else:
                            index_next = i
                            fmm_next = pop.opt_list[i].fmm_3
        # Add the next best Option to k_list3
        k_list3 = np.append(k_list3,index_next)
        fmm_min = fmm_next
        min_index = index_next
    
    # Return the filled k_lists
    return k_list1,k_list2,k_list3


def generate_k_lists(pop):
    """
    Inputs:
        pop = Population with current 'opt_list' array of Options
    Outputs:
        k_list1 = Array of 'opt_list' indices of single-system core loop designs to keep
        k_list2 = Array of 'opt_list' indices of single-system PERCS designs to keep
        k_list3 = Array of 'opt_list' indices of single-system PCS designs to keep
    """
    # Initiate the k_lists
    k_list1 = np.zeros(0,dtype=int)
    k_list2 = np.zeros(0,dtype=int)
    k_list3 = np.zeros(0,dtype=int)
    # Start by adding all non-dominated single-systems to respective k_lists
    for i in range(len(pop.opt_list)):
        if pop.opt_list[i].fmm_1 < 0.0:
            k_list1 = np.append(k_list1,i)
        if pop.opt_list[i].fmm_2 < 0.0:
            k_list2 = np.append(k_list2,i)
        if pop.opt_list[i].fmm_3 < 0.0:
            k_list3 = np.append(k_list3,i)
    
    proceed = 'no'
    while proceed == 'no':
        #----------------
        """ Core Loop """
        #----------------
        # Show initial core loop graph
        graph1(pop,k_list1)
        # Enter the AOP choices loop
        done1 = False
        while done1 == False:
            print("What would you like to do?")
            print("(1) add, (2) remove, (3) see graph, (4) move on, (5) see zoomed-in graph")
            action = eval(input("Action #: "))
            # Action 1: Add to keep list #1
            if action == 1:
                print("Give the range you want to look at")
                lb = eval(input("W_rcp lower bound = "))
                ub = eval(input("W_rcp upper bound = "))
                print(" ID         W_rcp            Cost")
                print("------------------------------------------")
                for i in range(len(pop.opt_list)):
                    # If opt.W_rcp is in the range
                    if pop.opt_list[i].W_rcp > lb and pop.opt_list[i].W_rcp < ub:
                        # Make sure it's not already in k_list1
                        not_kept = True
                        for v in k_list1:
                            if i == v: not_kept = False
                        # If it is not already kept, then print option data
                        if not_kept == True:
                            print(" ",i," ",pop.opt_list[i].W_rcp," ",pop.opt_list[i].cost_1)
                print()
                print("Which do you want to keep? [Ex: 3 6 12 19]")
                add_list = input()
                nums = add_list.split()
                for num in nums:
                    a = int(num)
                    added = False
                    for i in range(len(k_list1)):
                        if a < k_list1[i]:
                            k_list1 = np.insert(k_list1,i,a)
                            added = True
                            break
                    if added == False:
                        k_list1 = np.append(k_list1,a)
            # Action 2: Remove from keep list #1
            if action == 2:
                print(" ID         W_rcp            Cost")
                print("------------------------------------------")
                for k in k_list1:
                    print(" ",k," ",pop.opt_list[k].W_rcp," ",pop.opt_list[k].cost_1)
                print()
                print("Which do you want to remove? [Ex: 0 2 5 13]")
                remove_list = input()
                nums = remove_list.split()
                for num in nums:
                    r = int(num)
                    for i in range(len(k_list1)):
                        if r == k_list1[i]:
                            k_list1 = np.delete(k_list1,i)
                            break
            # Action 3: View core loop graph
            if action == 3:
                graph1(pop,k_list1)
            # Action 4: Move on from Core Loop
            if action == 4:
                done1 = True
            # Action 5: View zoomed-in core loop graph
            if action == 5:
                maxes = np.zeros(2,dtype=float)
                maxes[0] = eval(input("W_rcp upper bound = "))
                maxes[1] = eval(input("cost upper bound ="))
                graph1_zoom(pop,k_list1,maxes)
        #-----------------
        """ PERCS Loop """
        #-----------------
        # Show initial PERCS graph
        graph2(pop,k_list2)
        # Enter the AOP choices loop
        done2 = False
        while done2 == False:
            print("What would you like to do?")
            print("(1) add, (2) remove, (3) see graph, (4) move on, (5) see zoomed-in graph")
            action = eval(input("Action #: "))
            # Action 1: Add to keep list #2
            if action == 1:
                print("Give the range you want to look at")
                lb = eval(input("cost lower bound = "))
                ub = eval(input("cost upper bound = "))
                print(" ID         Cost            dT_int            Alpha")
                print("-----------------------------------------------------------------")
                for i in range(len(pop.opt_list)):
                    # If opt.zdt2_f1 is in the range
                    if pop.opt_list[i].cost_2 > lb and pop.opt_list[i].cost_2 < ub:
                        # Make sure it's not already in k_list2
                        not_kept = True
                        for v in k_list2:
                            if i == v: not_kept = False
                        # If it is not already kept, then print option data
                        if not_kept == True:
                            print(" ",i," ",pop.opt_list[i].cost_2," ",pop.opt_list[i].dT_int," ",pop.opt_list[i].alpha)
                print()
                print("Which do you want to keep? [Ex: 3 6 12 19]")
                add_list = input()
                nums = add_list.split()
                for num in nums:
                    a = int(num)
                    added = False
                    for i in range(len(k_list2)):
                        if a < k_list2[i]:
                            k_list2 = np.insert(k_list2,i,a)
                            added = True
                            break
                    if added == False:
                        k_list2 = np.append(k_list2,a)
            # Action 2: Remove from keep list #2
            if action == 2:
                print(" ID         Cost            dT_int            Alpha")
                print("-----------------------------------------------------------------")
                for k in k_list2:
                    print(" ",k," ",pop.opt_list[k].cost_2," ",pop.opt_list[k].dT_int," ",pop.opt_list[k].alpha)
                print()
                print("Which do you want to remove? [Ex: 0 2 5 13]")
                remove_list = input()
                nums = remove_list.split()
                for num in nums:
                    r = int(num)
                    for i in range(len(k_list2)):
                        if r == k_list2[i]:
                            k_list2 = np.delete(k_list2,i)
                            break
            # Action 3: View PERCS graph
            if action == 3:
                graph2(pop,k_list2)
            # Action 4: Move on from PERCS
            if action == 4:
                done2 = True
            # Action 5: View zoomed-in PERCS graph
            if action == 5:
                maxes = np.zeros(3,dtype=float)
                maxes[0] = eval(input("cost upper bound = "))
                maxes[1] = eval(input("dT_int upper bound ="))
                maxes[2] = eval(input("alpha upper bound ="))
                graph2_zoom(pop,k_list2,maxes)
        #---------------
        """ PCS Loop """
        #---------------
        # Show initial PERCS graph
        graph3(pop,k_list3)
        # Enter the AOP choices loop
        done3 = False
        while done3 == False:
            print("What would you like to do?")
            print("(1) add, (2) remove, (3) see graph, (4) move on")
            action = eval(input("Action #: "))
            # Action 1: Add to keep list #3
            if action == 1:
                print("Give the range you want to look at")
                lb = eval(input("eff lower bound = "))
                ub = eval(input("eff upper bound = "))
                print(" ID          eff                  Cost               Color")
                print("----------------------------------------------------------------")
                for i in range(len(pop.opt_list)):
                    # If opt.eff is in the range
                    if pop.opt_list[i].eff > lb and pop.opt_list[i].eff < ub:
                        # Make sure it's not already in k_list2
                        not_kept = True
                        for v in k_list3:
                            if i == v: not_kept = False
                        # If it is not already kept, then print option data
                        if not_kept == True:
                            print(" ",i," ",pop.opt_list[i].eff," ",pop.opt_list[i].cost_3,
                                  " ",pop.opt_list[i].color)
                print()
                print("Which do you want to keep? [Ex: 3 6 12 19]")
                add_list = input()
                nums = add_list.split()
                for num in nums:
                    a = int(num)
                    added = False
                    for i in range(len(k_list3)):
                        if a < k_list3[i]:
                            k_list3 = np.insert(k_list3,i,a)
                            added = True
                            break
                    if added == False:
                        k_list3 = np.append(k_list3,a)
            # Action 2: Remove from keep list #3
            if action == 2:
                print(" ID          eff                  Cost               Color")
                print("----------------------------------------------------------------")
                for k in k_list3:
                    print(" ",k," ",pop.opt_list[k].eff," ",pop.opt_list[k].cost_3,
                          " ",pop.opt_list[k].color)
                print()
                print("Which do you want to remove? [Ex: 0 2 5 13]")
                remove_list = input()
                nums = remove_list.split()
                for num in nums:
                    r = int(num)
                    for i in range(len(k_list3)):
                        if r == k_list3[i]:
                            k_list3 = np.delete(k_list3,i)
                            break
            # Action 3: View PCS graph
            if action == 3:
                graph3(pop,k_list3)
            # Action 4: Move on from PCS
            if action == 4:
                done3 = True
        # Is the DM ready to move on to the Crossover phase?
        proceed = input("Ready to proceed with Crossover phase? [yes/no]  ")
    # Return the filled k_lists
    return k_list1,k_list2,k_list3

# Create a graph for the Core
def graph1(p,k_list):
    print("Core Graph")
    x_max = 0.0
    y_max = 0.0
    # Find the maximums
    for opt in p.opt_list:
        if opt.W_rcp > x_max: x_max = opt.W_rcp
        if opt.cost_1 > y_max: y_max = opt.cost_1
    x_max = x_max + 1
    y_max = y_max + 1
    fig = plt.figure(figsize=(6,4))
    ax = fig.gca()
    plt.xlim([11,x_max]) # f1
    plt.ylim([0.5,y_max]) # f2
    ax.set_xticks(np.arange(11,x_max,(x_max-11)/10.0))
    ax.set_yticks(np.arange(0.5,y_max,(y_max-0.5)/10.0))
    for i in range(len(p.opt_list)):
        kept = False
        for v in k_list:
            if i == v:
                plt.scatter(p.opt_list[i].W_rcp,p.opt_list[i].cost_1,s=10,c='w',edgecolors='r')
                kept = True
        if kept == False:
            plt.scatter(p.opt_list[i].W_rcp,p.opt_list[i].cost_1,s=10,c='w',edgecolors='k')
    plt.grid()
    plt.xlabel('W_rcp (MW)')
    plt.ylabel('Cost ($1x10^9)')
    plt.show()

# Create a zoomed graph for the Core
def graph1_zoom(p,k_list,maxes):
    print("Zoomed-in Core Graph")
    x_min = 9.9e99
    y_min = 9.9e99
    # Find the minimums
    for opt in p.opt_list:
        if opt.W_rcp < x_min: x_min = opt.W_rcp
        if opt.cost_1 < y_min: y_min = opt.cost_1
    fig = plt.figure(figsize=(6,4))
    ax = fig.gca()
    plt.xlim([x_min-0.02,maxes[0]])
    plt.ylim([y_min-0.1,maxes[1]])
    ax.set_xticks(np.arange(x_min-0.02,maxes[0],(maxes[0]-x_min)/10.0))
    ax.set_yticks(np.arange(y_min-0.1,maxes[1],(maxes[1]-y_min)/10.0))
    for i in range(len(p.opt_list)):
        kept = False
        for v in k_list:
            if i == v:
                plt.scatter(p.opt_list[i].W_rcp,p.opt_list[i].cost_1,s=10,c='w',edgecolors='r')
                kept = True
        if kept == False:
            plt.scatter(p.opt_list[i].W_rcp,p.opt_list[i].cost_1,s=10,c='w',edgecolors='k')
    plt.grid()
    plt.xlabel('W_rcp (MW)')
    plt.ylabel('Cost ($1x10^9)')
    plt.show()

# Create a graph for the PERCS
def graph2(p,k_list):
    print("PERCS Graph")
    fig = plt.figure(figsize=(6,5))
    ax = fig.gca(projection='3d')
    for i in range(len(p.opt_list)):
        kept = False
        for v in k_list:
            if i == v:
                ax.scatter(p.opt_list[i].cost_2,p.opt_list[i].dT_int,p.opt_list[i].alpha,
                            c='w',edgecolors='r',marker='o')
                kept = True
        if kept == False:
            ax.scatter(p.opt_list[i].cost_2,p.opt_list[i].dT_int,p.opt_list[i].alpha,
                        c='w',edgecolors='k',marker='o',alpha=0.5)
    ax.view_init(azim=-110,elev=45)
    ax.set_xlabel('Cost ($1x10^8)',linespacing=3.2)
    ax.set_ylabel('\ndT_int',linespacing=3.1)
    ax.set_zlabel('\nAlpha',linespacing=3.1)
    plt.show()

# Create a zoomed graph for the PERCS
def graph2_zoom(p,k_list,maxes):
    print("Zoomed-in PERCS Graph")
    x_min = 9.9e99
    y_min = 9.9e99
    z_min = 9.9e99
    # Find the minimums
    for opt in p.opt_list:
        if opt.cost_2 < x_min: x_min = opt.cost_2
        if opt.dT_int < y_min: y_min = opt.dT_int
        if opt.alpha < z_min: z_min = opt.alpha
    fig = plt.figure(figsize=(6,5))
    ax = fig.gca(projection='3d')
    ax.set_xlim3d(x_min,maxes[0])
    ax.set_ylim3d(y_min,maxes[1])
    ax.set_zlim3d(z_min,maxes[2])
    for i in range(len(p.opt_list)):
        kept = False
        for v in k_list:
            if i == v:
                ax.scatter(p.opt_list[i].cost_2,p.opt_list[i].dT_int,p.opt_list[i].alpha,
                            c='w',edgecolors='r',marker='o')
                kept = True
        if kept == False:
            ax.scatter(p.opt_list[i].cost_2,p.opt_list[i].dT_int,p.opt_list[i].alpha,
                        c='w',edgecolors='k',marker='o',alpha=0.5)
    ax.view_init(azim=-110,elev=45)
    ax.set_xlabel('Cost ($1x10^8)',linespacing=3.2)
    ax.set_ylabel('\ndT_int',linespacing=3.1)
    ax.set_zlabel('\nAlpha',linespacing=3.1)
    plt.show()

# Create a graph for the PCS
def graph3(p,k_list):
    print("PCS Graph")
    x_max = 0.0
    y_max = 0.0
    # Find the maximums
    for opt in p.opt_list:
        if opt.eff > x_max: x_max = opt.eff
        if opt.cost_3 > y_max: y_max = opt.cost_3
    x_max = x_max + 0.005
    y_max = y_max + 1
    fig = plt.figure(figsize=(6,4))
    ax = fig.gca()
    plt.xlim([0.33,x_max]) # f1
    plt.ylim([0,y_max]) # f2
    ax.set_xticks(np.arange(0.33,x_max,0.005))
    ax.set_yticks(np.arange(0,y_max,2))
    for i in range(len(p.opt_list)):
        kept = False
        for v in k_list:
            if i == v:
                plt.scatter(p.opt_list[i].eff,p.opt_list[i].cost_3,s=12,
                            c=p.opt_list[i].color,edgecolors=p.opt_list[i].color)
                kept = True
        if kept == False:
            plt.scatter(p.opt_list[i].eff,p.opt_list[i].cost_3,s=12,
                        c='w',edgecolors=p.opt_list[i].color)
    plt.grid()
    plt.xlabel('Efficiency')
    plt.ylabel('Cost ($1x10^9)')
    plt.show()


############################################
"""""""""   Tough Love Functions   """""""""
############################################

def tough_love(opt_list):
    """
    Inputs:
        opt_list = Population's 'opt_list' array of tri-system Options
    Output:
        opt_list = Revised 'opt_list' array w/ TL'ed Options deleted
        save_x1s = Matrix of x1-arrays to save after TL
        save_x2s = Matrix of x2-arrays to save after TL
        save_x3s = Matrix of x3-arrays to save after TL
        save_ys = Matrix of y-arrays (linked with saved x3-arrays)
    """
    # Initiate the deleted_list
    deleted_list = np.zeros(0)
    # Initiate the non-dominated lists
    nd_list1 = np.zeros(0,dtype=int)
    nd_list2 = np.zeros(0,dtype=int)
    nd_list3 = np.zeros(0,dtype=int)
    for i in range(len(opt_list)):
        if opt_list[i].fmm_1 < 0.0:
            nd_list1 = np.append(nd_list1,i)
        if opt_list[i].fmm_2 < 0.0:
            nd_list2 = np.append(nd_list2,i)
        if opt_list[i].fmm_3 < 0.0:
            nd_list3 = np.append(nd_list3,i)
    # Show 3 initial graphs
    graph1_2(opt_list,nd_list1)
    graph2_2(opt_list,nd_list2)
    graph3_2(opt_list,nd_list3)
    #----------------------------
    """ Option Deletion Phase """
    #----------------------------
    done = False
    while done == False:
        # Prompt DM to make an action choice
        print("Tough Love Step: What would you like to do?")
        print("(1) see tables, (2) see graphs, (3) delete, (4) undelete, (5) move on")
        action = eval(input("Action #: "))
        # Action 1: View tables of objective function data
        if action == 1:
            print(" ID  fmm_o   fmm_1   fmm_2   fmm_3")
            print("----------------------------------")
            for i in range(len(opt_list)):
                print(" ",i," ",np.round(opt_list[i].fmm_o,5)," ",np.round(opt_list[i].fmm_1,5),
                      " ",np.round(opt_list[i].fmm_2,5)," ",np.round(opt_list[i].fmm_3,5))
            print()
            for i in range(len(opt_list)):
                print("-----------------------------------------------")
                print("opt_list index =",i)
                print("W_rcp =",np.round(opt_list[i].W_rcp,4),", Cost_1 =",np.round(opt_list[i].cost_1,4))
                print("Cost_2 =",np.round(opt_list[i].cost_2,4),", dT_int =",
                      np.round(opt_list[i].dT_int,2),", Alpha =",np.round(opt_list[i].alpha,7))
                print("Eff =",np.round(opt_list[i].eff,4),", Cost_3 =",
                      np.round(opt_list[i].cost_3,4),", color =",opt_list[i].color)
            print("-----------------------------------------------")
        # Action 2: View 3 graphs again
        if action == 2:
            # Show the three graphs with normal zooming
            graph1_2(opt_list,[])
            graph2_2(opt_list,[])
            graph3_2(opt_list,[])
            # Give the option to create a zoomed-in PERCS graph
            zoom_in = input("Would you like to see a zoomed-in PERCS graph? [yes/no] ")
            while zoom_in == 'yes':
                maxes = np.zeros(3,dtype=float)
                maxes[0] = eval(input("cost_2 upper bound = "))
                maxes[1] = eval(input("dT_int upper bound ="))
                maxes[2] = eval(input("alpha upper bound ="))
                graph2_2_zoom(opt_list,[],maxes)
                # Allow for another zoom-in if DM is not satisfied with the graph
                zoom_in = input("Zoom in further? [yes/no] ")
        # Action 3: Name Options to delete from the opt_list
        if action == 3:
            print("Which do you want to delete? [Ex: 0 2 5 13]")
            remove_list = input()
            nums = remove_list.split()
            # Go backwards through the list so the opt_list indeces match up correctly
            for d in range(len(nums)-1,-1,-1):
                if np.size(deleted_list) == 0:
                    deleted_list = np.array((opt_list[int(nums[d])]))
                else:
                    deleted_list = np.append(deleted_list,opt_list[int(nums[d])])
                opt_list = np.delete(opt_list,int(nums[d]))
            # Remind the DM that the opt_list got smaller and IDs may have shifted
            n_less = len(nums)
            print("REMEMBER: The list of options just got smaller by ",n_less," options.")
            print("          This may have shifted the ID values of the remaining options.")
            print("          To see the new ID values, select action #1.")
        # Action 4: Un-Delete some Options and return them to opt_list
        if action == 4:
            if len(deleted_list) > 0:
                # Print the obj. func. table for the current deleted_list
                print(" ID  fmm_o   fmm_1   fmm_2   fmm_3")
                print("----------------------------------")
                for i in range(len(deleted_list)):
                    print(" ",i," ",np.round(deleted_list[i].fmm_o,5)," ",np.round(deleted_list[i].fmm_1,5),
                          " ",np.round(deleted_list[i].fmm_2,5)," ",np.round(deleted_list[i].fmm_3,5))
                print()
                print("Which do you want to undelete? [Ex: 0 2 5 13]")
                add_list = input()
                nums = add_list.split()
                # Go backwards through the list so the deleted_list indeces match up correctly
                for a in range(len(nums)-1,-1,-1):
                    opt_list = np.append(opt_list,deleted_list[a])
                    deleted_list = np.delete(deleted_list,a)
            else:
                print("The deleted_list is currently empty.")
        # Action 5: Move on from TL step
        if action == 5:
            done = True
            break
    #---------------------------------
    """ Single-System Saving Phase """
    #---------------------------------
    # Phase used to save offspring data from parts of deleted Options
    # Create exmpty matrices for offspring data to save
    save_x1s = np.zeros(0)
    save_x2s = np.zeros(0)
    save_x3s = np.zeros(0)
    save_ys = np.zeros(0)
    # Ask DM if they want to save any deleted_list data to generate offspring
    gen_off = input("Do you want to save any data to generate offspring with? [yes/no] ")
    if gen_off == 'yes':
        # Print off the current deleted list data
        print(" ID  fmm_o   fmm_1   fmm_2   fmm_3")
        print("----------------------------------")
        for i in range(len(deleted_list)):
            print(" ",i," ",np.round(deleted_list[i].fmm_o,5)," ",np.round(deleted_list[i].fmm_1,5),
                  " ",np.round(deleted_list[i].fmm_2,5)," ",np.round(deleted_list[i].fmm_3,5))
        for i in range(len(deleted_list)):
            print("-----------------------------------------------")
            print("deleted_list index =",i)
            print("W_rcp =",np.round(deleted_list[i].W_rcp,4),", Cost_1 =",np.round(deleted_list[i].cost_1,4))
            print("Cost_2 =",np.round(deleted_list[i].cost_2,4),", dT_int =",
                  np.round(deleted_list[i].dT_int,2),", Alpha =",np.round(deleted_list[i].alpha,7))
            print("Eff =",np.round(deleted_list[i].eff,4),", Cost_3 =",
                  np.round(deleted_list[i].cost_3,4),", color =",deleted_list[i].color)
        print("-----------------------------------------------")
        print()
        # Ask DM if they want to save any Core loop data
        save_1 = input("Do you want to save any Core loop options? [yes/no] ")
        if save_1 == 'yes':
            # Select which Core loop options to save
            print("Which do you want to save? [Ex: 0 2 5 13]")
            save_list = input()
            nums = save_list.split()
            # Rotate through IDs to save and add to save_x1s array
            for num in nums:
                n = int(num)
                if len(save_x1s)==0:
                    save_x1s = np.zeros((1,5))
                    save_x1s[0,:] = deleted_list[n].x1
                else:
                    save_x1s = np.vstack((save_x1s,deleted_list[n].x1))
        # Ask DM if they want to save any PERCS loop data
        save_2 = input("Do you want to save any PERCS loop options? [yes/no] ")
        if save_2 == 'yes':
            # Select which PERCS loop options to save
            print("Which do you want to save? [Ex: 0 2 5 13]")
            save_list = input()
            nums = save_list.split()
            # Rotate through IDs to save and add to save_x2s array
            for num in nums:
                n = int(num)
                if len(save_x2s)==0:
                    save_x2s = np.zeros((1,6))
                    save_x2s[0,:] = deleted_list[n].x2
                else:
                    save_x2s = np.vstack((save_x2s,deleted_list[n].x2))
        # Ask DM if they want to save any PCS loop data
        save_3 = input("Do you want to save any PCS loop options? [yes/no] ")
        if save_3 == 'yes':
            # Select which PCS loop options to save
            print("Which do you want to save? [Ex: 0 2 5 13]")
            save_list = input()
            nums = save_list.split()
            # Rotate through IDs to save and add to save_x3s and save_ys arrays
            for num in nums:
                n = int(num)
                if len(save_x3s)==0:
                    save_x3s = np.zeros((1,9))
                    save_x3s[0,:] = deleted_list[n].x3
                    save_ys = np.zeros((1,6))
                    save_ys[0,:] = deleted_list[n].y
                else:
                    save_x3s = np.vstack((save_x3s,deleted_list[n].x3))
                    save_ys = np.vstack((save_ys,deleted_list[n].y))
    # Return the opt_list and the deleted_list data desired for generating offspring
    # Forget the deleted_list by not returning it
    return opt_list, save_x1s,save_x2s,save_x3s,save_ys

# Create a graph for the Core
def graph1_2(opt_list,nd_list):
    print("Core Graph")
    x_max = 0.0
    y_max = 0.0
    # Find the maximums
    for opt in opt_list:
        if opt.W_rcp > x_max: x_max = opt.W_rcp
        if opt.cost_1 > y_max: y_max = opt.cost_1
    x_max = x_max + 1
    y_max = y_max + 1
    fig = plt.figure(figsize=(6,4))
    ax = fig.gca()
    plt.xlim([11,x_max]) # W_rcp
    plt.ylim([0.5,y_max]) # cost_1
    ax.set_xticks(np.arange(11,x_max,(x_max-11)/10.0))
    ax.set_yticks(np.arange(0.5,y_max,(y_max-0.5)/10.0))
    for i in range(len(opt_list)):
        non_dom = False
        for v in nd_list:
            if i == v:
                plt.scatter(opt_list[i].W_rcp,opt_list[i].cost_1,s=10,c='w',edgecolors='r')
                non_dom = True
        if non_dom == False:
            plt.scatter(opt_list[i].W_rcp,opt_list[i].cost_1,s=10,c='w',edgecolors='k')
    plt.grid()
    plt.xlabel('W_rcp (MW)')
    plt.ylabel('Cost ($1x10^9)')
    plt.show()

# Create a graph for the PERCS
def graph2_2(opt_list,nd_list):
    print("PERCS Graph")
    fig = plt.figure(figsize=(6,5))
    ax = fig.gca(projection='3d')
    for i in range(len(opt_list)):
        non_dom = False
        for v in nd_list:
            if i == v:
                ax.scatter(opt_list[i].cost_2,opt_list[i].dT_int,opt_list[i].alpha,
                            c='w',edgecolors='r',marker='o')
                non_dom = True
        if non_dom == False:
            ax.scatter(opt_list[i].cost_2,opt_list[i].dT_int,opt_list[i].alpha,
                        c='w',edgecolors='k',marker='o',alpha=0.5)
    ax.view_init(azim=-110,elev=45)
    ax.set_xlabel('Cost ($1x10^8)',linespacing=3.2)
    ax.set_ylabel('\ndT_int',linespacing=3.1)
    ax.set_zlabel('\nAlpha',linespacing=3.1)
    plt.show()

# Create a zoomed graph for the PERCS
def graph2_2_zoom(opt_list,nd_list,maxes):
    print("Zoomed-in PERCS Graph")
    x_min = 9.9e99
    y_min = 9.9e99
    z_min = 9.9e99
    # Find the minimums
    for opt in opt_list:
        if opt.cost_2 < x_min: x_min = opt.cost_2
        if opt.dT_int < y_min: y_min = opt.dT_int
        if opt.alpha < z_min: z_min = opt.alpha
    fig = plt.figure(figsize=(6,5))
    ax = fig.gca(projection='3d')
    ax.set_xlim3d(x_min,maxes[0]) # cost_2
    ax.set_ylim3d(y_min,maxes[1]) # dT_int
    ax.set_zlim3d(z_min,maxes[2]) # alpha
    for i in range(len(opt_list)):
        non_dom = False
        for v in nd_list:
            if i == v:
                ax.scatter(opt_list[i].cost_2,opt_list[i].dT_int,opt_list[i].alpha,
                            c='w',edgecolors='r',marker='o')
                non_dom = True
        if non_dom == False:
            ax.scatter(opt_list[i].cost_2,opt_list[i].dT_int,opt_list[i].alpha,
                        c='w',edgecolors='k',marker='o',alpha=0.5)
    ax.view_init(azim=-110,elev=45)
    ax.set_xlabel('Cost ($1x10^8)',linespacing=3.2)
    ax.set_ylabel('\ndT_int',linespacing=3.1)
    ax.set_zlabel('\nAlpha',linespacing=3.1)
    plt.show()

# Create a graph for the PCS
def graph3_2(opt_list,nd_list):
    print("PCS Graph")
    x_max = 0.0
    y_max = 0.0
    # Find the maximums
    for opt in opt_list:
        if opt.eff > x_max: x_max = opt.eff
        if opt.cost_3 > y_max: y_max = opt.cost_3
    x_max = x_max + 0.005
    y_max = y_max + 1
    fig = plt.figure(figsize=(6,4))
    ax = fig.gca()
    plt.xlim([0.33,x_max]) # eff
    plt.ylim([0,y_max]) # cost_3
    ax.set_xticks(np.arange(0.33,x_max,0.005))
    ax.set_yticks(np.arange(0,y_max,2))
    for i in range(len(opt_list)):
        non_dom = False
        for v in nd_list:
            if i == v:
                plt.scatter(opt_list[i].eff,opt_list[i].cost_3,s=12,
                            c=opt_list[i].color,edgecolors=opt_list[i].color)
                non_dom = True
        if non_dom == False:
            plt.scatter(opt_list[i].eff,opt_list[i].cost_3,s=12,
                        c='w',edgecolors=opt_list[i].color)
    plt.grid()
    plt.xlabel('Efficiency')
    plt.ylabel('Cost ($1x10^9)')
    plt.show()
