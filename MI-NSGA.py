"""
Created on Mon Nov 05 03:52:36 2018
@author: Paul
"""

from Class2_Pop import *
from Func import *


"""""""""   MI-NSGA Code   """""""""

# Copy Excel template and create a "Results" document for this test run
wb = Book(r'C:\Users\pwild\Research\OPTIONS 3.0\Template.50.xlsm')
wb.save(r'C:\Users\pwild\Research\OPTIONS 3.0\Results.xlsm')
#wb = Book(r'C:\Users\pwild\Research\OPTIONS 3.0\Results.xlsm') #**Use if repopulating from Excel**#

# Start the Clock
tot_time = 0
time0 = time.perf_counter()
time1 = 0

# Initialize the Population
Pop = Population(80,40) #** First value should be divisible by 4 **#
Pop.Initiate(wb)
#Pop.Populate(wb,15,40) #**Use if repopulating from Excel**#

# Specify which iterations to create graphs and record current population data
graph_it = np.array((0,1,2,3,4,5,6,8,10,12,15,20,25,30,35,40,45,50))
tot_it = graph_it[np.size(graph_it)-1]
# Specify which iterations to perform Specialized Breeding step
sp_breed = np.array((1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,
                     28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50))
sb_type = 1 # Indicator for Automatic/Manual
dub_check = False # True if DM specifies desire to be offered chance of using SB next iteration
# Specify which iterations to enact "Tough Love" protocol during Competition step
tf_love = np.array((19,24,29,34,39,44,49)) #** Cannot happen on iteration 1 **#
tl_yes = False # Indicator for TL in current iteration

# Pause the clock for graphing
time1 = time.perf_counter()
dt = time1-time0
tot_time = tot_time + dt
#tot_time = wb.sheets['15'].range(2,2).value #**Use if repopulating from Excel**#

# Graph the initial population
Graph_Data(wb,Pop,0,dt,tot_time)
# Restart the Clock
time0 = time.perf_counter()

# Continue iterations until tot_it is met
for num_it in range(1,tot_it+1):
#for num_it in range(16,tot_it+1): #**Use if repopulating from Excel**#
    
    # Declare the start of each iteration
    m = time.localtime()
    if m[3]<=12:
        hr = m[3]
        if m[3]==0: hr = 12
        ap = "AM"
    if m[3]>12: 
        hr = m[3]-12
        ap = "PM"
    print ("Start Iteration",num_it,"  /  Time =", hr,":","%02.0f"%m[4],":","%02.0f"%m[5],ap)
    
    ############################################################
    """""""""  SELECTION :  Determine the Mating Pool  """""""""
    ############################################################
    
    #---------------------------
    """ Specialized Breeding """
    #---------------------------
    # Restart the AOP time for the new iteration
    dt_AOP = 0. # sec
    # Find out if this iteration will be graphed
    g_yes = False
    for g in graph_it:
        if num_it == g:
            g_yes = True
    # Find out if this iteration should include specialized breeding
    sb_yes = False
    for sb in sp_breed:
        if num_it == sb:
            sb_yes = True
    # If yes to specialized breeding and graphing, then proceed...
    # Or if yes to specialized breeding and the DM wants to do it again this iteration...
    if ((sb_yes == True) and (g_yes == True)) or ((sb_yes== True) and (dub_check==True)):
        # Start the AOP timer (so it can be subtracted from the dt for the iteration)
        t0_AOP = time.perf_counter()
        while True:
            # Ask if the DM wants to do Auto or Manual Specialized Breeding
            print("Which type of Specialized Breeding would you like to perform?")
            print("1 - Automatic // 2 - Manual")
            sb_type = eval(input())
            # Call the autogen_k_lists() or generate_k_lists() function to create
            # mating pools for specialized breeding step
            if sb_type == 1:
                dt_AOP = time.perf_counter() - t0_AOP
                k_list1,k_list2,k_list3 = autogen_k_lists(Pop,15)
                break
            elif sb_type == 2:
                k_list1,k_list2,k_list3 = generate_k_lists(Pop)
                # Ask if the DM wants to have the chance to perform Manual next iteration
                dc = input("Do you want to be prompted a choice for Auto or Manual next iteration? [y/n]: ")
                if dc=='y':
                    dub_check = True
                else:
                    dub_check = False
                break
    elif sb_yes == True:
        # Call the autogen_k_lists() function to create mating pools for specialized breeding step
        k_list1,k_list2,k_list3 = autogen_k_lists(Pop,15)
    else:
        #----------------------
        """ Normal Breeding """
        #----------------------
        # Initialize the Mating Pool
        mating_pool = np.empty(0,dtype=int)
        # Run through Selection 2x to fill the Mating Pool
        for o in range(2):
            # Meant to randomly pair the Options into 2's
            rand_pair = np.arange(np.size(Pop.opt_list))
            np.random.shuffle(rand_pair)
            # Rotate through the paired Options and undergo Selection
            for k in range(0,np.size(Pop.opt_list),2):
                # Compare fmm-values, and add winner to mating pool
                if Pop.opt_list[rand_pair[k]].fmm_o < Pop.opt_list[rand_pair[k+1]].fmm_o:
                    mating_pool = np.append(mating_pool,rand_pair[k])
                if Pop.opt_list[rand_pair[k]].fmm_o > Pop.opt_list[rand_pair[k+1]].fmm_o:
                    mating_pool = np.append(mating_pool,rand_pair[k+1])
                # No dominance by either Option
                if Pop.opt_list[rand_pair[k]].fmm_o == Pop.opt_list[rand_pair[k+1]].fmm_o:
                    if random.randint(0,2) == 0:
                        mating_pool = np.append(mating_pool,rand_pair[k])
                    else:
                        mating_pool = np.append(mating_pool,rand_pair[k+1])
        # Shuffle the Mating Pool... Make sure no duplicates are next to each other
        good = False
        while(good==False):
            for y in range(0,np.size(Pop.opt_list)-1,2):
                if mating_pool[y] == mating_pool[y+1]:
                    np.random.shuffle(mating_pool)
                    break
                if y == (np.size(Pop.opt_list)-2):
                    good = True
    
    ###########################################################
    """""""""  REPRODUCTION :  Double the Population  """""""""
    ###########################################################
    
    #---------------------------
    """ Specialized Breeding """
    #---------------------------
    # If yes to specialized breeding, then proceed...
    if sb_yes == True:
        # Specify how many specialized offspring to create
        num_gen_off = Pop.num_opt
        if sb_type == 2:
            num_gen_off = eval(input("How many specialized offspring do you want to create? [even int] "))
            dt_AOP = time.perf_counter() - t0_AOP
            # Return sb_type to 1, so that next time it won't ask if it's automaticly autogen_k_lists()
            sb_type = 1
        # Create empty matrices for specialized offspring
        n_rows = num_gen_off
        n_x1 = len(Pop.opt_list[0].x1)
        n_x2 = len(Pop.opt_list[0].x2)
        n_x3 = len(Pop.opt_list[0].x3)
        n_y = len(Pop.opt_list[0].y)
        offspring_x1s = np.zeros((n_rows,n_x1))
        offspring_x2s = np.zeros((n_rows,n_x2))
        offspring_x3s = np.zeros((n_rows,n_x3))
        offspring_ys = np.zeros((n_rows,n_y))
        # Create x-array combos and add to matrices
        for i in range(0,num_gen_off,2):
            # Reference random members of the k_lists
            i1_k1 = random.randint(0,len(k_list1))
            i2_k1 = random.randint(0,len(k_list1))
            while i1_k1 == i2_k1:
                i2_k1 = random.randint(0,len(k_list1))
            i1_k2 = random.randint(0,len(k_list2))
            i2_k2 = random.randint(0,len(k_list2))
            while i1_k2 == i2_k2:
                i2_k2 = random.randint(0,len(k_list2))
            i1_k3 = random.randint(0,len(k_list3))
            i2_k3 = random.randint(0,len(k_list3))
            while i1_k3 == i2_k3:
                i2_k3 = random.randint(0,len(k_list3))
            # Create empty x-arrays for both offspring
            o1_x1_k = np.zeros(n_x1)
            o1_x2_k = np.zeros(n_x2)
            o1_x3_k = np.zeros(n_x3)
            o1_y_k = np.zeros(n_y)
            o2_x1_k = np.zeros(n_x1)
            o2_x2_k = np.zeros(n_x2)
            o2_x3_k = np.zeros(n_x3)
            o2_y_k = np.zeros(n_y)
            # Rotate through the x1 parameters
            for j in range(n_x1):
                # Randomly assign parameters from parents to offspring
                if random.randint(0,2) == 0:
                    o1_x1_k[j] = Pop.opt_list[k_list1[i1_k1]].x1[j]
                    o2_x1_k[j] = Pop.opt_list[k_list1[i2_k1]].x1[j]
                else:
                    o1_x1_k[j] = Pop.opt_list[k_list1[i2_k1]].x1[j]
                    o2_x1_k[j] = Pop.opt_list[k_list1[i1_k1]].x1[j]
            # Rotate through the x2 parameters
            for j in range(n_x2):
                # Randomly assign parameters from parent to offspring
                if random.randint(0,2) == 0:
                    o1_x2_k[j] = Pop.opt_list[k_list2[i1_k2]].x2[j]
                    o2_x2_k[j] = Pop.opt_list[k_list2[i2_k2]].x2[j]
                else:
                    o1_x2_k[j] = Pop.opt_list[k_list2[i2_k2]].x2[j]
                    o2_x2_k[j] = Pop.opt_list[k_list2[i1_k2]].x2[j]
            # Rotate through the x3 parameters
            for j in range(n_x3):
                # Randomly assign parameters from parent to offspring
                if random.randint(0,2) == 0:
                    o1_x3_k[j] = Pop.opt_list[k_list3[i1_k3]].x3[j]
                    o2_x3_k[j] = Pop.opt_list[k_list3[i2_k3]].x3[j]
                else:
                    o1_x3_k[j] = Pop.opt_list[k_list3[i2_k3]].x3[j]
                    o2_x3_k[j] = Pop.opt_list[k_list3[i1_k3]].x3[j]
            # Randomly assign y parameters from parent to offspring
            if random.randint(0,2) == 0:
                o1_y_k = Pop.opt_list[k_list3[i1_k3]].y
                o2_y_k = Pop.opt_list[k_list3[i2_k3]].y
            else:
                o1_y_k = Pop.opt_list[k_list3[i2_k3]].y
                o2_y_k = Pop.opt_list[k_list3[i1_k3]].y
            # Fill offspring matrices with x-arrays
            offspring_x1s[i,:] = o1_x1_k
            offspring_x1s[i+1,:] = o2_x1_k
            offspring_x2s[i,:] = o1_x2_k
            offspring_x2s[i+1,:] = o2_x2_k
            offspring_x3s[i,:] = o1_x3_k
            offspring_x3s[i+1,:] = o2_x3_k
            offspring_ys[i,:] = o1_y_k
            offspring_ys[i+1,:] = o2_y_k
    else:
        #----------------------
        """ Normal Breeding """
        #----------------------
        # Create empty matrices for x- and y-array values of offspring
        n_rows = np.size(Pop.opt_list)
        n_x1 = np.size(Pop.opt_list[0].x1)
        n_x2 = np.size(Pop.opt_list[0].x2)
        n_x3 = np.size(Pop.opt_list[0].x3)
        n_y = np.size(Pop.opt_list[0].y)
        offspring_x1s = np.zeros((n_rows,n_x1))
        offspring_x2s = np.zeros((n_rows,n_x2))
        offspring_x3s = np.zeros((n_rows,n_x3))
        offspring_ys = np.zeros((n_rows,n_y))
        # Rotate through mating pairs in the mating pool
        for k in range(0,np.size(mating_pool),2):
            # Create the empty offspring x- and y-arrays
            o1_x1 = np.zeros(len(Pop.opt_list[0].x1))
            o2_x1 = np.zeros(len(Pop.opt_list[0].x1))
            o1_x2 = np.zeros(len(Pop.opt_list[0].x2))
            o2_x2 = np.zeros(len(Pop.opt_list[0].x2))
            o1_x3 = np.zeros(len(Pop.opt_list[0].x3))
            o2_x3 = np.zeros(len(Pop.opt_list[0].x3))
            o1_y = np.zeros(len(Pop.opt_list[0].y))
            o2_y = np.zeros(len(Pop.opt_list[0].y))
            # Rotate through the Core x-parameters
            for i in range(len(Pop.opt_list[0].x1)):
                # Randomly assign parameters from parents to offspring
                if random.randint(0,2) == 0:
                    o1_x1[i] = Pop.opt_list[mating_pool[k]].x1[i]
                    o2_x1[i] = Pop.opt_list[mating_pool[k+1]].x1[i]
                else:
                    o1_x1[i] = Pop.opt_list[mating_pool[k+1]].x1[i]
                    o2_x1[i] = Pop.opt_list[mating_pool[k]].x1[i]
            # Rotate through the PERCS x-parameters
            for i in range(len(Pop.opt_list[0].x2)):
                # Randomly assign parameters from parents to offspring
                if random.randint(0,2) == 0:
                    o1_x2[i] = Pop.opt_list[mating_pool[k]].x2[i]
                    o2_x2[i] = Pop.opt_list[mating_pool[k+1]].x2[i]
                else:
                    o1_x2[i] = Pop.opt_list[mating_pool[k+1]].x2[i]
                    o2_x2[i] = Pop.opt_list[mating_pool[k]].x2[i]
            # Rotate through the PCS x-parameters
            for i in range(len(Pop.opt_list[0].x3)):
                # Randomly assign parameters from parents to offspring
                if random.randint(0,2) == 0:
                    o1_x3[i] = Pop.opt_list[mating_pool[k]].x3[i]
                    o2_x3[i] = Pop.opt_list[mating_pool[k+1]].x3[i]
                else:
                    o1_x3[i] = Pop.opt_list[mating_pool[k+1]].x3[i]
                    o2_x3[i] = Pop.opt_list[mating_pool[k]].x3[i]
            # Randomly assign y-arrays from parents to offspring
            if random.randint(0,2) == 0:
                o1_y = Pop.opt_list[mating_pool[k]].y
                o2_y = Pop.opt_list[mating_pool[k+1]].y
            else:
                o1_y = Pop.opt_list[mating_pool[k+1]].y
                o2_y = Pop.opt_list[mating_pool[k]].y
            # Add the offspring x- and y-arrays to offspring matrices
            offspring_x1s[k,:] = o1_x1
            offspring_x1s[k+1,:] = o2_x1
            offspring_x2s[k,:] = o1_x2
            offspring_x2s[k+1,:] = o2_x2
            offspring_x3s[k,:] = o1_x3
            offspring_x3s[k+1,:] = o2_x3
            offspring_ys[k,:] = o1_y
            offspring_ys[k+1,:] = o2_y
    #----------------------------------------
    """ Random Offspring after Tough Love """
    #----------------------------------------
    # If tl_yes was switched to True during the previous iteration
    if tl_yes == True:
        #--------------------------------------------
        # Generate offspring from the saved x1-arrays
        #--------------------------------------------
        if len(TLo_x1s) > 0:
            # Rotate through the saved x1-arrays
            for i in range(len(TLo_x1s)):
                # Add the x1-array
                offspring_x1s = np.vstack((offspring_x1s,TLo_x1s[i,:]))
                # Pick other offspring data to generate new x2-, x3-, and y-arrays to accompany it
                o1 = random.randint(0,n_rows)
                o2 = random.randint(0,n_rows)
                while o1 == o2:
                    o2 = random.randint(0,n_rows)
                # Rotate through x2-values to create new x2-arrays
                tlo_x2 = np.zeros(len(Pop.opt_list[0].x2))
                for j in range(len(tlo_x2)):
                    # Randomly assign x2-values from the 2 random offspring x2-arrays
                    if random.randint(0,2) == 0:
                        tlo_x2[j] = offspring_x2s[o1,j]
                    else:
                        tlo_x2[j] = offspring_x2s[o2,j]
                # Rotate through x3-values to create new x3-arrays
                tlo_x3 = np.zeros(len(Pop.opt_list[0].x3))
                for j in range(len(tlo_x3)):
                    # Randomly assign x3-values from the 2 random offspring x3-arrays
                    if random.randint(0,2) == 0:
                        tlo_x3[j] = offspring_x3s[o1,j]
                    else:
                        tlo_x3[j] = offspring_x3s[o2,j]
                # Randomly assign one of the 2 random offsprings' y-array
                tlo_y = np.zeros(len(Pop.opt_list[0].y))
                if random.randint(0,2) == 0:
                    tlo_y = offspring_ys[o1,:]
                else:
                    tlo_y = offspring_ys[o2,:]
                # Add the new x2-, x3-, and y-arrays to offspring_x2s, offspring_x3s, and offspring_ys
                offspring_x2s = np.vstack((offspring_x2s,tlo_x2))
                offspring_x3s = np.vstack((offspring_x3s,tlo_x3))
                offspring_ys = np.vstack((offspring_ys,tlo_y))
        #--------------------------------------------
        # Generate offspring from the saved x2-arrays
        #--------------------------------------------
        if len(TLo_x2s) > 0:
            # Rotate through the saved x2-arrays
            for i in range(len(TLo_x2s)):
                # Add the x2-array
                offspring_x2s = np.vstack((offspring_x2s,TLo_x2s[i,:]))
                # Pick other offspring data to generate new x1-, x3-, and y-arrays to accompany it
                o1 = random.randint(0,n_rows)
                o2 = random.randint(0,n_rows)
                while o1 == o2:
                    o2 = random.randint(0,n_rows)
                # Rotate through x1-values to create new x1-arrays
                tlo_x1 = np.zeros(len(Pop.opt_list[0].x1))
                for j in range(len(tlo_x1)):
                    # Randomly assign x1-values from the 2 random offspring x1-arrays
                    if random.randint(0,2) == 0:
                        tlo_x1[j] = offspring_x1s[o1,j]
                    else:
                        tlo_x1[j] = offspring_x1s[o2,j]
                # Rotate through x3-values to create new x3-arrays
                tlo_x3 = np.zeros(len(Pop.opt_list[0].x3))
                for j in range(len(tlo_x3)):
                    # Randomly assign x3-values from the 2 random offspring x3-arrays
                    if random.randint(0,2) == 0:
                        tlo_x3[j] = offspring_x3s[o1,j]
                    else:
                        tlo_x3[j] = offspring_x3s[o2,j]
                # Randomly assign one of the 2 random offsprings' y-array
                tlo_y = np.zeros(len(Pop.opt_list[0].y))
                if random.randint(0,2) == 0:
                    tlo_y = offspring_ys[o1,:]
                else:
                    tlo_y = offspring_ys[o2,:]
                # Add the new x1-, x3-, and y-arrays to offspring_x1s, offspring_x3s, and offspring_ys
                offspring_x1s = np.vstack((offspring_x1s,tlo_x1))
                offspring_x3s = np.vstack((offspring_x3s,tlo_x3))
                offspring_ys = np.vstack((offspring_ys,tlo_y))
        #--------------------------------------------
        # Generate offspring from the saved x3-arrays
        #--------------------------------------------
        if len(TLo_x3s) > 0:
            # Rotate through the saved x3-arrays
            for i in range(len(TLo_x3s)):
                # Add the x3- and y-arrays
                offspring_x3s = np.vstack((offspring_x3s,TLo_x3s[i,:]))
                offspring_ys = np.vstack((offspring_ys,TLo_ys[i,:]))
                # Pick other offspring data to generate new x1- and x2-arrays to accompany it
                o1 = random.randint(0,n_rows)
                o2 = random.randint(0,n_rows)
                while o1 == o2:
                    o2 = random.randint(0,n_rows)
                # Rotate through x1-values to create new x1-arrays
                tlo_x1 = np.zeros(len(Pop.opt_list[0].x1))
                for j in range(len(tlo_x1)):
                    # Randomly assign x1-values from the 2 random offspring x1-arrays
                    if random.randint(0,2) == 0:
                        tlo_x1[j] = offspring_x1s[o1,j]
                    else:
                        tlo_x1[j] = offspring_x1s[o2,j]
                # Rotate through x2-values to create new x2-arrays
                tlo_x2 = np.zeros(len(Pop.opt_list[0].x2))
                for j in range(len(tlo_x2)):
                    # Randomly assign x2-values from the 2 random offspring x2-arrays
                    if random.randint(0,2) == 0:
                        tlo_x2[j] = offspring_x2s[o1,j]
                    else:
                        tlo_x2[j] = offspring_x2s[o2,j]
                # Add the new x1- and x2-arrays to offspring_x1s and offspring_x2s
                offspring_x1s = np.vstack((offspring_x1s,tlo_x1))
                offspring_x2s = np.vstack((offspring_x2s,tlo_x2))
    
    ###################################################################
    """""""""  MUTATION :  Introduce random variety into x's  """""""""
    ###################################################################
    
    # Create an empty Offspring population
    offspring_pop = Population(0,len(offspring_x1s))
    # Rotate through the rows of the offspring x- and y-matrices
    for k in range(0,len(offspring_x1s)):
        #-------------------------------
        """ Mutation of Core x-array """
        #-------------------------------
        x1_old = offspring_x1s[k,:]
        x1_new = x1_old
        # Rotate through each x1-value
        for m in range(0,np.size(x1_old)):
            # If one of the Dh_'s
            if m==2 or m==3 or m==4:
                # 40% Possibility of Mutation
                r = random.random()
                if r < 0.4:
                    dx_mut = 0.1 * random.randn()
                    x1_new[m] = np.round(x1_old[m] + dx_mut,3)
            # If other x's
            if m==0 or m==1:
                # 10% Possibility of Mutation
                r = random.random()
                if r < 0.1:
                    # Mutation for R_fuel
                    if m == 0:
                        dx_mut = 0.0005*random.randn()
                    # Mutation for H_fuel
                    if m == 1:
                        dx_mut = 1.0*random.randn()
                    # Add Mutation
                    if m == 0:
                        x1_new[m] = np.round(x1_old[m] + dx_mut,4)
                    if m == 1:
                        x1_new[m] = x1_old[m] + dx_mut
        #--------------------------------
        """ Mutation of PERCS x-array """
        #--------------------------------
        x2_old = offspring_x2s[k,:]
        x2_new = x2_old
        # Rotate through each x2-value
        for m in range(0,np.size(x2_old)):
            # 5% Possibility of Mutation
            r = random.random()
            if r < 0.05:
                # Normal distribution around original x-value
                if m == 0:
                    # Mutations for R_tank
                    dx_mut = 1.25*random.randn()
                if m == 1:
                    # Mutations for pitch
                    dx_mut = 0.1*random.randn()
                if m == 2:
                    # Mutations for D_h
                    dx_mut = 0.005*random.randn()
                if m == 3:
                    # Mutations for th
                    dx_mut = 0.001*random.randn()
                if m == 4:
                    # Mutations for Len
                    dx_mut = 2.0*random.randn()
                if m == 5:
                    # Mutations for elev
                    dx_mut = 4.0*random.randn()
                # Add mutation
                if m==0 or m==4 or m==5:
                    x2_new[m] = np.round(x2_old[m] + dx_mut,5)
                if m==1 or m==2 or m==3:
                    x2_new[m] = '%.5g'%(x2_old[m] + dx_mut)
        #------------------------------
        """ Mutation of PCS x-array """
        #------------------------------
        x3_old = offspring_x3s[k,:]
        x3_new = x3_old
        # Rotate through each x3-value
        for m in range(0,np.size(x3_old)):
            # 5% Possibility of Mutation
            r = random.random()
            if r < 0.05:
                # Normal distribution around original x-value
                if m == 2 or m == 4 or m == 7:
                    # Mutations for mf flow fractions
                    dx_mut = 0.01*random.randn()
                else:
                    # Mutations for T's and P's
                    dx_mut = random.randn()
                x3_new[m] = x3_old[m] + dx_mut
        
        # Conform new and/or mutated x- and y-values using constraints()
        x1_new,x2_new,x3_new,y_new = constraints(x1_new,x2_new,x3_new,offspring_ys[k,:])
        # Assign the new x- abd y-values back to the offspring x- and y-matrices
        offspring_x1s[k,:] = x1_new
        offspring_x2s[k,:] = x2_new
        offspring_x3s[k,:] = x3_new
        offspring_ys[k,:] = y_new
    
    # Calc all the Options in the offspring population
    offspring_pop.Breed(offspring_x1s,offspring_x2s,offspring_x3s,offspring_ys)
    # Combine the population and offspring into a "double_population"
    double_pop = np.append(Pop.opt_list,offspring_pop.opt_list)    
    
    ############################################################
    """""""""  COMPETITION :  Survival of the Fittest  """""""""
    ############################################################
    
    # Sort the double population using the Maximin function
    sorted_pop = maximin(double_pop)
    
    #------------------
    """ Kill Method """
    #------------------
    # "Survival of the Fittest" cannot allow failures, even if paired with non-dominated individual systems
    #    Forcefully kill off all tri-system options with any individual system failures
    kill_list = np.zeros(0,dtype=int)
    # Rotate through the tri-system options in sorted_pop
    for i in range(len(sorted_pop)):
        Core_failure = sorted_pop[i].failed
        PERCS_failure = sorted_pop[i].PERCS_failed
        PCS_failure = sorted_pop[i].pinch_point
        other_failure = sorted_pop[i].last_sec_penalty
        # If any sort of failure happened to one of the systems
        if Core_failure==True or PERCS_failure==True or PCS_failure==True or other_failure==True:
            # Add the sorted_list index to the kill_list
            kill_list = np.append(kill_list,i)
    print("kill_list =",kill_list)
    print("len(kill_list) =",len(kill_list))
    # Rotate backwards through the kill_list
    for j in range(len(kill_list)-1,-1,-1):
        # Delete the tri-system options from the sorted_list
        sorted_pop = np.delete(sorted_pop,kill_list[j])
    #-----------------------
    """ Tough Love stuff """
    #-----------------------
    # Find out if this iteration should include "tough love"
    tl_yes = False
    for tl in tf_love:
        if num_it == tl:
            tl_yes = True
    # If yes to tough love, then proceed...
    if tl_yes == True:
        # Start another AOP timer (so it can be subtracted from the dt for the iteration)
        t2_AOP = time.perf_counter()
        # Call the tough_love() function... which returns TL offspring data
        sorted_pop, TLo_x1s,TLo_x2s,TLo_x3s,TLo_ys = tough_love(sorted_pop)
        # Update the total for dt_AOP
        dt_AOP = dt_AOP + (time.perf_counter() - t2_AOP)
    #-----------------------------
    """ Survival of the Fittest"""
    #-----------------------------
    # Take only the best of the population to form the next population,
    #    and cut the population size down to the Pop.num_opt
    half = Pop.num_opt
    next_pop = np.array(sorted_pop[0:half])
    Pop.opt_list = next_pop
    
    ########################################################
    """""""""""""""  Graph the current data  """""""""""""""
    ########################################################
    
    # Rotate through the integers provided in 'graph_it' array
    for i in range(0,np.size(graph_it)):
        # If current iteration was specified, then call Graph_Data()
        if num_it == graph_it[i]:
            # Pause the clock for graphing
            time1 = time.perf_counter()
            dt = (time1-time0) - dt_AOP
            tot_time = tot_time + dt
            # Call the graphing function
            Graph_Data(wb,Pop,num_it,dt,tot_time)
            wb.save()
            # Restart the clock
            time0 = time.perf_counter()
    
    # REPEAT


# Calc the Options one last time, so that the final Options' data will be in the opt folders
print ("Calling Pop.calc_Options() one last time...")
global_Last_time = True
Pop.calc_Options(Pop.num_opt)
# Print the total computational run-time
print
print ('Time elapsed = ', tot_time/60.0, 'minutes')

# Rotate through the final Pareto front Options and make T-profiles, save as .png files
for i in range(len(Pop.opt_list)):
    # First make sure T_335_6 isn't too long
    if len(Pop.opt_list[i].T_335_6) > 10001:
        Pop.opt_list[i].t = Pop.opt_list[i].t[0:10001]
        Pop.opt_list[i].T_335_6 = Pop.opt_list[i].T_335_6[0:10001]
    fig = plt.figure(figsize=(7,5))
    plt.xlim([0,t_final])
    plt.ylim([480,640])
    plt.grid()
    plt.plot(Pop.opt_list[i].t,Pop.opt_list[i].T_335_6,'r-')
    name = r'C:\Users\pwild\Research\OPTIONS 3.0\Results - T_Profiles\Opt_'+repr(i+1)
    fig.savefig(name)
    plt.close(fig)
# Rotate through Options and print all T-profiles onto the same graph, save as .png file
fig = plt.figure(figsize=(7,5))
plt.xlim([0,t_final])
plt.ylim([480,640])
plt.grid()
for i in range(1,len(Pop.opt_list)):
    plt.plot(Pop.opt_list[i].t,Pop.opt_list[i].T_335_6,'r-')
plt.plot(Pop.opt_list[0].t,Pop.opt_list[0].T_335_6,'b-')
fig.savefig(r'C:\Users\pwild\Research\OPTIONS 3.0\Results - T_Profiles\All')
plt.close(fig)



