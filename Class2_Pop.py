"""
Created on Mon Nov 05 03:52:36 2018
@author: Paul
"""

### Boiler-Plate ###
import matplotlib.pylab as plt
import numpy as np
import scipy as sp
from numpy import random
import time
import threading
from queue import Queue
import os
import subprocess
from subprocess import Popen
import csv
import shutil

from Class0_Opt import *
from Class0_SS import *


q = Queue()

###############################################################################
"""""""""   Tri-System Population Class   """""""""############################
###############################################################################

global global_last_time
global_last_time = False
""" Change this value when changed in restart .i files """
global t_final
t_final = 10000 # seconds
global ss_fail_penalty
ss_fail_penalty = 700

class Population:
    """
    Inputs:
        init_num_opt = Initial population size of tri-system options
        num_opt = Population size after iteration #1
    Parameters:
        opt_list = Array of tri-system options (Option class)
        PCS_pop_list = Array of seed population of PCS options (PCS_Option class)
        PCS_keep_list = Array of best PCS seed options (PCS_Option class)
    Functions:
        Initiate() - Populates opt_list with initial tri-system options
        Breed() - (for Offspring populations) Populates opt_list using matrices of parent data
        grab() - Used to grab Excel data, if code is stopped and DM wants to pick up where it left off
        Populate() - Used after grab() to populate opt_list w/ Excel data
        RELAP_Job() - Run RELAP5 and PyPost files for specific tri-system option of opt_list
        threader() - Tells each thread to keep calling RELAP_Job() until queue is empty
        calc_Options() - Open threads, have them call threader(), then join them
        final_Option_calcs() - Run each tri-system option's remaining calculations
    """
    def __init__(self,init_num_Options,num_Options):
        self.init_num_opt = init_num_Options
        self.num_opt = num_Options
        self.opt_list = np.empty(0)
    
    """"""""" Initiate Population """""""""
    
    # Population will not initiate until this function is called
    def Initiate(self,wb):
        """
        Inputs:
            wb = Output Excel workbook to paste the PCS seed population data
        Functions:
            gen_rand_PCS_xy() - Generates random x,y-values for PCS superstructure seed options
            gen_rand_xy() - Generates random x,y-values for tri-system options
        """
        def gen_rand_PCS_xy(a):
            """
            Inputs:
                a = Integer representing a PCS superstructure configuration
            """
            # PCS x-parameters
            x_t3 = np.zeros(9)
            x_t3[0] = random.uniform(295.0,307.7) # mchx.Tout
            x_t3[1] = random.uniform(8.0,25.0)    # t1a.Pout
            x_t3[2] = random.uniform(0.05,0.19)   # mf_t1a
            x_t3[3] = random.uniform(6.0,20.0)    # t1b.Pout
            x_t3[4] = random.uniform(0.05,0.19)   # mf_t1b
            x_t3[5] = random.uniform(4.0,16.0)    # t1c.Pout
            x_t3[6] = random.uniform(3.0,13.0)    # t2a.Pout
            x_t3[7] = random.uniform(0.05,0.19)   # mf_t2a
            x_t3[8] = random.uniform(2.0,11.0)    # t2b.Pout
            # PCS y-parameters
            if a==1: y_t = np.array((0,1,0,0,0,0)) # red
            if a==2: y_t = np.array((0,0,0,0,0,0)) # firebrick
            if a==3: y_t = np.array((1,1,1,0,1,0)) # darkgreen
            if a==4: y_t = np.array((1,1,1,0,0,1)) # purple
            if a==5: y_t = np.array((1,0,1,0,0,1)) # deeppink
            if a==6: y_t = np.array((1,1,1,1,0,0)) # blue
            if a==7: y_t = np.array((1,0,1,1,0,0)) # cyan
            if a==8: y_t = np.array((1,0,0,0,0,0)) # orange
            if a==9: y_t = np.array((1,1,0,0,0,0)) # yellow
            # Make sure parameters conform with PCS_constraints()
            x_t3,y_t = PCS_constraints(x_t3,y_t)
            # Return x_t3 and y_t
            return x_t3,y_t
        
        ################################################################
        """ Initial Creation, Competition & Crossover of PCS Options """
        ################################################################
        t_one = time.perf_counter()
        # Calculate an appropriate init_PCS_copies value
        num_PCS_keep = int(self.init_num_opt / 2)
        init_PCS_copies = int(num_PCS_keep / 2)
        """ Create mini population list of PCS SS Options """
        # Create first PCS_Option and add to PCS_pop_list
        x_t3,y_t = gen_rand_PCS_xy(1)
        self.PCS_pop_list = np.array((PCS_Option(x_t3,y_t)))
        # Finish initiating PCS_Options and add to PCS_pop_list
        for i_a in range(1,10):
            if i_a == 1: i_b_start = 2
            if i_a != 1: i_b_start = 1
            for i_b in range(i_b_start,init_PCS_copies+1):
                x_t3,y_t = gen_rand_PCS_xy(i_a)
                self.PCS_pop_list = np.append(self.PCS_pop_list,PCS_Option(x_t3,y_t))
        # Calc fmm-values and sort the PCS_pop_list
        self.PCS_pop_list = PCS_maximin(self.PCS_pop_list)
        # Graph the new PCS_pop_list
        dt_pcs = time.perf_counter() - t_one
        PCS_Graph_Data(wb,self.PCS_pop_list,'0_pcs',dt_pcs)
        # Cut down the PCS_pop_list to half of init_num_opt
        self.PCS_keep_list = self.PCS_pop_list[0:num_PCS_keep]
        """ Double the PCS_keep_list size by creating new x- and y-arrays """
        # Initialize matrices to hold all x- and y-data for PCSs
        pcs_x_list = np.zeros((self.init_num_opt,9))
        pcs_y_list = np.zeros((self.init_num_opt,6))
        # Rotate through the PCS_keep_list
        for i in range(len(self.PCS_keep_list)):
            # Fill the pcs lists with the x- and y-data from PCS_keep_list
            pcs_x_list[i,:] = self.PCS_keep_list[i].x
            pcs_y_list[i,:] = self.PCS_keep_list[i].y
            # Every other time, do the following
            if i%2 == 0:
                # Create arrays for fake PCS offspring to add to the lists
                o1x,o1y = np.zeros(9),np.zeros(6)
                o2x,o2y = np.zeros(9),np.zeros(6)
                # Randomly assign x-parameters from parents to offspring
                for j in range(9):
                    if random.randint(0,2) == 0:
                        o1x[j] = self.PCS_keep_list[i].x[j]
                        o2x[j] = self.PCS_keep_list[i+1].x[j]
                    else:
                        o1x[j] = self.PCS_keep_list[i+1].x[j]
                        o2x[j] = self.PCS_keep_list[i].x[j]
                # Randomly assign y-parameters from parets to offspring
                for j in range(6):
                    if random.randint(0,2) == 0:
                        o1y[j] = self.PCS_keep_list[i].y[j]
                        o2y[j] = self.PCS_keep_list[i+1].y[j]
                    else:
                        o1y[j] = self.PCS_keep_list[i+1].y[j]
                        o2y[j] = self.PCS_keep_list[i].y[j]
                # Pass offspring x- and y-arrays through PCS_constraints
                o1x,o1y = PCS_constraints(o1x,o1y)
                o2x,o2y = PCS_constraints(o2x,o2y)
                # Add offspring to the x- and y-data lists
                ip = int(i+num_PCS_keep)
                pcs_x_list[ip,:] = o1x
                pcs_x_list[ip+1,:] = o2x
                pcs_y_list[ip,:] = o1y
                pcs_y_list[ip+1,:] = o2y
        
        def gen_rand_xy(v):
            # Core x-parameters
            x_t1 = np.zeros(5)
            x_t1[0] = np.round(random.uniform(0.0120,0.0159),4)  # ft (R_tuel)
            x_t1[1] = random.uniform(10.0,14.0)                  # ft (H_fuel)
            x_t1[2] = np.round(random.uniform(1.5,3.5),3)        # ft (Dh_00)
            x_t1[3] = np.round(random.uniform(1.5,3.5),3)        # ft (Dh_12)
            x_t1[4] = np.round(random.uniform(1.5,3.5),3)        # ft (Dh_14)
            # PERCS x-parameters
            x_t2 = np.zeros(6)
            x_t2[0] = np.round(random.uniform(3.28084,16.4042),5)  # ft (R_tank)
            x_t2[2] = '%.5g'%random.uniform(3.28084e-2,1.64042e-1) # ft (D_h)
            x_t2[1] = '%.5g'%random.uniform(1.25*x_t2[2],0.65617)  # ft (pitch)
            x_t2[3] = '%.5g'%random.uniform(3.28084e-3,8.2021e-3)  # ft (th)
            x_t2[4] = np.round(random.uniform(3.28084,49.2126),5)  # ft (Len)
            x_t2[5] = np.round(random.uniform(16.4042,131.234),5)  # ft (elev)
            # PCS x-parameters
            x_t3 = pcs_x_list[v,:]
            # PCS y-parameters
            y_t = pcs_y_list[v,:]
            #Subject new optimization parameters to the constraints function
            x_t1,x_t2,x_t3,y_t = constraints(x_t1,x_t2,x_t3,y_t)
            return x_t1,x_t2,x_t3,y_t
        
        #################################################
        """ Create a Population of Tri-System Options """
        #################################################
        # Creat first Option and add to opt_list
        x_t1,x_t2,x_t3,y_t = gen_rand_xy(0)
        self.opt_list = np.array((Option(x_t1,x_t2,x_t3,y_t)))
        # Finish initiating Options and add to opt_list
        for i in range(1,self.init_num_opt):
            x_t1,x_t2,x_t3,y_t = gen_rand_xy(i)
            self.opt_list = np.append(self.opt_list,Option(x_t1,x_t2,x_t3,y_t))
        
        # Calculate the newly initiated Options, assign fmm-values, and sort
        self.calc_Options(self.init_num_opt)
        self.final_Option_calcs()
        time.sleep(5)
        self.opt_list = maximin(self.opt_list)
    
    """"""""" Breed Population """""""""
    
    # Offspring population will not initiate until this function is called and
    def Breed(self,x1_matrix,x2_matrix,x3_matrix,y_matrix):
        # Add the first Option using the top row of the x- and y-matrices
        self.opt_list = np.array((Option(x1_matrix[0,:],x2_matrix[0,:],x3_matrix[0,:],y_matrix[0,:])))
        # Rotate through k to add the remaining Options using rows from the x- and y-matrices
        for k in range(1,self.num_opt):
            self.opt_list = np.append(self.opt_list,Option(x1_matrix[k,:],x2_matrix[k,:],x3_matrix[k,:],y_matrix[k,:]))
        # Calculate the newly initiated Options
        self.calc_Options(self.num_opt)
        self.final_Option_calcs()
        time.sleep(5)
    
    """"""""" Populate from Excel """""""""
    
    # Grab optimization parameters from Excel tab
    def grab(self,wb,tab,col):
        # Core loop x-parameters
        x1 = np.zeros(5)
        for j in range(len(x1)):
            x1[j] = wb.sheets[tab].range(j+14,col).value
        # PERCS x-parameters
        x2 = np.zeros(6)
        for j in range(len(x2)):
            x2[j] = wb.sheets[tab].range(j+19,col).value
        # PCS x-parameters
        x3 = np.zeros(9)
        for j in range(len(x3)):
            x3[j] = wb.sheets[tab].range(j+25,col).value
        # PCS y-parameters
        y = np.zeros(6)
        for j in range(len(y)):
            y[j] = wb.sheets[tab].range(j+35,col).value
        return x1,x2,x3,y
    
    # Generate the population using an Excel tab filled with options
    def Populate(self,wb,tb,num_opts):
        tab = repr(tb)
        col = 3
        # Create first Option and add to opt_list
        x_t1,x_t2,x_t3,y_t = self.grab(wb,tab,col)
        self.opt_list = np.array((Option(x_t1,x_t2,x_t3,y_t)))
        # Finish initiating Options and add to opt_list
        for i in range(1,num_opts):
            col = col + 1
            x_t1,x_t2,x_t3,y_t = self.grab(wb,tab,col)
            self.opt_list = np.append(self.opt_list,Option(x_t1,x_t2,x_t3,y_t))
        # Calculate the newly initiated Options, assign fmm-values, and sort
        self.calc_Options(num_opts)
        self.final_Option_calcs()
        time.sleep(5)
        self.opt_list = maximin(self.opt_list)
    
    """"""""" Calculate Options """""""""
    
    # Read .i file, Make changes to .i file, Run .bat file, 
    def RELAP_Job(self,job,th):
        """
        Inputs:
            job = Integer representing the index in the population's opt_list
            th = Local thread variable (not used here, but useful if needed)
        """
        # Assign the job number to the Option
        self.opt_list[job].opt_ID = job
        
        ########################
        """ Read the .i file """
        ########################
        # Find the correct input file, open it, and read lines
        first_part = r"C:\Users\pwild\Research\OPTIONS 3.0\Options\opt_"
        second_part = "\\zionpwr_PERCS_SBO_6in.i"
        i_file_name = first_part + str(job+1) + second_part
        with open(i_file_name,'r') as f:
            data = f.readlines()
        
        ########################################
        """ Make changes to and save .i file """
        ########################################
        #-----------------------------
        """ Make Zion Core changes """
        #-----------------------------
        # Initiate counters
        i_line = 0
        i_card = 0
        # Grab cards, i_vals, & vals lists from opt_list[job]
        card = self.opt_list[job].cards
        i_change = self.opt_list[job].i_vals
        change = self.opt_list[job].vals
        # Loop through the lines of the file
        done = False
        for line in data:
            # Split all the words of each line
            words = line.split()
            # Make sure the line has enough words to check for card number
            if np.size(words) >= 1:
                # If the line starts with the next card number
                while words[0] == card[i_card]:
                    # Turn the replacement value into a string of characters
                    replacement = str(change[i_card])
                    # Make replacement within the line of data
                    line = line.replace(words[i_change[i_card]-1],replacement)
                    data[i_line] = line
                    # Increase i_card count
                    i_card = i_card + 1
                    # Exit loop if last card change was made
                    if i_card+1 > np.size(card):
                        done = True
                        break
            # Previous break only affects the 'while', so here's another one
            if done == True:
                break
            # Increment the line index
            i_line = i_line + 1
        #-------------------------
        """ Make PERCS changes """
        #-------------------------
        # Initiate counters
        i_line = 0
        i_card = 0
        # Grab list_card, list_i_change, list_change lists from opt_list[job]
        card = self.opt_list[job].list_card
        i_change = self.opt_list[job].list_i_change
        change = self.opt_list[job].list_change
        # Loop through the lines of the file
        done = False
        for line in data:
            # Split all the words of each line
            words = line.split()
            # Make sure the line has enough words to check for card number
            if np.size(words) >= 1:
                # If the line starts with the next card number
                while words[0] == card[i_card]:
                    # Turn the replacement value into a string of characters
                    replacement = str(change[i_card])
                    # Make replacement within the line of data
                    line = line.replace(words[i_change[i_card]-1],replacement)
                    data[i_line] = line
                    # Increase i_card count
                    i_card = i_card + 1
                    # Exit loop if last card change was made
                    if i_card+1 > np.size(card):
                        done = True
                        break
            # Previous break only affects the 'while', so here's another one
            if done == True:
                break
            # Increment the line index
            i_line = i_line + 1
        #--------------------------------
        """ Write the updated .i file """
        #--------------------------------
        with open(i_file_name,'w') as f:
            for line in data:
                f.write(line)
        
        #################################################
        """ Run Core.bat file to generate new .r file """
        #################################################
        # Find the correct batch file location and run Core.bat
        first_part = r"C:\Users\pwild\Research\OPTIONS 3.0\Options\opt_"
        batchFileLocation = first_part + str(job+1)
        batchFileFullPath = os.path.join(batchFileLocation,"Core.bat")
        p = Popen(os.path.abspath(batchFileFullPath),stdin=subprocess.PIPE,cwd=batchFileLocation)
        stdout, stderr = p.communicate()
        # Place a copy of the new .r file in the "Core PyPost Data" folder
        original = first_part + str(job+1) + "\\zionpwr_PERCS_SBO_6in.r"
        copy = first_part + str(job+1) + "\\Core PyPost Data\\zionpwr_PERCS_SBO_6in.r"
        shutil.copyfile(original,copy)
        time.sleep(20)
        
        ############################################################################
        """ Check for failure due to 'Errors detected during input processing.'  """
        ############################################################################
        # 0******** Errors detected during input processing. #
        # Find the correct text file, open it, and read lines
        second_part2 = "\\zionpwr_PERCS_SBO_6in.txt"
        txt_file_name = first_part + str(job+1) + second_part2
        with open(txt_file_name,'r') as g:
            data2 = g.readlines()
        for line in data2:
            # Split all the words of each line
            words = line.split()
            # Make sure the line has enough words to check for 0*** error
            if np.size(words) >= 1:
                # If the line starts with the 0*** error
                if words[0] == '0********':
                    print("Checkpoint 1:",self.opt_list[job].opt_ID+1,"failed!")
                    self.opt_list[job].failed = True
        # Wait for 20 seconds while it checks for failure
        time.sleep(20)
        
        ################################
        """ Run Zion Core pypost.bat """
        ################################
        if self.opt_list[job].failed == False:
            # Find the correct batch fil location and run PyPost.bat
            # PyPost.bat runs a script that grabs .r data and puts it into .csv files
            second_part3 = "\\Core PyPost Data"
            batchFileLocation3 = first_part + str(job+1) + second_part3
            batchFileFullPath3 = os.path.join(batchFileLocation3,"PyPost.bat")
            p = Popen(os.path.abspath(batchFileFullPath3),stdin=subprocess.PIPE,cwd=batchFileLocation3)
            stdout, stderr = p.communicate()
            # Wait for 20 seconds while PyPost creates the .csv files
            time.sleep(20)
            # Save the batch file location so the Option can find the .csv files
            self.opt_list[job].csvFileLocation = batchFileLocation3
        else:
            # Warm DM that there was failure in running the Core PyPost file
            print("Checkpoint 2:",self.opt_list[job].opt_ID+1,"did not enter to run Core PyPost")
        
        ##################################################
        """ Run PERCS.bat file to generate new .r file """
        ##################################################
        if self.opt_list[job].failed == False:
            # Find the correct batch file location and run PERCS.bat
            batchFileFullPath4 = os.path.join(batchFileLocation,"PERCS.bat")
            p = Popen(os.path.abspath(batchFileFullPath4),stdin=subprocess.PIPE,cwd=batchFileLocation)
            stdout, stderr = p.communicate()
            # Place a copy of the new .r file in the "PERCS PyPost Data" folder
            original2 = first_part + str(job+1) + "\\zionpwr_PERCS_SBO_6in.r"
            copy2 = first_part + str(job+1) + "\\PERCS PyPost Data\\zionpwr_PERCS_SBO_6in.r"
            shutil.copyfile(original2,copy2)
        
        ############################
        """ Run PERCS pypost.bat """
        ############################
        if self.opt_list[job].failed == False:
            # Find the correct batch fil location and run PyPost.bat
            # PyPost.bat runs a script that grabs .r data and puts it into .csv files
            second_part4 = "\\PERCS PyPost Data"
            batchFileLocation5 = first_part + str(job+1) + second_part4
            batchFileFullPath5 = os.path.join(batchFileLocation5,"PyPost.bat")
            p = Popen(os.path.abspath(batchFileFullPath5),stdin=subprocess.PIPE,cwd=batchFileLocation5)
            stdout, stderr = p.communicate()
            # Wait for 20 seconds while PyPost creates the .csv files
            time.sleep(20)
            # Save the batch file location so the Option can find the .csv files
            self.opt_list[job].csvFileLocation2 = batchFileLocation5
    
    # The threader function grabs jobs off the Queue and calls RELAP_Job to solve 
    # the job, passing it the corresponding Option from the opt_list.
    def threader(self):
        # Create a local thread variable that is specific to the thread
        th = threading.local()
        th.count = True
        # Have thread loop forever
        while th.count == True:
            # Grab a job off the queue
            job = q.get()
            # Call RELAP_Job() to run RELAP5/PyPost files for specific tri-system option
            self.RELAP_Job(job,th)
            # Mark job as done
            q.task_done()
            # If the queue is empty, then stop looping
            if q.empty() == True:
                th.count = False
    
    def calc_Options(self,pp):
        """
        Inputs:
            pp = number of threads to open
        """
        # Open all the necessary threads
        for worker in range(pp):
            t = threading.Thread(target = self.threader)
            t.daemon = True
            t.start()
        # Add jobs to the Queue, 1 for each Option
        for job in range(len(self.opt_list)):
            q.put(job)
        # Join the threads and close the Queue
        q.join()
        
        # Declare the start of each iteration
        m = time.localtime()
        if m[3]<=12:
            hr = m[3]
            if m[3]==0: hr = 12
            ap = "AM"
        if m[3]>12: 
            hr = m[3]-12
            ap = "PM"
        print ("Threads closed  /  Time =", hr,":","%02.0f"%m[4],":","%02.0f"%m[5],ap)
    
    # Do the final calcs for the Options separately from calc_Options()
    def final_Option_calcs(self):
        ##############################################
        """ Run each Option's Zion and PERCS Calcs """
        ##############################################
        # Make sure this isn't the final run (which is just for the purpose of updating the RELAP5 files)
        if global_last_time == False:
            for j in range(len(self.opt_list)):
                # Say we started Zion/PERCS calcs for each Option
                m = time.localtime()
                if m[3]<=12:
                    hr = m[3]
                    if m[3]==0: hr = 12
                    ap = "AM"
                if m[3]>12: 
                    hr = m[3]-12
                    ap = "PM"
                print(self.opt_list[j].opt_ID+1," final Zion/PERCS calcs", hr,":","%02.0f"%m[4],":","%02.0f"%m[5],ap)
                # Only do if the Core loop did not fail
                if self.opt_list[j].failed == False:
                    time.sleep(60)
                    # Perform all final calculations
                    self.opt_list[j].final_ZION_calcs()
                    self.opt_list[j].final_PERCS_calcs()
                    self.opt_list[j].Alpha_calcs()
                # In the event of a fail, just define penalized obj. func. values
                else:
                    # Zion obj. funcs penalized
                    self.opt_list[j].cost_1 = 15 # $1x10^8
                    self.opt_list[j].W_rcp = 50 # MW
                    # PERCS obj. funcs penalized
                    self.opt_list[j].cost_2 = 300. # $1x10^8
                    self.opt_list[j].dT_int = 8*10**5.0 # K*s
                    self.opt_list[j].alpha = 0.015 # frac
        # Add penalty to any option with any crazy obj. function error
        for j in range(len(self.opt_list)):
            if self.opt_list[j].W_rcp < 0.0:
                self.opt_list[j].W_rcp = 50 # MW
                self.opt_list[j].cost_1 = 15 # $1x10^8
                print(self.opt_list[j].opt_ID," had a -W_rcp!")
                self.opt_list[j].last_sec_penalty = True
            if self.opt_list[j].alpha < 0.001:
                self.opt_list[j].cost_2 = 300. # $1x10^8
                self.opt_list[j].dT_int = 8*10**5.0 # K*s
                self.opt_list[j].alpha = 0.015 # frac
                print(self.opt_list[j].opt_ID," had a low Alpha!")
                self.opt_list[j].last_sec_penalty = True
                    
        ###################################
        """ Run each Option's PCS calcs """
        ###################################
        # Make sure this isn't the final run (which is just for the purpose of updating the RELAP5 files)
        if global_last_time == False:
            print("Doing final PCS calcs")
            for j in range(len(self.opt_list)):
                self.opt_list[j].PCS_SS_calcs()
                # Add penalty to fmm_3 if PERCS was penalized
                if self.opt_list[j].PERCS_failed == True:
                    self.opt_list[j].fmm_3 = abs(self.opt_list[j].fmm_3)*10.0

