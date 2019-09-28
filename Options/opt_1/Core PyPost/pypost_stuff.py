import csv
import os

#######################
#  Perliminary Stuff  #
#######################
# Set the units to SI
setUseSIUnits()
# Get the current directory
dir_path = os.path.dirname(os.path.realpath('__file__'))
first_part = str(dir_path)
# Loop through directory and change '\' to '/'
for k in range(len(first_part)):
    if first_part[k] == '\\':
        first_part = first_part[:k] + '/' + first_part[k+1:]
# Finish full path to .r file
file_path = first_part + "/zionpwr_PERCS_SBO_6in.r"
RELAP.openPlotFile(file_path,1)

####################
#  All tempf Data  #
####################
# Create and open tempf_data.csv file
tempf_file_name = first_part + "/tempf_data.csv"
f1 = open(tempf_file_name,'wb')
writer1 = csv.writer(f1)
# Grab the T_data and write to .csv file
T_106 = RELAP.getData(1,'tempf-106010000')
i = 0
for row in T_106:
    if i == len(T_106)-1:
        writer1.writerow([row[1]])
    i = i + 1
T_110 = RELAP.getData(1,'tempf-110010000')
i = 0
for row in T_110:
    if i == len(T_110)-1:
        writer1.writerow([row[1]])
    i = i + 1
f1.close()
################
#  All p Data  #
################
# Create and open p_data.csv file
p_file_name = first_part + "/p_data.csv"
f2 = open(p_file_name,'wb')
writer2 = csv.writer(f2)
# Grab the p_data and write to .csv file
P_106 = RELAP.getData(1,'p-106010000')
i = 0
for row in P_106:
    if i == len(P_106)-1:
        writer2.writerow([row[1]])
    i = i + 1
P_110 = RELAP.getData(1,'p-110010000')
i = 0
for row in P_110:
    if i == len(P_110)-1:
        writer2.writerow([row[1]])
    i = i + 1
num = 10000
for k in range(6):
    var_name = 'p-3350'+str(num)
    p_data = RELAP.getData(1,var_name)
    i = 0
    for row in p_data:
        if i == len(p_data)-1:
            writer2.writerow([row[1]])
        i = i + 1
    num = num + 10000
P_114 = RELAP.getData(1,'p-114010000')
i = 0
for row in P_114:
    if i == len(P_114)-1:
        writer2.writerow([row[1]])
    i = i + 1
f2.close()

#####################
#  All mflowj Data  #
#####################
# Create and open m_dot_data.csv file
m_file_name = first_part + "/mflowj_data.csv"
f3 = open(m_file_name,'wb')
writer3 = csv.writer(f3)
# Grab the m_dot_data and write to .csv file
nums = [100010000,400010000,600010000,200010000]
for num in nums:
    var_name = 'mflowj-'+str(num)
    m_data = RELAP.getData(1,var_name)
    i = 0
    for row in m_data:
        if i >= len(m_data)-7:
            writer3.writerow([row[1]])
        i = i + 1
j = 0
m_data = RELAP.getData(1,'mflowj-335010000')
for row in m_data:
    if j == len(m_data)-1:
        writer3.writerow([row[1]])
    j = j + 1
f3.close()
####################
#  All hvmix Data  #
####################
# Create and open m_dot_data.csv file
h_file_name = first_part + "/hvmix_data.csv"
f4 = open(h_file_name,'wb')
writer4 = csv.writer(f4)
# Grab the hvmix_data and write to .csv file
nums = [106010000,110010000,335010000,112050000,114010000,412050000,414010000,612050000,614010000,212050000,214010000]
for num in nums:
    var_name = 'hvmix-'+str(num)
    h_data = RELAP.getData(1,var_name)
    i = 0
    for row in h_data:
        if i == len(h_data)-1:
            writer4.writerow([row[1]])
        i = i + 1
f4.close()
#####################
#  All httemp Data  #
#####################
# Create and open httemp_data.csv file
T_file_name = first_part + "/httemp_data.csv"
f5 = open(T_file_name,'wb')
writer5 = csv.writer(f5)
# Grab the httemp_data and write to .csv file
num = 101
for k in range(6):
    var_name = 'httemp-336000'+str(num)
    T_data = RELAP.getData(1,var_name)
    i = 0
    for row in T_data:
        if i == len(T_data)-1:
            writer5.writerow([row[1]])
        i = i + 1
    num = num + 100
num = 117
for k in range(6):
    var_name = 'httemp-336000'+str(num)
    T_data = RELAP.getData(1,var_name)
    i = 0
    for row in T_data:
        if i == len(T_data)-1:
            writer5.writerow([row[1]])
        i = i + 1
    num = num + 100
f5.close()
####################
#  All quale Data  #
####################
# Create and open quale_data.csv file
x_file_name = first_part + "/quale_data.csv"
f6 = open(x_file_name,'wb')
writer6 = csv.writer(f6)
# Grab the quale_data and write to .csv file
num = 10000
for k in range(6):
    var_name = 'quale-3350'+str(num)
    x_data = RELAP.getData(1,var_name)
    i = 0
    for row in x_data:
        if i == len(x_data)-1:
            writer6.writerow([row[1]])
        i = i + 1
    num = num + 10000
f6.close()

######################
#  Close RELAP file  #
######################
RELAP.closeAll()
