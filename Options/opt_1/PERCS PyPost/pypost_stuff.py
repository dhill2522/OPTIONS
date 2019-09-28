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
###################
# All httemp data #
###################
# Grab the T_data needed for obj. function
httemp1 = RELAP.getData(1,'tempf-335060000')
# Open and write T_data.csv file
T_file_name = first_part + "/T_data.csv"
f = open(T_file_name,'wb')
writer = csv.writer(f)
writer.writerows(httemp1)
f.close()
##################
# All alpha data #
##################
# Grab the a_frac data needed for obj. function
a_data = RELAP.getData(1,'cntrlvar-704')
# Open and write it to alpha_data.csv file
a_file_name = first_part + "/Alpha_data.csv"
g = open(a_file_name,'wb')
writer_2 = csv.writer(g)
i = 0
for row in a_data:
    if i == len(a_data)-1:
        writer_2.writerow([row[1]])
    i = i + 1
num = 713
for k in range(99):
    var_name = 'cntrlvar-'+str(num)
    a_data = RELAP.getData(1,var_name)
    i = 0
    for row in a_data:
        if i == len(a_data)-1:
            writer_2.writerow([row[1]])
        i = i + 1
    num = num + 7
g.close()
####################
# Close RELAP file #
####################
RELAP.closeAll()
