# coding: utf-8

import subprocess
import os
import sys
import glob
import os.path
import numpy as np

ECL_path='./AA222_EclipseFile' #Folder name that contains the Eclipse files, currently on desktop
cmd='eclrun eclipse AA222.data' # This is the cmd command to run eclipse
fname='AA222.RSM' # This file contains the output of the eclipse run e.g., time, production rate, etc.
DETACHED_PROCESS = 0x00000008 # Don't worry about this guy

N=1 # number of particles aka number of function evaluations per iteration
oil_price = 60 # $/stb
gas_cost = 5.3E-3 # $/MSCF
discount_rate = 0.1 # Annual discount rate

NPV_outfile='NPV.txt' # This is where the NPV output for all particles will be stored at currently
f= open(NPV_outfile,"w+")
f.close()

# Run Eclipse N number of times per iteration and evaluate N NPVs
for n in range(N):
    os.chdir(ECL_path)
    subprocess.call(cmd, creationflags=DETACHED_PROCESS) # Run Eclipse - ECLRUN console will open and will close automatically - it's kinda annoying but will worry about getting rid of it later
    data = np.loadtxt(fname,skiprows=10)
    time = data[:,0] # Get time steps in days
    oil_rate = data[:,3] # Get field oil production rate STB/day
    gas_rate = data[:,2] # Get field gas production rate MSCF/day

    # Compute NPV
    NPV=0
    for i in range(len(time)-1):
        profit = oil_price * oil_rate[i+1] - gas_cost * gas_rate[i+1] # This is the profit between time tn-1 and tn
        time_factor = (time[i+1]-time[i])/np.power((1+discount_rate),(time[i+1]/365)) # This is the time multiplier between time tn-1 and tn
        NPV += time_factor * profit # Compute NPV and the final NPV will be the output

    # Currently I am printing NPV out in a file called 'NPVs.txt', this file is at the same location as this .py script
    # This file will capture the NPV values for all particles (in this current example 2 particles as N=2)
    # The first column is the particle index and the second column is the corresponding NPV value for that particle in this iteration. The two values are separated by a comma 
    os.chdir('..')
    try:        
        output = open(NPV_outfile, 'a')
        output.write(str(n+1)+','+str(NPV)+'\n')
        output.close()
    except:
        print("Unable to append to " + NPV_outfile)

