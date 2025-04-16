#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 12 19:33:17 2025

@author: katyrr
"""

"""
Created on Fri Feb  7 10:54:23 2025

@author: katyrr

HOW TO USE:
    
    - This code file is stored in /Code/Executables
    - Modules "functions.py" and "structs.py" are also stored in this directory.
    
    - Must have the following file structure set up before running (e.g. for 207Pb):
        /Code/Pb207/config_file.txt     (file path to config file)
        /Code/Pb207/Inputs/             (.DAT files will be generated here)
        /Code/Pb207/Outputs/            (.OUT files will be generated here)
        /Code/Pb207/Scripts/            (bash scripts for running codes will be 
                                         generated here)
        /Code/Pb207/Run/Batch1/         (this is where the first batch of data 
                                         points will be run through each program)
        /Code/Pb207/Run/Batch2/         (where the second batch of data points 
                                         is run through the programs in parallel)
        /Code/Pb207/Run/Batch3/         (etc. up to batch 8)
    
    - Before running, check that config_file_path is set correctly in [1. SET UP] below.
    
    - Run from terminal in working directory /Code/Executables using command: python3 main.py
    

"""

#%% 
""" 1. SET UP 
    
- Import modules.
- Set figure resolution.
- Create timers (one to time the whole program, and one to time small sections).
- State the nucleus being studied (for use in navigating file directories).

"""


import numpy as np                                   # for np.arange(start, end, step)
import math                                          # for ceil(num)
import matplotlib.pyplot as plt                      # for plotting graphs

import functions as fn                               # my own module file of functions
import structs as st                                 # my own module file of structs (classes, and read-only dicts)
import optuna_functions as opfn

import optuna                                        # for parameter optimisation
import subprocess                               # for calling shell scripts to run 
from matplotlib.ticker import FuncFormatter     # for formatting axis ticks
import matplotlib.tri as tri                    # for manual triangulation before drawing a contour plot

                        # my own module file of structs (classes, and read-only dicts)


import optuna                                        # for parameter optimisation

plt.rcParams['figure.dpi'] = 150
timer = st.Timer()
sub_timer = st.Timer()

nucleus = "Au179" #  "Pb207" #  "Fixing_Discontinuities_Pb207" # ##!!!
verbose = True # whether to print a lot of info, or just the essentials


#%% 2. READ CONFIG FILE 

'''
- Ignores empty lines, and lines beginning with * (to mark a comment).
- Check that the format of each line is correct (var_name value),
  else raise a ValueError.
  
- Saves deformation input, converting range inputs into arrays. 
- Raises a RuntimeError if: 
    deformation is input more than once, 
    or a range of E2PLUS values are input with a non-constant deformation,
    or both eps and gamma are input as a linear range (rather than a mesh).
- Raises a ValueError if any eps = 0.0, as this is not allowed in asyrmo.

- Converts all other inputs to the relevant type.
- All inputs are saved as name-value pairs in a dictionary.

'''


print("running for nucleus: " + nucleus)
timer.start()

lines = fn.read_file("../"+ nucleus +"/config.txt")

inputs = {}                                                                     
data_points = {}
experimental = {}

for l in range(len(lines)):                       
    
    # remove comments and blank lines, and check formatting  
    line = lines[l].strip()
                                
    if line == "" : continue                  
    if line[0] == "*" : continue              
    
    split_string = line.split(" ")  # split into name and value
    split_string = fn.remove_inline_comments(split_string)

    fn.check_line_format(split_string, line, l)         
    
    # save deformation inputs
    if (split_string[0]=="eps" or split_string[0]=="gamma" or split_string[0]=="single" or split_string[0]=="mesh"):
        
        if ("deformation_input" in inputs):
            raise RuntimeError("Deformation has already been input: " + inputs["deformation_input"])
        
        inputs["deformation_input"] = split_string[0]
        
        if split_string[0]=="mesh":
            data_points["eps"], data_points["gamma_degrees"] = fn.arrange_mesh(split_string[1].split(","))
            
        else: 
            data_points["eps"], data_points["gamma_degrees"] = fn.arrange_data_line(split_string)
            
        if split_string[0]=="eps":
            inputs["step"] = fn.get_range_step(split_string[1])
        elif split_string[0] == "gamma":
            inputs["step"] = fn.get_range_step(split_string[2])
                   
       
    # save E2PLUS input
    elif split_string[0]=="e2plus":
        
        if "eps" not in data_points:
            raise RuntimeError("missing deformation input (or perhaps deformation was input below e2plus?)")
            
        if split_string[1]=="0":
            # if e2plus has been input with value = "0", then it will later be 
            # calculated dynamically based on the deformation of each data point.
            data_points["e2plus"] = [0]*len(data_points["eps"])
        
        else:
            data_points["e2plus"], _ = fn.range_to_list(split_string[1])
            
            if (len(data_points["e2plus"])>1 
                and not(inputs["deformation_input"] == "single")):
                
                raise RuntimeError("Testing a range of e2plus is only " +
                               "supported for a single deformation input.")
            elif (len(data_points["e2plus"])==1
                  and len(data_points["eps"])>1):
                data_points["e2plus"] = [data_points["e2plus"][0]]*len(data_points["eps"])
                  
            else: 
                data_points["eps"] = data_points["eps"] * len(data_points["e2plus"])
                data_points["gamma_degrees"] = data_points["gamma_degrees"] * len(data_points["e2plus"])
                    
    # save other inputs
    elif split_string[0][:3]=="jp_":
        experimental[split_string[0]] = [float(n) for n in split_string[1].split(',')]
          
    elif split_string[0] in st.get_variable_list("int"):
        inputs[split_string[0]] = int(split_string[1])
        
    elif (split_string[0] in st.get_variable_list("float")
          or split_string[0][:3] in st.get_variable_list("float")):
        experimental[split_string[0]] = float(split_string[1])
                                      
    elif split_string[0] in st.get_variable_list("bool"):                                    
        inputs[split_string[0]] = bool(int(split_string[1]))
        
    else: # string
        inputs[split_string[0]] = split_string[1]


    if verbose: print(split_string[0] + ": \t" + split_string[1])
    


# check that there are no eps=0 values to test - this will cause the code to hang.
if 0.0 in data_points["eps"]:
    raise ValueError("Cannot test at eps=0.0.")


print("\n********** Finished reading lines: "+str(l)+ " **********\n")

# deallocate variables that are no longer needed
del [l, line, lines, split_string]



#%% 3. PROCESS INPUTS 
'''
- Convert gamma points to radians and save separately.
- Convert input fractional spins to floats and save separately.
- convert the nantj, noutj, ipout inputs to the correct format

- Using input A and Z, work out which is odd, to determine the nneupr input.
- Halve and ceiling for the overall index of the fermi level orbital.
- Generate the orbitals input for gampn (e.g. orbitals_string = "+4 19 20 21 22").

- Raises ValueError if an even-mass nucleus is input, or if A≠Z+N.

'''

# additionally save gamma in radians
data_points["gamma_radians"] = [n*np.pi/180 for n in data_points["gamma_degrees"]]


# add float versions of the spin parameters
if "gs_spin" in inputs:
    experimental["gs_spin_float"] = fn.spin_string_to_float(inputs["gs_spin"])

if "x1_spin" in inputs:
    experimental["x1_spin_float"] = fn.spin_string_to_float(inputs["x1_spin"])

if "x2_spin" in inputs:
    experimental["x2_spin_float"] = fn.spin_string_to_float(inputs["x2_spin"])

if "x3_spin" in inputs:
    experimental["x3_spin_float"] = fn.spin_string_to_float(inputs["x3_spin"])


# convert the nantj, noutj, ipout inputs to the correct format
inputs["nantj"] = inputs["nantj"].replace(",", " ")
inputs["noutj"] = inputs["noutj"].replace(",", " ")
inputs["ipout"] = inputs["ipout"].replace(",", " ")

# determine nneupr and calculate fermi level
inputs["N"] = inputs["A"]-inputs["Z"]


if inputs["A"]%2 == 0:
    raise ValueError("Input nucleus is even-A. Only odd-mass nuclei accepted.")
elif inputs["Z"]%2 == 0: 
    inputs["nneupr"] = "-1" 
    inputs["fermi_level"] = math.ceil(inputs["N"]/2)
    print("Calculating for odd NEUTRONS...")                 
elif inputs["N"]%2 == 0:
    inputs["nneupr"] = "1"
    inputs["fermi_level"] = math.ceil(inputs["Z"]/2)
    print("Calculating for odd PROTONS...")
else:
    raise ValueError("A ≠ N + Z; check inputs.")


# Generate a string that contains the number of orbitals and their indices, 
# in the correct format for input to gampn.
inputs["initial_orbitals"] = fn.write_orbitals(inputs["fermi_level"]//2, inputs["num_orbs"], inputs["par"])
# inputs["initial_orbitals"] = "-15 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38" #fn.write_orbitals(28, inputs["num_orbs"], inputs["par"]) #!!! hard coded version of the above




#%%

data_point = {}

t = 0

inputs["current_e2plus"] = data_points["e2plus"][t]

data_point["current_eps"] = data_points["eps"][t]
data_point["current_gamma"] = data_points["gamma_degrees"][t]

data_point["file_tag"] = "e%.3f_g%.1f_p%.3f_%s" % (data_point["current_eps"],  data_point["current_gamma"],
                                    inputs["current_e2plus"], inputs["nucleus"])

data_point["current_f002"] = "f002_"+data_point["file_tag"]+".dat"
data_point["current_f016"] = "f016_"+data_point["file_tag"]+".dat"
data_point["current_f017"] = "f017_"+data_point["file_tag"]+".dat"
data_point["current_f018"] = "f018_"+data_point["file_tag"]+".dat"


# test_output = opfn.run_model(inputs, data_point, verbose, nucleus)



#%%

trial = optuna.Trial

study = optuna.create_study(direction="minimize")
study.optimize(opfn.objective(trial, inputs, experimental), n_trials=5)


print("Best loss = ", study.best_value)
print("Best inputs: \n", study.best_params)

# note how long it took
timer.stop()
print("\n****************************************************************************************")
print("total runtime = %.2f seconds" % (timer.get_lapsed_time()))
print("****************************************************************************************\n")
