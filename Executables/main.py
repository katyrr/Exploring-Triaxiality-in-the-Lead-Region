#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  7 10:54:23 2025

@author: Katy Robson

HOW TO USE:
    
    - This code file is stored in /Code/Executables
    - Modules "functions.py", "structs.py", and "graph_plotting.py" are also stored in this directory.
    
    - Must have the following file structure set up before running (e.g. for 207Pb):
        
        /Code/Pb207/config.txt          (file path to config file)
        /Code/Pb207/Inputs/             (.DAT files will be generated here)
        /Code/Pb207/Outputs/            (.OUT files will be generated here)
        /Code/Pb207/Scripts/            (bash scripts for running codes will be 
                                         generated here)
        /Code/Pb207/Run/Batch1/         (this is where the first batch of data 
                                         points will be run through each program)
        /Code/Pb207/Run/Batch2/         (where the second batch of data points 
                                         is run through the programs in parallel)
        /Code/Pb207/Run/Batch3/         (etc. up to batch 8)
    
    - Before running, check that config.txt has been filled in with the desired input parameters.
    - Also set which graphs to plot in section 11 of this file
    
    - Run from terminal in working directory /Code/Executables using command: python3 main.py
    

"""

#%% 
""" 1. SET UP 
    
- Import modules.
- Set figure resolution.
- Create timers (one to time the whole program, and one to time small sections).
- State the nucleus being studied (for navigating file directories).

"""


import numpy as np                                  
import math                                        
import matplotlib.pyplot as plt                      

import functions as fn                              
import structs as st
import graph_plotting as gr


plt.rcParams['figure.dpi'] = 150  # set figure resolution

_timer = st.Timer()
_sub_timer = st.Timer()

nucleus = "Au179" #  "Pb207" # "Pt177" #  "Au179_test10000" # "Fixing_Discontinuities_Pb207" #!!!
verbose = False # whether to print a lot of info, or just the essentials


#%%
''' 2. READ CONFIG FILE 

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

print("\nrunning for nucleus: " + nucleus)
_timer.start()

_lines = fn.read_file("../"+ nucleus +"/config.txt")

inputs = {}                                                                     
data_points = {}
experimental = {}

for i in range(len(_lines)):                       
    
    # remove comments and blank lines, and check formatting  
    _line = _lines[i].strip()
                                
    if _line == "" : continue                  
    if _line[0] == "*" : continue              
    
    _split_string = _line.split(" ")  # split into name and value
    _split_string = fn.remove_inline_comments(_split_string, i)

    fn.check_line_format(_split_string, _line, i)         
    
    # save deformation inputs
    if (_split_string[0]=="eps" or _split_string[0]=="gamma" or _split_string[0]=="single" or _split_string[0]=="mesh"):
        
        if ("deformation_input" in inputs):
            raise RuntimeError("Deformation has already been input: " + inputs["deformation_input"])
        
        inputs["deformation_input"] = _split_string[0]
        
        if _split_string[0]=="mesh":
            data_points["eps"], data_points["gamma_degrees"] = fn.arrange_mesh(_split_string[1].split(","))
            
        else: 
            data_points["eps"], data_points["gamma_degrees"] = fn.arrange_data_line(_split_string)
            
        if _split_string[0]=="eps":
            inputs["step"] = fn.get_range_step(_split_string[1])
        elif _split_string[0] == "gamma":
            inputs["step"] = fn.get_range_step(_split_string[2])
                   
       
    # save E2PLUS input
    elif _split_string[0]=="e2plus":
        
        if "eps" not in data_points:
            raise RuntimeError("missing deformation input (or perhaps deformation was input below e2plus?)")
            
        if _split_string[1]=="0":
            # if e2plus has been input with value = "0", then it will later be 
            # calculated dynamically based on the deformation of each data point.
            data_points["e2plus"] = [0]*len(data_points["eps"])
        
        else:
            data_points["e2plus"], i = fn.range_to_list(_split_string[1])
            
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
    elif _split_string[0][:3]=="jp_":
        experimental[_split_string[0]] = [float(n) for n in _split_string[1].split(',')]
          
    elif _split_string[0] in st.get_variable_list("int"):
        inputs[_split_string[0]] = int(_split_string[1])
        
    elif (_split_string[0] in st.get_variable_list("float")
          or _split_string[0][:3] in st.get_variable_list("float")):
        experimental[_split_string[0]] = float(_split_string[1])
                                      
    elif _split_string[0] in st.get_variable_list("bool"):                                    
        inputs[_split_string[0]] = bool(int(_split_string[1]))
        
    else: # string
        inputs[_split_string[0]] = _split_string[1]


    if verbose: print(_split_string[0] + ": \t" + _split_string[1])
    


# check that there are no eps=0 values to test - this will cause the code to hang.
if 0.0 in data_points["eps"]:
    raise ValueError("Cannot test at eps=0.0.")

#%%
''' 3. PROCESS INPUTS 

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
    experimental["gs_spin_string"] = inputs["gs_spin"]

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
inputs["current_orbitals"] = fn.write_orbitals(inputs["fermi_level"]//2, inputs["num_orbs"], inputs["par"])
# inputs["current_orbitals"] = "-15 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38" #fn.write_orbitals(28, inputs["num_orbs"], inputs["par"]) #!! hard coded version of the above


#%%   
''' 4. WRITE GAMPN.DAT FILES 

- Loops through the data points to be tested:
    - Creates a tag referencing the nucleus, and the values of eps, gamma, and e2plus.
    - Uses string formatting with the inputs dictionary to write the .DAT file.
    - Creates/overwrites a .DAT file in the Inputs directory folder, named using the file tag. 

 '''

data_points["file_tags"] = []

for i in range(len(data_points["eps"])):
    
    inputs["current_eps"] = data_points["eps"][i]
    inputs["current_gamma"] = data_points["gamma_degrees"][i]
    
    if not(data_points["e2plus"][i]):
        # if the value of e2plus has been input as 0, calculate dynamically
        data_points["e2plus"][i] = fn.est_e2plus(data_points["eps"][i], inputs["A"])
    
    inputs["current_e2plus"] = data_points["e2plus"][i]
    
    _file_tag = "e%.3f_g%.1f_p%.3f_%s" % (inputs["current_eps"],  inputs["current_gamma"],
                                        inputs["current_e2plus"], inputs["nucleus"])
    
    inputs["current_f002"] = "f002_"+_file_tag+".dat"
    inputs["current_f016"] = "f016_"+_file_tag+".dat"
    inputs["current_f017"] = "f017_"+_file_tag+".dat"
    
    if data_points["gamma_degrees"][i] > 30:
        inputs["current_eps"] = data_points["eps"][i] * -1
        inputs["current_gamma"] = 60 - data_points["gamma_degrees"][i]
    
    fn.write_file("../"+nucleus+"/Inputs/GAM_"+_file_tag+".DAT", st.get_template("gampn") % inputs)
 
    data_points["file_tags"].append(_file_tag)

print("Number of data points = " + str(len(data_points["file_tags"])))
print("\nDeformation range being tested: \n\teps = [%.3f, %.3f], \n\tgamma = [%.1f, %.1f]."
      % (data_points["eps"][0], data_points["eps"][-1], data_points["gamma_degrees"][0], data_points["gamma_degrees"][-1]))



#%%
''' 5. WRITE AND RUN BASH SCRIPT TO EXECUTE GAMPN 

- Calculate batch settings:
    - How to divide up the data points into batches.
    - The maximum allowed time for the batch to run before assuming that it is hanging.
- Configure a script writer.
- Run the batches. The .OUT files are generated in the Outputs directory folder.

'''

# calculate batch settings

batch_settings = {}

batch_settings["num_batches"] = 8 # = number of cores for maximum efficiency with large data sets
batch_settings["num_per_batch"] = math.ceil(len(data_points["file_tags"]) / batch_settings["num_batches"])

# if the data set is small then use fewer cores for a minimum batch size of 20 to make the overhead worthwhile.
if batch_settings["num_per_batch"] < 20:
    batch_settings["num_per_batch"] = 20
    batch_settings["num_batches"] = math.ceil(len(data_points["eps"])/batch_settings["num_per_batch"])  

# each file takes ~ 0.1 seconds to run;
# allow double time plus an overhead of 2 s
batch_settings["allowed_time"] = 0.2*batch_settings["num_per_batch"]+50         


# configure a batch script writer for gampn, and run the batches
run_program = fn.configure_script_writer(data_points["file_tags"], nucleus, batch_settings["num_batches"], 
                                         batch_settings["num_per_batch"], batch_settings["allowed_time"], verbose)


print("Starting gampn...")
_sub_timer.start()

run_program("gampn")

_sub_timer.stop()
print("\n***** Finished running gampn in time = %.2f seconds. *****" % _sub_timer.get_lapsed_time())



#%%
''' 5. READ GAMPN.OUT FILE 

- For each data point (i.e. each GAMPN.OUT file):
    - Read the value of EFAC (the conversion factor from hw to eV).
    - Determine the line number of the fermi level orbital in the GAMPN.OUT file.
    - Read that line from GAMPN.OUT and get the energy, parity, and level index.
    - Generate the orbitals input for asyrmo (e.g. orbitals_string = "+11 19 20 21 22 23 24 25 26 27 28 29")
    
'''

# set up arrays to store data 
output_data = {}
output_data["fermi_parities"] = []
output_data["fermi_energies_hw"] = []
output_data["fermi_energies_mev"] = []
output_data["fermi_indices"] = []
data_points["asyrmo_orbitals"] = []


for i in data_points["file_tags"] :

    _lines = fn.read_file("../"+nucleus+"/Outputs/GAM_"+i+".OUT")
    
    # get the EFAC value (conversion factor from hw to MeV)
    inputs["efac"] = fn.get_efac(_lines)
    
    # locate the fermi level and extract data
    _fermi_level_line = fn.get_sp_level(_lines, inputs["fermi_level"], '0')
    
    _f_parity, _f_energy_hw, _f_index = fn.get_info(_fermi_level_line)
    
    output_data["fermi_parities"].append(_f_parity)                              
    output_data["fermi_energies_hw"].append(_f_energy_hw)   
    output_data["fermi_energies_mev"].append(_f_energy_hw * inputs["efac"])
    output_data["fermi_indices"].append(_f_index)
    
    # generate the orbitals input for asyrmo
    # version 5: dynamically finds the orbitals nearest to the fermi level in energy:
    data_points["asyrmo_orbitals"].append(fn.find_orbitals(_f_index, inputs["nu"], 
                                   inputs["par"], _f_energy_hw, _f_parity, _lines))
    



# Re-run gampn with the new set of orbitals, so that the strong-coupling basis can be maximised. 
# No need to re-read outputs, as the values read manually (energy levels etc) will be unaffected,
# only the matrix elements (read automatically by asyrmo) will be affected.
    
for i in range(len(data_points["file_tags"])):
    
    if data_points["gamma_degrees"][i] > 30:
        inputs["current_eps"] = data_points["eps"][i] * -1
        inputs["current_gamma"] = 60 - data_points["gamma_degrees"][i]
    else:
        inputs["current_eps"] = data_points["eps"][i]
        inputs["current_gamma"] = data_points["gamma_degrees"][i]
        
    inputs["current_e2plus"] = data_points["e2plus"][i]
    inputs["current_orbitals"] = data_points["asyrmo_orbitals"][i]
    
    inputs["current_f002"] = "f002_"+data_points["file_tags"][i]+".dat"
    inputs["current_f016"] = "f016_"+data_points["file_tags"][i]+".dat"
    inputs["current_f017"] = "f017_"+data_points["file_tags"][i]+".dat"
   
    
    fn.write_file("../"+nucleus+"/Inputs/GAM_"+data_points["file_tags"][i]+".DAT", st.get_template("gampn") % inputs)
 

_sub_timer.start()
run_program("gampn")
_sub_timer.stop()
print("***** Finished running gampn (again) in time = %.2f seconds. *****" % _sub_timer.get_lapsed_time())


_output_data = output_data # save a copy of the original before it's overwritten




#%%
''' 6. WRITE AND RUN ASYRMO 

- Use the existing list of file tags to write a .DAT file for each data point.
- Use the existing script writer to write and run asyrmo; dividing up the batches (as for gampn). 

'''

for i in range(len(data_points["file_tags"])):
        
    inputs["current_e2plus"] = data_points["e2plus"][i]
    inputs["current_orbitals"] = data_points["asyrmo_orbitals"][i]
    inputs["current_f016"] = "f016_"+data_points["file_tags"][i]+".dat"
    inputs["current_f017"] = "f017_"+data_points["file_tags"][i]+".dat"
    inputs["current_f018"] = "f018_"+data_points["file_tags"][i]+".dat"
    
    
    fn.write_file("../"+nucleus+"/Inputs/ASY_"+data_points["file_tags"][i]+".DAT", st.get_template("asyrmo") % inputs)

_sub_timer.start()
run_program("asyrmo")
_sub_timer.stop()

print("***** Finished running asyrmo in time = %.2f seconds. *****" % _sub_timer.get_lapsed_time())

#%%

''' 7. READ ASYRMO 

- Read the output files and check for the "NO DECOUPLING PARAMETERS CALCULATED" error.
- Record the value of the DELTA parameter.
- Check for the "SORRY I FOUND NO SOLUTIONS" error, and exclude those files from future analysis.

'''

output_data["delta"] = []

for i in data_points["file_tags"] :

    _lines = fn.read_file("../"+nucleus+"/Outputs/ASY_"+i+".OUT")
    
    if not "PARTICLE-ROTOR  MODEL" in _lines[0]: # then something has gone wrong
        raise RuntimeError("File " + i + " raised error in ASYRMO output: \n" + _lines[0] )
        
    output_data["delta"].append(fn.get_delta(_lines)) # this also checks for the "SORRY I FOUND NO SOLUTIONS" error.
        

#%%
''' 8. WRITE AND RUN PROBAMO 

- Use the existing list of file tags to write a .DAT file for each data point.
- Use the existing script writer to write and run probamo; dividing up the batches as for gampn. 

'''

for i in range(len(data_points["file_tags"])):

    inputs["current_e2plus"] = data_points["e2plus"][i]
    inputs["current_f017"] = "f017_"+data_points["file_tags"][i]+".dat"
    inputs["current_f018"] = "f018_"+data_points["file_tags"][i]+".dat"

    fn.write_file("../"+nucleus+"/Inputs/PROB_"+data_points["file_tags"][i]+".DAT", st.get_template("probamo") % inputs)

_sub_timer.start()
run_program("probamo")
_sub_timer.stop()

print("***** Finished running probamo in time = %.2f seconds. *****\n" % _sub_timer.get_lapsed_time())


#%%

''' 9. READ PROBAMO.OUT FILE

For each file:
    
- Read each line:
    - If it is a static transition, read the spin, energy, and magnetic moment.
    - Else ignore this line and move to the next.
    
- Sort the file data into categories:
    - Group lines by spin (e.g. spin 1/2 energies, spin 1/2 magnetic dipole moments, etc).
    - Additionally (separately) group lines by association with experimental data,
      (e.g. ground state energies, first excited state energies, etc).
      (this assumes that each input experimental excited state is the yrast state of that spin).
    - Fill missing gaps with NaN values, such that the same set of properties has 
      been recorded for every data point 
      
      (and if e.g. one data point found three spin 1/2 states, then all data points 
      should have a list of three spin 1/2 states, even if some of them are NaN)
    
    - Restructure the data set and save separately (now each property is recorded
      as a list of values for all data points, rather than each data point having
      a list of properties associated with it). 
      

'''

data_points["property_data"] = []

for i in range(len(data_points["file_tags"])):
    
    _file_data = {}
    
    _lines = fn.read_file("../"+nucleus+"/Outputs/Prob_"+data_points["file_tags"][i]+".OUT")

    for _line in _lines:
        
        _line_data = fn.read_data(_line)  # get the spin, energy, and magnetic moment from this line if it is a static moment, else return False
        
        if not(_line_data):
            continue # to next line in file
        
        # sort line_data into file_data according to its spin
        _file_data = fn.sort_by_spin(_line_data, _file_data)
        
        # additionally save data that corresponds to the experimental ground state and input excited states
        _file_data = fn.sort_by_expectation(_line_data, _file_data, inputs)
    
    _file_data = fn.missing_data(_file_data, inputs)
    data_points["property_data"].append(_file_data)

output_data = _output_data | fn.restructure_data(data_points["property_data"], inputs["ispin"], verbose)

# get energy gap between 9/2 and 13/2
output_data["gap_9_13"] = fn.find_gaps(output_data["spin_9/2_energies"], 3, output_data["spin_13/2_energies"], 1, 20) #!!!


# ensure all data sets have the same size and shape
# and mask ill-defined data poitns with reference to the DELTA data (any NaN values are ill-defined).

_mask = np.array([0 if np.isnan(x) else 1 for x in output_data["delta"]])

for i in output_data:
    if isinstance(output_data[i][0], list):
        output_data[i] = fn.fill_gaps(output_data[i])
        _list_mask = np.transpose(np.tile(_mask, (np.size(output_data[i][0]),1)))
        
        output_data[i] = np.where(_list_mask == 0, np.NaN, output_data[i])
    
    else: 
        output_data[i] = np.where(_mask == 0, np.NaN, output_data[i])

_output_data_dict = output_data # save a copy of the original before it's overwritten, so that the code can be run cell-by-cell without errors.


#%%
''' 10. PREPARE TO PLOT GRAPHS 

- Record each data set in an instance of class PropertyData
- Calculate graph plotting attributes and store within the class

- Raise a ValueError if the property isn't recognised 
  (i.e. if more data sets are read in the future, they cannot be plotted without
   first hard-coding the calculation of things like axis labels, contour levels,
   colour bar ticks, etc).

'''

# convert output_data from a dictionary of lists to a dictionary of PropertyData objects 
output_data = {}
for i in _output_data_dict:
    
    #print(i)
    
    output_data[i] = st.PropertyData(_output_data_dict[i], i)
    
    # calculate contour levels, colour bar ticks and labels, 
    # and assign experimental values and error tolerance if available.
    
    if output_data[i].prop == "energies" and not(output_data[i].sort=="Fermi"): 
        
        output_data[i].contour_levels = 10
        output_data[i].cbar_tick_labels = 0
        
        if output_data[i].sort == "gap":
            
            output_data[i].experimental_data, output_data[i].error_tolerance = gr.try_experimental(experimental, i, experimental["gap_en_tol"])
        
        elif output_data[i].sort == "Excited State ":
            
            output_data[i].experimental_data, output_data[i].error_tolerance = gr.try_experimental(experimental, "x" + output_data[i].num + "_energy", experimental["abs_en_tol"])
        
        elif output_data[i].sort == "Spin ":
            
            output_data[i].experimental_data, output_data[i].error_tolerance = gr.try_experimental(experimental, "jp_"+ output_data[i].num + inputs["par"], experimental["abs_en_tol"])
            
        else: raise ValueError("property not recognised: " + i)
        
    elif output_data[i].prop == "mag_moments":
        
        output_data[i].contour_levels = 6
        output_data[i].cbar_tick_labels = 0
            
        if output_data[i].sort == "Excited State ":
            
            output_data[i].experimental_data, output_data[i].error_tolerance = gr.try_experimental(experimental, "x" + output_data[i].num + "_mu", experimental["mu_tol"])
    
        elif output_data[i].sort == "Spin ":
            
            output_data[i].experimental_data = np.NaN
            output_data[i].error_tolerance = np.NaN
            
        elif output_data[i].sort == "Ground":
            
            output_data[i].experimental_data, output_data[i].error_tolerance = gr.try_experimental(experimental, "gs_mu", experimental["mu_tol"])
            
        else: raise ValueError("property not recognised: " + i)
        
    
    elif output_data[i].sort == "Ground":
        
        if output_data[i].prop == "spin_floats": 
        
            output_data[i].contour_levels = gr.calc_contour_levels(output_data[i].data)
            output_data[i].experimental_data, output_data[i].error_tolerance = gr.try_experimental(experimental, "gs_spin_float", experimental["mu_tol"])
            output_data[i].cbar_tick_labels = gr.calc_cbar_tick_labels(output_data[i].data, "half")
            output_data[i].cbar_ticks = gr.calc_cbar_ticks(output_data[i].contour_levels)
            
        elif output_data[i].prop == "spin_strings": continue
        else: raise ValueError("property not recognised: " + i)
    
    elif output_data[i].sort == "Fermi":
    
        if output_data[i].prop == "indices": 
            
            output_data[i].contour_levels = gr.calc_contour_levels(output_data[i].data)
            output_data[i].cbar_tick_labels = gr.calc_cbar_tick_labels(output_data[i].data, "int")
            output_data[i].cbar_ticks = gr.calc_cbar_ticks(output_data[i].contour_levels)
        
        elif output_data[i].prop == "energies":
        
            output_data[i].contour_levels = 10
            output_data[i].cbar_tick_labels = 0
        
        elif output_data[i].prop == "parities":
        
            output_data[i].contour_levels = 2
            output_data[i].cbar_tick_labels = 0
        
        else: raise ValueError("property not recognised: " + i)
            
        output_data[i].experimental_data = np.NaN
        output_data[i].error_tolerance = np.NaN
    
    elif output_data[i].prop == "delta":
        output_data[i].contour_levels = 10
        output_data[i].cbar_ticks = 0
        output_data[i].cbar_tick_labels = 0
        output_data[i].experimental_data = np.NaN
        output_data[i].error_tolerance = np.NaN
        
    else: raise ValueError("property not recognised: " + i)
    


output_data["all_energies"] = fn.collate_energy_data(output_data, len(data_points["file_tags"]), experimental["gs_spin_string"], experimental)

output_data["shifted_energies"] = fn.shift_energy_levels(output_data["all_energies"]) # recalculate all energies relative to the spin entered into fn.collate_energy_data() above

output_data["rms"] = fn.calc_rms_err(10, output_data["spin_1/2_energies"],output_data["spin_3/2_energies"], output_data["spin_5/2_energies"], output_data["spin_7/2_energies"], output_data["spin_9/2_energies"], output_data["spin_11/2_energies"], output_data["spin_13/2_energies"])




#%%


''' 11. PLOT GRAPHS

- Set which graphs should be plotted. Any not listed are False by default.

- If deformation was input as a mesh:
    - Plot filled contours in polar coordinates.
    - Draw a contour line to indicate the perimeter of the region where 
      the ground state spin was correctly reproduced, if requested.
    - Plot data point markers.
        - If experimental data is available, points that agree with experiment 
          (within tolerance) are marked in red.
        - Non-matching points are not plotted (unless there are fewer than 100 data points.)
    
- If only one of eps/gamma/e2plus is varied:
    - Plot line graphs.
    - Draw a green box around regions that have the correct ground state spin, if requested.
    - Plot a red line to indicate the experimental value, if available.
    
'''

# override settings
# inputs["mark_exp"] = 1
# inputs["mark_exp_tol"] = 0
# inputs["mark_points"] = 1
# inputs["mark_spin"] = 0

#!!! set graph subtitle:
subtitle = r'$E(2^+)$ = ' + str(inputs["current_e2plus"])

#!!! set which graphs to plot:

output_data["fermi_indices"].plot = 0
output_data["delta"].plot = 0
output_data["fermi_energies_mev"].plot = 0
output_data["fermi_energies_hw"].plot = 0

output_data["gs_mag_moments"].plot = 1
output_data["gs_spin_floats"].plot = 1

output_data["spin_1/2_energies"].plot = 0
output_data["spin_3/2_energies"].plot = 0
output_data["spin_5/2_energies"].plot = 0
output_data["spin_7/2_energies"].plot = 0
output_data["spin_9/2_energies"].plot = 0
output_data["spin_11/2_energies"].plot = 0
output_data["spin_13/2_energies"].plot = 0

output_data["spin_1/2_mag_moments"].plot = 0
output_data["spin_3/2_mag_moments"].plot = 0

output_data["x1_energies"].plot = 0
output_data["x2_energies"].plot = 0
output_data["x3_energies"].plot = 0

output_data["x1_mag_moments"].plot = 0

output_data["rms"].plot = 1

output_data["all_energies"].plot = 0
output_data["shifted_energies"].plot = 0

output_data["gap_9_13"].plot = 0


data_points["agreed"] = [0]*len(data_points["eps"])
num_comparisons = 0 

_sub_timer.start()



# start plotting graphs:
for i in output_data:
    
    prop = output_data[i]
    
    if not(prop.plot):
        continue
    inputs["current_graph"] = prop.title # makes several later inputs more efficient
    print("plotting graph: %(current_graph)s" % inputs) 
    
    if inputs["deformation_input"] == "mesh":  
        
        _fig, _ax = plt.subplots(subplot_kw=dict(projection='polar'))
        cax, cbar = gr.draw_contour_plot(_ax, prop, data_points)
        
        
        legend_handles = []
        
        if inputs["mark_spin"]:
            
            legend_handles = gr.mark_spin(inputs, data_points, output_data["gs_spin_floats"].data, legend_handles, _ax)
            
        # plot the data point markers, with comparison to experiment if possible
             
        legend_handles = gr.plot_points(data_points, prop, legend_handles, cbar, inputs)
        if np.isfinite(prop.experimental_data).all() and inputs["mark_exp"]: 
            num_comparisons += 1
        
        
        gr.format_fig('polar', _ax, legend_handles, '%(current_graph)s of %(nucleus)s' % inputs, subtitle)
        
        plt.show()
        
        
       
    elif (inputs["deformation_input"] ==  "gamma" 
          or inputs["deformation_input"] == "eps"
          or len(data_points["e2plus"]) > 1):
        
        # set which paramters are varied and which are constant
        var_sym, var, fix_sym, fix = gr.assign_parameters(inputs, data_points)
        
        _fig, _ax = plt.subplots() 
        
        legend_handles = []
        legend_handles, legend_title = gr.plot_line_data(data_points, prop, var, fix_sym, fix, legend_handles)

          
        # if experimental data is available, plot it in red for easy comparison
        if np.isfinite(prop.experimental_data).all() and not prop.num == "all": 
            legend_handles = gr.plot_exp_line(prop, inputs, var, legend_handles)

            
        # mark the range in which the correct ground state spin was calculated
        if inputs["mark_spin"]==1:
            
            correct_spin_range = gr.find_correct_spin(output_data["gs_spin_floats"].data, experimental["gs_spin_float"])
            if len(correct_spin_range) > 0:
                spin = gr.plot_correct_spin(correct_spin_range, var, inputs["step"], prop)
                legend_handles.append(spin)
                    
        gr.format_fig('linear', _ax, list(reversed(legend_handles)), 
                       '%(current_graph)s in %(nucleus)s' % inputs, subtitle, 
                       varied=var, x_label=var_sym, y_label=prop.axis_label, 
                       legend_title=legend_title)
        
        if prop.prop == "delta":
            _ax.set_ylim([0.2,1]) 
            
        if prop.cbar_tick_labels:        # then format for discrete values
            _ax.set_yticks(prop.cbar_ticks)
            _ax.set_yticklabels(prop.cbar_tick_labels)
        
        plt.show()
        
        
       

#%%
''' 12. ASSESS AGREEMENT OF CALCULATIONS WITH EXPERIMENT

- Print information to the console about the best agreement and its location.
- Plot a graph to show data point agreement across all data points.
- Print the mean energies of each level and the mean gs mag moment


'''
    
gr.check_agreement(verbose, data_points, num_comparisons)

agreement = st.PropertyData(data_points["agreed"], "Agreement of Data Points With Experimental Data")
agreement.contour_levels = np.arange(0, num_comparisons+2, dtype=int) #fn.calc_contour_levels(agreement.data)
agreement.cbar_ticks = gr.calc_cbar_ticks(agreement.contour_levels)
agreement.cbar_tick_labels = list(np.arange(0, num_comparisons+1, dtype=int)) #fn.calc_cbar_tick_labels(agreement.data, "int")
agreement.experimental_data = np.NaN
agreement.error_tolerance = np.NaN

agreement.plot = 0


if inputs["deformation_input"] == "mesh" and agreement.plot:  
    
    gr.plot_agreement(inputs, agreement, data_points, output_data, subtitle)
    
print("\n******** mean and standard error in the mean ******")

fn.report_mean(output_data["spin_1/2_energies"], verbose)
fn.report_mean(output_data["spin_3/2_energies"], verbose)
fn.report_mean(output_data["spin_5/2_energies"], verbose)
fn.report_mean(output_data["spin_7/2_energies"], verbose)
fn.report_mean(output_data["spin_9/2_energies"], verbose)
fn.report_mean(output_data["spin_11/2_energies"], verbose)
fn.report_mean(output_data["spin_13/2_energies"], verbose)
fn.report_mean(output_data["gs_mag_moments"], verbose)




# note how long it took
_sub_timer.stop()
_timer.stop()
print("\n****************************************************************************************")
print("finished plotting graphs in time = %.2f seconds" % (_sub_timer.get_lapsed_time()))
print("total runtime = %.2f seconds" % (_timer.get_lapsed_time()))
print("****************************************************************************************\n")





