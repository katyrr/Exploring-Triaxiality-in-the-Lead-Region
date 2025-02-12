#!/usr/bin/env python3
# -*- coding: utf-8 -*-
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

plt.rcParams['figure.dpi'] = 150
timer = st.Timer()
sub_timer = st.Timer()

nucleus = "Pb207"
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

timer.start()

lines = fn.read_file("../"+ nucleus +"/config.txt")

inputs = {}                                                                     
data_points = {}

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
            data_points["eps"], data_points["gamma_degrees"], inputs["step"] = fn.arrange_data_line(split_string)
                   
       
    # save E2PLUS input #!!! implement grodzins option?
    elif split_string[0]=="e2plus":
        
        e2plus_to_test, _ = fn.range_to_list(split_string[1])
        
        if (len(e2plus_to_test)>1 and not inputs["deformation_input"] == "single"):
            raise RuntimeError('''Testing a range of e2plus is only supported 
                             for a single deformation input.''')
    
    # save other inputs
    elif split_string[0][:3]=="jp_":
        inputs[split_string[0]] = [float(n) for n in split_string[1].split(',')]
          
    elif split_string[0] in st.get_variable_list("int"):
        inputs[split_string[0]] = int(split_string[1])
        
    elif split_string[0] in st.get_variable_list("float"):
        inputs[split_string[0]] = float(split_string[1])
                                      
    elif split_string[0] in st.get_variable_list("bool"):                                    
        inputs[split_string[0]] = bool(split_string[1])
        
    else: # string
        inputs[split_string[0]] = split_string[1]


    if verbose: print(split_string[0] + ": \t" + split_string[1])
    


# check that there are no eps=0 values to test - this will cause the code to hang.
if 0.0 in data_points["eps"]:
    raise ValueError("Cannot test at eps=0.0.")


print("\n********** Finished reading lines: "+str(l)+ " **********\n")

# deallocate variables that are no longer needed
del [l, line, lines, split_string]




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
    inputs["gs_spin_float"] = fn.spin_string_to_float(inputs["gs_spin"])

if "x1_spin" in inputs:
    inputs["x1_spin_float"] = fn.spin_string_to_float(inputs["x1_spin"])

if "x2_spin" in inputs:
    inputs["x2_spin_float"] = fn.spin_string_to_float(inputs["x2_spin"])

if "x3_spin" in inputs:
    inputs["x3_spin_float"] = fn.spin_string_to_float(inputs["x3_spin"])


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
inputs["gampn_orbitals"] = fn.write_orbitals(inputs["fermi_level"]//2, inputs["num_orbs"], inputs["par"])


#%%   
''' 4. WRITE GAMPN.DAT FILES 

- Loops through the data points to be tested:
    - Creates a tag referencing the nucleus, and the values of eps, gamma, and e2plus.
    - Uses string formatting with the inputs dictionary to write the .DAT file.
    - Creates/overwrites a .DAT file in the Inputs directory folder, named using the file tag. 

 '''

data_points["file_tags"] = []

for p in range(len(data_points["eps"])):
    for etp in e2plus_to_test:
        
        inputs["current_e2plus"] = etp
        inputs["current_eps"] = data_points["eps"][p]
        inputs["current_gamma"] = data_points["gamma_degrees"][p]
        
        file_tag = "e%.3f_g%.1f_p%.3f_%s" % (inputs["current_eps"],  inputs["current_gamma"],
                                            inputs["current_e2plus"], inputs["nucleus"])
        
        inputs["current_f002"] = "f002_"+file_tag+".dat"
        inputs["current_f016"] = "f016_"+file_tag+".dat"
        inputs["current_f017"] = "f017_"+file_tag+".dat"
        
        gampn_dat_text = st.get_template("gampn") % inputs
               
        fn.write_file("../"+nucleus+"/Inputs/GAM_"+file_tag+".DAT", gampn_dat_text)
 
        data_points["file_tags"].append(file_tag)

print("Number of data points = " + str(len(data_points["file_tags"])))
print("\nDeformation range being tested: \n\teps = [%.3f, %.3f], \n\tgamma = [%.1f, %.1f]."
      % (data_points["eps"][0], data_points["eps"][-1], data_points["gamma_degrees"][0], data_points["gamma_degrees"][-1]))


del [etp, p, file_tag, gampn_dat_text]



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
    batch_settings["num_batches"] = math.ceil(len(e2plus_to_test)*len(data_points["eps"])/batch_settings["num_per_batch"])  

batch_settings["allowed_time"] = 0.5*batch_settings["num_per_batch"]+5         #!!! each file takes ~ 0.1 seconds to run, and allow an overhead of 1s


# configure a batch script writer for gampn, and run the batches
run_program = fn.configure_script_writer(data_points["file_tags"], nucleus, batch_settings["num_batches"], 
                                         batch_settings["num_per_batch"], batch_settings["allowed_time"], verbose)

sub_timer.start()

run_program("gampn")

sub_timer.stop()
print("\n***** Finished running gampn in time = %.2f seconds. *****" % sub_timer.get_lapsed_time())



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


for file in data_points["file_tags"] :

    lines = fn.read_file("../"+nucleus+"/Outputs/GAM_"+file+".OUT")
    
    # get the EFAC value (conversion factor from hw to MeV)
    inputs["efac"] = fn.get_efac(lines)
    
    # locate the fermi level
    fermi_level_line = fn.get_fermi_level(lines, inputs["fermi_level"])
    
    # get data about the fermi level using the '#' as a reference point 
    hash_index = fermi_level_line.index("#")
    output_data["fermi_parities"].append(fermi_level_line[hash_index-2])                              
    output_data["fermi_energies_hw"].append(float(fermi_level_line[hash_index-10 : hash_index-4]))   
    output_data["fermi_energies_mev"].append(output_data["fermi_energies_hw"][-1]*inputs["efac"])
    output_data["fermi_indices"].append(fn.get_fermi_index(fermi_level_line, hash_index))
    
    # generate the orbitals input for asyrmo
    data_points["asyrmo_orbitals"].append(fn.write_orbitals(output_data["fermi_indices"][-1], inputs["nu"], inputs["par"]))
    
del [fermi_level_line, file, hash_index, lines]

_output_data = output_data # save a copy of the original before it's overwritten


#%%
''' 6. WRITE AND RUN ASYRMO 

- Use the existing list of file tags to write a .DAT file for each data point.
- Use the existing script writer to write and run asyrmo; dividing up the batches as for gampn. 

'''

for t in range(len(data_points["file_tags"])):
    
    if len(e2plus_to_test)>1:
        inputs["current_e2plus"] = e2plus_to_test[t]
    else:
        inputs["current_e2plus"] = e2plus_to_test[0]

    inputs["current_orbitals"] = data_points["asyrmo_orbitals"][t]
    inputs["current_f016"] = "f016_"+data_points["file_tags"][t]+".dat"
    inputs["current_f017"] = "f017_"+data_points["file_tags"][t]+".dat"
    inputs["current_f018"] = "f018_"+data_points["file_tags"][t]+".dat"
    
    asyrmo_dat_text = st.get_template("asyrmo") % inputs
    
    fn.write_file("../"+nucleus+"/Inputs/ASY_"+data_points["file_tags"][t]+".DAT", asyrmo_dat_text)

del [t, asyrmo_dat_text]


sub_timer.start()

run_program("asyrmo")

sub_timer.stop()

print("***** Finished running asyrmo in time = %.2f seconds *****" % sub_timer.get_lapsed_time())



#%%
''' 7. WRITE AND RUN PROBAMO 

- Use the existing list of file tags to write a .DAT file for each data point.
- Use the existing script writer to write and run probamo; dividing up the batches as for gampn. 

'''

for t in range(len(data_points["file_tags"])):
    
    if len(e2plus_to_test)>1:
        inputs["current_e2plus"] = e2plus_to_test[t]
    else:
        inputs["current_e2plus"] = e2plus_to_test[0]
    
    inputs["current_f017"] = "f017_"+data_points["file_tags"][t]+".dat"
    inputs["current_f018"] = "f018_"+data_points["file_tags"][t]+".dat"

    probamo_dat_text = st.get_template("probamo") % inputs
    
    fn.write_file("../"+nucleus+"/Inputs/PROB_"+data_points["file_tags"][t]+".DAT", probamo_dat_text)


del [probamo_dat_text, t]

sub_timer.start()

run_program("probamo")

sub_timer.stop()

print("***** Finished running probamo in time = %.2f seconds *****\n" % sub_timer.get_lapsed_time())


#%%

''' 8. READ PROBAMO.OUT FILE

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


# start reading files 

data_points["property_data"] = []

for t in range(len(data_points["file_tags"])):
    
    file_data = {}
    
    lines = fn.read_file("../"+nucleus+"/Outputs/Prob_"+data_points["file_tags"][t]+".OUT")

    for line in lines:
        
        line_data = fn.read_data(line)  # get the spin, energy, and magnetic moment from this line if it is a static moment, else return False
        
        if not(line_data):
            continue # to next line in file
        
        # sort line_data into file_data according to its spin
        file_data = fn.sort_by_spin(line_data, file_data)
        
        # additionally save data that corresponds to the experimental ground state and input excited states
        file_data = fn.sort_by_expectation(line_data, file_data, inputs)
    
    file_data = fn.missing_data(file_data, inputs)
    data_points["property_data"].append(file_data)

    
output_data = _output_data | fn.restructure_data(data_points["property_data"], inputs["ispin"], verbose)

# ensure all data sets have the same size and shape
for d in output_data:
    if isinstance(output_data[d][0], list):
        output_data[d] = fn.fill_gaps(output_data[d])

    
del [file_data, line, line_data, lines, t]

_output_data_dict = output_data # save a copy of the original before it's overwritten



#%%
''' 9. PREPARE TO PLOT GRAPHS 

- Record each data set in an instance of class PropertyData
- Calculate graph plotting attributes and store within the class

- Raise a ValueError if the property isn't recognised 
  (i.e. if more data sets are read in the future, they cannot be plotted without
   first hard-coding the calculation of things like axis labels, contour levels,
   colour bar ticks, etc).

'''

# convert output_data from a dictionary of lists to a dictionary of PropertyData objects 
output_data = {}
for p in _output_data_dict:
    
    output_data[p] = st.PropertyData(_output_data_dict[p], p)
    
    # calculate contour levels, colour bar ticks and labels, 
    # and assign experimental values and error tolerance if available.
    
    if output_data[p].prop == "energies" and not(output_data[p].sort=="Fermi"): 
        
        output_data[p].contour_levels = 10
        output_data[p].cbar_tick_labels = 0
        
        if output_data[p].sort == "Excited State ":
            
            output_data[p].experimental_data, output_data[p].error_tolerance = fn.try_experimental(inputs, "x" + output_data[p].num + "_energy", 50)
        
        elif output_data[p].sort == "Spin ":
            
            output_data[p].experimental_data, output_data[p].error_tolerance = fn.try_experimental(inputs, "jp_"+ output_data[p].num + inputs["par"], 50)
            
        else: raise ValueError("property not recognised: " + p)
        
    elif output_data[p].prop == "mag_moments":
            
        if output_data[p].sort == "Excited State ":
            
            output_data[p].experimental_data, output_data[p].error_tolerance = fn.try_experimental(inputs, "x" + output_data[p].num + "_mu", 0.2)
    
        elif output_data[p].sort == "Spin ":
            
            output_data[p].experimental_data = np.NaN
            output_data[p].error_tolerance = np.NaN
            
        elif output_data[p].sort == "Ground":
            
            output_data[p].experimental_data, output_data[p].error_tolerance = fn.try_experimental(inputs, "gs_mu", 0.2)
            
        else: raise ValueError("property not recognised: " + p)
        
        output_data[p].contour_levels = 8
        output_data[p].cbar_tick_labels = 0
    
    elif output_data[p].sort == "Ground":
        
        if output_data[p].prop == "spin_floats": 
        
            output_data[p].contour_levels = fn.calc_contour_levels(output_data[p].data)
            output_data[p].experimental_data, output_data[p].error_tolerance = fn.try_experimental(inputs, "gs_spin_float", 0.2)
            output_data[p].cbar_tick_labels = fn.calc_cbar_tick_labels(output_data[p].data)
            output_data[p].cbar_ticks = fn.calc_cbar_ticks(output_data[p].contour_levels)
            
        elif output_data[p].prop == "spin_strings": continue
    
    elif output_data[p].sort == "Fermi":
    
        if output_data[p].prop == "indices": 
            
            output_data[p].contour_levels = fn.calc_contour_levels(output_data[p].data)
            output_data[p].cbar_tick_labels = fn.calc_cbar_tick_labels(output_data[p].data)
            output_data[p].cbar_ticks = fn.calc_cbar_ticks(output_data[p].contour_levels)
        
        elif output_data[p].prop == "energies":
        
            output_data[p].contour_levels = 10
            output_data[p].cbar_tick_labels = 0
        
        elif output_data[p].prop == "parities":
        
            output_data[p].contour_levels = 2
            output_data[p].cbar_tick_labels = 0
        
        else: raise ValueError("property not recognised: " + p)
            
        output_data[p].experimental_data = np.NaN
        output_data[p].error_tolerance = np.NaN
    
    else: raise ValueError("property not recognised: " + p)
    
    
del [p]


#%%


''' 10. PLOT GRAPHS

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

# set which graphs to plot:
output_data["gs_mag_moments"].plot = False
output_data["gs_spin_floats"].plot = False
output_data["spin_1/2_energies"].plot = False
output_data["spin_3/2_energies"].plot = True
output_data["spin_5/2_energies"].plot = True
output_data["x1_energies"].plot = True
output_data["x2_energies"].plot = True

data_points["agreed"] = [0]*len(data_points["eps"])
num_comparisons = 0 

sub_timer.start()


# start plotting graphs:
for g in output_data:
    
    prop = output_data[g]
    if not(prop.plot):
        continue
    inputs["current_graph"] = prop.title # makes several later inputs more efficient
    print("plotting graph: %(current_graph)s" % inputs) 
    
    if inputs["deformation_input"] == "mesh":  
        
        fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
        cax, cbar = fn.draw_contour_plot(ax, prop, data_points)
        
        
        legend_handles = []
        
        if inputs["mark_spin"]:
            
            legend_handles = fn.mark_spin(inputs, data_points, output_data["gs_spin_floats"].data, legend_handles, ax)
            
        # plot the data point markers, with comparison to experiment if possible 
        if np.isfinite(prop.experimental_data).all(): 
            num_comparisons += 1
            legend_handles = fn.plot_points_with_experiment(data_points, prop, legend_handles, cbar)
            
        else: 
            legend_handles =  fn.plot_points_without_experiment(data_points, legend_handles)
               
        
        fn.format_fig('polar', ax, legend_handles, '%(current_graph)s of %(nucleus)s' % inputs)
        
        plt.show()
        
        
       
    elif (inputs["deformation_input"] ==  "gamma" 
          or inputs["deformation_input"] == "eps"
          or len(e2plus_to_test) > 1):
        
        # set which paramters are varied and which are constant
        var_sym, var, fix_sym, fix = fn.assign_parameters(inputs, e2plus_to_test, data_points)
        
        fig, ax = plt.subplots() 
        
        legend_handles = []
        
        # now plot the actual data
        if len(data_points["file_tags"]) < 100: marker_size = 5
        else: marker_size = 1 # use smaller markers if the data set is large
        
        if prop.sort == "Spin ":
        
            legend_handles, legend_title = fn.plot_multi_lines(prop, var, legend_handles, marker_size, fix_sym, fix[0])
        
        else:
            data, = plt.plot(var, prop.data, 'k-x', markersize=marker_size, label="%s = %s" % (fix_sym, fix[0]))
            legend_handles.append(data)
            legend_title = ""
          
        # if experimental data is available, plot it in red for easy comparison
        if np.isfinite(prop.experimental_data).all():  
            
            if prop.sort == "Spin ":
                for e in range(len(prop.experimental_data)):
                    exp, = plt.plot(var, np.full(len(var), prop.experimental_data[e]), 
                       'r-', label="experimental value")
            else:
                exp, = plt.plot(var, np.full(len(var), prop.experimental_data), 
                       'r-', label="experimental value")
                
            
            legend_handles.append(exp)
            
    
        # mark the range in which the correct ground state spin was calculated
        if inputs["mark_spin"]==1:
            
            correct_spin_range = fn.find_correct_spin(output_data["gs_spin_floats"].data, inputs["gs_spin_float"])
            if len(correct_spin_range) > 0:
                spin = fn.plot_correct_spin(correct_spin_range, var, inputs["step"], prop)
                legend_handles.append(spin)
            
                    
        fn.format_fig('linear', ax, list(reversed(legend_handles)), 
                       '%(current_graph)s in %(nucleus)s' % inputs, 
                       varied=var, x_label=var_sym, y_label=prop.axis_label, 
                       legend_title=legend_title)
        
        plt.show()
        
        del [var, var_sym, fix, fix_sym, legend_title, legend_handles]
        
       
del [g]

#%%

if len(data_points["file_tags"]) > 1:
    
    fn.check_agreement(verbose, data_points, num_comparisons)

# note how long it took
sub_timer.stop()
timer.stop()
print("\n****************************************************************************************")
print("finished plotting graphs in time = %.2f seconds" % (sub_timer.get_lapsed_time()))
print("total runtime = %.2f seconds" % (timer.get_lapsed_time()))
print("****************************************************************************************\n")





