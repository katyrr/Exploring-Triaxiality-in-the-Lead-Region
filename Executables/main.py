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
        
        (download the Examples folder from Github to make sure you have the right file structure)
    
    - Before running, check that config.txt has been filled in with the desired input parameters.
    - Also set which graphs to plot in section 11 of this file
    
    - Run from terminal in working directory "/Code/Executables" using command: "python3 main.py"
    

"""

#%% 
""" 1. SET UP 

- State the nucleus being studied (for navigating file directories).
- Import modules.
- Create timers (one to time the whole program, and one to time small sections).

"""

folder = "Pt177" #  "Examples" # "Au179" #  "Pb207" #!!! this is the name of the directory folder that contains your config file


import numpy as np                                  
import math                                        
import matplotlib.pyplot as plt                      

import functions as fn                              
import structs as st
import graph_plotting as gr

_timer = st.Timer()
_sub_timer = st.Timer()


#%%
''' 2. READ CONFIG FILE 

- Ignores empty lines, and lines beginning with * (to mark a comment).
- Checks that the format of each line is correct (var_name value),
  else raises a ValueError.
  
- Saves deformation input (converting any range inputs into arrays). 
- Raises a RuntimeError if: 
    deformation is input more than once, 
    or a range of E2PLUS values are input with a non-constant deformation,
    or both eps and gamma are input as a linear range (rather than a mesh).
- Raises a ValueError if any eps = 0.0, because this is not allowed in asyrmo.

- Converts all other inputs to the relevant type.
- All inputs are saved as name-value pairs in a dictionary.

''' 

_timer.start()


print("\nFetching config from folder: " + folder)
_lines = fn.read_file("../"+ folder +"/config.txt")

inputs = {}                                                                     
data_points = {}
experimental = {}
plot_props = {}

for i in range(len(_lines)):                       
    
    # remove comments and blank lines, and check formatting  
    _line = _lines[i].strip()
                                
    if _line == "" : continue                  
    if _line[0] == "*" : continue              
    
    _split_string = _line.split(" ")  # split into name and value
    _split_string = fn.remove_inline_comments(_split_string, i)

    fn.check_line_format(_split_string, _line, i)         
    
    # check what kind of input is stored in this line, and save it accordinly
    if _split_string[0] in ["eps", "gamma", "single","mesh"]:
        inputs, data_points = fn.save_deformation_input(inputs, data_points, _split_string)

    elif _split_string[0]=="e2plus":
        inputs, data_points = fn.save_e2plus_input(inputs, data_points, _split_string)
               
    elif _split_string[0] == "gs_spin":
        experimental = fn.save_gs_spin_input(experimental, _split_string)

    elif _split_string[0][:3]=="jp_":
        experimental[_split_string[0]] = [float(n) for n in _split_string[1].split(',')]
        
    elif _split_string[0][:5]=="plot_":
        plot_props[_split_string[0][5:]] = bool(int(_split_string[1]))
    
    elif _split_string[0][:6]=="engap_":
        experimental[_split_string[0]] = float(_split_string[1])
          
    else:
        inputs, experimental = fn.validate_input(inputs, experimental, _split_string)

    # sometimes helpful when debugging (change False to True, to check that config file is reading as expected) 
    if False: print(_split_string[0] + ": \t" + _split_string[1])


# check that we have all the inputs we need
for i in st.get_required_inputs():
    
    if not i in inputs and not i in experimental and not i in data_points:
       raise RuntimeError("missing input: " + i)

# check that there are no eps=0 values to test - that would cause the code to hang.
if 0.0 in data_points["eps"]:
    raise ValueError("Cannot test at eps=0.0.")
    

plt.rcParams['figure.dpi'] = inputs["figure_res"]  # set figure resolution

fn.setup_directory(folder, inputs["num_cores"], inputs["OS"])

#%%
''' 3. PROCESS INPUTS 

- Convert gamma points to radians and save separately.
- Convert input fractional spins to floats and save separately.
- convert the nantj, noutj, ipout inputs to the correct format

- Using input A and Z, work out which particle is odd, to determine the nneupr input.
- Halve and ceiling for the overall index of the fermi level orbital.
- Generate the orbitals input for gampn (e.g. "+4 19 20 21 22").

- Raises ValueError if an even-mass nucleus is input.

'''

# additionally save gamma in radians
data_points["gamma_radians"] = [n*np.pi/180 for n in data_points["gamma_degrees"]]


#!!! obsolete:
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
    raise RuntimeError("Check inputs of A and Z.")


# Generate a string containing the parity to calculate, the number of orbitals, 
#   and their indices, in the correct format for input to gampn.
inputs["current_orbitals"] = fn.write_orbitals(inputs["fermi_level"]//2, inputs["num_orbs"], inputs["par"])

# (useful for debugging) hard coded versions of the above:
#inputs["current_orbitals"] = "-15 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38" 
# inputs["current_orbitals"] = fn.write_orbitals(28, inputs["num_orbs"], inputs["par"])


#%%   
''' 4. WRITE GAMPN.DAT FILES 

- Loops through the data points to be tested:
    - Creates a tag referencing the nucleus, and the values of eps, gamma, and e2plus.
    - Uses string formatting with the inputs dictionary to write the .DAT file.
    - Creates/overwrites a .DAT file in the Inputs directory folder, named using the file tag. 

 '''

inputs["num"] = len(data_points["eps"])
print("Number of data points = ", inputs["num"])

data_points["file_tags"] = []

for i in range(inputs["num"]):
    
    if data_points["e2plus"][i] == 0:
        # if the value of e2plus has been input as 0, calculate dynamically
        data_points["e2plus"][i] = fn.est_e2plus(data_points["eps"][i], inputs["A"])
    
    inputs, data_points = fn.set_current(inputs, data_points, i)

    _file_tag = "e%.3f_g%.1f_p%.3f_%s" % (inputs["current_eps"],  inputs["current_gamma"],
                                        inputs["current_e2plus"], inputs["nucleus"])
    data_points["file_tags"].append(_file_tag)
    
    inputs["current_f002"] = "f002_"+_file_tag+".dat"
    inputs["current_f016"] = "f016_"+_file_tag+".dat"
    inputs["current_f017"] = "f017_"+_file_tag+".dat"
    

    fn.write_file("../"+folder+"/Inputs/GAM_"+_file_tag+".DAT", st.get_template("gampn") % inputs)
 
    

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

batch_settings["num_batches"] = inputs["num_cores"] # = number of cores for maximum efficiency with large data sets
batch_settings["num_per_batch"] = math.ceil(len(data_points["file_tags"]) / batch_settings["num_batches"])

# if the data set is small then use fewer cores for a minimum batch size of 20 to make the overhead worthwhile.
if batch_settings["num_per_batch"] < 20:
    batch_settings["num_per_batch"] = 20
    batch_settings["num_batches"] = math.ceil(len(data_points["eps"])/batch_settings["num_per_batch"])  

# Each file takes ~ 0.1 seconds to run;
#   Allow double time plus an overhead/extra of 10 seconds to ensure that the batch 
#   will finish even if the computer is running a bit slow today! If it takes longer 
#   than this, it is probably hanging, but you could increase 10 -> 60 seconds just to be sure.
batch_settings["allowed_time"] = 0.2*batch_settings["num_per_batch"]+10        


# configure a batch script writer, and run the batches.
run_program = fn.configure_script_writer(data_points["file_tags"], folder, batch_settings["num_batches"], 
                                         batch_settings["num_per_batch"], batch_settings["allowed_time"], 
                                         inputs["detailed_print"], inputs["OS"])

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
    - Dynamically locate the orbitals nearest to the fermi level in energy, for 
      future input.
      
- Re-run gampn for all data points, using the orbitals located in the step above. 
    
'''

# set up arrays to store data 
output_data = {"fermi_parities": [0]*inputs["num"], "fermi_energies_hw": [0]*inputs["num"], 
               "fermi_energies_mev": [0]*inputs["num"], "fermi_indices": [0]*inputs["num"]}
data_points["asyrmo_orbitals"] = []

for i in range(inputs["num"]):

    _lines = fn.read_file("../"+folder+"/Outputs/GAM_"+data_points["file_tags"][i]+".OUT")
    
    inputs["efac"] = fn.get_efac(_lines)
    _fermi_level_line = fn.get_sp_level(_lines, inputs["fermi_level"], '0')
    _f_parity, _f_energy_hw, _f_index = fn.get_info(_fermi_level_line)
    
    output_data["fermi_parities"][i] = _f_parity                             
    output_data["fermi_energies_hw"][i] = _f_energy_hw
    output_data["fermi_energies_mev"][i] = _f_energy_hw * inputs["efac"]
    output_data["fermi_indices"][i] = _f_index
    
    # dynamically finds the orbitals nearest to the fermi level in energy:
    data_points["asyrmo_orbitals"].append(fn.find_orbitals(_f_index, inputs["nu"], 
                                   inputs["par"], _f_energy_hw, _f_parity, _lines))


#%% 
'''6. RE-RUN GAMPN

- Re-run gampn with the new set of orbitals, so that the strong-coupling basis 
  can be maximised (to 15 orbitals) when calculating matrix elements.
- No need to re-read the outputs, because the properties we read earlier are 
  not affected, and the recalculated matrix elements will be passed to the next
  program automatically.


'''


for i in range(len(data_points["file_tags"])):
    
    inputs, data_points = fn.set_current(inputs, data_points, i, file_tag=data_points["file_tags"][i])
    inputs["current_orbitals"] = data_points["asyrmo_orbitals"][i]
   
    fn.write_file("../"+folder+"/Inputs/GAM_"+data_points["file_tags"][i]+".DAT", st.get_template("gampn") % inputs)
 
_sub_timer.start()
run_program("gampn")
_sub_timer.stop()
print("***** Finished running gampn (again) in time = %.2f seconds. *****" % _sub_timer.get_lapsed_time())

_output_data = output_data # save a copy of the original before it's overwritten (useful when running cell by cell)



#%%
''' 7. WRITE AND RUN ASYRMO 

- Use the existing list of file tags to write a .DAT file for each data point.
- Use the existing script writer to write and run asyrmo; dividing up the batches (as for gampn). 

'''

for i in range(inputs["num"]):
        
    inputs["current_e2plus"] = data_points["e2plus"][i]
    inputs["current_orbitals"] = data_points["asyrmo_orbitals"][i]
    
    inputs["current_f016"] = "f016_"+data_points["file_tags"][i]+".dat"
    inputs["current_f017"] = "f017_"+data_points["file_tags"][i]+".dat"
    inputs["current_f018"] = "f018_"+data_points["file_tags"][i]+".dat"
    
    fn.write_file("../"+folder+"/Inputs/ASY_"+data_points["file_tags"][i]+".DAT", st.get_template("asyrmo") % inputs)

_sub_timer.start()
run_program("asyrmo")
_sub_timer.stop()

print("***** Finished running asyrmo in time = %.2f seconds. *****" % _sub_timer.get_lapsed_time())

#%%

''' 8. READ ASYRMO 

- Read the output files and check for the "NO DECOUPLING PARAMETERS CALCULATED" error.
- Record the value of the DELTA parameter.
- Check for the "SORRY I FOUND NO SOLUTIONS" error, and exclude those files from future analysis.

'''

output_data["delta"] = []

for i in data_points["file_tags"]:

    _lines = fn.read_file("../"+folder+"/Outputs/ASY_"+i+".OUT")
    
    if not "PARTICLE-ROTOR  MODEL" in _lines[0]: # then something has gone wrong
        raise RuntimeError("File " + i + " raised error in ASYRMO output: \n" + _lines[0] )
        
    output_data["delta"].append(fn.get_delta(_lines)) # this also checks for the "SORRY I FOUND NO SOLUTIONS" error.
        

#%%
''' 9. WRITE AND RUN PROBAMO 

- Use the existing list of file tags to write a .DAT file for each data point.
- Use the existing script writer to write and run probamo; dividing up the batches as for gampn. 

'''

for i in range(inputs["num"]):

    inputs["current_e2plus"] = data_points["e2plus"][i]
    inputs["current_f017"] = "f017_"+data_points["file_tags"][i]+".dat"
    inputs["current_f018"] = "f018_"+data_points["file_tags"][i]+".dat"

    fn.write_file("../"+folder+"/Inputs/PROB_"+data_points["file_tags"][i]+".DAT", st.get_template("probamo") % inputs)

_sub_timer.start()
run_program("probamo")
_sub_timer.stop()

print("***** Finished running probamo in time = %.2f seconds. *****\n" % _sub_timer.get_lapsed_time())


#%%

''' 10. READ PROBAMO.OUT FILE

For each file:
    
- Read each line:
    - If it is a static transition, read the spin, energy, and magnetic moment.
    - Else ignore this line and move to the next.
    
- Sort the file data into categories:
    - Group lines by spin (e.g. spin 1/2 energies, spin 1/2 magnetic dipole moments, etc).
    - Additionally (separately) record the expected ground state (the lowest state with the 
      same spin as the experimental gs) properties as a group.
    - Fill missing gaps with NaN values, such that the same set of properties has 
      been recorded for every data point (and if e.g. one data point found three 
      spin 1/2 states, then all data points should have a list of three spin 1/2 states, 
      even if some of them are NaN).
    - Restructure the data set and save separately (now each property is recorded
      as a list of values for all data points, rather than each data point having
      a list of properties associated with it). 

- Calculate energy gaps between levels specified in the experimental section of the config input.

- Ensure all data sets have the same size and shape.
- Mask bad data points with reference to the DELTA data 
  (any DETLA = NaN values are bad, caused by some kind of convergence issue with BCS pairing).

'''

data_points["property_data"] = []

for i in range(inputs["num"]):
    
    _lines = fn.read_file("../"+folder+"/Outputs/Prob_"+data_points["file_tags"][i]+".OUT")
    
    _file_data = {}
    for _line in _lines:
        
        _line_data = fn.read_data(_line)  # get the spin, energy, and magnetic moment from this line if it is a static moment, else return False
        if not(_line_data):
            continue # to next line in file
        
        # sort line_data into file_data according to its spin
        _file_data = fn.sort_by_spin(_line_data, _file_data)
        # additionally save data that corresponds to the expected experimental ground state
        _file_data = fn.sort_by_expectation(_line_data, _file_data, inputs)
    
    _file_data = fn.missing_data(_file_data, inputs)
    data_points["property_data"].append(_file_data)

output_data = _output_data | fn.restructure_data(data_points["property_data"], inputs["ispin"], inputs["detailed_print"])

# get energy gap between third 9/2 and first 13/2 states
for i in experimental:
    if not "engap_" in i:
        continue
    
    _spin1, _idx1, _spin2, _idx2 = fn.parse_engap_input(i)
    
    output_data[i] = fn.find_gaps(output_data["spin_"+_spin1+"/2_energies"], _idx1, output_data["spin_"+_spin2+"/2_energies"], _idx2, experimental[i])

# output_data["engap_9.3_13.1"] = fn.find_gaps(output_data["spin_9/2_energies"], 3, output_data["spin_13/2_energies"], 1, 20) #!!!


# ensure all data sets have the same size and shape, and mask bad points

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
''' 11. PREPARE TO PLOT GRAPHS 

- Record each data set in an instance of class PropertyData.
- Calculate graph plotting attributes and store within the class.

- Raise a ValueError if the property isn't recognised 
  (i.e. if more data sets are recorded in the future, they cannot be plotted without
   first hard-coding the calculation of things like axis labels, contour levels,
   colour bar ticks, etc).
  
- Create a new data set containing all energies (of all spins) to plot together.
- Create a new data set with all energies shifted to be relative to the expected 
  ground state (not necessarily the same as the calculated ground state at all points).
  This makes the output lines look smoother (no sharp bends when the ground state changes).
- Create a new data set containing root mean squared error (i.e. discrepancy) between 
  the calculated lowest energy states of each spin and the exeperimental values (where available).

'''

# convert output_data from a dictionary of lists to a dictionary of PropertyData objects 
output_data = {}
for i in _output_data_dict:
    # print(i)
    output_data[i] = st.PropertyData(_output_data_dict[i], i)
    
    # calculate contour levels, colour bar ticks and labels, 
    # and assign experimental values and error tolerance if available.
    
    output_data[i] = gr.calculate_format_data(output_data[i], i, experimental)
    

output_data["all_energies"] = fn.collate_energy_data(output_data, len(data_points["file_tags"]), 
                                                     experimental["gs_spin_string"], experimental)

# recalculate all energies relative to the spin entered into fn.collate_energy_data() above
output_data["shifted_energies"] = fn.shift_energy_levels(output_data["all_energies"]) 

output_data["rms"] = fn.calc_rms_err(10, output_data["spin_1/2_energies"],
                     output_data["spin_3/2_energies"], output_data["spin_5/2_energies"], 
                     output_data["spin_7/2_energies"], output_data["spin_9/2_energies"], 
                     output_data["spin_11/2_energies"], output_data["spin_13/2_energies"])




#%%


''' 12. PLOT GRAPHS

- Set a subtitle containing the values of E2PLUS and GSFAC input, if requested.
- Set which graphs should be plotted (from config, or overwritten below). 
  Any not listed are False by default.

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


#!!! set graph subtitle:
if inputs["include_subtitle"]:
    subtitle = r'$E(2^+)$ = ' + str(inputs["current_e2plus"]) + '; gsfac = ' + str(inputs["gsfac"])
else:
    subtitle = ''
    
    
# set which graphs to plot:
for i in plot_props:
    if i in output_data:
        output_data[i].plot = plot_props[i]
    else:
        print("property not recorded, check that it is included in experimental data inputs:\n\t", i)

# override graph plotting options (useful when running cell by cell): #!!!
    
# output_data["fermi_indices"].plot = 0
# output_data["delta"].plot = 0
# output_data["fermi_energies_mev"].plot = 0
# output_data["fermi_energies_hw"].plot = 0

# output_data["gs_mag_moments"].plot = 0
# output_data["gs_quad_moments"].plot = 0
# output_data["gs_spin_floats"].plot = 1

# output_data["spin_1/2_energies"].plot = 0
# output_data["spin_3/2_energies"].plot = 0
# output_data["spin_5/2_energies"].plot = 0
# output_data["spin_7/2_energies"].plot = 0
# output_data["spin_9/2_energies"].plot = 0
# output_data["spin_11/2_energies"].plot = 0
# output_data["spin_13/2_energies"].plot = 0

# output_data["spin_1/2_mag_moments"].plot = 0
# output_data["spin_3/2_mag_moments"].plot = 0

# output_data["rms"].plot = 0

# output_data["all_energies"].plot = 0
# output_data["shifted_energies"].plot = 0

# output_data["gap_9_13"].plot = 0

# override settings
# inputs["mark_exp"] = 1
# inputs["mark_exp_tol"] = 0
# inputs["mark_points"] = 1
# inputs["mark_spin"] = 0


_sub_timer.start()

data_points["agreed"] = [0]*len(data_points["eps"])
num_comparisons = 0 

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
''' 13. ASSESS AGREEMENT OF CALCULATIONS WITH EXPERIMENT

- Print information about the best agreement and its location.
- Plot a graph to show data point agreement across all data points.
- Print the mean energies of each level and the mean gs moments, with standard error.
- Print the total runtime.

'''
    
gr.check_agreement(inputs["detailed_print"], data_points, num_comparisons)

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

fn.report_mean(output_data["spin_1/2_energies"], inputs["detailed_print"])
fn.report_mean(output_data["spin_3/2_energies"], inputs["detailed_print"])
fn.report_mean(output_data["spin_5/2_energies"], inputs["detailed_print"])
fn.report_mean(output_data["spin_7/2_energies"], inputs["detailed_print"])
fn.report_mean(output_data["spin_9/2_energies"], inputs["detailed_print"])
fn.report_mean(output_data["spin_11/2_energies"], inputs["detailed_print"])
fn.report_mean(output_data["spin_13/2_energies"], inputs["detailed_print"])
fn.report_mean(output_data["gs_mag_moments"], inputs["detailed_print"])
fn.report_mean(output_data["gs_quad_moments"], inputs["detailed_print"])

# note how long it took
_sub_timer.stop()
_timer.stop()
print("\n****************************************************************************************")
print("finished plotting graphs in time = %.2f seconds" % (_sub_timer.get_lapsed_time()))
print("total runtime = %.2f seconds" % (_timer.get_lapsed_time()))
print("****************************************************************************************\n")





