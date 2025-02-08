#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  7 10:53:52 2025

@author: katyrr


A separate module file to contain function definitions.

"""

import numpy as np                                                              # for np.arange(start, end, step)
import subprocess                                                               # for calling shell scripts to run 
import matplotlib.pyplot as plt                                                 # for plotting graphs
from matplotlib.ticker import FuncFormatter                                     # for formatting axis ticks
import matplotlib.tri as tri                                                    # for manual triangulation before drawing a contour plot


''' GENERAL FUNCTIONS '''

def spin_string_to_float(spin_string):
    if spin_string[1] == "/":
        spin_float = float(spin_string[0])/2
    elif spin_string[2] == "/":
        spin_float = float(spin_string[0:2])/2 
    else:
        raise ValueError("Cannot parse spin. \nCheck it has been input in " +
                         "the format 'n/2' where n is an int.")
    return spin_float

def spin_float_to_string(spin_float):
    numerator = float(spin_float)*2
    if str(numerator)[-1] != "0":
        raise ValueError("Cannot convert float to fraction. \nCheck it has " + 
                         "been input as a half-integer float.")
    else:
        return str(int(numerator))+"/2"
    
def read_file(path):
    with open(path, 'r') as f:
        lines = f.readlines()
    
    return lines

def write_file(path, text):
    
    with open(path, 'w') as f:
        f.write(text)






''' FUNCTIONS FOR READING CONFIG '''

def remove_inline_comments(split_string):
    for n in range(len(split_string)):                                          
        word = split_string[n]
        if word[0] == '*':
            return split_string[:n]
    return split_string

def check_line_format(split_string, line, l):
    
    if (split_string[0]=="eps" or split_string[0]=="gamma" or split_string[0]=="single"): expected_num_words = 3       
    else: expected_num_words = 2                                                # all other variables have only one value associated with the variable name
    
    if len(split_string)>expected_num_words : 
        raise ValueError("Line "+str(l)+ " is too long: " + line)                                
        
    elif len(split_string)<expected_num_words :
        raise ValueError("Line "+str(l)+ " is too short:" + line)
            
    
def range_to_list(range_string):
    split_range = [float(n) for n in range_string.split(",")]
    
    if len(split_range) > 1 :  
        step = split_range[2]           
        range_list =  np.arange(split_range[0],                          
                                 split_range[1]+step, 
                                 step)
    else : 
        step = None
        range_list = [split_range[0]]
    
    return range_list, step

def get_range_step(range_string):
    split_range = [float(n) for n in range_string.split(",")]
    return split_range[2]     

def arrange_data_line(in_val):
    
    # split_val = in_val.split(",")
    
    eps_to_test, eps_step = range_to_list(in_val[1])
    gamma_to_test, gamma_step = range_to_list(in_val[2])
        
    if (len(eps_to_test) > 1) and (len(gamma_to_test) > 1): 
        raise RuntimeError("Cannot vary both eps and gamma linearly.")
    
    # create arrays with repeated values, to cover every data point.
    # e.g. if eps_to_test = [0.01, 0.02, 0.03] and gamma_to_test = [15]
    # then eps_points = [0.01, 0.02, 0.03] and gamma_points = [15, 15, 15]
    # because all three eps values will be tested at gamma = 15º.
    eps_points = []
    gamma_points = []
    for e in range(len(eps_to_test)):
        for g in range(len(gamma_to_test)):
            eps_points.append(eps_to_test[e])
            gamma_points.append(gamma_to_test[g])
    
    if eps_step != None:
        return eps_points, gamma_points, eps_step
    elif gamma_step != None:
        return eps_points, gamma_points, gamma_step
    else:
        return eps_points, gamma_points, None

def arrange_mesh(in_val):
    
    eps_max = float(in_val[0])
    # the total number of data points = 0.5*outer_points*(outer_points + 1)
    outer_points = int(in_val[1])                                           
    
    # arrange a mesh of evenly distributed (eps,gamma) points to test over, 
    # such that eps_max is tested at (outer_points) gamma values. 
    # start from eps=0.001 because eps=0 is not an allowed input in asyrmo.
    eps_to_test = np.linspace(0.001, eps_max, num=outer_points)   
    
    eps_points = []
    gamma_points = []
    for e in range(len(eps_to_test)):
        for g in range(e+1):
            eps = np.round(eps_to_test[e], 3)
            eps_points.append(eps)
            if e==0: gamma = 0
            else: gamma = np.round(g*(60/e), 3) 
            gamma_points.append(gamma)
    eps_points = np.array(eps_points)
    gamma_points = np.array(gamma_points)
    
    return np.array(eps_points), np.array(gamma_points)






''' FUNCTIONS FOR RUNNING CODES '''

def write_orbitals(fermi_level, number, parity):
    
    first_index = fermi_level - number//2                  # select [num_orbs] orbitals (the fermi level plus [num_orbs//2] either side)
    last_index = fermi_level + number//2
    if number%2 == 1:                                                   # if an even number is requested, we need to add one to the final index to ensure the correct number are included
        last_index += 1
    orbitals = np.r_[first_index:last_index]                                        # generate a list of orbitals in unit steps inside this range with list slicing
    
    orbitals_string = parity+ str(number)                                       # e.g. orbitals_string = "4 19 20 21 22"
    
    for i in orbitals:
        orbitals_string += " "
        orbitals_string += str(i)
        
    return orbitals_string

def configure_script_writer(file_tags, nucleus, num_batches, 
                            num_per_batch, allowed_time):                       # these arguments will be the same throughout the execution, so we can define a family of script writers (one for each of gam/asy/prob)
    
    def run_script_batches(program):                                            # after configuring the script writer family, this is what we actually call (it's a "closure" function)
        
        if program == "gampn": abr = "GAM"
        elif program == "asyrmo": abr = "ASY"
        elif program == "probamo": abr = "PROB"
        else: raise ValueError("Input program ["+ program +"] not supported.")
        
        def write_script_batch(file_tag_batch, batch_index, nucleus, file_path):         # a sub-function, called in the loop below
            
            batch_folder = "Batch"+str(batch_index)
            script_text = ""
            
            for file in file_tag_batch :
                script_text += ("\n(cd ../"+nucleus+"/Run/"+batch_folder+
                                "; ./../../../Executables/MO/" + program + 
                                " < ../../Inputs/" +abr + "_" +file 
                                + ".DAT; cp "+program.upper()+ ".out ../../Outputs/"+
                                abr+"_"+file+".OUT)")
                                  # copy GAMPN.out to a new txt file with a more descriptive name
            
            script_text += ("\n\necho message from terminal: " +
                            "finished running " + program + " batch " + str(batch_index))
            
            script_file = open(file_path, 'w')
            script_file.write(script_text)
            script_file.close() 


        subprocesses = {}

        for b in range(num_batches):
            file_path = "../"+nucleus+"/Scripts/Run"+program.upper()+"_"+str(b+1)+".sh"
            batch_file_tags = file_tags[(b*num_per_batch):((b+1)*num_per_batch)]
            write_script_batch(batch_file_tags, b+1, nucleus, file_path)
            subprocesses[(program+"_"+str(b+1))] = subprocess.Popen(["sh", file_path])     # asynchronous call to start gampn as a subprocess

        for b in range(num_batches):
            subprocesses[(program+"_"+str(b+1))].wait(allowed_time)                        # wait to ensure it has finished before starting to read outputs, if it takes longer than the time limit seconds, throw an error to catch hangs.
            
    return run_script_batches



''' FUNCTIONS FOR READING GAMPN.OUT '''

def get_efac(lines):
     efac_line = lines.index("     KAPPA    MY     EPS   GAMMA    EPS4     EPS6     W0/W00   NMAX  COUPL     OMROT      EFAC      QFAC\n")
     return float(lines[efac_line+1][85:95].strip())
 
    
def get_fermi_level(lines, fermi_level_index):
    
    levels_header_line = lines.index("   #   ENERGY +/-(#)    <Q20>    <Q22>     <R2>     <JZ>      #   ENERGY +/-(#)    <Q20>    <Q22>     <R2>     <JZ>\n")
    fermi_level_line = fermi_level_index+levels_header_line+1               # calculate the line number of the fermi level in the GAMPN.OUT file (indexed from zero!)
    
    if fermi_level_index > 40:
        fermi_level_line -= 40
        whole_line = lines[fermi_level_line]
        half_line = whole_line[60:-1].strip()                                   # get only the second half of the line
    else:
        whole_line = lines[fermi_level_line]
        half_line = whole_line[0:60].strip()                                    # get only the first half of the line
    return half_line
 
def get_orbital_index(fermi_level_line, hash_index):
    single_parity_index = fermi_level_line[hash_index+1 : hash_index+3]
    if single_parity_index[1] == ")":                                           # in case the index is only a single digit, ignore the ")" that will have been caught
        single_parity_index = single_parity_index[0]
    return int(single_parity_index)




''' FUNCTIONS FOR READING PROBAMO.OUT '''

def read_data(line):
    line = line.strip()
    
    try:                                                                        # determine whether this line is a data row of the table (if ' - ' is present then it is)
        dash_index = line.index(" - ")
    except ValueError:
        return False # continue to next line
    
    spin_string = line[dash_index-4:dash_index].strip()
    final_spin_string = line[dash_index+11:dash_index+16].strip()               # get the spin of the final state (after transition)
    
    if not(spin_string == final_spin_string):
        return False

    spin_float = spin_string_to_float(spin_string)
    
    try:
        this_energy = float(line[:6].strip())
    except ValueError: # could not convert string to float: '0.0  1'
        this_energy = float(line[:5].strip())
    
    final_energy = float(line[dash_index+3:dash_index+10].strip())
                         
    if not(this_energy == final_energy):
        return False
    
    mag_moment = float(line[-8:].strip())
    line_data = {'spin_string': spin_string, 'spin_float':spin_float, 
                 'energy':this_energy, 'mag_moment':mag_moment}
    
    return line_data

def sort_by_spin(line_data, file_data):
    
    spin = "spin_"+line_data["spin_string"]
    
    if (spin+"_energies") in file_data:
        file_data[(spin+"_energies")].append(line_data["energy"])
        file_data[(spin+"_mag_moments")].append(line_data["mag_moment"])
        
    else:
        file_data[(spin+"_energies")] = [line_data["energy"]]
        file_data[(spin+"_mag_moments")] = [line_data["mag_moment"]]
    
    return file_data

def sort_by_expectation(line_data, file_data, inputs):
    
    if line_data["energy"] == 0.0: # record the ground state separately
        file_data["gs_spin_strings"] = line_data["spin_string"]
        file_data["gs_spin_floats"] = line_data["spin_float"]
        file_data["gs_mag_moments"] = line_data["mag_moment"]
    
    if "x1_spin" in inputs:
        if line_data["spin_string"] == inputs["x1_spin"]:
            file_data["x1_energies"] = line_data["energy"]
            file_data["x1_mag_moments"] = line_data["mag_moment"]
    
    if "x2_spin" in inputs:
        if line_data["spin_string"] == inputs["x2_spin"]:
            file_data["x2_energies"] = line_data["energy"]
            file_data["x2_mag_moments"] = line_data["mag_moment"]

    if "x3_spin" in inputs:
        if line_data["spin_string"] == inputs["x3_spin"]:
            file_data["x3_energies"] = line_data["energy"]
            file_data["x3_mag_moments"] = line_data["mag_moment"]
    
    return file_data
            
def missing_data(file_data, max_spin): 
    
    for i in range(max_spin):
        
        if i%2 == 0:
            continue # only half-int spins are calculated
        
        spin = "spin_"+str(i)+"/2"
        
        if not(spin+"_energies" in file_data):
            file_data[(spin+"_energies")] = [np.NaN]
            file_data[(spin+"_mag_moments")] = [np.NaN]

        
    return file_data
            
def restructure_data(old_data, max_spin):
# Take input data structured as a list of dictionaries. 
#       Each dictionary represents one deformation data-point.
#       Therefore the number of dictionaries is equal to the number of data points.
#       Each dictionary contains lists of states, organised by spin and property (energy or magnetic moment).
#       Therefore the number of lists in each dictionary is (1 + ((ISPIN+1)/2)*2) for gs magnetic moment and two properties.
# 
# Outputs the same data, reorganised into a new dictionary of lists.
#       Each list represents one property (energy or magnetic moment) and spin.
#       Therefore the number of lists is (1 + ((ISPIN+1)/2)*2) for gs magnetic moment and two properties.
#       Each list contains sub-lists of states, organised by deformation.
#       Therefore all the lists have the same length (equal to the number of data points).
#       This is more useful for plotting graphs.
    
    new_data = {}
    new_data["gs_spin_strings"] = []
    new_data["gs_spin_floats"] = []
    new_data["gs_mag_moments"] = []
    new_data["x1_mag_moments"] = []
    new_data["x2_mag_moments"] = []
    new_data["x3_mag_moments"] = []
    new_data["x1_energies"] = []
    new_data["x2_energies"] = []
    new_data["x3_energies"] = []
    

    for d in range(len(old_data)):
        new_data["gs_spin_strings"].append(old_data[d]["gs_spin_strings"])
        new_data["gs_spin_floats"].append(old_data[d]["gs_spin_floats"])
        new_data["gs_mag_moments"].append(old_data[d]["gs_mag_moments"])
        
        try:
            new_data["x1_energies"].append(old_data[d]["x1_energies"])
            new_data["x1_mag_moments"].append(old_data[d]["x1_mag_moments"])
        except(KeyError):
            print("Could not find any states with first excited spin in file " + str(d))
            new_data["x1_energies"].append(np.NaN)
            new_data["x1_mag_moments"].append(np.NaN)
        
        try:
            new_data["x2_energies"].append(old_data[d]["x2_energies"])
            new_data["x2_mag_moments"].append(old_data[d]["x2_mag_moments"])
        except(KeyError):
            print("Could not find any states with second excited spin in file " + str(d))
            new_data["x2_energies"].append(np.NaN)
            new_data["x2_mag_moments"].append(np.NaN)
        
        try:
            new_data["x3_mag_moments"].append(old_data[d]["x3_mag_moments"])
            new_data["x3_energies"].append(old_data[d]["x3_energies"])
        except(KeyError):
            print("Could not find any states with third excited spin in file " + str(d))
            new_data["x3_energies"].append(np.NaN)
            new_data["x3_mag_moments"].append(np.NaN)
    
    for i in range(max_spin): 
        
        if i%2 == 0:
            continue # only half-int spins are calculated
    
        spin = "spin_"+str(i)+"/2"
        
        new_data[spin+"_energies"] = []
        new_data[spin+"_mag_moments"] = []
        
        for d in range(len(old_data)):
           
            new_data[spin+"_energies"].append(old_data[d][spin+"_energies"])
            new_data[spin+"_mag_moments"].append(old_data[d][spin+"_mag_moments"])
            
    return new_data

def fill_gaps(multi_level_data):
    
    num_levels = [len(n) for n in multi_level_data]
    max_num = max(num_levels)
    
    for d in range(len(multi_level_data)):
        while len(multi_level_data[d]) < max_num:
            multi_level_data[d].append(np.NaN)
        
    return multi_level_data
    




''' FUNCTIONS FOR PLOTTING GRAPHS '''

def calc_contour_levels(data):
    
    min_contour = min(data)-0.5
    max_contour = max(data)+1.5
    return np.arange(min_contour, max_contour, 1.0)

def calc_cbar_tick_labels(data):
    cbar_ticks = np.arange(min(data), max(data)+1.0, 1.0)
    return [spin_float_to_string(n) for n in cbar_ticks]

def calc_cbar_ticks(contours):
    num = len(contours) - 1
    return [(contours[i] + contours[i+1]) / 2 for i in range(num)]                                # Calculate midpoints of levels for tick placement

def try_experimental(inputs, key, tolerance):
    try:
        return (inputs[key], tolerance)
    except KeyError:
        return (np.NaN, np.NaN)

def format_fig(polar_or_linear, ax, legend_handles, title, **kwargs):
    
    if polar_or_linear == 'polar':
        ax.set_thetamin(0)   
        ax.set_thetamax(60)  
        
        theta_ticks = np.arange(0, 70, 10)  
        ax.set_xticks(np.radians(theta_ticks))
        ax.yaxis.set_major_formatter(FuncFormatter(lambda x, _: f'{x:.2f}'))    # set the number of decimal places 
        
        plt.xlabel("ε")
        ax.text(45*np.pi/180, ax.get_rmax()*1.2, "γ", ha='center', va='center') # gamma axis label
        
        ax.legend(handles=legend_handles, loc="upper left", 
                           facecolor = 'lightblue', bbox_to_anchor=(-0.5, 1.1))
        
        ax.set_title(title, va='bottom', y=1.1)  
        
    elif polar_or_linear == 'linear':
        
        varied = kwargs.get("varied", None)
        x_label = kwargs.get("x_label", None)
        y_label = kwargs.get("y_label", None)
        legend_title = kwargs.get("legend_title", "")
    
        pad = 0.05*(varied[-1]-varied[0])
        ax.set_xlim([varied[0]-pad, varied[-1]+pad]) 
        
        plt.xlabel(x_label)
        plt.ylabel(y_label)
        
        ax.set_title(title, va='bottom', y=1.1) 
        
        legend = ax.legend(handles = legend_handles)
        
        if legend_title != "x":
            legend.set_title(legend_title)
    
    else: raise ValueError("unrecognised graph type: " + polar_or_linear + "; must be either 'polar' or 'linear'.")

def draw_contour_plot(ax, prop, data_points):
    
    if prop.sort == "Spin ":
        
        plot_data = np.transpose(prop.data)[0]
        
        print("multiple levels can't be plotted on a contour plot; \nplotting only the yrast state of spin " + prop.num)
        
        
    else:
        plot_data = np.array(prop.data)
        
    # some of the data points may be NaN, so manually create the triangulation and mask NaN triangles
    triang = tri.Triangulation(data_points["gamma_radians"], data_points["eps"])
    mask = np.any(np.isnan(plot_data[triang.triangles]), axis=1)
    triang.set_mask(mask)
    
    cax = ax.tricontourf(triang, plot_data, levels=prop.contour_levels)
    cbar = plt.colorbar(cax, pad=0.1, label=prop.axis_label)


    
    if prop.cbar_tick_labels:                                   # then format for discrete values
        
        cbar.set_ticks(prop.cbar_ticks)
        cbar.set_ticklabels(prop.cbar_tick_labels)
    
        
    return cax, cbar

def find_correct_spin(gs_spins, experimental_value):
    correct_spin_range = []                                             
    start_flag = False
    
    for i in range(len(gs_spins)):
        if gs_spins[i] == experimental_value and not start_flag:  # correct, and range hasn't started yet
            start_flag = True
            correct_spin_range.append(i)
            
        elif gs_spins[i] != experimental_value and start_flag:    # incorrect, and range has started
            start_flag = False
            correct_spin_range.append(i)
    
    return correct_spin_range

def plot_correct_spin(correct_spin_range, var, step, prop):
    
    for r in range(len(correct_spin_range)):
        
        if prop.sort == "Spin ":
            data = prop.data[0] 
            for i in range(1, len(prop.data)):
                data += prop.data[i]
        else:
            data = prop.data
            
        # front edge
        if correct_spin_range[r] == 0:                                  # the first value of eps in the range has the correct spin
            start_range = (var[correct_spin_range[r]]- step/2)          
        else:
            start_range = np.mean([var[correct_spin_range[r]-1], 
                                   var[correct_spin_range[r]]])
        
        correct_spin, = plt.plot([start_range, start_range], 
                                [min(data)-0.05*max(data), max(data)*1.05], 
                                'g-', label="range of correct spin")   
        # end edge
        if r%2==0:
            if r+1 == len(correct_spin_range):                          # the last value of eps in the range has the correct spin
                end_range = (var[-1]+step/2)
            else:
                end_range = np.mean([var[correct_spin_range[r+1]-1], 
                                     var[correct_spin_range[r+1]]])
       
        correct_spin, = plt.plot([end_range, end_range], 
                                [min(data)-0.05*max(data), max(data)*1.05], 
                                'g-', label="range of correct spin")  
    
        # top and bottom edges
        plt.plot([start_range, end_range],                          
                 [min(data)-0.05*max(data), 
                  min(data)-0.05*max(data)], 'g-')
        plt.plot([start_range, end_range], 
                 [max(data)*1.05, 
                  max(data)*1.05], 'g-')                  
        
    return correct_spin

def plot_points_with_experiment(data_points, prop, legend_handles, cbar):
    
    legend_hit = False
    legend_miss = False
    
    for r in range(len(data_points["eps"])):
        
        # use a smaller marker size for large data sets
        if len(data_points["file_tags"]) < 100: marker_size = 20
        else: marker_size = 5
        
        if prop.sort == "Spin ":
            error = [abs(prop.data[r][0] - exp) for exp in prop.experimental_data]
            match = [err < prop.error_tolerance for err in error]
        else:
            error = abs(prop.data[r] - prop.experimental_data)
            match = [error < prop.error_tolerance]
            
        if any(match): 
            data_points["agreed"][r] += 1                                       # record how many of the tested properties agree
            hit = plt.scatter(data_points["gamma_radians"][r], 
                      data_points["eps"][r], s=marker_size, edgecolor='red', 
                      facecolor='None', label="matches experiment")
            legend_hit = True
            
        # for data points that don't match experimental data, only plot them when the data set is quite small, to avoid cluttering the graph 
        elif len(data_points["file_tags"]) < 100: 
            miss, = plt.polar(data_points["gamma_radians"][r], 
                      data_points["eps"][r], 'wx', label="does not match experiment")
            legend_miss = True
    
    if legend_hit:
        legend_handles.append(hit)
    if legend_miss:
        legend_handles.append(miss)
    
    if prop.sort == "Spin ":
        for e in range(len(prop.experimental_data)):
            exp = cbar.ax.plot([0, 1], 
                               [prop.experimental_data[e], prop.experimental_data[e]], 
                               'r-', label = "experimental value")
    else:
        exp = cbar.ax.plot([0, 1], [prop.experimental_data, prop.experimental_data], 
                           'r-', label = "experimental value")
        
   
    legend_handles.append(exp[0])
    
    return legend_handles

def plot_points_without_experiment(data_points, legend_handles):
    
    # use a smaller marker size for large data sets
    if len(data_points["file_tags"]) < 100: marker_size = 5
    else: marker_size = 1
    
    all_points = plt.scatter(data_points["gamma_radians"], 
              data_points["eps"], s=marker_size, c='w', label="data point")
    legend_handles.append(all_points) 
    
    return legend_handles

def assign_parameters(inputs, e2plus_to_test, data_points):

    if inputs["deformation_input"] == "eps" :
        
        var_sym = "ε"
        var = data_points["eps"]
        fix_sym = "γ"
        fix = data_points["gamma_degrees"]
        
    elif inputs["deformation_input"] == "gamma" :
        
        var_sym = "γ / º"
        var = data_points["gamma_degrees"]
        fix_sym = "ε"
        fix = data_points["eps"]
        
    elif len(e2plus_to_test) > 1:
        
        var_sym = "E2PLUS / MeV"
        var = e2plus_to_test
        fix_sym = "(ε, γ)"
        fix = []
        
    else: raise ValueError("unrecognised graph request")

    return (var_sym, var, fix_sym, fix)

def plot_multi_lines(prop, var, legend_handles, marker_size, fix_sym, fix_val):
    data_by_line = np.transpose(prop.data)
    line_colours = ['k-x', 'b-x', 'y-x']
    line_labels = ["lowest energy", "second lowest energy", "third lowest energy"]
    
    for s in range(min(len(line_labels), np.size(data_by_line,0))):
        
        data, = plt.plot(var, data_by_line[s], line_colours[s], label=line_labels[s], markersize=marker_size)
        legend_handles.append(data)
        
    legend_title = "%s = %s" % (fix_sym, fix_val)
    
    return legend_handles, legend_title

def mark_spin(inputs, data_points, output_data, legend_handles, ax):

    correct_range = [inputs["gs_spin_float"]-0.5, 
                     inputs["gs_spin_float"]+0.5]
    spin_colour = (0,0,0) #(213/255,1,0)
    
    ax.tricontour(data_points["gamma_radians"], data_points["eps"], 
                  output_data["gs_spin_floats"].data, levels=correct_range,  
                  colors=[spin_colour], linewidths=1.0)
    spin_legend_proxy = plt.Line2D([], [], color=spin_colour, linewidth=1.0, label="region of correct g.s. spin") 
    legend_handles.append(spin_legend_proxy)

    return legend_handles





