#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 21 11:43 2025

@author: katyrr

Functions for reading config file and processing contents.

"""

import numpy as np
import math
import os

import src.functions.structs as st
import src.functions.file_handling as fh
import src.functions.run_ptrm as ptrm

from src.functions.spin_processing import spin_string_to_float


#--------------------------------------------------------------------------------------------------

def read_config(data_subfolder_path, code_settings, ptrm_inputs, data_points, experimental_data, graphs_to_plot):
    '''
    Locate and read the config file (if not found, create a new one from template).
    Read each line, and save inputs (correctly type-cast) as name-value pairs in dictionaries.
    Inputs are checked for validity and any missing required inputs.
    Format and (re)calculate any additional required parameters from the config inputs.
    
    
    Parameters
    ----------
    data_subfolder_path: string
        The directory path to the subfolder where the config file is 
        located and calculations will be done.

    code_settings : dictionary
        A dictionary containing settings for automating the running of the codes.

    ptrm_inputs : dictionary
        A dictionary that contains name-value pairs for ptrm inputs in the config file.

    data_points : dictionary
        A dictionary that contains lists of the variable inputs (deformations, etc)
        
    experimental_data : dictionary 
        A dictionary containing experimental data input via config.

    graphs_to_plot: dictionary
        A dictionary containing settings for graph plotting, input via config. 

    Returns
    -------
    None (dictionaries are edited in-place)
    '''

    config_path = fh.locate_config(data_subfolder_path)
    config_lines = fh.read_file(config_path)

    read_lines(config_lines, code_settings, ptrm_inputs, data_points, experimental_data, graphs_to_plot)
    validate_inputs(code_settings, ptrm_inputs, data_points, experimental_data, graphs_to_plot, data_subfolder_path)
    process_inputs(data_points, ptrm_inputs)

    code_settings["num_points"] = len(data_points["eps"])
    print("Number of data points = ", code_settings["num_points"])
    print("Deformation range:")
    print(f"\teps = [{data_points["eps"][0]:.3f}, {data_points["eps"][-1]:.3f}]")
    print(f"\tgamma = [{data_points["gamma_degrees"][0]:.1f}, {data_points["gamma_degrees"][-1]:.1f}] degrees")


    
#--------------------------------------------------------------------------------------------------

def read_lines(lines, code_settings, ptrm_inputs, data_points, experimental_data, graphs_to_plot):
    '''
    - Ignores empty lines, and lines beginning with * (to mark a comment).
    - Checks that the format of each line is correct ("var_name value").
    - Saves inputs via helper functions depending on the type of the input (inferred from the name)
    - All inputs are saved as name-value pairs in dictionaries.

    Parameters
    ----------
    lines: list of strings
        The lines read from the config file.
    
    code_settings : dictionary
        A dictionary containing settings for automating the running of the codes.

    ptrm_inputs : dictionary
        A dictionary that contains name-value pairs for ptrm inputs in the config file.

    data_points : dictionary
        A dictionary that contains lists of the variable inputs (deformations, etc)
        
    experimental_data : dictionary 
        A dictionary containing experimental data input via config.

    graphs_to_plot: dictionary
        A dictionary containing settings for graph plotting, input via config. 

    Returns
    -------
    None (dictionaries are mutated in place)

    Raises
    ------
    ValueError
        Occurs if the number of 'words' in a line does not match the expected number.
        Reports the original (unedited) line and the line number.
    
    '''

    for i in range(len(lines)):                       
        
        # remove/ignore comments and blank lines, and check formatting  
        line = lines[i].strip()
                                    
        if line == "" : continue                  
        if line[0] == "*" : continue              
        
        split_string = line.split(" ")  # split into name and value
        split_string = remove_inline_comments(split_string)

        if not check_line_format(split_string):
            raise ValueError(f"Incorrect number of values given for the parameter on line {i+1}: {line}")

        name = split_string[0]       
        
        match name:
            case "eps"|"gamma"|"single"|"mesh": 
                save_deformation_input(ptrm_inputs, data_points, split_string)
            case "e2plus": 
                save_e2plus_input(ptrm_inputs, data_points, split_string)
            case "gs_spin":
                save_gs_spin_input(experimental_data, split_string)
            case name if name[:3]=="jp_": 
                experimental_data[name] = [float(n) for n in split_string[1].split(',')]
            case name if name[:5]=="plot_":
                graphs_to_plot[name[5:]] = bool(int(split_string[1]))
            case name if name[:6]=="engap_":
                experimental_data[name] = float(split_string[1]) 
            case _:
                save_typecast_input(code_settings, ptrm_inputs, experimental_data, split_string)

        # print(f"{name}: \t {split_string[1:]}")


def validate_inputs(code_settings, ptrm_inputs, data_points, experimental_data, graphs_to_plot, data_subfolder_path):
    '''
    Some inputs can only take certain values (e.g. OS = "MacOS" or "64bit").
    This function checks that those inputs have a valid value.
    It also checks that all required settings have been input.

    Parameters
    ----------
    ptrm_inputs : dictionary
        A dictionary that contains name-value pairs for every input in the config file.

    data_points : dictionary
        A dictionary that contains lists of the variable inputs (deformations, etc)
        
    experimental_data : dictionary 
        A dictionary containing experimental data input via config.

    Raises
    ------
    ValueError
        Occurs if the input is not one of the valid values.
    
    RuntimeError
        Occurs if a required input is missing.

    Returns
    -------
    None

    '''
    all_inputs = {**ptrm_inputs, **experimental_data, **code_settings, **graphs_to_plot, **data_points}

    for i in st.get_required_inputs():
        if not i in all_inputs:
            raise RuntimeError(f"missing input: {i}")

    if 0.0 in data_points["eps"]:
        raise ValueError("Cannot test at eps=0.0.")
    
    restricted_inputs = st.get_restricted_inputs()
    
    for name in all_inputs:
        if name in restricted_inputs:
            allowed_values = restricted_inputs[name]

            if not all_inputs[name] in allowed_values:
                raise ValueError(f"Invalid input: \t {name} = {all_inputs[name]}.\nPlease choose from allowed values: {str(allowed_values)}")
    
    if str(ptrm_inputs["A"]) not in ptrm_inputs["nucleus"]:
        raise ValueError(f"The input value of A ({ptrm_inputs["A"]}) does not match the value in the name of the nucleus ({ptrm_inputs["nucleus"]}). Please fix the incorrect one in the config file.")

    if ptrm_inputs['nucleus'] not in data_subfolder_path:
        print(f"\nWARNING: the name of the nucleus being studied ({ptrm_inputs['nucleus']}) does NOT")
        print(f"\tappear in the name of the data subfolder ({os.path.basename(data_subfolder_path)}).")
        print("\tDid you make a typo?\n")

def process_inputs(data_points, ptrm_inputs):
    '''
    - Calculate and store gamma values in radians as well as degrees.
    - Convert the nantj, noutj, ipout inputs to the correct format.

    - Using input A and Z, work out which particle is odd, to determine the nneupr input.
    - Halve and ceiling for the overall index of the fermi level orbital.


    - Raises ValueError if an even-mass nucleus is input.

    Parameters
    ----------
    ptrm_inputs : dictionary
        A dictionary that contains name-value pairs for ptrm inputs in the config file.

    data_points : dictionary
        A dictionary that contains lists of the variable inputs (deformations, etc)
        
    
    Returns
    -------
    None (dictionaries are edited in-place)

    '''
    
    data_points["gamma_radians"] = [n*np.pi/180 for n in data_points["gamma_degrees"]]

    ptrm_inputs["nantj"] = ptrm_inputs["nantj"].replace(",", " ")
    ptrm_inputs["noutj"] = ptrm_inputs["noutj"].replace(",", " ")
    ptrm_inputs["ipout"] = ptrm_inputs["ipout"].replace(",", " ")

    ptrm_inputs["N"] = ptrm_inputs["A"]-ptrm_inputs["Z"]

    print(f"\nNucleus Data:")
    print(f"\tName: {ptrm_inputs['nucleus']}")
    print(f"\tA: {ptrm_inputs['A']}")
    print(f"\tZ: {ptrm_inputs['Z']}")
    print(f"\tN: {ptrm_inputs['N']}")

    if ptrm_inputs["A"]%2 == 0:
        raise ValueError("Input nucleus is even-A. Only odd-mass nuclei accepted.")
    elif ptrm_inputs["Z"]%2 == 0: 
        ptrm_inputs["nneupr"] = "-1" 
        ptrm_inputs["fermi_level"] = math.ceil(ptrm_inputs["N"]/2)
        print("\tOdd Neutrons\n")                 
    elif ptrm_inputs["N"]%2 == 0:
        ptrm_inputs["nneupr"] = "1"
        ptrm_inputs["fermi_level"] = math.ceil(ptrm_inputs["Z"]/2)
        print("\tOdd Protons\n")
    else:
        raise RuntimeError("Check inputs of A and Z.")
    
    ptrm_inputs["current_orbitals"] = ptrm.write_orbitals(
        ptrm_inputs["fermi_level"]//2, 
        ptrm_inputs["num_orbs"], 
        ptrm_inputs["par"]
    )
    
    
    
#--------------------------------------------------------------------------------------------------


def remove_inline_comments(split_string):
    """
    A function to remove inline comments from a line of text read from the config file.
    If the character "*" appears in the line, it marks the remainder of the line
    as a comment which should be ignored.

    Parameters
    ----------
    split_string : list of strings
        A single line of text, split into words by delimeter " ".

    Returns
    -------
    split_string : list of strings
        The same as the input split_string, minus any inline comments.

    """
    for n in range(len(split_string)):                                          
        word = split_string[n]
        if word[0] == '*':
            return split_string[:n]
        
    return split_string

def check_line_format(split_string):
    """
    A function to check the how many values are associated with this config input,
    compare it to the expected number, and raise an error if it is unexpected.
    
    Most inputs have expected length = 2 (the name of the variable and its value)
    but some have expected length = 3 (the name of the variable and two values).
    
    e.g. "Z" should have only one value (the proton number of the nucleus)  
    whereas "single" should have two values (the eps value and the gamma value 
    of the data point).

    Parameters
    ----------
    split_string : list of strings
        A single line of text, split into words by delimeter " ", 
        with any inline comments already removed.

    Returns
    -------
    False if the format was incorrect, otherwise True.

    """
    
    if (split_string[0]=="single" 
        or split_string[0]=="eps" 
        or split_string[0]=="gamma"): 
        
        expected_num_words = 3     
        
    else: expected_num_words = 2                                                
    
    if len(split_string) != expected_num_words : 
        return False
                                      
    return True

def save_deformation_input(ptrm_inputs, data_points, split_string):
    '''
    A function that reads the deformation input line of the config file, to determine
    what kind of deformation input has been made, and to arrange lists of eps and gamma
    values to test. Also records the step size, for linear inputs.
    

    Parameters
    ----------
    ptrm_inputs : dictionary
        A dictionary that contains name-value pairs for every input in the config file.

    data_points : dictionary
        A dictionary that contains lists of the variable inputs (deformations, etc)
        
    split_string : list of strings
        A single line of text, split into words by delimeter " ".

    Raises
    ------
    RuntimeError
        Occurs if multiple deformation inputs are made.

    Returns
    -------
    None (dictionaries modified in-place)

    '''
    if ("deformation_input" in ptrm_inputs):
        raise RuntimeError("Deformation has already been input: " + ptrm_inputs["deformation_input"])
    
    ptrm_inputs["deformation_input"] = split_string[0]

    # save data point arrays
    if split_string[0]=="mesh":
        data_points["eps"], data_points["gamma_degrees"] = arrange_mesh(split_string[1].split(","))
    else: 
        data_points["eps"], data_points["gamma_degrees"] = arrange_data_line(split_string)
        
    # save steps
    if split_string[0]=="eps":
        ptrm_inputs["step"] = get_range_step(split_string[1])
    elif split_string[0] == "gamma":
        ptrm_inputs["step"] = get_range_step(split_string[2])
        
def save_e2plus_input(ptrm_inputs, data_points, split_string):
    '''
    Reads E2PLUS input, arranges linear array if required, and stores in dictionary.

    Parameters
    ----------
    ptrm_inputs : dictionary
        A dictionary that contains name-value pairs for every input in the config file.

    data_points : dictionary
        A dictionary that contains lists of the variable inputs (deformations, etc)
        
    split_string : list of strings
        A single line of text, split into words by delimeter " ".

    Raises
    ------
    RuntimeError
        Occurs if deformation is not input (or if e2plus is input before deformation).
        Also if a range of E2PLUS is input with a non-point deformation.

    Returns
    -------
    None (dictionaries are modified in-place)

    '''
    
    if "eps" not in data_points:
        raise RuntimeError("missing deformation input (or perhaps deformation was input below e2plus? make sure deformation is input first.)")
    
    num_eps = len(data_points["eps"])

    if split_string[1]=="0":
        # if e2plus has been input with value = "0", then it will later be 
        # calculated dynamically based on the deformation of each data point.
        data_points["e2plus"] = np.zeros((num_eps,), dtype=int)
        return
    
    data_points["e2plus"], _ = range_to_list(split_string[1])
    num_e2plus = len(data_points["e2plus"])
    
    if (num_e2plus>1 and ptrm_inputs["deformation_input"] != "single"):
        raise RuntimeError("Testing a range of e2plus is only supported for a single deformation input.")
    
    elif (num_e2plus==1 and num_eps>1):
        data_points["e2plus"] = [data_points["e2plus"][0]]*num_eps
            
    else: 
        data_points["eps"] = data_points["eps"] * num_e2plus
        data_points["gamma_degrees"] = data_points["gamma_degrees"] * num_e2plus

def save_gs_spin_input(experimental_data, split_string):
    '''
    A function that reads the config line which states the experimental ground state spin.
    The formatting is checked, converted to float, and both string and float versions are recorded.

    Parameters
    ----------
    experimental_data : dictionary 
        A dictionary containing experimental data input via config.
        
    split_string : list of strings
        A single line of text, split into words by delimeter " ".

    Raises
    ------
    ValueError
        Occurs if the spin is not input in the format 'n/2' where n is an (odd) integer.

    Returns
    -------
    None (dictionary is modified in-place)

    '''
    
    try:
        experimental_data["gs_spin_float"] = spin_string_to_float(split_string[1])
    except ValueError:
        raise ValueError("wrong format for input of gs_spin, please input in the format '1/2' or '13/2', etc.")

    experimental_data["gs_spin_string"] = split_string[1]

def save_typecast_input(code_settings, ptrm_inputs, experimental_data, split_string):
    '''
    Converts input settings to the correct type before saving them in the relevant dictionary.

    Parameters
    ----------
    code_settings : dictionary
        A dictionary containing settings for automating the running of the codes.

    ptrm_inputs : dictionary
        A dictionary that contains name-value pairs for every input in the config file.
        
    experimental_data : dictionary 
        A dictionary containing experimental data input via config.

    split_string : list of strings
        A single line of text, split into words by delimeter " ".

    Raises
    ------
    ValueError
        Occurs if the type of the input has not been specified.

    Returns
    -------
    None

    '''

    name = split_string[0]
    
    match name:
        case "OS":
            code_settings[name] = split_string[1]

        case "figure_res" | "num_cores" : 
            code_settings[name] = int(split_string[1])

        case name if name in st.get_variable_list("bool"):  
            code_settings[name] = bool(int(split_string[1]))

        case name if name in st.get_variable_list("int"):
            ptrm_inputs[name] = int(split_string[1])

        case name if name in st.get_variable_list("experimental_float"):
            experimental_data[name] = float(split_string[1])

        case name if name[:3] in st.get_variable_list("experimental_float"):
            experimental_data[name] = float(split_string[1])

        case name if name in st.get_variable_list("settings_float"):
            ptrm_inputs[name] = float(split_string[1])

        case name if name in st.get_variable_list("string"): 
            ptrm_inputs[name] = split_string[1]

        case _ : raise ValueError(f"unrecognised input: {name}")

    
#--------------------------------------------------------------------------------------------------


def range_to_list(range_string):
    """
    A function to convert a string input range-and-step to a list of explicit float values.

    Parameters
    ----------
    range_string : string
        The input range-and-step, in the format "start,end,step".
        The start and end values are inclusive.
        A single value can also be input, in the format "value".

    Returns
    -------
    range_list : np.array of floats.
        A list of values in the input range, with the input step.
        If a single value was input, then this list will have length 1.
        
    step : float
        The step applied to the range.
        If a single value was input, then step=None.

    """
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
    """
    A function to get a float step value out of a string input range-and-step.

    Parameters
    ----------
    range_string :  string
        The input range-and-step, in the format "start,end,step".

    Returns
    -------
    step : float
        The step to be applied to the range.

    """
    
    split_range = [float(n) for n in range_string.split(",")]
    step = split_range[2]
    return step   


def arrange_data_line(split_string):
    """
    A function to convert ranges of eps and gamma values-to-be-tested into arrays 
    of data points, one of which has repeated values such that both arrays have 
    the same length = number of data points to be tested.
    
    e.g. if eps_to_test = [0.01, 0.02, 0.03] and gamma_to_test = [15]
    then eps_points = [0.01, 0.02, 0.03] and gamma_points = [15, 15, 15]
    because all three eps values will be tested at gamma = 15º.

    Parameters
    ----------
    split_string : list of strings
        A single line of text, split into words by delimeter " ", 
        with any inline comments already removed.
        
        Has length = 3, with the first value being the name 
        ("eps", or "gamma", or "single") to indicate which if any is varied,
        the second value being eps (a range or a single value), 
        and the third value being gamma (a range or a single value). 

    Raises
    ------
    RuntimeError
        Occurs if both eps and gamma are input as ranges. Only one at a time 
        can be varied linearly.

    Returns
    -------
    eps_points : np.array of floats
        Contains the eps value to test at each data point.
        Has length = number of data points.
        
    gamma_points : np.array of floats
        Contains the gamma value to test at each data point.
        Has length = number of data points.

    """
    
    eps_to_test, _ = range_to_list(split_string[1])
    gamma_to_test, _ = range_to_list(split_string[2])
        
    if (len(eps_to_test) > 1) and (len(gamma_to_test) > 1): 
        raise RuntimeError("Cannot vary both eps and gamma linearly.")

    eps_points = []
    gamma_points = []
    for e in range(len(eps_to_test)):
        for g in range(len(gamma_to_test)):
            eps_points.append(eps_to_test[e])
            gamma_points.append(gamma_to_test[g])
    
    return eps_points, gamma_points


def arrange_mesh(split_string):
    """
    A function to convert a given number of data points and maximum eps value  
    into arrays of data points, evenly distributed in polar coordinates,
    such that both arrays have the same length = number of data points to be tested.
    
    The total number of data points is determined by the input outer_points as 
    = 0.5*outer_points*(outer_points + 1).
    
    gamma is tested in the range [0, 60]º so that both oblate and prolate 
    deformations are covered, with only positive values of eps in the range 
    [0.001, eps_max] (cannot start at 0.0 because the asyrmo code will hang for 
    perfectly spherical inputs).

    Parameters
    ----------
    split_string : list of strings
        A single line of text, split into words by delimeter " ", 
        with any inline comments already removed.
        
        Has length = 3, with the first value being the name ("mesh"), 
        the second value being eps_max, and the third value being outer_points.
        
        outer_points is defined as the number of data points along one outer edge 
        of the wedge in polar coordinates (i.e. the number of gamma points 
        tested at eps_max, or the number of unique eps values tested).

    Returns
    -------
    eps_points : np.array of floats
        Contains the eps value to test at each data point.
        Has length = number of data points.
        
    gamma_points : np.array of floats
        Contains the gamma value to test at each data point.
        Has length = number of data points.

    """
    
    eps_max = float(split_string[0])
    outer_points = int(split_string[1])                                           
    
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
    
    return eps_points, gamma_points


    
