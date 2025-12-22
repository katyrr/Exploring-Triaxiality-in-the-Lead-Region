#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 21 11:43 2025

@author: katyrr

Functions for reading config file and processing contents.

"""
import os
import shutil
import numpy as np

import functions.structs as st

from functions.spin_processing import spin_string_to_float

def remove_inline_comments(split_string, line_index):
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
        
        if len(word) == 0:
            #raise ValueError("extra space on line = " + str(line_index+1))
            continue # this error doesn't actually matter, it just means there's accidentally a double space somewhere, we can ingore it
        
        if word[0] == '*':
            return split_string[:n]
    return split_string



def check_line_format(split_string, line, l):
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
        
    line : string
        The full text that was read from this line (unedited).
    
    l : int
        The line number of this line in the config file.

    Raises
    ------
    ValueError
        Occurs if the number of words does not match the expected number.
        Reports the original (unedited) line and the line number.

    Returns
    -------
    None.

    """
    
    if (split_string[0]=="single" 
        or split_string[0]=="eps" 
        or split_string[0]=="gamma"): 
        
        expected_num_words = 3     
        
    else: expected_num_words = 2                                                
    
    if len(split_string)>expected_num_words : 
        raise ValueError("Line "+str(l+1)+ " is too long: " + line)                                
        
    elif len(split_string)<expected_num_words :
        raise ValueError("Line "+str(l+1)+ " is too short:" + line)
            
    
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
    
    eps_to_test, eps_step = range_to_list(split_string[1])
    gamma_to_test, gamma_step = range_to_list(split_string[2])
        
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



def save_deformation_input(inputs, data_points, split_string):
    '''
    A function that reads the deformation input line of the config file, to determine
    what kind of deformation input has been made, and to arrange lists of eps and gamma
    values to test. Also records the step size, for linear inputs.
    

    Parameters
    ----------
    inputs : dictionary
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
    inputs : dictionary
        A dictionary that contains name-value pairs for every input in the config file.
        Now contains a new entry.
        
    data_points : dictionary
        A dictionary that contains lists of the variable inputs (deformations, etc).
        Now contains a new entry.

    '''
    if ("deformation_input" in inputs):
        raise RuntimeError("Deformation has already been input: " + inputs["deformation_input"])
    
    inputs["deformation_input"] = split_string[0]
    
    if split_string[0]=="mesh":
        data_points["eps"], data_points["gamma_degrees"] = arrange_mesh(split_string[1].split(","))
        
    else: 
        data_points["eps"], data_points["gamma_degrees"] = arrange_data_line(split_string)
        
    if split_string[0]=="eps":
        inputs["step"] = get_range_step(split_string[1])
    elif split_string[0] == "gamma":
        inputs["step"] = get_range_step(split_string[2])
        
    return inputs, data_points

def save_e2plus_input(inputs, data_points, split_string):
    '''
    

    Parameters
    ----------
    inputs : dictionary
        A dictionary that contains name-value pairs for every input in the config file.

    data_points : dictionary
        A dictionary that contains lists of the variable inputs (deformations, etc)
        
    split_string : list of strings
        A single line of text, split into words by delimeter " ".

    Raises
    ------
    RuntimeError
        Occurs if deformation is not input (or if e2plus is input before deformation).
        Also if the input deformation is not just a single point.

    Returns
    -------
    inputs : dictionary
        A dictionary that contains name-value pairs for every input in the config file.
        Now contains a new entry.
        
    data_points : dictionary
        A dictionary that contains lists of the variable inputs (deformations, etc).
        Now contains a new entry.

    '''
    
    if "eps" not in data_points:
        raise RuntimeError("missing deformation input (or perhaps deformation was input below e2plus?)")
        
    if split_string[1]=="0":
        # if e2plus has been input with value = "0", then it will later be 
        # calculated dynamically based on the deformation of each data point.
        data_points["e2plus"] = np.zeros((len(data_points["eps"]),), dtype=int)
    
    else:
        data_points["e2plus"], i = range_to_list(split_string[1])
        
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

    return inputs, data_points

def save_gs_spin_input(experimental, split_string):
    '''
    A function that reads the config line which states the experimental ground state spin.
    The formatting is checked, converted to float, and both string and float versions are recorded.

    Parameters
    ----------
    experimental : dictionary 
        A dictionary containing experimental data input via config.
        
    split_string : list of strings
        A single line of text, split into words by delimeter " ".

    Raises
    ------
    ValueError
        Occurs if the spin is not input in the format 'n/2' where n is an (odd) integer.

    Returns
    -------
    experimental : dictionary 
        A dictionary containing experimental data input via config.

    '''
    
    try:
        experimental["gs_spin_float"] = spin_string_to_float(split_string[1])
    except ValueError:
        raise ValueError("wrong format for input of gs_spin, please input in the format '1/2' or '13/2', etc.")

    experimental["gs_spin_string"] = split_string[1]
        
    return experimental


def validate_input(inputs, experimental, split_string):
    '''
    Some inputs can only take certain values (e.g. OS = "MacOS" or "64bit").
    This function checks that those inputs have a valid value.

    Parameters
    ----------
    inputs : dictionary
        A dictionary that contains name-value pairs for every input in the config file.
        
    experimental : dictionary 
        A dictionary containing experimental data input via config.

    split_string : list of strings
        A single line of text, split into words by delimeter " ".

    Raises
    ------
    ValueError
        Occurs if the input is not one of the valid values.

    Returns
    -------
    inputs : dictionary
        A dictionary that contains name-value pairs for every input in the config file.
        
    experimental : dictionary 
        A dictionary containing experimental data input via config.

    '''
    
    if split_string[0] in st.get_variable_list("int"):
        
        split_string[1] = int(split_string[1])
        dictionary = inputs
        
    elif (split_string[0] in st.get_variable_list("experimental_float")
          or split_string[0][:3] in st.get_variable_list("experimental_float")):
        
        split_string[1] = float(split_string[1])
        dictionary = experimental
        
    elif split_string[0] in st.get_variable_list("settings_float"):
        
        split_string[1] = float(split_string[1])
        dictionary = inputs
                                  
    elif split_string[0] in st.get_variable_list("bool"):  
                                  
        split_string[1] = bool(int(split_string[1]))
        dictionary = inputs
    
    elif split_string[0] in st.get_variable_list("string"): 
        
        dictionary = inputs                                  
        
    else: raise ValueError("unrecognised input: " + split_string[0])

    restricted_inputs = st.get_restricted_inputs()
    
    if split_string[0] in restricted_inputs:
        allowed_values = restricted_inputs[split_string[0]]

        if not split_string[1] in allowed_values:
            raise ValueError("Invalid input: \t" + split_string[0] + " = " + split_string[1] + ".\nPlease choose from allowed values: " + str(allowed_values))
            
    dictionary[split_string[0]] = split_string[1]
    
    
    return inputs, experimental
