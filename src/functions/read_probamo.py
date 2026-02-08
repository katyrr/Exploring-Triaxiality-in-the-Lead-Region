#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 21 11:47 2025

@author: katyrr

Functions for reading PROBAMO.OUT and processing contents.

"""

import numpy as np
import os

from src.functions.spin_processing import spin_string_to_float
from src.functions.parse_args import args

import src.functions.file_handling as fh



def read_probamo(num_points, data_subfolder_path, data_points, output_data, ptrm_inputs):
    '''
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
    (any DELTA==NaN values are bad, caused by some kind of convergence issue with BCS pairing).
    
    '''

    data_points["property_data"] = []
    
    for i in range(num_points):
        
        output_file_path = os.path.join(data_subfolder_path, "Outputs", f"PROB_{data_points["file_tags"][i]}.OUT")
        lines = fh.read_file(output_file_path)

        file_data = {}
        for l in lines:
            
            line_data = read_data(l)  # get the spin, energy, and magnetic moment from this line if it is a static moment, else return False
            if not(line_data):
                continue # to next line in file
            
            # sort line_data into file_data according to its spin
            file_data = sort_by_spin(line_data, file_data)
            # additionally save data that corresponds to the expected experimental ground state
            file_data = sort_by_expectation(line_data, file_data, ptrm_inputs)
        
        file_data = missing_data(file_data, ptrm_inputs)
        data_points["property_data"].append(file_data)
    
    return output_data

def process_data(output_data, data_points, experimental_data, ispin):
    #output_data = output_data_copy | rprob.restructure_data(data_points["property_data"], ptrm_inputs["ispin"], code_settings["print_details"])
    restructured_output_data = {**output_data, **restructure_data(data_points["property_data"], ispin)}

    # get energy gap between third 9/2 and first 13/2 states
    for i in experimental_data:
        if not "engap_" in i:
            continue
        
        spin1, idx1, spin2, idx2 = parse_engap_input(i)
        
        restructured_output_data[i] = find_gaps(restructured_output_data[f"spin_{spin1}/2_energies"], idx1, restructured_output_data[f"spin_{spin2}/2_energies"], idx2, experimental_data[i])

    # ensure all data sets have the same size and shape, and mask bad points

    mask = np.array([0 if np.isnan(x) else 1 for x in restructured_output_data["delta"]])

    for i in restructured_output_data:
        if isinstance(restructured_output_data[i][0], list):
            restructured_output_data[i] = fill_gaps(restructured_output_data[i])
            list_mask = np.transpose(np.tile(mask, (np.size(restructured_output_data[i][0]),1)))
            
            restructured_output_data[i] = np.where(list_mask == 0, np.NaN, restructured_output_data[i])
        
        else: 
            restructured_output_data[i] = np.where(mask == 0, np.NaN, restructured_output_data[i])
    
    return restructured_output_data



def find_gaps(spin_b_energies, index_b , spin_t_energies, index_t, exp):
    '''
    A function which calculates the energy gap between two energy levels of specified spins.
    Can specify which state of each spin to use (e.g. the lowest, or the second lowest, etc).
    Or can look at all states of the given spin, and use the one that returns an energy gap
    closest to the experimental value.
    

    Parameters
    ----------
    spin_b_energies : list of floats
        The calculated energies of the first specified spin (b).
        
    index_b : int
        Which calculated state of spin b to use (e.g 1 for the lowest state, 2 for the 
        second lowest, etc) - or 0 to look at all states of spin b.
        
    spin_t_energies : list of floats
        The calculated energies of the second specified spin (t).
        
    index_t : int
        Which calculated state of spin t to use (e.g 1 for the lowest state, 2 for the 
        second lowest, etc) - or 0 to look at all states of this spin t.
        
    exp : float
        Experimental value for the energy gap, in keV.

    Returns
    -------
    gaps : list of floats
        The calculated energy gap between the two states, at each deformation.

    '''

    # algorithm when index_b and index_t = 0:
    #   find the lowest spin b state
    #   check whether there are any spin t states at lower energy
    #   if not, NaN
    #   otherwise, get the energy gap between this spin t level and the closest spin b level below
    #   check whether this gap is smaller than the experimental gap 
    #   if so, AND there is another spin b level further down, get the gap to that level and compare with experiment 
    #   (if this is closer, keep it, otherwise return to the first value)
    #   otherwise, stop
    #   find the next spin t state
    #   get the best energy gap to its spin b below 
    #   if this gap is closer to experiment, keep 

    gaps = []
    
    for i in range(len(spin_t_energies)): # for each data point...
        
        gap = np.NaN
        gap_ref = np.inf
        
        if index_t == 0:
            all_t = spin_t_energies[i] # all the energy levels of spin t at this deformation
        else:
            try:
                all_t = [spin_t_energies[i][index_t-1]]
            except IndexError:
                gaps.append(np.NaN)
                continue
            
        if index_b == 0:
            all_b = spin_b_energies[i]
        else:
            try:
                all_b = [spin_b_energies[i][index_b-1]]
            except IndexError:
                gaps.append(np.NaN)
                continue
        
        for t in all_t:
            for b in all_b:
                
                if abs(t-b-exp) < gap_ref:
                    gap = t - b
                    gap_ref = abs(gap-exp)
                else:
                    continue
        gaps.append(gap)
        
    return gaps

def read_data(line):
    """
    A function which takes one line from the PROBAMO.OUT file, and checks 
    whether it describes a static state (internal transition). 
    
    If so, it reads the spin, energy, and magnetic dipole moment of the state.
    
    Otherwise it returns False and the line will be ignored.

    Parameters
    ----------
    line : string
        A single line read from the PROBAMO.OUT file.

    Returns
    -------
    line_data : dictionary
        A dictionary containing the state's spin (both as a fraction string 
        and as a float), energy, and magnetic dipole moment.
    
        If the line does not correspond to a static state, then returns False.
    
    """
    line = line.strip()
    
    # determine whether this line is a data row of the table 
    # (if ' - ' is present then it is).

    try:                                                                        
        dash_index = line.index(" - ")
    except ValueError:
        return False # continue to next line
    
    if not len(line)<100:
        return False # continue to next line
    
    # get the spins of the inital and final states of the transition.
    spin_string = line[dash_index-4:dash_index].strip()

    '''final_spin_string = line[dash_index+11:dash_index+16].strip()               
    if not(spin_string == final_spin_string):
        return False'''

    spin_float = spin_string_to_float(spin_string)
    
    #  get the energies of the inital and final states of the transition.
    try:
        this_energy = float(line[:6].strip())
    except ValueError: # could not convert string to float: '0.0  1'
        this_energy = float(line[:5].strip())
    
    '''final_energy = float(line[dash_index+3:dash_index+10].strip())    
    # determine whether the initial and final states are the same.
    if not(this_energy == final_energy):
        return False # continue to next line
    else:
        print(line)'''
    
    # read data.
    mag_moment = float(line[-8:].strip())
    quad_moment = float(line[-26:-16].strip())
    
    line_data = {'spin_string': spin_string, 'spin_float':spin_float, 
                 'energy':this_energy, 'mag_moment':mag_moment, 'quad_moment':quad_moment}
    
    return line_data



def sort_by_spin(line_data, file_data):
    """
    A function that sorts recorded data into groups based on spin.

    Parameters
    ----------
    line_data : dictionary
        Data about a single state of a single data point.
        Contains spin (both as a fractional string, and as a float), 
        energy in keV, and magnetic dipole moment in nuclear magnetons.
        
    file_data : dictionary
        A dictionary that will contain data about ALL the states of a single data point.
        May be empty or half-full when this function is called.
        Each entry in the dictionary will contain a list of values for a named property,
        e.g. "spin_1/2_energies" will be a list of energies for all states with 
        spin=1/2 in this file.

    Returns
    -------
    file_data : dictionary
        The same dictionary as was input, with data for one additional state now appended.

    """
    
    spin = "spin_"+line_data["spin_string"]
    
    if (spin+"_energies") in file_data:
        file_data[(spin+"_energies")].append(line_data["energy"])
        file_data[(spin+"_mag_moments")].append(line_data["mag_moment"])
        file_data[(spin+"_quad_moments")].append(line_data["quad_moment"])
        
    else:
        file_data[(spin+"_energies")] = [line_data["energy"]]
        file_data[(spin+"_mag_moments")] = [line_data["mag_moment"]]
        file_data[(spin+"_quad_moments")] = [line_data["quad_moment"]]
    
    return file_data


def sort_by_expectation(line_data, file_data, inputs):
    """
    A function that sorts recorded data into groups by association with input experimental data.
    If a calculated value for this experimental data has been recorded already,
    check whether the new value is closer, and only overwwrite if it is closer.

    Parameters
    ----------
    line_data : dictionary
        Data about a single state of a single data point.
        Contains spin (both as a fractional string, and as a float), 
        energy in keV, and magnetic dipole moment in nuclear magnetons.
        
    file_data : dictionary
        A dictionary that will contain data about ALL the states of a single data point.
        Will be half-full when this function is called.
        Each entry in the dictionary will contain the value of (/a list of values of) 
        a named property,
        e.g. "x1_energies" will be the energy for the states with (spin
        = experimental first excited state) in this file.
        
    inputs : dictionary
        A dictionary that contains name-value pairs for every input in the config file.
        May include some data about the experimental spectrum of the nucleus being studied,
        such as the ground state spin and magnetic moment, and the eneriges of the first
        few excited states.

    Returns
    -------
    file_data : dictionary
        The same dictionary as was input, with data for one additional state now appended.

    """
    
    if line_data["energy"] == 0.0: # record the ground state separately 
        file_data["gs_spin_strings"] = line_data["spin_string"]
        file_data["gs_spin_floats"] = line_data["spin_float"]
        file_data["gs_mag_moments"] = line_data["mag_moment"]
        file_data["gs_quad_moments"] = line_data["quad_moment"]
    
    
    
    # the rest of this function is obsolete
    
    if "x1_spin" in inputs: # assume that if x1_spin is input, then x1_energy will also have been input
        if line_data["spin_string"] == inputs["x1_spin"]:
            
            if ("x1_energies" in file_data): 
                # a lower energy state with the first excited spin has already been recorded
                existing_energy_error = abs(file_data["x1_energies"] - inputs["x1_energy"])
            else: existing_energy_error = np.inf
            
            new_energy_error = abs(line_data["energy"] - inputs["x1_energy"])
            
            if new_energy_error < existing_energy_error: # then overwrite
                
                file_data["x1_energies"] = line_data["energy"]
                file_data["x1_mag_moments"] = line_data["mag_moment"]
                file_data["x1_quad_moments"] = line_data["quad_moment"]
    
    if "x2_spin" in inputs:
        if line_data["spin_string"] == inputs["x2_spin"]:
        
            if ("x2_energies" in file_data): 
                # a lower energy state with the second excited spin has already been recorded
                existing_energy_error = abs(file_data["x2_energies"] - inputs["x2_energy"])
            else: existing_energy_error = np.inf
            
            new_energy_error = abs(line_data["energy"] - inputs["x2_energy"])
            
            if new_energy_error < existing_energy_error: # then overwrite
    
                file_data["x2_energies"] = line_data["energy"]
                file_data["x2_mag_moments"] = line_data["mag_moment"]
                file_data["x2_quad_moments"] = line_data["quad_moment"]

    if "x3_spin" in inputs:
        if line_data["spin_string"] == inputs["x3_spin"]:
            
            if ("x3_energies" in file_data): 
                # a lower energy state with the third excited spin has already been recorded
                existing_energy_error = abs(file_data["x3_energies"] - inputs["x3_energy"])
            else: existing_energy_error = np.inf
            
            new_energy_error = abs(line_data["energy"] - inputs["x3_energy"])
            
            if new_energy_error < existing_energy_error: # then overwrite
                file_data["x3_energies"] = line_data["energy"]
                file_data["x3_mag_moments"] = line_data["mag_moment"]
                file_data["x3_quad_moments"] = line_data["quad_moment"]
    
    return file_data
            

def missing_data(file_data, inputs): 
    """
    A function to fill missing data (i.e. data that was not calculated or not 
    found when reading PROBAMO.OUT).
    
    Ensures at least one state of each spin (from 1/2 up to ISPIN/2) has been 
    recorded, and fills gaps in energy and magnetic moment records with np.NaN.
    
    If experimental data has been provided, ensure that the closest energy match
    of the correct spin has been recorded, and if no states of that spin were
    found, fills the gap in energy and magnetic moment records with np.NaN.
    
    Assumes that a ground state will have been found (makes no checks).
    
    Parameters
    ----------
    file_data : dictionary
        A dictionary that will contain data about ALL the states calculated at a single data point.
        May be half-full when this function is called.
        Each entry in the dictionary will contain the value of (/a list of values of) 
        a named property.
        e.g. "x1_energies" will be the energy for the states with (spin 
        = experimental first excited state) in this file.
     
    inputs : dictionary
        Input settings from a config file, including the asyrmo input ISPIN.
        May contain some experimental data.

    Returns
    -------
    file_data : dictionary
        The same dictionary as was input, now with a full data set (any missing 
        data has been filled with np.NaN).


    """
    
    max_val = int(inputs["ispin"])+1
    
    for i in range(max_val): 
        
        if i%2 == 0:
            continue # only half-int spins are calculated
        
        spin = "spin_"+str(i)+"/2"
        
        if not(spin+"_energies" in file_data):
            file_data[(spin+"_energies")] = [np.NaN]
            # if the energy hasn't been recorded, then neither will the mag 
            # moment, and vice versa, because the code always outputs both.
            # Same for quadrupole moments
            file_data[(spin+"_mag_moments")] = [np.NaN]
            file_data[(spin+"_quad_moments")] = [np.NaN]

    if "x1_spin" in inputs: # assume that if x1_spin is input, then x1_energy will also have been input
        if not("x1_energies" in file_data): 
            file_data["x1_energies"] = np.NaN
            file_data["x1_mag_moments"] = np.NaN
    
    if "x2_spin" in inputs:
        if not("x2_energies" in file_data): 
            file_data["x2_energies"] = np.NaN
            file_data["x2_mag_moments"] = np.NaN
    
    if "x3_spin" in inputs:
        if not("x3_energies" in file_data): 
            file_data["x3_energies"] = np.NaN
            file_data["x3_mag_moments"] = np.NaN
        
    return file_data
            

def restructure_data(old_data, ispin):
    """
    Take input data structured as a list of dictionaries. 
    
    Outputs the same data, reorganised into a new dictionary of lists.
    This is more useful for plotting graphs.     

    Parameters
    ----------
    old_data : list of dictionaries
        Each dictionary represents one deformation data-point.
        Therefore the number of dictionaries is equal to the number of data points.
        Each dictionary contains lists of states, organised by spin and property 
        (energy or magnetic moment).
        Therefore the number of lists in each dictionary is (1 + ((ISPIN+1)/2)*2) 
        for gs magnetic moment and two properties.
        
    ispin : string
        The value of the asyrmo input ISPIN, from config.

    Returns
    -------
    new_data : dictionary of lists
        Contains the same data as the input, restructured into a new format.
        Each list represents one property (energy or magnetic moment) and spin.
        Therefore the number of lists is (1 + ((ISPIN+1)/2)*2) for gs magnetic moment and two properties.
        Each list contains sub-lists of states, organised by deformation.
        Therefore all the lists have the same length (equal to the number of data points).
        
    """

    new_data = {}
    new_data["gs_spin_strings"] = []
    new_data["gs_spin_floats"] = []
    new_data["gs_mag_moments"] = []
    new_data["x1_mag_moments"] = []
    new_data["x2_mag_moments"] = []
    new_data["x3_mag_moments"] = []
    new_data["gs_quad_moments"] = []
    new_data["x1_quad_moments"] = []
    new_data["x2_quad_moments"] = []
    new_data["x3_quad_moments"] = []
    new_data["x1_energies"] = []
    new_data["x2_energies"] = []
    new_data["x3_energies"] = []
    

    for d in range(len(old_data)):
        new_data["gs_spin_strings"].append(old_data[d]["gs_spin_strings"])
        new_data["gs_spin_floats"].append(old_data[d]["gs_spin_floats"])
        new_data["gs_mag_moments"].append(old_data[d]["gs_mag_moments"])
        new_data["gs_quad_moments"].append(old_data[d]["gs_quad_moments"])
        
        try:
            new_data["x1_energies"].append(old_data[d]["x1_energies"])
            new_data["x1_mag_moments"].append(old_data[d]["x1_mag_moments"])
            new_data["x1_quad_moments"].append(old_data[d]["x1_quad_moments"])
        except(KeyError):
            if args.verbose: 
                print("Could not find any states with first excited spin in file " + str(d))
            new_data["x1_energies"].append(np.NaN)
            new_data["x1_mag_moments"].append(np.NaN)
            new_data["x1_quad_moments"].append(np.NaN)
        
        try:
            new_data["x2_energies"].append(old_data[d]["x2_energies"])
            new_data["x2_mag_moments"].append(old_data[d]["x2_mag_moments"])
            new_data["x2_quad_moments"].append(old_data[d]["x2_quad_moments"])
        except(KeyError):
            if args.verbose: 
                print("Could not find any states with second excited spin in file " + str(d))
            new_data["x2_energies"].append(np.NaN)
            new_data["x2_mag_moments"].append(np.NaN)
            new_data["x2_quad_moments"].append(np.NaN)
        
        try:
            new_data["x3_mag_moments"].append(old_data[d]["x3_mag_moments"])
            new_data["x3_quad_moments"].append(old_data[d]["x3_quad_moments"])
            new_data["x3_energies"].append(old_data[d]["x3_energies"])
        except(KeyError):
            if args.verbose: 
                print("Could not find any states with third excited spin in file " + str(d))
            new_data["x3_energies"].append(np.NaN)
            new_data["x3_mag_moments"].append(np.NaN)
            new_data["x3_quad_moments"].append(np.NaN)
    
    
    max_val = int(ispin)+1

    for i in range(1, max_val, 2): 
    
        spin = "spin_"+str(i)+"/2"
        
        new_data[spin+"_energies"] = []
        new_data[spin+"_mag_moments"] = []
        new_data[spin+"_quad_moments"] = []
        
        for d in range(len(old_data)):
            
            new_data[spin+"_energies"].append(old_data[d][spin+"_energies"])
            new_data[spin+"_mag_moments"].append(old_data[d][spin+"_mag_moments"])
            new_data[spin+"_quad_moments"].append(old_data[d][spin+"_quad_moments"])
            
    return new_data


def fill_gaps(multi_level_data):
    """
    A function that fills gaps in multi-level data. 
    (e.g. spin_1/2_energies: there may be more than one spin 1/2 level calculated 
     for each data point, but not all data points will necessarily calculate the 
     same number of spin 1/2 levels. This function fills any gaps with np.NaN).

    Parameters
    ----------
    multi_level_data : a matrix (2D np.array or list) of floats
        Number of columns = number of data points (i.e. number of deformations being tested).
        Number of rows = number of levels of a certain spin calculated in each file (may be inconsistent).

    Returns
    -------
    multi_level_data : a matrix (2D np.array) of floats
        The same as input, now with any gaps filled with np.NaN, such that 
        every column has the same number of rows (a regular rectangular matrix).

    """
    
    num_levels = [len(n) for n in multi_level_data]
    max_num = max(num_levels)
    
    for d in range(len(multi_level_data)):
        while len(multi_level_data[d]) < max_num:
            multi_level_data[d].append(np.NaN)
        
    return multi_level_data
    

def parse_engap_input(engap_input):
    '''
    Get the spins and indices being requested in energy gap inputs.
    e.g. the config input "engap_9.3_13.1" is asking for the third 9/2 state
    and the first 13/2 state.

    Parameters
    ----------
    engap_input : string
        The input variable name for energy gap experimental data.

    Returns
    -------
    first_spin : string
        The numerator of the fractional spin, e.g. for 9/2 this is "9".
        
    first_index : int
        The index of the state, e.g. for the third state this is 3.
        
    second_spin : string
        The numerator of the fractional spin, e.g. for 13/2 this is "13".
        
    second_index : int
        The index of the state, e.g. for the first state this is 1.

    '''
    first_dot = engap_input.index(".")
    first_spin = engap_input[first_dot-2:first_dot]
    if first_spin[0] == "_":
        first_spin = first_spin[1]
    first_index = int(engap_input[first_dot+1:first_dot+2])
    
    second_dot = engap_input.index(".", first_dot+1)
    second_spin = engap_input[second_dot-2:second_dot]
    if second_spin[0] == "_":
        second_spin = second_spin[1]
    second_index = int(engap_input[second_dot+1:second_dot+2])
    
    return first_spin, first_index, second_spin, second_index
            