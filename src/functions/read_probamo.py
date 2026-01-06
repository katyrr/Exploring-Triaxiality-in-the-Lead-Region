#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 21 11:47 2025

@author: katyrr

Functions for reading PROBAMO.OUT and processing contents.

"""

import numpy as np
from functions.spin_processing import spin_string_to_float

import functions.structs as st

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
            

def restructure_data(old_data, ispin, verbose):
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
        
    verbose : bool
        True to print high detail messages to console.
        False to print only essential information to console.

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
            if verbose: 
                print("Could not find any states with first excited spin in file " + str(d))
            new_data["x1_energies"].append(np.NaN)
            new_data["x1_mag_moments"].append(np.NaN)
            new_data["x1_quad_moments"].append(np.NaN)
        
        try:
            new_data["x2_energies"].append(old_data[d]["x2_energies"])
            new_data["x2_mag_moments"].append(old_data[d]["x2_mag_moments"])
            new_data["x2_quad_moments"].append(old_data[d]["x2_quad_moments"])
        except(KeyError):
            if verbose: 
                print("Could not find any states with second excited spin in file " + str(d))
            new_data["x2_energies"].append(np.NaN)
            new_data["x2_mag_moments"].append(np.NaN)
            new_data["x2_quad_moments"].append(np.NaN)
        
        try:
            new_data["x3_mag_moments"].append(old_data[d]["x3_mag_moments"])
            new_data["x3_quad_moments"].append(old_data[d]["x3_quad_moments"])
            new_data["x3_energies"].append(old_data[d]["x3_energies"])
        except(KeyError):
            if verbose: 
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
    



def collate_energy_data(output_data, num_points, expected_gs_spin_string, experimental):
    """
    A function that takes all of the separate spin_n/2_energies data sets
    and combines them into a single collection. 
    
    Additionally makes a note of the energies of the lowest state of the expected 
    ground state spin, for future use. Saved as a list in the "shifts" property.
    

    Parameters
    ----------
    output_data : list of PropertyData objects
        Contains all the data read from the output files, packaged with graph
        plotting information.
        
    num_points : int
        The number of data points in the set (i.e. the number of files or 
        the number of deformations).

    Returns
    -------
    all_energy_levels : PropertyData object
        A collection of energy level data for all spins, packaged with
        graph plotting information.

    """
    all_level_data = np.zeros((num_points,1))
    spins = []
    
    for d in output_data:
        if output_data[d].prop == "energies" and output_data[d].sort == "Spin ":
            this_level_data = output_data[d].data
            all_level_data = np.hstack((all_level_data, this_level_data))
            spins += [output_data[d].num]*np.size(this_level_data,1)
            
            if output_data[d].num == expected_gs_spin_string:
                gs_spin_energies = this_level_data[:,0]  
                # the energies of the lowest state of the expected gs spin
            
            
    all_level_data = np.delete(all_level_data, [0], axis=1)
    
    indices = np.argsort(all_level_data[0,:])
    all_level_data = all_level_data[:, indices]
    spins = np.array(spins)
    spins = spins[indices]
    
    experimental_data = []
    expspins = []
    
    for e in experimental:
        if e[0:3] == "jp_":
            for i in experimental[e]:
                experimental_data.append(i)
                expspins.append(e[3:-1])
            
    
    all_energy_levels = st.PropertyData(all_level_data, "All Energy Levels")
    all_energy_levels.contour_levels = 10
    all_energy_levels.cbar_ticks = 0
    all_energy_levels.cbar_tick_labels = 0
    all_energy_levels.experimental_data = experimental_data
    all_energy_levels.explabels = expspins
    all_energy_levels.error_tolerance = np.NaN
    all_energy_levels.spins = [spin_string_to_float(n) for n in spins]
    all_energy_levels.shifts = gs_spin_energies
            
    return all_energy_levels

def shift_energy_levels(default_energy_levels):
    '''
    A function which takes a collated list of all energy levels, and recalculates
    them relative to the lowest state with the expected ground state spin.

    Parameters
    ----------
    default_energy_levels : PropertyData object
        All of the calculated energy levels, wrapped up in a PropertyData object.
        Contains a "shifts" property, which is a list of the energies of the lowest
        state with the expected ground state spin. These are used to perform the shift.

    Returns
    -------
    shifted_energies :  PropertyData object
        All energies, calculated relative to the expected ground state, and wrapped
        in a PropertyData object.

    '''
    
    shifted_energy_levels = []
    
    # for each deformation:
    for d in range(np.size(default_energy_levels.data, 0)):
        
        this_deformation_data = default_energy_levels.data[d,:]
        this_deformation_gs_energy = default_energy_levels.shifts[d]
        
        # shift energies to be relative to this
        shifted_energy_levels.append([(e - this_deformation_gs_energy) for e in this_deformation_data])
    
    shifted_energies = st.PropertyData(np.array(shifted_energy_levels), "All Energy Levels (Relative)")
    shifted_energies.contour_levels = 10
    shifted_energies.cbar_ticks = 0
    shifted_energies.cbar_tick_labels = 0
    shifted_energies.experimental_data = default_energy_levels.experimental_data
    shifted_energies.explabels = default_energy_levels.explabels
    shifted_energies.error_tolerance = np.NaN
    shifted_energies.spins = default_energy_levels.spins
            
    return shifted_energies


   
def calc_rms_err(range_min, *props):
    '''
    Calculate the root mean square error in the input properties. Any number of
    properties can be input. Technically they should all have the same units though.
    This function works on the assumption that they are energy levels in keV.
    
    For energy level data with more than one state, the "correct" energy level 
    is assumed to be the lowest - this is the one that is compared with the 
    experimental value.
    
    If any of the input properties are missing experimental data, they are 
    ignored.
    
    
    Assumes

    Parameters
    ----------
    range_min : int
        The lowest rms error in keV to put on the colour bar. It will use a log
        scale, so to highlight the best values you should put approx the minimum
        calculated RMS as the lowest value. Usually 10 keV is good, but for very 
        good calculations you could go down to 2 keV.
    
    *props : PropertyData objects
        Any number of properties, for which the rms error will be calculated.

    Returns
    -------
    rms : PropertyData object
        An object containing an array of the calculated rms errors at each deformation,
        as well as graph plotting settings.

    '''
    
    rms_data = np.zeros(np.size(props[0].data, axis=0)) # start a record of the rms error (for each data point) at zero
    num_nan = np.zeros(np.size(props[0].data, axis=0), dtype=int) # start a counter for the number of NaN values (for each data point) at zero
    
    
    for a in props:
        
        data = a.data
        exp = a.experimental_data
        
        if np.isnan(exp):
            continue   # ignore any properties with no experimental data
        
        if np.size(data, axis=1)>1:
            data = np.transpose(data)[0]
        
        if isinstance(exp, list):
            exp = exp[0] # use the first experimental value if more than one was input
        
        
        for d in range(len(data)):
            
            # count NaN values
            if np.isnan(data[d]):
                num_nan[d] += 1
                continue
        
            rms_data[d] += ((data[d]-exp)/1000)**2 # calculate in MeV to avoid overflow (the numbers get big quickly!)
    
    for d in range(len(rms_data)):
        rms_data[d] = np.sqrt( rms_data[d] / (np.size(props[0].data, axis=1)-num_nan[d])) * 1000 # convert back to keV at the end
    
    if range_min == 10:
        contour_levels = [10, 20, 30, 40, 50, 100, 200, 300, 400, 500,1000]
        cbar_ticks = [10, 20, 50, 100, 200, 500, 1000]
        cbar_tick_labels = [str(x) for x in cbar_ticks] 
    else:
        contour_levels = [2,5, 10, 20, 30, 40, 50, 100, 200, 300, 400, 500,1000]
        cbar_ticks = [2, 5, 10, 20, 50, 100, 200, 500, 1000]
        cbar_tick_labels = [str(x) for x in cbar_ticks] 
        
        
    rms = st.PropertyData(np.array(rms_data), "RMS energies")
    rms.range_min = range_min
    rms.contour_levels = contour_levels
    rms.cbar_ticks = cbar_ticks
    rms.cbar_tick_labels = cbar_tick_labels
    rms.experimental_data = np.NaN
    rms.error_tolerance = np.NaN
    
    return rms
