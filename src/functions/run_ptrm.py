#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 21 11:47 2025

@author: katyrr

Functions for preparing to run the PTRM codes.

"""

import os
import numpy as np
import subprocess
import math

import functions.file_handling as fh
import functions.structs as st

from functions.read_gampn import get_sp_level, get_info


def write_input_files(num_points, data_subfolder_path, program, ptrm_inputs, data_points, first_run=False):
    '''
    Create .DAT input files for each data point, specialised to the requested program.

    If this is the first program to be run (gampn for the first time) then some extra
    calculations must be done first:
        - Dynamically calculate e2plus for each data point (unless a fixed value has been specified)
        - Generate a file tag for each data point, referencing the nucleus, deformation, and e2plus 
          of that calculation
    
    Then loops through the data points:
        - Determine which properties are required for the requested program.
        - Select the correct value of those properties for the current data point.
        - Uses string formatting with the ptrm_inputs dictionary to write the .DAT file.
        - Creates/overwrites a .DAT file in the Inputs directory folder, named using the file tag.

    inputs
    ------
    num_points : int
        The number of data points being calculated 
    
    data_subfolder_path : string
        The directory path to the data subfolder where these calculations are to be done.
    
    program : string
        Either "gampn", "asyrmo", or "probamo"
    
    ptrm_inputs : dict
        A dictionary of values to be input to the ptrm programs.
    
    data_points : dict
        A dictionary of lists of eps, gamma, and e2plus values.
        (The length of each list equals num_points)
    
    first_run : bool
        Whether or not this is the first program to be run (gampn for the first time).
        Default = False. 
        If True then some additional setup calculations are performed.
    ''' 

    if first_run:

        for i in range(num_points):

            if data_points["e2plus"][i] == 0:
                # if the value of e2plus has been input as 0, calculate dynamically
                data_points["e2plus"][i] = est_e2plus(data_points["eps"][i], ptrm_inputs["A"])

            file_tag = set_current(i, ptrm_inputs, data_points, set_deformations=True, set_e2plus=True, create_file_tag=True)
            data_points["file_tags"].append(file_tag)


    for i in range(num_points):

        match program:
            case "gampn": 
                file_tag = set_current(i, ptrm_inputs, data_points, "002", "016", "017", 
                               set_deformations=True, set_orbitals= not first_run)
            case "asyrmo":
                file_tag = set_current(i, ptrm_inputs, data_points, "016", "017", "018", 
                               set_e2plus=True, set_orbitals=True)
            case "probamo":
                file_tag = set_current(i, ptrm_inputs, data_points, "017", "018")
            case _:
                raise ValueError(f"unrecognised program: {program}")

        
        abr = get_program_abr(program)
        file_path = os.path.join(data_subfolder_path, "Inputs", f"{abr}_{file_tag}.DAT")
        fh.write_file(file_path, st.get_template(program) % ptrm_inputs)


def set_current(i, ptrm_inputs, data_points, *args, 
                set_deformations=False, set_e2plus=False, create_file_tag=False, set_orbitals=False):
    '''
    A function which sets the current values of the deformation parameters, E2PLUS input, 
    and binary file names. For use in the "run gampn" loops.

    Parameters
    ----------
    i : int
        The index for the iteration of the loop.

    inputs : dictionary
        A dictionary that contains name-value pairs for every input in the config file.
        
    data_points : dictionary
        A dictionary that contains lists of the variable inputs (deformations, etc)
        
    *args : strings
        Optionally input any number of ptrm in/out binary file codes, which will be 
        named using the current file tag.
        Select from [002, 016, 017, 018]

    set_deformations : bool 
        (default False)
        Whether or not to set the current deformation parameters (current_eps and current_gamma_degrees) 
        Only necessary for gampn
        
    set_e2plus : bool
        (default False)
        Whether or not to set the current e2plus value
        Only necessary for asyrmos
    
    create_file_tag : bool
        (default False)
        Whether or not to generate the file tag from scratch 
        Only necessary on the first run (otherwise the existing tag is fetched)
    
    set_orbitals : bool
        (default False)
        Whether or not to set the current orbital input string 
        Only necessary for gampn and asyrmo

    Returns
    -------
    file_tag : string
        format "e[eps]_g[gamma]_p[e2plus]_[nucleus]",
        e.g. "e0.001_g10.0_p0.730_Pb207"
        for use in naming files.

    '''
    
    if set_deformations:
        if data_points["gamma_degrees"][i] > 30:
            ptrm_inputs["current_eps"] = data_points["eps"][i] * -1
            ptrm_inputs["current_gamma"] = 60 - data_points["gamma_degrees"][i]
        else:
            ptrm_inputs["current_eps"] = data_points["eps"][i]
            ptrm_inputs["current_gamma"] = data_points["gamma_degrees"][i]

    if set_e2plus:
        ptrm_inputs["current_e2plus"] = data_points["e2plus"][i]

    if set_orbitals:
        ptrm_inputs["current_orbitals"] = data_points["asyrmo_orbitals"][i]

    if create_file_tag:
        file_tag = f"e{ptrm_inputs['current_eps']:.3f}_g{ptrm_inputs['current_gamma']:.1f}_p{ptrm_inputs['current_e2plus']:.3f}_{ptrm_inputs["nucleus"]}"
    else:
        file_tag = data_points["file_tags"][i]

    for a in args:
        ptrm_inputs[f"current_f{a}"] = f"f{a}_{file_tag}.dat"

    return file_tag

    

def write_orbitals(fermi_level, number, parity):
    """
    A function to generate a parity + orbital string for input into gampn or asyrmo.
    e.g. orbitals_string = "+4 19 20 21 22" for 4 orbitals [19, 20, 21, 22] of positive parity.
    

    Parameters
    ----------
    fermi_level : int
        The index of the Fermi level orbital.
        Assumes that orbitals are indexed separately for positive and negative parities.
        This will be the central level in the orbital string.

        
    number : int
        The number of orbitals to include in the string.
        
    parity : string
        "+" for positive parity;
        "-" for negative parity;

    Returns
    -------
    orbitals_string : string
        The parity, number of orbitals, and list of orbitals, 
        formatted ready for input into gampn or asyrmo.
        e.g. orbitals_string = "+4 19 20 21 22".

    """
    
    # select [num_orbs] orbitals (the fermi level plus [num_orbs//2] either side)
    first_index = fermi_level - number//2                                       
    last_index = fermi_level + number//2
    
    # if an even number is requested, we need to add one to the final index to 
    # ensure the correct number are included.
    if number%2 == 1: last_index += 1                                                      
    
    # generate a list of orbitals in unit steps inside this range with list slicing
    orbitals = np.r_[first_index:last_index]                                    
   
    
    orbitals_string = parity+ str(number)                                       
    
    for i in orbitals:
        orbitals_string += " "
        orbitals_string += str(i)
        
    return orbitals_string


def find_orbitals(fermi_level, number, parity, fermi_energy, fermi_parity, lines):
    """
    A function to dynamically generate a parity + orbital string for input into asyrmo.
    e.g. orbitals_string = "+4 19 20 21 22" for 4 orbitals [19, 20, 21, 22] of positive parity.
    
    Chooses the orbitals nearest in energy to the fermi level, which may be unbalanced
    (e.g. may not be the fermi level ± 5, but could be the fermi level + 2 - 8).
    

    Parameters
    ----------
    fermi_level : int
        The index of the Fermi level orbital.
        Assumes that orbitals are indexed separately for positive and negative parities.
        
    number : int
        The number of orbitals to include in the string.
        
    parity : string
        "+" for positive parity;
        "-" for negative parity;

    Returns
    -------
    orbitals_string : string
        The parity, number of orbitals, and list of orbitals, 
        formatted ready for input into asyrmo.
        e.g. orbitals_string = "+4 19 20 21 22".

    """
    
    # this works, but... it could probably be more efficient. #!!!
    # also still need to implement checking that the requested orbital was included in the gampn input.
    
    if fermi_parity == parity:
        orbitals_list = [fermi_level]
        
        # next above/below is fermi_level ± 1
        
        
    else:
        orbitals_list = []
        
    # next above/below:
        
    fermi_line = get_sp_level(lines, fermi_level, fermi_parity)
    overall_index = int(fermi_line[0:2].strip())
    
    index_below = overall_index - 1
    index_above = overall_index + 1
    
    while len(orbitals_list) < number :
        if index_below < 1:
        
            line_above = get_sp_level(lines, index_above, '0')
            parity_above, energy_above, level_above = get_info(line_above)
            
            # can't add any more orbitals below, so fill up with orbitals above
            while parity_above != parity:
               index_above += 1
               if index_above > 80:
                   break
               
               line_above = get_sp_level(lines, index_above, '0')
               parity_above, energy_above, level_above = get_info(line_above)
               
            orbitals_list.append(level_above)
            index_above += 1
               
        
        elif index_above > 80:
            # can't add any more orbitals above, so fill up with orbitals below
            line_below = get_sp_level(lines, index_below, '0')
            parity_below, energy_below, level_below = get_info(line_below)
            
            while parity_below != parity:
               index_below -= 1
               if index_below < 1:
                   break
               
               line_below = get_sp_level(lines, index_below, '0')
               parity_below, energy_below, level_below = get_info(line_below)
               
            orbitals_list.append(level_below)
            index_below -= 1
            
        else:
            line_below = get_sp_level(lines, index_below, '0')
            parity_below, energy_below, level_below = get_info(line_below)
            
            line_above = get_sp_level(lines, index_above, '0')
            parity_above, energy_above, level_above = get_info(line_above)
            
            while parity_below != parity:
               index_below -= 1
               if index_below < 1:
                   break
               
               line_below = get_sp_level(lines, index_below, '0')
               parity_below, energy_below, level_below = get_info(line_below)
               
            while parity_above != parity:
               index_above += 1
               if index_above > 80:
                   break
               
               line_above = get_sp_level(lines, index_above, '0')
               parity_above, energy_above, level_above = get_info(line_above)
                
            lower_energy_gap = fermi_energy - energy_below
            upper_energy_gap = energy_above - fermi_energy
            
            if lower_energy_gap < upper_energy_gap:
                orbitals_list.append(level_below)
                index_below -= 1
            else:
                orbitals_list.append(level_above)
                index_above += 1
    
    orbitals_string = parity + str(number) + ' ' + ' '.join(str(x) for x in sorted(orbitals_list))
    
    return orbitals_string

def est_e2plus(eps, A):
    """
    A function to dynamically calculate the effective E2PLUS value of a single 
    data point, using Grodzin's relation. 
    
    eps is assumed to be interchangable with beta, to the precision of this estimate,
    as reccomended on p6 of the manual, from which the formula is taken:
        
    E2PLUS approx = 1100 / [pow(beta, 2) * pow(A, 7/3)] MeV
    
    #!!! the same p6 of the manual mentions a reference to a more sophisticated model...?
    

    Parameters
    ----------
    eps : float
        The value of eps with which to calculate E2PLUS.
        
    A : int
        The mass number of the nucleus.

    Returns
    -------
    e2plus : float
        The estimated effective value of e2plus, in MeV.

    """
    denom = pow(eps, 2)*pow(A, 7/3)
    e2plus = np.round(1225 / denom, 3)  # MeV
    
    return e2plus


def get_batch_settings(num_cores, num_points):
    '''
    Calculate batch settings (how to divide up the data points into batches).
        - If the data set is small then use fewer cores for a minimum batch size of 20, 
          to make the overhead worthwhile.
        - The maximum allowed time for the batch to run before assuming that it is hanging.
            - Each file takes ~ 0.1 seconds to run;
            - Allow double time plus an overhead/extra of 10 seconds to ensure that the batch 
              will finish even if the computer is running a bit slow today! If it takes longer 
              than this, it is probably hanging, but it does occasionallyl take longer... 
              You could increase the overhead allowance from 10 to 30 seconds, or more, depending
              on the machine.
    '''

    batch_settings = {}

    batch_settings["num_batches"] = num_cores # = number of cores for maximum efficiency with large data sets
    batch_settings["num_per_batch"] = math.ceil(num_points/batch_settings["num_batches"])

    if batch_settings["num_per_batch"] < 20:
        batch_settings["num_per_batch"] = 20
        batch_settings["num_batches"] = math.ceil(num_points/batch_settings["num_per_batch"])  

    batch_settings["allowed_time"] = 0.2*batch_settings["num_per_batch"]+10   

    return batch_settings

def configure_script_writer(folder_path, OS, batch_settings, file_tags):
    """
    A closure for configuring a general script writer, which can then be customised to 
    each program (gampn, asyrmo, probamo) while maintaining a consistent strategy 
    for dividing up file batches.

    Parameters
    ----------
    folder_path : string
        The absolute path to the folder containing the config file.

    OS : either "MacOS" or "64bit"
        Which version of the pre-compiled PTRM Fortran codes to use, depending 
        on your computer's operating system (Mac, or Windows/Linux)

    batch_settings : dict
        A dictionary containing information on how to divide batches.
        Includes keys: num_batches, num_per_batch, and allowed_time.

    file_tags : list of strings
        One tag for each data point, with format "e[eps]_g[gamma]_p[e2plus]_[nucleus]",
        e.g. "e0.001_g10.0_p0.730_Pb207".
        The list has length = number of data points.

    Returns
    -------
    run_script_batches(program) : function
        A function that runs the input "program", with data points divided up into batches.

        For each batch of data points:
            - Moves the working directory to the batch folder.
            - Gets the list of file tags corresponding to that batch.
            - Writes a bash script to execute the program for all the files in the batch.
            - Starts the bash script as a subprocess.
        
        When all the scripts have been set running, the function waits for them 
        all to complete before returning.
        
        Parameters
        ----------
        program : string
            "gampn" to run gampn
            "asyrmo" to run asyrmo
            "probamo" to run probamo

        Raises
        ------
        ValueError
            Occurs if a program other than gampn, asyrmo, or probamo is requested.

        Returns
        -------
        None.
    
    """

    def run_script_batches(program):

        print(f"running {program} in folder {folder_path}")
        match program:
            case "gampn": abr = "GAM"
            case "asyrmo": abr = "ASY"
            case "probamo": abr = "PROB"
            case _: raise ValueError(f"unrecognised program: {program}")

        if OS == "64bit":
            program = program.upper() + ".exe"

        subprocesses = {}
        num_per_batch = batch_settings["num_per_batch"]

        program_path = os.path.join("src", "ptrm", OS, "MO", program)
        abs_program_path = os.path.abspath(program_path)
        
        for b in range(batch_settings["num_batches"]):
            batch_file_tags = file_tags[(b*num_per_batch):((b+1)*num_per_batch)]
            
            file_path = os.path.join(folder_path, "Scripts", f"Run{program.upper()}_{b+1}.sh")
            run_folder_path = os.path.join(folder_path, "Run", f"Batch{b+1}")

            script_text = "" # "pwd" to print the working directory for debugging

            for file in batch_file_tags:
                
                # define input/output file paths relative to what WILL BE the working directory at runtime (run_folder_path)
                input_file_path = os.path.join(os.pardir, os.pardir, "Inputs", f"{abr}_{file}.DAT")
                output_file_path = os.path.join(os.pardir, os.pardir, "Outputs", f"{abr}_{file}.OUT")

                # move the working directory to run_folder_path and start a subprocess in that folder
                # call the program with the input file
                # copy the default output file to a new .OUT file with a more descriptive name in the Outputs folder.
                script_text += f"\n(cd {run_folder_path};{abs_program_path} < {input_file_path}; cp {program.upper()}.out {output_file_path})"

            #script_text += f"\n\necho message from terminal: finished running {program} batch {b+1}"
            
            script_file = open(file_path, 'w')
            script_file.write(script_text)
            script_file.close() 
            
            subprocesses[f"{program}_{b+1}"] = subprocess.Popen(["sh", file_path])   

            # asynchronous call to start the program as a subprocess
            
        for b in range(batch_settings["num_batches"]):
            # wait to ensure it has finished (before starting to read outputs!), 
            # if it takes longer than the time limit seconds, throw an error to catch hangs.
            subprocesses[f"{program}_{b+1}"].wait(batch_settings["allowed_time"])


    return run_script_batches


def get_program_abr(program):

    match program:
            case "gampn": return "GAM"
            case "asyrmo": return "ASY"
            case "probamo": return "PROB"
            case _: raise ValueError(f"unrecognised program: {program}")