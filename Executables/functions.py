#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  7 10:53:52 2025

@author: katyrr

A separate module file to contain function definitions.

CONTENTS:
--------

- GENERAL FUNCTIONS

- FUNCTIONS FOR READING CONFIG

- FUNCTIONS FOR RUNNING CODES

- FUNCTIONS FOR READING GAMPN.OUT
        
- FUNCTIONS FOR READING PROBAMO.OUT

- FUNCTIONS FOR ANALYSING RESULTS

"""

import numpy as np                              # for np.arrays
import subprocess                               # for calling shell scripts to run 
import structs as st                            # my own module file of structs (classes, and read-only dicts)

#%%

''' GENERAL FUNCTIONS '''

def spin_string_to_float(spin_string):
    """
    A function to covert the fractional (string) representation of a nuclear 
    spin to its float representation.

    Parameters
    ----------
    spin_string : string
        A fractional represention of a half-int spin.
        e.g. "1/2", "3/2", etc.

    Raises
    ------
    ValueError
        Occurs if the spin_string is not input in the correct format "n/2" where n is an integer.

    Returns
    -------
    spin_float : float
        The float representation of the spin. 
        e.g. 0.5, 1.5, etc.

    """
    if spin_string[1] == "/":
        spin_float = float(spin_string[0])/2
    elif spin_string[2] == "/":
        spin_float = float(spin_string[0:2])/2 
    else:
        raise ValueError("Cannot parse spin. \nCheck it has been input in " +
                         "the format 'n/2' where n is an int.")
    return spin_float


def spin_float_to_string(spin_float):
    """
    A function to convert a float representation of a nuclear spin 
    to its fractional representation, as a string.

    Parameters
    ----------
    spin_float : float
        The float representation of the spin. 
        e.g. 0.5, 1.5, etc.

    Raises
    ------
    ValueError
        Occurs if the ipnut spin_float cannot be converted to a half-integer fraction.
        i.e. must be n.5 where n is an integer.

    Returns
    -------
    spin_string : string
        A fractional represention of a half-int spin.
        e.g. "1/2", "3/2", etc.

    """
    numerator = float(spin_float)*2
    if str(numerator)[-1] != "0":
        raise ValueError("Cannot convert float to fraction. \nCheck it has " + 
                         "been input as a half-integer float.")
    else:
        spin_string = str(int(numerator))+"/2"
        return spin_string
    
    
    
def read_file(path):
    """
    A function to open a file, then read and return its full contents. 

    Parameters
    ----------
    path : string
        The relative file path from the current working directory to the file being opened.

    Returns
    -------
    lines : list of strings
        The contents of the file, where each element of the list is one line of the file.

    """
    with open(path, 'r') as f:
        lines = f.readlines()
    
    return lines



def write_file(path, text):
    """
    A function to open a file, and write some text to it. If the file already 
    exists, it will be overwritten. If it does not exist yet, it will be created.

    Parameters
    ----------
    path : string
        The relative file path from the current working directory to the file being opened.
        Includes the name of the file itself at the end of the file path.

    text : string
        The text that will be written to the file. 

    Returns
    -------
    None.

    """
    
    with open(path, 'w') as f:
        f.write(text)



#%%


''' FUNCTIONS FOR READING CONFIG '''

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




#%%

''' FUNCTIONS FOR RUNNING CODES '''

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

def configure_script_writer(file_tags, nucleus, num_batches, num_per_batch, allowed_time, verbose): 
    """
    A closure to configure a general script writer, which can then be customised to 
    each program (gampn, asyrmo, probamo) while maintaining a consistent strategy 
    for dividing up file batches.
    

    Parameters
    ----------
    file_tags : list of strings
        One tag for each data point, with format "e[eps]_g[gamma]_p[e2plus]_[nucleus]",
        e.g. "e0.001_g10.0_p0.730_Pb207".
        The list has length = number of data points.

    nucleus : string
        The nucleus being tested, with format e.g. "Pb207". 
        Used to navigate the file directory.
        
    num_batches : int
        The number of batches to divide data points into.
        
    num_per_batch : int
        The maximum number of data points per batch. 
        If the total number of data points is not exactly divisible by the 
        number of batches, then the last batch may contain fewer data points.
        
    allowed_time : float
        The maximum time in seconds to allow for the batch to run. 
        If the runtime exceeds this then the program is assumed to be hanging, 
        and execution is halted.
    
    verbose : bool
        True to print high detail messages to console
        False to print only essential information to console

    Returns
    -------
    run_script_batches(program) : function
        A function that runs the input "program", with data points divided up into batches.
    
    """
    def run_script_batches(program):
        """
        The inner function of the configure_script_writer closure. This is the 
        function that is returned from the closure, which can then be called.
        
        Dynamically defines a script writer customised to the program requested.
        
        For each batch of data points:
            - Determines the file path to the folder from which the batch will be run.
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
        if program == "gampn": abr = "GAM"
        elif program == "asyrmo": abr = "ASY"
        elif program == "probamo": abr = "PROB"
        else: raise ValueError("Input program ["+ program +"] not supported.")
        
        def write_script_batch(file_tag_batch, batch_index, nucleus, file_path): 
            """
            A sub-function dynamically defined to write bash scripts customised 
            to the program requested.
            
            The bash script executes the following method for each file in the batch:
                - Starts a subprocess in the [file_path] folder.
                - From that folder, calls the requested program.
                - Inputs the corresponding .DAT file for this data point into the program.
                - Copies the default .OUT file from the program to the Outputs 
                  folder with a more descriptive name.
            
            Then outputs a message to console to confirm that the batch script finished running. 

            Parameters
            ----------
            file_tag_batch : list of strings
                One tag for each data point in the batch, with format 
                e.g. "e0.001_g10.0_p0.730_Pb207".
                The list has length = number of data points in this batch.
                
            batch_index : int
                The index of this batch (numbered from 0 to num_batches-1).
                
            nucleus : string
                The nucleus being tested, with format e.g. "Pb207". 
                Used to navigate the file directory.
                
            file_path : string
                The relative file path from the current working directory to 
                the folder from which this batch should be run.

            Returns
            -------
            None.

            """
            
            
            batch_folder = "Batch"+str(batch_index)
            script_text = ""
            
            for file in file_tag_batch :
                script_text += ("\n(cd ../"+nucleus+"/Run/"+batch_folder+";" +  # start a subprocess and move to the Run folder for this batch.
                                " ./../../../Executables/MO/" + program +       # call the program from the Run folder.
                                " < ../../Inputs/"+ abr +"_"+ file +".DAT;" +   # input the relevant .DAT file to the program.
                                " cp "+program.upper()+ ".out"+                 # copy the default .OUT file... 
                                " ../../Outputs/"+ abr +"_" + file + ".OUT)")   # ...to a new .OUT file with a more descriptive name, in the Outputs folder.
            
            if verbose:
                script_text += ("\n\necho message from terminal: " +
                            "finished running "+ program +" batch "+ str(batch_index))
            
            script_file = open(file_path, 'w')
            script_file.write(script_text)
            script_file.close() 


        subprocesses = {}

        for b in range(num_batches):
            file_path = "../"+nucleus+"/Scripts/Run"+program.upper()+"_"+str(b+1)+".sh"
            batch_file_tags = file_tags[(b*num_per_batch):((b+1)*num_per_batch)]
            write_script_batch(batch_file_tags, b+1, nucleus, file_path)
            subprocesses[(program+"_"+str(b+1))] = subprocess.Popen(["sh", file_path])     
            # asynchronous call to start the program as a subprocess
            
        for b in range(num_batches):
            # wait to ensure it has finished (before starting to read outputs!), 
            # if it takes longer than the time limit seconds, throw an error to catch hangs.
            subprocesses[(program+"_"+str(b+1))].wait(allowed_time)             
            
    return run_script_batches





#%%

''' FUNCTIONS FOR READING GAMPN.OUT '''

def get_efac(lines):
     """ 
     A function which reads the full contents of the GAMPN.OUT file, 
     and returns the value of EFAC.
     
     Parameters
     ----------
     lines : list of strings
         The full contents of the GAMPN.OUT file.
         Each element of the list is a line read from the file.
    
     Returns:
     -------
     efac : float
         The value of the conversion factor from energy units hw to MeV.
     
     """
     ref = "     KAPPA    MY     EPS   GAMMA    EPS4     EPS6     W0/W00   NMAX  COUPL     OMROT      EFAC      QFAC\n"
     efac_line = lines.index(ref)
     efac = float(lines[efac_line+1][85:95].strip())
     return efac
 
    
def get_sp_level(lines, sp_index, parity):
    """ 
    A function which reads the full contents of the GAMPN.OUT file, 
    and returns the half-line containing data about the single particle level requested.
    
    Only half the line is required because the data for all 80 calculated levels 
    is output in two columns of 40 lines each.
    
    Parameters
    ----------
    lines : list of strings
        The full contents of the GAMPN.OUT file.
        Each element of the list is a line read from the file.
    
    index : int
        The index of the single particle level
        
    parity : string
        Defines how the single particle levels are indexed. 
        If '0', all the levels are counted (upwards from 1), regardless of parity.
        If '-', only the negative parity levels are counted (and positive parity levels are ignored.)
        If '+', only the positive parity levels are counted (and negative parity levels are ignored.)
        
    
    Returns:
    -------
    half_line : string
        A string containing data about the single particle level.
    
    """
    
    ref = "   #   ENERGY +/-(#)    <Q20>    <Q22>     <R2>     <JZ>      #   ENERGY +/-(#)    <Q20>    <Q22>     <R2>     <JZ>\n"
    levels_header_line = lines.index(ref)
    
    if parity == '0':
        # calculate the line number of the single particle level in the GAMPN.OUT file (indexed from zero!)
        sp_line = sp_index+levels_header_line+1  
        
        if sp_index > 40:
            # get only the second half of the line
            sp_line -= 40
            whole_line = lines[sp_line]
            half_line = whole_line[60:-1].strip()                                 
            
        else:
            # get only the first half of the line
            whole_line = lines[sp_line]
            half_line = whole_line[0:60].strip()   
            
    elif parity == '-' or parity == '+':
        
        reduced_lines = lines[levels_header_line+2 : levels_header_line+42]
        search = parity + '(#' + str(sp_index) + ')'
        
        found = False
        for l in reduced_lines:
            if search in l:
                whole_line = l
                found = True
                break
        
        if found:
            position = whole_line.index(search)
        else: raise RuntimeError("The requested single particle orbital was not calculated: " + search)
        
        if position > 60:
            half_line = whole_line[60:-1].strip()   
        else:
            half_line = whole_line[0:60].strip()   
        
    else: raise ValueError("Unrecognised parity. Allowed inputs are only '+', '-', or '0'.")
        
    
    return half_line
 
    
def get_sp_index(line, hash_index):
    """ 
    A function which takes a line of data about a single particle level, 
    and returns its orbital index.
    
    Parameters
    ----------
    line : string
        A string containing data about the single particle level.
    
    hash_index : int
        The index in the string that locates the "#" character.
        (Used as a reference point in the line).
        
    Returns:
    -------
    index : int
        The orbital index of the single particle level. 
        (numbered separately for positive and negative parity orbitals).
    
    """
    index_string = line[hash_index+1 : hash_index+3]
    if index_string[1] == ")":                                                  
        # in case the index is only a single digit, 
        # ignore the ")" that will have been caught:
        index_string = index_string[0]
        
    index = int(index_string)
    
    return index

def get_info(line):
    '''

    Parameters
    ----------
    line : string
        The line of text read from the GAMPN.OUT file which contains data 
        about a single particle (i.e. Nillson) level, including its energy,
        its index, its parity, etc.

    Returns
    -------
    parity : string
        The parity of the single particle level. Either '+' or '-'.
        
    energy_hw : float
        The energy of the single particle level, in oscillator (hbar omega) units.
        
    index : int
        The index of the single particle level, numbered separately for 
        positive and negative parity orbitals.

    '''
    hash_index = line.index("#")
    parity = line[hash_index-2]                            
    energy_hw = float(line[hash_index-10 : hash_index-4])  
    # energy_mev = energy_hw*inputs["efac"] # can't do this conversion without passing efac, and not really necessary in all cases anyway
    index = get_sp_index(line, hash_index)
    
    return (parity, energy_hw, index)



#%%

''' FUNCTIONS FOR READING ASYRMO.OUT '''

def get_delta(lines):
     """ 
     A function which reads the full contents of the ASYRMO.OUT file, 
     and returns the value of DELTA if it is well-defined, otherwise returns NaN.
     
     The value is deemed ill-defined if the error "SORRY I FOUND NO SOLUTION" appears.
     
     Parameters
     ----------
     lines : list of strings
         The full contents of the ASYRMO.OUT file.
         Each element of the list is a line read from the file.
    
     Returns:
     -------
     delta : float
         The value of the pairing gap energy DELTA in MeV.
     
     """
     for l in lines[10:20]:
         
         if "SORRY I FOUND NO SOLUTION" in l:
             delta = np.nan
             break
         
         if "DELTA=" in l:
             delta_string = l[7:13].strip()
             if delta_string == '*****':
                 delta = np.nan
             else:
                 delta = float(delta_string)
             break
     
     return delta
 
#%%

''' FUNCTIONS FOR READING PROBAMO.OUT '''

def find_gaps(spin_b_energies, index_b , spin_t_energies, index_t, exp):

    # find the lowest 13/2 state
    # check whether there are any 9/2 states at lower energy
    # if not, NaN
    # otherwise, get the energy gap between this 13/2 level and the closest 9/2 level below
    # check whether this gap is smaller than the experimental gap 
    # if so, AND there is another 9/2 level further down, get the gap to that level and compare with experiment (if this is closer, keep it, otherwise return to the first value)
    # otherwise, stop
    # find the next 13/2 state
    # get the best energy gap to its 9/2 below 
    # if this gap is closer to experiment, keep 

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
    
    # get the spins of the inital and final states of the transition.
    spin_string = line[dash_index-4:dash_index].strip()
    final_spin_string = line[dash_index+11:dash_index+16].strip()               
    
    if not(spin_string == final_spin_string):
        return False

    spin_float = spin_string_to_float(spin_string)
    
    #  get the energies of the inital and final states of the transition.
    try:
        this_energy = float(line[:6].strip())
    except ValueError: # could not convert string to float: '0.0  1'
        this_energy = float(line[:5].strip())
    
    final_energy = float(line[dash_index+3:dash_index+10].strip())
                         
    # determine whether the initial and final states are the same.
    if not(this_energy == final_energy):
        return False
    
    # read data.
    mag_moment = float(line[-8:].strip())
    line_data = {'spin_string': spin_string, 'spin_float':spin_float, 
                 'energy':this_energy, 'mag_moment':mag_moment}
    
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
        
    else:
        file_data[(spin+"_energies")] = [line_data["energy"]]
        file_data[(spin+"_mag_moments")] = [line_data["mag_moment"]]
    
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
            file_data[(spin+"_mag_moments")] = [np.NaN]

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
            if verbose: 
                print("Could not find any states with first excited spin in file " + str(d))
            new_data["x1_energies"].append(np.NaN)
            new_data["x1_mag_moments"].append(np.NaN)
        
        try:
            new_data["x2_energies"].append(old_data[d]["x2_energies"])
            new_data["x2_mag_moments"].append(old_data[d]["x2_mag_moments"])
        except(KeyError):
            if verbose: 
                print("Could not find any states with second excited spin in file " + str(d))
            new_data["x2_energies"].append(np.NaN)
            new_data["x2_mag_moments"].append(np.NaN)
        
        try:
            new_data["x3_mag_moments"].append(old_data[d]["x3_mag_moments"])
            new_data["x3_energies"].append(old_data[d]["x3_energies"])
        except(KeyError):
            if verbose: 
                print("Could not find any states with third excited spin in file " + str(d))
            new_data["x3_energies"].append(np.NaN)
            new_data["x3_mag_moments"].append(np.NaN)
    
    
    max_val = int(ispin)+1

    for i in range(max_val): 
        
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
                gs_spin_energies = this_level_data[:,0]  # the energies of the lowest state of the expected gs spin
            
            
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



#%%


''' FUNCTIONS FOR ANALYSING RESULTS '''

def report_mean(prop, verbose):
    '''
    Calculate the mean, standard deviation, and standard error in the mean
    of the input property.
    
    If the property is a set of energy levels then the mean/stdv/sterr will be 
    calculated separately for each level.
    
    Any NaN values are ignored (i.e. the size of the data set is assumed to decrease)

    Parameters
    ----------
    prop : PropertyData objetc
        The property for which the mean/stdv/sterr will be calculated
        
    verbose : bool
        Whether to print lots of info (i.e. all energy levels if more than one, and stdv as well as sterr)
        or just the basics (i.e the lowest energy level, and only stdv)

    Returns
    -------
    None.

    '''
    
    data = prop.data
    mask = ~np.isnan(data)
    not_nan = np.sum(mask, axis=0)
    
    mean = np.mean(data, axis=0, where=mask)
    stdv = np.std(data, axis=0, where=mask)
    
    print("\n" + prop.axis_label + ":")
    if verbose:
        print("\t # \t  mean ± sterr (stdv)")
    
    if isinstance(stdv, np.ndarray):
        sterr = [stdv[s]/np.sqrt(not_nan[s]) for s in range(len(stdv))]
        
        if verbose:
            itr = len(mean)
            
        else:
            itr = 1
        
        for i in range(itr):
            if verbose:
                print("\t", i+1, "\t %.1f ± %.1f \t (%.1f) " % (mean[i], sterr[i], stdv[i]))
            else:
                print("\t %.1f ± %.1f " % (mean[i], sterr[i]))
            
        
    else:
        sterr = stdv/np.sqrt(not_nan)
        
        if verbose:
            print("\t %.1f ± %.1f \t (%.1f) " % (mean, sterr, stdv))
        else:
            print("\t %.1f ± %.1f " % (mean, sterr))
    

    
            