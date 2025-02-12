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

- FUNCTIONS FOR PLOTTING GRAPHS

"""

import numpy as np                              # for np.arrays
import subprocess                               # for calling shell scripts to run 
import matplotlib.pyplot as plt                 # for plotting graphs
from matplotlib.ticker import FuncFormatter     # for formatting axis ticks
import matplotlib.tri as tri                    # for manual triangulation before drawing a contour plot



#%%

''' GENERAL FUNCTIONS 

- spin_float = spin_string_to_float(spin_string)

- spin_string = spin_float_to_string(spin_float)
    
- lines = read_file(path)

- void = write_file(path, text)
    
'''

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


''' FUNCTIONS FOR READING CONFIG 

- split_string = remove_inline_comments(split_string)

- void = check_line_format(split_string, line, l)
               
- range_list, step = range_to_list(range_string)

- step = get_range_step(range_string):

- eps_points, gamma_points = arrange_data_line(split_string)

- eps_points, gamma_points = arrange_mesh(split_string)

'''

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
        raise ValueError("Line "+str(l)+ " is too long: " + line)                                
        
    elif len(split_string)<expected_num_words :
        raise ValueError("Line "+str(l)+ " is too short:" + line)
            
    
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

''' FUNCTIONS FOR RUNNING CODES 

- orbitals_string = write_orbitals(fermi_level, number, parity):

- run_script_batches(program) = configure_script_writer(file_tags, nucleus, num_batches, num_per_batch, allowed_time)
    void = run_script_batches(program) # closure                              
        void = write_script_batch(file_tag_batch, batch_index, nucleus, file_path) # sub-function    

'''

def write_orbitals(fermi_level, number, parity):
    """
    A function to generate an parity + orbital string for input into gampn or asyrmo.
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

def est_e2plus(eps, A):
    """
    A function to dynamically calculate the effective E2PLUS value of a single 
    data point, using Grodzin's relation. 
    
    eps is assumed to be interchangable with beta, to the precision of this estimate,
    as reccomended on p6 of the manual, from which the formula is taken:
        
    E2PLUS approx = 1225 / [pow(beta, 2) * pow(A, 7/3)] MeV
    
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

''' FUNCTIONS FOR READING GAMPN.OUT

- efac = get_efac(lines)

- half_line = get_fermi_level(lines, fermi_level_index)
    
- fermi_index = get_fermi_index(fermi_level_line, hash_index)

'''

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
 
    
def get_fermi_level(lines, fermi_level_index):
    """ 
    A function which reads the full contents of the GAMPN.OUT file, 
    and returns the half-line containing data about the fermi level.
    
    Only half the line is required because the data for all 80 calculated levels 
    is output in two columns of 40 lines each.
    
    Parameters
    ----------
    lines : list of strings
        The full contents of the GAMPN.OUT file.
        Each element of the list is a line read from the file.
    
    fermi_level_index : int
        The index of the highest filled orbital, numbered from 1 (all parities together).
        
    Returns:
    -------
    half_line : string
        A string containing data about the fermi level.
    
    """
    
    ref = "   #   ENERGY +/-(#)    <Q20>    <Q22>     <R2>     <JZ>      #   ENERGY +/-(#)    <Q20>    <Q22>     <R2>     <JZ>\n"
    levels_header_line = lines.index(ref)
    fermi_level_line = fermi_level_index+levels_header_line+1                  
    # calculate the line number of the fermi level in the GAMPN.OUT file (indexed from zero!)
    
    if fermi_level_index > 40:
        fermi_level_line -= 40
        whole_line = lines[fermi_level_line]
        half_line = whole_line[60:-1].strip()                                 
        # get only the second half of the line
    else:
        whole_line = lines[fermi_level_line]
        half_line = whole_line[0:60].strip()                                   
        # get only the first half of the line
    
    return half_line
 
    
def get_fermi_index(fermi_level_line, hash_index):
    """ 
    A function which takes a line of data about the fermi level, 
    and returns its orbital index.
    
    Parameters
    ----------
    fermi_level_line : string
        A string containing data about the fermi level.
    
    hash_index : int
        The index in the string that locates the "#" character.
        (Used as a reference point in the line).
        
    Returns:
    -------
    fermi_index : int
        The orbital index of the Fermi level 
        (numbered separately for positive and negative parity orbitals).
    
    """
    fermi_index_string = fermi_level_line[hash_index+1 : hash_index+3]
    if fermi_index_string[1] == ")":                                                  
        # in case the index is only a single digit, 
        # ignore the ")" that will have been caught.
        fermi_index_string = fermi_index_string[0]
        
    fermi_index = int(fermi_index_string)
    
    return fermi_index



#%%

''' FUNCTIONS FOR READING PROBAMO.OUT 

- line_data = read_data(line)

- file_data = sort_by_spin(line_data, file_data)

- file_data = sort_by_expectation(line_data, file_data, inputs)

- file_data = missing_data(file_data, max_spin)

- new_data = restructure_data(old_data, max_spin)

- multi_level_data = fill_gaps(multi_level_data)

'''

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
    multi_level_data : a matrix (2D np.array) of floats
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
    


#%%

''' FUNCTIONS FOR PLOTTING GRAPHS 

- contour_levels = calc_contour_levels(data)

- cbar_tick_labels = calc_cbar_tick_labels(data)

- cbar_ticks = calc_cbar_ticks(contours)

- exp, tol = try_experimental(inputs, key, tolerance)
    
- void = format_fig(polar_or_linear, ax, legend_handles, title, **kwargs)

- cax, cbar = draw_contour_plot(ax, prop, data_points)

- correct_spin_range = find_correct_spin(gs_spins, experimental_value)

- correct_spin_handle = plot_correct_spin(correct_spin_range, var, step, prop)

- legend_handles = plot_points_with_experiment(data_points, prop, legend_handles, cbar)

- legend_handles = plot_points_without_experiment(data_points, legend_handles)

- var_sym, var, fix_sym, fix = assign_parameters(inputs, data_points)
    
- legend_handles, legend_title = plot_multi_lines(prop, var, legend_handles, marker_size, fix_sym, fix_val)

- legend_handles = mark_spin(inputs, data_points, output_data, legend_handles, ax)

'''

def calc_contour_levels(data):
    """
    A function to calculate a list of custom contour levels for a dataset. 

    Parameters
    ----------
    data : np.array of floats
        The data set to which the contour levels will apply.

    Returns
    -------
    contour_levels : np.array of floats
        The contour levels for a contour plot of this data.

    """
    
    min_contour = min(data)-0.5
    max_contour = max(data)+1.5
    
    contour_levels = np.arange(min_contour, max_contour, 1.0)
    return contour_levels


def calc_cbar_tick_labels(data):
    """
    A function to calculate a list of custom colour bar tick labels for a dataset. 

    Parameters
    ----------
    data : np.array of floats
        The data set to which these colour bar tick labels will apply.

    Returns
    -------
    cbar_tick_labels : list of strings
        The colour bar tick labels for this data.

    """
    cbar_ticks = np.arange(min(data), max(data)+1.0, 1.0)
    cbar_tick_labels = [spin_float_to_string(n) for n in cbar_ticks]
    return cbar_tick_labels


def calc_cbar_ticks(contours):
    """
    A function to calculate a list of custom colour bar ticks for a dataset. 

    Parameters
    ----------
    contour_levels : np.array of floats
        The contour levels for a contour plot of this data.

    Returns
    -------
    cbar_ticks : list of floats
        The colour bar ticks for this data.

    """
    num = len(contours) - 1
    
    # Calculate midpoints of levels for tick placement
    cbar_ticks = [(contours[i] + contours[i+1]) / 2 for i in range(num)]   
     
    return cbar_ticks


def try_experimental(inputs, key, tolerance):
    """
    A function which gets the experimental value of a nuclear property, and 
    returns it with the absolute error tolerance.
    
    If the experimental data has not been input, returns both the experimental 
    value and the tolerance as np.NaN.

    Parameters
    ----------
    inputs : dictionary
        Input settings from a config file. 
        May contain some experimental data.
        
    key : string
        The key for the requested data in the "inputs" dictionary.
       
    tolerance : float or int
        The absolute error tolerance when comparing calculated values to the 
        experimental value.

    Returns
    -------
    exp
        The experimental value.
    TYPE
        The absolute error tolerance.

    """
    try:
        exp = inputs[key]
        tol = tolerance
        
    except KeyError:
        exp = np.NaN
        tol = np.NaN
        
    return exp, tol


def format_fig(polar_or_linear, ax, legend_handles, title, **kwargs):
    """
    A function that handles formatting of a graph:
        - axis ranges
        - axis ticks and tick labels
        - axis labels
        - legend
        - graph title

    Parameters
    ----------
    polar_or_linear : string
        "polar" for a polar plot.
        "linear" for a linear plot.
        
    ax : Axes object
        The axes on which the graph is being plotted.
        
    legend_handles : list of object handles.
        A list of all the objects plotted on the graph which should be included in the legend.
        
    title : string
        The title for the graph.
        
    **kwargs : 
        
        - "varied" : list or npp.array of floats
            The independent variable.
            Only applies to linear graphs.
            Either eps, gamma, or e2plus.
            Default value = None.
            
        - "x_label" : string
            The label for the x-axis variable.
            Only applies to linear graphs.
            Default value = None.
        
        - "y_label" : string
            The label for the y-axis variable. 
            Only applies to linear graphs.
            Default value = None.
            
        - "legend_title" : string
            The title for the legend.
            Only applies to linear graphs.
            Default value = None.

    Raises
    ------
    ValueError
        Occurs in "polar_or_linear" is not one of "polar" or "linear".

    Returns
    -------
    None.

    """
    
    if polar_or_linear == 'polar':
        ax.set_thetamin(0)   
        ax.set_thetamax(60)  
        
        theta_ticks = np.arange(0, 70, 10)  
        ax.set_xticks(np.radians(theta_ticks))
        # set the number of decimal places 
        ax.yaxis.set_major_formatter(FuncFormatter(lambda x, _: f'{x:.2f}'))    
        
        plt.xlabel("ε")
        # gamma axis label
        ax.text(45*np.pi/180, ax.get_rmax()*1.2, "γ", ha='center', va='center') 
        
        
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
    
    else: 
        raise ValueError("unrecognised graph type: " + polar_or_linear + 
                         "; must be either 'polar' or 'linear'.")


def draw_contour_plot(ax, prop, data_points):
    """
    A function to draw a filled contour plot in polar coordinates.
    
    If the property is energy levels grouped by spin (e.g. spin_1/2_energies) 
    then the data may be 2D (i.e. contains multiple values for each data point).
    This cannot be represented on this kind of graph, so plot only the first 
    value (i.e. the yrast state of that spin).
    
    The data points are triangulated over the deformation space 
    (i.e. with respect to gamma and eps; no relation to the data values.)
    
    The triangulation is masked for any np.NaN values (i.e. missing data points
    are not included in the plot, and will appear in white.)
    
    The filled contours are plotted, with a colour bar.
    
    If custom colour bar ticks and labels have been calculated, then set those values.
    

    Parameters
    ----------
    ax : Axes object
        The axes on which the graph is being plotted.
        
    prop : PropertyData object
        The nuclear property that is being plotted, collected into a class with 
        information about its graph plotting features.
        
    data_points : dictionary
        Contains lists of data point values, including "gamma_radians" and "eps".
        All the lists in this dictionary have the same length = total number of data points.

    Returns
    -------
    cax : TriContourSet object
        The set of contour lines / regions for the plot.
        
    cbar : Colorbar object
        The color bar of the filled contour plot.

    """
    
    if prop.sort == "Spin ":
        
        plot_data = np.transpose(prop.data)[0]
        
        print('''multiple levels can't be plotted on a contour plot; 
              plotting only the yrast state of spin ''' + prop.num)
        
        
    else:
        plot_data = np.array(prop.data)
        
    # some of the data points may be NaN, 
    # so manually create the triangulation and mask NaN triangles.
    triang = tri.Triangulation(data_points["gamma_radians"], data_points["eps"])
    mask = np.any(np.isnan(plot_data[triang.triangles]), axis=1)
    triang.set_mask(mask)
    
    cax = ax.tricontourf(triang, plot_data, levels=prop.contour_levels)
    cbar = plt.colorbar(cax, pad=0.1, label=prop.axis_label)


    
    if prop.cbar_tick_labels:        # then format for discrete values
        
        cbar.set_ticks(prop.cbar_ticks)
        cbar.set_ticklabels(prop.cbar_tick_labels)
    
        
    return cax, cbar


def find_correct_spin(gs_spins, experimental_value):
    """
    A function to locate the 1D region(s) of a line graph in which the ground 
    state spin has been correctly reproduced.
    
    The boundaries of the region(s) are assumed to fall at data points where:
        - either the previous point was incorrect and this point is correct,
        - or the previous point was correct and this point is incorrect.

    Parameters
    ----------
    gs_spins : np.array of floats
        Contains the ground state spin calculated at each data point.
        
    experimental_value : float
        The experimental ground state spin.

    Returns
    -------
    correct_spin_range : list of ints
        The indices of the data points that lie at a boundary between correct 
        and incorrect ground state spins.

    """
    correct_spin_range = []                                             
    start_flag = False
    
    for i in range(len(gs_spins)):
        if gs_spins[i] == experimental_value and not start_flag:  
            # this point is correct and the previous point was incorrect
            start_flag = True
            correct_spin_range.append(i)
            
        elif gs_spins[i] != experimental_value and start_flag:    
            # this point is incorrect, and previous point was correct
            start_flag = False
            correct_spin_range.append(i)
    
    return correct_spin_range


def plot_correct_spin(correct_spin_range, var, step, prop):
    """
    A function to plot the region(s) of correct ground state spin onto a line graph,
    as a green box.

    Parameters
    ----------
    correct_spin_range : list of ints
        The indices of the data points that lie at a boundary between correct 
        and incorrect ground state spins.
        
    var : np.array of floats
        The independent variable.
        
    step : float
        The step size of the independent variable.
        
    prop : PropertyData object
        The nuclear property that is being plotted, collected into a class with 
        information about its graph plotting features.

    Returns
    -------
    correct_spin_handle : Line2D object
        The handle of the line that marks the correct spin boundary.

    """
    
    for r in range(len(correct_spin_range)):
        
        if prop.sort == "Spin ":
            data = prop.data[0] 
            for i in range(1, len(prop.data)):
                data += prop.data[i]
        else:
            data = prop.data
            
        # front edge
        if correct_spin_range[r] == 0:                                  
            # the first value of eps in the range has the correct spin
            start_range = (var[correct_spin_range[r]]- step/2)          
        else:
            start_range = np.mean([var[correct_spin_range[r]-1], 
                                   var[correct_spin_range[r]]])
        
        correct_spin, = plt.plot([start_range, start_range], 
                                [min(data)-0.05*max(data), max(data)*1.05], 
                                'g-', label="range of correct spin")   
        # end edge
        if r%2==0:
            if r+1 == len(correct_spin_range):                          
                # the last value of eps in the range has the correct spin
                end_range = (var[-1]+step/2)
            else:
                end_range = np.mean([var[correct_spin_range[r+1]-1], 
                                     var[correct_spin_range[r+1]]])
       
        correct_spin_handle, = plt.plot([end_range, end_range], 
                                [min(data)-0.05*max(data), max(data)*1.05], 
                                'g-', label="range of correct spin")  
    
        # top and bottom edges
        plt.plot([start_range, end_range],                          
                 [min(data)-0.05*max(data), 
                  min(data)-0.05*max(data)], 'g-')
        plt.plot([start_range, end_range], 
                 [max(data)*1.05, 
                  max(data)*1.05], 'g-')                  
        
    return correct_spin_handle


def plot_points_with_experiment(data_points, prop, legend_handles, cbar):
    """
    A function to plot data points on a graph in polar coordinates. 
    Additionally compares the value of each data point to an experimental value,
    and marks it in red if it agrees within tolerance.
    
    If the data set is large, the data point markers are smaller, and ONLY points 
    that match experiment are plotted, to avoid overly cluttering the graph.
    
    The experimental value is also marked with a red line on the colour bar.

    Parameters
    ----------
    data_points : dictionary
        Contains lists of data point values, including "gamma_radians" and "eps".
        All the lists in this dictionary have the same length = total number of data points.
        
    prop : PropertyData object
        The nuclear property that is being plotted, collected into a class with 
        information about its graph plotting features.
        
    legend_handles : list of object handles.
        A list of all the objects so-far plotted on the graph which should be 
        included in the legend.
        
    cbar : Colorbar object
        The color bar of the filled contour plot.
        
    Returns
    -------
    legend_handles : list of object handles.
        A list of all the objects so-far plotted on the graph which should be 
        included in the legend, now including the objects plotted in this function.
        
    
    """
    
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
            # update the record of how many of the tested properties agree
            data_points["agreed"][r] += 1                                       
            hit = plt.scatter(data_points["gamma_radians"][r], 
                      data_points["eps"][r], s=marker_size, edgecolor='red', 
                      facecolor='None', label="matches experiment")
            legend_hit = True
            
        # for data points that don't match experimental data, only plot them 
        # when the data set is quite small, to avoid cluttering the graph.
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
    """
    A function to plot data points on a graph in polar coordinates. 
    No comparison to any experimental data is made.
    
    If the data set is large, the data point markers are smaller, to avoid cluttering the graph.

    Parameters
    ----------
    data_points : dictionary
        Contains lists of data point values, including "gamma_radians" and "eps".
        All the lists in this dictionary have the same length = total number of data points.

    legend_handles : list of object handles.
        A list of all the objects so-far plotted on the graph which should be 
        included in the legend.
        
    Returns
    -------
    legend_handles : list of object handles.
        A list of all the objects so-far plotted on the graph which should be 
        included in the legend, now including the data points plotted in this function.
        
    

    """
    
    # use a smaller marker size for large data sets
    if len(data_points["file_tags"]) < 100: marker_size = 5
    else: marker_size = 1
    
    all_points = plt.scatter(data_points["gamma_radians"], 
              data_points["eps"], s=marker_size, c='w', label="data point")
    legend_handles.append(all_points) 
    
    return legend_handles


def assign_parameters(inputs, data_points):
    """
    For a linear plot, determine which is the independent variable, and which 
    are held constant.

    Parameters
    ----------
    inputs : dictionary
        Input settings from a config file. 
        May contain some experimental data.
        
    data_points : dictionary
        Contains lists of data point values, including "gamma_radians" and "eps".
        All the lists in this dictionary have the same length = total number of data points.


    Raises
    ------
    ValueError
        Occurs if none of eps, gamma, or e2plus are varied.

    Returns
    -------
    var_sym : string
        An axis label for the independent variable.
        "ε" for eps,
        "γ / º" for gamma,
        "E2PLUS / MeV" for e2plus.
        
        
    var : np.array of floats
        The independent variable values at each data point.
        
    fix_sym : string
        A symbol to represet the fixed variable(s).
        "ε" for eps (and e2plus implied),
        "γ" for gamma (and e2plus implied),
        "(ε, γ)" for eps and gamma.
        
    fix : no.array of floats
        The value of the fixed deformation variable.
        Has length 1.
        Empty if both eps and gamma are fixed.
       

    """

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
        
    elif len(data_points["e2plus"]) > 1:
        
        var_sym = "E2PLUS / MeV"
        var = data_points["e2plus"]
        fix_sym = "(ε, γ)"
        fix = []
        
    else: raise ValueError("unrecognised graph request")

    return (var_sym, var, fix_sym, fix)


def plot_multi_lines(prop, var, legend_handles, marker_size, fix_sym, fix_val):
    """
    A function to plot multiple lines for one property on a line graph.
    e.g. for spin_1/2_energies, which may contain multiple energy levels at 
    that spin for each data point.
    
    

    Parameters
    ----------
    prop : PropertyData object
        The nuclear property that is being plotted, collected into a class with 
        information about its graph plotting features.
        
    var : np.array of floats
        The independent variable.
        
    legend_handles : list of object handles.
        A list of all the objects so-far plotted on the graph which should be 
        included in the legend.
        
    marker_size : int
        The size of the data point markers.
        
    fix_sym : string
        A symbol representing the fixed deformation parameter.
        
    fix_val : np.array
        The value of the fixed deformation parameter.
        Has length 1.

    Returns
    -------
    legend_handles : list of object handles.
        A list of all the objects so-far plotted on the graph which should be 
        included in the legend, now including the lines plotted in this function.
    
    legend_title : string
        A title for the legend.

    """
    data_by_line = np.transpose(prop.data)
    line_colours = ['k-x', 'b-x', 'y-x']
    line_labels = ["lowest energy", "second lowest energy", "third lowest energy"]
    
    for s in range(min(len(line_labels), np.size(data_by_line,0))):
        
        data, = plt.plot(var, data_by_line[s], line_colours[s], label=line_labels[s], markersize=marker_size)
        legend_handles.append(data)
        
    legend_title = "%s = %s" % (fix_sym, fix_val)
    
    return legend_handles, legend_title


def mark_spin(inputs, data_points, spin_data, legend_handles, ax):
    """
    A function to plot the region(s) of correct ground state spin onto a polar plot,
    as a black contour line.

    Parameters
    ----------
    inputs : dictionary
        Input settings from a config file. 
        May contain some experimental data.
        
    data_points : dictionary
        Contains lists of data point values, including "gamma_radians" and "eps".
        All the lists in this dictionary have the same length = total number of data points.

    spin_data : np.array of floats
        The ground state spins calculated at each data point.
    
    legend_handles : list of object handles.
        A list of all the objects so-far plotted on the graph which should be 
        included in the legend.
    
    ax : Axes object
        The axes on which the graph is being plotted.

    Returns
    -------
    legend_handles : list of object handles.
        A list of all the objects so-far plotted on the graph which should be 
        included in the legend, now including the spin region plotted in this function.
    
    """

    correct_range = [inputs["gs_spin_float"]-0.5, 
                     inputs["gs_spin_float"]+0.5]
    spin_colour = (0,0,0) #(213/255,1,0)
    
    ax.tricontour(data_points["gamma_radians"], data_points["eps"], 
                  spin_data, levels=correct_range,  
                  colors=[spin_colour], linewidths=1.0)
    spin_legend_proxy = plt.Line2D([], [], color=spin_colour, linewidth=1.0, label="region of correct g.s. spin") 
    legend_handles.append(spin_legend_proxy)

    return legend_handles


def check_agreement(verbose, data_points, num_comparisons):
    """
    A function to check how well data points agreed with experimental values,
    and which data point(s) had the highest agreement.
    
    Parameters
    ----------
    

    Returns
    -------
    None.

    """
    
    sorted_indices = np.argsort(data_points["agreed"])
    sorted_eps = [data_points["eps"][i] for i in sorted_indices]
    sorted_gamma = [data_points["gamma_degrees"][i] for i in sorted_indices]
    
    max_agreement = data_points["agreed"][sorted_indices[-1]]
    
    unique_values, counts = np.unique(data_points["agreed"], return_counts=True)
    
    print("\n\n***** Agreement of each data point with experimental data: *****")
    if verbose:
        print(data_points["agreed"])
        
        print(dict(zip(unique_values, counts)))
        
        print("Number of data points with each level of agreement:")
        for i in range(len(unique_values)): # 0, 1
        
            print("\n\tAgreement = " + str(unique_values[i]) + ":")
            
            # when i = 0
            # we want the first 2 values of deformation, (2 = counts[0] = counts[i])
            # because they all have agreement = 0, (0 = unique_values[0] = unique_values[i])
            # so we need slice range [0:2], (0 = i, 2 = counts[i]) 
            
            # when i = 1
            # we want the next 8 values of deformation, (8 = counts[1] = counts[i])
            # because they all have agreement = 3, (3 = unique_values[1] = unique_values[i])
            # and account for existing values (2), (2 = counts[0] = counts[i-1] = sum(counts[0:1]) = sum(counts[0:i]))
            # so we need slice range [2, 10], (2 = sum(counts[0:i]), 10 = sum(counts[0:i+1])
        
            lower = sum(counts[0:i])
            upper = sum(counts[0:i+1])
            
            these_eps = sorted_eps[lower:upper]
            these_gamma = sorted_gamma[lower:upper]
            
            for j in range(len(these_eps)):
                print("\t\t(ε, γ) = (" + str(these_eps[j]) + ",\t" + str(these_gamma[j])+"º)")
        
       
    
    else:
        
        print("Highest agreement = " + str(max_agreement) + " / " + str(num_comparisons))
                
        print("\nPoints with agreement = " + str(max_agreement) + ":")
        
        i = len(unique_values)-1
        
        print("\n\tAgreement = " + str(unique_values[i]) + ":")
    
        lower = sum(counts[0:i])
        upper = sum(counts[0:i+1])
        
        these_eps = sorted_eps[lower:upper]
        these_gamma = sorted_gamma[lower:upper]
        
        for j in range(len(these_eps)):
            print("\t\t(ε, γ) = (" + str(these_eps[j]) + ",\t" + str(these_gamma[j])+"º)")
    
        
       
