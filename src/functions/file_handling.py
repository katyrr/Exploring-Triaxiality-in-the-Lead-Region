#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 21 11:37 2025

@author: katyrr

Functions for locating/creating/reading/writing files and folders.

"""

import os
import shutil

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


def locate_data_subfolder(argv):
    '''
    Use given command line argument to locate the data subfolder to be used for calculations.
    Typically named for the nucleus being studied, e.g. Pt177, but can be given any name.
    
    Searches for this folder in the data directory. If not found, it is created.
    If no command line argument is given, raises an error.
    
    Parameters
    ----------
    argv : list of strings
        A list of the input command line arguments.

    Returns
    -------
    abs_path : string
        The absolute directory path to the data subfolder to be used for calculations.

    Errors
    ------
    ValueError if there is no command line input.

    '''
    if len(argv) <= 1:
        raise ValueError("missing argument: name of folder containing config file")

    folder_name = argv[1]

    rel_path = os.path.join("data", folder_name)
    abs_path = os.path.abspath(rel_path)

    if not os.path.isdir(abs_path):
        print(f"folder not found at: {abs_path}, creating new folder")
        os.mkdir(abs_path)
    
    return abs_path

def locate_config(folder_path):
    '''
    Search for config file in given directory. 
    If not found, it is created from a template.
    
    Parameters
    ----------
    folder_path : string
        The absolute directory path to the data subfolder to be used for calculations.

    Returns
    -------
    config_path : string
        The absolute directory path to the config file.

    '''

    print(f"\nSearching for config file at: {folder_path}")

    folder_name = os.path.basename(folder_path)
    config_path = os.path.join(folder_path, f"config_{folder_name}.txt")

    if not os.path.isfile(config_path):
        template_rel_path = os.path.join("static", "config_template.txt")
        template_abs_path = os.path.abspath(template_rel_path)
        print(f"config file not found, generating new from template at: {template_abs_path}")
        shutil.copy(template_abs_path, config_path)

    return config_path

def setup_directory(folder, num_batches, OS):
    '''
    Attempts to give execute permissions to ptrm codes (may not work on Windows).
    Checks for the required directory structure, and creates any missing folders:
        
    Code/data/[folder]
                |-- inputs 
                |-- scripts
                |-- run
                    |-- batch1
                    |-- batch2 
                    |-- [etc. up to "batchN" for N=num_batches]
                |-- outputs
                |-- figures


    Parameters
    ----------
    folder : string
        The name of the folder that contains the config file.
        
    num_batches : int
        The number of batches to run the calculations in (i.e. the number of 
        computer processors to utilise, typically between 4-8 for most PCs).

    OS : either "MacOS" or "64bit"
        Which version of the pre-compiled PTRM Fortran codes to use, depending 
        on your computer's operating system (Mac, or Windows/Linux)

    Returns
    -------
    None.

    '''
    
    programs = ["gampn", "asyrmo", "probamo"]
    
    for i in programs:
            
        if OS == "64bit":
            i = i.upper() + ".exe"
        path_to_program = os.path.join("src", "ptrm", OS, "MO", i)

        permissions = oct(os.stat(path_to_program).st_mode)[-3:]
        if permissions != "775":
            print(f"{i} does not have execute permissions... attempting to turn on")
            os.chmod(path_to_program, 0o775)

            permissions = oct(os.stat(path_to_program).st_mode)[-3:]
            if permissions != "775":
                raise RuntimeError("Failed to turn on execute permissions: please do manually")
            else:
                print(f"\tSuccessfully turned on execute permissions for {i}")
                

              
    required_folders = ["inputs", "scripts", "run", "outputs", "figures"]
    for i in range(1, num_batches+1):
        batch_folder = os.path.join("run", f"batch{i}")
        required_folders.append(batch_folder)
        
    for i in required_folders:
        path_to_i = os.path.join(folder, i)
        abs_path_to_i = os.path.abspath(path_to_i)
        if not os.path.isdir(abs_path_to_i):
            os.mkdir(abs_path_to_i)
            print(f"created directory: {abs_path_to_i}") 

    return 

