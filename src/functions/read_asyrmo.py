#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 21 11:47 2025

@author: katyrr

Functions for reading ASYRMO.OUT and processing contents.

"""

import numpy as np
import os

import functions.file_handling as fh

def read_asyrmo(file_tags, data_subfolder_path, output_data):

    for i in file_tags:

        output_file_path = os.path.join(data_subfolder_path, "Outputs", f"ASY_{i}.OUT")
        lines = fh.read_file(output_file_path)
        
        if not "PARTICLE-ROTOR  MODEL" in lines[0]: # then something has gone wrong
            raise RuntimeError("File " + i + " raised error in ASYRMO output: \n" + lines[0] )
            
        output_data["delta"].append(get_delta(lines)) # this also checks for the "SORRY I FOUND NO SOLUTIONS" error.
            

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