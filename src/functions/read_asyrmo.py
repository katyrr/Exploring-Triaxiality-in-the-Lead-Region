#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 21 11:47 2025

@author: katyrr

Functions for reading ASYRMO.OUT and processing contents.

"""

import numpy as np

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