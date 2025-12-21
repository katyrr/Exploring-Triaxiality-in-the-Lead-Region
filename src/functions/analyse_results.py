#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 21 11:54 2025

@author: katyrr

Functions for reading ASYRMO.OUT and processing contents.

"""

import numpy as np


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
            