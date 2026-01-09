#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 21 11:54 2025

@author: katyrr

Functions for analysing the results of the full PTRM calculation.

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
    