#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 21 11:54 2025

@author: katyrr

Functions for analysing the results of the full PTRM calculation.

"""

import numpy as np

import src.functions.structs as st
import src.functions.graph_plotting as gr


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
    
    
    if isinstance(stdv, np.ndarray):
        sterr = [stdv[s]/np.sqrt(not_nan[s]) for s in range(len(stdv))]
    else:
        sterr = [stdv/np.sqrt(not_nan)]
        mean = [mean]
        stdv = [stdv]
        
    if verbose:
        itr = len(mean)
        print(f"\n{prop.axis_label}:")
        print("\t# \tmean \t± \tsterr \t(stdv)")
    else:
        itr = 1
    
    for i in range(itr): 
        if verbose:
            #print("\t", i+1, "\t %.1f ± %.1f \t (%.1f) " % (mean[i], sterr[i], stdv[i]))
            print(f"\t{i+1}\t{mean[i]:<5.1f}\t± \t{sterr[i]:<5.1f}\t({stdv[i]:^5.1f})")
        else:
            #print("\t %.1f ± %.1f " % (mean[i], sterr[i]))
            print(f"{prop.axis_label:<60}{mean[i]:<5.1f} ± {sterr[i]:<5.1f}")


def check_agreement(verbose, data_points, num_comparisons):
    """
    A function to check how well data points agreed with experimental values,
    and which data point(s) had the highest agreement.
    

    """
    
    sorted_indices = np.argsort(data_points["agreed"])
    sorted_eps = [data_points["eps"][i] for i in sorted_indices]
    sorted_gamma = [data_points["gamma_degrees"][i] for i in sorted_indices]
    
    max_agreement = data_points["agreed"][sorted_indices[-1]]
    print(f"Highest agreement = {max_agreement} / {num_comparisons}")
    if max_agreement == 0:
        return
    
    unique_values, counts = np.unique(data_points["agreed"], return_counts=True)
    
    print("\n\n***** Agreement of each data point with experimental data: *****")
    if verbose:
        print(data_points["agreed"])
        
        print(dict(zip(unique_values, counts)))
        
        print("Number of data points with each level of agreement:")
        for i in range(len(unique_values)): # e.g. 0, 1
        
            print(f"\n\tAgreement = {unique_values[i]}:")
            
            # e.g. logic: 

            # when i = 0
            # we want the first 2 values of deformation, (because counts[i] = counts[0] = 2)
            # because they all have agreement = 0, (because unique_values[i] = unique_values[0] = 0)
            # so we need slice range [0:2], (because i=0, and counts[i]=2) 
            
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
                print(f"\t\t(ε, γ) = ({sorted_eps[j]:.3f}, \t{sorted_gamma[j]:.1f}º)")
        
    
    else:
        print(f"\nPoints with agreement = {max_agreement}:")
        
        i = len(unique_values)-1
        
        print(f"\n\tAgreement = {unique_values[i]}:")
    
        lower = sum(counts[0:i])
        upper = sum(counts[0:i+1])
        
        these_eps = sorted_eps[lower:upper]
        these_gamma = sorted_gamma[lower:upper]
        
        r = min(len(these_eps), 5)
        for j in range(r):
            print(f"\t\t(ε, γ) = ({these_eps[j]:.3f},\t{these_gamma[j]:.1f}º)")
        
        if r==5 and len(these_eps) > 5:
            print(f"... etc, {len(these_eps)}")
      
def plot_agreement(data_points, num_comparisons, code_settings, ptrm_inputs, gs_spin_floats, subtitle, data_subfolder_path):
    agreement = st.PropertyData(data_points["agreed"], "Agreement of Data Points With Experimental Data")
    agreement.contour_levels = np.arange(0, num_comparisons+2, dtype=int) #fn.calc_contour_levels(agreement.data)
    agreement.cbar_ticks = gr.calc_cbar_ticks(agreement.contour_levels)
    agreement.cbar_tick_labels = list(np.arange(0, num_comparisons+1, dtype=int)) #fn.calc_cbar_tick_labels(agreement.data, "int")
    agreement.experimental_data = np.NaN
    agreement.error_tolerance = np.NaN

    agreement.plot = 0

    if ptrm_inputs["deformation_input"] == "mesh" and agreement.plot: 
        
        code_settings["current_graph"] = agreement.title
        print("plotting graph: %(current_graph)s" % code_settings) 
        gr.plot_mesh_graph(agreement, data_points, code_settings, ptrm_inputs, gs_spin_floats, subtitle, data_subfolder_path)
     