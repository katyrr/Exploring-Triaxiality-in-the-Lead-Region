#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  8 23:07:07 2025

@author: katyrr

"""


import matplotlib.pyplot as plt                 # for plotting graphs
from matplotlib.ticker import FuncFormatter     # for formatting axis ticks
import matplotlib.tri as tri                    # for manual triangulation before drawing a contour plot
import matplotlib.colors as colors
import numpy as np                              # for np.arrays
import os
import datetime

from src.functions.spin_processing import spin_string_to_float, spin_float_to_string
import src.functions.structs as st 

def plot_line_graph(prop, ptrm_inputs, data_points, code_settings, experimental_data, subtitle, gs_spin_floats, fig_path):
    '''
    - Draw a green box around regions that have the correct ground state spin, if requested.
    - Plot a red line to indicate the experimental value, if available.
    '''

    # set which paramters are varied and which are constant
    var_sym, var, fix_sym, fix = assign_parameters(ptrm_inputs, data_points)
    
    _, ax = plt.subplots() 
    
    legend_handles = []
    legend_handles, legend_title = plot_line_data(data_points, prop, var, fix_sym, fix, legend_handles)

    
    # if experimental data is available, plot it in red for easy comparison
    if np.isfinite(prop.experimental_data).all() and not prop.num == "all": 
        legend_handles = plot_exp_line(prop, code_settings, var, legend_handles)

        
    # mark the range in which the correct ground state spin was calculated
    if code_settings["mark_spin"]==1:
        
        correct_spin_range = find_correct_spin(gs_spin_floats.data, experimental_data["gs_spin_float"])
        if len(correct_spin_range) > 0:
            spin = plot_correct_spin(correct_spin_range, var, ptrm_inputs["step"], prop)
            legend_handles.append(spin)
                
    format_fig('linear', ax, list(reversed(legend_handles)), 
                '%(current_graph)s in %(nucleus)s' % ptrm_inputs, subtitle, 
                varied=var, x_label=var_sym, y_label=prop.axis_label, 
                legend_title=legend_title)
    
    if prop.prop == "delta":
        ax.set_ylim([0.2,1]) 
        
    if prop.cbar_tick_labels:        # then format for discrete values
        ax.set_yticks(prop.cbar_ticks)
        ax.set_yticklabels(prop.cbar_tick_labels)
    
    
    plt.savefig(fig_path, bbox_inches = 'tight')

    if code_settings["display_figures"]:
        plt.show()


def plot_mesh_graph(prop, data_points, code_settings, ptrm_inputs, gs_spin_floats, subtitle, fig_path):
    """
    - Draw a contour line to indicate the perimeter of the region where 
        the ground state spin was correctly reproduced, if requested.
        - Plot data point markers.
            - If experimental data is available, points that agree with experiment 
            (within tolerance) are marked in red.
            - Non-matching points are not plotted (unless there are fewer than 100 data points.)
    """

    _, ax = plt.subplots(subplot_kw=dict(projection='polar'))
    _, cbar = draw_contour_plot(ax, prop, data_points)
    
    legend_handles = []
    
    if code_settings["mark_spin"]:
        
        legend_handles = mark_spin(ptrm_inputs, data_points, gs_spin_floats.data, legend_handles, ax)
        
    # plot the data point markers, with comparison to experiment if possible
        
    legend_handles = plot_points(data_points, prop, legend_handles, cbar, code_settings)
    
    format_fig('polar', ax, legend_handles, '%(current_graph)s of %(nucleus)s' % ptrm_inputs, subtitle)

    plt.savefig(fig_path, bbox_inches = 'tight')


    if code_settings["display_figures"]:
        plt.show()


def prepare_data_to_plot(experimental_data, file_tags, restructured_output_data):
    ''' 6. PREPARE TO PLOT GRAPHS 

    - Record each data set in an instance of class PropertyData.
    - Calculate graph plotting attributes and store within the class.

    - Raise a ValueError if the property isn't recognised 
    (i.e. if more data sets are recorded in the future, they cannot be plotted without
    first hard-coding the calculation of things like axis labels, contour levels,
    colour bar ticks, etc).
    
    - Create a new data set containing all energies (of all spins) to plot together.
    - Create a new data set with all energies shifted to be relative to the expected 
    ground state (not necessarily the same as the calculated ground state at all points).
    This makes the output lines look smoother (no sharp bends when the ground state changes).
    - Create a new data set containing root mean squared error (i.e. discrepancy) between 
    the calculated lowest energy states of each spin and the exeperimental values (where available).

    '''

    # convert restructured_output_data from a dictionary of lists to a dictionary of PropertyData objects 
    data_to_plot = {}
    for i in restructured_output_data:
        # print(i)
        data_to_plot[i] = st.PropertyData(restructured_output_data[i], i)
        
        # calculate contour levels, colour bar ticks and labels, 
        # and assign experimental values and error tolerance if available.
        
        data_to_plot[i] = calculate_format_data(data_to_plot[i], i, experimental_data)
        

    data_to_plot["all_energies"] = collate_energy_data(data_to_plot, len(file_tags), 
                                                        experimental_data["gs_spin_string"], experimental_data)

    # recalculate all energies relative to the spin entered into fn.collate_energy_data() above
    data_to_plot["shifted_energies"] = shift_energy_levels(data_to_plot["all_energies"]) 

    data_to_plot["rms"] = calc_rms_err(10, data_to_plot["spin_1/2_energies"],
                        data_to_plot["spin_3/2_energies"], data_to_plot["spin_5/2_energies"], 
                        data_to_plot["spin_7/2_energies"], data_to_plot["spin_9/2_energies"], 
                        data_to_plot["spin_11/2_energies"], data_to_plot["spin_13/2_energies"])
    
    return data_to_plot


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


def calc_cbar_tick_labels(data, style):
    """
    A function to calculate a list of custom colour bar tick labels for a dataset. 

    Parameters
    ----------
    data : np.array of floats
        The data set to which these colour bar tick labels will apply.
    
    style : string
        "int" to format each label as an int
        "half" to format each label as a fraction of 2

    Returns
    -------
    cbar_tick_labels : list of strings
        The colour bar tick labels for this data.

    """
    
    cbar_ticks = np.arange(min(data), max(data)+1.0, 1.0)
    
    if style=="half":
        cbar_tick_labels = [spin_float_to_string(n) for n in cbar_ticks]
        
    elif style=="int":
        cbar_tick_labels = [int(n) for n in cbar_ticks]
        
    return cbar_tick_labels


def calc_cbar_ticks(contours):
    """
    A function to calculate a list of custom colour bar ticks for a dataset. 
    #!!! there is some kind of conflict between this and the func above...?

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
    
    if key[:2] == "jp":
        key = key[:-1] # remove the parity
    
    try:
        exp = inputs[key]
        tol = tolerance 
        
    except KeyError:
        exp = np.NaN
        tol = np.NaN
    
    return exp, tol


def calculate_format_data(output_item, name, experimental_data):
    '''
    Take a PropertyData object which has been initialised with data, and categorised
    with "num", "prop", and "sort" properties. Uses this info to determine/calculate
    graph plotting properties:
        - contour levels
        - colour bar ticks and tick labels
        - experimental data and tolerance
        

    Parameters
    ----------
    output_item : PropertyData object
        One property, which has been initialised as a PropertyData object but
        does not yet contain any information about graph plotting properties.
    
    name : string
        The name of the property.
        
    experimental_data : dictionary
        A dictionary of experimental data about the nucleus, input in the config file.

    Raises
    ------
    ValueError
        If an unregonsied property is input. This may be triggered when making extentsions
        to the code, to record more properties. Any new property will need to have its graph
        plotting properties hard coded in this function.

    Returns
    -------
    output_item : PropertyData object
        The same as the input objectm but now also including graph plotting information.

    '''

    if output_item.prop == "energies" and not(output_item.sort=="Fermi"): 
        
        output_item.contour_levels = 10
        output_item.cbar_tick_labels = 0
        
        if output_item.sort == "gap":
            
            output_item.experimental_data, output_item.error_tolerance = try_experimental(experimental_data, name, experimental_data["gap_en_tol"])
        
        elif output_item.sort == "Excited State ":
            
            output_item.experimental_data, output_item.error_tolerance = try_experimental(experimental_data, "x" + output_item.num + "_energy", experimental_data["abs_en_tol"])
        
        elif output_item.sort == "Spin ":
            
            output_item.experimental_data, output_item.error_tolerance = try_experimental(experimental_data, "jp_"+ output_item.num, experimental_data["abs_en_tol"])
            
        else: raise ValueError("property not recognised: " + name)
        
    elif output_item.prop == "mag_moments":
        
        output_item.contour_levels = 6
        output_item.cbar_tick_labels = 0
            
        if output_item.sort == "Excited State ":
            
            output_item.experimental_data, output_item.error_tolerance = try_experimental(experimental_data, "x" + output_item.num + "_mu", experimental_data["mu_tol"])
    
        elif output_item.sort == "Spin ":
            
            output_item.experimental_data = np.NaN
            output_item.error_tolerance = np.NaN
            
        elif output_item.sort == "Ground":
            
            output_item.experimental_data, output_item.error_tolerance = try_experimental(experimental_data, "gs_mu", experimental_data["mu_tol"])
            
        else: raise ValueError("property not recognised: " + name)
        
    elif output_item.prop == "quad_moments":
        
        output_item.contour_levels = 6
        output_item.cbar_tick_labels = 0
            
        if output_item.sort == "Excited State ":
            
            output_item.experimental_data, output_item.error_tolerance = try_experimental(experimental_data, "x" + output_item.num + "_mu", experimental_data["mu_tol"])
    
        elif output_item.sort == "Spin ":
            
            output_item.experimental_data = np.NaN
            output_item.error_tolerance = np.NaN
            
        elif output_item.sort == "Ground":
            
            output_item.experimental_data, output_item.error_tolerance = try_experimental(experimental_data, "gs_mu", experimental_data["mu_tol"])
            
        else: raise ValueError("property not recognised: " + name)
        
    
    elif output_item.sort == "Ground":
        
        if output_item.prop == "spin_floats": 
        
            output_item.contour_levels = calc_contour_levels(output_item.data)
            output_item.experimental_data, output_item.error_tolerance = try_experimental(experimental_data, "gs_spin_float", experimental_data["mu_tol"])
            output_item.cbar_tick_labels = calc_cbar_tick_labels(output_item.data, "half")
            output_item.cbar_ticks = calc_cbar_ticks(output_item.contour_levels)
            
        elif output_item.prop == "spin_strings": 
            output_item.contour_levels = 10
            output_item.cbar_ticks = 0
            output_item.cbar_tick_labels = 0
            output_item.experimental_data = np.NaN
            output_item.error_tolerance = np.NaN
            
        else: raise ValueError("property not recognised: " + name)
    
    elif output_item.sort == "Fermi":
    
        if output_item.prop == "indices": 
            
            output_item.contour_levels = calc_contour_levels(output_item.data)
            output_item.cbar_tick_labels = calc_cbar_tick_labels(output_item.data, "int")
            output_item.cbar_ticks = calc_cbar_ticks(output_item.contour_levels)
        
        elif output_item.prop == "energies":
        
            output_item.contour_levels = 10
            output_item.cbar_tick_labels = 0
        
        elif output_item.prop == "parities":
        
            output_item.contour_levels = 2
            output_item.cbar_tick_labels = 0
        
        else: raise ValueError("property not recognised: " + name)
            
        output_item.experimental_data = np.NaN
        output_item.error_tolerance = np.NaN
    
    elif output_item.prop == "delta":
        output_item.contour_levels = 10
        output_item.cbar_ticks = 0
        output_item.cbar_tick_labels = 0
        output_item.experimental_data = np.NaN
        output_item.error_tolerance = np.NaN
        
    else: raise ValueError("property not recognised: " + name)
    
    return output_item



def format_fig(polar_or_linear, ax, legend_handles, title, subtitle, **kwargs):
    """
    A function that handles formatting of a graph:
        - axis ranges
        - axis ticks and tick labels
        - axis labels
        - legend
        - graph title
        - subtitle

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
        
    subtitle : string
        The subtitle for the graph. May be empty.
        
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
        
        plt.xlabel("ε", size="xx-large")
        
        plt.xticks(size="large")
        #plt.yticks(size="xx-large") # doesn't seem to work for polar plots

        
        # gamma axis label
        ax.text(45*np.pi/180, ax.get_rmax()*1.2, "γ", ha='center', va='center', fontsize="xx-large") 
        
        if len(legend_handles) > 0:
            ax.legend(handles=legend_handles, loc="upper left", fontsize="x-large", 
                           bbox_to_anchor=(-0.5, 1.0))#, facecolor = '#D2D2D2', framealpha=0.9, )
        
        ax.set_title(title, va='bottom', y=1.1, fontsize="xx-large")  
        
        if subtitle != "":
            ax.text(0.05, 0.95, subtitle, transform=ax.transAxes, fontsize=12, #verticalalignment='top')
                    horizontalalignment = 'center', position=(0.5,1.06))
        
    elif polar_or_linear == 'linear':
        
        varied = kwargs.get("varied", None)
        x_label = kwargs.get("x_label", None)
        y_label = kwargs.get("y_label", None)
        legend_title = kwargs.get("legend_title", "")
    
        pad = 0.00*(varied[-1]-varied[0])
        ax.set_xlim([varied[0]-pad, varied[-1]+pad]) 
        
        #ax.margins(0) 
        
        plt.xlabel(x_label, size="xx-large")
        plt.ylabel(y_label, size="x-large")
        
        
        ax.set_title(title, va='bottom', y=1.1, fontsize="xx-large") 
        
        plt.xticks(size="large")
        plt.yticks(size="large")
        
        if "All " in title:
            legend1 = plt.legend(handles = legend_handles[len(legend_handles)//2:], loc="center left", bbox_to_anchor=(1.0, 0.75), fontsize="large", title_fontsize="x-large")
            legend2 = plt.legend(handles = legend_handles[:len(legend_handles)//2], loc="center left", bbox_to_anchor=(1.0, 0.25), fontsize="large", title_fontsize="large")
            plt.gca().add_artist(legend1)
            
            legend1.set_title(legend_title)
            legend2.set_title("experiment")
        else: 
            legend = ax.legend(handles = legend_handles, loc="upper center", bbox_to_anchor=(0.5, -0.15), fontsize="large", title_fontsize="x-large")
            
            if legend_title != "x":
                legend.set_title(legend_title)
        
        if subtitle != "":
            ax.text(0.05, 0.95, subtitle, transform=ax.transAxes, fontsize=10, horizontalalignment = 'center', #verticalalignment='top')
                    position=(0.5,1.03))
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
    
    # get the data to plot
    if prop.sort == "Spin ":
        
        plot_data = np.transpose(prop.data)[0]
        
        print('''multiple levels can't be plotted on a contour plot; 
              plotting only the yrast state of spin ''' + prop.num)
    else:
        plot_data = np.array(prop.data)
        
    # now start plotting
    
    # some of the data points may be NaN, 
    # so manually create the triangulation and mask NaN triangles.
    triang = tri.Triangulation(data_points["gamma_radians"], data_points["eps"])
    mask = np.any(np.isnan(plot_data[triang.triangles]), axis=1)
    triang.set_mask(mask)
    
    if prop.sort == "rms":
        mycmap = plt.get_cmap("viridis_r").copy()
        #mycmap.set_extremes(under='yellow', over='black', bad = 'red')
        cax = ax.tricontourf(triang, plot_data, levels=prop.contour_levels, cmap=mycmap, norm=colors.LogNorm(vmin=prop.range_min, vmax=500))
        cbar = plt.colorbar(cax, pad=0.1, extend='both')
    else:
        mycmap = plt.get_cmap("viridis").copy()
        cax = ax.tricontourf(triang, plot_data, levels=prop.contour_levels, cmap=mycmap)
        cbar = plt.colorbar(cax, pad=0.1)
        
    

    cbar.set_label(prop.axis_label, fontsize="x-large") 
    
    if prop.cbar_tick_labels:        # then format for discrete values
        
        cbar.set_ticks(prop.cbar_ticks)
        cbar.set_ticklabels(prop.cbar_tick_labels)
    
    cbar.ax.tick_params(labelsize="large")
        
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
        if len(data_points["file_tags"]) < 100: marker_size = 2
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
                      facecolor='None', label="data point that agrees with experiment")
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
            # exp_tol = cbar.ax.axhspan(prop.experimental_data[e]-prop.error_tolerance, prop.experimental_data[e]+prop.error_tolerance, facecolor='r', alpha=0.3, label="experimental tolerance")
           
    else:
        exp = cbar.ax.plot([0, 1], [prop.experimental_data, prop.experimental_data], 
                           'r-', label = "experimental value")
        # exp_tol = cbar.ax.axhspan(prop.experimental_data[e]-prop.error_tolerance, prop.experimental_data[e]+prop.error_tolerance, facecolor='r', alpha=0.3, label="experimental tolerance")

   
    legend_handles.append(exp[0])
    
    return legend_handles
    

def plot_points(data_points, prop, legend_handles, cbar, code_settings):
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

    code_settings : dictionary
        
    Returns
    -------
    legend_handles : list of object handles.
        A list of all the objects so-far plotted on the graph which should be 
        included in the legend, now including the objects plotted in this function.
        
    
    """
    
    if np.isfinite(prop.experimental_data).all() and code_settings["mark_exp"]==1:
    
        legend_hit = False
        legend_miss = False
        
        for r in range(len(data_points["eps"])):
            
            # use a smaller marker size for large data sets
            if len(data_points["file_tags"]) < 100: marker_size = 10
            else: marker_size = 1
            
            if prop.sort == "Spin ":
                error = [abs(prop.data[r][0] - exp) for exp in prop.experimental_data]
                match = [err < prop.error_tolerance for err in error]
            else:
                error = abs(prop.data[r] - prop.experimental_data)
                match = [error < prop.error_tolerance]
                
            if any(match): 
                # update the record of how many of the tested properties agree
                data_points["agreed"][r] += 1      

                if code_settings["mark_points"]==1:
                    
                    hit = plt.scatter(data_points["gamma_radians"][r], 
                          data_points["eps"][r], s=marker_size, #edgecolor='red', 
                          facecolor='red', label="data point that agrees \nwith experiment")
                    
                    legend_hit = True
                
            # for data points that don't match experimental data, only plot them 
            # when the data set is quite small, to avoid cluttering the graph.
            elif len(data_points["file_tags"]) < 100 and code_settings["mark_points"]==1: 
                miss, = plt.polar(data_points["gamma_radians"][r], 
                          data_points["eps"][r], 'wx', label="does not match experiment")
                legend_miss = True
        
        if legend_hit:
            legend_handles.append(hit)
        if legend_miss:
            legend_handles.append(miss)
        
        if code_settings["mark_exp"]==1:
            if prop.sort == "Spin ":
                for e in range(len(prop.experimental_data)):
                    exp = cbar.ax.plot([0, 1], 
                                       [prop.experimental_data[e], prop.experimental_data[e]], 
                                       'r-', label = "experimental value", linewidth=5)
                    # exp_tol = cbar.ax.axhspan(prop.experimental_data[e]-prop.error_tolerance, prop.experimental_data[e]+prop.error_tolerance, facecolor='r', alpha=0.3, label="experimental tolerance")
                    
            else:
                exp = cbar.ax.plot([0, 1], [prop.experimental_data, prop.experimental_data], 
                                   'r-', label = "experimental value", linewidth=5)
                # exp_tol = cbar.ax.axhspan(prop.experimental_data[e]-prop.error_tolerance, prop.experimental_data[e]+prop.error_tolerance, facecolor='r', alpha=0.3, label="experimental tolerance")
                
       
        legend_handles.append(exp[0])
    '''
    else: 
        # use a smaller marker size for large data sets
        if len(data_points["file_tags"]) < 100: marker_size = 5
        else: marker_size = 1
        
        if inputs["mark_points"]==1:
            all_points = plt.scatter(data_points["gamma_radians"], 
                      data_points["eps"], s=marker_size, c='w', label="data point")
            legend_handles.append(all_points) 
    '''
    
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


def assign_parameters(ptrm_inputs, data_points):
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

    if ptrm_inputs["deformation_input"] == "eps" :
        
        var_sym = "ε"
        var = data_points["eps"]
        fix_sym = "γ"
        fix = data_points["gamma_degrees"]
        
    elif ptrm_inputs["deformation_input"] == "gamma" :
        
        var_sym = "γ / º"
        var = data_points["gamma_degrees"]
        fix_sym = "ε"
        fix = data_points["eps"]
        
    elif len(data_points["e2plus"]) > 1:
        
        var_sym = r"$E(2^+)$ / MeV"
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
    line_colours = ['k-x', 'b-x', 'y-x', 'c-x', 'm-x']
    line_labels = ["lowest energy", "second lowest energy", "third lowest energy", "fourth lowest energy", "fifth lowest energy"]
    
    for s in range(min(len(line_labels), np.size(data_by_line,0))):
        
        data, = plt.plot(var, data_by_line[s], line_colours[s], label=line_labels[s], markersize=marker_size)
        legend_handles.append(data)
        
    legend_title = "%s = %s" % (fix_sym, fix_val)
    
    return legend_handles, legend_title

def plot_all_energies(prop, var, legend_handles, marker_size, fix_sym, fix_val):
    """
    A function to plot energy levels of all spins states calculated.
    
    Parameters
    ----------
    prop : PropertyData object
        The energy levels, collected into a class with 
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
    num_lines = np.size(prop.data, axis=1)
    
    
    line_colours = ['y', 'b', 'm', 'k', 'c', 'g', 'r']
    line_labels = ["1/2", "3/2", "5/2", "7/2", "9/2", "11/2", "13/2"]
    
    line_styles = ['', '-', '', '', '-', '', '-']
    
    idx = [int(((2*a)-1)/2) for a in prop.spins]
    used = []
    
    for s in range(num_lines): 
        if idx[s] in used:
            l = "_hidden"
            #continue #!!!
        else: 
            l = line_labels[idx[s]]
            used.append(idx[s])
            
        if line_styles[idx[s]] == "":
            continue
        
        data, = plt.plot(var, data_by_line[s], line_colours[idx[s]]+line_styles[idx[s]], label=l, markersize=marker_size, linewidth=0.75) # line_colours[s] // 'kx'
        legend_handles.append(data)
        
    legend_title = "%s = %s" % (fix_sym, fix_val)
    
    
    ''
    xrange = [min(var),max(var)]
    experimental_data = prop.experimental_data
    explabs = prop.explabels
    
    expidx = [int(((2*spin_string_to_float(a))-1)/2) for a in explabs]
    
    for e in range(len(experimental_data)):
        expl = plt.plot(xrange, [experimental_data[e],experimental_data[e]], line_colours[expidx[e]]+"--", label=explabs[e], linewidth=1.25)
        legend_handles.append(expl[0])
    
    '''
    exp_tol_1 = plt.plot(xrange, [experimental_data[0],experimental_data[0]], label="1/2")
    exp_tol_3  = plt.plot(xrange, [experimental_data[0],experimental_data[0]], label="3/2")
    exp_tol_5 = plt.plot(xrange, [experimental_data[0],experimental_data[0]], label="5/2")
    exp_tol_7 = plt.axhspan(xrange, [experimental_data[1],experimental_data[1]], label="7/2")
    exp_tol_9 = plt.axhspan(xrange, [experimental_data[2],experimental_data[2]], label="9/2")
    exp_tol_11 = plt.axhspan(xrange, [experimental_data[3],experimental_data[3]],label="11/2")
    exp_tol_13 = plt.axhspan(xrange, [experimental_data[4],experimental_data[4]],  label="13/2")

    legend_handles.append(exp_tol_1)
    legend_handles.append(exp_tol_3)
    legend_handles.append(exp_tol_5)
    legend_handles.append(exp_tol_7)
    legend_handles.append(exp_tol_9)
    legend_handles.append(exp_tol_11)
    legend_handles.append(exp_tol_13)
    '''
    
    
    return legend_handles, legend_title


def mark_spin(ptrm_inputs, data_points, spin_data, legend_handles, ax):
    """
    A function to plot the region(s) of correct ground state spin onto a polar plot,
    as a black contour line.

    Parameters
    ----------
    ptrm_inputs : dictionary
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

    correct_range = [ptrm_inputs["gs_spin_float"]-0.5, 
                     ptrm_inputs["gs_spin_float"]+0.5]
    spin_colour = (0,0,0) #(213/255,1,0)
    
    ax.tricontour(data_points["gamma_radians"], data_points["eps"], 
                  spin_data, levels=correct_range,  
                  colors=[spin_colour], linewidths=1.0)
    spin_legend_proxy = plt.Line2D([], [], color=spin_colour, linewidth=1.0, label="region of correct g.s. spin") 
    legend_handles.append(spin_legend_proxy)

    return legend_handles
  
       
def plot_line_data(data_points, prop, var, fix_sym, fix, legend_handles):
    '''
    Plot the line and "x" data point markers on a line graph. May be multiple lines,
    if the property is for all calculated states (of a given spin).

    Parameters
    ----------
    data_points : dictionary
        
    prop : PropertyData object
        The property being plotted.
        
    var : list of floats
        The independent variable.
        
    fix_sym : string
        The name/symbol of the controlled variable. Usually either ε or γ.
        
    fix : float
        The value of the controlled variable.
        
    legend_handles : list of handles


    Returns
    -------
    legend_handles : list of handles

    legend_title : string
        

    '''

    # now plot the actual data
    if len(data_points["file_tags"]) < 100: marker_size = 5
    else: marker_size = 1 # use smaller markers if the data set is large
    
    if len(fix)==0:
        fix= ["(%.3f, %.1f)"  % (data_points["eps"][0],data_points["gamma_degrees"][1] )]
        
    if prop.sort == "Spin ":
        
        if prop.num == "all":
            legend_handles, legend_title = plot_all_energies(prop, var, legend_handles, marker_size, fix_sym, fix[0])
        else:
            legend_handles, legend_title = plot_multi_lines(prop, var, legend_handles, marker_size, fix_sym, fix[0])
        
    else:
        
        data, = plt.plot(var, prop.data, 'k-x', markersize=marker_size, label="%s = %s" % (fix_sym, fix[0]))
        legend_handles.append(data)
        legend_title = ""
        
    return legend_handles, legend_title

def plot_exp_line(prop, code_settings, var, legend_handles):
    '''
    A function for plotting a red line on a line graph to indicate the experimental value.
    Can include a shaded region to indicate tolerance if requested.

    Parameters
    ----------
    prop : PropertyData object
        The property being plotted.
        
    code_settings : dictionary
        Dictionary of inputs from config file.
        
    var : list of floats
        The independent variable (usually either eps or gamma)
        
    legend_handles : list of handles
        A list of handles of objects already drawn on the graph, to be included in the legend.

    Returns
    -------
    legend_handles : list of handles
        Same as input, now including the experimental line.

    '''
    if prop.sort == "Spin ":
        for _ in range(len(prop.experimental_data)):
            if code_settings["mark_exp"]:
                exp, = plt.plot(var, np.full(len(var), prop.experimental_data[_]), 'r-', label="experimental value")
            if code_settings["mark_exp_tol"]:
                exp_tol = plt.axhspan(prop.experimental_data[_]-prop.error_tolerance, prop.experimental_data[_]+prop.error_tolerance, facecolor='r', alpha=0.2, label="experimental value")
    else:
        if code_settings["mark_exp"]:
            exp, = plt.plot(var, np.full(len(var), prop.experimental_data), 'r-', label="experimental value")
        if code_settings["mark_exp_tol"]:
            exp_tol = plt.axhspan(prop.experimental_data-prop.error_tolerance, prop.experimental_data+prop.error_tolerance, facecolor='r', alpha=0.2, label="experimental value")
    
    if code_settings["mark_exp"]:
        legend_handles.append(exp)
    if code_settings["mark_exp_tol"]:
        legend_handles.append(exp_tol)
        
    return legend_handles
