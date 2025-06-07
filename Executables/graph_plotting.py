#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  8 23:07:07 2025

@author: katyrr


CONTENTS:
--------

- contour_levels = calc_contour_levels(data)

- cbar_tick_labels = calc_cbar_tick_labels(data, style)

- cbar_ticks = calc_cbar_ticks(contours)

- exp, tol = try_experimental(inputs, key, tolerance)
    
- void = format_fig(polar_or_linear, ax, legend_handles, title, subtitle, **kwargs)

- cax, cbar = draw_contour_plot(ax, prop, data_points)

- correct_spin_range = find_correct_spin(gs_spins, experimental_value)

- correct_spin_handle = plot_correct_spin(correct_spin_range, var, step, prop)

- legend_handles = plot_points_with_experiment(data_points, prop, legend_handles, cbar)

- legend_handles = plot_points_without_experiment(data_points, legend_handles)

- var_sym, var, fix_sym, fix = assign_parameters(inputs, data_points)
    
- legend_handles, legend_title = plot_multi_lines(prop, var, legend_handles, marker_size, fix_sym, fix_val)

- legend_handles = mark_spin(inputs, data_points, output_data, legend_handles, ax)

"""




import matplotlib.pyplot as plt                 # for plotting graphs
from matplotlib.ticker import FuncFormatter     # for formatting axis ticks
import matplotlib.tri as tri                    # for manual triangulation before drawing a contour plot
import matplotlib.colors as colors
import numpy as np                              # for np.arrays
import structs as st                            # my own module file of structs (classes, and read-only dicts)
import functions as fn                          # my own module file of functions




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
        cbar_tick_labels = [fn.spin_float_to_string(n) for n in cbar_ticks]
        
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


def calculate_format_data(output_item, name, experimental):

    if output_item.prop == "energies" and not(output_item.sort=="Fermi"): 
        
        output_item.contour_levels = 10
        output_item.cbar_tick_labels = 0
        
        if output_item.sort == "gap":
            
            output_item.experimental_data, output_item.error_tolerance = try_experimental(experimental, name, experimental["gap_en_tol"])
        
        elif output_item.sort == "Excited State ":
            
            output_item.experimental_data, output_item.error_tolerance = try_experimental(experimental, "x" + output_item.num + "_energy", experimental["abs_en_tol"])
        
        elif output_item.sort == "Spin ":
            
            output_item.experimental_data, output_item.error_tolerance = try_experimental(experimental, "jp_"+ output_item.num, experimental["abs_en_tol"])
            
        else: raise ValueError("property not recognised: " + name)
        
    elif output_item.prop == "mag_moments":
        
        output_item.contour_levels = 6
        output_item.cbar_tick_labels = 0
            
        if output_item.sort == "Excited State ":
            
            output_item.experimental_data, output_item.error_tolerance = try_experimental(experimental, "x" + output_item.num + "_mu", experimental["mu_tol"])
    
        elif output_item.sort == "Spin ":
            
            output_item.experimental_data = np.NaN
            output_item.error_tolerance = np.NaN
            
        elif output_item.sort == "Ground":
            
            output_item.experimental_data, output_item.error_tolerance = try_experimental(experimental, "gs_mu", experimental["mu_tol"])
            
        else: raise ValueError("property not recognised: " + name)
        
    elif output_item.prop == "quad_moments":
        
        output_item.contour_levels = 6
        output_item.cbar_tick_labels = 0
            
        if output_item.sort == "Excited State ":
            
            output_item.experimental_data, output_item.error_tolerance = try_experimental(experimental, "x" + output_item.num + "_mu", experimental["mu_tol"])
    
        elif output_item.sort == "Spin ":
            
            output_item.experimental_data = np.NaN
            output_item.error_tolerance = np.NaN
            
        elif output_item.sort == "Ground":
            
            output_item.experimental_data, output_item.error_tolerance = try_experimental(experimental, "gs_mu", experimental["mu_tol"])
            
        else: raise ValueError("property not recognised: " + name)
        
    
    elif output_item.sort == "Ground":
        
        if output_item.prop == "spin_floats": 
        
            output_item.contour_levels = calc_contour_levels(output_item.data)
            output_item.experimental_data, output_item.error_tolerance = try_experimental(experimental, "gs_spin_float", experimental["mu_tol"])
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
    

def plot_points(data_points, prop, legend_handles, cbar, inputs):
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
    
    if np.isfinite(prop.experimental_data).all() and inputs["mark_exp"]==1:
    
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

                if inputs["mark_points"]==1:
                    
                    hit = plt.scatter(data_points["gamma_radians"][r], 
                          data_points["eps"][r], s=marker_size, #edgecolor='red', 
                          facecolor='red', label="data point that agrees \nwith experiment")
                    
                    legend_hit = True
                
            # for data points that don't match experimental data, only plot them 
            # when the data set is quite small, to avoid cluttering the graph.
            elif len(data_points["file_tags"]) < 100 and inputs["mark_points"]==1: 
                miss, = plt.polar(data_points["gamma_radians"][r], 
                          data_points["eps"][r], 'wx', label="does not match experiment")
                legend_miss = True
        
        if legend_hit:
            legend_handles.append(hit)
        if legend_miss:
            legend_handles.append(miss)
        
        if inputs["mark_exp"]==1:
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
    
    expidx = [int(((2*fn.spin_string_to_float(a))-1)/2) for a in explabs]
    
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
    #!!! missing docs

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
        
        if max_agreement > 0:
            print("\nPoints with agreement = " + str(max_agreement) + ":")
            
            i = len(unique_values)-1
            
            print("\n\tAgreement = " + str(unique_values[i]) + ":")
        
            lower = sum(counts[0:i])
            upper = sum(counts[0:i+1])
            
            these_eps = sorted_eps[lower:upper]
            these_gamma = sorted_gamma[lower:upper]
            
            r = min(len(these_eps), 5)
            for j in range(r):
                print("\t\t(ε, γ) = (" + str(these_eps[j]) + ",\t" + str(these_gamma[j])+"º)")
            
            if r==5 and len(these_eps) != 5:
                print("... etc, " + str(len(these_eps)))
        
       
def plot_line_data(data_points, prop, var, fix_sym, fix, legend_handles):

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

def plot_exp_line(prop, inputs, var, legend_handles):
    if prop.sort == "Spin ":
        for _ in range(len(prop.experimental_data)):
            if inputs["mark_exp"]:
                exp, = plt.plot(var, np.full(len(var), prop.experimental_data[_]), 'r-', label="experimental value")
            if inputs["mark_exp_tol"]:
                exp_tol = plt.axhspan(prop.experimental_data[_]-prop.error_tolerance, prop.experimental_data[_]+prop.error_tolerance, facecolor='r', alpha=0.2, label="experimental value")
    else:
        if inputs["mark_exp"]:
            exp, = plt.plot(var, np.full(len(var), prop.experimental_data), 'r-', label="experimental value")
        if inputs["mark_exp_tol"]:
            exp_tol = plt.axhspan(prop.experimental_data-prop.error_tolerance, prop.experimental_data+prop.error_tolerance, facecolor='r', alpha=0.2, label="experimental value")
    
    if inputs["mark_exp"]:
        legend_handles.append(exp)
    if inputs["mark_exp_tol"]:
        legend_handles.append(exp_tol)
        
    return legend_handles


def plot_agreement(inputs, agreement, data_points, output_data, subtitle):
    
    inputs["current_graph"] = agreement.title
    print("plotting graph: %(current_graph)s" % inputs) 
    
    _fig, _ax = plt.subplots(subplot_kw=dict(projection='polar'))
    cax, cbar = draw_contour_plot(_ax, agreement, data_points)
    
    _legend_handles = []
    
    if inputs["mark_spin"]:
        _legend_handles = mark_spin(inputs, data_points, output_data["gs_spin_floats"].data, _legend_handles, _ax)
    
    if inputs["mark_points"]:
        _legend_handles =  plot_points_without_experiment(data_points, _legend_handles)
           
    
    format_fig('polar', _ax, _legend_handles, '%(current_graph)s of %(nucleus)s' % inputs, subtitle)
    
    plt.show()
    
    