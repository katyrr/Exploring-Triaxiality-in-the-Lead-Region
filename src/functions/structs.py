#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  7 10:58:19 2025

@author: katyrr

A separate module file to contain structs such as classes and read-only data.


"""

import time                                                                    

class Timer():
    """
    A class that acts as a timer/stopwatch.

    Attributes
    ----------
    is_running : bool
        A record of whether the timer is currently running.
        
    start_time : float
        The time at which the timer was constructed.
        
    end_time : float
        The time at which end_timer() was called.
        
    Methods
    -------
    Timer():
        Constructor for the timer object.
        
    start():
        Starts the timer.
        
    stop():
        Ends the timer.
        
    get_lapsed_time():
        Returns the timeframe (in seconds) most recently recorded on the timer.
        
        
    Raises
    ------
    RuntimeError:
        Occurs if the timer is stopped before it is started, 
        or started while it is already running,
        or if the lapsed time is requested while it is still running.
    
    """
    
    def __init__(self):
        self.is_running = False
        
    def start(self):
        #if self.is_running:
        #    raise RuntimeError("The timer is already running.")
        self.start_time = time.time()
        self.is_running = True
        
        
    def stop(self):
        #if not(self.is_running):
        #    raise RuntimeError("The timer has not been started yet.")
        self.end_time = time.time()
        self.is_running = False
        
    def get_lapsed_time(self):
        #if self.is_running:
        #    raise RuntimeError("The timer is still running.")
        total_time= self.end_time - self.start_time
        return total_time
    


def get_restricted_inputs():

    
    restricted_input_values = {
        'par': ['+', '-'],
        'irec': [0,1],
        'icorr': [0,1],
        'istrch': [0,1],
        'iq': [0,1],
        'vmi': 	[0,1],
        'OS': ['MacOS', '64bit'],
        'detailed_print': [0,1],
        'mark_spin': [0,1],
        'mark_exp': [0,1],
        'mark_exp_tol': [0,1],
        'mark_points': [0,1],
        'include_subtitle': [0,1]
        }
    
    return restricted_input_values


def get_required_inputs():
    
    required_inputs = ['par', 'nucleus', 'Z','A', 'cutoff', 'gsfac', 'irec', 
                       'icorr', 'emin', 'emax', 'chsi', 'eta', 'nantj', 'noutj', 
                       'ipout', 'nu', 'imin', 'ispin', 'kmax', 'istrch', 'num_orbs', 
                       'nuu', 'nprot', 'nneutr', 'e2plur', 'ispec', 'iq', 'gr', 
                       'vmi', 'nmin', 'nmax', 'OS', 'num_cores', 'figure_res', 
                       'detailed_print', 'mu_tol', 'abs_en_tol', 'gap_en_tol', 
                       'mark_spin', 'mark_exp', 'mark_exp_tol', 'mark_points',
                       'e2plus', 'include_subtitle']
    
    return required_inputs


def get_variable_list(var_type):
    """

    Parameters
    ----------
    var_type : string
        The name of the type of variable being sought.
        Allowed values: "bool", "int", "float".

    Raises
    ------
    ValueError
        Occurs when the 'var_type' parameter is not one of bool/int/float.

    Returns
    -------
    var_list : list of strings
        A list of the variable inputs whose values are expected to have the
        requested type, and should be converted to this type when being read.

    """
    
    if var_type == "bool":
        var_list = ["mark_spin", "mark_exp", "mark_exp_tol", "mark_points", 
                    "detailed_print", "include_subtitle"]

    elif var_type == "int":
        var_list = ["A", "Z", "num_to_record", "num_orbs", "nu", "num_cores", "figure_res",
                    "irec", "icorr", "imin", "ispin", "kmax", "istrch", "nuu", "nprot", 
                    "nneutr", "ispec", "iq", "gr", "vmi", "nmin", "nmax"]

    elif var_type == "experimental_float":
        var_list = ["gs_energy", "gs_mu", "jp_", "engap", "mu_tol", "abs_en_tol", "gap_en_tol"]
        
    elif var_type == "settings_float":
        var_list = ["cutoff", "gsfac", "emin", "emax", "chsi", "eta", "e2plur"]
        
    elif var_type == "string":
        var_list = ["nucleus", "nantj", "noutj", "ipout", "OS", "par"]
        
    else: raise ValueError("unrecognised type: " + var_type)
    
    return var_list



def get_template(program):
    """
    A function to get a string that lays out the format of a .DAT file for the
    requested program, with string formatting inserts for filling values from
    a dictionary with the required keys.

    Parameters
    ----------
    program : string 
        The name of the program whose template should be returned.
        Allowed values: "gampn", "asyrmo", "probamo".

    Raises
    ------
    ValueError
        Occurs if the 'program' parameter is not one of gampn/aysrmo/probamo.

    Returns
    -------
    template : string
        The format for the .DAT file of the requested program.
        Contains string formatting inserts %(key)s in place of values.
        
    """
    
    if program == "gampn":
        
        template = '''
'%(current_f002)s' '%(current_f016)s' '%(current_f017)s'     file2, file16, file17
%(istrch)s,%(icorr)s,%(irec)s                          ISTRCH,ICORR,irec
9,%(nneupr)s                            NKAMY,NNEUPR
0.120,0.00,0.120,0.00
0.120,0.00,0.120,0.00
0.105,0.00,0.105,0.00
0.090,0.30,0.090,0.25
0.065,0.57,0.070,0.39
0.060,0.65,0.062,0.43        'STND'  KAPPA, MU (PROTON)
0.054,0.69,0.062,0.34        'STND'         MU (PROTON)
0.054,0.69,0.062,0.26        'STND'         MU (PROTON)
0.054,0.69,0.062,0.26        'STND'         MU (PROTON)
0
0.054,0.69,0.062,0.26
%(nuu)s,1,1,0                       NUU,IPKT,NOYES,ITRANS
%(emin)s,%(emax)s
%(current_orbitals)s                 fermi_parityP, NORBITP, LEVELP
%(current_orbitals)s                  fermi_parityN, NORBITN, LEVELN
%(Z)s,%(A)s                                                Z,A
%(current_eps)s,%(current_gamma)s,0.00,0.0,0.0000,8,8,0,0
(LAST CARD: EPS,GAMMA,EPS4,EPS6,OMROT,%(nprot)s,%(nneutr)s,NSHELP,NSHELN)
'''
        
    elif program == "asyrmo":
        
        template = '''
'%(current_f016)s' '%(current_f017)s' '%(current_f018)s' FILE16,FILE17,FILE18
1,0                                        IPKT,ISKIP
%(istrch)s,%(irec)s                                        ISTRCH,IREC
%(vmi)s,%(nmin)s,%(nmax)s,8,0.0188,100.00                      VMI,NMIN,NMAX,IARCUT,A00,STIFF
%(Z)s,%(A)s,%(imin)s,%(ispin)s,%(kmax)s,%(current_e2plus)s,%(e2plur)s                   Z,AA,IMIN,ISPIN,KMAX,E2PLUS,E2PLUR
19.2,7.4,15,%(chsi)s,%(eta)s                     GN0,GN1,IPAIR,CHSI,ETA
%(current_orbitals)s  
  %(nantj)s  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
  %(noutj)s  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
  %(ipout)s  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  IPOUT(I)
'''
        
    elif program ==  "probamo":
        
        template = '''
'%(current_f017)s' '%(current_f018)s'              FILE17,FILE18
1,0                               ipkt,iskip
%(istrch)s                                 isrtch
%(Z)s,%(A)s                           Z,AA
%(ispec)s,%(cutoff)s,%(iq)s,%(gsfac)s,%(gr)s                ISPEC,CUTOFF,IQ,GSFAC,GR
0.0000, 0.000,0.000, 0.000        BS2,BS4 (FOR S-STATE), BS2,BS4(P-STATE)
'''
        
    else: raise ValueError("unrecognised program: " + program)
    
    return template



class PropertyData:
    """
    A class collects information about a data set, for plotting a graph.

    Attributes
    ----------
    axis_label : string
        The text that should label the axis that this data is plotted on.
        Includes units (if applicable).
        
    cbar_tick_labels : list of strings, or 0
        Define custom tick labels on the colour bar of a filled contour plot.
        If 0, use pyplot's default tick labels. 
        
    cbar_ticks : list of floats, or 0
        Define custom placement of ticks on the colour bar of a filled contour plot.
        If 0, use pyplot's default tick placement. 
        
    contour_levels : np.array of floats, or int
        Define custom boundaries of contour levels.
        If int, use that number of equally spaced contour levels between the 
        max and min data value.
        
    data : list of floats, or np.array of floats
        The value of this nuclear property at each data point.
        Must have length = number of data points.
        Properties grouped by spin (e.g. spin_1/2_energies) have 
        depth = max number of levels found at that spin.
    
    error_tolerance : float
        The absolute error tolerance to allow when comparing to an experimental value.
        np.NaN if no experimental value has been input.
        
    eperimental_data : float
        The experimental value of this nuclear property.
        np.NaN if no experimental value has been input.
        
    num : string
        Used to categorise this nucleaar property.
        For properties that are grouped by spin, this is the spin of the group (e.g. "1/2").
        For properties that are grouped by association with experimental data, 
        this is the index of the energy level (e.g. "1" for the first excited state).
        Otherwise "0".
        
    plot : bool
        True if this data should be plotted.
        False if this data should not be plotted.
        
    prop : string
        Used to categorise this nuclear property.
        Describes what kind of data it is (e.g. "energies", "mag_moments", "spin_floats").
        
    sort : string
        Used to categorise this nuclear property.
        Describes which group the data belongs to (e.g. "Ground", "Excited State ", "Spin ").
    
        
    Methods
    -------
    PropertyData(data, name):
        Constructor for the PropertyData object.
        Saves "data" in the "data" attribute.
        Uses the "name" parameter to categories the data with "prop", "sort", "num" attributes.
        Sets the "title" and "axis_label" attributes.
        
        Parameters
        ----------
        
        data : list of floats, or np.array of floats
            The value of this nuclear property at each data point.
            Must have length = number of data points.
            Properties grouped by spin (e.g. spin_1/2_energies) have 
            depth = max number of levels found at that spin.
        
        name : string
            The key that was previously holding this data in a dictionary.
            Encodes information about "prop" and "sort", and may additionally contain "num".
            
      
        Raises
        ------
        ValueError:
            Occurs if the property cannot be categorised with the current implementation.
            Any future new properties must be hardcoded so that they can be 
            correctly categorised and plotted.
    
    """
    
    plot = False # don't plot by default
    
    def __init__(self, data, name): 
        
        self.data = data
        
        self.num, self.prop, self.sort = identify(name)
        
        self.title, self.axis_label = set_labels(self.num, self.prop, self.sort, name)
        


def identify(name):
    
    # work out what kind of property it is from its name
    if name[0:4] == "spin":
        num = name[5:9]
        if num[-1] == "_": num = num[0:3]
        prop = name[9:] 
        if prop[0] == "_": prop = prop[1:]
        sort = "Spin "
        
    elif name[0] == "x": # obsolete
        num = name[1]
        prop = name[3:]
        sort = "Excited State "
        
    elif name[0:2] == "gs":
        num = ""
        prop = name[3:]
        sort = "Ground"
    
    elif name[0:5] == "fermi":
        num = ""
        prop = name[6:14]
        sort = "Fermi"
        
    elif "Agreement" in name:
        num = ""
        prop = ""
        sort = ""
        
    elif name[0:5] == "engap":
        num = name[6:].replace("_", " ").strip().replace(" ", ") and ").replace(".", " (#") + ")"
        prop = "energies"
        sort = "gap"
        
    elif name == "RMS energies":
        num = ""
        prop = ""
        sort = "rms"
        
    elif "All E" in name:
        sort = "Spin "
        num = "all"
        prop = ""
    
    elif name == "delta":
        prop = "delta"
        num = ""
        sort = ""
    
    else: 
        raise ValueError("property not regonised: " + name)
    
    return num, prop, sort

        
def set_labels(num, prop, sort, name):

    # set graph title and axis label strings
    if sort == "gap":
        
        title = "Energy Gaps Between Spin States "+ num
        axis_label = "Energy Gap / keV"
        
    elif prop == "energies" and not(sort=="Fermi"): 
        
        title = "Energies of "+ sort + num + " States"
        axis_label = title + " / keV"
    
    elif sort == "rms":
         
        title = "RMS error in energies"
        axis_label = title + " / keV"
       
    elif prop == "mag_moments":
            
        if sort == "Excited State ": 
            title = "Magnetic Dipole Moment of Excited State " + num 
            
        elif sort == "Spin ": 
            title = "Magnetic Dipole Moment of Spin " + num + " States"
            
        elif sort == "Ground": 
            title = "Ground State Magnetic Dipole Moment"
            
        else: raise ValueError("property not recognised: " + name)
        
        axis_label = title + r' / $μ_{N}$'
    
    elif prop == "quad_moments":
            
        if sort == "Excited State ": 
            title = "Electric Quadrupole Moment of Excited State " + num 
            
        elif sort == "Spin ": 
            title = "Electric Quadrupole Moment of Spin " + num + " States"
            
        elif sort == "Ground": 
            title = "Ground State Electric Quadrupole Moment"
            
        else: raise ValueError("property not recognised: " + name)
        
        axis_label = title + r' / $eb$'
        
    elif sort == "Ground":
        
        title = "Ground State Spin"
        axis_label = title
    
    elif sort == "Fermi":
    
        if prop == "indices": 
            title = "Fermi Level Index"
            axis_label = title
            
        elif prop == "energies":
        
            title = "Fermi Energy"
            axis_label = title + name[-3:-1]
            
        elif prop == "parities":
        
            title = "Fermi Parity"
            axis_label = title
            
        else: raise ValueError("property not recognised: " + name)
    
    elif name == "Agreement of Data Points With Experimental Data":
        title = name
        axis_label = "Number of Matches"
        
    elif "All E" in name:
        title = name
        axis_label = "Energy / keV"
    
    elif name == "delta":
        title = "Pairing Gap Energy Δ"
        axis_label = "Δ / MeV"
        
    else:
        raise ValueError("property not recognised: " + name)
        
    return title, axis_label
