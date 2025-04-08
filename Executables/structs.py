#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  7 10:58:19 2025

@author: katyrr

A separate module file to contain structs such as classes and read-only data.

CONTENTS:
--------

- class Timer()

- template = get_template(program)

- var_list = get_variable_list(var_type)

- class PropertyData(data, name)

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
%(current_gampn_orbitals)s                 fermi_parityP, NORBITP, LEVELP
%(current_gampn_orbitals)s                  fermi_parityN, NORBITN, LEVELN
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
        var_list = ["mark_spin", "mark_exp"]

    elif var_type == "int":
        var_list = ["A", "Z", "num_to_record", "num_orbs", "nu"]

    elif var_type == "float":
        var_list = ["gs_energy", "x1_energy", "x2_energy", "x3_energy", "gs_mu", "x1_mu"]
        
    else: raise ValueError("unrecognised type: " + var_type)
    
    return var_list


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
        
        # work out what kind of property it is from its name
        if name[0:4] == "spin":
            num = name[5:9]
            if num[-1] == "_": num = num[0:3]
            prop = name[9:] 
            if prop[0] == "_": prop = prop[1:]
            sort = "Spin "
            
        elif name[0] == "x":
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
            
        elif (name == "Agreement of Data Points With Experimental Data"
              or name == "All Energy Levels"
              or name == "delta"):
            num = ""
            prop = ""
            sort = ""
        
        else: 
            raise ValueError("property not regonised: " + name)
        
        self.num = num
        self.prop = prop
        self.sort = sort
        
        # set graph title and axis label strings
        if prop == "energies" and not(sort=="Fermi"): 
            
            self.title = "Energies of "+ sort + num + " States"
            self.axis_label = self.title + " / keV"
           
        elif prop == "mag_moments":
                
            if sort == "Excited State ": 
                self.title = "Magnetic Dipole Moment of Excited State " + num 
                
            elif sort == "Spin ": 
                self.title = "Magnetic Dipole Moment of Spin " + num + " States"
                
            elif sort == "Ground": 
                self.title = "Ground State Magnetic Dipole Moment"
                
            else: raise ValueError("property not recognised: " + self.name)
            
            self.axis_label = self.title + r' / $μ_{N}$'
            
        elif sort == "Ground":
            
            self.title = "Ground State Spin"
            self.axis_label = self.title
        
        elif sort == "Fermi":
        
            if prop == "indices": 
                self.title = "Fermi Level Index"
                self.axis_label = self.title
                
            elif prop == "energies":
            
                self.title = "Fermi Energy"
                self.axis_label = self.title + name[-3:-1]
                
            elif prop == "parities":
            
                self.title = "Fermi Parity"
                self.axis_label = self.title
                
            else: raise ValueError("property not recognised: " + name)
        
        elif name == "Agreement of Data Points With Experimental Data":
            self.title = name
            self.axis_label = "Number of Matches"
            
        elif name == "All Energy Levels":
            self.title = name
            self.axis_label = "Energy / keV"
            self.sort = "Spin "
            self.num = "all"
        
        elif name == "delta":
            self.title = "Pairing Gap Energy Δ"
            self.prop = "delta"
            self.axis_label = "Δ / MeV"
            
            
        else:
            raise ValueError("property not recognised: " + name)
        
        
        
        