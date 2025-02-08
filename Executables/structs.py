#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  7 10:58:19 2025

@author: katyrr

A separate module file to contain structs such as classes and read-only dictionaries.


"""

import time                                                                     # for checking how long it took to run

class timer():
    """
    A class that acts as a timer/stopwatch.
    Can time two processes at once: 
      - An "outer" process which lasts for the whole duration of the timer.
        Starts and stops once.
      - An "inner" process which lasts for a fraction of the total duration 
        of the timer. Can be started and stopped multiple times.


    Attributes
    ----------
    outer_is_running : bool
        A record of whether the outer timer is currently running.
        
    inner_is_running : bool
        A record of whether the inner timer is currently running.
        Default value = False, until the inner timer is started.
        
    start_time : float
        The time at which the timer was constructed.
        
    end_time : float
        The time at which end_timer() was called.
        
    start_lapse_time : float
        The time at which start_timer_lapse() was most recently called.
        
    end_lapse_tme : float
        The time at which end_timer_lapse() was most recently called.

    Methods
    -------
    timer():
        Constructor for the timer object, which starts the outer timer.
        
    start_lapse():
        Starts the inner timer.
        
    end_lapse():
        Ends the inner timer.
        
    end_timer():
        Ends the outer timer.
        
    get_total_time():
        Returns the time in seconds recorded on the outer timer, 
        after it has been stopped.
        
    get_lapsed_time():
        Returns the time in seconds most recently recorded on the inner timer, 
        after a given lapse has been stopped.
        
        
    Raises
    ------
    
    RuntimeError:
        Occurs if the lapse timer is ended before it is started, 
        or started while it is already running.
        
        Additionally occurs if the recorded time of the outer/inner timer is 
        requested while it is still running.
    """
    
    def __init__(self):
        self.is_running = False
        
    def start(self):
        if self.is_running:
            raise RuntimeError("The timer is already running.")
        self.start_time = time.time()
        self.is_running = True
        
        
    def stop(self):
        if not(self.is_running):
            raise RuntimeError("The timer has not been started yet.")
        self.end_time = time.time()
        self.is_running = False
        
    def get_lapsed_time(self):
        if self.is_running:
            raise RuntimeError("The timer is still running.")
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
0,1,1,0                       NUU,IPKT,NOYES,ITRANS
%(emin)s,%(emax)s
%(gampn_orbitals)s                 fermi_parityP, NORBITP, LEVELP
%(gampn_orbitals)s                  fermi_parityN, NORBITN, LEVELN
%(Z)s,%(A)s                                                Z,A
%(current_eps)s,%(current_gamma)s,0.00,0.0,0.0000,8,8,0,0
(LAST CARD: EPS,GAMMA,EPS4,EPS6,OMROT,NPROT,NNEUTR,NSHELP,NSHELN)
'''
        
    elif program == "asyrmo":
        
        template = '''
'%(current_f016)s' '%(current_f017)s' '%(current_f018)s' FILE16,FILE17,FILE18
1,0                                        IPKT,ISKIP
%(istrch)s,%(irec)s                                        ISTRCH,IREC
%(vmi)s,4,4,8,0.0188,100.00                      VMI,NMIN,NMAX,IARCUT,A00,STIFF
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
0,%(cutoff)s,1,0.75,-1                ISPEC,CUTOFF,IQ,GSFAC,GR
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
        var_list = ["mark_spin"]

    elif var_type == "int":
        var_list = ["A", "Z", "num_to_record", "num_orbs", "nu"]

    elif var_type == "float":
        var_list = ["gs_energy", "x1_energy", "x2_energy", "x3_energy", "gs_mu", "x1_mu"]
        
    else: raise ValueError("unrecognised type: " + var_type)
    
    return var_list


class PropertyData:
    
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
            
        else: raise ValueError("property not regonised: " + name)
        
        self.num = num
        self.prop = prop
        self.sort = sort
        
        # set graph title and axis label strings
        if prop == "energies" and not(sort=="Fermi"): 
            
            self.title = "Energies of "+ sort + num + " States"
            self.axis_label = self.title + " / keV"
           
        elif prop == "mag_moments":
                
            if sort == "Excited State ": self.title = "Magnetic Dipole Moment of Excited State " + num 
            elif sort == "Spin ": self.title = "Magnetic Dipole Moment of Spin " + num + " States"
            elif sort == "Ground": self.title = "Ground State Magnetic Dipole Moment"
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
        
        else: raise ValueError("property not recognised: " + name)
        
        
        
        