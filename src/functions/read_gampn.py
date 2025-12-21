#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 21 11:47 2025

@author: katyrr

Functions for reading GAMPN.OUT and processing contents.

"""

def get_efac(lines):
     """ 
     A function which reads the full contents of the GAMPN.OUT file, 
     and returns the value of EFAC.
     
     Parameters
     ----------
     lines : list of strings
         The full contents of the GAMPN.OUT file.
         Each element of the list is a line read from the file.
    
     Returns:
     -------
     efac : float
         The value of the conversion factor from energy units hw to MeV.
     
     """
     ref = "     KAPPA    MY     EPS   GAMMA    EPS4     EPS6     W0/W00   NMAX  COUPL     OMROT      EFAC      QFAC\n"
     efac_line = lines.index(ref)
     efac = float(lines[efac_line+1][85:95].strip())
     return efac
 
    
def get_sp_level(lines, sp_index, parity):
    """ 
    A function which reads the full contents of the GAMPN.OUT file, 
    and returns the half-line containing data about the single particle level requested.
    
    Only half the line is required because the data for all 80 calculated levels 
    is output in two columns of 40 lines each.
    
    Parameters
    ----------
    lines : list of strings
        The full contents of the GAMPN.OUT file.
        Each element of the list is a line read from the file.
    
    index : int
        The index of the single particle level
        
    parity : string
        Defines how the single particle levels are indexed. 
        If '0', all the levels are counted (upwards from 1), regardless of parity.
        If '-', only the negative parity levels are counted (and positive parity levels are ignored.)
        If '+', only the positive parity levels are counted (and negative parity levels are ignored.)
        
    
    Returns:
    -------
    half_line : string
        A string containing data about the single particle level.
    
    """
    
    ref = "   #   ENERGY +/-(#)    <Q20>    <Q22>     <R2>     <JZ>      #   ENERGY +/-(#)    <Q20>    <Q22>     <R2>     <JZ>\n"
    levels_header_line = lines.index(ref)
    
    if parity == '0':
        # calculate the line number of the single particle level in the GAMPN.OUT file (indexed from zero!)
        sp_line = sp_index+levels_header_line+1  
        
        if sp_index > 40:
            # get only the second half of the line
            sp_line -= 40
            whole_line = lines[sp_line]
            half_line = whole_line[60:-1].strip()                                 
            
        else:
            # get only the first half of the line
            whole_line = lines[sp_line]
            half_line = whole_line[0:60].strip()   
            
    elif parity == '-' or parity == '+':
        
        reduced_lines = lines[levels_header_line+2 : levels_header_line+42]
        search = parity + '(#' + str(sp_index) + ')'
        
        found = False
        for l in reduced_lines:
            if search in l:
                whole_line = l
                found = True
                break
        
        if found:
            position = whole_line.index(search)
        else: raise RuntimeError("The requested single particle orbital was not calculated: " + search)
        
        if position > 60:
            half_line = whole_line[60:-1].strip()   
        else:
            half_line = whole_line[0:60].strip()   
        
    else: raise ValueError("Unrecognised parity. Allowed inputs are only '+', '-', or '0'.")
        
    
    return half_line
 
    
def get_sp_index(line, hash_index):
    """ 
    A function which takes a line of data about a single particle level, 
    and returns its orbital index.
    
    Parameters
    ----------
    line : string
        A string containing data about the single particle level.
    
    hash_index : int
        The index in the string that locates the "#" character.
        (Used as a reference point in the line).
        
    Returns:
    -------
    index : int
        The orbital index of the single particle level. 
        (numbered separately for positive and negative parity orbitals).
    
    """
    index_string = line[hash_index+1 : hash_index+3]
    if index_string[1] == ")":                                                  
        # in case the index is only a single digit, 
        # ignore the ")" that will have been caught:
        index_string = index_string[0]
        
    index = int(index_string)
    
    return index

def get_info(line):
    '''

    Parameters
    ----------
    line : string
        The line of text read from the GAMPN.OUT file which contains data 
        about a single particle (i.e. Nillson) level, including its energy,
        its index, its parity, etc.

    Returns
    -------
    parity : string
        The parity of the single particle level. Either '+' or '-'.
        
    energy_hw : float
        The energy of the single particle level, in oscillator (hbar omega) units.
        
    index : int
        The index of the single particle level, numbered separately for 
        positive and negative parity orbitals.

    '''
    hash_index = line.index("#")
    parity = line[hash_index-2]                            
    energy_hw = float(line[hash_index-10 : hash_index-4])  
    # energy_mev = energy_hw*inputs["efac"] # can't do this conversion without passing efac, and not really necessary in all cases anyway
    index = get_sp_index(line, hash_index)
    
    return (parity, energy_hw, index)


