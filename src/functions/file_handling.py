#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 21 11:37 2025

@author: katyrr

Functions for reading/writing files.

"""

def read_file(path):
    """
    A function to open a file, then read and return its full contents. 

    Parameters
    ----------
    path : string
        The relative file path from the current working directory to the file being opened.

    Returns
    -------
    lines : list of strings
        The contents of the file, where each element of the list is one line of the file.

    """
    with open(path, 'r') as f:
        lines = f.readlines()
    
    return lines



def write_file(path, text):
    """
    A function to open a file, and write some text to it. If the file already 
    exists, it will be overwritten. If it does not exist yet, it will be created.

    Parameters
    ----------
    path : string
        The relative file path from the current working directory to the file being opened.
        Includes the name of the file itself at the end of the file path.

    text : string
        The text that will be written to the file. 

    Returns
    -------
    None.

    """
    
    with open(path, 'w') as f:
        f.write(text)

