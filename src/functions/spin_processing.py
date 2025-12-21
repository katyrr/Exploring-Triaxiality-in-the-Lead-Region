#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 21 11:36 2025

@author: katyrr

Functions for processing spin data.

"""


def spin_string_to_float(spin_string):
    """
    A function to covert the fractional (string) representation of a nuclear 
    spin to its float representation.

    Parameters
    ----------
    spin_string : string
        A fractional represention of a half-int spin.
        e.g. "1/2", "3/2", etc.

    Raises
    ------
    ValueError
        Occurs if the spin_string is not input in the correct format "n/2" where n is an integer.

    Returns
    -------
    spin_float : float
        The float representation of the spin. 
        e.g. 0.5, 1.5, etc.

    """
    if spin_string[1] == "/":
        spin_float = float(spin_string[0])/2
    elif spin_string[2] == "/":
        spin_float = float(spin_string[0:2])/2 
    else:
        raise ValueError("Cannot parse spin. \nCheck it has been input in " +
                         "the format 'n/2' where n is an int.")
    return spin_float


def spin_float_to_string(spin_float):
    """
    A function to convert a float representation of a nuclear spin 
    to its fractional representation, as a string.

    Parameters
    ----------
    spin_float : float
        The float representation of the spin. 
        e.g. 0.5, 1.5, etc.

    Raises
    ------
    ValueError
        Occurs if the ipnut spin_float cannot be converted to a half-integer fraction.
        i.e. must be n.5 where n is an integer.

    Returns
    -------
    spin_string : string
        A fractional represention of a half-int spin.
        e.g. "1/2", "3/2", etc.

    """
    numerator = float(spin_float)*2
    if str(numerator)[-1] != "0":
        raise ValueError("Cannot convert float to fraction. \nCheck it has " + 
                         "been input as a half-integer float.")
    else:
        spin_string = str(int(numerator))+"/2"
        return spin_string
    