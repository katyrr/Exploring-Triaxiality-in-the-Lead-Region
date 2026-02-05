#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 22 10:24 2025

@author: katyrr

A timer class.

"""

import time                                                                  

class Timer():
    """
    A stopwatch for timing how long it takes to execute a portion of code.
    In the event of unexpected behaviour, the program does not crash, 
    but warnings are printed to console and timer data may not be available.

    Attributes
    ----------
    is_running : bool
        A record of whether the timer is currently running.
        
    start_time : float
        The time at which the .start() method was most recently called.
        
    end_time : float
        The time at which end_timer() was most recently called.
        
    Methods
    -------
    Timer():
        Constructor for the timer object.
        
    start():
        Starts the timer.
        Fails (prints warning) when is_running==True.
        
    stop():
        Stops the timer.
        Fails (prints warning) when is_running==False.
        
    time_elapsed = get_lapsed_time():
        Returns the time interval (in seconds) most recently recorded on the timer.
        Fails (returns None and prints warning) when is_running==True.
        Fails (returns None and prints warning) start_time or end_time are None.
    
    """
    
    def __init__(self):
        self.is_running = False
        self.start_time = None
        self.end_time = None
        
    def start(self):
        if self.is_running:
            print("WARNING: Could not start timer (already running).")
            return
        self.start_time = time.time()
        self.is_running = True
        
        
    def stop(self):
        if not(self.is_running):
            print("WARNING: Could not stop timer (not started yet).")
            return
        self.end_time = time.time()
        self.is_running = False
        
    def get_lapsed_time(self):
        if self.is_running:
            print("WARNING: Could not get lapsed time (the timer is still running).")
            return None
        if self.start_time is None or self.end_time is None:
            print("WARNING: Could not get lapsed time (recorded start and/or end time are None).")
            return None
        
        time_elapsed = self.end_time - self.start_time
        return time_elapsed
    
