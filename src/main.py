#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  7 10:54:23 2025

@author: Katy Robson




HOW TO USE:

- This code file is stored in /Code/src
- Modules "functions.py", "structs.py", and "graph_plotting.py" are also stored in this directory.

- How to use for the first time:

    1. Download the repository from GitHub: 
            git clone https://github.com/katyrr/Exploring-Triaxiality-in-the-Lead-Region
            cd Exploring-Triaxiality-in-the-Lead-Region
    
       If you have uv installed, you can then run:
            uv venv
            uv sync
       
       There is no need to activate the venv before running. Using uv is optional, but
       might reduce the risk of errors caused by a mismatch in python or framework versions.
    
    2. Check that the pre-compiled original PTRM codes in the /Code/src/ptrm/MacOS/MO/ 
       or /Code/src/ptrm/64bit/MO/ folder have execute permissions turned on. 
       
       On Mac you can fix this by navigating to the /Code/src/ptrm/MacOS/MO folder 
       in Terminal, and then using the commands "chmod 775 gampn", "chmod 775 asyrmo", 
       "chmod 775 probamo". Then use "ls -ltr" to see that the permissions (in the 
       leftmost column) now have 'x's (execute). You might also have to grant permission
       in settings (after attempting to run for the first time, the permission request 
       will appear in Privacy).
       
       On Windows you can fix this by navigating to the /Code/src/ptrm/64bit/MO 
       folder in PowerShell, and then using the commands "Unblock-File -Path GAMPN.exe", 
       "Unblock-File -Path ASYRMO.exe", "Unblock-File -Path PROBAMO.exe".
       
    3. Make a new folder in /Code, typically named as the nuclide you're calculating,
       e.g. /Code/Pt177.
    
    4. Copy the config file from /Code/Examples to your new folder. (All the other
       necessary folders and files will be created automatically when you run the
       code for the first time).
    
    5. Make any necessary changes to the settings in your config file in /Code/<folder>.

    6. Run the codes with command: "python src/main.py <folder>" 
                                or "Python3 src/main.py <folder>" 
                                or "uv run python src/main.py <folder>"
       
       Optional command line argument:

       "--display-figures" or "-d" to display figs on the screen as they are plotted
       (whether or not this option is used, the figures will be saved in the "figures" folder)


- How to use after the first time:

    1. Make changes to your config file in /Code/<folder>.
    
    2. Run the codes with command: "python src/main.py <folder>" 
                                or "Python3 src/main.py <folder>" 
                                or "uv run python src/main.py <folder>"
       
       Optional command line argument:

       "--display-figures" or "-d" to display figs on the screen as they are plotted
       (whether or not this option is used, the figures will be saved in the "figures" folder)
     
    

"""

import numpy as np                                  
import matplotlib.pyplot as plt    
import sys     

import functions.file_handling as fh 
import functions.read_config as rc     
import functions.run_ptrm as ptrm
import functions.read_gampn as rgam
import functions.read_asyrmo as rasy
import functions.read_probamo as rprob
import functions.structs as st
import functions.graph_plotting as gr
import functions.analyse_results as anyl

from classes.timer import Timer


def print_div():
    # Use to help organise console output into easy-to-read sections
    print("========================================================================================") 



def main():

    #%% 
    ''' 1. SET UP ---------------------------------------------------------------------------------

    - Create timers (one to time the whole program, and one to time small sections).
    - Read command line arguments and locate/create data subfolder.
    
    '''

    main_timer, sub_timer = Timer(), Timer()
    main_timer.start()

    print_div()
    data_subfolder_path = fh.locate_data_subfolder(sys.argv)
    
    
    #%%
    ''' 2. READ CONFIG FILE -----------------------------------------------------------------------

    - Create empty dictionaries for storing input settings and data.
    - Read the config file and save settings in dictionaries.
    - Generate the orbitals input for gampn (e.g. "+4 19 20 21 22").
    - Set figure resolution.
    - Create any missing directory subfolders.
    - Count the number of data points being calculated.
    - Print some reports to the console.

    ''' 

    code_settings, ptrm_inputs, data_points, experimental_data, graphs_to_plot = {}, {}, {}, {}, {}
    rc.read_config(data_subfolder_path, code_settings, ptrm_inputs, data_points, experimental_data, graphs_to_plot)

    ptrm_inputs["current_orbitals"] = ptrm.write_orbitals(ptrm_inputs["fermi_level"]//2, ptrm_inputs["num_orbs"], ptrm_inputs["par"])
    # (useful for debugging) hard coded versions of the above:
    # inputs["current_orbitals"] = "-15 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38" 
    # inputs["current_orbitals"] = fn.write_orbitals(28, inputs["num_orbs"], inputs["par"])

    plt.rcParams['figure.dpi'] = code_settings["figure_res"]  # set figure resolution
    if "-d" in sys.argv or "--display-figures" in sys.argv:
        print("DEBUG: display figs")
        code_settings["display_figures"] = True
    else:
        code_settings["display_figures"] = False

    fh.setup_directory(data_subfolder_path, code_settings["num_cores"], code_settings["OS"])
    
    code_settings["num_points"] = len(data_points["eps"])
    print("Number of data points = ", code_settings["num_points"])
    print("Deformation range:")
    print(f"\teps = [{data_points["eps"][0]:.3f}, {data_points["eps"][-1]:.3f}]")
    print(f"\tgamma = [{data_points["gamma_degrees"][0]:.1f}, {data_points["gamma_degrees"][-1]:.1f}] degrees")


    #%%   
    ''' 3. RUN GAMPN ------------------------------------------------------------------------------

    - Write the input .DAT files for the gampn code.
    - Calculate batch settings
    - Configure a script writer.
    - Run the batches. The .OUT files are generated in the 'outputs' directory folder.
    - Read the output files for the fermi level index, energy, and parity, so that new sets 
      of orbitals can be calculated dynamically.
    - Re-run gampn with the new set of orbitals, so that the strong-coupling basis 
      can be maximised (to 15 orbitals) when calculating matrix elements.
    - No need to re-read the outputs, because the properties we read earlier are 
      not affected, and the recalculated matrix elements will be passed to the next
      program automatically.

    '''

    data_points["file_tags"] = []
    ptrm.write_input_files(code_settings["num_points"], data_subfolder_path, "gampn", ptrm_inputs, data_points, first_run=True)
    
    batch_settings = ptrm.get_batch_settings(code_settings["num_cores"], code_settings["num_points"])
    run_program = ptrm.configure_script_writer(data_subfolder_path, code_settings["OS"], batch_settings, data_points["file_tags"])
    
    sub_timer.start()
    run_program("gampn")
    sub_timer.stop()

    print(f"\n***** Returned from gampn (first run) after {sub_timer.get_lapsed_time():.2f} seconds. *****\n")

    # set up arrays to store data 
    data_points["asyrmo_orbitals"] = []
    output_data = {"fermi_parities": [0]*code_settings["num_points"], 
                   "fermi_energies_hw": [0]*code_settings["num_points"], 
                   "fermi_energies_mev": [0]*code_settings["num_points"], 
                   "fermi_indices": [0]*code_settings["num_points"]}
    
    rgam.read_gampn(code_settings["num_points"], data_subfolder_path, data_points, output_data, ptrm_inputs)

    ptrm.write_input_files(code_settings["num_points"], data_subfolder_path, "gampn", ptrm_inputs, data_points)

    sub_timer.start()
    run_program("gampn")
    sub_timer.stop()
    print(f"***** Returned from gampn (second run) after {sub_timer.get_lapsed_time():.2f} seconds. *****\n")

    #%%
    ''' 4. RUN ASYRMO -----------------------------------------------------------------------------

    - Use the existing list of file tags to write a .DAT file for each data point.
    - Use the existing script writer to write and run asyrmo; dividing up the batches (as for gampn). 
    - Read the output files and check for the "NO DECOUPLING PARAMETERS CALCULATED" error.
    - Record the value of the DELTA parameter.
    - Check for the "SORRY I FOUND NO SOLUTIONS" error, and exclude those files from future analysis.

    '''

    ptrm.write_input_files(code_settings["num_points"], data_subfolder_path, "asyrmo", ptrm_inputs, data_points)
    
    sub_timer.start()
    run_program("asyrmo")
    sub_timer.stop()

    print(f"***** Returned from asyrmo after {sub_timer.get_lapsed_time():.2f} seconds. *****\n")

    output_data["delta"] = []

    rasy.read_asyrmo(data_points["file_tags"], data_subfolder_path, output_data)

    #%%
    ''' 5. RUN PROBAMO ----------------------------------------------------------------------------

    - Use the existing list of file tags to write a .DAT file for each data point.
    - Use the existing script writer to write and run probamo; dividing up the batches as for gampn. 
    - Read the probamo output files and store energy level data (with magnetic dipole and electric quadrupole moments)

    '''

    ptrm.write_input_files(code_settings["num_points"], data_subfolder_path, "probamo", ptrm_inputs, data_points)

    sub_timer.start()
    run_program("probamo")
    sub_timer.stop()

    print(f"***** Returned from probamo after {sub_timer.get_lapsed_time():.2f} seconds. *****\n")
    
    restructured_output_data = rprob.read_probamo(code_settings["num_points"], data_subfolder_path, data_points, output_data, experimental_data, ptrm_inputs, code_settings["print_details"])
    

    ''' 6. PLOT GRAPHS ----------------------------------------------------------------------------

    - Convert data into PropertyData class instances, for easy plotting later


    - Set a subtitle containing the values of E2PLUS and GSFAC input, if requested.
    - Set which graphs should be plotted (from config, or overwritten below). 
      Any not listed are False by default.

    - If deformation was input as a mesh, plot filled contours in polar coordinates.
    - If only one of eps/gamma/e2plus is varied, plot line graphs.
        
        
    '''

    data_to_plot = gr.prepare_data_to_plot(experimental_data, data_points["file_tags"], restructured_output_data)

    #!!! set graph subtitle:
    if code_settings["include_subtitle"]:
        subtitle = r'$E(2^+)$ = ' + str(ptrm_inputs["current_e2plus"]) + '; gsfac = ' + str(ptrm_inputs["gsfac"])
    else:
        subtitle = ''
        
        
    # set which graphs to plot:
    for i in graphs_to_plot:
        if i in data_to_plot:
            data_to_plot[i].plot = graphs_to_plot[i]
        else:
            print("property not recorded, check that it is included in experimental data inputs:\n\t", i)

    #!!! override graph plotting options (useful when running cell by cell):
        
    # data_to_plot["fermi_indices"].plot = 0
    # data_to_plot["delta"].plot = 0
    # data_to_plot["fermi_energies_mev"].plot = 0
    # data_to_plot["fermi_energies_hw"].plot = 0

    # data_to_plot["gs_mag_moments"].plot = 0
    # data_to_plot["gs_quad_moments"].plot = 0
    # data_to_plot["gs_spin_floats"].plot = 1

    # data_to_plot["spin_1/2_energies"].plot = 0
    # data_to_plot["spin_3/2_energies"].plot = 0
    # data_to_plot["spin_5/2_energies"].plot = 0
    # data_to_plot["spin_7/2_energies"].plot = 0
    # data_to_plot["spin_9/2_energies"].plot = 0
    # data_to_plot["spin_11/2_energies"].plot = 0
    # data_to_plot["spin_13/2_energies"].plot = 0

    # data_to_plot["spin_1/2_mag_moments"].plot = 0
    # data_to_plot["spin_3/2_mag_moments"].plot = 0

    # data_to_plot["rms"].plot = 0

    # data_to_plot["all_energies"].plot = 0
    # data_to_plot["shifted_energies"].plot = 0

    # data_to_plot["gap_9_13"].plot = 0

    # override settings
    # code_settings["mark_exp"] = 1
    # code_settings["mark_exp_tol"] = 0
    # code_settings["mark_points"] = 1
    # code_settings["mark_spin"] = 0


    sub_timer.start()

    data_points["agreed"] = [0]*len(data_points["eps"])
    num_comparisons = 0 

    # start plotting graphs:
    for i in data_to_plot:
        
        prop = data_to_plot[i]
        
        if not(prop.plot):
            continue
        
        ptrm_inputs["current_graph"] = prop.title # makes several later inputs more efficient
        print("plotting graph: %(current_graph)s" % ptrm_inputs) 
        
        if ptrm_inputs["deformation_input"] == "mesh":  
            
            gr.plot_mesh_graph(prop, data_points, code_settings, ptrm_inputs, data_to_plot["gs_spin_floats"], subtitle, data_subfolder_path)
            if np.isfinite(prop.experimental_data).all() and code_settings["mark_exp"]: 
                num_comparisons += 1
    

        elif (ptrm_inputs["deformation_input"] ==  "gamma" 
            or ptrm_inputs["deformation_input"] == "eps"
            or len(data_points["e2plus"]) > 1):
            
            gr.plot_line_graph(prop, ptrm_inputs, data_points, code_settings, experimental_data, subtitle, data_to_plot["gs_spin_floats"], data_subfolder_path)
            
            
        

    #%%
    ''' 7. ASSESS AGREEMENT OF CALCULATIONS WITH EXPERIMENT ---------------------------------------

    - Print information about the best agreement and its location.
    - Plot a graph to show data point agreement across all data points.
    - Print the mean energies of each level and the mean gs moments, with standard error.
    - Print the total runtime.

    '''
        
    anyl.check_agreement(code_settings["print_details"], data_points, num_comparisons)
    anyl.plot_agreement(data_points, num_comparisons, code_settings, ptrm_inputs, data_to_plot["gs_spin_floats"], subtitle, data_subfolder_path)
    
    print_div()
    print("\n******** mean and standard error in the mean *********")
    report_means=["spin_1/2_energies", "spin_3/2_energies", "spin_5/2_energies", "spin_7/2_energies", "spin_9/2_energies", "spin_11/2_energies", "spin_13/2_energies",
                  "gs_mag_moments", "spin_1/2_mag_moments", "spin_3/2_mag_moments", "spin_5/2_mag_moments", "spin_7/2_mag_moments", "spin_9/2_mag_moments", "spin_11/2_mag_moments", "spin_13/2_mag_moments",
                  "gs_quad_moments", "spin_1/2_quad_moments", "spin_3/2_quad_moments", "spin_5/2_quad_moments", "spin_7/2_quad_moments", "spin_9/2_quad_moments", "spin_11/2_quad_moments", "spin_13/2_quad_moments"]
    for i in report_means:
        anyl.report_mean(data_to_plot[i], code_settings["print_details"])

    sub_timer.stop()
    main_timer.stop()
    print("\n****************************************************************************************")
    print("Finished plotting graphs in time = %.2f seconds" % (sub_timer.get_lapsed_time()))
    print("total runtime = %.2f seconds" % (main_timer.get_lapsed_time()))
    print("****************************************************************************************\n")


if __name__ == "__main__":
    main()