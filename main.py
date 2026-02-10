#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  7 10:54:23 2025

@author: Katy Robson


To-do:
- clean up docs 
- write unit tests for individual funcs
- make a constants file for magic numbers
- add subtitle option to CLAs


================================= HOW TO USE FOR THE FIRST TIME: =================================

1. Download the repository from GitHub: 
        git clone https://github.com/katyrr/Exploring-Triaxiality-in-the-Lead-Region
        cd Exploring-Triaxiality-in-the-Lead-Region

    Optional:
    If you have uv installed, you can then run:
        uv venv
        uv sync
    
    There is no need to "activate" the venv before running. 
    Using uv is optional, but might reduce the risk of version errors.

2. Check that your current working directory is /Code. 
    If you run 'ls' you should see 'main.py' in the list of files.

3. Do a test run of the code:
        uv run main.py example_folder

    or (without uv):
        python main.py example_folder
    
    or (depending on your OS and python version): 
        Python3 main.py example_folder

    (a) It will automatically generate a folder called 'example_folder' in the 
        directory /Code/data where all the input files, output files, binaries, figures, 
        and helper scripts will be stored. You can give this folder any name (typically it
        might be the name of the nucleus you're studying, e.g. Pt177).

    (b) It will automatically generate a config.txt file in example_folder, from 
        the template at /Code/static/config_template.txt.

        You may wish to edit the template beforehand, to configure settings such as your
        operating system, how many CPU cores you want to utilise, your desired figure 
        resolution and formatting, which graphs you want plotted by default, and which nucleus 
        to calculate for. 
        
        The default is 177Pt at (ε, γ) = (0.26, 24º), on MacOS with 8 cores.

        The template will be used whenever you run the code in a new data folder for the 
        first time, so you may find it convenient to adjust your default settings now!

    (c) It will attempt to enable 'execute' permissions on the ptrm codes, which are stored
        in /Code/src/ptrm/MacOS/MO or /Code/src/ptrm/64bit/MO. Depending on your operating 
        system, it may or may not be successful. If it doesn't work, it will throw an error, 
        and you will have to enable these permissions manually:
    
        On Mac: navigate to the /Code/src/ptrm/MacOS/MO folder in Terminal, and then 
        use the commands "chmod 775 gampn", "chmod 775 asyrmo", "chmod 775 probamo". 
        Then use "ls -ltr" to see that the permissions (in the leftmost column) now 
        have 'x's (execute). You might also have to grant permission in Settings 
        (after attempting to run for the first time, the permission request will appear in Privacy).
    
        On Windows: navigate to the /Code/src/ptrm/64bit/MO folder in PowerShell, and then 
        use the commands "Unblock-File -Path GAMPN.exe", "Unblock-File -Path ASYRMO.exe", 
        "Unblock-File -Path PROBAMO.exe".

    (d) It will print progress and key results to the command line, and you can find more
        detailed results in the example_folder/outputs, and figures in 
        example_folder/figures.

        Optional:
        You can reroute the stdout to a text file with:
            touch path/to/saved_stdout.txt                             # creates the file
            uv run main.py example_folder > path/to/saved_stdout.txt   # sends stdout to the file



================================= HOW TO USE AFTER THE FIRST TIME: =================================

1. (a) If you just want to make some small adjustments to the config, you can edit the (existing)
        config file in /Code/data/example_folder, and move to step 2.
    
    (b) If you want to start a new folder (e.g. to calculate for a new nucleus), you can either 
        manually create a new folder in /Code/data, copy a config file from the template (or elsewhere), 
        and make the desired changes before moving to step 2. 

        Or you can skip straight to step 2 and let the code automatically generate a new folder 
        and config from the template.
    
    (i.e. if example_folder already exists with a config file, it will use that. 
    If it doesn't already exist, it will make one.)

2. Run the codes with command: "uv run python src/main.py example_folder"
                                or "python src/main.py example_folder" 
                                or "Python3 src/main.py example_folder" 
    
    Optional command line arguments:

    "--display-figures" or "-d":
        display figures on the screen as they are plotted
        (whether or not this option is used, the figures will be saved in the "figures" folder)

    "--verbose" or "-v":
        print more detailed information in stdout
     
===================================================================================================

"""

import numpy as np                                  
import matplotlib.pyplot as plt    
import sys     
import os
import datetime

import src.functions.file_handling as fh 
import src.functions.read_config as rc     
import src.functions.run_ptrm as ptrm
import src.functions.read_gampn as rgam
import src.functions.read_asyrmo as rasy
import src.functions.read_probamo as rprob
import src.functions.graph_plotting as gr
import src.functions.analyse_results as anyl

from src.classes.timer import Timer
from src.functions.parse_args import args



def print_div():
    # For organising console output into easy-to-read sections
    print("========================================================================================") 



def main():

    #%% 
    ''' 1. SET UP ---------------------------------------------------------------------------------

    - Create two timers (one to time the whole program, and one to time small sections).
    - Locate/create data subfolder, where inputs/outputs will be stored, and the ptrm will be run.
    - Read config file and sort contents into dictionaries.
    - Set figure resolution and display settings.
    - Ensure the data subfolder has the required directories and structure.
    
    '''
    timestamp = '{:%Y-%m-%d %H;%M;%S}'.format(datetime.datetime.now())
    main_timer, sub_timer = Timer(), Timer()
    main_timer.start()

    print_div()
    print(f"STARTING at {timestamp}")
    data_subfolder_path = fh.locate_data_subfolder()
    
    code_settings, ptrm_inputs, data_points, experimental_data, graphs_to_plot = {}, {}, {}, {}, {}
    rc.read_config(data_subfolder_path, code_settings, ptrm_inputs, 
                   data_points, experimental_data, graphs_to_plot)
    
    plt.rcParams['figure.dpi'] = code_settings["figure_res"]
    plt.rcParams.update({'figure.autolayout': True})
    code_settings["display_figures"] = args.display_figures

    fh.setup_directory(data_subfolder_path, code_settings["num_cores"], code_settings["OS"])
    
    print_div()

    # (useful for debugging) hard coded versions gampn input orbitals:
    # inputs["current_orbitals"] = "-15 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38" 
    # inputs["current_orbitals"] = fn.write_orbitals(28, inputs["num_orbs"], inputs["par"])


    #%%   
    ''' 2. RUN GAMPN ------------------------------------------------------------------------------

    - Write the input .DAT files for the gampn code.
    - Calculate batch settings and configure a script writer.
    - Run the batches. The .OUT files are generated in the 'outputs' directory folder.
    - Read the output files and store results in arrays within dictionaries. Data includes info 
      about the fermi level, and dynamically calculated orbital sets for the strong-coupling basis.
    - Re-run gampn with the new set of orbitals. There is no need to re-read the outputs, 
      because the properties we read earlier are not affected, and the recalculated matrix 
      elements will be passed to the next program automatically.

    '''

    data_points["file_tags"] = []
    ptrm.write_input_files(code_settings["num_points"], data_subfolder_path, "gampn", ptrm_inputs, data_points, first_run=True)
    
    batch_settings = ptrm.get_batch_settings(code_settings["num_cores"], code_settings["num_points"])
    run_program = ptrm.configure_script_writer(data_subfolder_path, code_settings["OS"], batch_settings, data_points["file_tags"])
    
    sub_timer.start()
    run_program("gampn")
    sub_timer.stop()
    print(f"***** Returned from gampn (first run) after {sub_timer.get_lapsed_time():.2f} seconds. *****\n")

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
    ''' 3. RUN ASYRMO -----------------------------------------------------------------------------

    - Write the input .DAT files for the asyrmo code.
    - Use the existing script writer to write and run asyrmo. 
    - Read the output files, record the DELTA parameter, and check for errors: 
        - "NO DECOUPLING PARAMETERS CALCULATED" error (fatal)
        - "SORRY I FOUND NO SOLUTIONS" error (exclude those files from further analysis)

    '''

    ptrm.write_input_files(code_settings["num_points"], data_subfolder_path, "asyrmo", ptrm_inputs, data_points)
    
    sub_timer.start()
    run_program("asyrmo")
    sub_timer.stop()
    print(f"***** Returned from asyrmo after {sub_timer.get_lapsed_time():.2f} seconds. *****\n")

    output_data["delta"] = []
    rasy.read_asyrmo(data_points["file_tags"], data_subfolder_path, output_data)



    #%%
    ''' 4. RUN PROBAMO ----------------------------------------------------------------------------

    - Write the input .DAT files for the probamo code.
    - Use the existing script writer to write and run probamo. 
    - Read the output files, and record energy level data, with magnetic dipole and electric 
      quadrupole moments. 
    - Restructure output data from dict of arrays to array of dicts, with bad data points masked, 
      and gaps filled with np.NaN.

    '''

    ptrm.write_input_files(code_settings["num_points"], data_subfolder_path, "probamo", ptrm_inputs, data_points)

    sub_timer.start()
    run_program("probamo")
    sub_timer.stop()
    print(f"***** Returned from probamo after {sub_timer.get_lapsed_time():.2f} seconds. *****\n")
    
    print_div()
    rprob.read_probamo(code_settings["num_points"], data_subfolder_path, data_points, output_data, ptrm_inputs)
    restructured_output_data = rprob.process_data(output_data, data_points, experimental_data, ptrm_inputs["ispin"])



    ''' 5. PLOT GRAPHS ----------------------------------------------------------------------------

    - Convert data into PropertyData class instances, for easy plotting later
    - Set a subtitle containing the values of E2PLUS and GSFAC input, if requested.
    - Set which graphs should be plotted (from config, or overwritten below). 
      Any not listed are False by default.

    - If deformation was input as a mesh, plot filled contours in polar coordinates.
    - If only one of eps/gamma/e2plus is varied, plot line graphs.
        
    '''

    data_to_plot = gr.prepare_data_to_plot(experimental_data, data_points["file_tags"], restructured_output_data)

    #!!! set graph subtitle:
    subtitle = ''
    if code_settings["include_subtitle"]:
        subtitle = r'$E(2^+)$ = ' + str(ptrm_inputs["current_e2plus"]) + '; gsfac = ' + str(ptrm_inputs["gsfac"])
        
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

    #!!! override settings
    # code_settings["mark_exp"] = 1
    # code_settings["mark_exp_tol"] = 0
    # code_settings["mark_points"] = 1
    # code_settings["mark_spin"] = 0


    sub_timer.start()

    data_points["agreed"] = [0]*len(data_points["eps"])
    num_comparisons = 0 

    print_div()
    
    # start plotting graphs:
    for i in data_to_plot:
        
        prop = data_to_plot[i]
        
        if not(prop.plot):
            continue
        
        ptrm_inputs["current_graph"] = prop.title # makes several later inputs more efficient
        print("plotting graph: %(current_graph)s" % ptrm_inputs) 
        
        file_name = f"{timestamp} {ptrm_inputs["current_graph"].replace("/", "_")}"
        fig_path = os.path.join(data_subfolder_path, "figures", file_name)

        if ptrm_inputs["deformation_input"] == "mesh":  
            
            gr.plot_mesh_graph(prop, data_points, code_settings, ptrm_inputs, data_to_plot["gs_spin_floats"], subtitle, fig_path)
            if np.isfinite(prop.experimental_data).all() and code_settings["mark_exp"]: 
                num_comparisons += 1
    

        elif (ptrm_inputs["deformation_input"] ==  "gamma" 
            or ptrm_inputs["deformation_input"] == "eps"
            or len(data_points["e2plus"]) > 1):
            
            gr.plot_line_graph(prop, ptrm_inputs, data_points, code_settings, experimental_data, subtitle, data_to_plot["gs_spin_floats"], fig_path)
            
    print_div()
    
    #%%
    ''' 6. ASSESS AGREEMENT OF CALCULATIONS WITH EXPERIMENT ---------------------------------------

    - Print information about the best agreement and its location.
    - Plot a graph to show data point agreement across all data points.
    - Print the mean energies of each level and the mean gs moments, with standard error.
    - Print the total runtime.

    '''
    
    anyl.check_agreement(data_points, num_comparisons)
    
    plot_agreement = False
    if plot_agreement:
        anyl.plot_agreement(data_points, num_comparisons, code_settings, ptrm_inputs, data_to_plot["gs_spin_floats"], subtitle, data_subfolder_path)
    print_div()

    print("\nMean and standard error in the mean:")
    if not args.verbose:
        print("[only printing lowest energy state of each spin; for all states use -v or --verbose]\n")
    report_means=["spin_1/2_energies", "spin_3/2_energies", "spin_5/2_energies", "spin_7/2_energies", "spin_9/2_energies", "spin_11/2_energies", "spin_13/2_energies",
                  "gs_mag_moments", "spin_1/2_mag_moments", "spin_3/2_mag_moments", "spin_5/2_mag_moments", "spin_7/2_mag_moments", "spin_9/2_mag_moments", "spin_11/2_mag_moments", "spin_13/2_mag_moments",
                  "gs_quad_moments", "spin_1/2_quad_moments", "spin_3/2_quad_moments", "spin_5/2_quad_moments", "spin_7/2_quad_moments", "spin_9/2_quad_moments", "spin_11/2_quad_moments", "spin_13/2_quad_moments"]
    for i in report_means:
        anyl.report_mean(data_to_plot[i])

    sub_timer.stop()
    main_timer.stop()
    print_div()
    print("Finished plotting graphs in time = %.2f seconds" % (sub_timer.get_lapsed_time()))
    print("Total runtime = %.2f seconds" % (main_timer.get_lapsed_time()))
    print_div()

    return


if __name__ == "__main__":
    main()