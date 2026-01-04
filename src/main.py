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


- How to use after the first time:

    1. Make changes to your config file in /Code/<folder>.
    
    2. Run the codes with command: "python src/main.py <folder>" 
                                or "Python3 src/main.py <folder>" 
                                or "uv run python src/main.py <folder>"
     
    

"""

import numpy as np                                   
import math                                        
import matplotlib.pyplot as plt   
import os   
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

from functions.spin_processing import spin_string_to_float

from classes.timer import Timer





def main():

    #%% 
    """ 1. SET UP ---------------------------------------------------------------------------------

    - Create timers (one to time the whole program, and one to time small sections).
    - Read command line arguments and locate/create data subfolder.
    
    """

    main_timer, sub_timer = Timer(), Timer()
    main_timer.start()

    print("**********************************************************") # helps find the start of the calculation in the console output
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
    fh.setup_directory(data_subfolder_path, code_settings["num_cores"], code_settings["OS"])
    
    code_settings["num_points"] = len(data_points["eps"])
    print("Number of data points = ", code_settings["num_points"])
    print("Deformation range:")
    print(f"\teps = [{data_points["eps"][0]:.3f}, {data_points["eps"][-1]:.3f}]")
    print(f"\tgamma = [{data_points["gamma_degrees"][0]:.1f}, {data_points["gamma_degrees"][-1]:.1f}] degrees")


    #%%   
    ''' 3. RUN GAMPN ------------------------------------------------------------------

    - Write the input .DAT files for the gampn code.
    - Calculate batch settings
    - Configure a script writer.
    - Run the batches. The .OUT files are generated in the 'outputs' directory folder.

    '''

    data_points["file_tags"] = []
    ptrm.write_input_files(code_settings["num_points"], data_subfolder_path, "gampn", ptrm_inputs, data_points, first_run=True)
    
    batch_settings = ptrm.get_batch_settings(code_settings["num_cores"], code_settings["num_points"])
    run_program = ptrm.configure_script_writer(data_subfolder_path, code_settings["OS"], batch_settings, data_points["file_tags"])
    
    sub_timer.start()
    run_program("gampn")
    sub_timer.stop()

    print(f"***** Started running gampn, returned after {sub_timer.get_lapsed_time():.2f} seconds. *****\n")


    #%%
    ''' 4. READ GAMPN.OUT FILE 

    - For each data point (i.e. each GAMPN.OUT file):
        - Read the value of EFAC (the conversion factor from hw to eV).
        - Determine the line number of the fermi level orbital in the GAMPN.OUT file.
        - Read that line from GAMPN.OUT and get the energy, parity, and level index.
        - Dynamically locate the orbitals nearest to the fermi level in energy, for 
        future input.
        
    - Re-run gampn for all data points, using the orbitals located in the step above. 
        
    '''

    # set up arrays to store data 
    output_data = {"fermi_parities": [0]*code_settings["num_points"], "fermi_energies_hw": [0]*code_settings["num_points"], 
                "fermi_energies_mev": [0]*code_settings["num_points"], "fermi_indices": [0]*code_settings["num_points"]}
    data_points["asyrmo_orbitals"] = []

    for i in range(code_settings["num_points"]):

        output_file_path = os.path.join(data_subfolder_path, "Outputs", f"GAM_{data_points["file_tags"][i]}.OUT")
        lines = fh.read_file(output_file_path)
        
        ptrm_inputs["efac"] = rgam.get_efac(lines)
        _fermi_level_line = rgam.get_sp_level(lines, ptrm_inputs["fermi_level"], '0')
        _f_parity, _f_energy_hw, _f_index = rgam.get_info(_fermi_level_line)
        
        output_data["fermi_parities"][i] = _f_parity                             
        output_data["fermi_energies_hw"][i] = _f_energy_hw
        output_data["fermi_energies_mev"][i] = _f_energy_hw * ptrm_inputs["efac"]
        output_data["fermi_indices"][i] = _f_index
        
        # dynamically finds the orbitals nearest to the fermi level in energy:
        data_points["asyrmo_orbitals"].append(ptrm.find_orbitals(_f_index, ptrm_inputs["nu"], 
                                    ptrm_inputs["par"], _f_energy_hw, _f_parity, lines))


    #%% 
    '''5. RE-RUN GAMPN

    - Re-run gampn with the new set of orbitals, so that the strong-coupling basis 
    can be maximised (to 15 orbitals) when calculating matrix elements.
    - No need to re-read the outputs, because the properties we read earlier are 
    not affected, and the recalculated matrix elements will be passed to the next
    program automatically.


    '''

    ptrm.write_input_files(code_settings["num_points"], data_subfolder_path, "gampn", ptrm_inputs, data_points)

    sub_timer.start()
    run_program("gampn")
    sub_timer.stop()
    print(f"***** Started running gampn (again), returned after {sub_timer.get_lapsed_time():.2f} seconds. *****\n")

    _output_data = output_data # save a copy of the original before it's overwritten (useful when running cell by cell)



    #%%
    ''' 6. RUN ASYRMO 

    - Use the existing list of file tags to write a .DAT file for each data point.
    - Use the existing script writer to write and run asyrmo; dividing up the batches (as for gampn). 

    '''

    ptrm.write_input_files(code_settings["num_points"], data_subfolder_path, "asyrmo", ptrm_inputs, data_points)
    
    sub_timer.start()
    run_program("asyrmo")
    sub_timer.stop()

    print(f"***** Started running asyrmo, returned after {sub_timer.get_lapsed_time():.2f} seconds. *****\n")

    #%%

    ''' 8. READ ASYRMO 

    - Read the output files and check for the "NO DECOUPLING PARAMETERS CALCULATED" error.
    - Record the value of the DELTA parameter.
    - Check for the "SORRY I FOUND NO SOLUTIONS" error, and exclude those files from future analysis.

    '''

    output_data["delta"] = []

    for i in data_points["file_tags"]:

        output_file_path = os.path.join(data_subfolder_path, "Outputs", f"ASY_{i}.OUT")
        lines = fh.read_file(output_file_path)
        
        if not "PARTICLE-ROTOR  MODEL" in lines[0]: # then something has gone wrong
            raise RuntimeError("File " + i + " raised error in ASYRMO output: \n" + lines[0] )
            
        output_data["delta"].append(rasy.get_delta(lines)) # this also checks for the "SORRY I FOUND NO SOLUTIONS" error.
            

    #%%
    ''' 9. WRITE AND RUN PROBAMO 

    - Use the existing list of file tags to write a .DAT file for each data point.
    - Use the existing script writer to write and run probamo; dividing up the batches as for gampn. 

    '''

    ptrm.write_input_files(code_settings["num_points"], data_subfolder_path, "probamo", ptrm_inputs, data_points)

    sub_timer.start()
    run_program("probamo")
    sub_timer.stop()

    print(f"***** Started running probamo, returned after {sub_timer.get_lapsed_time():.2f} seconds. *****\n")


    #%%

    ''' 10. READ PROBAMO.OUT FILE

    For each file:
        
    - Read each line:
        - If it is a static transition, read the spin, energy, and magnetic moment.
        - Else ignore this line and move to the next.
        
    - Sort the file data into categories:
        - Group lines by spin (e.g. spin 1/2 energies, spin 1/2 magnetic dipole moments, etc).
        - Additionally (separately) record the expected ground state (the lowest state with the 
        same spin as the experimental gs) properties as a group.
        - Fill missing gaps with NaN values, such that the same set of properties has 
        been recorded for every data point (and if e.g. one data point found three 
        spin 1/2 states, then all data points should have a list of three spin 1/2 states, 
        even if some of them are NaN).
        - Restructure the data set and save separately (now each property is recorded
        as a list of values for all data points, rather than each data point having
        a list of properties associated with it). 

    - Calculate energy gaps between levels specified in the experimental section of the config input.

    - Ensure all data sets have the same size and shape.
    - Mask bad data points with reference to the DELTA data 
    (any DETLA = NaN values are bad, caused by some kind of convergence issue with BCS pairing).

    '''

    data_points["property_data"] = []

    for i in range(code_settings["num_points"]):
        
        output_file_path = os.path.join(data_subfolder_path, "Outputs", f"PROB_{data_points["file_tags"][i]}.OUT")
        lines = fh.read_file(output_file_path)

        _file_data = {}
        for _line in lines:
            
            _line_data = rprob.read_data(_line)  # get the spin, energy, and magnetic moment from this line if it is a static moment, else return False
            if not(_line_data):
                continue # to next line in file
            
            # sort line_data into file_data according to its spin
            _file_data = rprob.sort_by_spin(_line_data, _file_data)
            # additionally save data that corresponds to the expected experimental ground state
            _file_data = rprob.sort_by_expectation(_line_data, _file_data, ptrm_inputs)
        
        _file_data = rprob.missing_data(_file_data, ptrm_inputs)
        data_points["property_data"].append(_file_data)

    output_data = _output_data | rprob.restructure_data(data_points["property_data"], ptrm_inputs["ispin"], code_settings["print_details"])

    # get energy gap between third 9/2 and first 13/2 states
    for i in experimental_data:
        if not "engap_" in i:
            continue
        
        _spin1, _idx1, _spin2, _idx2 = rprob.parse_engap_input(i)
        
        output_data[i] = rprob.find_gaps(output_data["spin_"+_spin1+"/2_energies"], _idx1, output_data["spin_"+_spin2+"/2_energies"], _idx2, experimental_data[i])

    # output_data["engap_9.3_13.1"] = fn.find_gaps(output_data["spin_9/2_energies"], 3, output_data["spin_13/2_energies"], 1, 20) #!!!


    # ensure all data sets have the same size and shape, and mask bad points

    _mask = np.array([0 if np.isnan(x) else 1 for x in output_data["delta"]])

    for i in output_data:
        if isinstance(output_data[i][0], list):
            output_data[i] = rprob.fill_gaps(output_data[i])
            _list_mask = np.transpose(np.tile(_mask, (np.size(output_data[i][0]),1)))
            
            output_data[i] = np.where(_list_mask == 0, np.NaN, output_data[i])
        
        else: 
            output_data[i] = np.where(_mask == 0, np.NaN, output_data[i])

    _output_data_dict = output_data # save a copy of the original before it's overwritten, so that the code can be run cell-by-cell without errors.


    #%%
    ''' 11. PREPARE TO PLOT GRAPHS 

    - Record each data set in an instance of class PropertyData.
    - Calculate graph plotting attributes and store within the class.

    - Raise a ValueError if the property isn't recognised 
    (i.e. if more data sets are recorded in the future, they cannot be plotted without
    first hard-coding the calculation of things like axis labels, contour levels,
    colour bar ticks, etc).
    
    - Create a new data set containing all energies (of all spins) to plot together.
    - Create a new data set with all energies shifted to be relative to the expected 
    ground state (not necessarily the same as the calculated ground state at all points).
    This makes the output lines look smoother (no sharp bends when the ground state changes).
    - Create a new data set containing root mean squared error (i.e. discrepancy) between 
    the calculated lowest energy states of each spin and the exeperimental values (where available).

    '''

    # convert output_data from a dictionary of lists to a dictionary of PropertyData objects 
    output_data = {}
    for i in _output_data_dict:
        # print(i)
        output_data[i] = st.PropertyData(_output_data_dict[i], i)
        
        # calculate contour levels, colour bar ticks and labels, 
        # and assign experimental values and error tolerance if available.
        
        output_data[i] = gr.calculate_format_data(output_data[i], i, experimental_data)
        

    output_data["all_energies"] = rprob.collate_energy_data(output_data, len(data_points["file_tags"]), 
                                                        experimental_data["gs_spin_string"], experimental_data)

    # recalculate all energies relative to the spin entered into fn.collate_energy_data() above
    output_data["shifted_energies"] = rprob.shift_energy_levels(output_data["all_energies"]) 

    output_data["rms"] = rprob.calc_rms_err(10, output_data["spin_1/2_energies"],
                        output_data["spin_3/2_energies"], output_data["spin_5/2_energies"], 
                        output_data["spin_7/2_energies"], output_data["spin_9/2_energies"], 
                        output_data["spin_11/2_energies"], output_data["spin_13/2_energies"])




    #%%


    ''' 12. PLOT GRAPHS

    - Set a subtitle containing the values of E2PLUS and GSFAC input, if requested.
    - Set which graphs should be plotted (from config, or overwritten below). 
    Any not listed are False by default.

    - If deformation was input as a mesh:
        - Plot filled contours in polar coordinates.
        - Draw a contour line to indicate the perimeter of the region where 
        the ground state spin was correctly reproduced, if requested.
        - Plot data point markers.
            - If experimental data is available, points that agree with experiment 
            (within tolerance) are marked in red.
            - Non-matching points are not plotted (unless there are fewer than 100 data points.)
        
    - If only one of eps/gamma/e2plus is varied:
        - Plot line graphs.
        - Draw a green box around regions that have the correct ground state spin, if requested.
        - Plot a red line to indicate the experimental value, if available.
        
    '''


    #!!! set graph subtitle:
    if code_settings["include_subtitle"]:
        subtitle = r'$E(2^+)$ = ' + str(ptrm_inputs["current_e2plus"]) + '; gsfac = ' + str(ptrm_inputs["gsfac"])
    else:
        subtitle = ''
        
        
    # set which graphs to plot:
    for i in graphs_to_plot:
        if i in output_data:
            output_data[i].plot = graphs_to_plot[i]
        else:
            print("property not recorded, check that it is included in experimental data inputs:\n\t", i)

    # override graph plotting options (useful when running cell by cell): #!!!
        
    # output_data["fermi_indices"].plot = 0
    # output_data["delta"].plot = 0
    # output_data["fermi_energies_mev"].plot = 0
    # output_data["fermi_energies_hw"].plot = 0

    # output_data["gs_mag_moments"].plot = 0
    # output_data["gs_quad_moments"].plot = 0
    # output_data["gs_spin_floats"].plot = 1

    # output_data["spin_1/2_energies"].plot = 0
    # output_data["spin_3/2_energies"].plot = 0
    # output_data["spin_5/2_energies"].plot = 0
    # output_data["spin_7/2_energies"].plot = 0
    # output_data["spin_9/2_energies"].plot = 0
    # output_data["spin_11/2_energies"].plot = 0
    # output_data["spin_13/2_energies"].plot = 0

    # output_data["spin_1/2_mag_moments"].plot = 0
    # output_data["spin_3/2_mag_moments"].plot = 0

    # output_data["rms"].plot = 0

    # output_data["all_energies"].plot = 0
    # output_data["shifted_energies"].plot = 0

    # output_data["gap_9_13"].plot = 0

    # override settings
    # code_settings["mark_exp"] = 1
    # code_settings["mark_exp_tol"] = 0
    # code_settings["mark_points"] = 1
    # code_settings["mark_spin"] = 0


    sub_timer.start()

    data_points["agreed"] = [0]*len(data_points["eps"])
    num_comparisons = 0 

    # start plotting graphs:
    for i in output_data:
        
        prop = output_data[i]
        
        if not(prop.plot):
            continue
        
        ptrm_inputs["current_graph"] = prop.title # makes several later inputs more efficient
        print("plotting graph: %(current_graph)s" % ptrm_inputs) 
        
        if ptrm_inputs["deformation_input"] == "mesh":  
            
            _fig, _ax = plt.subplots(subplot_kw=dict(projection='polar'))
            cax, cbar = gr.draw_contour_plot(_ax, prop, data_points)
            
            legend_handles = []
            
            if code_settings["mark_spin"]:
                
                legend_handles = gr.mark_spin(ptrm_inputs, data_points, output_data["gs_spin_floats"].data, legend_handles, _ax)
                
            # plot the data point markers, with comparison to experiment if possible
                
            legend_handles = gr.plot_points(data_points, prop, legend_handles, cbar, code_settings)
            if np.isfinite(prop.experimental_data).all() and code_settings["mark_exp"]: 
                num_comparisons += 1
            
            
            gr.format_fig('polar', _ax, legend_handles, '%(current_graph)s of %(nucleus)s' % ptrm_inputs, subtitle)
            
            plt.show()
            
            
        
        elif (ptrm_inputs["deformation_input"] ==  "gamma" 
            or ptrm_inputs["deformation_input"] == "eps"
            or len(data_points["e2plus"]) > 1):
            
            # set which paramters are varied and which are constant
            var_sym, var, fix_sym, fix = gr.assign_parameters(ptrm_inputs, data_points)
            
            _fig, _ax = plt.subplots() 
            
            legend_handles = []
            legend_handles, legend_title = gr.plot_line_data(data_points, prop, var, fix_sym, fix, legend_handles)

            
            # if experimental data is available, plot it in red for easy comparison
            if np.isfinite(prop.experimental_data).all() and not prop.num == "all": 
                legend_handles = gr.plot_exp_line(prop, code_settings, var, legend_handles)

                
            # mark the range in which the correct ground state spin was calculated
            if code_settings["mark_spin"]==1:
                
                correct_spin_range = gr.find_correct_spin(output_data["gs_spin_floats"].data, experimental_data["gs_spin_float"])
                if len(correct_spin_range) > 0:
                    spin = gr.plot_correct_spin(correct_spin_range, var, ptrm_inputs["step"], prop)
                    legend_handles.append(spin)
                        
            gr.format_fig('linear', _ax, list(reversed(legend_handles)), 
                        '%(current_graph)s in %(nucleus)s' % ptrm_inputs, subtitle, 
                        varied=var, x_label=var_sym, y_label=prop.axis_label, 
                        legend_title=legend_title)
            
            if prop.prop == "delta":
                _ax.set_ylim([0.2,1]) 
                
            if prop.cbar_tick_labels:        # then format for discrete values
                _ax.set_yticks(prop.cbar_ticks)
                _ax.set_yticklabels(prop.cbar_tick_labels)
            
            plt.show()
            
            
        

    #%%
    ''' 13. ASSESS AGREEMENT OF CALCULATIONS WITH EXPERIMENT

    - Print information about the best agreement and its location.
    - Plot a graph to show data point agreement across all data points.
    - Print the mean energies of each level and the mean gs moments, with standard error.
    - Print the total runtime.

    '''
        
    gr.check_agreement(code_settings["print_details"], data_points, num_comparisons)

    agreement = st.PropertyData(data_points["agreed"], "Agreement of Data Points With Experimental Data")
    agreement.contour_levels = np.arange(0, num_comparisons+2, dtype=int) #fn.calc_contour_levels(agreement.data)
    agreement.cbar_ticks = gr.calc_cbar_ticks(agreement.contour_levels)
    agreement.cbar_tick_labels = list(np.arange(0, num_comparisons+1, dtype=int)) #fn.calc_cbar_tick_labels(agreement.data, "int")
    agreement.experimental_data = np.NaN
    agreement.error_tolerance = np.NaN

    agreement.plot = 0
    if ptrm_inputs["deformation_input"] == "mesh" and agreement.plot:  
        gr.plot_agreement(code_settings, agreement, data_points, output_data, subtitle)
        
    print("\n******** mean and standard error in the mean ******")

    anyl.report_mean(output_data["spin_1/2_energies"], code_settings["print_details"])
    anyl.report_mean(output_data["spin_3/2_energies"], code_settings["print_details"])
    anyl.report_mean(output_data["spin_5/2_energies"], code_settings["print_details"])
    anyl.report_mean(output_data["spin_7/2_energies"], code_settings["print_details"])
    anyl.report_mean(output_data["spin_9/2_energies"], code_settings["print_details"])
    anyl.report_mean(output_data["spin_11/2_energies"], code_settings["print_details"])
    anyl.report_mean(output_data["spin_13/2_energies"], code_settings["print_details"])
    anyl.report_mean(output_data["gs_mag_moments"], code_settings["print_details"])
    anyl.report_mean(output_data["gs_quad_moments"], code_settings["print_details"])

    # note how long it took
    sub_timer.stop()
    main_timer.stop()
    print("\n****************************************************************************************")
    print("Finished plotting graphs in time = %.2f seconds" % (sub_timer.get_lapsed_time()))
    print("total runtime = %.2f seconds" % (main_timer.get_lapsed_time()))
    print("****************************************************************************************\n")





if __name__ == "__main__":
    main()