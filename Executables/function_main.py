#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import functions as fn                               # my own module file of functions
import structs as st                                 # my own module file of structs (classes, and read-only dicts)

import subprocess                               # for calling shell scripts to run 

#%%

def compute_loss(test_output, experimental):
    
    loss = 0
    
    
    # if the gs spin is wrong, that is a very bad sign, so assign this with heavy loss
    loss += 10 * int(test_output["gs_spin_float"] != experimental["gs_spin_float"])
    
    # the further the calculated magnetic moment is from the experimental value, the higher the loss should be 
    loss += (test_output["gs_mag_mom"] - experimental["gs_mu"])**2
    
    # take the lowest energy level of each spin for now
    for j in experimental:
        
        if "jp_" not in j: # ignore experimental values that aren't energy levels 
            continue
        
        spin = j[3:-1]
        
        lowest_calculated = test_output["spin_"+spin+"_energies"][0]
        lowest_experimental = experimental[j][0]
        
        loss += (lowest_calculated-lowest_experimental)**2
        
    return loss


#%%
    
    
def write_and_run(program, nucleus, file_tag, verbose):
     """
     The inner function of the configure_script_writer closure. This is the 
     function that is returned from the closure, which can then be called.
     
     Dynamically defines a script writer customised to the program requested.
     
     For each batch of data points:
         - Determines the file path to the folder from which the batch will be run.
         - Gets the list of file tags corresponding to that batch.
         - Writes a bash script to execute the program for all the files in the batch.
         - Starts the bash script as a subprocess.
     
     When all the scripts have been set running, the function waits for them 
     all to complete before returning.
     
     Parameters
     ----------
     program : string
         "gampn" to run gampn
         "asyrmo" to run asyrmo
         "probamo" to run probamo

     Raises
     ------
     ValueError
         Occurs if a program other than gampn, asyrmo, or probamo is requested.

     Returns
     -------
     None.

     """
     if program == "gampn": abr = "GAM"
     elif program == "asyrmo": abr = "ASY"
     elif program == "probamo": abr = "PROB"
     else: raise ValueError("Input program ["+ program +"] not supported.")

     script_file_path = "../"+nucleus+"/Scripts/no_batch/Run"+program.upper()+"_"+file_tag+".sh"
     
     script_text = ("\n(cd ../"+nucleus+"/Run/no_batch;" +  # start a subprocess and move to the Run folder for this batch.
                    " ./../../../Executables/MO/" + program +       # call the program from the Run folder.
                    " < ../../Inputs/"+ abr +"_"+ file_tag +".DAT;" +   # input the relevant .DAT file to the program.
                    " cp "+program.upper()+ ".out"+                 # copy the default .OUT file... 
                    " ../../Outputs/"+ abr +"_" + file_tag + ".OUT)")   # ...to a new .OUT file with a more descriptive name, in the Outputs folder.

     if verbose:
         script_text += ("\n\necho message from terminal: " +
                     "finished running "+ program +" for file "+ file_tag)
     
     with open(script_file_path, 'w') as f:
          f.write(script_text)

     subp = subprocess.Popen(["sh", script_file_path])     
     # asynchronous call to start the program as a subprocess
         
     allowed_time = 1 # wait a maximum of 1 second for it to finish running before assuming that it's hanging
     subp.wait(allowed_time)             


#%%

def run_model(inputs, file_tag, data_point):
 
    '''4. WRITE GAMPN.DAT FILES 
    
    - Loops through the data points to be tested:
        - Creates a tag referencing the nucleus, and the values of eps, gamma, and e2plus.
        - Uses string formatting with the inputs dictionary to write the .DAT file.
        - Creates/overwrites a .DAT file in the Inputs directory folder, named using the file tag. 
    
    '''
    
               
    fn.write_file("../"+nucleus+"/Inputs/GAM_"+file_tag+".DAT", st.get_template("gampn") % inputs)
    

    ''' 5. WRITE AND RUN BASH SCRIPT TO EXECUTE GAMPN 
    
    - Calculate batch settings:
        - How to divide up the data points into batches.
        - The maximum allowed time for the batch to run before assuming that it is hanging.
    - Configure a script writer.
    - Run the batches. The .OUT files are generated in the Outputs directory folder.
    
    '''
    

    write_and_run("gampn", nucleus, file_tag, verbose)
    
    
    ''' 5. READ GAMPN.OUT FILE 
    
    - For each data point (i.e. each GAMPN.OUT file):
        - Read the value of EFAC (the conversion factor from hw to eV).
        - Determine the line number of the fermi level orbital in the GAMPN.OUT file.
        - Read that line from GAMPN.OUT and get the energy, parity, and level index.
        - Generate the orbitals input for asyrmo (e.g. orbitals_string = "+11 19 20 21 22 23 24 25 26 27 28 29")
        
    '''
    
    # set up dictionary to store data 
    output_data = {}

    
    lines = fn.read_file("../"+nucleus+"/Outputs/GAM_"+file_tag+".OUT")
    
    # get the EFAC value (conversion factor from hw to MeV)
    inputs["efac"] = fn.get_efac(lines)
    
    # locate the fermi level and extract data
    fermi_level_line = fn.get_sp_level(lines, inputs["fermi_level"], '0')
    
    f_parity, f_energy_hw, f_index = fn.get_info(fermi_level_line)
    
    output_data["fermi_parities"] = f_parity                             
    output_data["fermi_energies_hw"] = f_energy_hw   
    output_data["fermi_energies_mev"]= f_energy_hw * inputs["efac"]
    output_data["fermi_indices"] = f_index
    
    # generate the orbitals input for asyrmo
    
    # version 5: dynamically finds the orbitals nearest to the fermi level in energy:
    data_point["asyrmo_orbitals"].append(fn.find_orbitals(f_index, inputs["nu"], 
                                   inputs["par"], f_energy_hw, f_parity, lines))
    

    # Re-run gampn with the new set of orbitals, so that the strong-coupling basis can be maximised. 
    # No need to re-read outputs, as the values read manually (energy levels etc) will be unaffected,
    # only the matrix elements (read automatically by asyrmo) will be affected.
        

    fn.write_file("../"+nucleus+"/Inputs/GAM_"+file_tag+".DAT", st.get_template("gampn") % inputs) # write DAT file
    write_and_run("gampn", nucleus, file_tag, verbose) # write bash script and run
    
    
    
    ''' 6. WRITE AND RUN ASYRMO 
    
    - Use the existing list of file tags to write a .DAT file for each data point.
    - Use the existing script writer to write and run asyrmo; dividing up the batches (as for gampn). 
    
    '''
        
    fn.write_file("../"+nucleus+"/Inputs/ASY_"+data_point["file_tags"]+".DAT", st.get_template("asyrmo") % inputs)
    write_and_run("asyrmo", nucleus, file_tag, verbose) # write bash script and run
    
    
    ''' 7. READ ASYRMO 
    
    - Read the output files and check for the "NO DECOUPLING PARAMETERS CALCULATED" error.
    - Record the value of the DELTA parameter.
    - Check for the "SORRY I FOUND NO SOLUTIONS" error, and exclude those files from future analysis.
    
    '''
    
    lines = fn.read_file("../"+nucleus+"/Outputs/ASY_"+file_tag+".OUT")
    
    if not "PARTICLE-ROTOR  MODEL" in lines[0]: # then something has gone wrong
        raise RuntimeError("File " + file_tag + " raised error in ASYRMO output: \n" + lines[0] )
        
    output_data["delta"] = fn.get_delta(lines) # this also checks for the "SORRY I FOUND NO SOLUTIONS" error.
            
    

    ''' 8. WRITE AND RUN PROBAMO 
    
    - Use the existing list of file tags to write a .DAT file for each data point.
    - Use the existing script writer to write and run probamo; dividing up the batches as for gampn. 
    
    '''
    
    fn.write_file("../"+nucleus+"/Inputs/PROB_"+data_point["file_tags"]+".DAT", st.get_template("probamo") % inputs)
    write_and_run("probamo", nucleus, file_tag, verbose) # write bash script and run
    
   
    
    ''' 9. READ PROBAMO.OUT FILE
    
    For each file:
        
    - Read each line:
        - If it is a static transition, read the spin, energy, and magnetic moment.
        - Else ignore this line and move to the next.
        
    - Sort the file data into categories:
        - Group lines by spin (e.g. spin 1/2 energies, spin 1/2 magnetic dipole moments, etc).
        - Additionally (separately) group lines by association with experimental data,
          (e.g. ground state energies, first excited state energies, etc).
          (this assumes that each input experimental excited state is the yrast state of that spin).
        - Fill missing gaps with NaN values, such that the same set of properties has 
          been recorded for every data point 
          
          (and if e.g. one data point found three spin 1/2 states, then all data points 
          should have a list of three spin 1/2 states, even if some of them are NaN)
        
        - Restructure the data set and save separately (now each property is recorded
          as a list of values for all data points, rather than each data point having
          a list of properties associated with it). 
          
    
    '''
    

    file_data = {}
    
    lines = fn.read_file("../"+nucleus+"/Outputs/Prob_"+data_point["file_tags"]+".OUT")

    for line in lines:
        
        line_data = fn.read_data(line)  # get the spin, energy, and magnetic moment from this line if it is a static moment, else return False
        
        if not(line_data):
            continue # to next line in file
        
        # sort line_data into file_data according to its spin
        file_data = fn.sort_by_spin(line_data, file_data)
        
        # additionally save data that corresponds to the experimental ground state and input excited states
        file_data = fn.sort_by_expectation(line_data, file_data, inputs)
    
    output_data["property_data"] = fn.missing_data(file_data, inputs)
    
    return output_data
    








'''
    
    output_data = fn.restructure_data(data_point["property_data"], inputs["ispin"], verbose)
    
    # get energy gap between 9/2 and 13/2
    output_data["gap_9_13"] = fn.find_gaps(output_data["spin_9/2_energies"], 2, output_data["spin_13/2_energies"], 1, 20) #!!!
    
    
    # ensure all data sets have the same size and shape
    # and mask ill-defined data poitns with reference to the DELTA data (any NaN values are ill-defined).
    
    mask = np.array([0 if np.isnan(x) else 1 for x in output_data["delta"]])
    
    for d in output_data:
        if isinstance(output_data[d][0], list):
            output_data[d] = fn.fill_gaps(output_data[d])
            list_mask = np.transpose(np.tile(mask, (np.size(output_data[d][0]),1)))
            
            output_data[d] = np.where(list_mask == 0, np.NaN, output_data[d])
        
        else: 
            output_data[d] = np.where(mask == 0, np.NaN, output_data[d])
    
    _output_data_dict = output_data # save a copy of the original before it's overwritten, so that the code can be run cell-by-cell without errors.
    
    
    10. PREPARE TO PLOT GRAPHS 
    
    - Record each data set in an instance of class PropertyData
    - Calculate graph plotting attributes and store within the class
    
    - Raise a ValueError if the property isn't recognised 
      (i.e. if more data sets are read in the future, they cannot be plotted without
       first hard-coding the calculation of things like axis labels, contour levels,
       colour bar ticks, etc).
    
    
    
    # convert output_data from a dictionary of lists to a dictionary of PropertyData objects 
    output_data = {}
    for p in _output_data_dict:
        
        output_data[p] = st.PropertyData(_output_data_dict[p], p)
        
        # calculate contour levels, colour bar ticks and labels, 
        # and assign experimental values and error tolerance if available.
        
        if output_data[p].prop == "energies" and not(output_data[p].sort=="Fermi"): 
            
            output_data[p].contour_levels = 10
            output_data[p].cbar_tick_labels = 0
            
            if output_data[p].sort == "gap":
                
                output_data[p].experimental_data, output_data[p].error_tolerance = fn.try_experimental(inputs, p, inputs["gap_en_tol"])
            
            elif output_data[p].sort == "Excited State ":
                
                output_data[p].experimental_data, output_data[p].error_tolerance = fn.try_experimental(inputs, "x" + output_data[p].num + "_energy", inputs["abs_en_tol"])
            
            elif output_data[p].sort == "Spin ":
                
                output_data[p].experimental_data, output_data[p].error_tolerance = fn.try_experimental(inputs, "jp_"+ output_data[p].num + inputs["par"], inputs["abs_en_tol"])
                
            else: raise ValueError("property not recognised: " + p)
            
        elif output_data[p].prop == "mag_moments":
            
            output_data[p].contour_levels = 6
            output_data[p].cbar_tick_labels = 0
                
            if output_data[p].sort == "Excited State ":
                
                output_data[p].experimental_data, output_data[p].error_tolerance = fn.try_experimental(inputs, "x" + output_data[p].num + "_mu", inputs["mu_tol"])
        
            elif output_data[p].sort == "Spin ":
                
                output_data[p].experimental_data = np.NaN
                output_data[p].error_tolerance = np.NaN
                
            elif output_data[p].sort == "Ground":
                
                output_data[p].experimental_data, output_data[p].error_tolerance = fn.try_experimental(inputs, "gs_mu", inputs["mu_tol"])
                
            else: raise ValueError("property not recognised: " + p)
            
        
        elif output_data[p].sort == "Ground":
            
            if output_data[p].prop == "spin_floats": 
            
                output_data[p].contour_levels = fn.calc_contour_levels(output_data[p].data)
                output_data[p].experimental_data, output_data[p].error_tolerance = fn.try_experimental(inputs, "gs_spin_float", inputs["mu_tol"])
                output_data[p].cbar_tick_labels = fn.calc_cbar_tick_labels(output_data[p].data, "half")
                output_data[p].cbar_ticks = fn.calc_cbar_ticks(output_data[p].contour_levels)
                
            elif output_data[p].prop == "spin_strings": continue
            else: raise ValueError("property not recognised: " + p)
        
        elif output_data[p].sort == "Fermi":
        
            if output_data[p].prop == "indices": 
                
                output_data[p].contour_levels = fn.calc_contour_levels(output_data[p].data)
                output_data[p].cbar_tick_labels = fn.calc_cbar_tick_labels(output_data[p].data, "int")
                output_data[p].cbar_ticks = fn.calc_cbar_ticks(output_data[p].contour_levels)
            
            elif output_data[p].prop == "energies":
            
                output_data[p].contour_levels = 10
                output_data[p].cbar_tick_labels = 0
            
            elif output_data[p].prop == "parities":
            
                output_data[p].contour_levels = 2
                output_data[p].cbar_tick_labels = 0
            
            else: raise ValueError("property not recognised: " + p)
                
            output_data[p].experimental_data = np.NaN
            output_data[p].error_tolerance = np.NaN
        
        elif output_data[p].prop == "delta":
            output_data[p].contour_levels = 10
            output_data[p].cbar_ticks = 0
            output_data[p].cbar_tick_labels = 0
            output_data[p].experimental_data = np.NaN
            output_data[p].error_tolerance = np.NaN
            
        else: raise ValueError("property not recognised: " + p)
        
    del [p]
    
    output_data["all_energies"] = fn.collate_energy_data(output_data, len(data_points["file_tags"]))
    
     11. PLOT GRAPHS
    
    - Set which graphs should be plotted. Any not listed are False by default.
    
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
        
    
    
    #!!! set graph subtitle:
    subtitle = "e2plus = 0.13, gsfac = 0.7 \n second 9/2 and first 13/2"# "no coriolis attenuation (chsi = eta = 1.0) \n 15 orbitals calculated dynamically"
    
    #!!! set which graphs to plot:
    
    output_data["fermi_indices"].plot = 0
    output_data["delta"].plot = 0
    output_data["fermi_energies_mev"].plot = 0
    output_data["fermi_energies_hw"].plot = 0
    
    output_data["gs_mag_moments"].plot = 1
    output_data["gs_spin_floats"].plot = 1
    
    output_data["spin_1/2_energies"].plot = 0 
    output_data["spin_3/2_energies"].plot = 1
    output_data["spin_5/2_energies"].plot = 0
    output_data["spin_7/2_energies"].plot = 0
    output_data["spin_9/2_energies"].plot = 1
    output_data["spin_11/2_energies"].plot = 0
    output_data["spin_13/2_energies"].plot = 1
    
    output_data["spin_1/2_mag_moments"].plot = 0
    output_data["spin_3/2_mag_moments"].plot = 0
    
    output_data["x1_energies"].plot = 0
    output_data["x2_energies"].plot = 0
    output_data["x3_energies"].plot = 0
    
    output_data["x1_mag_moments"].plot = 0
    
    output_data["all_energies"].plot = 1
    
    output_data["gap_9_13"].plot = 1
    
    
    data_points["agreed"] = [0]*len(data_points["eps"])
    num_comparisons = 0 
    
    sub_timer.start()
    
    
    
    
    
    # start plotting graphs:
    for g in output_data:
        
        prop = output_data[g]
        if not(prop.plot):
            continue
        inputs["current_graph"] = prop.title # makes several later inputs more efficient
        print("plotting graph: %(current_graph)s" % inputs) 
        
        if inputs["deformation_input"] == "mesh":  
            
            fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
            cax, cbar = fn.draw_contour_plot(ax, prop, data_points)
            
            
            legend_handles = []
            
            if inputs["mark_spin"]:
                
                legend_handles = fn.mark_spin(inputs, data_points, output_data["gs_spin_floats"].data, legend_handles, ax)
                
            # plot the data point markers, with comparison to experiment if possible
                 
            legend_handles = fn.plot_points(data_points, prop, legend_handles, cbar, inputs)
            if np.isfinite(prop.experimental_data).all() and inputs["mark_exp"]: 
                num_comparisons += 1
            
            
            fn.format_fig('polar', ax, legend_handles, '%(current_graph)s of %(nucleus)s' % inputs, subtitle)
            
            plt.show()
            
            
           
        elif (inputs["deformation_input"] ==  "gamma" 
              or inputs["deformation_input"] == "eps"
              or len(data_points["e2plus"]) > 1):
            
            # set which paramters are varied and which are constant
            var_sym, var, fix_sym, fix = fn.assign_parameters(inputs, data_points)
            
            fig, ax = plt.subplots() 
            
            legend_handles = []
            
            # now plot the actual data
            if len(data_points["file_tags"]) < 100: marker_size = 5
            else: marker_size = 1 # use smaller markers if the data set is large
            
            if prop.sort == "Spin ":
                
                if prop.num == "all":
                    legend_handles, legend_title = fn.plot_all_energies(prop, var, legend_handles, marker_size, fix_sym, fix[0])
                else:
                    legend_handles, legend_title = fn.plot_multi_lines(prop, var, legend_handles, marker_size, fix_sym, fix[0])
                
            else:
                data, = plt.plot(var, prop.data, 'k-x', markersize=marker_size, label="%s = %s" % (fix_sym, fix[0]))
                legend_handles.append(data)
                legend_title = ""
              
            # if experimental data is available, plot it in red for easy comparison
            if np.isfinite(prop.experimental_data).all():  
                
                if prop.sort == "Spin ":
                    for e in range(len(prop.experimental_data)):
                        if inputs["mark_exp"]:
                            exp, = plt.plot(var, np.full(len(var), prop.experimental_data[e]), 'r-', label="experimental value")
                        if inputs["mark_exp_tol"]:
                            exp_tol = ax.axhspan(prop.experimental_data[e]-prop.error_tolerance, prop.experimental_data[e]+prop.error_tolerance, facecolor='r', alpha=0.2, label="experimental range with tolerance")
                else:
                    if inputs["mark_exp"]:
                        exp, = plt.plot(var, np.full(len(var), prop.experimental_data), 'r-', label="experimental value")
                    if inputs["mark_exp_tol"]:
                        exp_tol = ax.axhspan(prop.experimental_data-prop.error_tolerance, prop.experimental_data+prop.error_tolerance, facecolor='r', alpha=0.2, label="experimental range with tolerance")
    
                if inputs["mark_exp"]:
                    legend_handles.append(exp)
                if inputs["mark_exp_tol"]:
                    legend_handles.append(exp_tol)
                
        
            # mark the range in which the correct ground state spin was calculated
            if inputs["mark_spin"]==1:
                
                correct_spin_range = fn.find_correct_spin(output_data["gs_spin_floats"].data, experimental["gs_spin_float"])
                if len(correct_spin_range) > 0:
                    spin = fn.plot_correct_spin(correct_spin_range, var, inputs["step"], prop)
                    legend_handles.append(spin)
                
                        
            fn.format_fig('linear', ax, list(reversed(legend_handles)), 
                           '%(current_graph)s in %(nucleus)s' % inputs, subtitle, 
                           varied=var, x_label=var_sym, y_label=prop.axis_label, 
                           legend_title=legend_title)
            
            if prop.prop == "delta":
                ax.set_ylim([0.2,1]) 
                
            if prop.cbar_tick_labels:        # then format for discrete values
                
                ax.set_yticks(prop.cbar_ticks)
                ax.set_yticklabels(prop.cbar_tick_labels)
            
            plt.show()
            
            del [var, var_sym, fix, fix_sym, legend_title, legend_handles]
            
           
    del [g]
    
    12. ASSESS AGREEMENT OF CALCULATIONS WITH EXPERIMENT
    
    - Print information to the console about the best agreement and its location.
    - Plot a graph to show data point agreement across all data points.
    

    
    inputs["mark_points"] = 0
        
    fn.check_agreement(verbose, data_points, num_comparisons)
    
    agreement = st.PropertyData(data_points["agreed"], "Agreement of Data Points With Experimental Data")
    agreement.contour_levels = np.arange(0, num_comparisons+2, dtype=int) #fn.calc_contour_levels(agreement.data)
    agreement.cbar_ticks = fn.calc_cbar_ticks(agreement.contour_levels)
    agreement.cbar_tick_labels = list(np.arange(0, num_comparisons+1, dtype=int)) #fn.calc_cbar_tick_labels(agreement.data, "int")
    agreement.experimental_data = np.NaN
    agreement.error_tolerance = np.NaN
    
    agreement.plot = 1
    
    
    if inputs["deformation_input"] == "mesh" and agreement.plot:  
        
        inputs["current_graph"] = agreement.title
        print("plotting graph: %(current_graph)s" % inputs) 
        
        fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
        cax, cbar = fn.draw_contour_plot(ax, agreement, data_points)
        
        legend_handles = []
        
        if inputs["mark_spin"]:
            legend_handles = fn.mark_spin(inputs, data_points, output_data["gs_spin_floats"].data, legend_handles, ax)
        
        if inputs["mark_points"]:
            legend_handles =  fn.plot_points_without_experiment(data_points, legend_handles)
               
        
        fn.format_fig('polar', ax, legend_handles, '%(current_graph)s of %(nucleus)s' % inputs, subtitle)
        
        plt.show()
        
    
    # note how long it took
    sub_timer.stop()
    timer.stop()
    print("\n****************************************************************************************")
    print("finished plotting graphs in time = %.2f seconds" % (sub_timer.get_lapsed_time()))
    print("total runtime = %.2f seconds" % (timer.get_lapsed_time()))
    print("****************************************************************************************\n")

    '''
    
    
    
    
