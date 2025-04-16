#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import functions as fn                               # my own module file of functions
import structs as st                                 # my own module file of structs (classes, and read-only dicts)
import optuna_functions as opfn

import optuna                                        # for parameter optimisation
import subprocess                               # for calling shell scripts to run 



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

def run_model(inputs, data_point, verbose, nucleus):
 
    '''4.  GAMPN    '''
    
    data_point["current_orbitals"] = inputs["initial_orbitals"]
    fn.write_file("../"+nucleus+"/Inputs/GAM_"+data_point["file_tag"]+".DAT", st.get_template("gampn") % (inputs | data_point ))
    

    write_and_run("gampn", nucleus, data_point["file_tag"], verbose)
    
    
    # set up dictionary to store data 
    output_data = {}

    
    lines = fn.read_file("../"+nucleus+"/Outputs/GAM_"+data_point["file_tag"]+".OUT")
    
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
    data_point["current_orbitals"] = (fn.find_orbitals(f_index, inputs["nu"], 
                                   inputs["par"], f_energy_hw, f_parity, lines))
    
    
    # Re-run gampn with the new set of orbitals, so that the strong-coupling basis can be maximised. 
        
    if data_point["current_orbitals"] != inputs["initial_orbitals"]:
        fn.write_file("../"+nucleus+"/Inputs/GAM_"+data_point["file_tag"]+".DAT", st.get_template("gampn") % (inputs | data_point )) # write DAT file
        write_and_run("gampn", nucleus, data_point["file_tag"], verbose) # write bash script and run
    
    
    
    ''' 6. ASYRMO    '''
        
    fn.write_file("../"+nucleus+"/Inputs/ASY_"+data_point["file_tag"]+".DAT", st.get_template("asyrmo") % (inputs | data_point ))
    write_and_run("asyrmo", nucleus, data_point["file_tag"], verbose) # write bash script and run
    
    
    lines = fn.read_file("../"+nucleus+"/Outputs/ASY_"+data_point["file_tag"]+".OUT")
    
    if not "PARTICLE-ROTOR  MODEL" in lines[0]: # then something has gone wrong
        raise RuntimeError("File " + data_point["file_tag"] + " raised error in ASYRMO output: \n" + lines[0] )
        
    output_data["delta"] = fn.get_delta(lines) # this also checks for the "SORRY I FOUND NO SOLUTIONS" error.
            
    ''' 8. PROBAMO    '''
    
    fn.write_file("../"+nucleus+"/Inputs/PROB_"+data_point["file_tag"]+".DAT", st.get_template("probamo") % (inputs | data_point ))
    write_and_run("probamo", nucleus, data_point["file_tag"], verbose) # write bash script and run
    
    
    ''' 9. READ PROBAMO.OUT FILE
    
        
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
    
    '''
    

    file_data = {}
    
    lines = fn.read_file("../"+nucleus+"/Outputs/Prob_"+data_point["file_tag"]+".OUT")

    for line in lines:
        
        line_data = fn.read_data(line)  # get the spin, energy, and magnetic moment from this line if it is a static moment, else return False
        
        if not(line_data):
            continue # to next line in file
        
        # sort line_data into file_data according to its spin
        file_data = fn.sort_by_spin(line_data, file_data)
        
        # additionally save data that corresponds to the experimental ground state and input excited states
        file_data = fn.sort_by_expectation(line_data, file_data, inputs)

    
    return output_data | fn.missing_data(file_data, inputs)


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


def objective(trial, inputs, experimental):
    
    data_point = {}

    data_point["current_eps"] = trial.suggest_float("EPS", 0.001, 0.5)
    data_point["current_gamma"] = trial.suggest_float("GAM", 0.0, 60.0)

    data_point["file_tag"] = "e%.3f_g%.1f_p%.3f_%s" % (data_point["current_eps"],  data_point["current_gamma"],
                                        inputs["current_e2plus"], inputs["nucleus"])

    data_point["current_f002"] = "f002_"+data_point["file_tag"]+".dat"
    data_point["current_f016"] = "f016_"+data_point["file_tag"]+".dat"
    data_point["current_f017"] = "f017_"+data_point["file_tag"]+".dat"
    data_point["current_f018"] = "f018_"+data_point["file_tag"]+".dat"
    
    output_data = run_model(inputs, data_point, True, inputs["nucleus"])
    
    loss = compute_loss(output_data, experimental)
    
    return loss
