#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  8 10:07:25 2024

@author: katyrr


run from Outputs directory folder in terminal, with: 
    python3 ../../../Executables/write_inputs.py

or run from Spyder console with: 
    runfile('/Users/katyrr/Downloads/MSci Project/Code/Executables/write_inputs.py', wdir='/Users/katyrr/Downloads/MSci Project/Code/Tl207/MO/Outputs')

to debug (use breakpoints) run from Spyder console with: 
    debugfile('/Users/katyrr/Downloads/MSci Project/Code/Executables/write_inputs.py', wdir='/Users/katyrr/Downloads/MSci Project/Code/tl207/MO/Outputs')

"""

#%%


import numpy as np                                                              # for np.arange(start, end, step)
import subprocess                                                               # for calling shell scripts to run 
import math                                                                     # for ceil(num)
import matplotlib.pyplot as plt                                                 # for plotting graphs
import time                                                                     # for checking how long it took to run
from matplotlib.ticker import FuncFormatter                                     # for formatting axis ticks
import matplotlib.tri as tri                                                    # for manual triangulation before drawing a contour plot

plt.rcParams['figure.dpi'] = 150
        
templates = {}

templates["gampn"] = '''
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

templates["asyrmo"] = '''
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

templates["probamo"] = '''
'%(current_f017)s' '%(current_f018)s'              FILE17,FILE18
1,0                               ipkt,iskip
%(istrch)s                                 isrtch
%(Z)s,%(A)s                           Z,AA
0,%(cutoff)s,1,0.75,-1                ISPEC,CUTOFF,IQ,GSFAC,GR
0.0000, 0.000,0.000, 0.000        BS2,BS4 (FOR S-STATE), BS2,BS4(P-STATE)
'''

variable_types = {}

variable_types["bool"] = ["mark_spin"]

variable_types["int"] = ["A", "Z", "num_to_record", "num_orbs", "nu"]

variable_types["float"] = ["gs_energy", "x1_energy", "x2_energy", "x3_energy",
                           "gs_mu", "x1_mu"]

timer_data = {}

#%%

''' FUNCTIONS FOR READING CONFIG AND RUNNING CODES '''



def spin_string_to_float(spin_string):
    if spin_string[1] == "/":
        spin_float = float(spin_string[0])/2
    elif spin_string[2] == "/":
        spin_float = float(spin_string[0:2])/2 
    else:
        raise ValueError("Cannot parse spin. \nCheck it has been input in " +
                         "the format 'n/2' where n is an int.")
    return spin_float

def spin_float_to_string(spin_float):
    numerator = float(spin_float)*2
    if str(numerator)[-1] != "0":
        raise ValueError("Cannot convert float to fraction. \nCheck it has " + 
                         "been input as a half-integer float.")
    else:
        return str(int(numerator))+"/2"
    
def range_to_list(range_string):
    split_range = [float(n) for n in range_string.split(",")]
    
    
    if len(split_range) > 1 :  
        step = split_range[2]           
        range_list =  np.arange(split_range[0],                          
                                 split_range[1]+step, 
                                 step)
    else : range_list = [split_range[0]]
    
    return range_list

def get_range_step(range_string):
    split_range = [float(n) for n in range_string.split(",")]
    return split_range[2]     
    


def configure_script_writer(file_tags,  num_batches, 
                            num_per_batch, allowed_time):                       # these arguments will be the same throughout the execution, so we can define a family of script writers (one for each of gam/asy/prob)
    
    def run_script_batches(program):                                            # after configuring the script writer family, this is what we actually call (it's a "closure" function)
        
        if program == "gampn": abr = "GAM"
        elif program == "asyrmo": abr = "ASY"
        elif program == "probamo": abr = "PROB"
        else: raise ValueError("Input program ["+ program +"] not supported.")
        
        def write_script_batch(file_tag_batch, batch_index, file_path):         # a sub-function, called in the loop below
            
            batch_folder = "Batch"+str(batch_index)
            script_text = ""
            
            for file in file_tag_batch :
                script_text += ("\n(cd "+batch_folder+
                                "; ./../../../../Executables/MO/" + program + 
                                " < ../../Inputs/" +abr + "_" +file + ".DAT)")
                script_text += ("\ncp "+batch_folder+"/"+program.upper()+
                                ".out "+abr+"_"+file+".OUT")                        # copy GAMPN.out to a new txt file with a more descriptive name
            
            script_text += ("\n\necho message from terminal: " +
                            "finished running " + program + " batch " + str(batch_index))
            
            script_file = open(file_path, 'w')
            script_file.write(script_text)
            script_file.close() 


        subprocesses = {}

        for b in range(num_batches):
            file_path = "../Run"+program.upper()+"_"+str(b+1)+".sh"
            batch_file_tags = file_tags[(b*num_per_batch):((b+1)*num_per_batch)]
            write_script_batch(batch_file_tags, b+1, file_path)
            subprocesses[(program+"_"+str(b+1))] = subprocess.Popen(["sh", file_path])     # asynchronous call to start gampn as a subprocess

        for b in range(num_batches):
            subprocesses[(program+"_"+str(b+1))].wait(allowed_time)                        # wait to ensure it has finished before starting to read outputs, if it takes longer than the time limit seconds, throw an error to catch hangs.
            
    return run_script_batches



''' OPEN AND READ CONFIG FILE '''
# ignores empty lines, and lines beginning with * (to mark a comment)
# counts and outputs the number of non-empty and non-comment lines (should be equal to the number of variables input in the config file)
# expects each line to have the format:    var_name value
# outputs a warning if the format is unexpected, otherwise splits the string at " " and assigns the variable to a dictionary
# tests whether deformation has been input as a range, mesh, or single value, and creates a list of values to test (which may contain just one value)

timer_data["start"] = time.time()

print("reading from config file...") 
with open("../Inputs/config.txt", 'r') as f:
    config_lines = f.readlines()
del f

line_count = 0                                                                  # start a counter for all lines  (including empty lines and comments)
inputs = {}                                                                     # create an empty dictionary to hold settings for input

for line in config_lines:                                                              # proccess each line
    line_count += 1
    
    # remove comments and blank lines
    line_string = line.strip()                                                  # extract text from line
    if line_string == "" : continue                                             # skip blank rows
    if line_string[0] == "*" : continue                                         # skip comment lines
    
    split_string = line_string.split(" ")                                       # split into name and value
        
    for n in range(len(split_string)):                                          # remove inline comments 
        word = split_string[n]
        if word[0] == '*':             
            del split_string[n:]
            break  
    
    # warn if the input looks unexpected
    if (split_string[0]=="eps" 
        or split_string[0]=="gamma" 
        or split_string[0]=="single"):
          expected_num_words = 3                                                # these variables have two values associated with the variable name
    else: expected_num_words = 2                                                # all other variables have only one value associated with the variable name
    
    if len(split_string)>expected_num_words : 
        raise ValueError("Line "+str(line_count)+ ''' has too many words! 
              Check var name contains no spaces, and separate list values 
              with commas not spaces. \nBad input: ''' + line_string)                                
        
    elif len(split_string)<expected_num_words :
        raise ValueError("Line "+str(line_count)+''' is missing either 
              variable name or value. Check they are both present
              and seperated with a space. \nBad input: ''' + line_string)         
    
    # save inputs that may have a range of values   
    if (split_string[0]=="eps" 
        or split_string[0]=="gamma" 
        or split_string[0]=="single"):                                          # each of these methods of inputting deformation contain two values: eps/range of eps, and gamma/range of gamma
        
        # warn if multiple sets of deformation parameters have been input
        if ("deformation_input" in inputs):
            raise ValueError("Deformation parameter sets have been input " + 
                  "multiple times. \nCheck that config file only contains " +
                  "one set of deformation inputs. \nAttempted: " + 
                  split_string[0] + "\nbut already had: " + 
                  inputs["deformation_input"]+".")
        
        inputs["deformation_input"] = split_string[0]
        
        # parse input
        eps_to_test = range_to_list(split_string[1])
        gamma_to_test = range_to_list(split_string[2])
        
        if len(eps_to_test) > 1:
            inputs["step"] = get_range_step(split_string[1])
            
        if len(gamma_to_test) > 1:
            inputs["step"] = get_range_step(split_string[2])
            
        if (len(eps_to_test) > 1) and (len(gamma_to_test) > 1): 
            raise ValueError("cannot vary both eps and gamma linearly; choose one or use mesh")
        
        
        # create arrays with repeated values, to cover every data point.
        # e.g. if eps_to_test = [0.01, 0.02, 0.03] and gamma_to_test = [15]
        # then eps_points = [0.01, 0.02, 0.03] and gamma_points = [15, 15, 15]
        # because all three eps values will be tested at gamma = 15º.
        eps_points = []
        gamma_points = []
        for e in range(len(eps_to_test)):
            for g in range(len(gamma_to_test)):
                eps_points.append(eps_to_test[e])
                gamma_points.append(gamma_to_test[g])
                
    elif split_string[0]=="mesh":
        
        if "deformation_input" in inputs:
            raise ValueError("Deformation parameter sets have been input " +
                        "multiple times. \nCheck that config file only " +
                        "contains one set of deformation inputs." +
                        "\nAttempted: " + split_string[0] + 
                        "\nbut already had: " + inputs["deformation_input"]+".")
        
        inputs["deformation_input"] = split_string[0]
        
        values = split_string[1].split(",")
        inputs["eps_max"] = float(values[0])
        outer_points = int(values[1])                                           # the total number of data points generated = 0.5*outer_points*(outer_points + 1)
        
        # arrange a mesh of evenly distributed (eps,gamma) points to test over, 
        # such that eps_max is tested at outer_points gamma values. 
        eps_to_test = np.linspace(0.001, inputs["eps_max"], num=outer_points)   # start from eps=0.001 because eps=0 is not an allowed input when it comes to asyrmo
                                   
        eps_points = []
        gamma_points = []
        for e in range(len(eps_to_test)):
            for g in range(e+1):
                eps = np.round(eps_to_test[e], 3)
                eps_points.append(eps)
                if e==0: gamma = 0
                else: gamma = np.round(g*(60/e), 3) 
                gamma_points.append(gamma)
        eps_points = np.array(eps_points)
        gamma_points = np.array(gamma_points)
        
        del [eps, gamma, outer_points, values]
        
        
    elif split_string[0]=="e2plus":
        
        e2plus_to_test = range_to_list(split_string[1])
        
        if (len(e2plus_to_test)>1 and not inputs["deformation_input"] == "single"):
            raise ValueError("testing a range of e2plus values is " +
                      "only supported \nfor a single deformation input. " +
                      '''Check that deformation is input as "single".''')
    
    elif split_string[0][:3]=="jp_":
        inputs[split_string[0]] = [float(n) for n in split_string[1].split(',')]
          
    elif split_string[0] in variable_types["int"]:
        inputs[split_string[0]] = int(split_string[1])
        
    elif split_string[0] in variable_types["float"]:
        inputs[split_string[0]] = float(split_string[1])
                                      
    elif split_string[0] in variable_types["bool"]:                                    
        inputs[split_string[0]] = bool(split_string[1])
        
    else: # string
        inputs[split_string[0]] = split_string[1]

    
    print(split_string[0] + ": \t" + split_string[1])
    


# check that there are no eps=0 values to test - this will cause the code to hang.
if 0.0 in eps_to_test:
    raise ValueError("Cannot test at eps=0.0 because the code will only work" +
                " for non-spherical deformations. \nCheck deformation inputs.")
    
data_points = {}
data_points["eps"] = eps_points
data_points["gamma_degrees"] = gamma_points
data_points["gamma_radians"] = [n*np.pi/180 for n in gamma_points]
del [eps_points, gamma_points]


#print(inputs)
print("********** Finished reading lines: "+str(line_count)+ " **********\n")

# deallocate variables that are no longer needed
del [e,g,n, line_count, expected_num_words] 
del [config_lines, line, line_string, split_string, word]


#%%
''' PROCESS INPUTS '''
# conver input fractional spins to floats and save separately
# using input A and Z, work out which is odd, to determine the nneupr input
# half and ceiling for the overall index of the fermi level orbital
# generate the orbitals input for gampn (e.g. orbitals_string = "+4 19 20 21 22")


# add float versions of the spin parameters
if "gs_spin" in inputs:
    inputs["gs_spin_float"] = spin_string_to_float(inputs["gs_spin"])

if "x1_spin" in inputs:
    inputs["x1_spin_float"] = spin_string_to_float(inputs["x1_spin"])

if "x2_spin" in inputs:
    inputs["x2_spin_float"] = spin_string_to_float(inputs["x2_spin"])

if "x3_spin" in inputs:
    inputs["x3_spin_float"] = spin_string_to_float(inputs["x3_spin"])


# convert the nantj, noutj, ipout inputs to the correct format
inputs["nantj"] = inputs["nantj"].replace(",", " ")
inputs["noutj"] = inputs["noutj"].replace(",", " ")
inputs["ipout"] = inputs["ipout"].replace(",", " ")

# determine nneupr and calculate fermi level
inputs["N"] = inputs["A"]-inputs["Z"]
    
if inputs["Z"]%2 == 0: # Z is even so N is odd so NNEUPR = -1
    inputs["nneupr"] = "-1" 
    inputs["fermi_level"] = math.ceil(inputs["N"]/2)
    print("calculating for odd neutrons...")                 
elif inputs["N"]%2 == 0:
    inputs["nneupr"] = "1"
    inputs["fermi_level"] = math.ceil(inputs["Z"]/2)
    print("calculating for odd protons...")
else:
    raise ValueError('''Input nucleus is even-even, 
                     but this code only applies to odd-mass nuclei.''')


# generate a string that contains the number of orbitals and their indices, in the correct format for input to gampn
gampn_orbitals = inputs["par"]
gampn_orbitals += str(inputs["num_orbs"])                                       # e.g. orbitals_string = "4 19 20 21 22"

first_index = inputs["fermi_level"]//2 - inputs["num_orbs"]//2                  # select [num_orbs] orbitals (the fermi level plus [num_orbs//2] either side)
last_index = inputs["fermi_level"]//2 + inputs["num_orbs"]//2
if inputs["num_orbs"]%2 == 1:                                                   # if an even number is requested, we need to add one to the final index to ensure the correct number are included
    last_index += 1
orbitals = np.r_[first_index:last_index]                                        # generate a list of orbitals in unit steps inside this range with list slicing


for i in orbitals:
    gampn_orbitals += " "
    gampn_orbitals += str(i)
    
inputs["gampn_orbitals"] = gampn_orbitals

del [i, first_index, last_index, gampn_orbitals, orbitals]

    
#%%   
''' WRITE GAMPN.DAT FILES USING DICT '''
# loops through the data points to be tested (specified by the arrays eps_points and gamma_points, which should have the same length) 
# each iteration:
#   creates a tag referencing the nucleus being tested, and the values of eps and gamma being tested in this iteration
#   creates/overwrites a .DAT file in the Inputs directory folder, named using the tag
#   uses string formatting with the dictionary to write the .DAT file 
#   (currently uses the provided example as a base, and changes as few settings as possible!)
# if anything went wrong, the existing .DAT file won't be overwritten 

print("writing GAMPN.DAT files...")

file_tags = []                                                          # create a list to record which files were written

for p in range(len(data_points["eps"])):
    for etp in e2plus_to_test:
         
        inputs["current_e2plus"] = etp
        inputs["current_eps"] = data_points["eps"][p]                                   # store eps in the dict for ease of use when string formatting below
        inputs["current_gamma"] = data_points["gamma_degrees"][p]
        
        file_tag = "e%.3f_g%.1f_p%.3f_%s" % (inputs["current_eps"], 
                                       inputs["current_gamma"],
                                       inputs["current_e2plus"],
                                       inputs["nucleus"])
        inputs["current_f002"] = "f002_"+file_tag+".dat"
        inputs["current_f016"] = "f016_"+file_tag+".dat"
        inputs["current_f017"] = "f017_"+file_tag+".dat"
        
        gampn_dat_text = templates["gampn"] % inputs
               
        with open(("../Inputs/GAM_"+file_tag+".DAT"), 'w') as f:
            f.write(gampn_dat_text)
        del f

        file_tags.append(file_tag)

print('''%d input files were written, 
      for eps in range [%.3f, %.3f],
      and gamma in range [%.1f, %.1f].'''
      % (len(file_tags), data_points["eps"][0], data_points["eps"][-1], data_points["gamma_degrees"][0], data_points["gamma_degrees"][-1]))

data_points["file_tags"] = file_tags
del [etp, p, file_tag, file_tags, gampn_dat_text]



#%%
''' WRITE AND RUN BASH SCRIPT TO EXECUTE GAMPN '''
# writes a bash shell script to run each of the GAMPN .DAT files written above
# after each one is run, the output file GAMPN.out is copied to a new text file with a more descriptive name (based on file tag),
# so that the data is not lost when the next iteration runs gampn again and GAMPN.out is overwritten

timer_data["lapse"] = time.time()
batch_settings = {}

batch_settings["allowed_time"] = 0.1*len(data_points["file_tags"]) + 10         #!!! each file takes ~ 0.06 seconds to run, as a rough average, so allow 0.1 seconds per file to be safe, with an overhead of 0.2                                                               # time in seconds to allow for running the bash script before timing out (assuming hanging code)

batch_settings["num_batches"] = 8                                               # = number of cores for maximum efficiency with large data sets
batch_settings["num_per_batch"] = math.ceil(len(data_points["file_tags"])/
                                            batch_settings["num_batches"])
if batch_settings["num_per_batch"] < 20:
    batch_settings["num_per_batch"] = 20
    batch_settings["num_batches"] = math.ceil(len(e2plus_to_test)*
                                              len(data_points["eps"])/
                                              batch_settings["num_per_batch"])             # if the data set is small then use fewer cores for a minimum batch size of 20 to make the overhead worth it

run_program = configure_script_writer(data_points["file_tags"], 
                                      batch_settings["num_batches"], 
                                      batch_settings["num_per_batch"], 
                                      batch_settings["allowed_time"])

run_program("gampn")

timer_data["lapse_end"] = time.time()
print("finished running gampn in time = %.2f seconds" % (timer_data["lapse_end"]-timer_data["lapse"]))
timer_data["lapse"] = timer_data["lapse_end"]




#%%
''' READ GAMPN.OUT FILE '''
# for each deformation (i.e. each GAMPN.OUT file):
#   determine the line number of the fermi level orbital in the GAMPN.OUT file
#   read that line from GAMPN.OUT and get the energy, parity, and single-party-index
#   generate the orbitals input for asyrmo (e.g. orbitals_string = "+11 19 20 21 22 23 24 25 26 27 28 29")

print("\nreading GAMPN.OUT files...")

asyrmo_orbitals    = []                                                         # an empty array to store inputs for asyrmo.dat

fermi_energies_hw  = []                                                         # for the fermi energy at each deformation
fermi_energies_mev = []
fermi_indices      = []                                                         # for the index and parity of the fermi level
fermi_parities     = [] 

for file in data_points["file_tags"] :

    with open(("GAM_"+file+".OUT"), 'r') as f:
        lines = f.readlines()
    del f
    
    # get the EFAC value (conversion factor from hw to MeV)
    efac_line = lines.index("     KAPPA    MY     EPS   GAMMA    EPS4     EPS6     W0/W00   NMAX  COUPL     OMROT      EFAC      QFAC\n")
    inputs["efac"] = float(lines[efac_line+1][85:95].strip())
    
    # locate the fermi level
    levels_header_line = lines.index("   #   ENERGY +/-(#)    <Q20>    <Q22>     <R2>     <JZ>      #   ENERGY +/-(#)    <Q20>    <Q22>     <R2>     <JZ>\n")
    fermi_level_line = inputs["fermi_level"]+levels_header_line+1               # calculate the line number of the fermi level in the GAMPN.OUT file (indexed from zero!)
    
    if inputs["fermi_level"] > 40:
        fermi_level_line -= 40
        whole_line = lines[fermi_level_line]
        half_line = whole_line[60:-1].strip()                                   # get only the second half of the line
    else:
        whole_line = lines[fermi_level_line]
        half_line = whole_line[0:60].strip()                                    # get only the first half of the line
        
    # get data from the half-line containing data on the fermi level orbital
    hash_index = half_line.index("#")                                           # use the index of '#' as a reference point 
    fermi_parities.append(half_line[hash_index-2])                              
    fermi_energies_hw.append(float(half_line[hash_index-10 : hash_index-4]))   
    fermi_energies_mev.append(fermi_energies_hw[-1]*inputs["efac"])

    single_parity_index = half_line[hash_index+1 : hash_index+3]
    if single_parity_index[1] == ")":                                           # in case the index is only a single digit, ignore the ")" that will have been caught
        single_parity_index = single_parity_index[0]
    fermi_indices.append(int(single_parity_index))
    
    # generate the orbitals input for asyrmo
    first_index = int(single_parity_index) - inputs["nu"]//2
    last_index = int(single_parity_index) + inputs["nu"]//2
    if inputs["nu"]%2 == 1:                                                     # then we need to add one to the final index to ensure the correct number are included
        last_index += 1
    orbitals = np.r_[first_index:last_index]
    
    orbitals_string = fermi_parities[-1] + str(inputs["nu"])                    # e.g. orbitals_string = "+11 19 20 21 22 23 24 25 26 27 28 29"
    for l in orbitals:
        orbitals_string += " "
        orbitals_string += str(l)
    asyrmo_orbitals.append(orbitals_string)


output_data = {}

output_data["fermi_parities"] = fermi_parities
output_data["fermi_energies_hw"] = fermi_energies_hw
output_data["fermi_energies_mev"] = fermi_energies_mev
output_data["fermi_indices"] = fermi_indices

data_points["asyrmo_orbitals"] = asyrmo_orbitals
print("finished reading %d files\n" % len(asyrmo_orbitals))

del [efac_line, levels_header_line, fermi_level_line, first_index, last_index]
del [file, half_line, hash_index, l]
del [lines, orbitals_string, single_parity_index, whole_line, orbitals]
del [fermi_energies_hw, fermi_energies_mev, fermi_indices, fermi_parities, asyrmo_orbitals]

_output_data = output_data # save a copy of the original before it's overwritten

#%%
''' WRITING ASYRMO.DAT FILE '''
# for each deformation, writes a .DAT file for ASYRMO using the file tag to name, and the provided example as a base


print("writing ASYRMO.DAT files...")


for t in range(len(data_points["file_tags"])):

    inputs["current_orbitals"] = data_points["asyrmo_orbitals"][t]
    
    if len(e2plus_to_test)>1:
        inputs["current_e2plus"] = e2plus_to_test[t]
    else:
        inputs["current_e2plus"] = e2plus_to_test[0]

    inputs["current_f016"] = "f016_"+data_points["file_tags"][t]+".dat"
    inputs["current_f017"] = "f017_"+data_points["file_tags"][t]+".dat"
    inputs["current_f018"] = "f018_"+data_points["file_tags"][t]+".dat"
    
    asyrmo_dat_text = templates["asyrmo"] % inputs
   
    with open(("../Inputs/ASY_"+data_points["file_tags"][t]+".DAT"), 'w') as f:
        f.write(asyrmo_dat_text)
    del f
    
print("finished writing %d ASY.DAT files" % (t+1))
del [t, asyrmo_dat_text]



#%%
''' WRITE AND RUN BASH SCRIPT TO EXECUTE ASYRMO '''
# writes a bash shell script to run each of the ASYRMO .DAT files written above
# after each one is run, the output file ASYRMO.out is copied to a new text file with a more descriptive name (based on file tag),
# so that the data is not lost when the next iteration runs asyrmo again and ASYRMO.out is overwritten

timer_data["lapse"] = time.time()

run_program("asyrmo")

timer_data["lapse_end"] = time.time()
print("finished running asyrmo in time = %.2f seconds" % (timer_data["lapse_end"]-timer_data["lapse"]))
timer_data["lapse"] = timer_data["lapse_end"]




#%%
''' WRITING PROBAMO.DAT FILE '''
# for each deformation, writes a .DAT file for PROBAMO, using the file tag to name, and the provided example as a base

for file in data_points["file_tags"]:
    
    inputs["current_f017"] = "f017_"+file+".dat"
    inputs["current_f018"] = "f018_"+file+".dat"

    probamo_dat_text = templates["probamo"] % inputs
    
    with open(("../Inputs/PROB_"+file+".DAT"), 'w') as f:
        f.write(probamo_dat_text)    
    del f


print("finished writing PROB.DAT files")

del [probamo_dat_text, file]


#%%
''' WRITE AND RUN BASH SCRIPT TO EXECUTE PROBAMO '''
# writes a bash shell script to run each of the PROBAMO .DAT files written above
# after each one is run, the output file PROBAMO.out is copied to a new text file with a more descriptive name (based on file tag),
# so that the data is not lost when the next iteration runs asyrmo again and PROBAMO.out is overwritten


timer_data["lapse"] = time.time()

run_program("probamo")

timer_data["lapse_end"] = time.time()
print("finished running probamo in time = %.2f seconds" 
      % (timer_data["lapse_end"]-timer_data["lapse"]))
timer_data["lapse"] = timer_data["lapse_end"]

#%%

''' READ PROBAMO.OUT FILE '''
# for each deformation (i.e. each PROBAMO.OUT file):
# the row containing the ground state is identified by energy=0.0
# the ground state spin and magnetic dipole moment are recorded
# the lowest energy (yrast) state with the same spin as the experimental first excited state is identified
# the "first excitation energy" is recorded
# similarly for the second and third experimental excited states


''' FUNCTIONS FOR READING PROBAMO '''

def read_data(line):
    line = line.strip()
    
    try:                                                                        # determine whether this line is a data row of the table (if ' - ' is present then it is)
        dash_index = line.index(" - ")
    except ValueError:
        return False # continue to next line
    
    spin_string = line[dash_index-4:dash_index].strip()
    final_spin_string = line[dash_index+11:dash_index+16].strip()               # get the spin of the final state (after transition)
    
    if not(spin_string == final_spin_string):
        return False

    spin_float = spin_string_to_float(spin_string)
    
    try:
        this_energy = float(line[:6].strip())
    except ValueError: # could not convert string to float: '0.0  1'
        this_energy = float(line[:5].strip())
    
    final_energy = float(line[dash_index+3:dash_index+10].strip())
                         
    if not(this_energy == final_energy):
        return False
    
    mag_moment = float(line[-8:].strip())
    line_data = {'spin_string': spin_string, 'spin_float':spin_float, 
                 'energy':this_energy, 'mag_moment':mag_moment}
    
    return line_data

def sort_by_spin(line_data, file_data):
    
    spin = "spin_"+line_data["spin_string"]
    
    if (spin+"_energies") in file_data:
        file_data[(spin+"_energies")].append(line_data["energy"])
        file_data[(spin+"_mag_moments")].append(line_data["mag_moment"])
        
    else:
        file_data[(spin+"_energies")] = [line_data["energy"]]
        file_data[(spin+"_mag_moments")] = [line_data["mag_moment"]]
    
    return file_data

def sort_by_expectation(line_data, file_data):
    
    if line_data["energy"] == 0.0: # record the ground state separately
        file_data["gs_spin_strings"] = line_data["spin_string"]
        file_data["gs_spin_floats"] = line_data["spin_float"]
        file_data["gs_mag_moments"] = line_data["mag_moment"]
    
    if "x1_spin" in inputs:
        if line_data["spin_string"] == inputs["x1_spin"]:
            file_data["x1_energies"] = line_data["energy"]
            file_data["x1_mag_moments"] = line_data["mag_moment"]
    
    if "x2_spin" in inputs:
        if line_data["spin_string"] == inputs["x2_spin"]:
            file_data["x2_energies"] = line_data["energy"]
            file_data["x2_mag_moments"] = line_data["mag_moment"]

    if "x3_spin" in inputs:
        if line_data["spin_string"] == inputs["x3_spin"]:
            file_data["x3_energies"] = line_data["energy"]
            file_data["x3_mag_moments"] = line_data["mag_moment"]
    
    return file_data
            
def missing_data(file_data): 
    
    for i in range(int(inputs["ispin"])+1):
        
        if i%2 == 0:
            continue # only half-int spins are calculated
        
        spin = "spin_"+str(i)+"/2"
        
        if not(spin+"_energies" in file_data):
            file_data[(spin+"_energies")] = [np.NaN]
            file_data[(spin+"_mag_moments")] = [np.NaN]

        
    return file_data
            
def restructure_data(old_data):
# Take input data structured as a list of dictionaries. 
#       Each dictionary represents one deformation data-point.
#       Therefore the number of dictionaries is equal to the number of data points.
#       Each dictionary contains lists of states, organised by spin and property (energy or magnetic moment).
#       Therefore the number of lists in each dictionary is (1 + ((ISPIN+1)/2)*2) for gs magnetic moment and two properties.
# 
# Outputs the same data, reorganised into a new dictionary of lists.
#       Each list represents one property (energy or magnetic moment) and spin.
#       Therefore the number of lists is (1 + ((ISPIN+1)/2)*2) for gs magnetic moment and two properties.
#       Each list contains sub-lists of states, organised by deformation.
#       Therefore all the lists have the same length (equal to the number of data points).
#       This is more useful for plotting graphs.
    
    new_data = {}
    new_data["gs_spin_strings"] = []
    new_data["gs_spin_floats"] = []
    new_data["gs_mag_moments"] = []
    new_data["x1_mag_moments"] = []
    new_data["x2_mag_moments"] = []
    new_data["x3_mag_moments"] = []
    new_data["x1_energies"] = []
    new_data["x2_energies"] = []
    new_data["x3_energies"] = []
    

    for d in range(len(old_data)):
        new_data["gs_spin_strings"].append(old_data[d]["gs_spin_strings"])
        new_data["gs_spin_floats"].append(old_data[d]["gs_spin_floats"])
        new_data["gs_mag_moments"].append(old_data[d]["gs_mag_moments"])
        
        try:
            new_data["x1_energies"].append(old_data[d]["x1_energies"])
            new_data["x1_mag_moments"].append(old_data[d]["x1_mag_moments"])
        except(KeyError):
            print("Could not find any states with first excited spin in file " + str(d))
            new_data["x1_energies"].append(np.NaN)
            new_data["x1_mag_moments"].append(np.NaN)
        
        try:
            new_data["x2_energies"].append(old_data[d]["x2_energies"])
            new_data["x2_mag_moments"].append(old_data[d]["x2_mag_moments"])
        except(KeyError):
            print("Could not find any states with second excited spin in file " + str(d))
            new_data["x2_energies"].append(np.NaN)
            new_data["x2_mag_moments"].append(np.NaN)
        
        try:
            new_data["x3_mag_moments"].append(old_data[d]["x3_mag_moments"])
            new_data["x3_energies"].append(old_data[d]["x3_energies"])
        except(KeyError):
            print("Could not find any states with third excited spin in file " + str(d))
            new_data["x3_energies"].append(np.NaN)
            new_data["x3_mag_moments"].append(np.NaN)
    
    for i in range(int(inputs["ispin"])+1): 
        
        if i%2 == 0:
            continue # only half-int spins are calculated
    
        spin = "spin_"+str(i)+"/2"
        
        new_data[spin+"_energies"] = []
        new_data[spin+"_mag_moments"] = []
        
        for d in range(len(old_data)):
           
            new_data[spin+"_energies"].append(old_data[d][spin+"_energies"])
            new_data[spin+"_mag_moments"].append(old_data[d][spin+"_mag_moments"])
            
    return new_data
    
print("\nreading PROBAMO.OUT files...")

def fill_gaps(multi_level_data):
    
    num_levels = [len(n) for n in multi_level_data]
    max_num = max(num_levels)
    
    for d in range(len(multi_level_data)):
        while len(multi_level_data[d]) < max_num:
            multi_level_data[d].append(np.NaN)
        
    return multi_level_data
    





# start reading files 

data_points["property_data"] = []

for t in range(len(data_points["file_tags"])):
    file = data_points["file_tags"][t]
    file_data = {}
    
    with open(("PROB_"+file+".OUT"), 'r') as f:
        lines = f.readlines()
    del f

    for line in lines:
        
        line_data = read_data(line)  # get the spin, energy, and magnetic moment from this line if it is a static moment, else return False
        
        if not(line_data):
            continue # to next line in file
        
        # sort line_data into file_data according to its spin
        file_data = sort_by_spin(line_data, file_data)
        
        # separately save data that corresponds to the experimental ground state and input excited states
        file_data = sort_by_expectation(line_data, file_data)
    
    file_data = missing_data(file_data)
    data_points["property_data"].append(file_data)

    
output_data = _output_data | restructure_data(data_points["property_data"])

transposed_data = {}
for d in output_data:
    if isinstance(output_data[d][0], list):
        output_data[d] = fill_gaps(output_data[d])
        transposed_data[d] = np.transpose(np.array(output_data[d]))
    
del [file, file_data, line, line_data, lines, t, d]

_output_data_dict = output_data # save a copy of the original before it's overwritten

#%%
''' PREPARE TO PLOT GRAPHS '''
# all arrays of data to be plotted are collected into a single matrix.
# lists of the titles to be used in each plot are constructed.
# similarly for contour levels, colour bar ticks, data axis labels and units, 
#       experimental data to compare results to, absolute error tolerance, 
#       and whether the data set is continuous or discrete.

''' FUNCTIONS FOR ORGANISING DATA '''

def calc_contour_levels(data):
    
    min_contour = min(data)-0.5
    max_contour = max(data)+1.5
    return np.arange(min_contour, max_contour, 1.0)

def calc_cbar_tick_labels(data):
    cbar_ticks = np.arange(min(data), max(data)+1.0, 1.0)
    return [spin_float_to_string(n) for n in cbar_ticks]

def calc_cbar_ticks(contours):
    num = len(contours) - 1
    return [(contours[i] + contours[i+1]) / 2 for i in range(num)]                                # Calculate midpoints of levels for tick placement

def try_experimental(inputs, key, tolerance):
    try:
        return (inputs[key], tolerance)
    except KeyError:
        return (np.NaN, np.NaN)

class PropertyData:
    
    plot = False # don't plot by default
    
    def __init__(self, data): 
        
        self.data = data
        
def sort_property(name):
    
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
    
    return num, prop, sort
    





# convert output_data from a dictionary of lists to a dictionary of PropertyData objects 
output_data = {}
for p in _output_data_dict:
    
    output_data[p] = PropertyData(_output_data_dict[p])
    
    # work out what kind of property it is from its name
    num, prop, sort = sort_property(p)
    
    
    # determine graph plotting attributes
    if prop == "energies" and not(sort=="Fermi"): 
        output_data[p].title = "Energies of "+ sort + num + " States"
        output_data[p].contour_levels = 10
        output_data[p].cbar_tick_labels = 0
        output_data[p].axis_label = output_data[p].title + " / keV"
        
        
        if sort == "Excited State ":
            
            output_data[p].experimental_data, output_data[p].error_tolerance = try_experimental(inputs, "x" + num + "_energy", 50)
        
        elif sort == "Spin ":
            
            try:
                output_data[p].experimental_data = inputs["jp_"+num+inputs["par"]]
                output_data[p].error_tolerance = 50
            except KeyError:
                output_data[p].experimental_data = np.NaN
                output_data[p].error_tolerance = np.NaN
            
            
        else: raise ValueError("property not recognised: " + p)
        
    elif prop == "mag_moments":
            
        if sort == "Excited State ":
            output_data[p].title = "Magnetic Dipole Moment of Excited State " + num 
            output_data[p].experimental_data, output_data[p].error_tolerance = try_experimental(inputs, "x" + num + "_mu", 0.2)
    
        elif sort == "Spin ":
            output_data[p].title = "Magnetic Dipole Moment of Spin " + num + " States"
            output_data[p].experimental_data = np.NaN
            output_data[p].error_tolerance = np.NaN
            
        elif sort == "Ground":
            output_data[p].title = "Ground State Magnetic Dipole Moment"
            output_data[p].experimental_data, output_data[p].error_tolerance = try_experimental(inputs, "gs_mu", 0.2)
            
        else: raise ValueError("property not recognised: " + p)
        
        output_data[p].contour_levels = 8
        output_data[p].axis_label = output_data[p].title + r' / $μ_{N}$'
        output_data[p].cbar_tick_labels = 0
    
    elif sort == "Ground":
        
        if prop == "spin_floats": 
        
            output_data[p].title = "Ground State Spin"
            output_data[p].contour_levels = calc_contour_levels(output_data[p].data)
            output_data[p].axis_label = output_data[p].title
            
            output_data[p].experimental_data, output_data[p].error_tolerance = try_experimental(inputs, "gs_spin_float", 0.2)
            output_data[p].cbar_tick_labels = calc_cbar_tick_labels(output_data[p].data)
            output_data[p].cbar_ticks = calc_cbar_ticks(output_data[p].contour_levels)
            
        elif prop == "spin_strings": continue
    
    elif sort == "Fermi":
    
        if prop == "indices": 
            output_data[p].title = "Fermi Level Index"
            output_data[p].contour_levels = calc_contour_levels(output_data[p].data)
            
            output_data[p].axis_label = output_data[p].title
            output_data[p].cbar_tick_labels = calc_cbar_tick_labels(output_data[p].data)
            output_data[p].cbar_ticks = calc_cbar_ticks(output_data[p].contour_levels)
            
        
        elif prop == "energies":
        
            output_data[p].title = "Fermi Energy"
            output_data[p].contour_levels = 10
            output_data[p].axis_label = output_data[p].title + p[-3:-1]
            output_data[p].cbar_tick_labels = 0
        
        elif prop == "parities":
        
            output_data[p].title = "Fermi Parity"
            output_data[p].contour_levels = 2
            output_data[p].axis_label = output_data[p].title
            output_data[p].cbar_tick_labels = 0
        
        else: raise ValueError("property not recognised: " + p)
            
        output_data[p].experimental_data = np.NaN
        output_data[p].error_tolerance = np.NaN
    
    else: raise ValueError("property not recognised: " + p)
    
    output_data[p].sort = sort
    output_data[p].prop = prop
    output_data[p].num = num
    
    
del [sort, prop, p, num]


#%%


''' PLOT GRAPHS'''
# plots variation of each data set with distortion
# if a mesh of eps and gamma has been tested, plot wedge contour plots
#   for each data point in each data set:
#       use '+' data markers for positive parity points, and '.' markers for negative parity points
#       compare with experimental value (if available)
#       if the absolute error is below a tolerance, consider this to be a "match"; plot the data marker in red
#       otherwise plot in white
#   keep a list of data points at which all the data sets match their corresponding experimental values (within tolerance) 
# if only one of eps and gamma is varying, plot line graphs 
#   mark the experimental value in red for easy comparison
#   mark the region of correct ground state spin in green for easy identification of the relevant (and meaningful) results

''' FUNCTIONS FOR GRAPH PLOTTING '''

def format_fig(polar_or_linear, ax, legend_handles, title, **kwargs):
    
    if polar_or_linear == 'polar':
        ax.set_thetamin(0)   
        ax.set_thetamax(60)  
        
        theta_ticks = np.arange(0, 70, 10)  
        ax.set_xticks(np.radians(theta_ticks))
        ax.yaxis.set_major_formatter(FuncFormatter(lambda x, _: f'{x:.2f}'))    # set the number of decimal places 
        
        plt.xlabel("ε")
        ax.text(45*np.pi/180, ax.get_rmax()*1.2, "γ", ha='center', va='center') # gamma axis label
        
        ax.legend(handles=legend_handles, loc="upper left", 
                           facecolor = 'lightblue', bbox_to_anchor=(-0.5, 1.1))
        
        ax.set_title(title, va='bottom', y=1.1)  
        
    elif polar_or_linear == 'linear':
        
        varied = kwargs.get("varied", None)
        x_label = kwargs.get("x_label", None)
        y_label = kwargs.get("y_label", None)
        legend_title = kwargs.get("legend_title", "")
    
        pad = 0.05*(varied[-1]-varied[0])
        ax.set_xlim([varied[0]-pad, varied[-1]+pad]) 
        
        plt.xlabel(x_label)
        plt.ylabel(y_label)
        
        ax.set_title(title, va='bottom', y=1.1) 
        
        legend = ax.legend(handles = legend_handles)
        
        if legend_title != "x":
            legend.set_title(legend_title)
    
    else: raise ValueError("unrecognised graph type: " + polar_or_linear + "; must be either 'polar' or 'linear'.")

def draw_contour_plot(ax, prop, data_points):
    
    if prop.sort == "Spin ":
        
        plot_data = np.transpose(prop.data)[0]
        
        print("multiple levels can't be plotted on a contour plot; \nplotting only the yrast state of spin " + prop.num)
        
        
    else:
        plot_data = np.array(prop.data)
        
    # some of the data points may be NaN, so manually create the triangulation and mask NaN triangles
    triang = tri.Triangulation(data_points["gamma_radians"], data_points["eps"])
    mask = np.any(np.isnan(plot_data[triang.triangles]), axis=1)
    triang.set_mask(mask)
    
    cax = ax.tricontourf(triang, plot_data, levels=prop.contour_levels)
    cbar = plt.colorbar(cax, pad=0.1, label=prop.axis_label)


    
    if prop.cbar_tick_labels:                                   # then format for discrete values
        
        cbar.set_ticks(prop.cbar_ticks)
        cbar.set_ticklabels(prop.cbar_tick_labels)
    
        
    return cax, cbar

def find_correct_spin(gs_spins, experimental_value):
    correct_spin_range = []                                             
    start_flag = False
    
    for i in range(len(gs_spins)):
        if gs_spins[i] == experimental_value and not start_flag:  # correct, and range hasn't started yet
            start_flag = True
            correct_spin_range.append(i)
            
        elif gs_spins[i] != experimental_value and start_flag:    # incorrect, and range has started
            start_flag = False
            correct_spin_range.append(i)
    
    return correct_spin_range

def plot_correct_spin(correct_spin_range, var, step, prop):
    
    for r in range(len(correct_spin_range)):
        if correct_spin_range[r] == 0:                                  # the first value of eps in the range has the correct spin
            start_range = (var[correct_spin_range[r]]- step/2)          
        else:
            start_range = np.mean([var[correct_spin_range[r]-1], 
                                   var[correct_spin_range[r]]])
        
        if prop.sort == "Spin ":
            data = prop.data[0]
        else:
            data = prop.data
            
        correct_spin, = plt.plot([start_range, start_range], 
                                [min(data)+0.05*max(data), max(data)*1.05], 
                                'g-', label="range of correct spin")   # this plots the front and end edges of the box
        if r%2==0:
            if r+1 == len(correct_spin_range):                          # the last value of eps in the range has the correct spin
                end_range = (var[-1]+step/2)
            else:
                end_range = np.mean([var[correct_spin_range[r+1]-1], 
                                     var[correct_spin_range[r+1]]])
            
            plt.plot([start_range, end_range],                          # this plots the bottom edge of the box
                     [min(data)+0.05*max(data), 
                      min(data)+0.05*max(data)], 'g-')
            plt.plot([start_range, end_range], 
                     [max(data)*1.05, 
                      max(data)*1.05], 'g-')                  # this plots the top edge of the box
    return correct_spin

def plot_points_with_experiment(data_points, prop, legend_handles, cbar):
    
    legend_hit = False
    legend_miss = False
    
    for r in range(len(data_points["eps"])):
        
        # use a smaller marker size for large data sets
        if len(data_points["file_tags"]) < 100: marker_size = 20
        else: marker_size = 5
        
        if prop.sort == "Spin ":
            error = [abs(prop.data[r][0] - exp) for exp in prop.experimental_data]
            match = [err < prop.error_tolerance for err in error]
        else:
            error = abs(prop.data[r] - prop.experimental_data)
            match = [error < prop.error_tolerance]
            
        if any(match): 
            data_points["agreed"][r] += 1                                       # record how many of the tested properties agree
            hit = plt.scatter(data_points["gamma_radians"][r], 
                      data_points["eps"][r], s=marker_size, edgecolor='red', 
                      facecolor='None', label="matches experiment")
            legend_hit = True
            
        # for data points that don't match experimental data, only plot them when the data set is quite small, to avoid cluttering the graph 
        elif len(data_points["file_tags"]) < 100: 
            miss, = plt.polar(data_points["gamma_radians"][r], 
                      data_points["eps"][r], 'wx', label="does not match experiment")
            legend_miss = True
    
    if legend_hit:
        legend_handles.append(hit)
    if legend_miss:
        legend_handles.append(miss)
    
    if prop.sort == "Spin ":
        for e in range(len(prop.experimental_data)):
            exp = cbar.ax.plot([0, 1], 
                               [prop.experimental_data[e], prop.experimental_data[e]], 
                               'r-', label = "experimental value")
    else:
        exp = cbar.ax.plot([0, 1], [prop.experimental_data, prop.experimental_data], 
                           'r-', label = "experimental value")
        
   
    legend_handles.append(exp[0])
    
    return legend_handles

def plot_points_without_experiment(data_points, legend_handles):
    
    # use a smaller marker size for large data sets
    if len(data_points["file_tags"]) < 100: marker_size = 5
    else: marker_size = 1
    
    all_points = plt.scatter(data_points["gamma_radians"], 
              data_points["eps"], s=marker_size, c='w', label="data point")
    legend_handles.append(all_points) 
    
    return legend_handles

def assign_parameters(inputs, e2plus_to_test, data_points):

    if inputs["deformation_input"] == "eps" :
        
        var_sym = "ε"
        var = data_points["eps"]
        fix_sym = "γ"
        fix = data_points["gamma_degrees"]
        
    elif inputs["deformation_input"] == "gamma" :
        
        var_sym = "γ / º"
        var = data_points["gamma_degrees"]
        fix_sym = "ε"
        fix = data_points["eps"]
        
    elif len(e2plus_to_test) > 1:
        
        var_sym = "E2PLUS / MeV"
        var = e2plus_to_test
        fix_sym = "(ε, γ)"
        fix = []
        
    else: raise ValueError("unrecognised graph request")

    return (var_sym, var, fix_sym, fix)

def plot_multi_lines(prop, var, legend_handles):
    data_by_line = np.transpose(prop.data)
    line_colours = ['k-x', 'b-x', 'y-x']
    line_labels = ["lowest energy", "second lowest energy", "third lowest energy"]
    
    for s in range(min(len(line_labels), np.size(data_by_line,0))):
        
        data, = plt.plot(var, data_by_line[s], line_colours[s], label=line_labels[s])
        legend_handles.append(data)
        
    legend_title = "%s = %s" % (fix_sym, fix[0])
    
    return legend_handles, legend_title

def mark_spin(inputs, data_points, output_data, legend_handles):

    correct_range = [inputs["gs_spin_float"]-0.5, 
                     inputs["gs_spin_float"]+0.5]
    spin_colour = (0,0,0) #(213/255,1,0)
    
    ax.tricontour(data_points["gamma_radians"], data_points["eps"], 
                  output_data["gs_spin_floats"].data, levels=correct_range,  
                  colors=[spin_colour], linewidths=1.0)
    spin_legend_proxy = plt.Line2D([], [], color=spin_colour, linewidth=1.0, label="region of correct g.s. spin") 
    legend_handles.append(spin_legend_proxy)

    return legend_handles






# set which graphs to plot:
output_data["gs_mag_moments"].plot = True
output_data["gs_spin_floats"].plot = True
output_data["spin_1/2_energies"].plot = True
output_data["spin_3/2_energies"].plot = True
output_data["spin_5/2_energies"].plot = True
output_data["x1_energies"].plot = True

data_points["agreed"] = [0]*len(data_points["eps"])


for g in output_data:
    
    prop = output_data[g]
    if not(prop.plot):
        continue
    inputs["current_graph"] = prop.title # makes several inputs more efficient
    print("plotting graph: %(current_graph)s" % inputs) 
    
    
    
    
    if inputs["deformation_input"] == "mesh":  
        
        fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
        cax, cbar = draw_contour_plot(ax, prop, data_points)
        
        
        legend_handles = []
        
        if inputs["mark_spin"]:
            
            legend_handles = mark_spin(inputs, data_points, output_data, legend_handles)
            
        # plot the data point markers, with comparison to experiment if possible 
        if np.isfinite(prop.experimental_data).all(): 
            
            legend_handles = plot_points_with_experiment(data_points, prop, legend_handles, cbar)
            
        else: 
            
            legend_handles =  plot_points_without_experiment(data_points, legend_handles)
               
        
        format_fig('polar', ax, legend_handles, '%(current_graph)s of %(nucleus)s' % inputs)
        
        plt.show()
        
        
       
    elif (inputs["deformation_input"] ==  "gamma" 
          or inputs["deformation_input"] == "eps"
          or len(e2plus_to_test) > 1):
        
        # set which paramters are varied and which are constant
        var_sym, var, fix_sym, fix = assign_parameters(inputs, e2plus_to_test, data_points)
        
        fig, ax = plt.subplots() 
        
        legend_handles = []
        
        # now plot the actual data
        if prop.sort == "Spin ":
        
            legend_handles, legend_title = plot_multi_lines(prop, var, legend_handles)
        
        else:
            data, = plt.plot(var, prop.data, 'k-x', label="%s = %s" % (fix_sym, fix[0]))
            legend_handles.append(data)
            legend_title = ""
          
        # if experimental data is available, plot it in red for easy comparison
        if np.isfinite(prop.experimental_data).all():  
            
            if prop.sort == "Spin ":
                for e in range(len(prop.experimental_data)):
                    exp, = plt.plot(var, np.full(len(var), prop.experimental_data[e]), 
                       'r-', label="experimental value")
            else:
                exp, = plt.plot(var, np.full(len(var), prop.experimental_data), 
                       'r-', label="experimental value")
                
            
            legend_handles.append(exp)
            
    
        # mark the range in which the correct ground state spin was calculated
        if inputs["mark_spin"]==1:
            
            correct_spin_range = find_correct_spin(output_data["gs_spin_floats"].data, inputs["gs_spin_float"])
            if len(correct_spin_range) > 0:
                spin = plot_correct_spin(correct_spin_range, var, inputs["step"], prop)
                legend_handles.append(spin)
            
                    
        format_fig('linear', ax, list(reversed(legend_handles)), 
                       '%(current_graph)s in %(nucleus)s' % inputs, 
                       varied=var, x_label=var_sym, y_label=prop.axis_label, 
                       legend_title=legend_title)
        
        plt.show()
        
        del [var, var_sym, fix, fix_sym, legend_title, legend_handles]
        
       
del [g]


print("\n\nAgreement of each data point with experimental data: ")
print(data_points["agreed"])

# note how long it took
timer_data["lapse_end"] = time.time()
print("finished plotting graphs in time = %.2f seconds" % (timer_data["lapse_end"]-timer_data["lapse"]))
print("total runtime = %.2f seconds" % (timer_data["lapse_end"]-timer_data["start"]))


