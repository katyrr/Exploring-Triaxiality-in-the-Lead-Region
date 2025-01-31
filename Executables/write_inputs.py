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

variable_types["float"] = ["gs_energy", "fx_energy", "sx_energy", "tx_energy",
                           "gs_mu", "fx_mu"]


#%%

''' DEFINE FUNCTIONS '''



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
        range_list =  np.arange(split_range[0],                          
                                 split_range[1]+split_range[2], 
                                 split_range[2])
    else : range_list = [split_range[0]]
    
    return range_list


def configure_script_writer(file_tags,  num_batches, 
                            num_per_batch, allowed_time):                       # these arguments will be the same throughout the execution, so we can define a family of script writers (one for each of gam/asy/prob)
    
    def run_script_batches(program):                                            # after configuring the script writer family, this is what we actually call (it's a "closure" function)
        
        if program == "gampn": abr = "GAM"
        elif program == "asyrmo": abr = "ASY"
        elif program == "probamo": abr = "PROB"
        else: raise ValueError("Input program ["+ program +"] not supported.")
        
        def write_script_batch(file_tag_batch, batch_index, file_path):         # a sub-function, called in the loop below
            
            script_text = ""
            
            for file in file_tag_batch :
                script_text += ("\n./../../../Executables/MO/" + program + 
                              " < ../Inputs/" + abr + "_" + file + ".DAT")
                script_text += ("\ncp " + program.upper() + ".out " + 
                              abr + "_" + file + ".OUT")                        # copy GAMPN.out to a new txt file with a more descriptive name
            
            #script_text += ("\nrm f002_*")                                     #!!! delete the f002 output, as it is empty and not used
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



#%%


''' OPEN AND READ CONFIG FILE '''
# ignores empty lines, and lines beginning with * (to mark a comment)
# counts and outputs the number of non-empty and non-comment lines (should be equal to the number of variables input in the config file)
# expects each line to have the format:    var_name value
# outputs a warning if the format is unexpected, otherwise splits the string at " " and assigns the variable to a dictionary
# tests whether deformation has been input as a range, mesh, or single value, and creates a list of values to test (which may contain just one value)

timer_start = time.time()

print("reading from config file...") 
with open("../Inputs/config.txt", 'r') as f:
    config_lines = f.readlines()

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

if "fx_spin" in inputs:
    inputs["fx_spin_float"] = spin_string_to_float(inputs["fx_spin"])

if "sx_spin" in inputs:
    inputs["sx_spin_float"] = spin_string_to_float(inputs["sx_spin"])

if "tx_spin" in inputs:
    inputs["tx_spin_float"] = spin_string_to_float(inputs["tx_spin"])


# convert the nantj, noutj, ipout inputs to the correct format
inputs["nantj"] = inputs["nantj"].replace(",", " ")
inputs["noutj"] = inputs["noutj"].replace(",", " ")
inputs["ipout"] = inputs["ipout"].replace(",", " ")

# determine nneupr and calculate fermi level
inputs["N"] = inputs["A"]-inputs["Z"]
    
if inputs["Z"]%2 == 0:
    inputs["nneupr"] = "1" 
    inputs["fermi_level"] = math.ceil(inputs["Z"]/2)
    print("calculating for odd protons...")                 
elif inputs["N"]%2 == 0:
    inputs["nneupr"] = "-1"
    inputs["fermi_level"] = math.ceil(inputs["N"]/2)
    print("calculating for odd neutrons...")
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

for p in range(len(gamma_points)):
    for etp in e2plus_to_test:
         
        inputs["current_e2plus"] = etp
        inputs["current_eps"] = eps_points[p]                                   # store eps in the dict for ease of use when string formatting below
        inputs["current_gamma"] = gamma_points[p]
        
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

        file_tags.append(file_tag)

print('''%d input files were written, 
      for eps in range [%.3f, %.3f],
      and gamma in range [%.1f, %.1f].'''
      % (len(file_tags), eps_points[0], eps_points[-1], gamma_points[0], gamma_points[-1]))

del [etp, p, file_tag, gampn_dat_text]



#%%
''' WRITE AND RUN BASH SCRIPT TO EXECUTE GAMPN '''
# writes a bash shell script to run each of the GAMPN .DAT files written above
# after each one is run, the output file GAMPN.out is copied to a new text file with a more descriptive name (based on file tag),
# so that the data is not lost when the next iteration runs gampn again and GAMPN.out is overwritten

timer_lapse = time.time()

allowed_time = 0.1*len(file_tags) + 10                                          #!!! each file takes ~ 0.06 seconds to run, as a rough average, so allow 0.1 seconds per file to be safe, with an overhead of 0.2                                                               # time in seconds to allow for running the bash script before timing out (assuming hanging code)

num_batches = 8                                                                   # = number of cores for maximum efficiency with large data sets
num_per_batch = math.ceil(len(file_tags)/num_batches)
if num_per_batch < 20:
    num_per_batch = 11
    num_batches = math.ceil(len(e2plus_to_test)*len(eps_points)/num_per_batch)             # if the data set is small then use fewer cores for a minimum batch size of 20 to make the overhead worth it

run_program = configure_script_writer(file_tags, num_batches, num_per_batch, allowed_time)

run_program("gampn")

timer_lapse_new = time.time()
print("finished running gampn in time = %.2f seconds" % (timer_lapse_new-timer_lapse))
timer_lapse = timer_lapse_new




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

for file in file_tags :
    
    # open the file and read
    with open(("GAM_"+file+".OUT"), 'r') as f:
        lines = f.readlines()
    
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
print("finished reading %d files\n" % len(asyrmo_orbitals))

del [efac_line, levels_header_line, fermi_level_line, first_index, last_index]
del [file, half_line, hash_index, l]
del [lines, orbitals_string, single_parity_index, whole_line, orbitals]
del [fermi_energies_hw, fermi_energies_mev, fermi_indices, fermi_parities]



#%%
''' WRITING ASYRMO.DAT FILE '''
# for each deformation, writes a .DAT file for ASYRMO using the file tag to name, and the provided example as a base


print("writing ASYRMO.DAT files...")


for t in range(len(file_tags)):

    inputs["current_orbitals"] = asyrmo_orbitals[t]
    
    if len(e2plus_to_test)>1:
        inputs["current_e2plus"] = e2plus_to_test[t]
    else:
        inputs["current_e2plus"] = e2plus_to_test[0]

    inputs["current_f016"] = "f016_"+file_tags[t]+".dat"
    inputs["current_f017"] = "f017_"+file_tags[t]+".dat"
    inputs["current_f018"] = "f018_"+file_tags[t]+".dat"
    
    asyrmo_dat_text = templates["asyrmo"] % inputs
   
    with open(("../Inputs/ASY_"+file_tags[t]+".DAT"), 'w') as f:
        f.write(asyrmo_dat_text)
    
print("finished writing %d ASY.DAT files" % (t+1))
del [t, asyrmo_dat_text]



#%%
''' WRITE AND RUN BASH SCRIPT TO EXECUTE ASYRMO '''
# writes a bash shell script to run each of the ASYRMO .DAT files written above
# after each one is run, the output file ASYRMO.out is copied to a new text file with a more descriptive name (based on file tag),
# so that the data is not lost when the next iteration runs asyrmo again and ASYRMO.out is overwritten

timer_lapse = time.time()

run_program("asyrmo")

timer_lapse_new = time.time()
print("finished running asyrmo in time = %.2f seconds" % (timer_lapse_new-timer_lapse))
timer_lapse = timer_lapse_new




#%%
''' WRITING PROBAMO.DAT FILE '''
# for each deformation, writes a .DAT file for PROBAMO, using the file tag to name, and the provided example as a base

for file in file_tags:
    
    inputs["current_f017"] = "f017_"+file+".dat"
    inputs["current_f018"] = "f018_"+file+".dat"

    probamo_dat_text = templates["probamo"] % inputs
    
    with open(("../Inputs/PROB_"+file+".DAT"), 'w') as f:
        f.write(probamo_dat_text)    


print("finished writing PROB.DAT files")

del [probamo_dat_text, file]


#%%
''' WRITE AND RUN BASH SCRIPT TO EXECUTE PROBAMO '''
# writes a bash shell script to run each of the PROBAMO .DAT files written above
# after each one is run, the output file PROBAMO.out is copied to a new text file with a more descriptive name (based on file tag),
# so that the data is not lost when the next iteration runs asyrmo again and PROBAMO.out is overwritten


timer_lapse = time.time()

run_program("probamo")

timer_lapse_new = time.time()
print("finished running probamo in time = %.2f seconds" 
      % (timer_lapse_new-timer_lapse))
timer_lapse = timer_lapse_new

#%%

''' READ PROBAMO.OUT FILE '''
# for each deformation (i.e. each PROBAMO.OUT file):
# the row containing the ground state is identified by energy=0.0
# the ground state spin and magnetic dipole moment are recorded
# the lowest energy (yrast) state with the same spin as the experimental first excited state is identified
# the "first excitation energy" is recorded
# similarly for the second and third experimental excited states

print("\nreading PROBAMO.OUT files...")

if inputs["spin_or_excitation"]=="excitation":
    gs_mag_mom = []
    gs_spins = []
    fx_energies = []
    sx_energies = []
    tx_energies = []
    
else: # inputs["spin_or_excitation"]=="spin"
    num_to_record = int(inputs["num_to_record"])-1 
    
    mag_mom_1 = []
    mag_mom_3 = []
    mag_mom_5 = []
    
    energies_1 = []
    energies_3 = []
    energies_5 = []
    
    for i in range(num_to_record+1):
        mag_mom_1.append([])
        mag_mom_3.append([])
        mag_mom_5.append([])
        
        energies_1.append([])
        energies_3.append([])
        energies_5.append([])
    
    


for f in range(len(file_tags)) :
    
    file = file_tags[f]

    probamo_out_file_name = "PROB_"+file+".OUT"
    
    if inputs["spin_or_excitation"]=="excitation":
        # locate the yrast (lowest energy) state of each spin
        fx_energy = math.inf                                                        # initialise to a large value so that the first '<' comparison will be true
        sx_energy = math.inf
        tx_energy = math.inf
    
    else: # inputs["spin_or_excitation"]=="spin"
        states_1 = 0 # for counting how many states of each spin have been calculated
        states_3 = 0
        states_5 = 0
    
    probamo_out_file = open(probamo_out_file_name, 'r')
    for line in probamo_out_file:
        line = line.strip()
        
        try:                                                                    # determine whether this line is a data row of the table (if a ' - ' is present (with a space either side!), then it is)
            dash_index = line.index(" - ")
        except ValueError:
            continue    
        else:
            spin_string = line[dash_index-4:dash_index].strip()
            if spin_string[1] == "/":
                spin = float(spin_string[0])/2                        # get the spin of this row
            else: 
                spin = float(spin_string[:2])/2
                
            final_spin_string = line[dash_index+11:dash_index+16].strip()              # get the spin of the final state (after transition)
            if final_spin_string[1] == "/":
                final_spin = float(final_spin_string[0])/2                        # get the spin of this row
            else: 
                final_spin = float(final_spin_string[:2])/2
            
            try:
                this_energy = float(line[:6].strip())
            except ValueError: # could not convert string to float: '0.0  1'
                this_energy = float(line[:5].strip())
            
            final_energy = float(line[dash_index+3:dash_index+10].strip())
                                                                # if this row doesn't contain data, ignore it
        
        if inputs["spin_or_excitation"]=="excitation":
            frac_spin = (str(int(float(spin)*2))+"/2")
            if line[0:3] == "0.0":                                                  # then this is the ground state
                gs_mag_mom.append(float(line[-8:].strip()))                                # get the ground state magnetic moment
                gs_spins.append(spin)                                   
            
            elif frac_spin == inputs["fx_spin"]:                                 # check to see if spin matches experimental first excited state
                if this_energy < fx_energy:                                         # this will eventually locate the yrast state of this spin (i.e. the first excited state)
                    fx_energy = this_energy
                    fx_line = line
            elif frac_spin == inputs["sx_spin"]:                                 # repeat for second excited state
                if this_energy < sx_energy:
                    sx_energy = this_energy
                    sx_line = line
            elif frac_spin == inputs["tx_spin"]:                                 # repeat for third excited state
                if this_energy < tx_energy:
                    tx_energy = this_energy
                    tx_line = line
                    
        else: # inputs["spin_or_excitation"]=="spin"
            if spin == 0.5  and final_spin == 0.5 and this_energy==final_energy:        # check that this line represents an internal transition (static moment)
                if states_1 > num_to_record: # only save the first three
                    continue
                mag_mom_1[states_1].append(float(line[-8:].strip()))                                 
                energies_1[states_1].append(float(line[0:6].strip())) 
                
                states_1 += 1
                                             
            
            elif spin == 1.5 and final_spin == 1.5 and this_energy==final_energy:
                if states_3 > num_to_record:
                    continue
                mag_mom_3[states_3].append(float(line[-8:].strip()))                                 
                energies_3[states_3].append(float(line[0:6].strip()))   
                
                states_3 += 1
                
            elif spin == 2.5 and final_spin == 2.5 and this_energy==final_energy:
                if states_5 > num_to_record:
                    continue
                mag_mom_5[states_5].append(float(line[-8:].strip()))                                 
                energies_5[states_5].append(float(line[0:6].strip()))   
                
                states_5 += 1
    
    if inputs["spin_or_excitation"]=="spin":
        for s in range(num_to_record+1):
            if len(mag_mom_1[s]) == f: # then this state couldn't be found in this file, so placehold with NaN
                mag_mom_1[s].append(np.NaN)
                    
            if len(mag_mom_3[s]) == f: # then this state couldn't be found in this file, so placehold with NaN
                mag_mom_3[s].append(np.NaN)
           
            if len(mag_mom_5[s]) == f: # then this state couldn't be found in this file, so placehold with NaN
                mag_mom_5[s].append(np.NaN)
            
            if len(energies_1[s]) == f: # then this state couldn't be found in this file, so placehold with NaN
                energies_1[s].append(np.NaN)
                    
            if len(energies_3[s]) == f: # then this state couldn't be found in this file, so placehold with NaN
                energies_3[s].append(np.NaN)
           
            if len(energies_5[s]) == f: # then this state couldn't be found in this file, so placehold with NaN
                energies_5[s].append(np.NaN)
    else:  #inputs["spin_or_excitation"]=="excitation":
        try:
            fx_energies.append(float(fx_line[0:5].strip()))
        except NameError:
            print("could not find an excited state with spin: " + inputs["fx_spin"] + " in file: " + probamo_out_file_name)
            # raise ValueError("could not find the excited states (no states of the appropriate spins)")
                
        try:
            sx_energies.append(float(sx_line[0:5].strip()))
        except NameError:
            print("could not find an excited state with spin: " + inputs["sx_spin"] + " in file: " + probamo_out_file_name)
            # raise ValueError("could not find the excited states (no states of the appropriate spins)")
        
        try:
            tx_energies.append(float(tx_line[0:5].strip()))
        except NameError:
            print("could not find an excited state with spin: " + inputs["tx_spin"] + " in file: " + probamo_out_file_name)
            # raise ValueError("could not find the excited states (no states of the appropriate spins)")
        
    probamo_out_file.close()

# convert collected data from string to float
#gs_mag_mom = [float(n) for n in gs_mag_mom]
#gs_spins = [float(n[0])/2 for n in gs_spins]

print("finished reading %d files\n" % len(asyrmo_orbitals))

#%%

''' PREPARE TO PLOT GRAPHS '''
# gamma is converted to radians.
# all arrays of data to be plotted are collected into a single matrix.
# lists of the titles to be used in each plot are constructed.
# similarly for contour levels, colour bar ticks, data axis labels and units, 
#       experimental data to compare results to, absolute error tolerance, 
#       and whether the data set is continuous or discrete.

for r in range(len(gamma_points)):
    gamma_points[r] *= np.pi/180

if inputs["spin_or_excitation"]=="excitation":
    # locate the yrast (lowest energy) state of each spin
    
    data_matrix = [gs_mag_mom, output_data["fermi_energies_mev"], output_data["fermi_indices"], 
                   gs_spins, fx_energies, sx_energies, tx_energies]

    graphs_to_print = ["Ground State Magnetic Dipole Moment", "Fermi Energy", "Fermi Level Parity And Index", 
                       "Ground State Spin", "First Excitation Energy", "Second Excitation Energy"]

    fermi_index_colour_levels = np.arange(min(output_data["fermi_indices"])-0.5, max(output_data["fermi_indices"])+1.5, 1.0) 
    gs_spin_colour_levels = np.arange(min(gs_spins)-0.5, max(gs_spins)+1.5, 1.0)
    contour_levels = [8, 10, fermi_index_colour_levels, 
                      gs_spin_colour_levels, 10, 10]

    fermi_index_cbar_ticks = np.arange(min(output_data["fermi_indices"]), max(output_data["fermi_indices"])+1.0, 1.0)
    fermi_index_cbar_ticks = [str(int(n)) for n in fermi_index_cbar_ticks]
    gs_spin_cbar_ticks = np.arange(min(gs_spins), max(gs_spins)+1.0, 1.0)
    gs_spin_cbar_ticks = [(str(int(n*2))+"/2") for n in gs_spin_cbar_ticks]
    cbar_ticks = [0,0, fermi_index_cbar_ticks, 
                  gs_spin_cbar_ticks,0,0]

    data_axis_labels = [r'μ / $μ_{N}$', 'fermi energy / MeV', 'fermi level parity and index',
                        'ground state spin I', "First Excitation Energy / keV", "Second Excitation Energy / keV"]

    experimental_data = [[], [], [], 
                         [], [], []]
    if "gs_mu" in inputs:
        experimental_data[0] = inputs["gs_mu"]
    if "gs_spin" in inputs:
        experimental_data[3] = inputs["gs_spin_float"]
    if "fx_energy" in inputs:
        experimental_data[4] = inputs["fx_energy"]
    if "sx_energy" in inputs:
        experimental_data[5] = inputs["sx_energy"]
        
    error_tolerance = [0.2, 0.0, 0.0, 0.1, 50, 50]                                #!!! these are a bit arbitrary...


else: # inputs["spin_or_excitation"]=="spin"

    data_matrix = [mag_mom_1, mag_mom_3, mag_mom_5, 
                   [output_data["fermi_energies_mev"]], [output_data["fermi_indices"]], 
                   energies_1, energies_3, energies_5]
    
    if "eps_max" in inputs or num_to_record == 0: 
        graphs_to_print = ["Magnetic Dipole Moment of Lowest Spin 1/2 State", "Magnetic Dipole Moment of Lowest Spin 3/2 State", "Magnetic Dipole Moment of Lowest Spin 5/2 State", 
                           "Fermi Energy", "Fermi Level Parity And Index", 
                           "Energy of Lowest Spin 1/2 State", "Energy of Lowest Spin 3/2 State", "Energy of Lowest Spin 5/2 State"]

    else:   
        graphs_to_print = ["Magnetic Dipole Moments of Spin 1/2 States", "Magnetic Dipole Moments of Spin 3/2 States", "Magnetic Dipole Moments of Spin 5/2 States", 
                           "Fermi Energy", "Fermi Level Parity And Index", 
                           "Energies of Spin 1/2 States", "Energies of Spin 3/2 States", "Energies of Spin 5/2 States"]

    fermi_index_colour_levels = np.arange(min(output_data["fermi_indices"])-0.5, max(output_data["fermi_indices"])+1.5, 1.0) 
    contour_levels = [8, 8, 8, 
                      10, fermi_index_colour_levels, 
                      10, 10, 10]

    fermi_index_cbar_ticks = np.arange(min(output_data["fermi_indices"]), max(output_data["fermi_indices"])+1.0, 1.0)
    fermi_index_cbar_ticks = [str(int(n)) for n in fermi_index_cbar_ticks]
    cbar_ticks = [0,0,0,
                  0,fermi_index_cbar_ticks, 
                  0,0,0]

    data_axis_labels = [r'μ / $μ_{N}$', r'μ / $μ_{N}$',r'μ / $μ_{N}$',
                        'fermi energy / MeV', 'fermi level parity and index',
                        "Energy / keV", "Energy / keV", "Energy / keV"]

    experimental_data = [[], [], [],
                         [], [],
                         [], [], []]
    if "gs_mu" in inputs:
        
        if inputs["gs_spin"] == "0.5":
            experimental_data[0] = float(inputs["gs_mu"])
        if inputs["gs_spin"] == "1.5":
            experimental_data[1] = float(inputs["gs_mu"])
        if inputs["gs_spin"] == "2.5":
            experimental_data[2] = float(inputs["gs_mu"])
        
    if "fx_energy" in inputs:
        if inputs["fx_spin"] == "1/2":
            experimental_data[5] = float(inputs["fx_energy"])
            if "fx_mu" in inputs:
                experimental_data[0] = float(inputs["fx_mu"])
        if inputs["fx_spin"] == "3/2":
            experimental_data[6] = float(inputs["fx_energy"])
            if "fx_mu" in inputs:
                experimental_data[1] = float(inputs["fx_mu"])
        if inputs["fx_spin"] == "5/2":
            experimental_data[7] = float(inputs["fx_energy"])
            if "fx_mu" in inputs:
                experimental_data[2] = float(inputs["fx_mu"])
    if "sx_energy" in inputs:
        if inputs["sx_spin"] == "1/2":
            experimental_data[5] = float(inputs["sx_energy"])
        if inputs["sx_spin"] == "3/2":
            experimental_data[6] = float(inputs["sx_energy"])
        if inputs["sx_spin"] == "5/2":
            experimental_data[7] = float(inputs["sx_energy"])
        
    error_tolerance = [0.2, 0.2, 0.2, 0.0, 0.0, 50, 50, 50]                                #!!! these are a bit arbitrary...


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

if inputs["deformation_input"] == "mesh":   
    
    
    agreed_points = [] # an array to store a list of points that have properties in agreement with the experimental data
    
    for g in range(len(graphs_to_print)):
        
        inputs["current_graph"] = graphs_to_print[g]
        print("\nplotting graph of %(current_graph)s variation..." % inputs)
        
        fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
        
        ax.set_thetamin(0)   
        ax.set_thetamax(60)  
        ax.set_rmax(1.0)
        
        theta_ticks = np.arange(0, 70, 10)  
        ax.set_xticks(np.radians(theta_ticks))
        ax.yaxis.set_major_formatter(FuncFormatter(lambda x, _: f'{x:.2f}'))    # set the number of decimal places 

        
        # create filled colour contour plots
        c_level_boundaries = contour_levels[g]
        if inputs["spin_or_excitation"] == "excitation":
            
            cax = ax.tricontourf(gamma_points, eps_points, data_matrix[g], levels=c_level_boundaries)
            cbar = plt.colorbar(cax, pad=0.1, label=data_axis_labels[g])
            if isinstance(c_level_boundaries, list):                                # then format for discrete values
                cax = ax.tricontourf(gamma_points, eps_points, data_matrix[g], levels=c_level_boundaries)
                cbar = plt.colorbar(cax, pad=0.1, label=data_axis_labels[g])
        
                ticks = [(c_level_boundaries[i] + c_level_boundaries[i+1]) / 2 
                         for i in range(len(c_level_boundaries) - 1)]               # Calculate midpoints of levels for tick placement
                cbar.set_ticks(ticks)
                cbar.set_ticklabels(cbar_ticks[g])
                
        else: # inputs["spin_or_excitation"] == "spin"
            this_data = data_matrix[g]
            
            
            cax = ax.tricontourf(gamma_points, eps_points, this_data[0], levels=c_level_boundaries)
            cbar = plt.colorbar(cax, pad=0.1, label=data_axis_labels[g])
            if isinstance(c_level_boundaries, list):                                # then format for discrete values
                cax = ax.tricontourf(gamma_points, eps_points, this_data[0], levels=c_level_boundaries)
                cbar = plt.colorbar(cax, pad=0.1, label=data_axis_labels[g])
        
                ticks = [(c_level_boundaries[i] + c_level_boundaries[i+1]) / 2 
                         for i in range(len(c_level_boundaries) - 1)]               # Calculate midpoints of levels for tick placement
                cbar.set_ticks(ticks)
                cbar.set_ticklabels(cbar_ticks[g])
                
            
        # mark the range in which the correct ground state spin was calculated
        if inputs["mark_spin"]==1:
            ax.tricontour(gamma_points, eps_points, data_matrix[3], 
                          levels=[float(inputs["gs_spin"])-0.5,float(inputs["gs_spin"])+0.5],  
                          colors=[(213/255,1,0)], linewidths=1.0)
        
        
        # create flags to record which ROIs should be included in the legend (each flag will only be turned on if that ROI is plotted)
        dot_hit_flag = False
        plus_hit_flag = False
        dot_miss_flag = False
        plus_miss_flag = False
        dot_flag = False
        plus_flag = False
        
        # plot the data point markers, and check which deformations give a result that matches experiment
        for r in range(len(gamma_points)):
            if experimental_data[g]:                                            # if the experimental data exists for comparison
                error = abs(data_matrix[g][r] - experimental_data[g])
                if error < error_tolerance[g]: # error/experimental_data[g] < 0.1: #!!! if they agree within 10%, plot the data point in red rather than white
                    if g==0 and output_data["fermi_parities"][r]==inputs["par"]:       # check that the ground state parity has been reproduced
                        agreed_points.append(r)
                    if output_data["fermi_parities"][r] == "-":
                        plt.polar(gamma_points[r], eps_points[r], 'r.')
                        dot_hit_flag = True
   
                    elif output_data["fermi_parities"][r] == "+":
                        plt.polar(gamma_points[r], eps_points[r], 'r+')
                        plus_hit_flag = True
                else:
                    if r in agreed_points:
                        del agreed_points[agreed_points.index(r)]               # this point is no longer in agreement, so remove it from the list
                    if output_data["fermi_parities"][r] == "-":
                        plt.polar(gamma_points[r], eps_points[r], 'w.')
                        dot_miss_flag = True
                    else: # fermi_parities[r] == "+":
                        plt.polar(gamma_points[r], eps_points[r], 'w+')
                        plus_miss_flag = True
                exp, = cbar.ax.plot([0, 1], [experimental_data[g], experimental_data[g]], 'r-')
                        
            # if there is no experimental data available for comparison, just plot all points in white            
            elif output_data["fermi_parities"][r] == "-":
                plt.polar(gamma_points[r], eps_points[r], 'w.')
                dot_flag = True
            else: # fermi_parities[r] == "+":
                plt.polar(gamma_points[r], eps_points[r],'w+')
                plus_flag = True
                
        # plot four data points (unseen outside the data range) with desired formatting to use for the legend
        dot, = plt.polar(3.0, 0.2, 'k.', label="negative parity")
        plus, = plt.polar(3.0, 0.2, 'k+', label="positive parity")
        dot_miss, = plt.polar(3.0, 0.2, 'k.', label="parity (-) miss")
        plus_miss, = plt.polar(3.0, 0.2, 'k+', label="parity (+) miss")
        dot_hit, = plt.polar(3.0, 0.2, 'r.', label="parity (-) match")
        plus_hit, = plt.polar(3.0, 0.2, 'r+', label="parity (+) match")
        exp, = plt.polar(3.0, 0.2, 'r-', label="experiment")
        spin, = plt.polar(3.0, 0.2, '-', color=(213/255,1,0), linewidth=2.0, 
                          label="boundary of correct\nground state spin")
        
        
        legend_handles=[dot, plus]
        '''
        if dot_hit_flag: legend_handles.append(dot_hit)
        if plus_hit_flag: legend_handles.append(plus_hit)
        if dot_miss_flag: legend_handles.append(dot_miss)
        if plus_miss_flag: legend_handles.append(plus_miss)
        if dot_flag: legend_handles.append(dot)
        if dot_flag: legend_handles.append(plus)
        '''
        
        annotation_handles = []
        if experimental_data[g]:
            annotation_handles.append(exp)
        if inputs["mark_spin"]==1:
            annotation_handles.append(spin)
        
        legend = ax.legend(handles=legend_handles, loc="upper left", bbox_to_anchor=(-0.3, 0.78))
        legend.set_title("data points\n---------------")
        
        if len(annotation_handles) > 0:
            annotation_legend = ax.legend(handles=annotation_handles,  loc="upper left", bbox_to_anchor=(-0.3, 1.12))
            annotation_legend.set_title('annotations\n-------------')
            plt.gca().add_artist(legend) # adding annotation_legend removed the original legend, so we have to manually add it back in
        
        
        # add title and axis labels
        ax.set_title('%(current_graph)s of %(nucleus)s' % inputs, va='bottom', y=1.1)                      
        plt.xlabel("ε")
        ax.text(65*np.pi/180, ax.get_rmax()*1.05, "γ", ha='center', va='center') # gamma axis label
        
        
        
        plt.show()
        
    print("\n\nnumber of points that agreed with experimental data: " +
          str(len(agreed_points)))
    print("\n")
    for p in agreed_points:
        print("point %d: ε=%.3f, γ=%d" 
              % (p, eps_points[p], gamma_points[p]*180/np.pi))
    print("\n")
        












elif inputs["deformation_input"] == "eps" or inputs["deformation_input"] == "gamma":                                              # then plot a line graph of data variation with eps (or gamma)

    # locate the range in which the excitation energy data is meaningful (the range in which the ground state spin is correct)
    correct_spin_range = []                                                     #!!! do this with a mask instead? would allow easy combination of conditions (parity) by multiplication...
    start_flag = False
    
    if inputs["mark_spin"]==1:
        for i in range(len(gs_spins)):
            if gs_spins[i] == experimental_data[3] and not start_flag:              # correct, and range hasn't started yet
                start_flag = True
                correct_spin_range.append(i)
            elif gs_spins[i] != experimental_data[3] and start_flag:                # incorrect, and range has started
                start_flag = False
                correct_spin_range.append(i)
        
    # plot graphs
    for g in range(len(graphs_to_print)):
        inputs["current_graph"] = graphs_to_print[g]
        print("plotting graph of %(current_graph)s variation..." % inputs)
        
        fig, ax = plt.subplots()                                                # create figure and axes
        
        # plot four data points (unseen outside the data range) with desired formatting to use for the legend
        dot, = plt.plot(3.0, 0.0, 'k.', label="negative parity")
        plus, = plt.plot(3.0, 0.0, 'k+', label="positive parity")
        
        # create flags to record which ROIs should be included in the legend (each flag will only be turned on if that ROI is plotted)
        dot_flag = False
        plus_flag = False
        
          
        
        
        
        # if eps was varied and gamma was fixed:
        if inputs["deformation_input"] == "eps" :
            inputs["line"] = "ε"
            inputs["fixed"] = "γ / º"
            # print("plotting line graph of variation with eps...")
            
            pad = 0.05*(eps_to_test[-1]-eps_to_test[0])
            ax.set_xlim([eps_to_test[0]-pad, eps_to_test[-1]+pad]) 
            
            # if experimental data is available, plot it in red for easy comparison
            if experimental_data[g]:
                exp, = plt.plot(eps_to_test, 
                       np.full(len(eps_to_test), float(experimental_data[g])), 
                       'r-', label="experimental value")
            
            # mark the range in which the correct ground state spin was calculated
            if inputs["mark_spin"]==1:
                for r in range(len(correct_spin_range)):
                    if correct_spin_range[r] == 0:                                  # the first value of eps in the range has the correct spin
                        start_range = (eps_to_test[correct_spin_range[r]]-
                                       float(inputs["eps"][2])/2)          
                    else:
                        start_range = np.mean([eps_to_test[correct_spin_range[r]-1], 
                                               eps_to_test[correct_spin_range[r]]])
                    
                    correct_spin, = plt.plot([start_range, start_range], 
                                             [min(data_matrix[g])+0.05*max(data_matrix[g]), 
                                              max(data_matrix[g])*1.05], 
                                             'g-', label="range of correct spin")   # this plots the front and end edges of the box
                    if r%2==0:
                        if r+1 == len(correct_spin_range):                          # the last value of eps in the range has the correct spin
                            end_range = (eps_to_test[-1]+
                                         float(inputs["eps"][2])/2)
                        else:
                            end_range = np.mean([eps_to_test[correct_spin_range[r+1]-1], 
                                                 eps_to_test[correct_spin_range[r+1]]])
                        
                        plt.plot([start_range, end_range],                          # this plots the bottom edge of the box
                                 [min(data_matrix[g])+0.05*max(data_matrix[g]), 
                                  min(data_matrix[g])+0.05*max(data_matrix[g])], 'g-')
                        plt.plot([start_range, end_range], 
                                 [max(data_matrix[g])*1.05, 
                                  max(data_matrix[g])*1.05], 'g-')                  # this plots the top edge of the box
                
                
            # now plot the actual data
            if inputs["spin_or_excitation"]=="excitation":
                data, = plt.plot(eps_to_test, data_matrix[g], 'k-', label="γ = %s" % gamma_to_test[0])
                for p in range(len(eps_to_test)):
                    if output_data["fermi_parities"][p] == "-":
                        plt.plot(eps_to_test[p], data_matrix[g][p], 'k.', label='negative parity')
                        dot_flag = True
                    elif output_data["fermi_parities"][p] == "+":
                        plt.plot(eps_to_test[p], data_matrix[g][p], 'k+', label='positive parity')
                        plus_flag = True
                        
            else: # inputs["spin_or_excitation"]=="spin":
                this_data = data_matrix[g]
                line_colours = ['k-', 'b-', 'y-']
                line_labels = ["lowest energy", "second lowest energy", "third lowest energy"]
                data_handles = []
                for s in range(len(this_data)):
                    
                    data, = plt.plot(eps_to_test, this_data[s], line_colours[s], label=line_labels[s])
                    data_handles.append(data)
                    for p in range(len(this_data[s])):
                        if output_data["fermi_parities"][p] == "-":
                            plt.plot(eps_to_test[p], this_data[s][p], 'k.', label='negative parity')
                            dot_flag = True
                        elif output_data["fermi_parities"][p] == "+":
                            plt.plot(eps_to_test[p], this_data[s][p], 'k+', label='positive parity')
                            plus_flag = True
                legend_title = "γ = %s" % gamma_to_test[0]
            
            
            
        # if gamma was varied and eps was fixed:   
        else:                                                                   # inputs["line"] == "gamma":
            inputs["line"] = "γ / º"
            inputs["fixed"] = "ε"
            # print("plotting line graph of variation with gamma...")
            
            pad = 0.05*(gamma_to_test[-1]-gamma_to_test[0])
            ax.set_xlim([gamma_to_test[0]-pad, gamma_to_test[-1]+pad])
            
            # if experimental data is available, plot it in red for easy comparison
            if experimental_data[g]:
                exp, = plt.plot(gamma_to_test, 
                       np.full(len(gamma_to_test), float(experimental_data[g])), 
                       'r-', linewidth=3, label="experimental value")
            
            # mark the range in which the correct ground state spin was calculated
            if inputs["mark_spin"]==1:
                for r in range(len(correct_spin_range)):
                    if correct_spin_range[r] == 0:                                  # the first value of gamma in the range has the correct spin
                        start_range = (gamma_to_test[correct_spin_range[r]]-
                                       float(inputs["gamma"][2])/2)
                    else:
                        start_range = np.mean([gamma_to_test[correct_spin_range[r]-1], 
                                               gamma_to_test[correct_spin_range[r]]])
                        
                    correct_spin, = plt.plot([start_range, start_range], 
                                             [min(data_matrix[g])-0.05*max(data_matrix[g]), 
                                              max(data_matrix[g])*1.05], 
                                             'g-', label="range of correct spin")
                    
                    if r%2==0:
                        if r+1 == len(correct_spin_range):                          # the last value of gamma in the range has the correct spin
                            end_range = (gamma_to_test[-1]+
                                         float(inputs["gamma"][2])/2)
                            
                            correct_spin, = plt.plot([end_range, end_range], 
                                            [min(data_matrix[g])-0.05*max(data_matrix[g]), 
                                             max(data_matrix[g])*1.05], 
                                            'g-', label="range of correct spin")    # this plots the front and end edges of the box
                           
                        else:
                            end_range = np.mean([gamma_to_test[correct_spin_range[r+1]-1], 
                                                 gamma_to_test[correct_spin_range[r+1]]])
                        
                        plt.plot([start_range, end_range],                          # this plots the bottom edge of the box
                                 [min(data_matrix[g])-0.05*max(data_matrix[g]),
                                  min(data_matrix[g])-0.05*max(data_matrix[g])], 'g-')
                        plt.plot([start_range, end_range], 
                                 [max(data_matrix[g])*1.05, 
                                  max(data_matrix[g])*1.05], 'g-')                  # this plots the top edge of the box
                
                
            # plot the actual data (first a line graph, then the points seperately so that the parities can be indicated by the markers)
            
            if inputs["spin_or_excitation"]=="excitation":
                data, = plt.plot(gamma_to_test, data_matrix[g], 'k-', label="ε = %s" % eps_to_test[0])
                for p in range(len(gamma_to_test)):
                    if output_data["fermi_parities"][p] == "-":
                        dot, = plt.plot(gamma_to_test[p], data_matrix[g][p], 'k.', label='negative parity')
    
                    elif output_data["fermi_parities"][p] == "+":
                        plus, = plt.plot(gamma_to_test[p], data_matrix[g][p], 'k+', label='positive parity')
            
            else: # inputs["spin_or_excitation"]=="spin":
                this_data = data_matrix[g]
                line_colours = ['k-', 'b-', 'y-', 'c-', 'm-']
                line_labels = ["lowest energy", "second lowest energy", "third lowest energy", "fourth lowest energy", "fifth lowest energy"]
                data_handles = []
                for s in range(len(this_data)):
                    
                    data, = plt.plot(gamma_to_test, this_data[s], line_colours[s], label=line_labels[s])
                    data_handles.append(data)
                    for p in range(len(this_data[s])):
                        if output_data["fermi_parities"][p] == "-":
                            dot, = plt.plot(gamma_to_test[p], this_data[s][p], 'k.', label='negative parity')
        
                        elif output_data["fermi_parities"][p] == "+":
                            plus, = plt.plot(gamma_to_test[p], this_data[s][p], 'k+', label='positive parity')
                legend_title = "ε = %s" % eps_to_test[0]
                
            
        
        # add title and axis labels
        ax.set_title('%(current_graph)s in %(nucleus)s' 
                     % inputs, va='bottom', y=1.1)                     
        plt.xlabel("%(line)s" % inputs)
        plt.ylabel(data_axis_labels[g])
        
        # add legend (depending on what ROIs have been drawn)
        if correct_spin_range and experimental_data[g] and inputs["mark_spin"]==1: 
            legend_handles=[exp, correct_spin, data, dot]
        elif experimental_data[g]:
            legend_handles=[exp, data]
        else: 
            legend_handles=[data]
            
        if dot_flag:
            legend_handles.append(dot)
        
        if plus_flag:
            legend_handles.append(plus)
        
        if inputs["spin_or_excitation"] == "spin":
             legend_handles.remove(data)
             legend_handles += data_handles
    
        legend = ax.legend(handles = list(reversed(legend_handles)))
        
        if inputs["spin_or_excitation"] == "spin":
            legend.set_title(legend_title)
            
        plt.show()










elif len(e2plus_to_test)>1:                                                  # then plot a line graph of data variation with e2plus

    # locate the range in which the excitation energy data is meaningful (the range in which the ground state spin is correct)
    correct_spin_range = []                                                     #!!! do this with a mask instead? would allow easy combination of conditions (parity) by multiplication...
    start_flag = False
    
    if inputs["mark_spin"]==1:
        for i in range(len(gs_spins)):
            if gs_spins[i] == experimental_data[3] and not start_flag:              # correct, and range hasn't started yet
                start_flag = True
                correct_spin_range.append(i)
            elif gs_spins[i] != experimental_data[3] and start_flag:                # incorrect, and range has started
                start_flag = False
                correct_spin_range.append(i)
        
    # plot graphs
    for g in range(len(graphs_to_print)):
        inputs["current_graph"] = graphs_to_print[g]
        print("plotting graph of %(current_graph)s variation..." % inputs)
        
        fig, ax = plt.subplots()                                                # create figure and axes
        
        # plot four data points (unseen outside the data range) with desired formatting to use for the legend
        dot, = plt.plot(3.0, 0.0, 'k.', label="negative parity")
        plus, = plt.plot(3.0, 0.0, 'k+', label="positive parity")
        
        # create flags to record which ROIs should be included in the legend (each flag will only be turned on if that ROI is plotted)
        dot_flag = False
        plus_flag = False
        
        
        pad = 0.05*(e2plus_to_test[-1]-e2plus_to_test[0])
        ax.set_xlim([e2plus_to_test[0]-pad, e2plus_to_test[-1]+pad]) 
        
        # if experimental data is available, plot it in red for easy comparison
        if experimental_data[g]:
            exp, = plt.plot(e2plus_to_test, 
                   np.full(len(e2plus_to_test), float(experimental_data[g])), 
                   'r-', label="experimental value")
        
        # mark the range in which the correct ground state spin was calculated
        if inputs["mark_spin"]==1:
            for r in range(len(correct_spin_range)):
                if correct_spin_range[r] == 0:                                  # the first value of eps in the range has the correct spin
                    start_range = (e2plus_to_test[correct_spin_range[r]]-
                                   float(inputs["eps"][2])/2)          
                else:
                    start_range = np.mean([e2plus_to_test[correct_spin_range[r]-1], 
                                           e2plus_to_test[correct_spin_range[r]]])
                
                correct_spin, = plt.plot([start_range, start_range], 
                                         [min(data_matrix[g])+0.05*max(data_matrix[g]), 
                                          max(data_matrix[g])*1.05], 
                                         'g-', label="range of correct spin")   # this plots the front and end edges of the box
                if r%2==0:
                    if r+1 == len(correct_spin_range):                          # the last value of eps in the range has the correct spin
                        end_range = (e2plus_to_test[-1]+
                                     float(inputs["e2plus"][2])/2)
                    else:
                        end_range = np.mean([e2plus_to_test[correct_spin_range[r+1]-1], 
                                             e2plus_to_test[correct_spin_range[r+1]]])
                    
                    plt.plot([start_range, end_range],                          # this plots the bottom edge of the box
                             [min(data_matrix[g])+0.05*max(data_matrix[g]), 
                              min(data_matrix[g])+0.05*max(data_matrix[g])], 'g-')
                    plt.plot([start_range, end_range], 
                             [max(data_matrix[g])*1.05, 
                              max(data_matrix[g])*1.05], 'g-')                  # this plots the top edge of the box
            
                
        # now plot the actual data
        if inputs["spin_or_excitation"]=="excitation":
            data, = plt.plot(e2plus_to_test, data_matrix[g], 'k-', label="ε = %s, γ = %s" % (eps_to_test[0], gamma_to_test[0]))
            for p in range(len(e2plus_to_test)):
                if output_data["fermi_parities"][p] == "-":
                    plt.plot(e2plus_to_test[p], data_matrix[g][p], 'k.', label='negative parity')
                    dot_flag = True
                elif output_data["fermi_parities"][p] == "+":
                    plt.plot(e2plus_to_test[p], data_matrix[g][p], 'k+', label='positive parity')
                    plus_flag = True
                    
        else: # inputs["spin_or_excitation"]=="spin":
            this_data = data_matrix[g]
            line_colours = ['k-', 'b-', 'y-']
            line_labels = ["lowest energy", "second lowest energy", "third lowest energy"]
            data_handles = []
            for s in range(len(this_data)):
                
                data, = plt.plot(e2plus_to_test, this_data[s], line_colours[s], label=line_labels[s])
                data_handles.append(data)
                for p in range(len(this_data[s])):
                    if output_data["fermi_parities"][p] == "-":
                        plt.plot(e2plus_to_test[p], this_data[s][p], 'k.', label='negative parity')
                        dot_flag = True
                    elif output_data["fermi_parities"][p] == "+":
                        plt.plot(e2plus_to_test[p], this_data[s][p], 'k+', label='positive parity')
                        plus_flag = True
            legend_title = "ε = %s, γ = %s" % (eps_to_test[0], gamma_to_test[0])
        
            
        
        # add title and axis labels
        ax.set_title('%(current_graph)s in %(nucleus)s' 
                     % inputs, va='bottom', y=1.1)                     
        plt.xlabel("E2PLUS / MeV")
        plt.ylabel(data_axis_labels[g])
        
        # add legend (depending on what ROIs have been drawn)
        if correct_spin_range and experimental_data[g] and inputs["mark_spin"]==1: 
            legend_handles=[exp, correct_spin, data, dot]
        elif experimental_data[g]:
            legend_handles=[exp, data]
        else: 
            legend_handles=[data]
            
        if dot_flag:
            legend_handles.append(dot)
        
        if plus_flag:
            legend_handles.append(plus)
        
        if inputs["spin_or_excitation"] == "spin":
             legend_handles.remove(data)
             legend_handles += data_handles
    
        legend = ax.legend(handles = list(reversed(legend_handles)))
        
        if inputs["spin_or_excitation"] == "spin":
            legend.set_title(legend_title)
            
        plt.show()

# note how long it took
timer_lapse_new = time.time()
print("finished plotting graphs in time = %.2f seconds" % (timer_lapse_new-timer_lapse))
print("total runtime = %.2f seconds" % (timer_lapse_new-timer_start))





