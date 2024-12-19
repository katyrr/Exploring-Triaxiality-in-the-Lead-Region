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
from matplotlib import rc                                                       # for TeX formatting
import time                                                                     # for checking how long it took to run

plt.rcParams['figure.dpi'] = 150
rc('text', usetex=False)

config_file_name = "../Inputs/config.txt"

                                                                 
start_wall_time = time.time()


#%%



''' OPEN AND READ CONFIG FILE '''
# ignores empty lines, and lines beginning with * (to mark a comment)
# counts and outputs the number of non-empty and non-comment lines (should be equal to the number of variables input in the config file)
# expects each line to have the format:    var_name value
# outputs a warning if the format is unexpected, otherwise splits the string at " " and assigns the variable to a dictionary
# tests whether deformation has been input as a range, mesh, or single value, and creates a list of values to test (which may contain just one value)


print("reading from config file...\n")

line_count = 0                                                                  # start a counter for lines read (not including ignored lines)
all_lines_count = 0                                                             # start a counter for all lines  (including empty lines and comments)
bad_lines = []                                                                  # create an empty list to store indices of bad lines, so that they can be reported

input_settings = {}                                                             # create an empty dictionary to hold settings for input

config_file = open(config_file_name, 'r')
for line in config_file:
    
    all_lines_count += 1
    
    line_string = line.strip()                                                  # extract text from line
    if line_string == "" : continue                                             # skip blank rows
    if line_string[0] == "*" : continue                                         # skip comment lines
    line_count += 1 
    
    split_string = line_string.split(" ")                                       # split into name and value
        
    for n in range(len(split_string)):
        word = split_string[n]
        if word[0] == '*':                                                      # remove inline comments 
            del split_string[n:]
            break  
    
    if split_string[0]=="eps" or split_string[0]=="gamma" or split_string[0]=="single":
          expected_num_words = 3                                                  # these variables have two values associated with the variable name
    else: expected_num_words = 2                                                  # all other variables have only one value associated with the variable name
    
    # warn if the input looks unexpected
    if len(split_string)>expected_num_words :                                 
        print("Line "+str(all_lines_count)+ ''' has too many words! 
              Check var name contains no spaces, and separate list 
              values with commas not spaces. Ignoring extra words.\n''')
        bad_lines.append(line_count)
        
    elif len(split_string)<expected_num_words :
        print("Line "+str(all_lines_count)+''' is missing either 
              variable name or value. Check they are both present 
              and seperated with a space. Skipping this line.\n''')
        continue
        bad_lines.append(line_count)
    
    
    # save deformation input in the input_settings dictionary
    if split_string[0]=="eps" or split_string[0]=="gamma" or split_string[0]=="single": 
        
        # warn if multiple sets of deformation parameters have been input
        if "eps" in input_settings or "gamma" in input_settings or "eps_max" in input_settings:
            print('''Deformation parameter sets have been input multiple times. 
                  Check that config file only contains one set of deformation inputs. 
                  Ignoring excess inputs.''')
            continue
        
        input_settings["eps"] = split_string[1]
        input_settings["eps"] = input_settings["eps"].split(",") 
        if len(input_settings["eps"]) > 1 :                                     # then a range of eps values has been input
            input_settings["line"] = "eps"
            eps_to_test =  np.arange(float(input_settings["eps"][0]),           # add a step to the end value (because arange function is non-inclusive)
                                    float(input_settings["eps"][1])+float(input_settings["eps"][2]), 
                                    float(input_settings["eps"][2]))
        else : eps_to_test = [float(input_settings["eps"][0])]                  # even though it only has one element, a list is easier to input to a for loop later
        
        input_settings["gamma"] = split_string[2]
        input_settings["gamma"] = input_settings["gamma"].split(",")            # same for gamma
        if len(input_settings["gamma"]) > 1 :                                             
            input_settings["line"] = "gamma"
            gamma_to_test = np.arange(float(input_settings["gamma"][0]),  
                                    float(input_settings["gamma"][1])+float(input_settings["gamma"][2]), 
                                    float(input_settings["gamma"][2]))
        else : gamma_to_test = [float(input_settings["gamma"][0])] 
        
        
        eps_points = []
        gamma_points = []
        for l in range(len(eps_to_test)):
            for p in range(len(gamma_to_test)):
                eps_points.append(eps_to_test[l])
                gamma_points.append(gamma_to_test[p])
                
        
    elif split_string[0]=="mesh":
        
        if "eps" in input_settings or "gamma" in input_settings or "eps_max" in input_settings:
            print('''Deformation parameter sets have been input multiple times. 
                  Check that config file only contains one set of deformation inputs. 
                  Ignoring excess inputs.''')
            continue
        
        input_settings["eps_max"] = split_string[1]
        input_settings["eps_max"] = float(input_settings["eps_max"])            # convert from string to float
        
        # arrange a mesh of evenly distributed (eps,gamma) points to test over, such that the largest value of eps is tested with outer_points gamma values 
        outer_points = int(input_settings["outer_points"])                      # the total number of data points generated = 0.5*outer_points*(outer_points + 1)
        eps_to_test = np.linspace(0.001,                                        # start from eps=0.001 because eps=0 is not an allowed input when it comes to asyrmo
                                  input_settings["eps_max"], num=outer_points)
        eps_points = []
        gamma_points = []
        for l in range(len(eps_to_test)):
            for p in range(l+1):
                eps = np.round(eps_to_test[l], 3)
                eps_points.append(eps)
                if l==0: gamma = 0
                else: gamma = np.round(p*(60/l)) 
                gamma_points.append(gamma)
        eps_points = np.array(eps_points)
        gamma_points = np.array(gamma_points)
        
    else:
        # save other (non-deformation) parameters
        input_settings[split_string[0]] = split_string[1]
config_file.close()
    
input_settings["A"] = int(input_settings["A"])                                  # convert A and Z from strings to ints
input_settings["Z"] = int(input_settings["Z"])

input_settings["num_orbs"] = int(input_settings["num_orbs"])       

if "fx_spin" in input_settings:
    input_settings["fx_spin"] = (str(int(float(input_settings["fx_spin"])*2))+"/2")
else: input_settings["fx_spin"] = "0"
if "sx_spin" in input_settings:
    input_settings["sx_spin"] = (str(int(float(input_settings["sx_spin"])*2))+"/2")
else: input_settings["sx_spin"] = "0"
if "tx_spin" in input_settings:
    input_settings["tx_spin"] = (str(int(float(input_settings["tx_spin"])*2))+"/2")
else: input_settings["tx_spin"] = "0"


print(input_settings)
print("\n\nFinished reading lines: "+str(line_count))
print(str(len(bad_lines))+" of these had unexpected formatting.\n\n")



#%%
''' CALCULATE FERMI LEVEL '''
# using input A and Z, work out which is odd
# check that the odd particle matches the input of nneupr (the input nneupr will take precedence, but a mismatch may indicate an input error)
# half and ceiling for the overall index of the fermi level orbital
# generate the orbitals input for gampn (e.g. orbitals_string = "+4 19 20 21 22")

print("calcualting fermi level...")


input_settings["N"] = input_settings["A"]-input_settings["Z"]                   # calculate N = A - Z

if input_settings["nneupr"] == "1":
    print("calculating for odd protons...")
    input_settings["fermi_level"] = math.ceil(input_settings["Z"]/2)            # calculate fermi level orbital index
    if input_settings["Z"]%2 == 0:
        print('''this nucleus has even Z but nneupr = 1; check inputs 
              (the code will treat this as the core and assume an odd proton on top)''')
    
elif input_settings["nneupr"]  == "-1":    
    print("calculating for odd neutrons...")
    input_settings["fermi_level"] = math.ceil(input_settings["Z"]/2)
    if input_settings["Z"]%2 == 0:
        print('''this nucleus has even N but nneupr = -1; check inputs 
              (the code will treat this as the core and assume an odd neutron on top)''')

else:
    print("missing nneupr input")
    raise ValueError


print("fermi level = " + str(input_settings["fermi_level"]) + "\n")

# generate a string that contains the number of orbitals and their indices, in the correct format for input to gampn
gampn_orbitals = input_settings["par"]
gampn_orbitals += str(input_settings["num_orbs"])                               # e.g. orbitals_string = "4 19 20 21 22"

first_index = input_settings["fermi_level"]//2 - input_settings["num_orbs"]//2  # select [num_orbs] orbitals (the fermi level plus [num_orbs//2] either side)
last_index = input_settings["fermi_level"]//2 + input_settings["num_orbs"]//2
if input_settings["num_orbs"]%2 == 1:                                           # if an even number is requested, we need to add one to the final index to ensure the correct number are included
    last_index += 1
orbitals = np.r_[first_index:last_index]                                        # generate a list of orbitals in unit steps inside this range with list slicing
                                                 
for l in orbitals:
    gampn_orbitals += " "
    gampn_orbitals += str(l)
    
input_settings["gampn_orbitals"] = gampn_orbitals
    
    
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

write_count = 0                                                                 # start counting how many files get written
written_file_tags = []                                                          # create a list to record which files were written

for p in range(len(gamma_points)):
    
    input_settings["current_eps"] = eps_points[p]                               # store eps in the dict for ease of use when string formatting below
    input_settings["current_gamma"] = gamma_points[p]
    
    file_tag = "%s_e%.3f_g%.1f" % (input_settings["nucleus"], 
                                   input_settings["current_eps"], 
                                   input_settings["current_gamma"])
    
    input_settings["current_f002"] = "f002_"+file_tag+".dat"
    input_settings["current_f016"] = "f016_"+file_tag+".dat"
    input_settings["current_f017"] = "f017_"+file_tag+".dat"
    
    try:                                                                        # only overwrite the existing input file if this code is successful
        new_input_text =     '''
        
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
    
    ''' % input_settings
            
    except KeyError:
        print("Could not write GAM_"+file_tag+".DAT because an input is missing."+
              " Check config file is correct (and saved!) Will not attempt to overwrite existing file.")
        raise
        
    else: 
        gampn_dat_file_path = "../Inputs/GAM_"+file_tag+".DAT" 
        gampn_dat_file = open(gampn_dat_file_path, 'w')
        gampn_dat_file.write(new_input_text)
        gampn_dat_file.close()     

    write_count += 1
    written_file_tags.append(file_tag)

print("%d input files were written, \nfor eps in range [%.3f, %.3f], \nand gamma in range [%d, %d].\n" 
      % (write_count, eps_points[0], eps_points[-1], gamma_points[0], gamma_points[-1]))




#%%
''' WRITE AND RUN BASH SCRIPT TO EXECUTE GAMPN '''
# writes a bash shell script to run each of the GAMPN .DAT files written above
# after each one is run, the output file GAMPN.out is copied to a new text file with a more descriptive name (based on file tag),
# so that the data is not lost when the next iteration runs gampn again and GAMPN.out is overwritten


shell_script_file_path = "../RunGAMPN.sh"                                       # this is where the shell script will be created
timer_lapse = time.time()
print("running gampn; time so far elapsed = %.2f seconds" % (timer_lapse-start_wall_time))
allowed_time = 0.1*write_count + 10                                            # each file takes ~ 0.06 seconds to run, as a rough average, so allow 0.1 seconds per file to be safe, with an overhead of 0.2                                                               # time in seconds to allow for running the bash script before timing out (assuming hanging code)

new_shell_script_text = ""
for file in written_file_tags :
    new_shell_script_text += ("\n./../../../Executables/MO/gampn < ../Inputs/GAM_"+file+".DAT")
    new_shell_script_text += ("\ncp GAMPN.out GAM_"+file+".OUT")                # copy GAMPN.out to a new txt file with a more descriptive name

new_shell_script_text += "\n\necho message from terminal: finished running gampn"

shell_script_file = open(shell_script_file_path, 'w')
shell_script_file.write(new_shell_script_text)
shell_script_file.close()                                                       

gampn = subprocess.Popen(["sh", "./../RunGAMPN.sh"])                            # start gampn as a subprocess (#!!! asynchronous - potential for parallelisation??)
gampn.wait(allowed_time)                                                        # wait to ensure it has finished before starting to read outputs, if it takes longer than 10 seconds, timeout (#!!! may need to allow a longer timeout for larger datasets)


new_timer_lapse = time.time()
print("finished running gampn in time = %.2f seconds" % (new_timer_lapse-timer_lapse))
timer_lapse = new_timer_lapse



                                                       





#%%
''' READ GAMPN.OUT FILE '''
# for each deformation (i.e. each GAMPN.OUT file):
#   determine the line number of the fermi level orbital in the GAMPN.OUT file
#   read that line from GAMPN.OUT and get the energy, parity, and single-party-index
#   generate the orbitals input for asyrmo (e.g. orbitals_string = "+11 19 20 21 22 23 24 25 26 27 28 29")

print("\nreading GAMPN.OUT files...")

asyrmo_orbitals    = []                                                         # an empty array to store inputs for asyrmo.dat
fermi_energies     = []                                                         # for the fermi energy at each deformation
fermi_indices      = []                                                         # for the index and parity of the fermi level
fermi_parities     = [] 

for file in written_file_tags :
    
    # open the file and read
    gampn_out_file_name = "GAM_"+file+".OUT"
    gampn_out_file = open(gampn_out_file_name, 'r')
    lines = gampn_out_file.readlines()
    gampn_out_file.close()
    
    # get the EFAC value (conversion factor from hw to MeV)
    efac_header = lines.index("     KAPPA    MY     EPS   GAMMA    EPS4     EPS6     W0/W00   NMAX  COUPL     OMROT      EFAC      QFAC\n")
    efac = float(lines[efac_header+1][85:95].strip())
    
    # locate the header line of the table of single particle levels, and calculate the location of the fermi level relative to the header line
    header = lines.index("   #   ENERGY +/-(#)    <Q20>    <Q22>     <R2>     <JZ>      #   ENERGY +/-(#)    <Q20>    <Q22>     <R2>     <JZ>\n")
    fermi_level_line = input_settings["fermi_level"]+header+1                   # calculate the line number of the fermi level in the GAMPN.OUT file (indexed from zero!)
    if input_settings["fermi_level"] > 40:
        fermi_level_line -= 40
        whole_line = lines[fermi_level_line]
        half_line = whole_line[60:-1].strip()                                   # get only the second half of the line
    else:
        whole_line = lines[fermi_level_line]
        half_line = whole_line[0:60].strip()                                    # get only the first half of the line
        
    # from the half-line containing data on the fermi level orbital, extract the parity, single-parity-index, and energy
    hash_index = half_line.index("#")                                           # use the index of '#' as a reference point 
    fermi_parities.append(half_line[hash_index-2])                              #!!! rather than reading the parity out like this, just set it as a parameter? 
    fermi_energies.append(float(half_line[hash_index-10 : hash_index-4])*efac)              #!!! convert to float here? and convert to MeV (read EFAC conversion factor from output)
    
    single_parity_index = half_line[hash_index+1 : hash_index+3]
    if single_parity_index[1] == ")":                                           # in case the index is only a single digit, ignore the ")" that will have been caught
        single_parity_index = single_parity_index[0]
    fermi_indices.append(int(single_parity_index))
    
    
    # generate a string that contains the number of orbitals and their indices, in the correct format for input to asyrmo
    nu = int(input_settings["nu"])
    first_index = int(single_parity_index) - nu//2
    last_index = int(single_parity_index) + nu//2
    if nu%2 == 1:                                                               # then we need to add one to the final index to ensure the correct number are included
        last_index += 1
    orbitals = np.r_[first_index:last_index]
    
    orbitals_string = fermi_parities[-1] + input_settings["nu"]                 # e.g. orbitals_string = "+11 19 20 21 22 23 24 25 26 27 28 29")
    for l in orbitals:
        orbitals_string += " "
        orbitals_string += str(l)
    asyrmo_orbitals.append(orbitals_string)

print("finished reading %d files\n" % len(asyrmo_orbitals))




#%%
''' WRITING ASYRMO.DAT FILE '''
# for each deformation, writes a .DAT file for ASYRMO using the file tag to name, and the provided example as a base


print("writing ASYRMO.DAT files...")

# convert the nantj, noutj, ipout inputs to the correct format
input_settings["nantj"] = input_settings["nantj"].replace(",", " ")
input_settings["noutj"] = input_settings["noutj"].replace(",", " ")
input_settings["ipout"] = input_settings["ipout"].replace(",", " ")


for f in range(len(written_file_tags)):
    
    input_settings["current_orbitals"] = asyrmo_orbitals[f]
    
    file = written_file_tags[f]
    input_settings["current_f016"] = "f016_"+file+".dat"
    input_settings["current_f017"] = "f017_"+file+".dat"
    input_settings["current_f018"] = "f018_"+file+".dat"
    
    try:                                                                        # only overwrite the existing input file if this code is successful
        new_input_text =     '''

'%(current_f016)s' '%(current_f017)s' '%(current_f018)s' FILE16,FILE17,FILE18
1,0                                        IPKT,ISKIP
%(istrch)s,%(irec)s                                        ISTRCH,IREC
%(vmi)s,4,4,8,0.0188,100.00                      VMI,NMIN,NMAX,IARCUT,A00,STIFF
%(Z)s,%(A)s,%(imin)s,%(ispin)s,%(kmax)s,%(e2plus)s,%(e2plur)s                   Z,AA,IMIN,ISPIN,KMAX,E2PLUS,E2PLUR
19.2,7.4,15,%(chsi)s,%(eta)s                     GN0,GN1,IPAIR,CHSI,ETA
%(current_orbitals)s  
  %(nantj)s  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
  %(noutj)s  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
  %(ipout)s  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  IPOUT(I)

    ''' % input_settings
    
    except KeyError:
        print("Could not write ASY_"+file+".DAT because an input is missing."+
              " Check config file is correct (and saved!) Will not attempt to overwrite existing file.")
        raise
        
    else: 
        asyrmo_dat_file_path = "../Inputs/ASY_"+file+".DAT" 
        asyrmo_dat_file = open(asyrmo_dat_file_path, 'w')
        asyrmo_dat_file.write(new_input_text)
        asyrmo_dat_file.close()     
    

print("finished writing %d ASY.DAT files" % (f+1))




#%%
''' WRITE AND RUN BASH SCRIPT TO EXECUTE ASYRMO '''
# writes a bash shell script to run each of the ASYRMO .DAT files written above
# after each one is run, the output file ASYRMO.out is copied to a new text file with a more descriptive name (based on file tag),
# so that the data is not lost when the next iteration runs asyrmo again and ASYRMO.out is overwritten


shell_script_file_path = "../RunASYRMO.sh"                                      # this is where the shell script will be created

new_shell_script_text = ""
for file in written_file_tags :
    
    new_gampn_out_file_name = "ASY_"+file+".OUT"
    
    new_shell_script_text += ("\n./../../../Executables/MO/asyrmo < ../Inputs/ASY_"+file+".DAT")
    new_shell_script_text += ("\ncp ASYRMO.out "+new_gampn_out_file_name)        # copy ASYRMO.out to a new txt file with a more descriptive name

new_shell_script_text += "\n\necho message from terminal: finished running asyrmo"

shell_script_file = open(shell_script_file_path, 'w')
shell_script_file.write(new_shell_script_text)
shell_script_file.close()                                                       

asyrmo = subprocess.Popen(["sh", "./../RunASYRMO.sh"])
asyrmo.wait(allowed_time)

new_timer_lapse = time.time()
print("finished running asyrmo in time = %.2f seconds" % (new_timer_lapse-timer_lapse))
timer_lapse = new_timer_lapse




#%%
''' WRITING PROBAMO.DAT FILE '''
# for each deformation, writes a .DAT file for PROBAMO, using the file tag to name, and the provided example as a base

for file in written_file_tags:
    
    input_settings["current_f017"] = "f017_"+file+".dat"
    input_settings["current_f018"] = "f018_"+file+".dat"

    try:
        new_input_text =     '''

'%(current_f017)s' '%(current_f018)s'              FILE17,FILE18
1,0                               ipkt,iskip
%(istrch)s                                 isrtch
%(Z)s,%(A)s                           Z,AA
0,%(cutoff)s,1,0.75,-1                ISPEC,CUTOFF,IQ,GSFAC,GR
0.0000, 0.000,0.000, 0.000        BS2,BS4 (FOR S-STATE), BS2,BS4(P-STATE)

    ''' % input_settings
    
    except KeyError:
        print("Could not write PROB_"+file+".DAT because an input is missing."+
              " Check config file is correct (and saved!) Will not attempt to overwrite existing file.")
        raise
        
    else: 
        probamo_dat_file_path = "../Inputs/PROB_"+file+".DAT" 
        probamo_dat_file = open(probamo_dat_file_path, 'w')
        probamo_dat_file.write(new_input_text)
        probamo_dat_file.close()     
    

print("finished writing PROB.DAT files")






#%%
''' WRITE AND RUN BASH SCRIPT TO EXECUTE PROBAMO '''
# writes a bash shell script to run each of the PROBAMO .DAT files written above
# after each one is run, the output file PROBAMO.out is copied to a new text file with a more descriptive name (based on file tag),
# so that the data is not lost when the next iteration runs asyrmo again and PROBAMO.out is overwritten


shell_script_file_path = "../RunPROBAMO.sh"                                     # this is where the shell script will be created
new_shell_script_text = ""

for file in written_file_tags :
    
    new_gampn_out_file_name = "PROB_"+file+".OUT"
    
    new_shell_script_text += ("\n./../../../Executables/MO/probamo < ../Inputs/PROB_"+file+".DAT")
    new_shell_script_text += ("\ncp PROBAMO.out "+new_gampn_out_file_name)        # copy ASYRMO.out to a new txt file with a more descriptive name

new_shell_script_text += "\n\necho messgae from terminal: finished running probamo"

shell_script_file = open(shell_script_file_path, 'w')
shell_script_file.write(new_shell_script_text)
shell_script_file.close()                                                       

probamo = subprocess.Popen(["sh", "./../RunPROBAMO.sh"])
probamo.wait(allowed_time)

new_timer_lapse = time.time()
print("finished running probamo in time = %.2f seconds" 
      % (new_timer_lapse-timer_lapse))
timer_lapse = new_timer_lapse

#%%

''' READ PROBAMO.OUT FILE '''
# for each deformation (i.e. each PROBAMO.OUT file):
# the row containing the ground state is identified by energy=0.0
# the ground state spin and magnetic dipole moment are recorded
# the lowest energy (yrast) state with the same spin as the experimental first excited state is identified
# the "first excitation energy" is recorded
# similarly for the second and third experimental excited states

print("\nreading PROBAMO.OUT files...")

gs_mag_mom = []
gs_spins = []
fx_energies = []
sx_energies = []
tx_energies = []


for file in written_file_tags :

    probamo_out_file_name = "PROB_"+file+".OUT"
    
    # locate the yrast (lowest energy) state of each spin
    fx_energy = math.inf                                                        # initialise to a large value so that the first '<' comparison will be true
    sx_energy = math.inf
    tx_energy = math.inf
    
    probamo_out_file = open(probamo_out_file_name, 'r')
    for line in probamo_out_file:
        line = line.strip()
        
        try:                                                                    # determine whether this line is a data row of the table (if a '-' is present, then it is)
            dash_index = line.index(" - ")
            spin = line[dash_index-3:dash_index]                                # get the spin of this row
        except ValueError:
            continue                                                            # if this row doesn't contain data, ignore it
        
        if line[0:3] == "0.0":                                                  # then this is the ground state
            gs_mag_mom.append(line[-8:].strip())                                # get the ground state magnetic moment
            gs_spins.append(spin)                                   
        
        elif spin == input_settings["fx_spin"]:                                 # check to see if spin matches experimental first excited state
            this_energy = float(line[0:5].strip())
            if this_energy < fx_energy:                                         # this will eventually locate the yrast state of this spin (i.e. the first excited state)
                fx_energy = this_energy
                fx_line = line
        elif spin == input_settings["sx_spin"]:                                 # repeat for second excited state
            this_energy = float(line[0:5].strip())
            if this_energy < sx_energy:
                sx_energy = this_energy
                sx_line = line
        elif spin == input_settings["tx_spin"]:                                 # repeat for third excited state
            this_energy = float(line[0:5].strip())
            if this_energy < tx_energy:
                tx_energy = this_energy
                tx_line = line
    
    try:
        fx_energies.append(float(fx_line[0:5].strip()))
        sx_energies.append(float(sx_line[0:5].strip()))
        tx_energies.append(float(sx_line[0:5].strip()))
    except NameError:
        print("could not find the excited states (no states of the appropriate spins)")
        # raise ValueError("could not find the excited states (no states of the appropriate spins)")
        
    probamo_out_file.close()

# convert collected data from string to float
gs_mag_mom = [float(n) for n in gs_mag_mom]
gs_spins = [float(n[0])/2 for n in gs_spins]

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


data_matrix = [gs_mag_mom, fermi_energies, fermi_indices, 
               gs_spins, fx_energies, sx_energies, tx_energies]

graphs_to_print = ["Ground State Magnetic Dipole Moment", "Fermi Energy", "Fermi Level Parity And Index", 
                   "Ground State Spin", "First Excitation Energy", "Second Excitation Energy"]

fermi_index_colour_levels = np.arange(min(fermi_indices)-0.5, max(fermi_indices)+1.5, 1.0) 
gs_spin_colour_levels = np.arange(min(gs_spins)-0.5, max(gs_spins)+1.5, 1.0)
contour_levels = [8, 10, fermi_index_colour_levels, 
                  gs_spin_colour_levels, 10, 10]

fermi_index_cbar_ticks = np.arange(min(fermi_indices), max(fermi_indices)+1.0, 1.0)
fermi_index_cbar_ticks = [str(int(n)) for n in fermi_index_cbar_ticks]
gs_spin_cbar_ticks = np.arange(min(gs_spins), max(gs_spins)+1.0, 1.0)
gs_spin_cbar_ticks = [(str(int(n*2))+"/2") for n in gs_spin_cbar_ticks]
cbar_ticks = [0,0, fermi_index_cbar_ticks, 
              gs_spin_cbar_ticks,0,0]

data_axis_labels = [r'μ / $μ_{N}$', 'fermi energy / MeV', 'fermi level parity and index',
                    'ground state spin I', "First Excitation Energy / keV", "Second Excitation Energy / keV"]

experimental_data = [[], [], [], 
                     [], [], []]
if "gs_mu" in input_settings:
    experimental_data[0] = float(input_settings["gs_mu"])
if "gs_spin" in input_settings:
    experimental_data[3] = float(input_settings["gs_spin"])
if "fx_energy" in input_settings:
    experimental_data[4] = float(input_settings["fx_energy"])
if "sx_energy" in input_settings:
    experimental_data[5] = float(input_settings["sx_energy"])
    
error_tolerance = [0.2, 0.0, 0.0, 0.1, 50, 50]                                #!!! these are a bit arbitrary...


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

if "eps_max" in input_settings:   
    
    
    agreed_points = [] # an array to store a list of points that have properties in agreement with the experimental data
    
    for g in range(len(graphs_to_print)):
        
        input_settings["current_graph"] = graphs_to_print[g]
        print("\nplotting graph of %(current_graph)s variation..." % input_settings)
        
        fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
        
        ax.set_thetamin(0)   
        ax.set_thetamax(60)  
        ax.set_rmax(1.0)
        
        theta_ticks = np.arange(0, 70, 10)  
        ax.set_xticks(np.radians(theta_ticks))  
        
        # create filled colour contour plots
        c_level_boundaries = contour_levels[g]
        cax = ax.tricontourf(gamma_points, eps_points, data_matrix[g], levels=c_level_boundaries)
        cbar = plt.colorbar(cax, pad=0.1, label=data_axis_labels[g])
        if isinstance(c_level_boundaries, list):                                # then format for discrete values
            cax = ax.tricontourf(gamma_points, eps_points, data_matrix[g], levels=c_level_boundaries)
            cbar = plt.colorbar(cax, pad=0.1, label=data_axis_labels[g])
    
            ticks = [(c_level_boundaries[i] + c_level_boundaries[i+1]) / 2 
                     for i in range(len(c_level_boundaries) - 1)]               # Calculate midpoints of levels for tick placement
            cbar.set_ticks(ticks)
            cbar.set_ticklabels(cbar_ticks[g])
        
        ax.tricontour(gamma_points, eps_points, data_matrix[3], 
                      levels=[float(input_settings["gs_spin"])-0.5,float(input_settings["gs_spin"])+0.5],  
                      colors=[(213/255,1,0)], linewidths=1.0)
        
        
        # create flags to record which ROIs should be included in the legend (each flag will only be turned on if that ROI is plotted)
        dot_hit_flag = False
        plus_hit_flag = False
        dot_miss_flag = False
        plus_miss_flag = False
        dot_flag = False
        dot_flag = False
        
        # plot the data point markers, and check which deformations give a result that matches experiment
        for r in range(len(gamma_points)):
            if experimental_data[g]:                                            # if the experimental data exists for comparison
                error = abs(data_matrix[g][r] - experimental_data[g])
                if error < error_tolerance[g]: # error/experimental_data[g] < 0.1: #!!! if they agree within 10%, plot the data point in red rather than white
                    if g==0 and fermi_parities[r]==input_settings["par"]:       # check that the ground state parity has been reproduced
                        agreed_points.append(r)
                    if fermi_parities[r] == "-":
                        plt.polar(gamma_points[r], eps_points[r], 'r.')
                        dot_hit_flag = True
   
                    elif fermi_parities[r] == "+":
                        plt.polar(gamma_points[r], eps_points[r], 'r+')
                        plus_hit_flag = True
                else:
                    if r in agreed_points:
                        del agreed_points[agreed_points.index(r)]               # this point is no longer in agreement, so remove it from the list
                    if fermi_parities[r] == "-":
                        plt.polar(gamma_points[r], eps_points[r], 'w.')
                        dot_miss_flag = True
                    else: # fermi_parities[r] == "+":
                        plt.polar(gamma_points[r], eps_points[r], 'w+')
                        plus_miss_flag = True
                exp, = cbar.ax.plot([0, 1], [experimental_data[g], experimental_data[g]], 'r-')
                        
            # if there is no experimental data available for comparison, just plot all points in white            
            elif fermi_parities[r] == "-":
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
        exp, = plt.polar(3.0, 0.2, 'r-', label="experimental value")
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
        
        
        annotation_legend = ax.legend(handles=[spin, exp],  loc="upper left", bbox_to_anchor=(-0.3, 1.12))
        annotation_legend.set_title('annotations\n----------------------')
        
        
        legend = ax.legend(handles=legend_handles, loc="upper left", bbox_to_anchor=(-0.3, 0.78))
        legend.set_title("data points\n---------------")
        
        plt.gca().add_artist(annotation_legend)
        
        
        # add title and axis labels
        ax.set_title('%(current_graph)s of %(nucleus)s' % input_settings, va='bottom', y=1.1)                      
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
        












elif "line" in input_settings:                                                  # then plot a line graph of data variation with eps (or gamma)

    # locate the range in which the excitation energy data is meaningful (the range in which the ground state spin is correct)
    correct_spin_range = []                                                     #!!! do this with a mask instead? would allow easy combination of conditions (parity) by multiplication...
    start_flag = False
    for i in range(len(gs_spins)):
        if gs_spins[i] == experimental_data[3] and not start_flag:              # correct, and range hasn't started yet
            start_flag = True
            correct_spin_range.append(i)
        elif gs_spins[i] != experimental_data[3] and start_flag:                # incorrect, and range has started
            start_flag = False
            correct_spin_range.append(i)
    
    # plot graphs
    for g in range(len(graphs_to_print)):
        input_settings["current_graph"] = graphs_to_print[g]
        print("plotting graph of %(current_graph)s variation..." % input_settings)
        
        fig, ax = plt.subplots()                                                # create figure and axes
        
        # if eps was varied and gamma was fixed:
        if input_settings["line"] == "eps" or input_settings["line"] == "ε":
            input_settings["line"] = "ε"
            input_settings["fixed"] = "γ / º"
            # print("plotting line graph of variation with eps...")
            
            # if experimental data is available, plot it in red for easy comparison
            if experimental_data[g]:
                exp, = plt.plot(eps_to_test, 
                       np.full(len(eps_to_test), float(experimental_data[g])), 
                       'r-', label="experimental value")
            
            # mark the range in which the correct ground state spin was calculated
            for r in range(len(correct_spin_range)):
                if correct_spin_range[r] == 0:                                  # the first value of eps in the range has the correct spin
                    start_range = (eps_to_test[correct_spin_range[r]]-
                                   float(input_settings["eps"][2])/2)          
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
                                     float(input_settings["eps"][2])/2)
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
            data, = plt.plot(eps_to_test, data_matrix[g], 'k-x', label="γ = %s" % gamma_to_test[0])

        # if gamma was varied and eps was fixed:   
        else:                                                                   # input_settings["line"] == "gamma":
            input_settings["line"] = "γ / º"
            input_settings["fixed"] = "ε"
            # print("plotting line graph of variation with gamma...")
            
            # if experimental data is available, plot it in red for easy comparison
            if experimental_data[g]:
                exp, = plt.plot(gamma_to_test, 
                       np.full(len(gamma_to_test), float(experimental_data[g])), 
                       'r-', label="experimental value")
            
            # mark the range in which the correct ground state spin was calculated
            for r in range(len(correct_spin_range)):
                if correct_spin_range[r] == 0:                                  # the first value of gamma in the range has the correct spin
                    start_range = (gamma_to_test[correct_spin_range[r]]-
                                   float(input_settings["gamma"][2])/2)
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
                                     float(input_settings["gamma"][2])/2)
                        
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
            data, = plt.plot(gamma_to_test, data_matrix[g], 'k-', label="ε = %s" % eps_to_test[0])
            for p in range(len(gamma_to_test)):
                if fermi_parities[p] == "-":
                    dot, = plt.polar(gamma_to_test[p], data_matrix[g][p], 'k.', label='negative parity')

                elif fermi_parities[p] == "+":
                    plus, = plt.polar(gamma_to_test[p], data_matrix[g][p], 'k+', label='positive parity')
                    
            
        
        
        # add title and axis labels
        ax.set_title('Variation of %(current_graph)s With %(line)s in %(nucleus)s' 
                     % input_settings, va='bottom', y=1.1)                     
        plt.xlabel("%(line)s" % input_settings)
        plt.ylabel(data_axis_labels[g])
        
        # add legend (depending on what ROIs have been drawn)
        if correct_spin_range: legend = ax.legend(handles=[exp, correct_spin, data, dot, plus])
        else:
            legend = ax.legend(handles=[exp, data, dot, plus])
            print("none of these deformations yeilded a correct ground state spin")
        
        plt.show()

# note how long it took
new_timer_lapse = time.time()
print("finished plotting graphs in time = %.2f seconds" % (new_timer_lapse-timer_lapse))
print("total runtime = %.2f seconds" % (new_timer_lapse-start_wall_time))






