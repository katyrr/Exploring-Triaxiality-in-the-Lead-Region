#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  8 10:07:25 2024

@author: katyrr


run from Outputs directory in terminal, with command: 
    python3 ../../../Executables/write_inputs.py

or run from console with: 
    runfile('/Users/katyrr/Downloads/MSci Project/Code/Executables/write_inputs.py', wdir='/Users/katyrr/Downloads/MSci Project/Code/Tl207/MO/Outputs')

to debug (use breakpoints) run from console with: 
    debugfile('/Users/katyrr/Downloads/MSci Project/Code/Executables/write_inputs.py', wdir='/Users/katyrr/Downloads/MSci Project/Code/tl207/MO/Outputs')

"""

#%%


import numpy as np                                                              # for np.arange(start, end, step)
import subprocess                                                               # for calling shell scripts to run 
import math                                                                     # for ceil(num)
import matplotlib.pyplot as plt                                                 # for plotting graphs

plt.rcParams['figure.dpi'] = 150

config_file_name = "../Inputs/config.txt"
                                                                 



#%%



''' OPEN AND READ CONFIG FILE '''
# ignores empty lines, and lines beginning with * (to mark a comment)
# counts and outputs the number of non-empty and non-comment lines (should be equal to the number of variables input in the config file)
# expects each line to have the format:    var_name value
# outputs a warning if the format is unexpected, otherwise splits the string at " " and assigns the variable to a dictionary
# tests whether deformation has been input as a range or as a single value, and creates a list of values to test (which may contain just one value)


config_file = open(config_file_name, 'r')
print("reading from config file...\n")


line_count = 0                                                                  # start a counter for lines read
all_lines_count = 0                                                             # start a counter for all lines (including empty lines and comments)
bad_lines = []                                                                  # create an empty list to store indices of bad lines

input_settings = {}                                                             # create an empty dictionary to hold settings for input


for line in config_file:
    
    all_lines_count += 1
    
    line_string = line.strip()                                                  # extract text from line
    if line_string == "" : continue                                             # skip blank rows
    if line_string[0] == "*" : continue                                         # skip comment lines
    
    split_string = line_string.split(" ")                                       # split into name and value (and potentially comments starting with *)
        
    line_count += 1                                                            
    
    for n in range(len(split_string)):
        word = split_string[n]
        if word[0] == '*':                                                      # remove inline comments 
            del split_string[n:]
            break  
    
    if split_string[0]=="eps" or split_string[0]=="gamma" or split_string[0]=="single":                      # these are the only cases where 3 input words are expected
        expected_num_words = 3
    else:
        expected_num_words = 2
    
    if len(split_string)>expected_num_words :                                                    # check whether the line has been formatted as expected
        print("Line "+str(all_lines_count)+ ''' has too many words! 
              Check var name contains no spaces, and seperate list 
              values with commas not spaces. Ignoring extra words.\n''')
        bad_lines.append(line_count)
        
    elif len(split_string)<expected_num_words :
        print("Line "+str(all_lines_count)+''' is missing either 
              variable name or value. Check they are both present 
              and seperated with a space. Skipping this line.\n''')
        bad_lines.append(line_count)
    
    
    if split_string[0]=="eps" or split_string[0]=="gamma" or split_string[0]=="single": 
        if "eps" in input_settings or "gamma" in input_settings or "eps_max" in input_settings:
            print("Deformation parameters have been input multiple times. Check that config file only contains one deformation input. Ignoring excess inputs.")
            continue
        input_settings["eps"] = split_string[1]
        input_settings["gamma"] = split_string[2]
        
    elif split_string[0]=="mesh":
        if "eps" in input_settings or "gamma" in input_settings or "eps_max" in input_settings:
            print("Deformation parameters have been input multiple times. Check that config file only contains one deformation input. Ignoring excess inputs.")
            continue
        input_settings["eps_max"] = split_string[1]
        
    else:
        input_settings[split_string[0]] = split_string[1]                       # store in dictionary







if "eps" in input_settings and "gamma" in input_settings:
    
    input_settings["eps"] = input_settings["eps"].split(",") 
    if len(input_settings["eps"]) > 1 :                                             # then a range of eps values has been input
        input_settings["line"] = "eps"
        eps_to_test =  np.arange(float(input_settings["eps"][0]),  
                                float(input_settings["eps"][1])+float(input_settings["eps"][2]), # add a step to the final value (because range is non-inclusive)
                                float(input_settings["eps"][2]))
    else : eps_to_test = [float(input_settings["eps"][0])]                          # even though it only has one element, a list is easier to input to a for loop later
    
    input_settings["gamma"] = input_settings["gamma"].split(",")                    # same for gamma
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
            eps = eps_to_test[l]
            eps_points.append(eps)
            
            gamma = int(np.round(gamma_to_test[p]))
            gamma_points.append(gamma)
            
    
    print(input_settings)
    print("\n\nFinished reading lines: "+str(line_count))
    print(str(len(bad_lines))+" of these had unexpected formatting.\n\n")






    
elif "eps_max" in input_settings:                                               # then arrange a mesh of 91 evenly distributed (eps,gamma) points to test over
    
    input_settings["eps_max"] = float(input_settings["eps_max"])                # convert from string to float
    eps_to_test = np.linspace(0.001, input_settings["eps_max"], num=20)         #!!!start from eps=0.001 because eps=0 is not an allowed input when it comes to asyrmo
    
    eps_points = []
    gamma_points = []
    for l in range(len(eps_to_test)):
        for p in range(l+1):
            eps = np.round(eps_to_test[l], 3)
            eps_points.append(eps)
            if l==0:
                gamma = 0
            else:
                gamma = np.round(p*(60/l)) 
            gamma_points.append(gamma)
            
    print(input_settings)
    print("\n\nFinished reading lines: "+str(line_count))
    print(str(len(bad_lines))+" of these had unexpected formatting.\n\n")
    
            
    



#%%


''' CALCULATE FERMI LEVEL '''
# using input A and Z
# work out which is odd
# half and ceiling for the overall index of the fermi level orbital

print("calcualting fermi level...")

input_settings["A"] = int(input_settings["A"])                                  # convert A and Z from strings to ints
input_settings["Z"] = int(input_settings["Z"])
input_settings["N"] = input_settings["A"]-input_settings["Z"]                   # calculate N = A - Z

if input_settings["nneupr"] == "1":
    print("calculating for odd protons...")
    input_settings["fermi_level"] = math.ceil(input_settings["Z"]/2)            # calculate fermi level
    if input_settings["Z"]%2 == 0:
        print("this nucleus has even Z; check inputs (the code will treat this as the core and assume an odd proton on top)")
    
elif input_settings["nneupr"]  == "-1":    
    print("calculating for odd neutrons...")
    input_settings["fermi_level"] = math.ceil(input_settings["Z"]/2)
    if input_settings["Z"]%2 == 0:
        print("this nucleus has even N; check inputs (the code will treat this as the core and assume an odd neutron on top)")
    

else:
    print("missing nneupr input")
    raise ValueError


print("fermi level = " + str(input_settings["fermi_level"]) + "\n")




num_to_calc = 11

first_index = input_settings["fermi_level"]//2 - num_to_calc//2
last_index = input_settings["fermi_level"]//2 + num_to_calc//2
if num_to_calc%2 == 1:                                                               # then we need to add one to the final index to ensure the correct number are included
    last_index += 1
orbitals = np.r_[first_index:last_index]

gampn_orbitals = str(num_to_calc)                                             # e.g. orbitals_string = "4 19 20 21 22" (set the parity later)
for l in orbitals:
    gampn_orbitals += " "
    gampn_orbitals += str(l)
    
input_settings["gampn_orbitals"] = gampn_orbitals




    
    
    #%%
    
    
    
''' WRITE GAMPN.DAT FILES USING DICT '''
# checks which mode to run (MO or WS)
# sets up a for loop through the list of eps values to test, nested inside a similar loop over gamma 
# each iteration:
#   creates a tag referencing the code being run, the nucleus being tested, and the value of eps being tested in this iteration
#   creates/overwrites a .DAT file in the Inputs directory folder, named using the tag
#   uses string formatting with the dictionary to write the .DAT file 
#   (currently uses the provided example as a base, and changes as few settings as possible!)
# if anything went wrong, the existing .DAT file won't be overwritten 

print("writing GAMPN.DAT files for a mesh...")

write_count = 0                                                                 # start counting how many files get written
written_file_tags = []                                                          # create a list to record which files were written

if (input_settings["mode"] == "MO") :
    
    '''
    porbs = ""
    porb_index = int(input_settings["forbitp"])
    for o in range(int(input_settings["norbitp"])):
        porbs += (" " + str(porb_index))
        porb_index += 1
    input_settings["porbs"] = porbs
    
    norbs = ""
    norb_index = int(input_settings["forbitn"])
    for o in range(int(input_settings["norbitn"])):
        norbs += (" " + str(norb_index))
        norb_index += 1
    input_settings["norbs"] = norbs
    '''
   # %(iparp)s%(norbitp)s%(porbs)s                  IPARP, NORBITP, LEVELP
   # %(iparn)s%(norbitn)s%(norbs)s                  IPARN, NORBITN, LEVELN
    
    
    for p in range(len(gamma_points)):
        
        input_settings["current_eps"] = eps_points[p]                                   # store eps in the dict for ease of use when string formatting below
        input_settings["current_gamma"] = gamma_points[p]
        
        file_tag = "e%.3f_g%d" % (input_settings["current_eps"], input_settings["current_gamma"])
        
        input_settings["current_f002"] = "f002_"+file_tag+".dat"
        input_settings["current_f016"] = "f016_"+file_tag+".dat"
        input_settings["current_f017"] = "f017_"+file_tag+".dat"
        
        try:                                                                    # only overwrite the existing input file if this code is successful
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
%(par)s%(gampn_orbitals)s                 IPARP, NORBITP, LEVELP
%(par)s%(gampn_orbitals)s                  IPARN, NORBITN, LEVELN
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

    print("%d input files were written, \nfor eps in range [%.3f, %.3f], \nand gamma in range [%d, %d].\n" % (write_count, eps_points[0], eps_points[-1], gamma_points[0], gamma_points[-1]))


else : 
    print("No deformation parameters were input. Check config file.")
    








#%%




''' WRITE AND RUN BASH SCRIPT TO EXECUTE GAMPN '''
# writes a bash shell script to run each of the GAMPN .DAT files written above
# after each one is run, the output file GAMPN.out is copied to a new text file with a more descriptive name (based on file tag),
# so that the data is not lost when the next iteration runs gampn again and GAMPN.out is overwritten


shell_script_file_path = "../RunGAMPN.sh"                                # this is where the shell script will be created


new_shell_script_text = "echo running GAMPN ..."

for file in written_file_tags :
    
    
    new_shell_script_text += ("\n./../../../Executables/MO/gampn < ../Inputs/GAM_"+file+".DAT")
    new_shell_script_text += ("\ncp GAMPN.out GAM_"+file+".OUT")                # copy GAMPN.out to a new txt file with a more descriptive name


new_shell_script_text += "\n\necho done"

shell_script_file = open(shell_script_file_path, 'w')
shell_script_file.write(new_shell_script_text)
shell_script_file.close()                                                       

subprocess.call(["sh", "./../RunGAMPN.sh"])






                                                       





#%%




''' READ GAMPN.OUT FILE '''
# for each deformation (i.e. each GAMPN.OUT file):
#   work out the line number of the fermi level orbital in the GAMPN.OUT file
#   read that line from GAMPN.OUT and get the energy, parity, and orbital index restricted to that parity
#   construct a string for ipar,nu,[level(j),j=1,nu]

print("\nreading GAMPN.OUT files...")

asyrmo_inputs      = []                                                         # an empty array to store inputs for asyrmo.dat
fermi_energies     = []                                                         # for the fermi energy at each deformation
fermi_indices      = []                                                         # for the index and parity of the fermi level
fermi_parities     = [] 

for file in written_file_tags :
    
    gampn_out_file_name = "GAM_"+file+".OUT"
    gampn_out_file = open(gampn_out_file_name, 'r')
    
    lines = gampn_out_file.readlines()
    
    header = lines.index("   #   ENERGY +/-(#)    <Q20>    <Q22>     <R2>     <JZ>      #   ENERGY +/-(#)    <Q20>    <Q22>     <R2>     <JZ>\n")
    
        
        
    fermi_level_line = input_settings["fermi_level"]+header+1                   # calculate the line number of the fermi level in the GAMPN.OUT file (indexed from zero!)
    if input_settings["fermi_level"] > 40:
        fermi_level_line -= 40
        whole_line = lines[fermi_level_line]
        half_line = whole_line[60:-1].strip()                                   # get only the second half of the line
        
    else:
        whole_line = lines[fermi_level_line]
        half_line = whole_line[0:60].strip()                                    # get only the first half of the line
    
    hash_index = half_line.index("#")                                           # use the index of # as a reference point 
    ipar = half_line[hash_index-2]                                              #!!! rather than reading the parity out like this, just set it as a parameter? 
    single_parity_index = half_line[hash_index+1 : hash_index+3]
    fermi_energy = half_line[hash_index-10 : hash_index-4]
    fermi_energies.append(fermi_energy)
    fermi_index = int(single_parity_index)
    #if (ipar=="-"): fermi_index *= 1
    fermi_parities.append(ipar)
    fermi_indices.append(fermi_index)
    
    if single_parity_index[1] == ")":                                           # in case the index is only a single digit, ignore the ")" that will have been caught
        single_parity_index = single_parity_index[0]
    
    nu = int(input_settings["nu"])
    
    first_index = int(single_parity_index) - nu//2
    last_index = int(single_parity_index) + nu//2
    if nu%2 == 1:                                                               # then we need to add one to the final index to ensure the correct number are included
        last_index += 1
    orbitals = np.r_[first_index:last_index]
    
    orbitals_string = ipar + input_settings["nu"]                               # e.g. orbitals_string = "+4 19 20 21 22"
    for l in orbitals:
        orbitals_string += " "
        orbitals_string += str(l)
    
    asyrmo_inputs.append(orbitals_string)
    

print("finished reading %d files\n" % len(asyrmo_inputs))




#%%


''' WRITING ASYRMO.DAT FILE '''
# for each deformation, writes a .DAT file for ASYRMO using the file tag to name, and the provided example as a base

input_settings["nantj"] = input_settings["nantj"].replace(",", " ")
input_settings["noutj"] = input_settings["noutj"].replace(",", " ")
input_settings["ipout"] = input_settings["ipout"].replace(",", " ")

count = 0
for file in written_file_tags:
    
    input_settings["current_orbitals"] = asyrmo_inputs[count]
    count += 1
    
    input_settings["current_f016"] = "f016_"+file+".dat"
    input_settings["current_f017"] = "f017_"+file+".dat"
    input_settings["current_f018"] = "f018_"+file+".dat"
    
    try:                                                                    # only overwrite the existing input file if this code is successful
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
    

print("finished writing ASY.DAT files")




#%%





''' WRITE AND RUN BASH SCRIPT TO EXECUTE ASYRMO '''
# writes a bash shell script to run each of the ASYRMO .DAT files written above
# after each one is run, the output file ASYRMO.out is copied to a new text file with a more descriptive name (based on file tag),
# so that the data is not lost when the next iteration runs asyrmo again and ASYRMO.out is overwritten


shell_script_file_path = "../RunASYRMO.sh"                                # this is where the shell script will be created


new_shell_script_text = "echo running ASYRMO ..."

for file in written_file_tags :
    
    new_gampn_out_file_name = "ASY_"+file+".OUT"
    
    new_shell_script_text += ("\n./../../../Executables/MO/asyrmo < ../Inputs/ASY_"+file+".DAT")
    new_shell_script_text += ("\ncp ASYRMO.out "+new_gampn_out_file_name)        # copy ASYRMO.out to a new txt file with a more descriptive name

new_shell_script_text += "\n\necho done"

shell_script_file = open(shell_script_file_path, 'w')
shell_script_file.write(new_shell_script_text)
shell_script_file.close()                                                       

subprocess.call(["sh", "./../RunASYRMO.sh"])







#%%




''' WRITING PROBAMO.DAT FILE '''
# for each deformation, writes a .DAT file for PROBAMO, using the file tag to name, and the provided example as a base

for file in written_file_tags:
    
    input_settings["current_f017"] = "f017_"+file+".dat"
    input_settings["current_f018"] = "f018_"+file+".dat"

    try:                                                                        # only overwrite the existing input file if this code is successful
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
        asyrmo_dat_file_path = "../Inputs/PROB_"+file+".DAT" 
        asyrmo_dat_file = open(asyrmo_dat_file_path, 'w')
        asyrmo_dat_file.write(new_input_text)
        asyrmo_dat_file.close()     
    

print("finished writing PROB.DAT files")






#%%



''' WRITE AND RUN BASH SCRIPT TO EXECUTE PROBAMO '''
# writes a bash shell script to run each of the PROBAMO .DAT files written above
# after each one is run, the output file PROBAMO.out is copied to a new text file with a more descriptive name (based on file tag),
# so that the data is not lost when the next iteration runs asyrmo again and PROBAMO.out is overwritten


shell_script_file_path = "../RunPROBAMO.sh"                                # this is where the shell script will be created


new_shell_script_text = "echo running PROBAMO ..."

for file in written_file_tags :
    
    new_gampn_out_file_name = "PROB_"+file+".OUT"
    
    new_shell_script_text += ("\n./../../../Executables/MO/probamo < ../Inputs/PROB_"+file+".DAT")
    new_shell_script_text += ("\ncp PROBAMO.out "+new_gampn_out_file_name)        # copy ASYRMO.out to a new txt file with a more descriptive name

new_shell_script_text += "\n\necho done"

shell_script_file = open(shell_script_file_path, 'w')
shell_script_file.write(new_shell_script_text)
shell_script_file.close()                                                       

subprocess.call(["sh", "./../RunPROBAMO.sh"])






#%%

''' READ PROBAMO.OUT FILE '''
# for each deformation (i.e. each PROBAMO.OUT file):
#   work out the line number of the fermi level orbital in the GAMPN.OUT file
#   read that line from GAMPN.OUT and get the energy, parity, and orbital index restricted to that parity
#   construct a string for ipar,nu,[level(j),j=1,nu]

print("\nreading PROBAMO.OUT files...")

#mag_moments     = []                                                            # an empty array to store the magnetic dipole moment at each deformation
gs_mag_mom = []
gs_spins = []

for file in written_file_tags :
    
    probamo_out_file_name = "PROB_"+file+".OUT"
    probamo_out_file = open(probamo_out_file_name, 'r')
    
    lines = probamo_out_file.readlines()
    
    ''' # this is assuming that the ground state will be spin 1/2!
    try:
        header = lines.index("       E1     II      EF     IF     B(E2)             B(M1)              GD        D       EREL    T(E2)      T(M1)\n")
    except ValueError: #!!! need to actually implement a backup plan for this - read error message and return to gampn with a new input
        print("Could not read probamo output at deformation " +file +"\nProbably because File18 was not generated from asyrmo, so probamo could not be run.\nTypically caused by incorrect input of orbitals into gampn.")
        raise
    spin1_line_index = header+3                                                 # calculate the line number of the spin 1/2 state internal transition
    spin1_line = lines[spin1_line_index]
    mag_moments.append(spin1_line[52:60].strip())                               # get the magnetic moment
    '''
    
    for l in range(len(lines)):
        this_line = lines[l].strip()
        if this_line[0:3] == "0.0":
            gs_line_index = l
            break
    gs_mag_mom.append(this_line[-8:].strip())                             # get the magnetic moment
    gs_spins.append(this_line[4:10].strip())


print("finished reading %d files\n" % len(asyrmo_inputs))

#%%

eps_points = np.array(eps_points)
gamma_points = np.array(gamma_points)
gamma_points = [float(n) for n in gamma_points]
for r in range(len(gamma_points)):
    gamma_points[r] *= np.pi/180

fermi_energies = [float(n) for n in fermi_energies]
gs_mag_mom = [float(n) for n in gs_mag_mom]
gs_spins = [float(n[0])/2 for n in gs_spins]
data_matrix = [gs_mag_mom, fermi_energies, fermi_indices, gs_spins]

fermi_index_colour_levels = np.arange(min(fermi_indices)-0.5, max(fermi_indices)+1.5, 1.0) 
gs_spin_colour_levels = np.arange(min(gs_spins)-0.5, max(gs_spins)+1.5, 1.0)

fermi_index_cbar_ticks = np.arange(min(fermi_indices), max(fermi_indices)+1.0, 1.0)
fermi_index_cbar_ticks = [str(int(n)) for n in fermi_index_cbar_ticks]
gs_spin_cbar_ticks = np.arange(min(gs_spins), max(gs_spins)+1.0, 1.0)
gs_spin_cbar_ticks = [(str(int(n*2))+"/2") for n in gs_spin_cbar_ticks]

graphs_to_print = ["Magnetic Dipole Moment the Ground State", "Fermi Energy", "Fermi Level Parity And Index", "Ground State Spin"] #, "First Excitation Energy"]
# for gs Jp, find EI=0.0 in probamo.out and read spin; for first exc en, input expected Jp, find in probamo.out, and read exc en 
data_axis_labels = [r'μ / $μ_{N}$', 'fermi energy / $\hbar\omega_{0}$', 'fermi level parity and index', 'ground state spin I']


    

#%%

''' PLOT GRAPHS'''
# plots variation in magnetic dipole moment with distortion (only if distortion was input as a mesh)
#     converts oblate input values (with gamma between 0-30º and eps<0) to the equivalent parameterisation: 
#     eps -> -eps (now positive) and gamma -> 60º-gamma (now in range 30-60º)
#     which allows all of the values to be plotted in polar coordinates.
#     
#     additionally plots the positions of the data points as white x marks.

if "eps_max" in input_settings:   
    
    
    
    contour_levels = [15, 10, fermi_index_colour_levels, gs_spin_colour_levels]
    cbar_ticks = [0,0, fermi_index_cbar_ticks, gs_spin_cbar_ticks]
    
    dc = ['c','c','d','d'] # record whether each data set is continuous or discrete
    
    
    for g in range(len(graphs_to_print)):
        
        input_settings["current_graph"] = graphs_to_print[g]
        print("\nplotting graph of %(current_graph)s variation..." % input_settings)
        
        fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
        
        ax.set_thetamin(0)   # Start angle in degrees
        ax.set_thetamax(60)  # End angle in degrees
        ax.set_rmax(1.0)
        
        theta_ticks = np.arange(0, 70, 10)  # Create tick locations in degrees
        ax.set_xticks(np.radians(theta_ticks))  # Convert to radians for set_xticks
        
        for r in range(len(gamma_points)):
            if fermi_parities[r] == "-":
                dot, = plt.polar(gamma_points[r], eps_points[r], 'k.', label="parity (-)")
            elif fermi_parities[r] == "+":
                plus, = plt.polar(gamma_points[r], eps_points[r], 'k+', label="parity (+)")
        
        
        legend = ax.legend(handles=[dot, plus], loc="upper left", bbox_to_anchor=(-0.1, 1.0))
        legend.set_title("data points")
        
        handles, labels = ax.get_legend_handles_labels()
        for h in handles:
            h.set_markerfacecolor('white')  # Set marker color to white on the graph, but don't update the legend afterwards (so that the legend markers remain black for visibility)
            h.set_markeredgecolor('white')
                
        if dc[g] == 'c':
            cax = ax.tricontourf(gamma_points, eps_points, data_matrix[g], levels=contour_levels[g])#, vmin=-0.8, vmax=5.4) # vmin/max are range of colourmap, -0.8/6.4 for mu, 4.5, 6.5 for fermi energy
            plt.colorbar(cax, pad=0.1, label=data_axis_labels[g])
        else:
            ''' # this isn't much good... the interpolation is too obvious. At least with the contour plot, each point has a correct and clearly-defined colour
            c_level_boundaries = contour_levels[g]
            cax = ax.tripcolor(gamma_points, eps_points, data_matrix[g], shading='flat', cmap='viridis')
            cbar = plt.colorbar(cax, pad=0.1, label=data_axis_labels[g])
            ticks = [(c_level_boundaries[i] + c_level_boundaries[i+1]) / 2 for i in range(len(c_level_boundaries) - 1)] # Calculate midpoints of levels for tick placement
            cbar.set_ticks(ticks)#c_level_boundaries[:-1])
            cbar.set_ticklabels(cbar_ticks[g])
            '''
            c_level_boundaries = contour_levels[g]
            cax = ax.tricontourf(gamma_points, eps_points, data_matrix[g], levels=c_level_boundaries)
            cbar = plt.colorbar(cax, pad=0.1, label=data_axis_labels[g])
            ticks = [(c_level_boundaries[i] + c_level_boundaries[i+1]) / 2 for i in range(len(c_level_boundaries) - 1)] # Calculate midpoints of levels for tick placement
            cbar.set_ticks(ticks)
            cbar.set_ticklabels(cbar_ticks[g])#['22-', '21-', '', '20-', '', '19-', '', '19+', '', '20+', '', '21+', '22+'])
            
           
            
            
        ax.set_title('%(current_graph)s of %(nucleus)s' % input_settings, va='bottom', y=1.1)  # Adjust y value as needed
        plt.xlabel("ε")
        ax.text(65*np.pi/180, ax.get_rmax()*1.05, "γ", ha='center', va='center')
       
        plt.show()

elif "line" in input_settings: # then plot a line graph of parammeter variation with eps (or gamma)
    
    
    for g in range(len(graphs_to_print)):

        
        input_settings["current_graph"] = graphs_to_print[g]
        print("\nplotting graph of %(current_graph)s variation..." % input_settings)
        
        data = np.array(data_matrix[g])
        fig, ax = plt.subplots()
        
        
        if input_settings["line"] == "eps":
            input_settings["line"] = "ε"
            input_settings["fixed"] = "γ / º"
            print("plotting line graph of variation with eps...")
            plt.plot(eps_to_test, data, '-x', label=gamma_to_test)
            
        else: # input_settings["line"] == "gamma":
            input_settings["line"] = "γ / º"
            input_settings["fixed"] = "ε"
            print("plotting line graph of variation with gamma...")
            plt.plot(gamma_to_test, data, '-x', label=eps_to_test)
        
        ax.set_title('variation of %(current_graph)s with %(line)s in %(nucleus)s' % input_settings, va='bottom', y=1.1)  # Adjust y value as needed
        plt.xlabel("%(line)s" % input_settings)
        plt.ylabel(data_axis_labels[g])
        
        legend = ax.legend()
        legend.set_title("%(fixed)s" % input_settings)
       
        plt.show()
