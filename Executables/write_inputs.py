#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  8 10:07:25 2024

@author: katyrr


run from Outputs directory in terminal, with command: 
    python3 ../../../Executables/write_inputs.py

or run from console with: 
    runfile('/Users/katyrr/Downloads/MSci Project/Code/Executables/write_inputs.py', wdir='/Users/katyrr/Downloads/MSci Project/Code/Ir191_test_new_inputs/MO/Outputs')

to debug (use breakpoints) run from console with: 
    debugfile('/Users/katyrr/Downloads/MSci Project/Code/Executables/write_inputs.py', wdir='/Users/katyrr/Downloads/MSci Project/Code/Ir191_test_new_inputs/MO/Outputs')

"""




import numpy as np                                                              # for np.arange(start, end, step)
import subprocess                                                               # for calling shell scripts to run 
import math                                                                     # for ceil(num)
import matplotlib.pyplot as plt                                                 # for plotting graphs
#from mpl_toolkits.axes_grid1 import make_axes_locatable                         # for adjusting the position of the colour bar on the plots

config_file_name = "../Inputs/config.txt"
                                                                 







''' OPEN AND READ CONFIG FILE '''
# ignores empty lines, and lines beginning with * (to mark a comment)
# counts and outputs the number of non-empty and non-comment lines (should be equal to the number of variables input in the config file)
# expects each line to have the format:    var_name value
# outputs a warning if the format is unexpected, otherwise splits the string at " " and assigns the variable to a dictionary
# tests whether eps has been input as a range or as a single value, and creates a list of values to test (which may contain just one value), for easy input to a for loop later
# similarly for gamma


config_file = open(config_file_name, 'r')
print("reading from config file...\n")


line_count = 0                                                                  # start a counter for lines
bad_lines = []                                                                  # create an empty list to store indices of bad lines

input_settings = {}                                                             # create an empty dictionary to hold settings for input


for line in config_file:
    
    line_string = line.strip()                                                  # extract text from line
    if line_string == "" : continue                                             # skip blank rows
    if line_string[0] == "*" : continue                                         # skip comment lines
    
    split_string = line_string.split(" ")                                       # split into name and value (and potentially comments starting with *)
        
    line_count += 1                                                             # increase line count now so that lines will be indexed from 1
    
    for n in range(len(split_string)):
        word = split_string[n]
        if word[0] == '*':
            del split_string[n:]
            break                                                # remove comments that start at the end of the line
    
    if len(split_string)>2 :                                                    # check whether the line has been formatted asa expected
        print("Line "+str(line_count)+ ''' has too many words! 
              Check var name contains no spaces, and don't seperate list 
              values with spaces. Ignoring extra words.\n''')
        bad_lines.append(line_count)
        
    elif len(split_string)<2 :
        print("Line "+str(line_count)+''' is missing either 
              variable name or value. Check they are both present 
              and seperated with a space. Skipping this line.\n''')
        bad_lines.append(line_count)
        
    input_settings[split_string[0]] = split_string[1]                           # store in dictionary


input_settings["eps"] = input_settings["eps"].split(",") 
if len(input_settings["eps"]) > 1 :                                             # then a range of eps values has been input
    eps_to_test = np.linspace(float(input_settings["eps"][0]),                  # create a list of 13 eps values to be tested
                              float(input_settings["eps"][1]), 
                              13)  

else : eps_to_test = [float(input_settings["eps"][0])]                          # even though it only has one element, a list is easier to input to a for loop later

eps_points = []
gamma_points = []                                                              # arrange a set of gamma values such that the eps/gamma plane is covered with an even distribution of 91 grid points
for l in range(len(eps_to_test)):
    for p in range(l+1):
        eps = eps_to_test[l]
        eps_points.append(eps)
        if l==0:
            gamma = 0
        else:
            gamma = p*(60/l)*np.pi/180
        gamma_points.append(gamma)
        



print(input_settings)
print("\n\nFinished reading lines: "+str(line_count))
print(str(len(bad_lines))+" of these had unexpected formatting.\n\n")













''' WRITE GAMPN.DAT FILES USING DICT '''
# checks which mode to run (MO or WS)
# sets up a for loop through the list of eps values to test, nested inside a similar loop over gamma 
# each iteration:
#   creates a tag referencing the code being run, the nucleus being tested, and the value of eps being tested in this iteration
#   creates/overwrites a .DAT file in the Inputs directory folder, named using the tag
#   uses string formatting with the dictionary to write the .DAT file 
#   (currently uses the provided example as a base, and changes as few settings as possible!)
# if anything went wrong, the existing .DAT file won't be overwritten 

print("writing GAMPN.DAT files...")

write_count = 0                                                                 # start counting how many files get written
written_file_tags = []                                                          # create a list to record which files were written

if (input_settings["mode"] == "MO") :
    
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
    
    for gamma in gamma_points:
        input_settings["current_gamma"] = gamma

        for eps in eps_points :
            input_settings["current_eps"] = eps                                     # store eps in the dict for ease of use when string formatting below
            
            #file_tag = "%s_eps_%.3f_gamma_%d" % (input_settings["nucleus"], eps, gamma)
            file_tag = "e%.3f_g%d" % (eps, gamma)
            
            input_settings["current_f002"] = "f002_"+file_tag+".dat"
            input_settings["current_f016"] = "f016_"+file_tag+".dat"
            input_settings["current_f017"] = "f017_"+file_tag+".dat"
            
            try:                                                                    # only overwrite the existing input file if this code is successful
                new_input_text =     '''
        
'%(current_f002)s' '%(current_f016)s' '%(current_f017)s'     file2, file16, file17
1,0,1                          ISTRCH,ICORR,irec
9,1                            NKAMY,NNEUPR
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
%(iparp)s%(norbitp)s%(porbs)s                  IPARP, NORBITP, LEVELP
%(iparn)s%(norbitn)s%(norbs)s                  IPARN, NORBITN, LEVELN
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

'''
print("waiting 5 seconds for gampn to finish running...")
time.sleep(5)
'''










''' CALCULATE FERMI LEVEL '''
# using input A and Z
# work out which is odd
# half and ceiling for the overall index of the fermi level orbital


input_settings["A"] = int(input_settings["A"])                                  # convert A and Z from strings to ints
input_settings["Z"] = int(input_settings["Z"])
input_settings["N"] = input_settings["A"]-input_settings["Z"]                   # calculate N = A - Z

if input_settings["A"]%2 == 0:
    print("A is even, but this program applies to odd-A nuclei. Check inputs! Proceeding with calculation for protons...")
    input_settings["P_or_N"] = "P"
    input_settings["fermi_level"] = math.ceil(input_settings["Z"]/2)            # calculate fermi level
    
elif input_settings["Z"]%2  == 1:    
    print("Z is odd; calculating for protons...")
    input_settings["P_or_N"] = "P"
    input_settings["fermi_level"] = math.ceil(input_settings["Z"]/2)
else:
    print("N is odd; calculating for neutrons...")
    input_settings["P_or_N"] = "N"
    input_settings["fermi_level"] = math.ceil(input_settings["N"]/2)

                                                       










''' READ GAMPN.OUT FILE '''
# for each deformation (i.e. each GAMPN.OUT file):
#   work out the line number of the fermi level orbital in the GAMPN.OUT file
#   read that line from GAMPN.OUT and get the energy, parity, and orbital index restricted to that parity
#   construct a string for ipar,nu,[level(j),j=1,nu]

print("\nreading GAMPN.OUT files...")

asyrmo_inputs      = []                                                         # an empty array to store inputs for asyrmo.dat
fermi_energies     = []                                                         # an empty array to store the fermi energy at each deformation

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
    ipar = half_line[hash_index-2]
    single_parity_index = half_line[hash_index+1 : hash_index+3]
    fermi_energy = half_line[hash_index-10 : hash_index-4]
    fermi_energies.append(fermi_energy)
    
    if single_parity_index[1] == ")":                                           # in case the index is only a single digit
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








''' WRITING ASYRMO.DAT FILE '''
# for each deformation, writes a .DAT file for ASYRMO using the file tag to name, and the provided example as a base

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
1,1                                        ISTRCH,IREC
0,4,4,8,0.0188,100.00                      VMI,NMIN,NMAX,IARCUT,A00,STIFF
%(Z)s,%(A)s,%(imin)s,%(ispin)s,%(kmax)s,0.210,0.0                   Z,AA,IMIN,ISPIN,KMAX,E2PLUS,E2PLUR
19.2,7.4,15,0.0,1.0                     GN0,GN1,IPAIR,CHSI,ETA
%(current_orbitals)s  
  3  3  3  3  3  3  3  0  0  0  0  0  0  0  0  0  0  0  0  0  0
  2  2  2  2  2  2  2  0  0  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  1  1  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  IPOUT(I)

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












''' WRITING PROBAMO.DAT FILE '''
# for each deformation, writes a .DAT file for PROBAMO, using the file tag to name, and the provided example as a base

for file in written_file_tags:
    
    input_settings["current_f017"] = "f017_"+file+".dat"
    input_settings["current_f018"] = "f018_"+file+".dat"

    try:                                                                        # only overwrite the existing input file if this code is successful
        new_input_text =     '''

'%(current_f017)s' '%(current_f018)s'              FILE17,FILE18
1,0                               ipkt,iskip
1                                 isrtch
%(Z)s,%(A)s                           Z,AA
0,1500.,1,0.75,-1                ISPEC,CUTOFF,IQ,GSFAC,GR
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








''' READ ASYRMO.OUT FILE '''
# for each deformation (i.e. each GAMPN.OUT file):
#   work out the line number of the fermi level orbital in the GAMPN.OUT file
#   read that line from GAMPN.OUT and get the energy, parity, and orbital index restricted to that parity
#   construct a string for ipar,nu,[level(j),j=1,nu]

print("\nreading PROBAMO.OUT files...")

mag_moments     = []                                                            # an empty array to store the magnetic dipole moment at each deformation

for file in written_file_tags :
    
    probamo_out_file_name = "PROB_"+file+".OUT"
    probamo_out_file = open(probamo_out_file_name, 'r')
    
    lines = probamo_out_file.readlines()
    header = lines.index("       E1     II      EF     IF     B(E2)             B(M1)              GD        D       EREL    T(E2)      T(M1)\n")
    
    spin1_line_index = header+3                                                 # calculate the line number of the spin 1/2 state internal transition
    spin1_line = lines[spin1_line_index]
    mag_moments.append(spin1_line[52:60].strip())                               # get only the first half of the line
    

    

print("finished reading %d files\n" % len(asyrmo_inputs))






''' PLOT GRAPHS'''
# plots variation in magnetic dipole moment with distortion
#     converts oblate input values (with gamma between 0-30º and eps<0) to the equivalent parameterisation: 
#     eps -> -eps (now positive) and gamma -> 60º-gamma (now in range 30-60º)
#     which allows all of the values to be plotted in polar coordinates.
#     
#     additionally plots the positions of the data points as white x marks.



print("\nplotting graph of magnetic dipole moment variation...")


gamma_points = np.array(gamma_points)
eps_points = np.array(eps_points)
#fermi_energies = [float(n) for n in fermi_energies]
mag_moments = [float(n) for n in mag_moments]


gamma_to_plot, eps_to_plot = np.meshgrid(gamma_points, eps_points)

fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))

ax.set_thetamin(0)   # Start angle in degrees
ax.set_thetamax(60)  # End angle in degrees

theta_ticks = np.arange(0, 70, 10)  # Create tick locations in degrees
ax.set_xticks(np.radians(theta_ticks))  # Convert to radians for set_xticks

'''
# (eps,gamma) and (-eps,60-gamma) correspnod to the same shape - covert to positive eps so that it can be plotted in polar coordinates
for r in range(len(gamma_to_test)):
    for c in range(len(eps_to_test)):
        if(eps_to_plot[c][r]<0):
            eps_to_plot[c][r] *= -1
            gamma_to_plot[c][r] = 60-gamma_to_plot[c][r]
        gamma_to_plot[c][r] *= np.pi/180
        plt.polar(gamma_to_plot[c][r], eps_to_plot[c][r], 'wx')
            
#fermi_energies_to_plot = np.array(fermi_energies).reshape((len(eps_to_test),len(gamma_to_test)))
mag_moments_to_plot = np.array(mag_moments).reshape((len(eps_to_test),len(gamma_to_test)))

#cax = ax.contourf(gamma_to_plot, eps_to_plot, fermi_energies_to_plot, 10)
cax = ax.contourf(gamma_to_plot, eps_to_plot, mag_moments_to_plot, 10)

# Add color bar to the right of the plot
plt.colorbar(cax)

plt.show()
'''




