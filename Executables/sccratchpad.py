#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 17:45:40 2024

@author: katyrr
"""




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
        eps_to_test =  np.arange(float(input_settings["eps"][0]),  
                                float(input_settings["eps"][1])+float(input_settings["eps"][2]), # add a step to the final value (because range is non-inclusive)
                                float(input_settings["eps"][2]))
    else : eps_to_test = [float(input_settings["eps"][0])]                          # even though it only has one element, a list is easier to input to a for loop later
    
    input_settings["gamma"] = input_settings["gamma"].split(",")                    # same for gamma
    if len(input_settings["gamma"]) > 1 :                                             
        gamma_to_test = np.arange(float(input_settings["gamma"][0]),  
                                float(input_settings["gamma"][1])+float(input_settings["gamma"][2]), 
                                float(input_settings["gamma"][2]))
    else : gamma_to_test = [float(input_settings["gamma"][0])] 
    
    
    eps_points = []
    gamma_points = []
    for l in range(len(eps_to_test)):
        for p in range(len(gamma_to_test)):
            eps = np.round(eps_to_test[l], 3)
            eps_points.append(eps)
            
            gamma = int(np.round(gamma_to_test[p]))  #*np.pi/180 #!!!don't convert to radians yet 
            gamma_points.append(gamma)
            
    
    print(input_settings)
    print("\n\nFinished reading lines: "+str(line_count))
    print(str(len(bad_lines))+" of these had unexpected formatting.\n\n")






    
elif "eps_max" in input_settings:                                               # then arrange a mesh of 91 evenly distributed (eps,gamma) points to test over
    
    input_settings["eps_max"] = float(input_settings["eps_max"])                # convert from string to float
    eps_to_test = np.linspace(0.001, input_settings["eps_max"], num=6)          # start from eps=0.001 because eps=0 is not an allowed input when it comes to asyrmo
    
    eps_points = []
    gamma_points = []
    for l in range(len(eps_to_test)):
        for p in range(l+1):
            eps = np.round(eps_to_test[l], 3)
            eps_points.append(eps)
            if l==0:
                gamma = 0
            else:
                gamma = int(np.round(p*(60/l)))  #*np.pi/180 #!!!don't convert to radians yet 
            gamma_points.append(gamma)
            
    print(input_settings)
    print("\n\nFinished reading lines: "+str(line_count))
    print(str(len(bad_lines))+" of these had unexpected formatting.\n\n")
    
            
    
    
    
    
    
    
    
    
''' WRITE GAMPN.DAT FILES USING DICT (CASE 2: mesh) '''
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


else : 
    print("No deformation parameters were input. Check config file.")
    













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
import numpy as np
import matplotlib.pyplot as plt


#eps_to_test = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2]
#eps_to_test = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2]
#eps_to_test = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]
eps_min = 0.0
eps_max = 0.83
eps_to_test = np.linspace(eps_min, eps_max, num=6)

eps_to_plot = [] 
gamma_to_plot = []
fermi_energies_to_plot = [] 

fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))

ax.set_thetamin(0)   # Start angle in degrees
ax.set_thetamax(60)  # End angle in degrees

theta_ticks = np.arange(0, 70, 10)  # Create tick locations in degrees
ax.set_xticks(np.radians(theta_ticks))  # Convert to radians for set_xticks

for l in range(len(eps_to_test)):
    for p in range(l+1):
        eps = eps_to_test[l]
        eps_to_plot.append(eps)
        if l==0:
            gamma = 0
        else:
            gamma = p*(60/l)*np.pi/180
        
        gamma_to_plot.append(gamma)
        fermi_energies_to_plot.append(eps)
        plt.polar(gamma, eps, 'kx')

        

cax = ax.tricontourf(gamma_to_plot, eps_to_plot, fermi_energies_to_plot, levels=10)

plt.colorbar(cax)
plt.show()
     
'''


'''
string_list = ["1.1", "2.2", "3.3"]
float_list = [float(num) for num in string_list]

gamma_to_test = [0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0]
eps_to_test = [-0.5, -0.3, -0.1, 0.1, 0.3, 0.5]
fermi_energies =  [5.9771, 6.1252, 6.1799, 6.1011, 6.1113, 6.0125, 5.9858, 6.1180, 6.1791, 6.1069, 6.1088, 6.0142, 6.0315, 6.1011, 6.1770, 6.1168, 6.1045, 6.0195, 6.0148, 6.0896, 6.1734, 6.1272, 6.1022, 6.0294, 6.0257, 6.1312, 6.1685, 6.1372, 6.1043, 6.0524, 6.0658, 6.1375, 6.1623, 6.1465, 6.1114, 6.0770, 6.1026, 6.1231, 6.1549, 6.1549, 6.1231, 6.1026]


gamma_to_test = np.array(gamma_to_test)
eps_to_test = np.array(eps_to_test)
fermi_energies = np.array(fermi_energies) 

gamma_to_plot, eps_to_plot = np.meshgrid(gamma_to_test, eps_to_test)

fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))

ax.set_thetamin(0)   # Start angle in degrees
ax.set_thetamax(60)  # End angle in degrees

theta_ticks = np.arange(0, 70, 10)  # Create tick locations in degrees
ax.set_xticks(np.radians(theta_ticks))  # Convert to radians for set_xticks

# (eps,gamma) and (-eps,60-gamma) correspnod to the same shape - covert to positive eps so that it can be plotted in polar coordinates
for r in range(len(gamma_to_test)):
    for c in range(len(eps_to_test)):
        if(eps_to_plot[c][r]<0):
            eps_to_plot[c][r] *= -1
            gamma_to_plot[c][r] = 60-gamma_to_plot[c][r]
        gamma_to_plot[c][r] *= np.pi/180
        plt.polar(gamma_to_plot[c][r], eps_to_plot[c][r], 'wx')
            
fermi_energies_to_plot = np.array(fermi_energies).reshape((len(eps_to_test),len(gamma_to_test)))

cax = ax.contourf(gamma_to_plot, eps_to_plot, fermi_energies_to_plot, 10)

plt.colorbar(cax)
plt.show()
'''
