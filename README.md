## INTRODUCTION

This repository contains a Python interface to the set of three Particle-plus-Triaxial-Rotor programmes, which are used to model the spectroscopic properties of an odd-mass deformed nucleus.

These were originally described in: S.E. Larsson, G. Leander, and I. Ragnarsson. “Nuclear core-quasiparticle coupling”. In: Nuclear Physics A 307.2 (Sept. 1978), pp. 189-223. issn: 0375-9474. DOI: 10.1016/0375-9474(78)90613-9. URL: http://dx.doi.org/10.1016/0375-9474(78)90613-9,

and updated in: Ingemar Ragnarsson and Paul B. Semmes. “Description of nuclear moments and nuclear spectra in the particle-rotor model”. In: Hyperfine Interactions 43.1–4 (Dec. 1988), pp. 423–440. issn: 1572-9540. DOI: 10 . 1007 / bf02398323. URL: http://dx.doi.org/10.1007/BF02398323.

The code in this repository automates calculations over a range of deformations (quadrupole ε and/or triaxial γ). It is designed to run in parallel on an 8-core machine. This can be adapted to a machine of any size by editing the variable "num_cores" in the config file.

Developed for an MSci dissertation project in 2024-25.



## CONTENTS OF THIS REPOSITORY

The /Code/Execuables folder contains a main.py file, plus three module files "functions.py", "structs.py", and "graph_plotting.py".

It also contains the original PTRM pre-compiled codes for MacOS and Windows (64bit).

(Additionally, the files labelled "optuna" can be used as an alternative machine learning style method for optimising input parameters, but this implementation is currently still very basic!)

The Pt177, Au179, Pb207 folders contain example config files for those three nuclei.

The Code/Examples folder contains example inputs, outputs, binary files, shell scripts, etc for Pt177. This folder indicates the required directory structure, but all of the subfolders will be created automatically by the python script when you run it for the first time.




## HOW TO USE 


**For the first time:**

1. Download the repository from GitHub. 
2. Check that the pre-compiled original PTRM codes in the /Code/Executables/MacOS/MO/ or /Code/Executables/64bit/MO/ folder have execute permissions turned on. See the troubleshooting section at the bottom for how to do this.
3. Make a new folder in /Code.
4. Copy the config file from /Code/Examples to your new folder.
5. Make any necessary changes to the settings in the config file.
6. Open main.py in a text editor. At the top of section "1. SET UP", change the 'folder' variable to the name of your new folder. This is so that it can find your config file.
7. Navigate to the /Code/Executables folder in you computer terminal. Run the codes with command "python main.py".


**After the first time:**

1. Make changes to your config file.
2. Check that the 'folder' variable at the top of section "1. SET UP" in main.py is correctly set to the name of the folder you want to run.
3. Navigate to the /Code/Executables folder in you computer terminal. Run the codes with command "python main.py".


**If you don't need to run the whole thing:**

If you open it in Spyder (or similar) then you can run it cell by cell. Especially useful when calculating large data sets, since you can fiddle with graph plotting easily and quickly *after* doing the calculations.




## TROUBLESHOOTING

**Exec format error**

Trying to run the 64bit codes in a Linux operating system (or via ssh)? It won't work; you have to use either Windows or MacOS.

**Permission denied**

The pre-compiled original PTRM codes in the /Code/Executables folder don't have execute permissions. On Mac you can fix this by navigating to the /Code/Executables/MacOS/MO folder in terminal, and then using the commands "chmod 775 gampn", "chmod 775 asyrmo", chmod 775 probamo". Then use "ls -ltr" to see that the permissions (in the leftmost column) now have 'x's (execute). 

On windows this doesn't seem to work. You can get around it by rebooting in linux (or using ssh), modifying the permissions with the same method as above, then returning to Windows to actually run the codes.


The "/Examples" folder shows the file structure, and provides some example input and output files. The three nuclei folders contain example config files for input to the code.

