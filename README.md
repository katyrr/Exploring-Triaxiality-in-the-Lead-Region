## INTRODUCTION

This repository contains a Python interface to the set of three Particle-plus-Triaxial-Rotor programmes, which are used to model the spectroscopic properties of an odd-mass deformed nucleus.

These were originally described in: S.E. Larsson, G. Leander, and I. Ragnarsson. “Nuclear core-quasiparticle coupling”. In: Nuclear Physics A 307.2 (Sept. 1978), pp. 189-223. issn: 0375-9474. DOI: 10.1016/0375-9474(78)90613-9. URL: http://dx.doi.org/10.1016/0375-9474(78)90613-9,

and updated in: Ingemar Ragnarsson and Paul B. Semmes. “Description of nuclear moments and nuclear spectra in the particle-rotor model”. In: Hyperfine Interactions 43.1–4 (Dec. 1988), pp. 423–440. issn: 1572-9540. DOI: 10 . 1007 / bf02398323. URL: http://dx.doi.org/10.1007/BF02398323.

The code in this repository automates calculations over a range of deformations (quadrupole ε and/or triaxial γ). It is designed to run in parallel on an 8-core machine. This can be adapted to a machine of any size by editing the variable "num_cores" in the config file.

Developed for an MSci dissertation project in 2024-25.
Some changes have been made to improve usability since the completion of the project.



## CONTENTS OF THIS REPOSITORY

The `/Code` root folder contains the main.py file. The code should always be run from the root of the repo.

The `/Code/data` folder contains example config files for three nuclei (Au179, Pt177, Pb207) and an example folder showing how input files, output files, binaries, figures, and helper scripts are organised when the code is run.

The `/Code/src` folder contains definitions of classes and functions used in main.py. It also has a ptrm folder which contains the original pre-compiled PTRM codes for both MacOS and Windows (64bit). 
(NOTE: The `optuna` folder contains experimental code which implemented a machine learning style method for optimising input parameters using python's Optuna framework, but this has not been maintained.)

The `/Code/static` folder contains a config.txt template.

The `/Code/tests` folder contains pytest unit tests, which can be run with "uv run pytest".



## HOW TO USE FOR THE FIRST TIME: 

1. Download the repository from GitHub: 
   ```
   git clone https://github.com/katyrr/Exploring-Triaxiality-in-the-Lead-Region
   cd Exploring-Triaxiality-in-the-Lead-Region
   ```  

   Optional:  
   If you have uv installed, you can then run:
   ```
   uv venv
   uv sync
   ```

   There is no need to "activate" the venv before running. 
   Using uv is optional, but might reduce the risk of version errors.

2. Check that your current working directory is `/Code`. 
   If you run `ls` you should see `main.py` in the list of files.

3. Do a test run of the code:
   `uv run main.py example_folder`

   or (without uv):
   `python main.py example_folder`
    
   or (depending on your OS and python version): 
   `Python3 main.py example_folder`

   1. It will automatically generate a folder called `example_folder` in the 
      directory `/Code/data` where all the input files, output files, binaries, figures, 
      and helper scripts will be stored. You can give this folder any name (typically it
      might be the name of the nucleus you're studying, e.g. Pt177). 

   2. It will automatically generate a config.txt file in `example_folder`, from 
      the template at `/Code/static/config_template.txt`.  

      You may wish to edit the template beforehand, to configure settings such as your
      operating system, how many CPU cores you want to utilise, your desired figure 
      resolution and formatting, which graphs you want plotted by default, and which nucleus 
      to calculate for.  
   
      The default is 177Pt at (ε, γ) = (0.26, 24º), on MacOS with 8 cores.  

      The template will be used whenever you run the code in a new data folder for the 
      first time, so you may find it convenient to adjust your default settings now!

   3. It will attempt to enable 'execute' permissions on the ptrm codes, which are stored
      in `/Code/src/ptrm/MacOS/MO` or `/Code/src/ptrm/64bit/MO`. Depending on your operating 
      system, it may or may not be successful. If it doesn't work, it will throw an error, 
      and you will have to enable these permissions manually:  

      On Mac: navigate to the `/Code/src/ptrm/MacOS/MO` folder in Terminal, and then 
      use the commands `chmod 775 gampn`, `chmod 775 asyrmo`, `chmod 775 probamo`. 
      Then use `ls -ltr` to see that the permissions (in the leftmost column) now 
      have `x`s (execute). You might also have to grant permission in Settings 
      (after attempting to run for the first time, the permission request will appear in Privacy).  

      On Windows: navigate to the `/Code/src/ptrm/64bit/MO` folder in PowerShell, and then 
      use the commands `Unblock-File -Path GAMPN.exe`, `Unblock-File -Path ASYRMO.exe`, 
      `Unblock-File -Path PROBAMO.exe`.  

   4. It will print progress and key results to the command line, and you can find more
      detailed results in the `example_folder/outputs`, and figures in 
      `example_folder/figures`.  
       
      Optional:
      You can reroute the stdout to a text file with:
      ```
      touch path/to/saved_stdout.txt                             # creates the file
      uv run main.py example_folder > path/to/saved_stdout.txt   # sends stdout to the file
      ```



## HOW TO USE AFTER THE FIRST TIME: 

1. 1. If you just want to make some small adjustments to the config, you can edit the (existing)
      config file in `/Code/data/example_folder`, and move to step 2.  
    
   2. If you want to start a new folder (e.g. to calculate for a new nucleus), you can either
      manually create a new folder in `/Code/data`, copy a config file from the template (or elsewhere), and make the desired changes before moving to step 2.  
       
      Or you can skip straight to step 2 and let the code automatically generate a new folder 
      and config from the template.  
    
      (i.e. if `example_folder` already exists with a config file, it will use that. 
      If it doesn't already exist, it will make one.)

2. Run the codes with command: `uv run python src/main.py example_folder`
                            or `python src/main.py example_folder`
                            or `Python3 src/main.py example_folder`
    
    Optional command line arguments:

    `--display-figures` or `-d`:
        Display figures on the screen as they are plotted. 
        (whether or not this option is used, the figures will be saved in the "figures" folder)

    `--verbose` or `-v`:
        Print more detailed information in stdout.

    `--replot` or `-rp`:
        Don't recalculate results, just replot graphs using most recent (existing) output files.
        
        This can be useful in large data sets, when you want some graphs that you didn't request to be plotted the first time (in the config file), or perhaps to change the formatting of the graphs.

        This option can't be used if no output files are found.
     
