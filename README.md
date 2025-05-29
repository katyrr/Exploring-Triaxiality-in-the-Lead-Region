This repository contains a Python interface to the set of three Particle-plus-Triaxial-Rotor programmes, which are used to model the spectroscopic properties of an odd-mass deformed nucleus.

These were originally described in: S.E. Larsson, G. Leander, and I. Ragnarsson. “Nuclear core-quasiparticle coupling”. In: Nuclear Physics A 307.2 (Sept. 1978), pp. 189-223. issn: 0375-9474. DOI: 10.1016/0375-9474(78)90613-9. URL: http://dx.doi.org/10.1016/0375-9474(78)90613-9,

and updated in: Ingemar Ragnarsson and Paul B. Semmes. “Description of nuclear moments and nuclear spectra in the particle-rotor model”. In: Hyperfine Interactions 43.1–4 (Dec. 1988), pp. 423–440. issn: 1572-9540. DOI: 10 . 1007 / bf02398323. URL: http://dx.doi.org/10.1007/BF02398323.

The code in this repository automates calculations over a range of deformations (quadrupole ε and triaxial γ). It is designed to run in parallel on an 8-core machine. This can be adapted to a machine of any size by editing the variable batch_settings["num_batches"] in section 5 of main.py.

Contains a main.py file, plus three module files "functions.py", "structs.py", and "graph_plotting.py", in the folder "/Executables".

The "/Examples" folder shows the file structure, and provides some example input and output files. The three nuclei folders contain example config files for input to the code.

Developed for an MSci dissertation project in 2024-25.
