# In development, updated version with examples coming soon

This is a packing defect script for membranes or lipid droplets. In the 'scripts' folder, we have 2 scripts: run_defect.py calculates the defects of the system, incuding those below the protein. On the other hand, run_radius.py modifies the output to then only take into account defects ignoring a given radius around the protein.

* The run_defect.py creates a .gro file for each frame in the trajectory.
* run_radius.py takes those gro files, and deletes defects within a given radius of the protein.
* GROMACS has to in python path becuase we use gmx to renumber

There is an included 10-step trajectoy and topology file in the "traj" folder

python scripts/run_defect.py
python scripts/run_radius.py

will make output file of "filtered_gro_paper" (just example name)


