import os
import subprocess

# User setup
atdynpath = '/home/okuwaki/GENESIS_Tutorials-2022/Programs/genesis-2.1.2/bin/atdyn'

# Read directory names from dirnames.in
with open("dirnames.in", "r") as f:
    directories = [line.strip() for line in f if line.strip()]

file_names = ["ligand", "receptor", "complex"]
content_template = """
[INPUT]
prmtopfile = {name}.prmtop   # AMBER parameter topology file
ambcrdfile = {name}.inpcrd      # AMBER coordinates file

[OUTPUT]
dcdfile = {name}.dcd          # DCD trajectory file
rstfile = {name}.rst          # restart file

[ENERGY]
forcefield       = AMBER   # [AMBER]
electrostatic    = CUTOFF  # [CUTOFF]
switchdist       = 99.9    # switch distance
cutoffdist       = 99.9    # cutoff distance
pairlistdist     = 100.0    # pair-list distance
implicit_solvent = GBSA    # [GBSA]
gbsa_salt_cons   = 0.15
gbsa_surf_tens   = 0.0072
gbsa_vdw_offset  = 0.09

[MINIMIZE]
method           = SD      # [SD] or [LBFGS]
nsteps           = 1       # number of minimization steps
eneout_period    = 1       # energy output period
crdout_period    = 1       # coordinates output period
rstout_period    = 1       # restart output period
nbupdate_period  = 1

[BOUNDARY]
type             = NOBC    # [NOBC]
"""

for directory in directories:
    if not os.path.exists(directory):
        os.makedirs(directory)
    for name in file_names:
        # Change to the directory
        os.chdir(directory)

        content = content_template.format(name=name)
        inp_file = f"{name}.inp"
        log_file = f"{name}.log"
        with open(inp_file, "w") as file:
            file.write(content)
        # Execute the command
        subprocess.run(f"mpirun -np 1 {atdynpath} {inp_file} | tee {log_file}", shell=True)
        # Change back to the original directory
        os.chdir("..")

