import os
import subprocess
import matplotlib.pyplot as plt

# Variables
group2_value = "rno:306"  # Change this value as needed
rmsd_analysis_path = "/home/okuwaki/GENESIS_Tutorials-2022/Programs/genesis-2.1.2/bin/rmsd_analysis"
combined_output_file = 'combined_output.txt'

# Template for rmsd.inp
rmsd_inp_template = """
[INPUT]
prmtopfile          = complex.prmtop
ambcrdfile          = complex.inpcrd

[OUTPUT]
rmsfile             = output.rms

[TRAJECTORY]
trjfile1      = step4_out.dcd    # trajectory file
md_step1      = 1000                  # number of MD steps
mdout_period1 = 1000                  # MD output period
ana_period1   = 1000                  # analysis period
repeat1       = 1
trj_format    = DCD                     # (PDB/DCD)
trj_type      = COOR+BOX                # (COOR/COOR+BOX)
trj_natom     = 0                       # (0:uses reference PDB atom count)

[SELECTION]
group1        = an:CA
group2        = {group2_value}

[FITTING]
fitting_method = TR+ROT           # [NO,TR,TR+ROT,TR+ZROT,XYTR,XYTR+ZROT]
fitting_atom   = 1                # atom group
mass_weight    = NO               # mass-weight is not applied
 
[OPTION]
check_only     = NO               # (YES/NO)
analysis_atom  = 2                # atom group
"""

# Read directory names from dirnames.in
with open('dirnames.in', 'r') as file:
    directories = file.read().splitlines()

# Function to create rmsd.inp
def create_rmsd_inp(group2_value):
    return rmsd_inp_template.format(group2_value=group2_value)

# Iterate through directories and run rmsd_analysis
for directory in directories:
    os.chdir(directory)
    rmsd_inp_content = create_rmsd_inp(group2_value)
    with open('rmsd.inp', 'w') as file:
        file.write(rmsd_inp_content)
    command = f"{rmsd_analysis_path} rmsd.inp > rmsd.log"
    subprocess.run(command, shell=True)
    os.chdir('..')

# Collect output.rms files and combine them
with open(combined_output_file, 'w') as outfile:
    for directory in directories:
        rms_file_path = os.path.join(directory, 'output.rms')
        if os.path.exists(rms_file_path):
            with open(rms_file_path, 'r') as infile:
                for line in infile:
                    if not line.startswith('#'):
                        columns = line.split()
                        if len(columns) > 1:
                            outfile.write(f"{directory} {columns[1]}\n")

# Plot the combined results
directories = []
values = []

with open(combined_output_file, 'r') as file:
    for line in file:
        columns = line.split()
        directories.append(columns[0])
        values.append(float(columns[1]))

plt.figure(figsize=(10, 5))
plt.bar(directories, values)
plt.xlabel('Directory')
plt.ylabel('RMSD Value')
plt.title('RMSD Analysis Results')
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig('rmsd_analysis_results.png')
plt.show()


