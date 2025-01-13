import os
import subprocess

# User setup
spdynpath = '/home/okuwaki/GENESIS_Tutorials-2022/Programs/genesis-2.1.2/bin/spdyn'
np = 1
OMP_NUM_THREADS = 4

minstep = 500 # 5000
heatstep = 500 # 50000
nvtstep = 500 # 100000
nptstep = 1000 # 1000000

# Read directory names from dirnames.in
with open("dirnames.in", "r") as f:
    directories = [line.strip() for line in f if line.strip()]

# Templates for each step
content_template = """
[INPUT]
prmtopfile          = complex.prmtop
ambcrdfile          = complex.inpcrd
ambreffile          = complex.inpcrd

[OUTPUT]
dcdfile             = step1_out.dcd
rstfile             = step1_out.rst

[ENERGY]
electrostatic       = PME
dielec_const        = 1.0
switchdist          = 8.0
cutoffdist          = 8.0
pairlistdist        = 10.0
dispersion_corr     = EPRESS
pme_alpha           = auto
pme_alpha_tol       = 1.0E-5
pme_nspline         = 4
pme_max_spacing     = 1.2
forcefield          = AMBER
nonb_limiter        = YES
minimum_contact     = 0.5

[MINIMIZE]
method              = SD
nsteps              = {minstep}
crdout_period       = 500
eneout_period       = 100
rstout_period       = {minstep}
nbupdate_period     = 10
force_scale_init    = 5.0E-5
force_scale_max     = 1.0E-4

[SELECTION]
group1              = an:CA or an:N or an:C or an:O

[RESTRAINTS]
nfunctions          = 1
function1           = POSI
select_index1       = 1
constant1           = 10.0
direction1          = ALL

[BOUNDARY]
type                = PBC
box_size_x          = {box_size_x}
box_size_y          = {box_size_y}
box_size_z          = {box_size_z}
"""

content_template2 = """
[INPUT]
prmtopfile          = complex.prmtop
ambcrdfile          = complex.inpcrd
ambreffile          = complex.inpcrd
rstfile             = step1_out.rst

[OUTPUT]
dcdfile             = step2_out.dcd
dcdvelfile          = step2_out_vel.dcd
rstfile             = step2_out.rst

[DYNAMICS]
integrator          = VVER
timestep            = 0.002
nsteps              = {heatstep}
crdout_period       = 500
iseed               = 31415
annealing           = YES
anneal_period       = 500
dtemperature        = 3.0
initial_time        = 0.0
eneout_period       = 500
velout_period       = 0
rstout_period       = {heatstep}
stoptr_period       = 10
nbupdate_period     = 10
verbose             = YES

[ENSEMBLE]
ensemble            = NVT
temperature         = 0.1
tpcontrol           = BUSSI
tau_t               = 1.0

[ENERGY]
electrostatic       = PME
dielec_const        = 1.0
switchdist          = 8.0
cutoffdist          = 8.0
pairlistdist        = 10.0
dispersion_corr     = EPRESS
pme_alpha           = auto
pme_alpha_tol       = 1.0E-5
pme_nspline         = 4
pme_max_spacing     = 1.2
nonb_limiter        = YES
minimum_contact     = 0.5
forcefield          = AMBER

[CONSTRAINTS]
rigid_bond          = YES
fast_water          = YES
shake_iteration     = 500
shake_tolerance     = 1.0E-10
water_model         = WAT
hydrogen_type       = NAME

[SELECTION]
group1              = an:CA or an:N or an:C or an:O

[RESTRAINTS]
nfunctions          = 1
function1           = POSI
select_index1       = 1
constant1           = 10.0
direction1          = ALL

[BOUNDARY]
type                = PBC
"""

content_template3 = """
[INPUT]
prmtopfile          = complex.prmtop
ambcrdfile          = complex.inpcrd
ambreffile          = complex.inpcrd
rstfile             = step2_out.rst

[OUTPUT]
dcdfile             = step3_out.dcd
dcdvelfile          = step3_out_vel.dcd
rstfile             = step3_out.rst

[DYNAMICS]
integrator          = VVER
timestep            = 0.002
nsteps              = {nvtstep}
crdout_period       = 500
iseed               = 31415
annealing           = NO
initial_time        = 0.0
eneout_period       = 500
velout_period       = 0
rstout_period       = {nvtstep}
stoptr_period       = 10
nbupdate_period     = 10
verbose             = YES

[ENSEMBLE]
ensemble            = NVT
temperature         = 300.0
tpcontrol           = BUSSI
tau_t               = 5.0

[ENERGY]
electrostatic       = PME
dielec_const        = 1.0
switchdist          = 8.0
cutoffdist          = 8.0
pairlistdist        = 10.0
dispersion_corr     = EPRESS
pme_alpha           = auto
pme_alpha_tol       = 1.0E-5
pme_nspline         = 4
pme_max_spacing     = 1.2
nonb_limiter        = YES
minimum_contact     = 0.5
forcefield          = AMBER

[CONSTRAINTS]
rigid_bond          = YES
fast_water          = YES
shake_iteration     = 500
shake_tolerance     = 1.0E-10
water_model         = WAT
hydrogen_type       = NAME

[SELECTION]
group1              = an:CA or an:N or an:C or an:O

[RESTRAINTS]
nfunctions          = 1
function1           = POSI
select_index1       = 1
constant1           = 10.0
direction1          = ALL

[BOUNDARY]
type                = PBC
"""

content_template4 = """
[INPUT]
prmtopfile          = complex.prmtop
ambcrdfile          = complex.inpcrd
ambreffile          = complex.inpcrd
rstfile             = step3_out.rst


[OUTPUT]
dcdfile             = step4_out.dcd
dcdvelfile          = step4_out_vel.dcd
rstfile             = step4_out.rst

[DYNAMICS]
integrator          = VVER
timestep            = 0.002
nsteps              = {nptstep}
crdout_period       = 1000
iseed               = 31415
annealing           = NO
initial_time        = 0.0
eneout_period       = 1000
velout_period       = 0
rstout_period       = {nptstep}
stoptr_period       = 10
nbupdate_period     = 10
verbose             = YES

[ENSEMBLE]
ensemble            = NPT
temperature         = 300.0
tpcontrol           = BUSSI
pressure            = 1.0
tau_p               = 5.0
tau_t               = 5.0
isotropy            = ISO

[ENERGY]
electrostatic       = PME
dielec_const        = 1.0
switchdist          = 8.0
cutoffdist          = 8.0
pairlistdist        = 10.0
dispersion_corr     = EPRESS
pme_alpha           = auto
pme_alpha_tol       = 1.0E-5
pme_nspline         = 4
pme_max_spacing     = 1.2
nonb_limiter        = YES
minimum_contact     = 0.5
forcefield          = AMBER

[CONSTRAINTS]
rigid_bond          = YES
fast_water          = YES
shake_iteration     = 500
shake_tolerance     = 1.0E-10
water_model         = WAT
hydrogen_type       = NAME

[SELECTION]

[BOUNDARY]
type                = PBC
"""

os.environ["OMP_NUM_THREADS"] = str(OMP_NUM_THREADS)

templates = [
    ("step1.inp", content_template),
    ("step2.inp", content_template2),
    ("step3.inp", content_template3),
    ("step4.inp", content_template4)
]

# Function to write template content to file and run MD calculation
def run_md_step(step_file, content):
    inp_file = step_file
    log_file = f"{step_file}.log"

    # Write the content to the .inp file
    with open(inp_file, "w") as f:
        f.write(content)

    # Run the MD calculation
    subprocess.run(f"mpirun -np {np} {spdynpath} {inp_file} | tee {log_file}", shell=True)

# Function to get box sizes from complex.inpcrd
def get_box_sizes():
    with open("complex.inpcrd", "r") as f:
        lines = f.readlines()
        last_line = lines[-1].split()
        box_size_x = last_line[0]
        box_size_y = last_line[1]
        box_size_z = last_line[2]
    return box_size_x, box_size_y, box_size_z

# Iterate over each directory and run the MD steps
for directory in directories:
    os.chdir(directory)
    box_size_x, box_size_y, box_size_z = get_box_sizes()

    # Replace placeholders in templates with actual box sizes
    for step_file, content in templates:
        content = content.format(box_size_x=box_size_x, box_size_y=box_size_y, box_size_z=box_size_z,
                                 minstep=minstep, heatstep=heatstep, nvtstep=nvtstep, nptstep=nptstep)
        run_md_step(step_file, content)

    os.chdir("..")


