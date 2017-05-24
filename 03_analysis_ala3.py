# Following https://www.htmd.org/docs/latest/tutorials/protein-folding-analysis.html


from htmd import *
from htmd.projections.metricplumed2 import *
htmd.projections.metricplumed2._getPlumedRoot()

# Trivial example

# Step 1 - Load the trajectories

ala3=Molecule("ala3.prmtop")
ala3.read("ala3_nvt_100ps.xtc")


# Step 2 - Decide the lower-dimensional projection

ala3_rama="""
TORSION ATOMS=5,7,9,15   LABEL=ALA1_PHI
TORSION ATOMS=7,9,15,17  LABEL=ALA1_PSI
TORSION ATOMS=15,17,19,25  LABEL=ALA2_PHI
TORSION ATOMS=17,19,25,27  LABEL=ALA2_PSI
TORSION ATOMS=25,27,29,35  LABEL=ALA3_PHI
TORSION ATOMS=27,29,35,37  LABEL=ALA3_PSI"""
ala3_rama_cv=MetricPlumed2(ala3_rama)
ala3_rama_data=ala3_rama_cv.project(ala3)

ala3_rama_data.shape
# (nframes,ncv)



# Step 3 - Perform the kinetic analysis. We follow the default workflow.

# Step 3a - Cluster in microstates

# Step 3b - Build MM

# Plot timescales

# Decide lag time

# Plot fes

# FES + microstates

# Visualize the states (to choose source and sink states)

# Kinetics on and off and dG

# Plot flux pathways




