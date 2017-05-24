# A series of examples of the PLUMED CV object-oriented model in
# HTMD. If you haven't done so yet, please review the simpler examples
# before. Here we try and follow the "standard" Markov workflow for
# the analysis of biomolecular simulations.

# We shall follow the steps used at
# https://www.htmd.org/docs/latest/tutorials/protein-folding-analysis.html
# for this simple Ala3 system, which will be analyzed using the 6
# phi-psi Ramachandran angles computed via Plumed. Note that this is
# for illustration only: HTMD has "native" backbone torsion functions
# for this.

# The analysis is based on a single 1-us long trajectory of Ala3
# solvated with ~8000 water molecules and simulated in the NPT
# ensemble with OpenMM. Integration timestep was 2 fs, writing a frame
# every 0.1 ps.

# The results are not guaranteed to be correct. In particular, the
# clustering function ignores the periodicity.


from htmd import *
from htmd.projections.metricplumed2 import *
htmd.projections.metricplumed2._getPlumedRoot()

# Trivial example - single trajectory. Done this way we can only
# evaluate CVs, but not build a Markov model

# Step 1 - Load the trajectories

ala3=Molecule("ala3.prmtop")
ala3.read("ala3_npt_1us_every_100ps.xtc")

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


# Since we are at it, strip out water
ala3.filter("not water")
ala3.write("ala3_nowater.pdb")
ala3.write("ala3_npt_1us_every_100ps_nowater.xtc")


# ------------------------

# Now use a simlist, i.e. (potentially) multiple trajectories. This is
# more involved, but enables the Model classes etc. We can build a
# simlist even with a single trajectory, as in this case.  Each
# trajectory has to go in a separate directory.

# Step 1 - Load the trajectories

sl=simlist(['ala3/traj1'],"ala3_nowater.pdb")


# Step 2 - Decide the lower-dimensional projection

ala3_rama="""
TORSION ATOMS=5,7,9,15   LABEL=ALA1_PHI
TORSION ATOMS=7,9,15,17  LABEL=ALA1_PSI
TORSION ATOMS=15,17,19,25  LABEL=ALA2_PHI
TORSION ATOMS=17,19,25,27  LABEL=ALA2_PSI
TORSION ATOMS=25,27,29,35  LABEL=ALA3_PHI
TORSION ATOMS=27,29,35,37  LABEL=ALA3_PSI"""
ala3_rama_cv=MetricPlumed2(ala3_rama)

# Step 2a - Create the Metric, bind it, and project
ala3_metric = Metric(sl)

ala3_metric.set(ala3_rama_cv)

ala3_data=ala3_metric.project()
ala3_data.fstep=.1

ala3_data.dat.shape             # 1 element only
ala3_data.dat[0].shape          # 10000 frames x 6 CVs as expected


# Step 3 - Perform the kinetic analysis. We follow the default workflow.

# Optional: TICA
# Optional: bootstrap

# Step 3a - Cluster in microstates
ala3_data.cluster(MiniBatchKMeans(n_clusters=100))

# Step 3b - Build MM
model=Model(ala3_data)

# Plot timescales
lags=numpy.linspace(1,stop=100,num=100)
model.plotTimescales(lags=lags,units="frames")

# Decide lag time is 2 ns, 4 macrostates
model.markovModel(2,4,units="ns")

# Plot fes
model.eqDistribution()

# FES + microstates

# Visualize the states (to choose source and sink states)
model.viewStates(protein=True)


# Kinetics on and off and dG
kin = Kinetics(model, temperature=300,  source=3, sink=0)

kin.getRates()
kin.plotRates()


# Plot flux pathways




