# Following https://www.htmd.org/docs/latest/tutorials/protein-folding-analysis.html

# Villin dataset (partitioned into 3 parts), 3GB:
# http://pub.htmd.org/tutorials/protein-folding-analysis/datasets.tar.gz

from htmd import *
from htmd.projections.metricplumed2 import *
htmd.projections.metricplumed2._getPlumedRoot()


# Trivial example
kidkix = Molecule("1kdx")
kidkix.numFrames                    # 17 frames

metric = MetricPlumed2(['d1: DISTANCE ATOMS=1,200',
                        'd2: DISTANCE ATOMS=5,6'])

metric.project(kidkix)


contacts = PlumedCV("COORDINATION", label="coo",
                    R_0=7,
                    GROUPA=PlumedGroup(kidkix,"kix","chain A and name CA"),
                    GROUPB=PlumedGroup(kidkix,"kid","chain B and name CA") )
metric2 = MetricPlumed2(contacts)
print(metric2)
metric2.project(kidkix)




# Show help
# 
htmd.projections.metricplumed2.genTemplate("COORDINATION")
htmd.projections.metricplumed2.genTemplate("COORDINATION",include_optional=True)



# Frames are saved every 20ps


ala3=Molecule("ala3.prmtop")
ala3.read("ala3_nvt_100ps.xtc")

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



m=Molecule("villin/datasets/1/filtered/filtered.pdb")

mca=m.copy()
mca.filter("name CA")

hydrophobic_resnames = ["ALA","VAL","ILE","LEU","MET","PHE","TYR","TRP","NLE"]
hyd=numpy.in1d(mca.resname,hydrophobic_resnames)
hyd_resid=mca.resid[hyd]

pg=[]
for resid in hyd_resid:
    pg.append(PlumedGroup(m,
                          "hg_{}".format(len(pg)),
                          "resid {} and not name N CA C O and noh".format(resid)))

ngroups=len(pg)

group_pairs=[]
for g1 in range(ngroups):
    for g2 in range(g1+1,ngroups):
        group_pairs.append(PlumedCV("COORDINATION",
                                    GROUPA=pg[g1],
                                    GROUPB=pg[g2],
                                    R_0="3.5",
                                    label="c_{}".format(len(group_pairs)+1)))



armsd=PlumedCV("ALPHARMSD",RESIDUES="all", R_0=1.0, label="armsd")
    
rg=PlumedCV("GYRATION",ATOMS=mca,label="rg")


    
# Step 1 - Load the trajectories
sets = glob('villin/datasets/*/')
sims = []
for s in sets:
    fsims = simlist(glob(s + '/filtered/*/'), 'datasets/1/filtered/filtered.pdb')
    sims = simmerge(sims, fsims)



# Step 2 - Decide the lower-dimensional projection


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




