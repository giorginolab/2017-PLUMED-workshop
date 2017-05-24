# A series of examples of the PLUMED CV object-oriented model in HTMD.

# First, install HTMD according to the instructions. To use the latest
# version, clone from https://github.com/Acellera/htmd and prepend the
# directory to PYTHONPATH. Note that as of June 2017 the features are
# under development in the "toni-devel-plumed" branch, which will be
# merged later.

# In the unlikely case that  plumed 2 is not in path
# import os
# os.environ["PATH"]="/usr/local/bin:/Users/toni/bin:"+os.environ["PATH"]


# Import HTMD
from htmd import *

# Until this is merged in master
from htmd.projections.metricplumed2 import *

# This is just to check if and where PLUMED2 is found.
htmd.projections.metricplumed2._getPlumedRoot()



# Trivial example: colvars as strings
kidkix = Molecule("1kdx")
kidkix.numFrames                    # 17 frames

metric = MetricPlumed2(['d1: DISTANCE ATOMS=1,200',
                        'd2: DISTANCE ATOMS=5,6'])

metric.project(kidkix)



# The "COORDINATION" variable, encapsulated in the PlumedCV object
contacts = PlumedCV("COORDINATION", label="coo",
                    R_0=7,
                    GROUPA=PlumedGroup(kidkix,"kix","chain A and name CA"),
                    GROUPB=PlumedGroup(kidkix,"kid","chain B and name CA") )
metric2 = MetricPlumed2(contacts)
print(metric2)
metric2.project(kidkix)



# ------------------------



# Show help
# 
htmd.projections.metricplumed2.genTemplate("COORDINATION")

htmd.projections.metricplumed2.genTemplate("COORDINATION",include_optional=True)

"""
 This is a very basic tutorial. For further examples, see
  * The presentation of this talk
  * Embedded docstrings and their online version at [1], to be updated
  * The other files, especially villin colvars
"""  


[1] https://www.htmd.org/docs/latest/htmd.projections.metricplumed2.html 
