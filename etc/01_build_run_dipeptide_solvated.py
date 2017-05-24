

# In[3]:

from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout


# In[5]:

tleap_in=b"""
source leaprc.protein.ff14SB
source leaprc.water.tip3p
m = sequence { ACE ALA NME }
solvatebox m TIP3PBOX 8
saveAmberParm m diala_s.prmtop diala_s.rst7
savepdb m diala_s.pdb
quit
"""

from subprocess import Popen, PIPE, STDOUT
tleap_out,tleap_err = Popen(['tleap','-f','-'], stdout=PIPE, stdin=PIPE, stderr=PIPE)    .communicate(input=tleap_in)

if tleap_err != b'':
    print("****** Something went wrong, STDERR follows ******")
    print(tleap_err.decode())
    print("****** STDOUT follows ******")

print(tleap_out.decode())


# In[6]:

prmtop = AmberPrmtopFile('diala_s.prmtop')
inpcrd = AmberInpcrdFile('diala_s.rst7')


# In[8]:

system = prmtop.createSystem(nonbondedMethod=CutoffPeriodic, 
                             nonbondedCutoff=1*nanometer,
                             constraints=HBonds)


# In[9]:

barostat = openmm.MonteCarloBarostat(1.0*bar, 300.0*kelvin, 25)
system.addForce(barostat)


# In[10]:

integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)


# In[14]:

try:
    platform = openmm.Platform.getPlatformByName('CUDA')
except:
    platform = openmm.Platform.getPlatformByName('CPU')       

simulation = Simulation(prmtop.topology, system, integrator, platform)
simulation.context.setPositions(inpcrd.positions)

if inpcrd.boxVectors is not None:
    simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)

PDBReporter("pre_s.pdb",1).report(simulation,simulation.context.getState(getPositions=True))


# In[17]:

simulation.minimizeEnergy()
PDBReporter("mini_s.pdb",1).report(simulation,simulation.context.getState(getPositions=True))


# In[19]:

simulation.reporters.append(DCDReporter('output_s.dcd', 1000))
simulation.reporters.append(StateDataReporter("output_s.log", 1000, step=True,
        potentialEnergy=True, temperature=True,volume=True))
simulation.step(500000000)
simulation.saveState("checkpoint_s.state")


# -------------

# In[ ]:



