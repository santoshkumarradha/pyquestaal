# PyQUESTAAL

[![DOI](https://zenodo.org/badge/258655688.svg)](https://zenodo.org/badge/latestdoi/258655688)

This is a calculator class that has been written to interface calculations done using [QUESTAAL](http://questaal.org) with highthroughput calculators like [ASE](https://wiki.fysik.dtu.dk/ase/). 

>Note: The module can be used as a stand alone python controller for questaal jobs as well.

### Features

  - make symmetry line files supported by QUESTAAL on the fly with a given structure
  - plot bands directly *(requires plotquestaal.py)*
  - relaxations included
  - can read forces from output files
  - Control number of processors



### Installation
 1. Copy the files lmf.py to the working directory
 2. Import them and use !
 3. Make sure "lmf" is recognizable in the path and also modify the type of mpi call. default "mpirun"
 
### Examples
Simple example usage shown in 	[simple_examples.ipynb](https://github.com/santoshkumarradha/pyquestaal/blob/pyquestaal/simple_example.ipynb) 
Comprehensive example shown in [ tutorial_example.ipynb](https://github.com/santoshkumarradha/pyquestaal/blob/pyquestaal/tutorial_example.ipynb)

To use it with ASE, one can do something like
```python
from lmf import lmf #load the lmf calculator
import numpy as np
from ase.io import read
from ase.units import kJ
from ase.eos import EquationOfState
from ase import Atoms
from ase.io.trajectory import Trajectory
def test():
    a = 4.0  # approximate lattice constant
    b = a / 2
    ag = Atoms('Ag',
               cell=[(0, b, b), (b, 0, b), (b, b, 0)],
               pbc=1,
               calculator=lmf()) # Use QUESTAAL's LMF as calculator
    cell = ag.get_cell()
    traj = Trajectory('Ag.traj', 'w')
    for x in np.linspace(0.95, 1.05, 5):
        ag.set_cell(cell * x, scale_atoms=True)
        ag.get_potential_energy()
        traj.write(ag)
    

    configs = read('Ag.traj@0:5')  # read 5 configurations
    # Extract volumes and energies:
    volumes = [ag.get_volume() for ag in configs]
    energies = [ag.get_potential_energy() for ag in configs]
    eos = EquationOfState(volumes, energies)
    v0, e0, B = eos.fit()
    print(B / kJ * 1.0e24, 'GPa')
    eos.plot('Ag-eos.png')
test()
```
### Citation

If you find this work useful, please cite using
```
@misc{https://doi.org/10.5281/zenodo.4292415,
  doi = {10.5281/ZENODO.4292415},
  url = {https://zenodo.org/record/4292415},
  author = {Radha, Santosh Kumar},
  title = {pyQuestaal: An python interface or Questaal quantum codes.},
  publisher = {Zenodo},
  year = {2020},
  copyright = {Open Access}
}
```

### To Do
   - control over all variables
   - Interface with pyamtgen
   - Add the developed lattice relaxation module which uses Genetic Algorithm
   - Make more standalone
   - More complex parser for questaal output (look for questaal-reader repo)


