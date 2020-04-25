# PyQUESTAAL


This is a calculator class that has been written to interface calculations done using [QUESTAAL](http://questaal.org) with highthroughput calculators like [ASE](https://wiki.fysik.dtu.dk/ase/). Note: The module can be used as a stad alone python controler for questaal jobs. 

### Features

  - make symmetry line files supported by QUESTAAL on the fly with a given structure
  - plot bands directly *(requires plotquestaal.py)*
  - relaxiations included
  - can read forces from output files
  - Control number of processors



### Installation
 1. Copy the files lmf.py to the working directory
 2. Import them and use !
 3. Make sure "lmf" is recogonizable in the path and also modify the type of mpi call. default "mpirun"
 
### Examples

Example usage shown in 	[examples.ipynb](https://github.com/santoshkumarradha/pyquestaal/blob/pyquestaal/example.ipynb)

### To Doooooooo......
   - control over all variables
   - Interface with pyamtgen
   - Add the developed lattice relaxation module which uses Genetic Algorithm
   - Make more standalone
   - More complex praser for questaal output
