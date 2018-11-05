# Python Implementation of the Ising Model Simulation

## Requirements
- Python 3
- NumPy
- Numba

## Usage
Adjust run parameters (number of sweeps, lattice size, etc.) in the marked section in the
beginning of ising.py.
Then run
```
python ising.py [datadir]
```
where `datadir` is an optional argument to specify a directory to write the output files to.
It defaults to `data`.
Note that the directory and all files inside it get deleted if it already exists!
