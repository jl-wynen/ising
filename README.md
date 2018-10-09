# Ising
This repository is a collection of a number of different implementations of the 2D Ising Model in various languages.
All implementations are simple and cover only the basics of the model.
This is meant to serve as a comparison of different programming languages with respect to ease of use and run efficiency.

## Model
The Ising Model is given by the Hamiltonian (assuming no external magnetic field)  
<img src="https://latex.codecogs.com/svg.latex?H(s)&space;=&space;-J&space;\sum_{\langle&space;i,j&space;\rangle}\,&space;s_i&space;s_j" title="H(s) = -J \sum_{\langle i,j \rangle}\, s_i s_j" />  
where s is a configuration of spins with s_i=+1,-1. The angle brackets denote nearest neighbours.

The programs measure the average magnetisation per spin m = mean(s) and the magnetic susceptibility χ = var(m)/(k_B T).
k_B is the Boltzmann constant and T the temperature.

The only dimensionless parameter of this model is J/(k_B T) = Jβ.
Hence we can set J=k_B=1 without loss of generality.
All calculations are therefore parmeterized through the temperature T.

## Implementation
The simulation is implemented in several different languages:
- C++
- Python
- Rust

Each implementation performs the same basic calculations.
They thermalise for a given number of steps and an initial temperature.
Then they iterate over a range of temperatures and re-thermalise for each one before running production MC evolution and
measuring the energy and average magnetisation.
For each measurement, a whole sweep of the lattice is performed, i.e. Nx*Ny single site updates.

For simplicity, all parameters are hard-coded into the individual programs.
Refer to the README.md files in each subdirectory for language specific information.

## File Format
All implementations write their output into files of the same format, a directory containing two kinds of files:
- `temperatures.dat`: Lists all temperatures of the run with associated indices.
                      Each line of that file has the form `index: temperature`.
- `<idx>.dat`: Observables measured for the temperature associated with index idx.
               First row is the energy, second row the average magnetisation.

The files can be read and analysed via the code provided in directory `ana`.
