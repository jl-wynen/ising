# Comparison of implementations of the 2D Ising model

This directory contains several implementations of the 2D model in the languages
- C++
- Python
- Rust

The comparison simulations are restricted to h=0 and parameterize the model through the temperature T, setting J=1.

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
