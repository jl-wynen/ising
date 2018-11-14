# N-Dimensional Ising Model

This directory contains a C++ implementation and Python analysis scripts for the N-dimensional Ising model.
It differs from the implementation in [comparison/cpp](/comparison/cpp) in many ways.
It allows for non zero external field h and does not use the temperature to parameterize the model but rahter the dimensionless J and h.
The implementation is much more structured and spread out over several files.

## Build
```
mkdir build
cd build
cmake .. [-DCMAKE_BUILD_TYPE=RELEASE]
cmake --build .
```

## Tests
You can build and run unit tests from the build directory via
```
make ising-test
ising-test
```

## Run
The program takes two arguments:
```
ising <infile> <outdir>
```
where
- `infile` is a YAML file describing the run, see the sample input file [input.yml](n-dimensional/input.yml).
- `outdir` is a directory to write the output files to.

## Analysis
Plotting and analysis scripts can be found in [ana](n-dimensional/ana).
