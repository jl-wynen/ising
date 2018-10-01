# C++ Implementation of the Ising Model Simulation

## Requirements
- CMake version >= 3.9.0
- C++17 compiler

## Usage
- Prepare build:
```
mkdir build
cd build
cmake .. [-DCMAKE_BUILD_TYPE=RELEASE]
```

- Adjust run parameters (number of sweeps, lattice size, etc.) in the marked section in the
  beginning of ising.cpp.

- Compile:
```
cmake --build .
```

- Run:
```
ising [datadir]
```
where `datadir` is an optional argument to specify a directory to write the output files to.
It defaults to `data`.
Note that the directory and all files inside it get deleted if it already exists!
