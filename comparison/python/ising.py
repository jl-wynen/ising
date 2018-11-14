#! /usr/bin/env python3

"""
Python implementation of the Ising model simulation.
"""

import sys
import dataclasses
from pathlib import Path
import math
import time

import numpy as np
import numba

#--------------------------
# Set run parameters here.

NTHERM_INIT = 1000  # number of thermalisation sweeps in the beginning
NTHERM = 1000  # number of thermalisation sweeps per temperature
NPROD = 10000  # number of production sweeps (with measurements) per temperature

NX = 5  # number of lattice sites in x direction
NY = 5  # number of lattice sites in y direction

# all the temperatures we want to measure at
# TEMPERATURES = np.concatenate((np.arange(6., 0.4, -0.5), np.arange(0.4, 6., 0.5)))
TEMPERATURES = np.arange(6., 0.4, -0.5)

# seed for the random number generator
SEED = 538

# End of run parameters.
#------------------------

def prepare_datadir():
    """
    Create the output data directory and write the temperature file.
    Deletes the directory and all its contents if it exists.

    Returns:
        Path to the directory.
    """

    if len(sys.argv) == 2:
        datadir = Path(__file__).resolve().parent/sys.argv[1]
    else:
        datadir = Path(__file__).resolve().parent/"data"
    print(f"Saving data in '{datadir}'")

    # delete if exists
    if datadir.exists():
        print(f"Data directory '{datadir}' exists, deleting!")
        for fil in datadir.iterdir():
            fil.unlink()
        datadir.rmdir()

    # create directory
    datadir.mkdir()

    # write temperatures
    with open(datadir/"temperatures.dat", "w") as tempf:
        for i, temp in enumerate(TEMPERATURES):
            tempf.write(f"{i:d}: {temp:f}\n")

    return datadir

@dataclasses.dataclass
class Observables:
    """
    Dataclass to store observables for each configuration on
    which measurements are performed.
    Each member is a list containing the Monte-Carlo history of that observable.
    """
    energy: list
    magnetisation: list

def hamiltonian(cfg):
    """
    Evaluate the Hamiltonian on a given configuration.
    """
    # Use numpy.roll to shift the array in x or y direction by one site.
    # This way, cfg*no.roll(...) multiplies nearest neighbours with the
    # speed of numpy.
    return - 2*(np.sum(cfg*np.roll(cfg, -1, axis=0))
                + np.sum(cfg*np.roll(cfg, -1, axis=1)))

def magnetisation(cfg):
    """
    Compute the magnetisation of a given configuration.
    """
    return np.mean(cfg)

# This function is simple but gets called a lot.
# So compile it to native code using Numba.
@numba.jit(nopython=True)
def delta_e(cfg, x, y):
    """
    Compute the change in energy if the spin at site (x, y) were flipped.
    """

    xp1 = (x+1)%NX
    xm1 = x-1 if x > 0 else NX-1
    yp1 = (y+1)%NY
    ym1 = y-1 if y > 0 else NY-1
    return 2*cfg[x, y]*(cfg[x, ym1] + cfg[x, yp1] + cfg[xm1, y] + cfg[xp1, y])

def evolve(cfg, energy, beta, nsweep, obs):
    """
    Evolve a spin configuration in Monte-Carlo time by flipping spins
    at random sites nsweep*NX*NY times and accepting or rejecting the change
    using the Metropolis-Hastings algroithm.
    Measures observables every NX*NY steps, i.e. once per sweep.

    Parameters:
        cfg: Initial configuration.
        energy: Initial energy.
        beta: Inverse temperature.
        nsweep: Number of spin flips sweeps over the lattice.
        obs: Instance of Observables to store energy and magnetisation.
             Can be None in which case nothings gets stored.

    Returns:
        - Final configuration.
        - Final energy.
        - Number of accepted spin flips.
    """

    # copy the config so we don't change it on the outside
    cfg = cfg.copy()
    # track number of accepted changes
    naccept = 0
    # store lattice size for convenience
    nx, ny = cfg.shape
    latsize = nx*ny

    for _ in range(nsweep):
        for _ in range(latsize):
            # pick a random lattice site
            x = np.random.randint(0, nx)
            y = np.random.randint(0, ny)

            # compute proposed change in energy
            delta = delta_e(cfg, x, y)

            # test if change is accepted
            # The first check is not necessary for this to be correct but avoids
            # evaluating the costly exponential and RNG.
            if delta <= 0 or math.exp(-beta*delta) > np.random.uniform(0, 1):
                cfg[x, y] = -cfg[x, y]  # apply the accepted change
                energy += delta
                naccept += 1

        # save energy and magnetisation is reqested
        if obs is not None:
            obs.energy.append(energy)
            # We don't need the magnetisation during the sweep, so don't keep track
            # of it but recompute in the end, it is cheap.
            obs.magnetisation.append(magnetisation(cfg))

    return cfg, energy, naccept

def main():
    """
    Simulate a spin 1/2 lattice using the Ising model.
    """

    np.random.seed(SEED)

    datadir = prepare_datadir()

    # hot start
    cfg = np.random.choice([-1, 1], (NX, NY), replace=True)
    energy = hamiltonian(cfg)

    # start timing now, the above doesn't really count
    start_time = time.time()

    # run initial thermalization
    cfg, energy, naccept = evolve(cfg, energy, 1/TEMPERATURES[0], NTHERM_INIT, None)
    print(f"Initial thermalisation acceptance rate: {naccept/NTHERM_INIT/NX/NY}")

    for itemp, temperature in enumerate(TEMPERATURES):
        print(f"Running T = {temperature}")
        beta = 1 / temperature

        # thermalise
        cfg, energy, naccept = evolve(cfg, energy, beta, NTHERM, None)
        print(f"  Thermalisation acceptance rate: {naccept/NTHERM/NX/NY}")

        # measure
        obs = Observables([], [])
        cfg, energy, naccept = evolve(cfg, energy, beta, NPROD, obs)
        print(f"  Production acceptance rate: {naccept/NPROD/NX/NY}")

        # save result
        np.savetxt(str(datadir/f"{itemp}.dat"), (obs.energy, obs.magnetisation))

    end_time = time.time()
    print(f"Duration in wall clock time: {end_time-start_time:.3f}s")

if __name__ == "__main__":
    main()
