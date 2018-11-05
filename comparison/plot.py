import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

def read_temperatures(datadir):
    """
    Read the temperatures file.

    Returns:
        - List of indices correponding to file names for respective.
        - List of temperatures.
    """

    with open(datadir/"temperatures.dat", "r") as tempf:
        idxs, temps = zip(*(line.split(":") for line in tempf))
    return list(map(int, idxs)), \
        list(map(float, temps))

def read_data(datadir, idx):
    """
    Read data for given temperature index.
    """
    return np.loadtxt(datadir/f"{idx}.dat")

def plot_energy(ax, temperatures, energy):
    ax.set_title("Energy")
    ax.plot(temperatures, energy, marker=".")
    ax.set_xlabel("$T$ / J")
    ax.set_ylabel("$E(T)$ / J")

def plot_magnetisation(ax, temperatures, magnetisation):
    ax.set_title("Magnetisation")
    ax.plot(temperatures, magnetisation, marker=".")
    ax.set_xlabel("$T$ / J")
    ax.set_ylabel("$|m(T)|$")

def plot_susceptibility(ax, temperatures, susceptibility):
    ax.set_title("Susceptibility")
    ax.plot(temperatures, susceptibility, marker=".")
    ax.set_xlabel(r"$T$ / J")
    ax.set_ylabel(r"$|\chi(T)|$")

def plot(temperatures, energy, magnetisation, susceptibility):
    fig = plt.figure(figsize=(16, 5))
    plot_energy(fig.add_subplot(131), temperatures, energy)
    plot_magnetisation(fig.add_subplot(132), temperatures, magnetisation)
    plot_susceptibility(fig.add_subplot(133), temperatures, susceptibility)
    fig.tight_layout()

def main():
    if len(sys.argv) != 2:
        print("You need to provide the data directory as a command line argument!")
        sys.exit(1)

    datadir = Path(sys.argv[1])

    idxs, temps = read_temperatures(datadir)

    energy = []
    magnetisation = []
    susceptibility = []
    for idx, temp in zip(idxs, temps):
        beta = 1/temp
        ener, mag = read_data(datadir, idx)

        energy.append(np.mean(ener))
        magnetisation.append(np.mean(np.abs(mag)))
        susceptibility.append(beta*np.var(np.abs(mag)))

    plot(temps, energy, magnetisation, susceptibility)

    plt.show()

if __name__ == "__main__":
    main()
