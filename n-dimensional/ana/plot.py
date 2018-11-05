import sys
from pathlib import Path
import re

import numpy as np
import matplotlib.pyplot as plt

class Metadata:
    def __init__(self, J, h, shape):
        self.J = J
        self.h = h
        self.shape = shape

def loadMetadata(fname):
    paramStr = open(fname, "r").readline()

    match = re.match(r"# J=([^ ]+) h=([^ ]+) shape=\[(\d+(, \d+)*)\]", paramStr)
    if not match or len(match.groups()) < 4:
        raise RuntimeError(f"Invalid first line in file {fname}: '{paramStr}'")

    return Metadata(float(match[1]),
                    float(match[2]),
                    [int(d) for d in match[3].split(",")])

def loadFile(fname):
    meta = loadMetadata(fname)
    energy, magn = np.loadtxt(fname, skiprows=1, delimiter=",")
    return meta, energy, magn


def main():
    if len(sys.argv) != 2:
        print("You need to provide the data directory as a command line argument!")
        sys.exit(1)

    datadir = Path(sys.argv[1])

    energy = []
    magnetisation = []
    h = []
    for fname in datadir.iterdir():
        meta, e, m = loadFile(fname)
        e = np.mean(e)
        m = np.mean(m)
        energy.append(e)
        magnetisation.append(np.abs(m))
        h.append(meta.h)

    plt.figure()
    plt.title("energy")
    plt.plot(h, energy, ls="", marker="x")
    plt.figure()
    plt.title("magn")
    plt.plot(h, magnetisation, ls="", marker="x")

    plt.show()

if __name__ == "__main__":
    main()
