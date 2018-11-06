import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

from fileio import loadFile

def main():
    if len(sys.argv) != 2:
        print("You need to provide the data directory as a command line argument!")
        sys.exit(1)

    datadir = Path(sys.argv[1])

    energy = []
    magnetisation = []
    h = []
    for fname in [d for d in datadir.iterdir() if d.suffix=='.dat']:
        meta, dat = loadFile(fname)
        e = np.mean(dat[0])
        m = np.mean(dat[1])
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
