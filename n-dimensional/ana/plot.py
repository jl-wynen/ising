import sys
from pathlib import Path
from collections import Counter
from operator import itemgetter

import numpy as np
import matplotlib.pyplot as plt

from fileio import loadFile

def findVarying(params):
    varying = [not np.all(p == p[0]) for p in params]
    if Counter(varying)[True] == 0:
        raise RuntimeError("No varying parameters")
    if Counter(varying)[True] > 1:
        raise RuntimeError("More than one parameter is varying")
    return varying.index(True)

def splitVarFixed(meta):
    # extrac individual params from meta object
    decomposed = list(zip(*((m.shape, m.J, m.h) for m in meta)))
    kwparams = dict(shape=np.array(decomposed[0]),
                    J=np.array(decomposed[1]),
                    h=np.array(decomposed[2]))

    # find parameter that varies
    names, params = zip(*kwparams.items())
    varying = findVarying(params)
    varName, varVal = names[varying], params[varying]

    # remove varying from dict
    del kwparams[varName]
    # keep only fist elem (the rest are the same anyway)
    fixed = {name: val[0] for name, val in kwparams.items()}

    return (varName, varVal), fixed

def xsorted(x, y):
    return tuple(zip(*sorted(zip(x, y), key=itemgetter(0))))

def main():
    if len(sys.argv) != 2:
        print("You need to provide the data directory as a command line argument!")
        sys.exit(1)

    datadir = Path(sys.argv[1])

    meta, energy, magn = zip(*((meta, np.mean(dat[0]), np.mean(dat[1]))
                               for meta, dat in map(loadFile, filter(lambda d: d.suffix == ".dat",
                                                                     datadir.iterdir()))))
    relEnergy = [e/m.latsize() for e, m in zip(energy, meta)]
    (xlabel, x), fixed = splitVarFixed(meta)

    fig = plt.figure(figsize=(11, 5))
    fig.suptitle(",   ".join("{}={}".format(n, v) for n, v in fixed.items()))

    axe = fig.add_subplot(121)
    axe.set_xlabel(xlabel)
    axe.set_ylabel(r"$E / \Lambda$")
    axe.plot(*xsorted(x, relEnergy), ls="-", marker=".")

    axm = fig.add_subplot(122)
    axm.set_xlabel(xlabel)
    axm.set_ylabel(r"$m$")
    axm.plot(*xsorted(x, magn), ls="-", marker=".")

    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()

if __name__ == "__main__":
    main()
