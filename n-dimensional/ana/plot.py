import sys
from pathlib import Path
from collections import Counter
from operator import itemgetter

import numpy as np
import matplotlib.pyplot as plt

from fileio import loadFile

def findVarying(params):
    "Return index of varying parameter. There must be exactly 1, the others must be fixed."

    varying = [not np.all(p == p[0]) for p in params]
    if Counter(varying)[True] == 0:
        raise RuntimeError("No varying parameters")
    if Counter(varying)[True] > 1:
        raise RuntimeError("More than one parameter is varying")
    return varying.index(True)

def splitVarFixed(meta):
    "Identify varying parameter and split it off from the rest."

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
    "Sort both x and y based on x."
    return tuple(zip(*sorted(zip(x, y), key=itemgetter(0))))

def meanAndErr(x, nbs=100, bslength=None):
    """
    Compute mean and error of x, the latter via bottstrap
    Parameters:
      x: 1D array of data.
      nbs: Number of bootstrap samples
      bslength: Size of each bootstrap sample. Defaults to len(x)
    """

    if not bslength:
        bslength = len(x)
    bootstrapIndices = np.random.randint(0, len(x), [nbs, bslength])

    return np.mean(x), np.std(np.mean(x[bootstrapIndices], axis=1), axis=0)

def loadData(datadir, skip=1):
    """
    Load parameters, energies, and magnetisations from all datafiles in directory.
    Parameters:
      datadir: Input directory. All files datadir/*.dat are loaded.
      skip: Stride for data arrays. Uses only every skip-th entry.
    """

    metas = []
    energies = []
    magns = []
    for m, dat in map(loadFile, filter(lambda d: d.suffix == ".dat", datadir.iterdir())):
        metas.append(m)
        energies.append(meanAndErr(dat[0, ::skip]/m.latsize()))
        magns.append(meanAndErr(np.abs(dat[1, ::skip])))

    (xlabel, x), fixed = splitVarFixed(metas)
    return (xlabel, x), fixed, energies, magns


def main():
    if len(sys.argv) != 2:
        print("You need to provide the data directory as a command line argument!")
        sys.exit(1)
    datadir = Path(sys.argv[1])

    np.random.seed(123)

    (xlabel, x), fixed, ener, magn = loadData(datadir, skip=20)

    fig = plt.figure(figsize=(11, 5))
    fig.suptitle(",   ".join("{}={}".format(n, v) for n, v in fixed.items()))

    axe = fig.add_subplot(121)
    axe.set_xlabel(xlabel)
    axe.set_ylabel(r"$E / \Lambda$")
    xs, es = xsorted(x, ener)
    axe.errorbar(xs, *zip(*es), ls="-", marker=".")

    axm = fig.add_subplot(122)
    axm.set_xlabel(xlabel)
    axm.set_ylabel(r"$|m|$")
    xs, ms = xsorted(x, magn)
    axm.errorbar(xs, *zip(*ms), ls="-", marker=".")

    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()

if __name__ == "__main__":
    main()
