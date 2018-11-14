import re
from operator import mul
from functools import reduce

import numpy as np

class Metadata:
    "Store metadata on ensembles."

    def __init__(self, J, h, shape):
        self.J = J
        self.h = h
        self.shape = shape

    def latsize(self):
        return reduce(mul, self.shape)

def loadMetadata(fname):
    "Load metadata from first line of file."

    paramStr = open(fname, "r").readline()

    match = re.match(r"# J=([^ ]+) h=([^ ]+) shape=\[(\d+(, \d+)*)\]", paramStr)
    if not match or len(match.groups()) < 4:
        raise RuntimeError(f"Invalid first line in file {fname}: '{paramStr}'")

    return Metadata(float(match[1]),
                    float(match[2]),
                    [int(d) for d in match[3].split(",")])

def loadFile(fname):
    "Load meta- and 'normal' data from file."
    return loadMetadata(fname), \
        np.loadtxt(fname, skiprows=1, delimiter=",")
