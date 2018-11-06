"""
Show an animation of configurations as they evolve.

Command line arguments are input-dir, [ensemble number],
where ensemble number is optional and defaults to 0.
"""

import sys
from pathlib import Path
from functools import reduce

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import matplotlib.gridspec as gridspec

from fileio import loadFile

def main():
    if len(sys.argv) < 2:
        print("You need to provide the data directory as a command line argument!")
        sys.exit(1)

    datadir = Path(sys.argv[1])
    ensemble = int(sys.argv[2]) if len(sys.argv) == 3 else 0


    # load data
    meta, cfgs = loadFile(datadir/f"{ensemble:04d}.cfg")
    _, data = loadFile(datadir/f"{ensemble:04d}.dat")
    energy = data[0]/reduce(lambda x, y: x*y, meta.shape, 1)
    magn = data[1]

    nframes = cfgs.shape[0]

    # make figure
    fig = plt.figure(figsize=(7.5, 5))
    gspec = gridspec.GridSpec(2, 3)

    # create image
    axim = fig.add_subplot(gspec[0:2, 0:2])
    image = axim.imshow(-np.ones(meta.shape), vmin=-1, vmax=+1,
                        cmap="Greys", interpolation="nearest")
    axim.set_title(fr"t = 000 / ?, J = {meta.J}, h = {meta.h}")

    # create plot for energy
    axe = fig.add_subplot(gspec[0, 2])
    epoints, = axe.plot(np.arange(nframes), energy, c="C0")
    epoints.set_data([], [])
    axe.set_xlabel(r"$N_\mathrm{tr}$")
    axe.set_ylabel(r"$H(s) / L^2$")

    # create plot for magnetisation
    axm = fig.add_subplot(gspec[1, 2])
    mpoints, = axm.plot(np.arange(nframes), magn, c="C1")
    mpoints.set_data([], [])
    axm.set_xlabel(r"$N_\mathrm{tr}$")
    axm.set_ylabel(r"$m(s)$")

    fig.tight_layout()

    # update image and plots for given frame
    def animate(frame):
        axim.set_title(fr"t = {frame:03d} / {nframes},   J = {meta.J}, h = {meta.h}")
        image.set_data(np.reshape(cfgs[frame], meta.shape))
        epoints.set_data(np.arange(frame), energy[:frame])
        mpoints.set_data(np.arange(frame), magn[:frame])

    # animate the whole shebang
    anim = animation.FuncAnimation(fig, animate, frames=nframes,
                                   interval=200)

    # save to file
    # anim.save("anim.mp4", fps=7, extra_args=["-vcodec", "libx264"])

    plt.show()

if __name__ == "__main__":
    main()
