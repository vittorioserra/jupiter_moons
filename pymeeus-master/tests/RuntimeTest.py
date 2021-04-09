import time
from statistics import mean

from pymeeus_optimized.JupiterMoons import JupiterMoons as JM_opt
from pymeeus.JupiterMoons import JupiterMoons as JM
from pymeeus.Epoch import Epoch as Epoch
from pymeeus_optimized.Epoch import Epoch as Epoch_opt

import numpy as np


def unoptimizedCalc(epchs):
    times = []
    starttime = time.time()
    for e in epchs:
        stime = time.time()
        JM.rectangular_positions_jovian_equatorial(e)
        times.append(time.time() - stime)

    times.pop(0)
    print("Unoptimized time [s]: ", time.time() - starttime, "\nAvg (without first): ", mean(times))


def optimizedCalc(epchs):
    times = []
    starttime = time.time()
    for e in epchs:
        stime = time.time()
        JM_opt.rectangular_positions_jovian_equatorial(e)
        times.append(time.time() - stime)

    times.pop(0)
    print("Optimized time [s]: ", time.time() - starttime, "\nAvg (without first): ", mean(times))


def vecCall(epchs):
    starttime = time.time()
    moon_coord = np.vectorize(JM_opt.rectangular_positions_jovian_equatorial, otypes=[tuple])
    arr = moon_coord(epchs)
    print("Optimized time [s]: ", time.time() - starttime)


def unvecCall(epchs):
    starttime = time.time()
    for e in epchs:
        JM_opt.rectangular_positions_jovian_equatorial(e)
    print("Unoptimized time [s]: ", time.time() - starttime)


if __name__ == "__main__":
    epochs = []
    epochs_opt = []
    for t in range(2460000, 2461000):
        epochs.append(Epoch(t))
        epochs_opt.append(Epoch_opt(t))

    #unoptimizedCalc(epochs)
    optimizedCalc(epochs_opt)

    # unvecCall(epochs_opt)
    # vecCall(epochs_opt)
