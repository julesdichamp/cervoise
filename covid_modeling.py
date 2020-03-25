#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 12 17:13:07 2019
@author: florian
"""
# import pyximport
# pyximport.install(pyimport=True)
import sys
import bisect
import fastrand
import multiprocessing
from multiprocessing import Queue
from queue import Empty
import time
import progressbar
import matplotlib.pyplot as plt
import seaborn as sns
import vapeplot

# from vmtk import vmtkscripts
import numpy as np
import scipy as sp

# import diffusion_analysis as da
import numexpr
import json
import scipy.stats as ss
import configparser
import re
from transonic import boost

# sys.path.append("/home/florian/PycharmProjects/drunkWalker/")
# import computeMSD

print("\nrunning 2D case diffusion on a reg grid")
t0 = time.time()

# =============================================================================
# read input file. Either specified by argument or generic inputParams.ini
# =============================================================================

Config = configparser.ConfigParser()
if len(sys.argv) > 1:
    inputFile = sys.argv[1]
else:
    inputFile = "inputParams.ini"

Config.read(inputFile)

# =============================================================================
# feed the dictionnary
# =============================================================================
def ConfigSectionMap():
    dict1 = {}
    sections = ["Domain", "Properties", "Others"]
    for section in sections:
        options = Config.options(section)

        for option in options:
            try:
                dict1[option] = Config.get(section, option)
                if dict1[option] == -1:
                    DebugPrint("skip: %s" % option)
            except:
                print("exception on %s!" % option)
                dict1[option] = None
    return dict1


params = {
    "Nx": 500,  # grid size
    "Ny": 500,  # 1 for 2D
    "lx": 1e-6 / 1,  # spacing
    "D": 3e-12,  # diffusion Coeff (m2/s)
    "simTime": 25,  # physical time
    "nParticles": np.int(500000),  # number of particls
    "source": "disk",  # source type manual or disk or pt
    "diskSourceRadius": 25, #source radius
    "refX": 128,  # position on x of line activation
    "obstacles": False,  # add obstacles
    "network": True,
    "obsLength": 5,
    "obstacleSpacing": 10,
    "nCpus": 1,  # nb of process
    "radiusROI": 100,  # radius of the ROI for t_1/2 in 2d
    #    "drift": [1, 2., 1, 1],
    "Cx": 250,
    "Cy": 250,
    "onlyDomain": 1,  # stop after generating the domain
    "dimMSD": 2,  # compute MSD in 1 or 2D
    "msdFileName": "msd_ny45_obstacles",
}

# params = ConfigSectionMap()
# =============================================================================
# recompute the number of time intervals for resampling at the end
# =============================================================================
params["dt"] = params["lx"] ** 2 / (2 * params["D"])
params["timeIntervals"] = int(params["simTime"] / params["dt"])

print("Simulation parameters are:\n")
print(json.dumps(params, indent=4, sort_keys=0))

print("Physical domain is :", params["lx"] * params["Nx"], params["lx"] * params["Ny"])


def logSpaced(L, pointsPerDecade=15):
    """Generate an array of log spaced integers smaller than L"""
    nbdecades = np.log10(L)
    return np.unique(
        np.logspace(
            start=0,
            stop=nbdecades,
            num=nbdecades * pointsPerDecade,
            base=10,
            endpoint=False,
        ).astype(int)
    )


@boost
def rr(i, n):
    """ generate list of random int of length i bwtween 0 and n"""
    return fastrand.pcg32bounded(n) * i


def pair_from_y(y, ROImean):
    x = np.arange(0, len(ROImean))
    y_lookup = sorted(map(tuple, map(reversed, zip(x, ROImean))))

    return tuple(
        reversed(y_lookup[min(bisect.bisect_left(y_lookup, (y,)), len(y_lookup) - 1)])
    )


@boost
def get_msdTraj(traj, dt, with_nan=False):

    shifts = np.arange(1, len(traj), dtype="int")
    shifts = logSpaced(len(traj), 20)
    msd = np.empty((len(shifts), 2), dtype="float")
    msd[:] = np.nan

    msd[:, 1] = shifts * dt

    for i, shift in enumerate(shifts):
        diffs = traj[:-shift] - traj[shift:]
        if with_nan:
            diffs = diffs[~np.isnan(diffs).any(axis=1)]
        diffs = np.square(diffs).sum(axis=0)

        if len(diffs) > 0:
            msd[i, 0] = np.mean(diffs)

    #    msd = pd.DataFrame(msd)
    #    msd.columns = ["msd", "delay"]
    #
    #    msd.set_index('delay', drop=True, inplace=True)
    #    msd.dropna(inplace=True)

    return msd[:, 0]


@boost
def get_msd(traj, dt, with_nan=False):

    #    shifts = np.arange(1, len(traj), dtype='int')
    shifts = logSpaced(len(traj) * params["dt"], 30)
    msd = np.empty((len(shifts), 2), dtype="float")
    msd[:] = np.nan

    msd[:, 1] = shifts * dt

    for i, shift in enumerate(shifts):
        diffs = traj[:-shift] - traj[shift:]
        if with_nan:
            diffs = diffs[~np.isnan(diffs).any(axis=1)]
        diffs = np.square(diffs).sum(axis=0)

        if len([diffs]) > 0:
            msd[i, 0] = np.mean(diffs)

    #    msd = pd.DataFrame(msd)
    #    msd.columns = ["msd", "delay"]
    #
    #    msd.set_index('delay', drop=True, inplace=True)
    #    msd.dropna(inplace=True)

    return msd[:, 1]


@boost
def findNeighbours(i, j):
    # =============================================================================
    #     first find bc list for the coords in bcGrid, called n
    #     from n, find possible directions, called possibilities
    #
    # =============================================================================
    # 1D case
    if params["Ny"] == 1:
        j = 0
        if bcGrid[i, j] == 1:
            possibilities = [0, 1, 0, 0]
        elif bcGrid[i, j] == 2:
            possibilities = [1, 0, 0, 0]
        else:
            possibilities = [1, 1, 0, 0]

    else:
        # find neighbours bc
        # 2D case
        if i == 0:  # bord
            if j == 0:  # coin 1
                n = [0, bcGrid[i + 1, j], bcGrid[i, j + 1], 0]
            elif j == params["Ny"] - 1:  # coin 1
                n = [0, bcGrid[i + 1, j], 0, bcGrid[i, j - 1]]
            else:  # bord sans coin
                n = [0, bcGrid[i + 1, j], bcGrid[i, j + 1], bcGrid[i, j - 1]]
        elif i == params["Nx"] - 1:
            if j == 0:  # coin 1
                n = [bcGrid[i - 1, j], 0, bcGrid[i, j + 1], 0]
            elif j == params["Ny"] - 1:  # coin 1
                n = [bcGrid[i - 1, j], 0, 0, bcGrid[i, j - 1]]
            else:  # bord sans coin
                n = [bcGrid[i - 1, j], 0, bcGrid[i, j + 1], bcGrid[i, j - 1]]
        elif j == params["Ny"] - 1:
            n = [bcGrid[i - 1, j], bcGrid[i + 1, j], 0, bcGrid[i, j - 1]]
        elif j == 0:
            n = [bcGrid[i - 1, j], bcGrid[i + 1, j], bcGrid[i, j + 1], 0]
        else:
            n = [bcGrid[i - 1, j], bcGrid[i + 1, j], bcGrid[i, j + 1], bcGrid[i, j - 1]]

        # find where I am
        if bcGrid[i, j] == 0:
            possibilities = [1, 1, 1, 1]
        elif bcGrid[i, j] == 1:
            possibilities = [1, 1, 1, 0]
        elif bcGrid[i, j] == 2:
            possibilities = [1, 1, 0, 1]
        elif bcGrid[i, j] == 4:
            possibilities = [1, 1, 0, 0]
        elif bcGrid[i, j] == 3:
            possibilities = [0, 0, 0, 0]
        #    elif bcGrid[i, j] == 5:
        #        possibilities = [0, 1, 1, 1]
        #    elif bcGrid[i, j] == 6:
        #        possibilities = [1, 0, 1, 1]
        # all positions of wall
        g = list(i for i, e in enumerate(n) if e == 3)
        for direction in g:
            possibilities[direction] = 0
        #    print(possibilities)

    return np.array(possibilities)


# =============================================================================
# # init progressbar
# =============================================================================
widgets = [
    " [",
    progressbar.Timer(),
    "] ",
    progressbar.Bar(),
    " (",
    progressbar.ETA(),
    ") ",
]
bar = progressbar.ProgressBar(
    min_value=0, max_value=params["nParticles"], widgets=widgets
)
# =============================================================================
# # grid of initial conditions
# =============================================================================
grid = np.zeros([params["Nx"], params["Ny"]], dtype=np.uint32)

if params["source"] == "disk":
    # =============================================================================
    #     # source position in grid
    # =============================================================================


    circleIni = np.zeros([params["Nx"], params["Ny"]], dtype=np.uint32)
    xIni = np.arange(0, params["Nx"])
    yIni = np.arange(0, params["Ny"])
    circleIni[
        (xIni[np.newaxis, :] - params["Cx"]) ** 2
        + (yIni[:, np.newaxis] - params["Cy"]) ** 2
        < params["diskSourceRadius"] ** 2
    ] = 1
    # =============================================================================
    #     # nb of activated cells
    # =============================================================================
    nIni = np.sum(circleIni)
    circleIni *= int(params["nParticles"] / nIni)
#    grid[params["Cx"], params["Cy"]] = params["nParticles"]
    grid = circleIni

elif params["source"] == "point":
    grid[params["Cx"], params["Cy"]] = params["nParticles"]
elif params["source"] == "line":
    grid[int(params["refX"]), :] = int(params["nParticles"] / params["Ny"])
#    grid[int(params["refX"]+1), :] = int(params["nParticles"] / params["Ny"])


if params["network"]:
    bcGrid = np.ones([params["Nx"], params["Ny"]], dtype=np.uint32) * 3
    for i in range(params["Nx"]):
        for j in range(params["Ny"]):
            if i % params["obsLength"] == 0 or j % params["obsLength"] == 0:
                bcGrid[i, j] = 0
    grid[bcGrid==3] = 0
else:
    bcGrid = np.zeros_like(grid)

# =============================================================================
# # grid of BC, 1=top,2=bottom, 3 = inner walls
# =============================================================================

if params["Ny"] == 1:
    bcGrid[:, 0] = 4
    bcGrid[0, 0] = 1
    bcGrid[-1, 0] = 2
else:
    bcGrid[:, 0], bcGrid[:, params["Ny"] - 1] = 1, 2
    bcGrid[0, :], bcGrid[-1, :] = 0, 0
#    bcGrid[0, :], bcGrid[params["Nx"] - 1, :] = 5, 6

if params["obstacles"]:
    for i in range(params["Nx"] - 6):
        for j in range(params["Ny"]):
            if i % int(params["obstacleSpacing"]) == 0 and (
                j < params["obsLength"] or j > params["Ny"] - params["obsLength"] - 1
            ):
                bcGrid[i, j] = 3
                bcGrid[i + 1, j] = 3
                bcGrid[i + 2, j] = 3
                bcGrid[i + 3, j] = 3
                bcGrid[i + 4, j] = 3
                bcGrid[i + 5, j] = 3
                bcGrid[i + 6, j] = 3


# =============================================================================
# # loop to find neighbours
# =============================================================================
neighbours = {}
for i in range(params["Nx"]):
    for j in range(params["Ny"]):
        neighbours[i, j] = findNeighbours(i, j)


@boost
def walkVart(grid):
    print("\nwalker texas ranger started\n")

    posList = []
    partList = []

    # X[i]
    X = np.zeros(params["timeIntervals"] + 50)
    X2 = np.zeros(params["timeIntervals"] + 50)
    #    Y2 = np.zeros(params["timeIntervals"]+50)

    Nevent = np.zeros(params["timeIntervals"] + 50)
    incrPart = 0
    for i in range(params["Nx"]):
        for j in range(params["Ny"]):
            nPartHere = grid[i, j]
            if nPartHere != 0:
                for part in range(nPartHere):
                    bar.update(incrPart)
                    t = 0
                    ii, jj = i, j
                    posList.append([t, [ii, jj]])
                    while t <= params["simTime"]:
                        if not 0 < ii < params["Nx"] - 1:
                            break
                        #                        avail = findNeighbours(ii, jj)
                        avail = neighbours[ii, jj]
                        posDirection = np.sum(avail)
                        if posDirection < 1:
                            break

                        w = 1 * np.sum(avail) * params["D"] / params["lx"] ** 2
                        dt = -(1 / w) * np.log(1 - np.random.rand())

                        t += dt

                        randList = fastrand.pcg32bounded(posDirection)
                        # =============================================================================
                        #                         # for the particle part, going to
                        # =============================================================================
                        goto = np.nonzero(avail)[0][randList]
                        if goto == 0:
                            posList.append([t, [ii - 1, jj]])
                            ii -= 1
                        if goto == 1:
                            posList.append([t, [ii + 1, jj]])
                            ii += 1
                        if goto == 2:
                            posList.append([t, [ii, jj + 1]])
                            jj += 1
                        if goto == 3:
                            posList.append([t, [ii, jj - 1]])
                            jj -= 1

                        #                        print('pos list', t, dt, t/params["dt"], '\n')
                        #                        print('pos list', t, dt, t/params["dt"], '\n')
                        #                        print(posList[-1][1])
                        timid = int(t / params["dt"])
                        if params["dimMSD"] == 2:
                            X[timid] += posList[-1][1][0] + posList[-1][1][1]
                            X2[timid] += (posList[-1][1][0] + posList[-1][1][1]) ** 2
                        else:
                            X[timid] += posList[-1][1][0]
                            X2[timid] += posList[-1][1][0] ** 2

                        Nevent[timid] += 1

                    partList.append(posList)
                    posList = []
                    incrPart += 1
    """
    interpt = np.linspace(0, params["simTime"], params["timeIntervals"])
    t, x, y = [], [], []
    partID = 0
# =============================================================================
#     # toutes les part, tous les temps
# =============================================================================
    interpPos = np.empty([len(partList), 2, len(interpt)])
    for idp, p in enumerate(partList):
        #    print("part", idp)
        for timeStep in p:
            t.append(timeStep[0])
            x.append(timeStep[1][0])
            y.append(timeStep[1][1])
        interpX = np.interp(interpt, t, x)
        interpY = np.interp(interpt, t, y)
        interpPos[idp, 0, :] = interpX
        interpPos[idp, 1, :] = interpY
        #    print("i")
        #    plt.plot(interpX, interpY)
        #    plt.show()
        partID += 1
        t, x, y = [], [], []
    return interpPos
    """
    #    if params["dimMSD"] == 2:
    #        return X*params["lx"], X2*params["lx"]**2,Y*params["lx"], Y2*params["lx"]**2, Nevent.astype(int)
    #    else:
    return X * params["lx"], X2 * params["lx"] ** 2, Nevent.astype(int)


def mp_factorizer(subGridList, nprocs):
    def worker(trajList, out_q):
        # =============================================================================
        #         # put into queue the result of simulations
        # =============================================================================
        returnMSDtmp = get_msd(trajList)
        out_q.put(returnMSDtmp)

    # init of the result queue
    out_q = Queue()
    # init of process list
    procs = []
    # loop over ncpus to start simulations
    for i in range(nprocs):
        p = multiprocessing.Process(
            target=worker, args=(trajList[i], out_q)  # function to be ran on cpus
        )
        procs.append(p)
        print("starting process", i, "on ", np.sum(trajList[i]), "particles")
        p.start()

    # init of arrays to gather all data
    # +1 because t=0 added at the begining of walk
    returnMSD = [np.zeros_like(trajList14 * 14[0])] * (params["simTime"] + 1)
    # process to gather data while the simu is running
    liveprocs = list(procs)
    while liveprocs:
        try:
            while 1:
                returnGridTmp = out_q.get(False)
                print("len tmp", len(returnGridTmp), "len dst", len(returnGrid))
                returnGrid = [x + y for x, y in zip(returnGridTmp, returnGrid)]

        except Empty:
            pass

        time.sleep(0.5)  # Give tasks a chance to put more data in
        if not out_q.empty():
            continue

        liveprocs = [p for p in liveprocs if p.is_alive()]
    # processes are ended
    for p in procs:
        p.join()

    return returnGrid  # , ROImeanList


# if params["onlyDomain"]:
#    sys.exit()

subGrids = [(grid / params["nCpus"]).astype(np.uint32)] * params["nCpus"]
print("starting run on ", params["nCpus"], " processes")
# returnGrid = mp_factorizer(subGrids, params["nCpus"])
interpt = np.linspace(0, params["simTime"], params["timeIntervals"])

X, X2, Nevent = walkVart(subGrids[0])
print("\n\ntime =", time.time() - t0)
x2minusx = (
    X2[: len(interpt)] / Nevent[: len(interpt)]
    - (X[: len(interpt)] / Nevent[: len(interpt)]) ** 2
)


plt.xlabel("Time (s)")
plt.ylabel("MSD (micrometer ** 2)")
deuxdD = ss.linregress(interpt, x2minusx * 1e12)
print("MSD(tau=1s) = {0:.2f} micrometer ** 2/s".format(deuxdD.slope))
if params["dimMSD"] == 2:
    print(
        "The effective diffusion coefficient is {0:.2f} micrometer ** 2/s".format(
            deuxdD.slope / 4
        )
    )
else:
    print(
        "The effective diffusion coefficient is {0:.2f} micrometer ** 2/s".format(
            deuxdD.slope / 2
        )
    )

title = "slope = {0:.2f} micrometer ** 2/s".format(deuxdD.slope)
plt.title(title)
sns.lineplot(x=interpt, y=x2minusx * 1e12)
imgname = "_".join("{!s}={!r}".format(key, val) for (key, val) in params.items())
imgname = re.sub("\W+", "", imgname) + ".png"
plt.savefig(imgname[:100])

plt.show()
# msd = get_msd(x2minusx, params["dt"])


# partList = mp_factorizer(subGrids, params["nCpus"])

# print("SUM", np.sum(returnGrid))
# np.savetxt("ROImean.csv", ROImeanList)

# compute t_1/2
"""
t12 = pair_from_y(np.max(ROImeanList)/2., ROImeanList)
ya, yb = t12[0], t12[0]+1
xa, xb = ROImeanList[ya], ROImeanList[yb]
ya, yb = np.float(ya), np.float(yb)
y12 = np.max(ROImeanList)/2.
t12_interp = ya + (t12 - xa) * (yb - ya) / (xb - xa)
print("t_1/2 = ", t12[0], 'or ', y12 * params["lx"]**2/(2 * params["D"]), ' s')
"""

# interpt = np.linspace(0, params["simTime"], params["timeIntervals"])
# t, x, y = [], [], []
# partID = 0
# toutes les part, tous les temps
# interpPos = np.empty([len(partList), 2, len(interpt)])
# for idp, p in enumerate(partList):
#    #    print("part", idp)
#    for timeStep in p:
#        t.append(timeStep[0])
#        x.append(timeStep[1][0])
#        y.append(timeStep[1][1])
#    interpX = np.interp(interpt, t, x)
#    interpY = np.interp(interpt, t, y)
#    interpPos[idp, 0, :] = interpX
#    interpPos[idp, 1, :] = interpY
#
#    #    print("i")
#    #    plt.plot(interpX, interpY)
#    #    plt.show()
#    partID += 1
#    t, x, y = [], [], []

# =============================================================================
# # compute number of parts in ROI
# =============================================================================
"""
if params["Ny"] == 1:
    nPartInRoi = [0] * params["timeIntervals"]
    for part in range(interpPos.shape[0]):
        posx = interpPos[part, 0, :]
        for it, x in enumerate(posx):
            if params["refX"] + params["radiusROI"]/2. > x > params["refX"] - params["radiusROI"]/2.:
                nPartInRoi[it] += 1
else:
    nPartInRoi = [0] * params["timeIntervals"]
    for part in range(interpPos.shape[0]):
        posx, posy = interpPos[part, 0, :], interpPos[part, 1, :]
        for it, x in enumerate(posx):
            #            if (posx[it]-params["Cx"])**2 + (posy[it]-params["Cy"])**2 < params["radiusROI"]**2:
            if params["refX"] + params["radiusROI"]/2 > x > params["refX"] - params["radiusROI"]/2:
                nPartInRoi[it] += 1
# =============================================================================
# # t_1/2
# =============================================================================
interpt = np.linspace(0, params["simTime"], params["timeIntervals"])
f = sp.interpolate.interp1d(nPartInRoi, interpt)
print("\n\nt_1/2 is ", f(nPartInRoi[0]/2), "s")
"""

"""
# plotting
# sns.set()
var = []
sns.set_style("ticks")
sns.despine()
#sns.set_palette(sns.dark_palette("purple"))
vapeplot.set_palette('sunset_r')
for t in range(params["timeIntervals"]):
    if params["Ny"] == 1:
        var.append(np.var(interpPos[:, 0, t]-125))
        var.append(msd1(interpPos[:, 0, t]))
    else:
#        var.append(np.var(interpPos[:, :, t]-125))
        var.append(msd(interpPos[:, 0, t],interpPos[:, 1, t]))
plt.ylabel("variance")
plt.xlabel("time step")
sns.lineplot(x=interpt, y=var)
m, b = np.polyfit(interpt, var, 1)
print("2dD =", m)
plt.show()
"""

# =============================================================================
# plot trajectories
# =============================================================================
plot_traj = 0
if plot_traj == 1:
    sns.set_style("whitegrid")
    for part in range(int(params["nParticles"] / 100)):
        #        sns.lineplot(x=interpPos[part,0,:], y=interpPos[part,1,:])
        sns.scatterplot(x=interpPos[part, 0, :], y=interpPos[part, 1, :])
    # sns.scatterplot(x=interpPos[10, 0, :], y=interpPos[10, 1, :])

    sns.scatterplot(x=interpPos[0, 0, :], y=interpPos[0, 1, :])
    # plt.xlim(210,230)
    # plt.ylim(210,230)
    plt.show()


"""
for t in range(params["timeIntervals"]):
    posx, posy = [], []
#    for part in range(interpPos.shape[0]):
#        posx.append(interpPos[part, 0, t])
#        posy.append(interpPos[part, 1, t])
    if t % 10 == 0:
#        posy = np.array(posy)
#        posx = np.array(posx)
#        plt.scatter(posx, posy)
#        for part in range(interpPos.shape[0]):
#            if part % 100 == 0:
#                    plt.scatter(interpPos[part,0,t], interpPos[part,1,t])
        sns.distplot(interpPos[:,1,t], hist=False)
#        plt.scatter()
        plt.xlim(0, params["Nx"])
        plt.ylim(0, params["Ny"])
#        sns.jointplot(x=posx, y=posy, kind="kde")
#        sns.distplot(posy, hist=False)
#        plt.show()
#        plt.ylim(0, 0.04)
plt.ylabel("density")
plt.xlabel("x")
plt.show()
"""
"""
# =============================================================================
# MSD
# =============================================================================
posx, posy = [], []
path = []
msd = np.zeros(logSpaced(len(traj), 20).shape)
dt = params["simTime"]/params["timeIntervals"]
tau = 1/dt
for part in range(interpPos.shape[0]):
    # for part in range(50):
    path = []
    for t in range(params["timeIntervals"]):
        path.append([interpPos[part, 0, t]])  # ,interpPos[part, 1, t]])
#    if part == 0:
#        msd = np.zeros(logSpaced(len(path), 20).shape)
    path = np.array(path)
    msd += get_msd(path, dt, with_nan=False)
msd /= params["nParticles"]
np.save(params["msdFileName"], msd)
# =============================================================================
# #msd scaling
# =============================================================================
msdscale = (1e-6/params["lx"])**2
sns.lineplot(x=interpt[:-1], y=msd/msdscale)
plt.xlabel('Time (s)')
plt.ylabel('MSD (micrometer ** 2)')
deuxdD = ss.linregress(interpt[:-1], msd/msdscale)
title = 'slope = ' + str(deuxdD[0]) + ' micrometer ** 2/s'
plt.title(title)
plt.show()
print("\n\ntime =", time.time() - t0)
"""
