#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wen Mar 25 12:35:45
@author: florian
"""
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

import numpy as np
import scipy as sp

import numexpr
import json
import scipy.stats as ss
import configparser
import re
from transonic import boost

# sys.path.append("/home/florian/PycharmProjects/drunkWalker/")
# import computeMSD

print("\ncovid19 epidemic modeling with master equation stochastic model")
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
