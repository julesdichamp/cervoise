#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wen Mar 25 12:35:45
@author: florian
"""
# pyximport.install(pyimport=True)
import sys
import os
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
import parser

import momepy
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import networkx as nx

from getFranceGeomDpts import *
from models import *
from plotModels import *

print(
    "\nCERVOISE = Covid19 : Estimation of Regional Virus Organization Infered using Stochastic Equation = Starting"
)
t0 = time.time()

# =============================================================================
# read input file. Either specified by argument or generic inputParams.ini
# =============================================================================
if len(sys.argv) > 1:
    inputFile = sys.argv[1]
else:
    inputFile = "config.ini"


def LoadRoads(road_names):
    """
	Test function to display main roads in France metropole
	In: file name with roads, should be 'road.shp'
	"""
    print("LoadRoads")
    shapefile = gpd.read_file(road_names)
    shapefile = shapefile.drop(shapefile[shapefile["geometry"] == None].index)
    f, ax = plt.subplots(figsize=(10, 10))
    shapefile.plot(ax=ax)
    ax.set_axis_off()
    plt.show()


def InitParams(config_name):
    """
	in: config file name 
	out: dictionnary with all needed parameters
	"""
    print("InitParams")
    config = configparser.ConfigParser()
    config.read(config_name)
    params = {}
    #  load input data infos
    params["graph_name"] = config.get("inputs", "graph_name")
    params["data_fname"] = config.get("inputs", "data_name")
    params["roads_name"] = config.get("inputs", "roads_name")
    params["shapefile_name"] = config.get("inputs", "shapefile_name")
    params["gareFile_name"] = config.get("inputs", "gareFile_name")
    params["ageFile_name"] = config.get("inputs", "ageFile_name")

    # load options infos
    params["location"] = config.get("options", "location")

    # load model infos
    params["model_name"] = config.get("model", "model_name")
    if params["model_name"] == "SIRWellMixed":  # load parameters needed for this model
        params["k_SI"] = config.getfloat("model", "k_SI")
        params["k_IR"] = config.getfloat("model", "k_IR")
        params["S_ini"] = config.getfloat("model", "S_ini")
        params["I_ini"] = config.getfloat("model", "I_ini")
        params["R_ini"] = config.getfloat("model", "R_ini")

    else:
        print("[ERROR]: model name not found. Exit here.")
        exit(0)

    # load settigns infos
    params["time_start"] = config.getfloat("settings", "time_start")
    params["time_end"] = config.getfloat("settings", "time_end")

    #  load outputs infos
    params["outputs_dir"] = config.get("outputs", "outputs_dir")
    if not (
        os.path.exists(params["outputs_dir"])
    ):  #  if directory does not exist, create it
        os.mkdir(params["outputs_dir"])

    return params


def InitGraph(params):
    print("InitGraph")
    """
	initialize a graph of regions with positions, and ids per regions
	in: params, dictionnary with all needed parameters
	out: a graph of type networkx.Graph() 
	"""
    LIST_METROPOLE = [
        "Auvergne-Rhône-Alpes",
        "Bourgogne-Franche-Comté",
        "Bretagne",
        "Centre-Val de Loire",
        "Corse",
        "Grand Est",
        "Hauts-de-France",
        "Ile-de-France",
        "Normandie",
        "Nouvelle-Aquitaine",
        "Occitanie",
        "Pays de la Loire",
        "Provence-Alpes-Côte d’Azur",
    ]

    LIST_OUTREMER = [
        "Guadeloupe",
        "Guyane",
        "Martinique",
        "Mayotte",
        "La Réunion",
        "Saint-Barthélémy",
        "Saint-Martin",
    ]

    # from Florian
    DIC_POSITIONS_METROPOLE = {
        "Auvergne-Rhône-Alpes": [772022.9490728247, 2059666.8371020914],
        "Bourgogne-Franche-Comté": [787177.8408282029, 2251335.950928572],
        "Bretagne": [215039.24516678863, 2365217.5071779615],
        "Centre-Val de Loire": [550917.5994336334, 2275882.307395058],
        "Corse": [1160476.2800258105, 1706885.7403624766],
        "Grand Est": [841720.2901229224, 2415159.5340459584],
        "Hauts-de-France": [631486.0019718618, 2552427.4478215473],
        "Ile-de-France": [612377.0469593473, 2412305.1064892523],
        "Normandie": [437037.0999347288, 2460471.8503777585],
        "Nouvelle-Aquitaine": [431896.2266105056, 2023453.7858529706],
        "Occitanie": [583863.6352416947, 1855575.2526474672],
        "Pays de la Loire": [361717.2237577763, 2279805.18747009],
        "Provence-Alpes-Côte d’Azur": [899251.0580184147, 1890939.1725853563],
    }

    DIC_POSITIONS_OUTREMER = {
        "Guadeloupe": [16.2490067, -61.5650444],
        "Guyane": [4.8417097, -58.6416891],
        "Martinique": [14.6113732, -60.9620777],
        "Mayotte": [-12.823048, 45.1520755],
        "La Réunion": [-21.1309332, 55.5265771],
        "Saint-Barthélémy": [17.9036287, -62.811568843006896],
        "Saint-Martin": [18.06481835, -63.049171129033],
    }

    # dictionnary that associates an id to a region
    dic_regions_metropole = {}
    dic_regions_metropole_positions = {}
    for i in range(0, len(LIST_METROPOLE)):
        dic_regions_metropole[i] = LIST_METROPOLE[i]
        dic_regions_metropole_positions[i] = DIC_POSITIONS_METROPOLE[LIST_METROPOLE[i]]

    dic_regions_outremer = {}
    dic_regions_outremer_positions = {}
    for i in range(0, len(LIST_OUTREMER)):
        dic_regions_outremer[len(dic_regions_metropole) + i] = LIST_OUTREMER[i]
        dic_regions_outremer_positions[
            len(dic_regions_metropole) + i
        ] = DIC_POSITIONS_OUTREMER[LIST_OUTREMER[i]]

    #  set full dictionnary depending on option
    dic_regions = {}
    dic_regions_positions = {}
    if params["location"] == "metropole":
        dic_regions = dic_regions_metropole
        dic_regions_positions = dic_regions_metropole_positions
    elif params["location"] == "outremer":
        dic_regions = dic_regions_outremer
        dic_regions_positions = dic_regions_outremer_positions
    else:
        for x in dic_regions_metropole.keys():
            dic_regions[x] = dic_regions_metropole[x]
            dic_regions_positions[x] = dic_regions_metropole_positions[x]
        for x in dic_regions_outremer.keys():
            dic_regions[x] = dic_regions_outremer[x]
            dic_regions_positions[x] = dic_regions_outremer_positions[x]

            # add nodes of graph
    G = nx.Graph()
    for x in dic_regions.keys():
        G.add_node(x)

    #  All neighborings connections ONLY FOR METROPOLE NOW
    #  Corse to be better dealt with?
    dic_all_connections = {}
    dic_all_connections["Bretagne"] = ["Normandie", "Pays de la Loire"]
    dic_all_connections["Normandie"] = [
        "Bretagne",
        "Pays de la Loire",
        "Centre-Val de Loire",
        "Ile-de-France",
        "Hauts-de-France",
    ]
    dic_all_connections["Pays de la Loire"] = [
        "Bretagne",
        "Normandie",
        "Centre-Val de Loire",
        "Nouvelle-Aquitaine",
    ]
    dic_all_connections["Centre-Val de Loire"] = [
        "Pays de la Loire",
        "Normandie",
        "Ile-de-France",
        "Bourgogne-Franche-Comté",
        "Auvergne-Rhône-Alpes",
        "Nouvelle-Aquitaine",
    ]
    dic_all_connections["Hauts-de-France"] = ["Grand Est", "Ile-de-France", "Normandie"]
    dic_all_connections["Nouvelle-Aquitaine"] = [
        "Pays de la Loire",
        "Centre-Val de Loire",
        "Auvergne-Rhône-Alpes",
        "Occitanie",
    ]
    dic_all_connections["Grand Est"] = [
        "Hauts-de-France",
        "Ile-de-France",
        "Bourgogne-Franche-Comté",
    ]
    dic_all_connections["Bourgogne-Franche-Comté"] = [
        "Ile-de-France",
        "Grand Est",
        "Centre-Val de Loire",
        "Auvergne-Rhône-Alpes",
    ]
    dic_all_connections["Auvergne-Rhône-Alpes"] = [
        "Occitanie",
        "Nouvelle-Aquitaine",
        "Centre-Val de Loire",
        "Bourgogne-Franche-Comté",
        "Provence-Alpes-Côte d’Azur",
    ]
    dic_all_connections["Provence-Alpes-Côte d’Azur"] = [
        "Auvergne-Rhône-Alpes",
        "Occitanie",
        "Corse",
    ]
    dic_all_connections["Occitanie"] = [
        "Nouvelle-Aquitaine",
        "Auvergne-Rhône-Alpes",
        "Provence-Alpes-Côte d’Azur",
    ]
    dic_all_connections["Corse"] = ["Provence-Alpes-Côte d’Azur"]
    dic_all_connections["Ile-de-France"] = LIST_METROPOLE

    for key in dic_all_connections.keys():
        # get id of region
        id1 = list(dic_regions.keys())[list(dic_regions.values()).index(key)]
        for value in dic_all_connections[key]:  # loop over all connections regions
            #  get id of conected region
            id2 = list(dic_regions.keys())[list(dic_regions.values()).index(value)]
            edge = (id1, id2)
            G.add_edge(*edge)  # unpack edge tuple*

    # fig, ax = plt.subplots()
    # nx.draw(G, with_labels=True, font_weight='bold', pos=dic_regions_positions, labels=dic_regions)
    # plt.show()

    return G


def LaunchModel(params):
    if params["model_name"] == "SIRWellMixed":
        outputs = SIRWellMixed(params)
        plot_SIRWellMixed(params, outputs)
    elif params["model_name"] == "SIROnRegions":
        SIRonRegions(params)
    elif params["model_name"] == "SIROnDepartments":
        SIROnDepartments(params)
    else:
        print("[ERROR]: invalid name in parser for model name. Exit here.\n")
        exit(0)


# =============================================================================
# LAUNCH WORKFLOW
# =============================================================================
params = InitParams(inputFile)
# franceGeom = getFranceGeomDpts(params)
# graph = InitGraph(params)
# LoadRoads(params['roads_name'])
LaunchModel(params)
