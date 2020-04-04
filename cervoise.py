#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wen Mar 25 12:35:45
@author: florian
"""
# pyximport.install(pyimport=True)
import sys
import os
import glob
import bisect
import fastrand
import multiprocessing
from multiprocessing import Queue
from queue import Empty
import time
import progressbar
import matplotlib.pyplot as plt
# import cartopy.crs as ccrs

import numpy as np
import scipy as sp

import numexpr
import json
import scipy.stats as ss
import configparser
import re
import parser

# import momepy
import pandas as pd
import geopandas as gpd
import networkx as nx

import getFranceGeom
# from getFranceGeom import *
from models import *
from plotModels import *
import distutils
from distutils.util import strtobool

# contains functions to perform fit of the model
import lmfit
from fit import InitFit, LaunchFit


print(
    "\nCERVOISE = Covid19 : Estimation of Regional Virus Organization Infered using Stochastic Equation: Starting"
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
    shapefile = shapefile.drop(shapefile[shapefile["geometry"] is None].index)
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
    #  load input data infos
    params["graph_name"] = config.get("inputs", "graph_name")
    params["data_dir"] = config.get("inputs", "data_dir")
    params["roads_name"] = config.get("inputs", "roads_name")
    params["shapefile_name"] = config.get("inputs", "shapefile_name")
    params["gareFile_name"] = config.get("inputs", "gareFile_name")
    params["ageFile_name"] = config.get("inputs", "ageFile_name")
    params['airplaneFile_name'] = config.get("inputs", "airplaneFile_name")
    params['citiesFile_name'] = config.get("inputs", 'citiesFile_name')

    params['covid19DepartmentsDir_name'] = config.get("inputs", "covid19DepartmentsDir_name")
    params['covid19RegionsDir_name'] = config.get("inputs", "covid19RegionsDir_name")
    # get all file names we are interested in in this directory (do not keep polygons files)s
    params['allCovid19DepartmentsFileNames'] = glob.glob(
        params['covid19DepartmentsDir_name'] + '/*departements-france-*.json')
    params['allCovid19RegionsFileNames'] = glob.glob(params['covid19RegionsDir_name'] + '/*regions-france-2020*.json')
    params['departments_name'] = config.get("inputs", "departments_name")
    params['citiesFile_name'] = config.get("inputs", "citiesFile_name")

    # contains all data of getFranceGeomDpts(params) or getFranceGeomRegs(params)
    # we can directly read it if it was saved so we do not compute it again
    params['franceGeom_name'] = config.get("inputs", "franceGeom_name")
    # location of saved file
    # params['graphFile'] = config.get("inputs", 'graphFile')
    # if not params['graphFile'].strip():
    #     params['graphFile'] = 'franceGeomRegions.txt'

    # load options infos
    # params["genAndQuit"] = distutils.util.strtobool(config.get("graph", "genAndQuit"))  # quit after graph is generated and saved
    # JD: it fails for me otherwise...
    params["genAndQuit"] = strtobool(config.get("graph", "genAndQuit"))  # quit after graph is generated and saved
    params["location"] = config.get("graph", "location")
    if params["location"] != "metropole":
        print('[ERROR]: location can only be metropole for now. Exit here.')
        exit(0)

    params["graph_type"] = config.get("graph", "graph_type")  # what type of graph is is

    # load model infos
    params["model_name"] = config.get("model", "model_name")
    if params["model_name"] == "SIRWellMixed":  # load parameters needed for this model
        params["k_SI"] = config.getfloat("model", "k_SI")
        params["k_IR"] = config.getfloat("model", "k_IR")
        params["S_ini"] = config.getint("model", "S_ini")
        params["I_ini"] = config.getint("model", "I_ini")
        params["R_ini"] = config.getint("model", "R_ini")
    elif params["model_name"] == "SIROnGraph":  # load parameters needed for this model
        # For now all parameters are the same wherever the region
        params["k_SI"] = config.getfloat("model", "k_SI")
        params["k_IR"] = config.getfloat("model", "k_IR")
        params["S_ini"] = config.getint("model", "S_ini")
        params["I_ini"] = config.getint("model", "I_ini")
        params["R_ini"] = config.getint("model", "R_ini")
    elif params["model_name"] == "SIRDWellMixed":  # load parameters needed for this model
        params["k_SI"] = config.getfloat("model", "k_SI")
        params["k_IR"] = config.getfloat("model", "k_IR")
        params["k_ID"] = config.getfloat("model", "k_ID")
        params["S_ini"] = config.getint("model", "S_ini")
        params["I_ini"] = config.getint("model", "I_ini")
        params["R_ini"] = config.getint("model", "R_ini")
        params["D_ini"] = config.getint("model", "D_ini")
    elif params["model_name"] == "SIRDOnGraph":  # load parameters needed for this model
        # For now all parameters are the same wherever the region
        params["k_SI"] = config.getfloat("model", "k_SI")
        params["k_IR"] = config.getfloat("model", "k_IR")
        params["k_ID"] = config.getfloat("model", "k_ID")
        params["S_ini"] = config.getint("model", "S_ini")
        params["I_ini"] = config.getint("model", "I_ini")
        params["R_ini"] = config.getint("model", "R_ini")
        params["D_ini"] = config.getint("model", "D_ini")
    elif params["model_name"] == "SIRerlangWellMixed":  # load parameters needed for this model
        params["k_SI"] = config.getfloat("model", "k_SI")
        params["k_IR"] = config.getfloat("model", "k_IR")
        params["S_ini"] = config.getint("model", "S_ini")
        params["I_ini"] = config.getint("model", "I_ini")
        params["R_ini"] = config.getint("model", "R_ini")
        params["nb_intermediate_states"] = config.getint("model", "nb_intermediate_states")
    elif params["model_name"] == "SIRerlangOnGraph":  # load parameters needed for this model
        # For now all parameters are the same wherever the region
        params["k_SI"] = config.getfloat("model", "k_SI")
        params["k_IR"] = config.getfloat("model", "k_IR")
        params["S_ini"] = config.getint("model", "S_ini")
        params["I_ini"] = config.getint("model", "I_ini")
        params["R_ini"] = config.getint("model", "R_ini")
        params["nb_intermediate_states"] = config.getint("model", "nb_intermediate_states")
    elif params["model_name"] == "RealSIRDOnGraph":  # load parameters needed for this model
        # For now all parameters are the same wherever the region: 
        # fix R_0 per region and addapt k_IS accordingly to population
        params["R_0"] = config.getfloat("model", "R_0")
        params["k_IR"] = config.getfloat("model", "k_IR")        
        params["k_ID"] = config.getfloat("model", "k_ID")
        params["time_start_day"] = config.get("model", "time_start_day")
    elif params["model_name"] == "RealSIRDWellMixed":  # load parameters needed for this model
        # For now all parameters are the same wherever the region: 
        # fix R_0 per region and addapt k_IS accordingly to population
        params["R_0"] = config.getfloat("model", "R_0")
        params["k_IR"] = config.getfloat("model", "k_IR")        
        params["k_ID"] = config.getfloat("model", "k_ID")
        params["time_start_day"] = config.get("model", "time_start_day")
    else:
        print("[ERROR]: model name not found. Exit here.")
        exit(0)        
    
    # load settigns infos
    params["time_start"] = config.getfloat("settings", "time_start")
    params["time_end"] = config.getfloat("settings", "time_end")
    
    # load fit infos
    # to know if we fit or not
    params["fit_model"] = strtobool(config.get("fit", "fit_model"))  
    dates_for_fit = config.get("fit", "dates_for_fit")        
    params["dates_for_fit"] = [p.replace(" ", "").replace('\n','') for p in
                              dates_for_fit.split(",")]
    method = config.get("fit", "method")
    
    #  load outputs infos
    params["outputs_dir"] = config.get("outputs", "outputs_dir")
    if not (
            os.path.exists(params["outputs_dir"])
    ):  # if directory does not exist, create it
        os.mkdir(params["outputs_dir"])
    plot_options = config.get("outputs", "plot_options")
    # remove spaces and return to line if there are different plot options
    params["plot_options"] = [p.replace(" ", "").replace('\n','') for p in
                              plot_options.split(",")]
    
    params["save_frequency"] = config.getint("outputs", "save_frequency")

    return params


def InitGraph(params):
    print("InitGraph")
    """
    initialize a graph of regions with positions, and ids per regions
    in: params, dictionnary with all needed parameters
    out: a graph of type networkx.Graph()
    """
    ### --- 1/ read data we need depending of what is written on params["graph_type"]
    franceGeom = getFranceGeom.readFranceGeom(params)
    G = nx.Graph()

    ### --- Construct graph
    # same routine for regions and departments as it is handled by readFranceGeom(params)
    if (params['graph_type'] == 'regions') or (params['graph_type'] == 'departments'):
        dic_total_population = {}
        dic_data_covid = {}
        dic_positions = {}
        dic_names = {}
        # add nodes to the graph
        for key in franceGeom.keys():
            id_node = franceGeom[key]["id"]  # get ID attributed from readFranceGeom(params)
            G.add_node(id_node)

            # add data on node_attributes, for now:
            # - Total population
            # - Confirmed cases
            # - Deaths
            # - Position
            dic_total_population[id_node] = franceGeom[key]['agesSex']['ftot'] + franceGeom[key]['agesSex']['htot']
            dic_data_covid[id_node] = franceGeom[key]['data_covid']            

            centroid = [franceGeom[key]['centroid'][0], franceGeom[key]['centroid'][1]]
            dic_positions[id_node] = centroid
            dic_names[id_node] = key
            id_node += 1

        # add attributes with positions to graph
        nx.set_node_attributes(G, name='positions', values=dic_positions)
        nx.set_node_attributes(G, name='names', values=dic_names)
        nx.set_node_attributes(G, name='total_pop', values=dic_total_population)
        nx.set_node_attributes(G, name='data_covid', values=dic_data_covid)

        # add edges to the graph
        for key in franceGeom.keys():
            id_node1 = franceGeom[key]["id"]  # get ID attributed from readFranceGeom(params)
            for keyneigh in franceGeom[key]['neighbours']:
                id_node2 = franceGeom[keyneigh]["id"]  # get ID attributed from readFranceGeom(params)
                edge12 = (id_node1, id_node2)
                edge21 = (id_node2, id_node1)
                if (not (edge12 in G.edges())) or (not (edge21 in G.edges())):
                    G.add_edge(*edge12)  # unpack edge tuple*
        # add edges based on airways connections
        for key in franceGeom.keys():
            id_node1 = franceGeom[key]["id"]  # get ID attributed from readFranceGeom(params)
            for connection in franceGeom[key]['airways']:
                keyneigh = connection[0]
                pax = connection[1]  # not used yet
                id_node2 = franceGeom[keyneigh]["id"]  # get ID attributed from readFranceGeom(params)
                # make sure to not add loop on itself
                if id_node1 != id_node2:
                    edge12 = (id_node1, id_node2)
                    edge21 = (id_node2, id_node1)
                    if (not (edge12 in G.edges())) or (not (edge21 in G.edges())):
                        G.add_edge(*edge12)  # unpack edge tuple*


    elif (params['graph_type'] == 'fullConnectedRegions') or (params['graph_type'] == 'fullConnectedDepartments'):

        dic_total_population = {}
        dic_data_covid = {}
        dic_positions = {}
        dic_names = {}
        # add nodes to the graph
        id_node = 0
        for key in franceGeom.keys():
            id_node = franceGeom[key]["id"]  # get ID attributed from readFranceGeom(params)
            G.add_node(id_node)

            # add data on node_attributes, for now:
            # - Total population
            # - Confirmed cases
            # - Deaths
            # - Position
            dic_total_population[id_node] = franceGeom[key]['agesSex']['ftot'] + franceGeom[key]['agesSex']['htot']
            dic_data_covid[id_node] = franceGeom[key]['data_covid']
            centroid = [franceGeom[key]['centroid'][0], franceGeom[key]['centroid'][1]]
            dic_positions[id_node] = centroid
            dic_names[id_node] = key
            id_node += 1
        # add attributes with positions to graph
        nx.set_node_attributes(G, name='positions', values=dic_positions)
        nx.set_node_attributes(G, name='names', values=dic_names)
        nx.set_node_attributes(G, name='total_pop', values=dic_total_population)
        nx.set_node_attributes(G, name='data_covid', values=dic_data_covid)

        # add edges of graph: all regions are connected to every other region
        for n1 in G.nodes():
            for n2 in G.nodes():
                if n1 != n2:
                    edge = (n1, n2)
                    G.add_edge(*edge)  # unpack edge tuple*

    print("num_edges = ", G.number_of_edges())
    print("num_nodes = ", G.number_of_nodes())

    print("nodes attributes = ", set([k for n in G.nodes for k in G.nodes[n].keys()]))
    print("edges attributes = ", set([k for n in G.edges for k in G.edges[n].keys()]))
    # ### --- plot graph
    # fig, ax = plt.subplots()    
    # dic_pos = nx.get_node_attributes(G, 'positions')
    # dic_labels = nx.get_node_attributes(G, 'names')            
    # nx.draw(G, with_labels=True, font_weight='bold', pos=dic_pos, labels=dic_labels)
    # plt.show()

    return G


def LaunchModel(params, graph):
    """
    Launch given model depending on name
    In: - params, dictionnary that contains all needed parameters
        - graph
    """
    print('LaunchModel')
    # initialize models
    params, graph = InitModels(params, graph)    
    if params["model_name"] == "SIRWellMixed":
        outputs = SolveOneCompartment(params)
        plot_SIRWellMixed(params, outputs)
    elif params["model_name"] == "SIRDWellMixed":
        outputs = SolveOneCompartment(params)
        plot_SIRDWellMixed(params, outputs)
    elif params["model_name"] == "SIRerlangWellMixed":
        outputs = SolveOneCompartment(params)
        plot_SIRerlangWellMixed(params, outputs)
    elif params["model_name"] == "SIROnGraph":
        outputs, graph = SolveOnGraph(params, graph)
        plot_SIROnGraph(params, outputs, graph)
    elif params["model_name"] == "SIRDOnGraph":
        outputs, graph = SolveOnGraph(params, graph)
        plot_SIRDOnGraph(params, outputs, graph)
    elif params["model_name"] == "SIRerlangOnGraph":
        outputs, graph = SolveOnGraph(params, graph)
        plot_SIRerlangOnGraph(params, outputs, graph)
    elif params["model_name"] == "RealSIRDOnGraph":
        outputs , graph = SolveOnGraph(params, graph)
        plot_TEST(params,outputs, graph)
    
    else:
        print("[ERROR]: invalid name in parser for model name. Exit here.\n")
        exit(0)


# =============================================================================
# LAUNCH WORKFLOW
# =============================================================================
params = InitParams(inputFile)
#  generate graph if file does not already exist
if not os.path.isfile(params['franceGeom_name']):
    print('geom files does not exist')
    params["franceFileExist"] = False
    G = InitGraph(params)
    if params['genAndQuit']:
        print("\ngraph has been generated - quitting now")
        print('file saved here :', params['franceGeom_name'])
        sys.exit(0)
else:
    if params['genAndQuit']:
        print("\ngraph already exist - nothing was done.")
        print('dumbass')
        sys.exit(0)
    else:
        params["franceFileExist"] = True
        G = InitGraph(params)

if (params['fit_model']):
    print("Start fitting procedure, Inch'Allah...")
    paramsToFit = InitFit()
    df = LaunchFit(paramsToFit, params, G)

LaunchModel(params, G)

t1 = time.time()
print("Code ran for: ", t1 - t0, " seconds")
