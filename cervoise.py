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
from numpy import * # in order to read numpy array rfom files

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
from plotData import *
from ODEs import *

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
    # if params["location"] != "metropole":
        # print('[ERROR]: location can only be metropole for now. Exit here.')
        # exit(0)

    params["graph_type"] = config.get("graph", "graph_type")  # what type of graph is is

    # load model infos
    params["model_name"] = config.get("model", "model_name")
    #__________________________________________________________________#
    #
    # --- WELL MIXED (= ONE COMPARTMENT) --- #
    #__________________________________________________________________#
    
    if params["model_name"] == "SIRWellMixed":  # load parameters needed for this model
        params["k_SI"] = config.getfloat("model", "k_SI")
        params["k_IR"] = config.getfloat("model", "k_IR")
        params["S_ini"] = config.getint("model", "S_ini")
        params["I_ini"] = config.getint("model", "I_ini")
        params["R_ini"] = config.getint("model", "R_ini")    
    elif params["model_name"] == "SIRHigherOrderWellMixed":  # load parameters needed for this model
        params["k_SI"] = config.getfloat("model", "k_SI")
        params["k_IR"] = config.getfloat("model", "k_IR")
        params["n_order"] = config.getint("model", "n_order")        
        params["S_ini"] = config.getint("model", "S_ini")
        params["I_ini"] = config.getint("model", "I_ini")
        params["R_ini"] = config.getint("model", "R_ini")
    elif params["model_name"] == "SIRHigherOrderSeveralReactionsWellMixed":  # load parameters needed for this model
        params["tilde_k_SI"] = config.getfloat("model", "tilde_k_SI")
        params["k_IR"] = config.getfloat("model", "k_IR")
        params["n_order"] = config.getint("model", "n_order")
        # get all probabilities of all reactions
        for i in range(0,params["n_order"]+1):
            str_to_search = 'p_'+str(i)
            params[str_to_search] = config.getfloat("model", str_to_search)
        params["S_ini"] = config.getint("model", "S_ini")
        params["I_ini"] = config.getint("model", "I_ini")
        params["R_ini"] = config.getint("model", "R_ini")            
    elif params["model_name"] == "SIHigherOrderSeveralReactionsWellMixed":  # load parameters needed for this model
        params["tilde_k_SI"] = config.getfloat("model", "tilde_k_SI")        
        params["n_order"] = config.getint("model", "n_order")
        # get all probabilities of all reactions
        for i in range(0,params["n_order"]+1):
            str_to_search = 'p_'+str(i)
            params[str_to_search] = config.getfloat("model", str_to_search)
        params["S_ini"] = config.getint("model", "S_ini")
        params["I_ini"] = config.getint("model", "I_ini")        
    elif params["model_name"] == "SIRDWellMixed":  # load parameters needed for this model
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
    elif params["model_name"] == "SIRerlang1WellMixed":  # load parameters needed for this model
        params["k_SI"] = config.getfloat("model", "k_SI")
        params["k_IR"] = config.getfloat("model", "k_IR")
        params["S_ini"] = config.getint("model", "S_ini")
        params["I_ini"] = config.getint("model", "I_ini")
        params["R_ini"] = config.getint("model", "R_ini")
        params["nb_intermediate_states"] = config.getint("model", "nb_intermediate_states")        
    elif params["model_name"] == "SIRerlang2WellMixed":  # load parameters needed for this model
        params["k_SI"] = config.getfloat("model", "k_SI")
        params["k_IR"] = config.getfloat("model", "k_IR")
        params["S_ini"] = config.getint("model", "S_ini")
        params["I_ini"] = config.getint("model", "I_ini")
        params["R_ini"] = config.getint("model", "R_ini")
        params["nb_intermediate_states"] = config.getint("model", "nb_intermediate_states")
    elif params["model_name"] == "SIWellMixed":  # load parameters needed for this model
        params["k_SI"] = config.getfloat("model", "k_SI")        
        params["S_ini"] = config.getint("model", "S_ini")
        params["I_ini"] = config.getint("model", "I_ini")
    elif params["model_name"] == "SIerlangWellMixed":  # load parameters needed for this model
        params["k_SI"] = config.getfloat("model", "k_SI")        
        params["S_ini"] = config.getint("model", "S_ini")
        params["I_ini"] = config.getint("model", "I_ini")
        params["nb_intermediate_states"] = config.getint("model", "nb_intermediate_states")
    elif params["model_name"] == "SIerlangWellMixedModified":  # load parameters needed for this model
        params["k_SI"] = config.getfloat("model", "k_SI")        
        params["S_ini"] = config.getint("model", "S_ini")
        params["I_ini"] = config.getint("model", "I_ini")
        params["nb_intermediate_states"] = config.getint("model", "nb_intermediate_states")
    elif params["model_name"] == "SIerlangOtherWellMixed":  # load parameters needed for this model
        params["k_SI"] = config.getfloat("model", "k_SI")        
        params["S_ini"] = config.getint("model", "S_ini")
        params["I_ini"] = config.getint("model", "I_ini")
        params["nb_intermediate_states"] = config.getint("model", "nb_intermediate_states")
    elif params["model_name"] == "SIAndIntermediateErlangWellMixed":  # load parameters needed for this model
        params["k_SL"] = config.getfloat("model", "k_SL")
        params["k_LI"] = config.getfloat("model", "k_LI")
        params["S_ini"] = config.getint("model", "S_ini")
        params["I_ini"] = config.getint("model", "I_ini")
        params["L_ini"] = config.getint("model", "L_ini")
        params["nb_intermediate_states"] = config.getint("model", "nb_intermediate_states")
    elif params["model_name"] == "SIAndIntermediateErlangWellMixedModified":  # load parameters needed for this model
        params["k_SL"] = config.getfloat("model", "k_SL")
        params["k_LI"] = config.getfloat("model", "k_LI")
        params["S_ini"] = config.getint("model", "S_ini")
        params["I_ini"] = config.getint("model", "I_ini")
        params["L_ini"] = config.getint("model", "L_ini")
        params["nb_intermediate_states"] = config.getint("model", "nb_intermediate_states")
    elif params["model_name"] == "SIerlangFixSWellMixed":  # load parameters needed for this model    
        params["k_SI"] = config.getfloat("model", "k_SI")        
        params["S_ini"] = config.getint("model", "S_ini")
        params["I_ini"] = config.getint("model", "I_ini")
        params["nb_intermediate_states"] = config.getint("model", "nb_intermediate_states")
    elif params["model_name"] == "SILinearErlangWellMixed":  # load parameters needed for this model
        params["k_SI"] = config.getfloat("model", "k_SI")        
        params["S_ini"] = config.getint("model", "S_ini")
        params["I_ini"] = config.getint("model", "I_ini")    
        params["nb_intermediate_states"] = config.getint("model", "nb_intermediate_states")        
    elif params["model_name"] == "SIRGroupMeetingsWellMixed":  # load parameters needed for this model
        params["k_SI"] = config.getfloat("model", "k_SI")
        params["k_IR"] = config.getfloat("model", "k_IR")
        params["n_max"] = config.getint("model", "n_max")
        params["S_ini"] = config.getint("model", "S_ini")
        params["I_ini"] = config.getint("model", "I_ini")
        params["R_ini"] = config.getint("model", "R_ini")            
    
    # --- deals with real data
    elif params["model_name"] == "RealSIRDWellMixed":  # load parameters needed for this model        
        # For now all parameters are the same wherever the region: 
        # fix R_0 per region and addapt k_IS accordingly to population
        params["R_0"] = config.getfloat("model", "R_0")
        params["k_IR"] = config.getfloat("model", "k_IR")        
        params["k_ID"] = config.getfloat("model", "k_ID")
        params["time_start_day"] = config.get("model", "time_start_day")
    elif params["model_name"] == "RealSIRLDWellMixed":  # load parameters needed for this model        
        params["k_SI"]  = config.getfloat("model", "k_SI")
        params["k_IR"]  = config.getfloat("model", "k_IR")        
        params["k_IL"]  = config.getfloat("model", "k_IL")
        params["k_LD"]  = config.getfloat("model", "k_LD")
        params["alpha"] = config.getfloat("model", "alpha")        
        params["time_start_day"] = config.get("model", "time_start_day")
        params["I_ini"] = config.getint("model", "I_ini")
        params["R_ini"] = config.getint("model", "R_ini")
        params["L_ini"] = config.getint("model", "L_ini")
    elif params["model_name"] == "RealSIRLDmodifiedWellMixed":  # load parameters needed for this model        
        params["k_SI"]  = config.getfloat("model", "k_SI")
        params["k_IR"]  = config.getfloat("model", "k_IR")        
        params["k_IL"]  = config.getfloat("model", "k_IL")
        params["k_LR"] = config.getfloat("model", "k_LR")
        params["k_LD"]  = config.getfloat("model", "k_LD")
        params["alpha"] = config.getfloat("model", "alpha")        
        params["time_start_day"] = config.get("model", "time_start_day")
        params["I_ini"] = config.getint("model", "I_ini")
        params["R_ini"] = config.getint("model", "R_ini")
        params["L_ini"] = config.getint("model", "L_ini")
    elif params["model_name"] == "RealSEIRLDWellMixed":  # load parameters needed for this model        
        params["time_start_day"] = config.get("model", "time_start_day")
        params["k_SE"] = config.getfloat("model", "k_SE")
        params["k_EI"] = config.getfloat("model", "k_EI")
        params["k_IR"] = config.getfloat("model", "k_IR")
        params["k_IL"] = config.getfloat("model", "k_IL")
        params["k_LR"] = config.getfloat("model", "k_LR")
        params["k_LD"] = config.getfloat("model", "k_LD")
        # fraction that relates confirmed and real infected
        # I_confirmed = alpha * I
        params["alpha"] = config.getfloat("model", "alpha")
        # initial number of population
        params["E_ini"] = config.getint("model", "E_ini")
        params["I_ini"] = config.getint("model", "I_ini")
        params["R_ini"] = config.getint("model", "R_ini")
        params["L_ini"] = config.getint("model", "L_ini")
    elif params["model_name"] == "RealSEIRLDWellMixed":  # load parameters needed for this model        
        params["time_start_day"] = config.get("model", "time_start_day")
        params["k_SE"] = config.getfloat("model", "k_SE")
        params["k_EI"] = config.getfloat("model", "k_EI")
        params["k_IR"] = config.getfloat("model", "k_IR")
        params["k_IL"] = config.getfloat("model", "k_IL")
        params["k_LR"] = config.getfloat("model", "k_LR")
        params["k_LD"] = config.getfloat("model", "k_LD")
        # fraction that relates confirmed and real infected
        # I_confirmed = alpha * I
        params["alpha"] = config.getfloat("model", "alpha")
        # initial number of population
        params["E_ini"] = config.getint("model", "E_ini")
        params["I_ini"] = config.getint("model", "I_ini")
        params["R_ini"] = config.getint("model", "R_ini")
        params["L_ini"] = config.getint("model", "L_ini")
    elif params["model_name"] == "RealSEIsIasHoRadRDWellMixed":  # load parameters needed for this model        
        params["time_start_day"] = config.get("model", "time_start_day")       
        params["k_SE"]    = config.getfloat("model", "k_SE")
        params["k_EIs"]   = config.getfloat("model", "k_EIs") 
        params["k_EIas"]  = config.getfloat("model", "k_EIas") 
        params["k_IsR"]   = config.getfloat("model", "k_IsR") 
        params["k_IsHo"]  = config.getfloat("model", "k_IsHo") 
        params["k_IasR"]  = config.getfloat("model", "k_IasR") 
        params["k_HoRad"] = config.getfloat("model", "k_HoRad") 
        params["k_HoD"]   = config.getfloat("model", "k_HoD") 
        # initial number of population
        # params["S_ini"]   = config.getint("model", "S_ini")
        params["E_ini"]   = config.getint("model", "E_ini")
        params["Is_ini"]  = config.getint("model", "Is_ini")
        params["Ias_ini"] = config.getint("model", "Ias_ini")
        # params["Ho_ini"]  = config.getint("model", "Ho_ini")
        params["R_ini"]   = config.getint("model", "R_ini")            
    
    # --- agent based tracking of individuals
    elif params["model_name"] == "SIRWellMixedAgentBased":  # load parameters needed for this model
        params["k_SI"] = config.getfloat("model", "k_SI")
        params["k_IR"] = config.getfloat("model", "k_IR")
        params["S_ini"] = config.getint("model", "S_ini")
        params["I_ini"] = config.getint("model", "I_ini")
        params["R_ini"] = config.getint("model", "R_ini")
    elif params["model_name"] == "SIerlangWellMixedAgentBased":  # load parameters needed for this model
        params["k_SI"] = config.getfloat("model", "k_SI")        
        params["S_ini"] = config.getint("model", "S_ini")
        params["I_ini"] = config.getint("model", "I_ini")
        params["nb_intermediate_states"] = config.getint("model", "nb_intermediate_states")
    elif params["model_name"] == "SIAndIntermediateErlangWellMixedAgentBased":  # load parameters needed for this model
        params["k_SL"] = config.getfloat("model", "k_SL")
        params["k_LI"] = config.getfloat("model", "k_LI")
        params["S_ini"] = config.getint("model", "S_ini")
        params["L_ini"] = config.getint("model", "L_ini")
        params["I_ini"] = config.getint("model", "I_ini")
        params["nb_intermediate_states"] = config.getint("model", "nb_intermediate_states")    
    elif params["model_name"] == "SIerlangModifiedWellMixedAgentBased":  # load parameters needed for this model        
        params["k_SI"] = config.getfloat("model", "k_SI")        
        params["S_ini"] = config.getint("model", "S_ini")
        params["I_ini"] = config.getint("model", "I_ini")
        params["nb_intermediate_states"] = config.getint("model", "nb_intermediate_states")
    elif params["model_name"] == "SIRerlangKrylovaWellMixedAgentBased":  # load parameters needed for this model
        params["k_SI"] = config.getfloat("model", "k_SI")        
        params["k_IR"] = config.getfloat("model", "k_IR")        
        params["S_ini"] = config.getint("model", "S_ini")
        params["I_ini"] = config.getint("model", "I_ini")
        params["R_ini"] = config.getint("model", "R_ini")
        params["nb_intermediate_states"] = config.getint("model", "nb_intermediate_states")

    #__________________________________________________________________#
    #
    # --- ON GRAPH = SEVERAL COMPARTMENTS --- #
    #__________________________________________________________________#
    
    elif params["model_name"] == "SIROnGraph":  # load parameters needed for this model
        # For now all parameters are the same wherever the region
        params["k_SI"] = config.getfloat("model", "k_SI")
        params["k_IR"] = config.getfloat("model", "k_IR")
        params["S_ini"] = config.getint("model", "S_ini")
        params["I_ini"] = config.getint("model", "I_ini")
        params["R_ini"] = config.getint("model", "R_ini")
    elif params["model_name"] == "SIRDOnGraph":  # load parameters needed for this model
        # For now all parameters are the same wherever the region
        params["k_SI"] = config.getfloat("model", "k_SI")
        params["k_IR"] = config.getfloat("model", "k_IR")
        params["k_ID"] = config.getfloat("model", "k_ID")
        params["S_ini"] = config.getint("model", "S_ini")
        params["I_ini"] = config.getint("model", "I_ini")
        params["R_ini"] = config.getint("model", "R_ini")
        params["D_ini"] = config.getint("model", "D_ini")
    elif params["model_name"] == "SIRerlangOnGraph":  # load parameters needed for this model
        # For now all parameters are the same wherever the region
        params["k_SI"] = config.getfloat("model", "k_SI")
        params["k_IR"] = config.getfloat("model", "k_IR")
        params["S_ini"] = config.getint("model", "S_ini")
        params["I_ini"] = config.getint("model", "I_ini")
        params["R_ini"] = config.getint("model", "R_ini")
        params["nb_intermediate_states"] = config.getint("model", "nb_intermediate_states")    
    elif params["model_name"] == "TestScalingHoppingRates":  # load parameters needed for this model
        ### --- MARWANE --- ###
        print('coucou')
        ### ...
        ### --------------- ###
    # --- deals with real data
    
    # --- deals with real data
    elif params["model_name"] == "RealSIRDOnGraph":  # load parameters needed for this model
        # For now all parameters are the same wherever the region: 
        # fix R_0 per region and addapt k_IS accordingly to population
        params["R_0"] = config.getfloat("model", "R_0")
        params["k_IR"] = config.getfloat("model", "k_IR")        
        params["k_ID"] = config.getfloat("model", "k_ID")
        params["time_start_day"] = config.get("model", "time_start_day")    
    elif params["model_name"] == "RealSEIRLDOnGraphInitTEST":  # load parameters needed for this model
        params["time_start_day"] = config.get("model", "time_start_day")
        params["k_SE"] = config.getfloat("model", "k_SE")
        params["k_EI"] = config.getfloat("model", "k_EI")
        params["k_IR"] = config.getfloat("model", "k_IR")
        params["k_IL"] = config.getfloat("model", "k_IL")
        params["k_LR"] = config.getfloat("model", "k_LR")
        params["k_LD"] = config.getfloat("model", "k_LD")        
        # fraction that relates confirmed and real infected
        # I_confirmed = alpha * I
        params["alpha"] = config.getfloat("model", "alpha")
        # initial number of population
        params["E_ini"] = config.getint("model", "E_ini")
        params["I_ini"] = config.getint("model", "I_ini")
        params["L_ini"] = config.getint("model", "L_ini")
    
    # --- agent based tracking of individuals
    elif params["model_name"] == "SIROnGraphAgentBased":  # load parameters needed for this model
        # For now all parameters are the same wherever the region
        params["k_SI"] = config.getfloat("model", "k_SI")
        params["k_IR"] = config.getfloat("model", "k_IR")
        params["S_ini"] = config.getint("model", "S_ini")
        params["I_ini"] = config.getint("model", "I_ini")
        params["R_ini"] = config.getint("model", "R_ini")
    elif params["model_name"] == "SEIIssROnGraphAgentBased":  # load parameters needed for this model
        # For now all parameters are the same wherever the region
        params["k_SE_I"] = config.getfloat("model", "k_SE_I")
        params["k_SE_Iss"] = config.getfloat("model", "k_SE_Iss")
        params["k_EI"]   = eval(config.get("model", "k_EI"))
        params["k_EIss"] = eval(config.get("model", "k_EIss"))                
        params["k_IR"]   = eval(config.get("model", "k_IR"))
        params["k_IssR"] = eval(config.get("model", "k_IssR"))        
        params["S_ini"] = config.getint("model", "S_ini")
        params["E_ini"] = config.getint("model", "E_ini")
        params["I_ini"] = config.getint("model", "I_ini")
        params["Iss_ini"] = config.getint("model", "Iss_ini")
        params["R_ini"] = config.getint("model", "R_ini")
    elif params["model_name"] == "TestSimpleReturnHomeModel":  # load parameters needed for this model
        ### --- MARWANE --- ###
        print('coucou')
        ### ...
        ### --------------- ###
    else:
        print("[ERROR]: model name: ",params['model_name']," not found. Exit here.")
        exit(0)        
    
    # load settigns infos
    params["time_start"] = config.getfloat("settings", "time_start")
    params["time_end"] = config.getfloat("settings", "time_end")
    
    if config.has_option("settings", "seed_num"):
        params["seed_num"] = config.getint("settings", "seed_num")
    else: # if seed is not given then default is 0
        params["seed_num"] = 0
        
    print(params["seed_num"])
    # load fit infos
    # to know if we fit or not
    params["fit_model"] = strtobool(config.get("fit", "fit_model"))
    params["fit_with_ODEs"] = strtobool(config.get("fit", "with_ODEs"))
    dates_for_fit = config.get("fit", "dates_for_fit")        
    params["dates_for_fit"] = [p.replace(" ", "").replace('\n','') for p in
                              dates_for_fit.split(",")]
    params['method'] = config.get("fit", "method")
    
    #  load outputs infos
    params["outputs_dir"] = config.get("outputs", "outputs_dir")
    if not (
            os.path.exists(params["outputs_dir"])
    ):  # if directory does not exist, create it
        params["outputs_dir"]
        os.makedirs(params["outputs_dir"])
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
    params['total_pop'] = 0 # feed dictionnary with aggregated data from graph
    params['total_covid_data'] = {} # feed dictionnary with aggregated data from graph
    
    list_outremer = ['Mayotte', 'La Réunion', 'Guyane', 'Guadeloupe', 'Martinique']
    
    ### --- Construct graph
    # same routine for regions and departments as it is handled by readFranceGeom(params)
    if (params['graph_type'] == 'regions') or (params['graph_type'] == 'departments'):
        dic_total_population = {}
        dic_data_covid = {}
        dic_data_hosp = {}
        dic_positions = {}
        dic_names = {}
        dic_areas = {}        
        # add nodes to the graph
        for key in franceGeom.keys():
            id_node = franceGeom[key]["id"]  # get ID attributed from readFranceGeom(params)
            
            if (params['location'] == 'metropole'):
                # only add departments/ergions which are in metropole
                if not(key in list_outremer):
                    G.add_node(id_node)                    
                    # add data on node_attributes, for now:
                    # - Total population
                    # - Confirmed cases (old source)
                    # - Deaths (old source)
                    # - Recovered (old source)
                    # - Hospital data
                    # - Position 
                    # - Area
                    dic_total_population[id_node] = franceGeom[key]['agesSex']['ftot'] + franceGeom[key]['agesSex']['htot']
                    dic_data_covid[id_node] = franceGeom[key]['data_covid']
                    dic_data_hosp[id_node] = franceGeom[key]['data_hosp']
                    centroid = [franceGeom[key]['centroid'][0], franceGeom[key]['centroid'][1]]
                    dic_positions[id_node] = centroid                    
                    dic_names[id_node] = key
                    dic_areas[id_node] = franceGeom[key]['area']                    
                    # print("id_node, dic_names[id_node] = " , id_node, dic_names[id_node])
                    params['total_pop'] += dic_total_population[id_node]
                    for _date in dic_data_covid[id_node].keys():
                        if not(_date in params['total_covid_data']):
                            params['total_covid_data'][_date] = {}
                            params['total_covid_data'][_date]['Confirmed'] = 0
                            params['total_covid_data'][_date]['Deaths']    = 0
                            params['total_covid_data'][_date]['Recovered'] = 0
                        params['total_covid_data'][_date]['Confirmed'] += dic_data_covid[id_node][_date]['Confirmed']
                        params['total_covid_data'][_date]['Deaths']    += dic_data_covid[id_node][_date]['Deaths']
                        params['total_covid_data'][_date]['Recovered'] += dic_data_covid[id_node][_date]['Recovered']
            else:
                G.add_node(id_node)
                # add data on node_attributes, for now:                
                # - Total population
                # - Confirmed cases
                # - Deaths
                # - Position
                # - Area
                dic_total_population[id_node] = franceGeom[key]['agesSex']['ftot'] + franceGeom[key]['agesSex']['htot']
                dic_data_covid[id_node] = franceGeom[key]['data_covid']
                dic_data_hosp[id_node] = franceGeom[key]['data_hosp']                

                centroid = [franceGeom[key]['centroid'][0], franceGeom[key]['centroid'][1]]
                dic_positions[id_node] = centroid
                dic_names[id_node] = key
                dic_areas[id_node] = franceGeom[key]['area']
                
                # print("id_node, dic_names[id_node] = " , id_node, dic_names[id_node])
                params['total_pop'] += dic_total_population[id_node]
                for _date in dic_data_covid[id_node].keys():
                    if not(_date in params['total_covid_data']):
                        params['total_covid_data'][_date] = {}
                        params['total_covid_data'][_date]['Confirmed'] = 0
                        params['total_covid_data'][_date]['Deaths']    = 0
                        params['total_covid_data'][_date]['Recovered'] = 0
                    params['total_covid_data'][_date]['Confirmed'] += dic_data_covid[id_node][_date]['Confirmed']
                    params['total_covid_data'][_date]['Deaths']    += dic_data_covid[id_node][_date]['Deaths']
                    params['total_covid_data'][_date]['Recovered'] += dic_data_covid[id_node][_date]['Recovered']

        # add attributes with positions to graph
        nx.set_node_attributes(G, name='positions', values=dic_positions)
        nx.set_node_attributes(G, name='names', values=dic_names)
        nx.set_node_attributes(G, name='total_pop', values=dic_total_population)
        nx.set_node_attributes(G, name='data_covid', values=dic_data_covid)
        nx.set_node_attributes(G, name='data_hosp', values=dic_data_hosp)
        nx.set_node_attributes(G, name='areas', values=dic_areas)        
        
        if (params['location'] == 'metropole'):
            # add edges to the graph
            for key in franceGeom.keys():
                if not(key in list_outremer):
                    id_node1 = franceGeom[key]["id"]  # get ID attributed from readFranceGeom(params)
                    for keyneigh in franceGeom[key]['neighbours']:
                        if not(keyneigh in list_outremer):
                            id_node2 = franceGeom[keyneigh]["id"]  # get ID attributed from readFranceGeom(params)
                            edge12 = (id_node1, id_node2)
                            edge21 = (id_node2, id_node1)
                            if (not (edge12 in G.edges())) or (not (edge21 in G.edges())):
                                G.add_edge(*edge12)  # unpack edge tuple*
            # add edges based on airways connections
            for key in franceGeom.keys():
                if not(key in list_outremer):
                    id_node1 = franceGeom[key]["id"]  # get ID attributed from readFranceGeom(params)
                    for connection in franceGeom[key]['airways']:
                        keyneigh = connection[0]
                        if not(keyneigh in list_outremer):
                            pax = connection[1]  # not used yet
                            id_node2 = franceGeom[keyneigh]["id"]  # get ID attributed from readFranceGeom(params)
                            # make sure to not add loop on itself
                            if id_node1 != id_node2:
                                edge12 = (id_node1, id_node2)
                                edge21 = (id_node2, id_node1)
                                if (not (edge12 in G.edges())) or (not (edge21 in G.edges())):
                                    G.add_edge(*edge12)  # unpack edge tuple*
        
        else:
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
        dic_data_hosp = {}
        dic_positions = {}
        dic_names = {}
        dic_areas = {}        
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
            # - Area
            dic_total_population[id_node] = franceGeom[key]['agesSex']['ftot'] + franceGeom[key]['agesSex']['htot']
            dic_data_covid[id_node] = franceGeom[key]['data_covid']
            dic_data_hosp[id_node] = franceGeom[key]['data_hosp']
            centroid = [franceGeom[key]['centroid'][0], franceGeom[key]['centroid'][1]]
            dic_positions[id_node] = centroid
            dic_names[id_node] = key
            dic_areas[id_node] = franceGeom[key]['area']
            
            params['total_pop'] += dic_total_population[id_node]
            for _date in dic_data_covid[id_node].keys():
                if not(_date in params['total_covid_data']):
                    params['total_covid_data'][_date] = {}
                    params['total_covid_data'][_date]['Confirmed'] = 0
                    params['total_covid_data'][_date]['Deaths']    = 0
                    params['total_covid_data'][_date]['Recovered'] = 0
                params['total_covid_data'][_date]['Confirmed'] += dic_data_covid[id_node][_date]['Confirmed']
                params['total_covid_data'][_date]['Deaths']    += dic_data_covid[id_node][_date]['Deaths']
                params['total_covid_data'][_date]['Recovered'] += dic_data_covid[id_node][_date]['Recovered']        
            id_node += 1
        # add attributes with positions to graph
        nx.set_node_attributes(G, name='positions', values=dic_positions)
        nx.set_node_attributes(G, name='names', values=dic_names)
        nx.set_node_attributes(G, name='total_pop', values=dic_total_population)
        nx.set_node_attributes(G, name='data_covid', values=dic_data_covid)
        nx.set_node_attributes(G, name='data_hosp', values=dic_data_hosp)
        nx.set_node_attributes(G, name='areas', values=dic_areas)

        # add edges of graph: all regions are connected to every other region
        for n1 in G.nodes():
            for n2 in G.nodes():
                if n1 != n2:
                    edge = (n1, n2)
                    G.add_edge(*edge)  # unpack edge tuple*


    # add  densities as  node_attribute so we don't have to compute it 
    # every time we need it
    dic_densities = {}
    dic_areas = nx.get_node_attributes(G, 'areas')    
    dic_total_pop = nx.get_node_attributes(G, 'total_pop')   
    for id_node in G.nodes():
        dic_densities[id_node] = dic_total_population[id_node] / dic_areas[id_node]
        
    nx.set_node_attributes(G, name='densities', values=dic_densities)
    # print("num_edges = ", G.number_of_edges())
    # print("num_nodes = ", G.number_of_nodes())
    # print("nodes attributes = ", set([k for n in G.nodes for k in G.nodes[n].keys()]))
    # print("edges attributes = ", set([k for n in G.edges for k in G.edges[n].keys()]))
    # dic_areas = nx.get_node_attributes(G, 'areas')    
    # dic_total_pop = nx.get_node_attributes(G, 'total_pop')    
    # total_area = 0.
    # total_pop = 0
    # for k in dic_areas.keys():
        # total_area += dic_areas[k]
        # total_pop += dic_total_pop[k]
        # print(k, dic_areas[k])
        
    # print('\ntotal_pop, total_area, density = ' , total_pop, total_area, total_pop / total_area)
    # exit(0)
    # ### --- plot graph
    # fig, ax = plt.subplots()
    # dic_pos = nx.get_node_attributes(G, 'positions')
    # dic_labels = nx.get_node_attributes(G, 'names')    
    # nx.draw_networkx_edges(G, font_weight='bold', pos=dic_pos, alpha = 0.5 )
    # nx.draw_networkx_nodes(G, font_weight='bold', pos=dic_pos, alpha = 0.5, node_size = 500)
    # nx.draw_networkx_labels(G, font_weight='bold', pos=dic_pos, labels=dic_labels, font_size=15)
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
    elif params["model_name"] == "SIRHigherOrderWellMixed":
        outputs = SolveOneCompartment(params)
        plot_SIRHigherOrderWellMixed(params, outputs)
    elif params["model_name"] == "SIRHigherOrderSeveralReactionsWellMixed":
        outputs = SolveOneCompartment(params)
        plot_SIRHigherOrderSeveralReactionsWellMixed(params, outputs)
    elif params["model_name"] == "SIHigherOrderSeveralReactionsWellMixed":
        outputs = SolveOneCompartment(params)
        plot_SIHigherOrderSeveralReactionsWellMixed(params, outputs)
    elif params["model_name"] == "SIRDWellMixed":
        outputs = SolveOneCompartment(params)
        plot_SIRDWellMixed(params, outputs)
    elif params["model_name"] == "SIRerlangWellMixed":
        outputs = SolveOneCompartment(params)
        plot_SIRerlangWellMixed(params, outputs)
    elif params["model_name"] == "SIRerlang1WellMixed":
        outputs = SolveOneCompartment(params)
        plot_SIRerlang1WellMixed(params, outputs)
    elif params["model_name"] == "SIRerlang2WellMixed":
        outputs = SolveOneCompartment(params)
        plot_SIRerlang2WellMixed(params, outputs)
    elif params["model_name"] == "SIWellMixed":
        outputs = SolveOneCompartment(params)
        plot_SIWellMixed(params, outputs)
    elif params["model_name"] == "SILinearErlangWellMixed":
        outputs = SolveOneCompartment(params)
        plot_SILinearErlangWellMixed(params, outputs)        
    elif params["model_name"] == "SIerlangWellMixed":
        outputs = SolveOneCompartment(params)
        plot_SIerlangWellMixed(params, outputs)
    elif params["model_name"] == "SIerlangWellMixedModified":
        outputs, waiting_times = SolveOneCompartment(params)
        plot_SIerlangWellMixedModified(params, outputs, waiting_times)
    elif params["model_name"] == "SIerlangOtherWellMixed":
        outputs = SolveOneCompartment(params)
        plot_SIerlangOtherWellMixed(params, outputs)
    elif params["model_name"] == "SIerlangFixSWellMixed":
        outputs = SolveOneCompartment(params)
        plot_SIerlangFixSWellMixed(params, outputs)
    elif params["model_name"] == "SIAndIntermediateErlangWellMixed":
        outputs = SolveOneCompartment(params)
        plot_SIAndIntermediateErlangWellMixed(params, outputs)     
    elif params["model_name"] == "SIAndIntermediateErlangWellMixedModified":
        outputs, waiting_times = SolveOneCompartment(params)
        plot_SIAndIntermediateErlangWellMixedModified(params, outputs, waiting_times)
    elif params["model_name"] == "SIRGroupMeetingsWellMixed":
        outputs = SolveOneCompartment(params)
        plot_SIRGroupMeetingsWellMixed(params, outputs)
    elif params["model_name"] == "RealSEIRLDWellMixed":
        outputs = SolveOneCompartment(params)
        dates_to_plot = []
        # copy list otherwize references are copied... and it messes up
        # with params['dates_for_fit'] list later on (e.g. on the fit function)
        dates_to_plot = params['dates_for_fit'].copy()
        dates_to_plot.insert(0, params['time_start_day'])
        plotDataAndSEIRLDResults(params, graph,outputs, dates_to_plot)    
    elif params["model_name"] == "RealSEIsIasHoRadRDWellMixed":                
        outputs = SolveOneCompartment(params)
        dates_to_plot = []
        # copy list otherwize references are copied... and it messes up
        # with params['dates_for_fit'] list later on (e.g. on the fit function)
        dates_to_plot = params['dates_for_fit'].copy()
        dates_to_plot.insert(0, params['time_start_day'])
        
        plotDataAndRealSEIsIasHoRadRDWellMixed(params, graph,outputs, dates_to_plot)            
    elif params["model_name"] == "SIRWellMixedAgentBased":
        outputs = SolveAgentBasedOneCompartment(params)
        plot_SIRWellMixedAgentBased(params, outputs)
    elif params["model_name"] == "SIerlangWellMixedAgentBased":
        outputs = SolveAgentBasedOneCompartment(params)
        plot_SIerlangWellMixedAgentBased(params, outputs)
    elif params["model_name"] == "SIAndIntermediateErlangWellMixedAgentBased":
        outputs = SolveAgentBasedOneCompartment(params)
        plot_SIAndIntermediateErlangWellMixedAgentBased(params, outputs)
    elif params["model_name"] == "SIerlangModifiedWellMixedAgentBased":
        outputs, waiting_times = SolveAgentBasedOneCompartment(params)
        plot_SIerlangModifiedWellMixedAgentBased(params, outputs, waiting_times)
    elif params["model_name"] == "SIRerlangKrylovaWellMixedAgentBased":
        outputs = SolveAgentBasedOneCompartment(params)
        plot_SIRerlangKrylovaWellMixedAgentBased(params, outputs)
    elif params["model_name"] == "SIROnGraph":
        outputs, graph = SolveOnGraph(params, graph)
        plot_SIROnGraph(params, outputs, graph)
    elif params["model_name"] == "SIRDOnGraph":
        outputs, graph = SolveOnGraph(params, graph)
        plot_SIRDOnGraph(params, outputs, graph)
    elif params["model_name"] == "SIRerlangOnGraph":
        outputs, graph = SolveOnGraph(params, graph)
        plot_SIRerlangOnGraph(params, outputs, graph)
    ### --- MARWANE START --- ###
    elif params["model_name"] == "TestScalingHoppingRates":
        outputs, graph = SolveOnGraph(params, graph)
        plot_TestScalingHoppingRates(params, outputs, graph)   
    ### --- MARWANE END   --- ###  
    elif params["model_name"] == "RealSIRDOnGraph":                 
        outputs , graph = SolveOnGraph(params, graph)
        plot_TEST(params,outputs, graph)
    elif params["model_name"] == "RealSEIRLDOnGraphInitTEST":
        outputs, G = SolveOnGraph(params, graph)
        dates_to_plot = []
        # copy list otherwize references are copied... and it messes up
        # with params['dates_for_fit'] list later on (e.g. on the fit function)
        dates_to_plot = params['dates_for_fit'].copy()
        dates_to_plot.insert(0, params['time_start_day'])        
        plotDataAndSEIRLDOnGraphResults(params, G, outputs, dates_to_plot)
    elif params["model_name"] == "SIROnGraphAgentBased":        
        outputs, graph = SolveAgentBasedOnGraph(params, graph)
        plot_SIROnGraphAgentBased(params, outputs, graph)
    elif params["model_name"] == "SEIIssROnGraphAgentBased":      
        outputs, graph = SolveAgentBasedOnGraph(params, graph)
        plot_SEIIssROnGraphAgentBased(params, outputs, graph)
    ### --- MARWANE START --- ###
    elif params["model_name"] == "TestSimpleReturnHomeModelAgentBased":
        outputs, graph = SolveAgentBasedOnGraph(params, graph)
        plot_TestSimpleReturnHomeModelAgentBasedl(params, outputs, graph)
    ### --- MARWANE END   --- ###  
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

# plot france data 
#TODO: kind of useless now or it should be adapted later
# to deal with new data
if ("PlotData" in params['plot_options']):
    dates_to_plot = []
    # copy list otherwize references are copied... and it messes up
    # with params['dates_for_fit'] list later on (e.g. on the fit function)
    dates_to_plot = params['dates_for_fit'].copy()
    dates_to_plot.insert(0, params['time_start_day'])
    # print("data at start ", params['time_start_day'] , ' = ' , params['total_covid_data'][params['time_start_day']])
    # for _date in params['dates_for_fit']:
        # print ('data at ' , _date , ' = ' , params['total_covid_data'][_date])
    # plotDataFrance(params, G, dates_to_plot)    
    plotDataHospitalFrance(params, G, dates_to_plot)    

# launch fit procedure 
# TODO: save df to some format or do something with it later on
if (params['fit_model']):
    print("Start fitting procedure, Inch'Allah...")
    paramsToFit = InitFit(params)     
    params['num_iter_fit'] = 0 # intialize number of iteration
    df = LaunchFit(paramsToFit, params, G)    

# reads params and outputs from previously generated files
# params.txt and outputs.txt located in direcotry params['outputs_dir']
# so far only work for specific models
if ("PlotResults" in params['plot_options']):
    if (params["model_name"] == "RealSEIRLDOnGraphInitTEST"):
        _dir=params['outputs_dir']
        fname=_dir+'/outputs.txt'
        s = open(fname, 'r').read()
        outputs = list(eval(s))
        
        fname=_dir+'/params.txt'
        s = open(fname, 'r').read()                
        params = dict(eval(s))
        
        # outputs, G = SolveOnGraph(params, graph)
        dates_to_plot = []
        # copy list otherwize references are copied... and it messes up
        # with params['dates_for_fit'] list later on (e.g. on the fit function)
        dates_to_plot = params['dates_for_fit'].copy()
        dates_to_plot.insert(0, params['time_start_day'])        
    else:
        print('[ERROR]: can not plot model ' , params["model_name"],  ' yet. Exit here')
        exit(0)

# launch model either with odes or not
# in case of ODEs only specific models are available
# for old tests with R0 dynamic
else:
    if (params['fit_with_ODEs']):
        params , G = InitModelODEs(params, G)
        print('coucou')
        
        outputs = LaunchModelODEs(params) 
        
        # save params inputs
        f = open(params['outputs_dir']+'/params.txt',"w")
        f.write(str(params))
        f.close()       
        
        # save outputs
        f = open(params['outputs_dir']+'/outputs.txt',"w")
        f.write(str(outputs))
        f.close()
        
        if params["model_name"] == "SIRWellMixed":
            R0 = params['k_SI'] * params['S_ini']/ (params['k_IR'])
            print('R0 = ' , R0)
            plot_SIRWellMixed(params, outputs)
            
        elif params["model_name"] == "RealSEIRLDWellMixed":      
            R0 = params['k_SE'] / (params['k_IR'] + params['k_IL'])
            print('R0 = ' , R0)
            dates_to_plot = []
            # copy list otherwize references are copied... and it messes up
            # with params['dates_for_fit'] list later on (e.g. on the fit function)
            dates_to_plot = params['dates_for_fit'].copy()
            dates_to_plot.insert(0, params['time_start_day'])
            plotDataAndSEIRLDResults(params, G,outputs, dates_to_plot)  
        elif params["model_name"] == "RealSIRLDWellMixed":      
            R0 = params['k_SI'] / (params['k_IR'] + params['k_IL'])
            print('R0 = ' , R0)
            dates_to_plot = []
            # copy list otherwize references are copied... and it messes up
            # with params['dates_for_fit'] list later on (e.g. on the fit function)
            dates_to_plot = params['dates_for_fit'].copy()
            dates_to_plot.insert(0, params['time_start_day'])
            plotDataAndSIRLDResults(params, G,outputs, dates_to_plot)  
        elif params["model_name"] == "RealSIRLDmodifiedWellMixed":      
            R0 = params['k_SI'] / (params['k_IR'] + params['k_IL'])
            print('R0 = ' , R0)
            dates_to_plot = []
            # copy list otherwize references are copied... and it messes up
            # with params['dates_for_fit'] list later on (e.g. on the fit function)
            dates_to_plot = params['dates_for_fit'].copy()
            dates_to_plot.insert(0, params['time_start_day'])
            plotDataAndSIRLDmodidifedResults(params, G,outputs, dates_to_plot)  
        elif params["model_name"] == "SIRWellMixed":
            R0 = params['k_SI'] / (params['k_IR'])
            print('R0 = ' , R0)
            plot_SIRWellMixed(params,outputs)
        elif params["model_name"] == "RealSEIsIasHoRadRDWellMixed":
            dates_to_plot = params['dates_for_fit'].copy()
            dates_to_plot.insert(0, params['time_start_day'])   
            
            plotDataAndRealSEIsIasHoRadRDWellMixed(params,G, outputs, dates_to_plot) 
            
        else:
            print('[ERROR]: can not plot model ' , params["model_name"],  ' yet. Exit here')
            exit(0)
    else:
        LaunchModel(params, G)

t1 = time.time()
print("Code ran for: ", t1 - t0, " seconds")
