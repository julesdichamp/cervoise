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
import parser

import momepy
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import networkx as nx

# from regions import 

# sys.path.append("/home/florian/PycharmProjects/drunkWalker/")
# import computeMSD

print("\nCERVOISE = Covid19 : Estimation of Regional Virus Organization Infered using Stochastic Equation = Starting")
t0 = time.time()

# =============================================================================
# read input file. Either specified by argument or generic inputParams.ini
# =============================================================================
if len(sys.argv) > 1:
    inputFile = sys.argv[1]
else:
    inputFile = "config.ini"



def LoadRoads(road_names):
	print('LoadRoads')
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
	print('InitParams')
	config = configparser.ConfigParser()
	config.read(config_name)		
	params={}
	# load input data infos
	params['graph_name'] = config.get('input_data','graph_name')
	params['data_fname'] = config.get('input_data','data_name')
	params['roads_name'] = config.get('input_data','roads_name')
	
	# load options infos
	params['location'] = config.get('options','location')
	return params


def InitGraph(params):
	print('InitGraph')
	LIST_METROPOLE = [
	'Auvergne-Rhône-Alpes',
	'Bourgogne-Franche-Comté',
	'Bretagne',
	'Centre-Val de Loire',
	'Corse',
	'Grand Est',
	'Hauts-de-France',
	'Ile-de-France',
	'Normandie',
	'Nouvelle-Aquitaine',
	'Occitanie',
	'Pays de la Loire',
	"Provence-Alpes-Côte d’Azur"]
	
	LIST_OUTREMER = [
	"Guadeloupe",
	"Guyane",
	"Martinique",
	"Mayotte",
	"La Réunion",
	"Saint-Barthélémy",
	"Saint-Martin" ]
	
	# from kaggle data on death/reported cases in France animated map
	DIC_POSITIONS_METROPOLE = {'Auvergne-Rhône-Alpes': [45.2968119,4.6604809],
							 'Bourgogne-Franche-Comté' : [47.0510946,5.0740568],							 
							 'Bretagne' : [48.2640845,-2.9202408],
							 'Centre-Val de Loire' : [47.5490251,1.7324062],
							 'Corse' : [42.188089649999995,9.068413771427695],
							 'Grand Est' : [48.4845157,6.113035],
							 'Hauts-de-France' : [50.1024606,2.7247515],
							 'Ile-de-France' : [48.6443057,2.7537863],
							 'Normandie' : [49.0677708,0.3138532],
							 'Nouvelle-Aquitaine' : [45.4039367,0.3756199],
							 'Occitanie' : [43.6487851,2.3435684],
							 'Pays de la Loire' : [47.6594864,-0.8186143],
							 'Provence-Alpes-Côte d’Azur' : [44.0580563,6.0638506]
							 }
							 
	DIC_POSITIONS_OUTREMER = {'Guadeloupe' : [16.2490067,-61.5650444],
							 'Guyane' : [4.8417097,-58.6416891],
							 'Martinique' : [14.6113732,-60.9620777],
							 'Mayotte' : [-12.823048,45.1520755],
							 'La Réunion' : [-21.1309332,55.5265771],
							 'Saint-Barthélémy' : [17.9036287,-62.811568843006896],
							 'Saint-Martin' : [18.06481835,-63.049171129033]
							 }

	
	# dictionnary that associates an id to a region
	dic_regions_metropole = {}
	dic_regions_metropole_positions = {}
	for i in range(0,len(LIST_METROPOLE)):
		dic_regions_metropole[i] = LIST_METROPOLE[i]
		dic_regions_metropole_positions[i] = DIC_POSITIONS_METROPOLE[LIST_METROPOLE[i]]
	
	dic_regions_outremer = {}
	dic_regions_outremer_positions = {}
	for i in range(0,len(LIST_OUTREMER)):
		dic_regions_outremer[len(dic_regions_metropole) + i] = LIST_OUTREMER[i]
		dic_regions_outremer_positions[len(dic_regions_metropole) + i] = DIC_POSITIONS_OUTREMER[LIST_OUTREMER[i]]
	
	# set full dictionnary depending on option
	dic_regions={}
	dic_regions_positions={}
	if (params['location']=="metropole"):
		dic_regions = dic_regions_metropole
		dic_regions_positions = dic_regions_metropole_positions
	elif (params['location']=="outremer"):
		dic_regions = dic_regions_outremer
		dic_regions_positions = dic_regions_outremer_positions
	else:
		for x in dic_regions_metropole.keys() :
			dic_regions[x] = dic_regions_metropole[x]
			dic_regions_positions[x] = dic_regions_metropole_positions[x]
		for x in dic_regions_outremer.keys() :
			dic_regions[x] = dic_regions_outremer[x]
			dic_regions_positions[x] = dic_regions_outremer_positions[x]

	G = nx.Graph()
	for x in dic_regions.keys():
		print(x, dic_regions[x])
		G.add_node(x)	
		
	# All neighborings connections ONLY FOR METROPOLE NOW
	# Corse to be better dealt with?	
	dic_all_connections = {}
	dic_all_connections["Bretagne"] = ["Normandie", "Pays de la Loire"]
	dic_all_connections["Normandie"] = ["Bretagne", "Pays de la Loire", "Centre-Val de Loire", "Ile-de-France", "Hauts-de-France"]
	dic_all_connections["Pays de la Loire"] = ["Bretagne", "Normandie", "Centre-Val de Loire", "Nouvelle-Aquitaine"]
	dic_all_connections["Centre-Val de Loire"] = ["Pays de la Loire", "Normandie", "Ile-de-France", "Bourgogne-Franche-Comté", "Auvergne-Rhône-Alpes", "Nouvelle-Aquitaine"]
	dic_all_connections["Hauts-de-France"] = ["Grand Est", "Ile-de-France", "Normandie"]
	dic_all_connections["Nouvelle-Aquitaine"] = ["Pays de la Loire", "Centre-Val de Loire", "Auvergne-Rhône-Alpes", "Occitanie"]
	dic_all_connections["Grand Est"] = ["Hauts-de-France", "Ile-de-France", "Bourgogne-Franche-Comté"]
	dic_all_connections["Bourgogne-Franche-Comté"] = ["Ile-de-France", "Grand Est", "Centre-Val de Loire", "Auvergne-Rhône-Alpes"]
	dic_all_connections["Auvergne-Rhône-Alpes"] = ["Occitanie", "Nouvelle-Aquitaine", "Centre-Val de Loire", "Bourgogne-Franche-Comté", "Provence-Alpes-Côte d’Azur"]
	dic_all_connections["Provence-Alpes-Côte d’Azur"] = ["Auvergne-Rhône-Alpes", "Occitanie", "Corse"]
	dic_all_connections["Occitanie"] = ["Nouvelle-Aquitaine", "Auvergne-Rhône-Alpes", "Provence-Alpes-Côte d’Azur"]
	dic_all_connections["Corse"] = ["Provence-Alpes-Côte d’Azur"]
	dic_all_connections["Ile-de-France"] = LIST_METROPOLE
	
	for key in dic_all_connections.keys() :
		# get id of region
		id1 = list(dic_regions.keys())[list(dic_regions.values()).index(key)]
		for value in dic_all_connections[key]: # loop over all connections regions
			# get id of conected region
			id2 = list(dic_regions.keys())[list(dic_regions.values()).index(value)]
			edge = (id1, id2)
			G.add_edge(*edge)  # unpack edge tuple*
	
	# print("edges = " , G.edges())	
	fig, ax = plt.subplots()
	nx.draw(G, with_labels=True, font_weight='bold', pos=dic_regions_positions)	
	ax.tick_params(left=True, bottom=True, labelleft=True, labelbottom=True)
	plt.show()

# =============================================================================
# LAUNCH WORKFLOW
# =============================================================================
params = InitParams(inputFile)
graph = InitGraph(params)

# LoadRoads(params['roads_name'])
