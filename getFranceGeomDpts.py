#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 25 22:04:48 2020

@author: florian
"""
import geopandas as gpd
import pandas as pd
import numpy as np


def getFranceGeomDpts(params):
    # %% data spatiales                
    name_shapefile = params['shapefile_name']     
    shapefile = gpd.read_file(name_shapefile)
    
    # %% data train    
    gareFile = params['gareFile_name'] 
    gares = pd.read_csv(gareFile, delimiter=";")
    # sum by dpt
    freqGares = np.zeros(shapefile.shape[0])

    # %% data pop/age/sex    
    ageFile = params['ageFile_name']    
    ages = pd.read_excel(ageFile)

    # %% sort and sum by dpt data train
    for i in range(gares.shape[0]):
        gare = gares.loc[i]
        cp = str(gare["Code postal"])
        if len(cp) < 5:
            cp = "0" + cp
        dpt = int(cp[:2]) - 1  # parceque les dpts commencent à 01 et les arrays a 0
        freq = gare["Total Voyageurs 2018"]
        freqGares[dpt] += freq

    # %% create dict
    franceGeom = {}
    dptID = 0
    for dpt in range(shapefile.shape[0]):
        dptID += 1
        dptGeom = shapefile["geometry"][dpt]
        dptname = shapefile["nom"][dpt]
        if type(dptname) == bytes:
            dptname = dptname.decode("latin-1", "replace")
        else:
            dptname = dptname

        franceGeom[dptname] = {}
        franceGeom[dptname]["agesSex"] = ages.loc[dpt]
        franceGeom[dptname]["freqGare"] = freqGares[dpt]
        franceGeom[dptname]["neighbours"] = {}
        franceGeom[dptname]["length"] = {}

        print("\n\n", dptname, " centroid is ", dptGeom.centroid.coords[0])
        print("\n\n", dptname, " ages/sex : ", franceGeom[dptname]["agesSex"])

        franceGeom[dptname]["centroid"] = dptGeom.centroid.coords[0]
        # look for neighbours
        neighList = []
        lenList = []

        for nb in range(shapefile.shape[0] - 1):
            if nb != dpt:
                neighbour = shapefile["geometry"][nb]
                if dptGeom.touches(neighbour):
                    dptneighname = shapefile["nom"][nb]
                    if type(dptneighname) == bytes:
                        dptneighname = dptneighname.decode("latin-1", "replace")
                    else:
                        dptneighname = dptneighname
                    neighList.append(dptneighname)
                    l = dptGeom.boundary.intersection(neighbour.boundary)
                    l = l.length
                    lenList.append(l)

                    print(
                        "\t",
                        dptname,
                        " touches ",
                        dptneighname,
                        " on ",
                        l,
                        " m",
                    )

        franceGeom[dptname]["neighbours"] = neighList
        franceGeom[dptname]["length"] = lenList
    return franceGeom    

if __name__ == "__main__":
    #  https://www.data.gouv.fr/fr/datasets/r/ab21d892-aa39-466a-a3ed-dc25ff076b78
    name_shapefile = "data/departments/dpts/departements-20140306-5m.shp"
    #  https://www.data.gouv.fr/fr/datasets/frequentation-en-gares/
    gareFile = "data/departments/frequentation-gares.csv"
    #  https://www.insee.fr/fr/statistiques/1893198
    ageFile = "data/departments/estim-pop-dep-sexe-gca-1975-2020.xls"
    params = {}
    params['shapefile_name'] = name_shapefile
    params['gareFile_name']  = gareFile
    params['ageFile_name']   = ageFile
    getFranceGeomDpts(params)



