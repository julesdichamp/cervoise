#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 25 22:04:48 2020

@author: florian
"""
import geopandas as gpd
import pandas as pd
import numpy as np
import glob

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

    # %% data airplanes
    airplaneFile = params['airplaneFile_name']
    airplanes = pd.read_excel(airplaneFile)
   
    # %% covid19 data
    covidDptsFileNames = params['allCovid19DepartmentsFileNames']
    covidRegionsFileNames = params['allCovid19RegionsFileNames']
    
    dic_data_covid = {}    
    for fname in covidDptsFileNames:
        _date = fname[-15:-5]
        df = pd.read_json(fname, orient='columns')
        for i in range(0,df.shape[0]):
            properties = df.iloc[i,0]['properties']
            dpt_name = properties['Province/State']
            if ('Confirmed'  in properties):
                num_confirmed = properties['Confirmed']
            else:
                num_confirmed = 0
            if ('Deaths' in properties):
                num_deaths = properties['Deaths']
            else:
                num_deaths = 0

            # rename it because it has a different name in other data set...
            if (dpt_name == "Territoire de Belfort"):
                dpt_name = "Territoire-de-Belfort"

            # check that department is not alerady in dictionnary
            if (not(dpt_name in dic_data_covid)):
                dic_data_covid[dpt_name] = {}
            dic_data_covid[dpt_name][_date] = {}
            dic_data_covid[dpt_name][_date]['Confirmed'] = num_confirmed
            dic_data_covid[dpt_name][_date]['Deaths'] = num_deaths
            # print(_date, dpt_name, num_confirmed, num_deaths)
    
    # the data set in departments does not include La Reunion for some reason
    # data for La Reunion is however available is directory regions-france
    # La Reunion is both a region and a department so we include it to
    # the data from regions folder
    for fname in covidRegionsFileNames:        
        _date = fname[-15:-5]
        df = pd.read_json(fname, orient='columns')
        for i in range(0,df.shape[0]):
            properties = df.iloc[i,0]['properties']
            dpt_name = properties['Province/State']
            if (dpt_name == "La Réunion" ):
                if ('Confirmed'  in properties):
                    num_confirmed = properties['Confirmed']
                else:
                    num_confirmed = 0
                if ('Deaths' in properties):
                    num_deaths = properties['Deaths']
                else:
                    num_deaths = 0
                if (not(dpt_name in dic_data_covid)):
                    dic_data_covid[dpt_name] = {}                    
                    
                dic_data_covid[dpt_name][_date] = {}
                dic_data_covid[dpt_name][_date]['Confirmed'] = num_confirmed
                dic_data_covid[dpt_name][_date]['Deaths'] = num_deaths
                # print(_date, dpt_name, num_confirmed, num_deaths)
    
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
        franceGeom[dptname]["agesSex"] = ages.loc[dpt].to_dict() # added to_dict() so it can be read easily from files        
        franceGeom[dptname]["freqGare"] = freqGares[dpt]
        franceGeom[dptname]["neighbours"] = {}
        franceGeom[dptname]["length"] = {}
        
        if (not(dptname in dic_data_covid)):
            print("[ERROR]: ", dptname , "not find in dic_data_covid. Exit here")
            exit(0)
        
        franceGeom[dptname]["data_covid"] = dic_data_covid[dptname]
        # print("\n\n", dptname, " centroid is ", dptGeom.centroid.coords[0])
        # print("\n\n", dptname, " ages/sex : ", franceGeom[dptname]["agesSex"])

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
        
    # add an ID that will be used in node 
    # IT IS NOT RELATED TO THE DEPARTMENT NUMBER 
    # BECAUSE IT STARTS FROM 0
    id_dpt = 0
    for dptname in franceGeom.keys():
        franceGeom[dptname]["id"] = id_dpt
        id_dpt+=1
    
    outname = params['data_dir'] + '/franceGeomDepartments.txt'    
    f = open(outname,"w")
    f.write( str(franceGeom) )
    
    return franceGeom
    

def getFranceGeomRegions(params):
    # we start from data at department scale    
    franceGeomDpts = getFranceGeomDpts(params)    
    df = pd.read_csv(params['departments_name'])
    dpts_names_from_list = df['NOM']
    all_correct = True
    for dpt_name in dpts_names_from_list:
        if not(dpt_name in franceGeomDpts):
            print(dpt_name, " not found in franceGeomDpts")
            all_correct = False
    
    if not(all_correct):
        print("we miss at least one department")
    else: 
        print("we found all departments")
        
    df = pd.read_csv(params['departments_name'])    
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

    LIST_DPTS_REGS_OUTREMER = [
        "Guadeloupe",
        "Guyane",
        "Martinique",
        "Mayotte",
        "La Réunion"        
    ]
    # LIST_OUTREMER = [
        # "Guadeloupe",
        # "Guyane",
        # "Martinique",
        # "Mayotte",
        # "La Réunion",
        # "Saint-Barthélémy",
        # "Saint-Martin",
    # ]
    LIST_REGIONS = LIST_METROPOLE + LIST_DPTS_REGS_OUTREMER
    
    dpts_regions_from_list = df['REGION']
    list_regions = []
    for r in dpts_regions_from_list:
        if not(r in list_regions):
            list_regions.append(r)
    
    all_correct = True
    for reg in LIST_REGIONS:
        if not(reg in list_regions):
            print(reg, " not found in LIST_REGIONS")
            all_correct = False
        
    if not(all_correct):
        print("we miss at least one region")
    else: 
        print("we found all regions")
        
    # for each  region get list of departments in that region
    franceGeomRegions = {}
    for regName in LIST_REGIONS:
        print(regName)
        # indices = [df['REGION'].index(regName)]    
        franceGeomRegions[regName] = {}
        indices = [i for i, r in enumerate(df['REGION']) if r == regName]
        franceGeomRegions[regName]['departments'] = df['NOM'][indices].to_list()
        x = 0.
        y = 0.
        franceGeomRegions[regName]["data_covid"] = {}
        franceGeomRegions[regName]["freqGare"] = 0.
        franceGeomRegions[regName]["agesSex"] = {}
                
        # loop over all departments that are part of the region
        for dptName in franceGeomRegions[regName]['departments']:
            # sum all confirmed and deaths from departments in region and 
            # set those values in region dictionnary
            # do it for each date there are data for
            for _date in franceGeomDpts[dptName]["data_covid"].keys():
                if not(_date in franceGeomRegions[regName]["data_covid"]):
                    franceGeomRegions[regName]["data_covid"][_date] = {}
                    franceGeomRegions[regName]["data_covid"][_date]['Confirmed'] = 0
                    franceGeomRegions[regName]["data_covid"][_date]['Deaths'] = 0
                franceGeomRegions[regName]["data_covid"][_date]['Confirmed'] += franceGeomDpts[dptName]["data_covid"][_date]['Confirmed']
                franceGeomRegions[regName]["data_covid"][_date]['Deaths'] += franceGeomDpts[dptName]["data_covid"][_date]['Deaths']
                
            # add train stations frequentations = sum over all train frequentations
            # in departments
            franceGeomRegions[regName]["freqGare"] += franceGeomDpts[dptName]["freqGare"]
            
            # add ages per sex data = sum over all data in each department
            for key in franceGeomDpts[dptName]["agesSex"].keys():
                # do not add those keys because they are useless now
                if ( (key != "dptname") and (key != "dpt") ):
                    # if the key is not in the franceGeomRegions dictionnary
                    # we initialize as 0 the number of people in respective age class
                    if not(key in franceGeomRegions[regName]["agesSex"].keys()):
                        franceGeomRegions[regName]["agesSex"][key] = 0.
                    # we sum the number of people in respective age class for all departments
                    else:
                        franceGeomRegions[regName]["agesSex"][key] += franceGeomDpts[dptName]["agesSex"][key]
            # get centroid from region as barycenter from all centroids departments
            x += franceGeomDpts[dptName]["centroid"][0] / len(franceGeomRegions[regName]['departments'])
            y += franceGeomDpts[dptName]["centroid"][1] / len(franceGeomRegions[regName]['departments'])                    
        
            
        # set centroid for region
        franceGeomRegions[regName]["centroid"] = (x,y)
        
        # prepare empty lists
        franceGeomRegions[regName]["neighbours"] = []
        
    # set neighbours of regions
    # loop over all departments that are part of the region
    for regName1 in franceGeomRegions.keys():        
        for regName2 in franceGeomRegions.keys():
            if (regName1 != regName2):
                for dptName1 in franceGeomRegions[regName1]['departments']:
                    for dptName2 in franceGeomRegions[regName2]['departments']:
                        if ((dptName1 in franceGeomDpts[dptName2]["neighbours"]) and
                            (dptName2 in franceGeomDpts[dptName1]["neighbours"])):
                                print(regName1, " is neighbour of ", regName2)
                                if not(regName1 in franceGeomRegions[regName2]['neighbours']):
                                    franceGeomRegions[regName2]['neighbours'].append(regName1)
                                if not(regName2 in franceGeomRegions[regName1]['neighbours']):
                                    franceGeomRegions[regName1]['neighbours'].append(regName2)                            

    
    # add an ID that will be used in node 
    # IT IS NOT RELATED TO THE DEPARTMENT NUMBER 
    # BECAUSE IT STARTS FROM 0
    id_reg = 0
    for regName in franceGeomRegions.keys():
        franceGeomRegions[regName]["id"] = id_reg
        id_reg+=1
    
    # print(franceGeomRegions)
    outname = params['data_dir'] + '/franceGeomRegions.txt'    
    f = open(outname,"w")
    f.write( str(franceGeomRegions))
    return franceGeomRegions
        
def readFranceGeom(params):
    if not(params['franceGeom_name']): # test if string is empty
        # in that case we need to process data files to get dictionnary
        # of data
        if (params["graph_type"] == "departments"):
            franceGeom = getFranceGeomDepartments(params)
        elif (params["graph_type"] == "regions") or (params["graph_type"] == "fullConnectedRegions"):
            franceGeom = getFranceGeomRegions(params)
        else:
            print("Not sure what to do with that graph_type = " , params["graph_type"] , ", exit here.")
    else:
        s = open(params['franceGeom_name'], 'r').read()
        franceGeom = eval(s)
        # print(franceGeom)
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
    
    airplaneFile_name = 'data/bulletin_stat_trafic_aerien_2018-converti.xlsx'
    # %% https://github.com/kalisio/covid-19/
    # %% full data is located in data/covid19_france_departments/
    # covid19DepartmentsDir_name = data/covid19_france_departments/departements-france/
    covid19DepartmentsDir_name = 'data/covid19_france_departments/departements-france/'
    covid19RegionsDir_name = 'data/covid19_france_departments/regions-france/'
    # list of regions and french departments
    #https://www.regions-et-departements.fr/departements-francais#departements_fichiers
    # some departments names were modified to fit other names
    departments_name= 'data/departments/departements-francais.csv'
    # contains a dictionnary with all franceGeom data
    params['airplaneFile_name'] = airplaneFile_name
    
    params['covid19DepartmentsDir_name'] = covid19DepartmentsDir_name
    params['covid19RegionsDir_name'] = covid19RegionsDir_name
    # get all file names we are interested in in this directory (do not keep polygons files)s
    params['allCovid19DepartmentsFileNames'] = glob.glob(params['covid19DepartmentsDir_name'] +'/*departements-france-*.json')
    params['allCovid19RegionsFileNames'] = glob.glob(params['covid19RegionsDir_name'] +'/*regions-france-2020*.json')
    params['departments_name'] = departments_name
    params['data_dir'] = 'data/'
    # franceGeomRegions = getFranceGeomDpts(params)
    franceGeomRegions = getFranceGeomRegions(params)

