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
import collections

# path to Vicopo.py must be added in PYTHONPATH
# from Vicopo import Vicopo


#WARNING: there exists in destination of airways between departments
# and outremer which are not departments:
# Polynésie-française
# there are also internal lines between outre-mer territories that are not covered yet
def getFranceGeomDpts(params):
    print('getFranceGeomDpts')
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
    airways = pd.read_excel(airplaneFile)

    citiesFile = params['citiesFile_name']
    cities = pd.read_csv(citiesFile)
    airwaysDict = {}
    for airway in airways.iterrows():
        # print(airway[1]["in"], airway[1]["out"])
        if not airway[1]["in"] in airwaysDict:
            airwaysDict[airway[1]["in"]] = [[airway[1]["out"], airway[1]["Pax"]]]
        else:
            airwaysDict[airway[1]["in"]].append([airway[1]["out"], airway[1]["Pax"]])
    #building the dict
    # keys = dpts numbers
    # values = list of [dpt number arrival, total pax]
    airwaysDictdpt = {}
    for flightSource in airwaysDict:
        # print(flightSource, cities[cities.eq(flightSource).any(1)].dpt.values[0])
        flightSourcedpt = cities[cities.eq(flightSource).any(1)].dpt.values[0]
        if not flightSourcedpt in airwaysDictdpt:
            if not flightSourcedpt in airwaysDictdpt:
                airwaysDictdpt[flightSourcedpt] = {}
            for flightDest in airwaysDict[flightSource]:
                
                # print(flightDest[0])
                # print(cities[cities.eq(flightDest[0]).any(1)].dpt.values[0])
                flightDestdpt = cities[cities.eq(flightDest[0]).any(1)].dpt.values[0]
                if not bool(airwaysDictdpt[flightSourcedpt]):
                    airwaysDictdpt[flightSourcedpt] = [[flightDestdpt, int(flightDest[1])]]
                elif not flightDestdpt in airwaysDictdpt[flightSourcedpt]:
                    airwaysDictdpt[flightSourcedpt].append([flightDestdpt, int(flightDest[1])])
                else:
                    airwaysDictdpt[flightSourcedpt][airwaysDictdpt[flightSourcedpt][0]==flightDestdpt] += int(flightDest[1])
    # print (airwaysDictdpt)    
    # cleaning dict as there are several flightSource dpt

    # %% covid19 data
    covidDptsFileNames = params['allCovid19DepartmentsFileNames']
    covidRegionsFileNames = params['allCovid19RegionsFileNames']

    dic_data_covid = {}
    for fname in covidDptsFileNames:
        _date = fname[-15:-5]
        df = pd.read_json(fname, orient='columns')
        for i in range(0,df.shape[0]):
            properties = df.iloc[i,1]['properties']
            dpt_name = properties['Province/State']
            if ('Confirmed'  in properties):
                num_confirmed = properties['Confirmed']
            else:
                num_confirmed = 0
            if ('Deaths' in properties):
                num_deaths = properties['Deaths']
            else:
                num_deaths = 0
            if ('Recovered'  in properties):
                num_recovered = properties['Recovered']
            else:
                num_recovered = 0
            # rename it because it has a different name in other data set...
            if (dpt_name == "Territoire de Belfort"):
                dpt_name = "Territoire-de-Belfort"

            # check that department is not alerady in dictionnary
            if (not(dpt_name in dic_data_covid)):
                dic_data_covid[dpt_name] = {}
            dic_data_covid[dpt_name][_date] = {}
            dic_data_covid[dpt_name][_date]['Confirmed'] = num_confirmed
            dic_data_covid[dpt_name][_date]['Deaths'] = num_deaths
            dic_data_covid[dpt_name][_date]['Recovered'] = num_recovered
            # print(_date, dpt_name, num_confirmed, num_deaths)

    # the data set in departments does not include La Reunion for some reason
    # data for La Reunion is however available is directory regions-france
    # La Reunion is both a region and a department so we include it to
    # the data from regions folder
    for fname in covidRegionsFileNames:
        _date = fname[-15:-5]
        df = pd.read_json(fname, orient='columns')
        for i in range(0,df.shape[0]):
            properties = df.iloc[i,1]['properties']
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
                if ('Recovered'  in properties):
                    num_recovered = properties['Recovered']
                else:
                    num_recovered = 0
                if (not(dpt_name in dic_data_covid)):
                    dic_data_covid[dpt_name] = {}

                dic_data_covid[dpt_name][_date] = {}
                dic_data_covid[dpt_name][_date]['Confirmed'] = num_confirmed
                dic_data_covid[dpt_name][_date]['Deaths'] = num_deaths
                dic_data_covid[dpt_name][_date]['Recovered'] = num_recovered
                # print(_date, dpt_name, num_confirmed, num_deaths)
    
    _list = [x for x in dic_data_covid.keys()]
    
    # %% sort and sum by dpt data train
    for i in range(gares.shape[0]):
        gare = gares.loc[i]
        cp = str(gare["Code postal"])
        if len(cp) < 5:
            cp = "0" + cp
        dpt = int(cp[:2]) - 1  # parceque les dpts commencent à 01 et les arrays a 0
        freq = gare["Total Voyageurs 2018"]
        freqGares[dpt] += freq

    # read departments name
    df = pd.read_csv(params['departments_name'], delimiter=',')
    nums  = df['NUMÉRO'].to_list()
    names = df['NOM'].to_list()
    dic_dptname_to_num = {}
    dic_dptnum_to_name = {}
    for i in range(0,len(nums)):
        dic_dptname_to_num[names[i]] = nums[i]
        dic_dptnum_to_name[nums[i]] = names[i]    
    
    # %% create dict
    franceGeom = {}
    dptID = 0
    for dpt in range(shapefile.shape[0]):        
        # dptID+=1 
        dptGeom = shapefile["geometry"][dpt]
        dptname = shapefile["nom"][dpt]
        if type(dptname) == bytes:
            dptname = dptname.decode("latin-1", "replace")
        else:
            dptname = dptname
        if not (dptname in dic_dptname_to_num):
            print("[ERROR]: we didn't find ",dptname, "in dic_dptname_to_num. Exit here.")
            exit(0)
        else:
            dptID = dic_dptname_to_num[dptname]
        # print(dpt, dptID, dptname)
        
        franceGeom[dptname] = {}
        franceGeom[dptname]["agesSex"] = ages.loc[dpt].to_dict() # added to_dict() so it can be read easily from files
        franceGeom[dptname]["freqGare"] = freqGares[dpt]
        franceGeom[dptname]["neighbours"] = {}
        franceGeom[dptname]["length"] = {}                
        
        if (not(dptname in dic_data_covid)):
            print("[ERROR]: ", dptname , "not find in dic_data_covid. Exit here")
            exit(0)
        
        # add airways data to dic (rename dptID into dpt names to match the keys
        franceGeom[dptname]["airways"] = {}
        if (str(dptID) in airwaysDictdpt):            
            franceGeom[dptname]["airways"] = [ [ dic_dptnum_to_name[x[0]] , x[1] ] for x in airwaysDictdpt[str(dptID)]]            
            # for x in franceGeom[dptname]["airways"]:
                # if (not(x[0] in _list)):
                    # print( x , " is not in franceGeom keys")
                    
        # add data_covid to dic
        franceGeom[dptname]["data_covid"] = dic_data_covid[dptname]

        # add centroid to dic
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


        franceGeom[dptname]["neighbours"] = neighList
        franceGeom[dptname]["length"] = lenList

    # add an ID that will be used in node
    # IT IS NOT RELATED TO THE DEPARTMENT NUMBER
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
    print('getFranceGeomRegions')
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

    # for each  region get list of departments in that region
    franceGeomRegions = {}
    for regName in LIST_REGIONS:
        # indices = [df['REGION'].index(regName)]
        franceGeomRegions[regName] = {}
        indices = [i for i, r in enumerate(df['REGION']) if r == regName]
        franceGeomRegions[regName]['departments'] = df['NOM'][indices].to_list()
        x = 0.
        y = 0.
        franceGeomRegions[regName]["data_covid"] = {}
        franceGeomRegions[regName]["freqGare"] = 0.
        franceGeomRegions[regName]["agesSex"] = {}
        franceGeomRegions[regName]["airways"] = [] 
        # loop over all departments that are part of the region
        for dptName in franceGeomRegions[regName]['departments']:
            # sum all confirmed and deaths from departments in region and
            # set those values in region dictionnary
            # do it for each date there are data for
            for _date in franceGeomDpts[dptName]["data_covid"].keys():
                if not(_date in franceGeomRegions[regName]["data_covid"]):
                    franceGeomRegions[regName]["data_covid"][_date] = {}
                    franceGeomRegions[regName]["data_covid"][_date]['Confirmed'] = 0
                    franceGeomRegions[regName]["data_covid"][_date]['Deaths'] = 0
                    franceGeomRegions[regName]["data_covid"][_date]['Recovered'] = 0
                franceGeomRegions[regName]["data_covid"][_date]['Confirmed'] += franceGeomDpts[dptName]["data_covid"][_date]['Confirmed']
                franceGeomRegions[regName]["data_covid"][_date]['Deaths'] += franceGeomDpts[dptName]["data_covid"][_date]['Deaths']
                franceGeomRegions[regName]["data_covid"][_date]['Recovered'] += franceGeomDpts[dptName]["data_covid"][_date]['Recovered']

            # add train stations frequentations = sum over all train frequentations
            # in departments
            franceGeomRegions[regName]["freqGare"] += franceGeomDpts[dptName]["freqGare"]
            
            # add ages per sex data = sum over all data in each department
            for key in franceGeomDpts[dptName]["agesSex"].keys():
                # do not add those keys because they are useless now
                if ( (key != "dptname") and (key != "dpt") ):
                    # if the key is not in the franceGeomRegions dictionnary
                    # we initialize as 0 the number of people in respective age class
                    if not(key in franceGeomRegions[regName]["agesSex"].keys()):
                        franceGeomRegions[regName]["agesSex"][key] = 0.
                    # we sum the number of people in respective age class for all departments
                    else:
                        franceGeomRegions[regName]["agesSex"][key] += franceGeomDpts[dptName]["agesSex"][key]
            # get centroid from region as barycenter from all centroids departments
            x += franceGeomDpts[dptName]["centroid"][0] / len(franceGeomRegions[regName]['departments'])
            y += franceGeomDpts[dptName]["centroid"][1] / len(franceGeomRegions[regName]['departments'])
            
        # set centroid for region
        franceGeomRegions[regName]["centroid"] = (x,y)

        # prepare empty lists
        franceGeomRegions[regName]["neighbours"] = []
        
    # set neighbours of regions
    # loop over all departments that are part of the region
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
    
    # add neighbours based on airways connections 
    # also iterate on pax to get it at region scale
    for regName1 in franceGeomRegions.keys():
        for dptName1 in franceGeomRegions[regName1]['departments']:
            list_airplanes_connections = franceGeomDpts[dptName1]['airways']
            for connection in list_airplanes_connections:
                dptname_connected = connection[0]
                pax = connection[1]
                # look at which region dpt_connected corresponds
                for regName2 in franceGeomRegions.keys():
                    # we found the corresponding region
                    if (dptname_connected in franceGeomRegions[regName2]['departments']):
                        regName_connected = regName2
                
                # check that this connection between region does not already exist
                _count = 0
                connection_exists = False
                for check_connection in franceGeomRegions[regName1]["airways"]:
                    # this already region already had a connection
                    if (regName_connected == check_connection[0]):                        
                        indice = _count
                        connection_exists = True
                    _count += 1
                # if a connection already exists we add pax value to already 
                # existing pax value
                if (connection_exists):
                    franceGeomRegions[regName1]["airways"][indice][1] += pax
                # we the pax and the new connection to the list ['airways']
                else:
                    franceGeomRegions[regName1]["airways"].append([regName_connected,pax])
    
    # add an ID that will be used in node
    # IT IS NOT RELATED TO THE DEPARTMENT NUMBER
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
    if not(params['franceGeom_name'].strip()): # test if string is empty (or filled with blanks)
        # in that case we need to process data files to get dictionnary
        # of data
        if (params["graph_type"] == "departments"):
            franceGeom = getFranceGeomDepartments(params)
        elif (params["graph_type"] == "regions") or (params["graph_type"] == "fullConnectedRegions"):
            franceGeom = getFranceGeomRegions(params)
        else:
            print("Not sure what to do with that graph_type = " , params["graph_type"] , ", exit here.")
    else:
        print('read franceGeom file: ' , params['franceGeom_name'])
        s = open(params['franceGeom_name'], 'r').read()
        franceGeom = eval(s)
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

    airplaneFile_name = 'data/bulletin_stat_trafic_aerien_2018-france_passager.xlsx'

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

    citiesFile_name = 'data/villes_france_correct.csv'

    params['citiesFile_name'] = citiesFile_name

    params['data_dir'] = 'data/'
    # franceGeomRegions = getFranceGeomDpts(params)
    franceGeomRegions = getFranceGeomRegions(params)

