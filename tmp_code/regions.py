########
# import
#

import shutil
import os
import re
import pandas as pd
import numpy as np
import datetime
import locale
import platform
# for geo lat/lon
from geopy.geocoders import Nominatim

import scrapy

import pdfminer_to_text


########
# PATCHS
#

# patch Month in french
#locale.setlocale(locale.LC_ALL, 'fr_FR')
# patch VSCODE issue : https://github.com/Microsoft/vscode-python/issues/5468

##################
# useful functions
#
# save before scraping
def clean_file(path_file_name):
    '''
    Clean file already traited : rename file with date
    '''
    try:
        d = datetime.datetime.now()
        str_date = '_' + d.strftime("%Y%m%d_%H_%M_%S")
       
        res_re = re.search('\.\w+$', path_file_name)
        
        path_file_name_saved = \
            path_file_name[0:res_re.start()] + str_date + res_re.group(0)
         
        shutil.move(path_file_name, path_file_name_saved) 
        print('File {} moved!'.format(path_file_name_saved))
    except:
        print('File {} does not exist!'.format(path_file_name))

# manage accents
def conv_noaccent(chaine):
    '''
    Convert string without accents

    ex: chaine = "mà chàîne âvèc dès càräctères spéciâüx"
        print(conv_noaccent(chaine))
        returns : 'ma chaine avec des caracteres speciaux'
    '''
    import unicodedata
    
    return unicodedata.normalize('NFKD', chaine).encode('ASCII', 'ignore').decode('ASCII')
    
#chaine = "mà chàîne âvèc dès càräctères spéciâüx"
#print(conv_noaccent(chaine))

def create_df_tables_fr():
    #
    # Create a dataFrame to load data
    #
    return pd.DataFrame(columns=LIST_COL_FRA, index=[0])

def update_df_fr(df_add):
    #
    # Update dataFrame save date by date for areas of France
    #
    # First : 
    # If save file doesn't exsit create it
    # else  load it.
    # After : 
    # if new date, add the row
    # else update the row.
    # 
    if not os.path.exists(PATH_DF_FR):
        if MODE_PLATFORM == "KAGGLE":
            df_fr = pd.read_csv(URL_DF_FR)
        else:
            df_fr = df_add.copy()
    else:
        df_fr = pd.read_csv(PATH_DF_FR)
        # add / update 
        for iloc_idx in range(df_add.shape[0]):
            date_curr = df_add.iloc[iloc_idx]["date"]
            # if date exist, replace
            if date_curr in df_fr["date"].values:
                # num row in df_fr 
                num_row_df_fr = df_fr["date"].isin([date_curr]).values.argmax()
                # update this row
                df_fr.iloc[num_row_df_fr] = df_add.iloc[iloc_idx]
            # if new date, add
            else: 
                df_fr = df_fr.append(df_add.iloc[iloc_idx])
    # sort by date
    df_fr.sort_values(by=['date'], inplace=True)
    # reindex
    df_fr.reset_index(drop=True, inplace=True)
    # save old csv file
    clean_file(PATH_DF_FR)
    # save csv
    df_fr.to_csv(PATH_DF_FR, index=False)

    return df_fr


#############
# DEFINITIONS
#

# for LOCAL mode : 
PATH_INPUT = 'data/covid19_france_sources/'
PATH_TO_SAVE_DATA = "."

PATH_INPUT_KAGGLE = '/kaggle/input/covid19-france-sources/sources' #for kaggle
PATH_TO_SAVE_DATA_KAGGLE = '/kaggle/working' #for kaggle
PATH_GOOGLE_COLAB = '' #IF YOU USE GOOGLE COLAB, ENTER YOUR WORKING FOLDER
# get computer name
COMPUTERNAME = platform.node()
# select current platform
MY_PLATFORM = platform.system()
# check if Google colab need Drive ?
if re.match("^/content", os.getcwd()):
    MODE_PLATFORM = "GOOGLE_COLAB"
    #print("GOOGLE COLAB MODE")
    from google.colab import drive
    drive.mount('/content/drive', force_remount=True)
    os.chdir(PATH_GOOGLE_COLAB)

elif re.match("^/kaggle", os.getcwd()):
    MODE_PLATFORM = "KAGGLE"
    #print("KAGGLE COLAB MODE")
    PATH_INPUT = PATH_INPUT_KAGGLE 
    PATH_TO_SAVE_DATA = PATH_TO_SAVE_DATA_KAGGLE 
    # !pip install pdfminer.six
else:
    MODE_PLATFORM = "LOCAL"
    
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

TOT_METROPOLE = "Total Métropole"

LIST_OUTREMER = [
"Guadeloupe",
"Guyane",
"Martinique",
"Mayotte",
"La Réunion",
"Saint-Barthélémy",
"Saint-Martin" ]

TOT_OUTREMER = "Total Outre Mer"

LIST_SUBTOT = [
TOT_METROPOLE,
TOT_OUTREMER]

LIST_AREA_FR = LIST_METROPOLE + LIST_OUTREMER

# list of formatted area for scrapy Items
LIST_ITEM = []
for area_curr in LIST_AREA_FR:
    LIST_ITEM.append(area_curr.replace("-","_").replace("ô","o").\
    replace("é","e").replace("’","_").replace(" ","_"))

PATH_CSV = PATH_INPUT + '/csv_fr'
LIST_FILE_CSV = ['spf_web_front_2020_02_28.csv', 'spf_web_front_2020_03_02.csv']

PATH_PDF = PATH_INPUT + '/pdf_fr'
PATH_SPF_FRONT_OUTPUT = PATH_TO_SAVE_DATA + '/' + 'table_spf_front_update.json'

URL_SPF_FRONT = 'https://www.santepubliquefrance.fr/maladies-et-traumatismes/' + \
    'maladies-et-infections-respiratoires/infection-a-coronavirus/articles/' + \
    'infection-au-nouveau-coronavirus-sars-cov-2-covid-19-france-et-monde'

LIST_COL_FRA = ["date", "url", "source", "country","iso_alpha", "total"] + LIST_SUBTOT + \
    LIST_METROPOLE + LIST_OUTREMER + ["tot_metrop_calc", "tot_outre_calc", 
        "check_metrop", "check_outre", "check_tot"]

PATH_DF_FR = PATH_TO_SAVE_DATA + '/' + 'df_fr.csv'
URL_DF_FR = "https://raw.githubusercontent.com/jeugregg/coronavirusModel/master/df_fr.csv"
PATH_CSV_DF_GEO_FRA = PATH_TO_SAVE_DATA + '/' + 'df_geo_fra.csv'

# RETRIEVED DATA FROM WORLDOMETERS
PATH_TABLES_WORLDO = PATH_TO_SAVE_DATA + '/' + 'table_worldo_fr.json'


###################################
# CSV file for old datas (from spf)
#
df_tables_csv_fr = create_df_tables_fr()

'/Users/gregory/Documents/CloudStationSinchon/Applications/python/CoronaVirus/code/coronavirusModel/sources/csv_fr'
df_tables_csv_fr.at[0, "date"] = '2020-02-28'
df_tables_csv_fr.at[0, "url"] = 'https://web.archive.org/web/20200229212809/https://www.santepubliquefrance.fr/maladies-et-traumatismes/maladies-et-infections-respiratoires/infection-a-coronavirus/articles/covid-19-situation-epidemiologique-en-france'

df_tables_csv_fr.at[1, "date"] = '2020-03-02'
df_tables_csv_fr.at[1, "url"] = 'https://web.archive.org/web/20200303124836/https://www.santepubliquefrance.fr/maladies-et-traumatismes/maladies-et-infections-respiratoires/infection-a-coronavirus/articles/covid-19-situation-epidemiologique-en-france'

df_tables_csv_fr["source"] = "santepubliquefrance.fr"
df_tables_csv_fr["country"] = "France"
df_tables_csv_fr["iso_alpha"] = "FRA"



for idx_csv in range(len(LIST_FILE_CSV)):
    df_csv = pd.read_csv(PATH_CSV + '/' + LIST_FILE_CSV[idx_csv])
    # total
    df_tables_csv_fr.at[idx_csv, 'total'] = df_csv[df_csv["area"] == "Total"]["nb_cases"].values[0]
    # total metropole
    df_tables_csv_fr.at[idx_csv, TOT_METROPOLE] =  \
        df_csv[df_csv["area"] == 'France métropolitaine']["nb_cases"].values[0]
    for area_curr in LIST_METROPOLE + LIST_OUTREMER:
        nb_cases_found = df_csv[df_csv["area"] == area_curr]["nb_cases"].values
        if nb_cases_found.shape[0] == 0:
            nb_cases_found = 0
        else:
            nb_cases_found = nb_cases_found[0]
        df_tables_csv_fr.at[idx_csv, area_curr] = nb_cases_found
    
    # tot outremer
    df_tables_csv_fr[TOT_OUTREMER] = df_tables_csv_fr.filter(items=LIST_OUTREMER).sum(axis=1)

# add csv data to final dataframe
df_fr = update_df_fr(df_tables_csv_fr)   


###################
# PDF SOURCES FILES
# 

# check DATA input folder
list_pdf = []
for dirname, _, filenames in os.walk(PATH_PDF):
    for filename in filenames:
        list_pdf.append(os.path.join(dirname, filename))
        
df_pdf_fr = pd.DataFrame(columns=["url", "source", "path"], index=range(len(list_pdf)))
df_pdf_fr["source"] = "santepubliquefrance.fr"
df_pdf_fr["path"] = list_pdf


##########
# READ PDFs
#
from pdfminer_to_text import *

def find_pdf_fr_total(str_text):
    #
    # Find french total into text of PDF
    #
    
    # "Cas confirmés Nombre de cas 949"
    sub_pre = "Cas confirmés Nombre de cas"
    str_pattern = "(?<={})\s+\d+".format(sub_pre)
    list_found = re.search(str_pattern, str_text)
    if list_found != None:
        return np.int(list_found[0].strip())
    else:
        return np.nan

def find_pdf_fr_date(str_text):
    #
    # Find french date into texte of PDF 
    #
    re_found = re.search("(?<=Point\sépidémiologique\s-\sSituation\sau\s)\d+\s\w+\s\d\d\d\d", str_text)    
    if re_found != None:
        print((re_found.span))
        # str_date_fr = re_found[0]
        str_date_fr = str(re_found)[0]
        str_date_fr = str_date_fr.replace("mars", "March") # DIRTY PATCH... TODO 
        return datetime.datetime.strptime(str_date_fr,'%d %B %Y').strftime("%Y-%m-%d")
    else:
        return None


for idx_pdf in range(df_pdf_fr.shape[0]):

    pdf_file = df_pdf_fr.at[idx_pdf, "path"]
    str_text = extract_text_from_pdf(pdf_file)
    if idx_pdf == 0:
        df_tables_pdf_fr = create_df_tables_fr()

    df_tables_pdf_fr.at[idx_pdf, "date"] = find_pdf_fr_date(str_text)
    df_tables_pdf_fr.at[idx_pdf, "country"] = "France"
    df_tables_pdf_fr.at[idx_pdf, "iso_alpha"] = "FRA"
    df_tables_pdf_fr.at[idx_pdf, "total"] = find_pdf_fr_total(str_text)
    df_tables_pdf_fr.at[idx_pdf, "source"] = 'santepubliquefrance.fr'
    #Région   Auvergne-Rhône-Alpes 102 Bourgogne-Franche-Comté 129 Bretagne 40 
    #Centre-Val de Loire  16 Corse 5 Grand Est 250 Hauts-de-France 173 
    #Ile-de-France 104 Normandie 11 Nouvelle-Aquitaine 17 Occitanie 36 
    #Pays de la Loire 18 Provence-Alpes-Côte d’Azur 38 Total Métropole  939 
    #Guadeloupe 0 Guyane 5 Martinique 2 Mayotte 0 La Réunion  0 Saint-Barthélemy
    #1 Saint-Martin 2 Total Outre Mer 10

    for reg_curr in LIST_METROPOLE + LIST_OUTREMER + LIST_SUBTOT:
        str_pattern = "(?<={})\s+\d+".format(conv_noaccent(reg_curr))
        list_found = re.search(str_pattern, conv_noaccent(str_text))
        if list_found != None:
            df_tables_pdf_fr.at[idx_pdf, reg_curr] = \
                np.int(list_found[0].strip())
        else:
            df_tables_pdf_fr.at[idx_pdf, reg_curr] = 0
    # check total
    if df_tables_pdf_fr.at[idx_pdf, "total"] is np.nan:
        df_tables_pdf_fr.at[idx_pdf, "total"] = \
            df_tables_pdf_fr.at[idx_pdf, TOT_METROPOLE] + \
            df_tables_pdf_fr.at[idx_pdf, TOT_OUTREMER]
    # check total metropole
    df_tables_pdf_fr.at[idx_pdf, "tot_metrop_calc"] = 0
    for reg_curr in LIST_METROPOLE:
        df_tables_pdf_fr.at[idx_pdf, "tot_metrop_calc"] = \
            df_tables_pdf_fr.at[idx_pdf, "tot_metrop_calc"] + \
            df_tables_pdf_fr.at[idx_pdf, reg_curr]
    df_tables_pdf_fr.at[idx_pdf, "check_metrop"] = \
        df_tables_pdf_fr.at[idx_pdf, "tot_metrop_calc"] == \
        df_tables_pdf_fr.at[idx_pdf, TOT_METROPOLE]
    # check total outre mer
    df_tables_pdf_fr.at[idx_pdf, "tot_outre_calc"] = 0
    for reg_curr in LIST_OUTREMER:
        df_tables_pdf_fr.at[idx_pdf, "tot_outre_calc"] = \
            df_tables_pdf_fr.at[idx_pdf, "tot_outre_calc"] + \
            df_tables_pdf_fr.at[idx_pdf, reg_curr]
    df_tables_pdf_fr.at[idx_pdf, "check_outre"] = \
        df_tables_pdf_fr.at[idx_pdf, "tot_outre_calc"] == \
        df_tables_pdf_fr.at[idx_pdf, TOT_OUTREMER]
    # check Total
    df_tables_pdf_fr.at[idx_pdf, "check_tot"] = \
        df_tables_pdf_fr.at[idx_pdf, "total"] == \
        df_tables_pdf_fr.at[idx_pdf, TOT_METROPOLE] + \
        df_tables_pdf_fr.at[idx_pdf, TOT_OUTREMER] 
    # final flag OK
    df_tables_pdf_fr.at[idx_pdf, "flag_ok"] = \
        df_tables_pdf_fr.at[idx_pdf, "check_outre"] & \
        df_tables_pdf_fr.at[idx_pdf, "check_metrop"] & \
        df_tables_pdf_fr.at[idx_pdf, "check_tot"]

# add csv data to final dataframe
df_fr = update_df_fr(df_tables_pdf_fr)
    
# if __name__ == '__main__':
    # print(extract_text_from_pdf('w9.pdf'))
