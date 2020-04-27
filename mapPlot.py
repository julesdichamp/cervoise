#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 11:44:44 2020

@author: florian
"""
from bokeh.models import ColumnDataSource, GMapOptions, GeoJSONDataSource
from pyproj import Proj, transform
import geopandas as gpd
from bokeh.plotting import figure, output_file, show
from bokeh.tile_providers import CARTODBPOSITRON, get_provider
import json
from bokeh.io import show, output_notebook, output_file
from bokeh.models import GeoJSONDataSource, HoverTool, LinearColorMapper
from bokeh.plotting import figure
from bokeh.palettes import Viridis6


import numpy as np, pandas as pd, cartopy.crs as ccrs
import shapefile as shp
import matplotlib.pyplot as plt
import seaborn as sns
import folium
from folium import plugins
from branca.colormap import linear
from folium.plugins import TimeSliderChoropleth
sns.set(style="whitegrid", palette="pastel", color_codes=True)

sns.mpl.rc("figure", figsize=(10, 6))

regContour ='/data/regions.geojson'
# Read the geojson file
map_df = gpd.read_file(regContour)
# Convert the IDs to int
map_df.ID = map_df.code
# Get info about data
map_df.info()

nb = G.number_of_nodes()
reg = pd.DataFrame(columns=["long", "lat", "name", "pop", 'death', 'NCCENR'])

timeList = list(G.nodes[i]['data_covid'].keys())

for i in range(nb):
    print(G.nodes[i]['names'])
    if G.nodes[i]['names'] == 'Provence-Alpes-Côte d’Azur':
        name = "Provence-Alpes-Côte d'Azur"
    elif G.nodes[i]['names'] == 'Ile-de-France':
        name = 'Île-de-France'
    else:
        name =G.nodes[i]['names']
    reg.loc[i] = [
        G.nodes[i]["positions"][0],
        G.nodes[i]["positions"][1],
        i,
        G.nodes[i]["total_pop"],
        int(G.nodes[i]['data_covid']['2020-03-22']['Deaths']),
        name,
    ]

#merging geo with data
map_df.set_index('nom', inplace=True)
reg.set_index('NCCENR', inplace=True)
merged = map_df.join(reg)
merged.set_index('name',inplace=True)
merged[[ 'code','geometry','pop','death']].to_file('prout.json', driver='GeoJSON')

styledata = {}

for country in reg.name:
    d = []
    for t in timeList:
        d.append(G.nodes[country]['data_covid'][t]['Deaths'])
    df = pd.DataFrame(
        {'color': d,
         'opacity': d},
        index=timeList
    )
    # df = df.cumsum()
    # df.sample(n_sample, replace=False).sort_index()
    # styledata[reg.iloc[country].NCCENR] = df
    styledata[country] = df

max_color, min_color, max_opacity, min_opacity = 0, 0, 0, 0

for country, data in styledata.items():
    max_color = max(max_color, data['color'].max())
    min_color = min(max_color, data['color'].min())
    max_opacity = max(max_color, data['opacity'].max())
    max_opacity = min(max_color, data['opacity'].max())


cmap = linear.Greys_09.scale(min_color, max_color)


def norm(x):
    return (x - x.min()) / (x.max() - x.min())


for country, data in styledata.items():
    data['color'] = data['color'].apply(cmap)
    data['opacity'] = norm(data['opacity'])*0.5 +0.4


styledict = {
    str(country): data.to_dict(orient='index') for
    country, data in styledata.items()
}

mapa = folium.Map(location=[48, 2],tiles='Stamen Toner', zoom_start=7)
g = TimeSliderChoropleth(
    merged.to_json(),
    styledict=styledict,
).add_to(mapa)

mapa.choropleth(
      geo_data='prout.json', # geoJSON file
      name='choropleth',
      data=merged, # Pandas dataframe
      columns=['code','death'], # key and value of interest from the dataframe
      key_on='feature.properties.code', # key to link the json file and the dataframe
      fill_color='Greys', # colormap
      fill_opacity=0.0,
      line_opacity=0.0,
      smooth_factor=10,
      legend_name='Death'
)

mapa.save("index.html")