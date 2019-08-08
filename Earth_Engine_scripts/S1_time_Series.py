# -*- coding: utf-8 -*-
"""
Load the Sentinel-1 ImageCollection from the Earth engine
Read shapefile and select Modulo/Turno/Lote to retrieve the time series
Get the time series of a manually defined polygon
Get the time series of a manually defined point (lat,lon)

Pending:
    comment the code (specially in the library)
    Enable customization of plots (Change limits, titles, location of legends, dates x-axis ...)
    Filter image collection by geometry (e.g. Doñana, Nazca, Agrisar, etc)
    Combine incidence angles for each polarization
    Visualice maps
    
@author: Cristian Silva
"""
import ee
ee.Initialize()
#import time
from datetime import datetime
#from pprint import pprint
#import geextract
#from geextract import ts_extract, get_date
#import numpy as np
#import matplotlib.pyplot as plt
#import pandas as pd
import geopandas
import EE_library
#import matplotlib.dates as mdates
#---------------------------- Functions -----------------------------------------
# Define function to map over imageCollection to perform spatial aggregation 
#def _reduce_region(image,geometria):
    #if radius is not None or feature is not None:
    # Define spatial aggregation function
#--------------------------------------------------------------------------------------------    
save_ts = "No"
path_out_img = ""
#----------------------------Parameters: DATES, DESCENDING OR DESCENDING, POLARISATION------- 
str2='2019-08-05'
str1='2016-12-01'
Direction='ASCENDING'
#Direction='DESCENDING'
#date2 = ee.Date(str2)
date2=ee.Date(datetime.today())
date1 = ee.Date(str1)
ds = date2.difference(date1, 'day')
orbit_No=91
#---------------------------------- Load the Sentinel-1 ImageCollection -----------------------
sentinel1 = ee.ImageCollection('COPERNICUS/S1_GRD').filterDate(date1, date2).filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV')).filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH')).filter(ee.Filter.eq('instrumentMode', 'IW')).filter(ee.Filter.eq('orbitProperties_pass', Direction)).filter(ee.Filter.eq('relativeOrbitNumber_start', orbit_No)); 
#sentinel1 = ee.ImageCollection('COPERNICUS/S1_GRD').filterDate(date1, date2).filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV')).filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH')).filter(ee.Filter.eq('instrumentMode', 'IW')); 

#------------------------------------ Read shapefile and select Modulo/Turno/Lote to retrieve the time series ------------
geo_df = geopandas.read_file("D:\\PhD Info\\EO4 Cultivar\\Ground data\\20180504_Danper_Ground_data\Shapefiles\\Danper_Campositan_parcels_boundaries_lotes_EPSG 4326 - WGS 84.shp")
# MT to link with pics
#List_of_interesting_modulos = ["M14","M14","M14","M14","M24","M12"]
#List_of_interesting_turnos =  ["T02","T04","T05","T06","T02","T05"]
List_of_interesting_modulos = ["M12"]
List_of_interesting_turnos =  ["T02"]
# MT with phenology
#List_of_interesting_modulos = ["M19","M11","M24","M10","M23","M06","M22","M12","M14",'M07']
#List_of_interesting_turnos =  ["T04","T01","T02","T05","T05","T02","T03","T02","T05",'T04'] # Turno representativo de ese modulo
filter_geo_df_by="ID_MT"
stats = 'mean'    #   median, max, min, minMax
out_path=""
dicc_of_geometries1={}
# If it was not Peru but a generic case, filter the geodataframe by the column desired and use ID_MT to locate the desired geometry
for o in range(len(List_of_interesting_modulos)):
    #ID_MT="M23T05"
    ID_MT=List_of_interesting_modulos[o]+List_of_interesting_turnos[o]
    geo_df_IDMT=geo_df.loc[geo_df[filter_geo_df_by] == ID_MT]  
    # Use automatic_ax_lim="No" if it is desired to do interparcel visual comparisons. All plot will have then same y limits               
    dicc_of_geometries1[ID_MT]=EE_library.time_series_from_shp(sentinel1,geo_df_IDMT,ID_MT,stats,save_ts,out_path,automatic_ax_lim="No")
    #dicc_of_geometries1[ID_MT].to_csv(path_or_buf="C:\\Users\\crisj\\Box\\Danper_Data_analysis\\Output_xlsx\\EE_"+ID_MT+"_"+Direction+"_"+str(orbit_No))
#-------------------------------------Define manually the point/polygon to get time series for ----------------------------------------------------------
#Geometries for test
M06T01L237 = ee.Geometry.Polygon([[[-78.81792738406443, -8.491412672228353],
                                   [-78.81786032883906, -8.492288097425805],
                                   [-78.81688668696665, -8.492158110174252],
                                   [-78.8170020219543, -8.491298601525253]]])

M24T03L085 = ee.Geometry.Polygon([[[-78.79227799971744, -8.532604673332454],
                                   [-78.79161281188175, -8.533113957068709],
                                   [-78.7907759626691, -8.532052948518402],
                                   [-78.79147333701297, -8.531554273481474]]])
                                   
M07T04L020 = ee.Geometry.Polygon([[[-78.82064779049574, -8.471746481621691],
                                   [-78.82111449486433, -8.472574198772755],
                                   [-78.82030446774183, -8.47303050362005],
                                   [-78.81983239895521, -8.47214442278201]]])  

Parcela_A_Doñana=ee.Geometry.Polygon([[     [-6.109310154489094,37.05954307806896 ],
                                            [-6.100791458657795,37.0607417086733  ],
                                            [-6.102529530099446,37.06158073882878 ],
                                            [-6.110297207406575,37.060536230486036],
                                            [-6.109310154489094,37.05954307806896 ]]])                               

geometria=['M06T01L237','M24T03L085','M07T04L020'] # Iterate over this list of geometries in case more than one region time series is needed
#---------------------------------------- Time series of a region (image.reduceRegion) -------------------
sentinel1 = ee.ImageCollection('COPERNICUS/S1_GRD').filterDate(date1, date2).filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV')).filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH')).filter(ee.Filter.eq('instrumentMode', 'IW')).filter(ee.Filter.eq('orbitProperties_pass', Direction)); 
geometry=Parcela_A_Doñana #M06T01L237    #  
sentinel1=sentinel1.filterBounds(geometry)
stats = 'mean'    #   median, max, min, minMax
title='M06T01L237'
# Use automatic_ax_lim="No" if it is desired to do interparcel visual comparisons. All plot will have then same y limits               
df_region=EE_library.Time_series_of_a_region(sentinel1,geometry,stats,title,automatic_ax_lim="Yes")
#---------------------------------------------------- Time series of a point -------------------------------------------
lon= -78.79694540336544
lat=-8.519651605063213
point = ee.Geometry.Point(lon, lat)
title='Point in M06T01L237'   
df_point=EE_library.Time_series_of_a_point(sentinel1,point,title)



