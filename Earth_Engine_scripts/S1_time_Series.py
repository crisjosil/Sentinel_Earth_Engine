# -*- coding: utf-8 -*-
"""
Created on Sun Aug  4 01:03:05 2019

@author: crisj
"""
import ee
ee.Initialize()
import time
from datetime import datetime
from pprint import pprint
import geextract
from geextract import ts_extract, get_date
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import geopandas

#Parameters: DATES, DESCENDING OR DESCENDING, POLARISATION 
str2='2019-08-04'
str1='2016-12-01'
#Direction='ASCENDING'
Direction='DESCENDING'
date2 = ee.Date(str2)
date1 = ee.Date(str1)
ds = date2.difference(date1, 'day')
#---------------- Load the Sentinel-1 ImageCollection -----------------------
sentinel1 = ee.ImageCollection('COPERNICUS/S1_GRD').filterDate(date1, date2).filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV')).filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH')).filter(ee.Filter.eq('instrumentMode', 'IW')).filter(ee.Filter.eq('orbitProperties_pass', Direction)).filter(ee.Filter.eq('relativeOrbitNumber_start', 142)); 

lon= -78.79694540336544
lat=-8.519651605063213
geometry = ee.Geometry.Point(lon, lat)
polygon = ee.Geometry.Polygon([[[-78.79737094642212,-8.520166448964439], 
                                [-78.79594669343521,-8.520017902912258], 
                                [-78.79635170699646,-8.519153151534463], 
                                [-78.79773840905716,-8.519322918830348], 
                                [-78.79737094642212,-8.520166448964439]]])

Read_Sentinel_data = sentinel1.filterBounds(geometry).getRegion(geometry, 30).getInfo()
out = geextract.dictify(Read_Sentinel_data)

def get_S1_date(out):
    """ Obtains the SAR image acquisition date from the product ID

    Args:
        out (List of dictionaries): Each dictionary corresponds to the sentinel1 info of an image

    Returns:
        Dates : A list of the acquisition dates
    
    """
    Dates=np.zeros(len(out)).tolist()
    VH=np.zeros(len(out)).tolist()
    VV=np.zeros(len(out)).tolist()
    angle=np.zeros(len(out)).tolist()
    df=pd.DataFrame()
    for i in range(len(out)):
        a=out[i]['id']
        b=a.split(sep='_')
        b1=b[4].split(sep='T')
        b2=b1[0]
        Dates[i]=datetime.strptime(b2, "%Y%m%d")
        VH[i]=out[i]['VH']
        VV[i]=out[i]['VV']
        angle[i]=out[i]['angle']
    
    
    df['Dates']=Dates
    df['VH']=VH
    df['VV']=VV
    df['angle']=angle
    df.set_index('Dates', inplace=True)
    
    return(df)
        
df=get_S1_date(out)

fig, (ax, ax1,) = plt.subplots(2, 1, sharex=True,figsize=(19,9)) 
ax.plot(df['VH'])
ax1.plot(df['VV'])

# Define function to map over imageCollection to perform spatial aggregation 
def _reduce_region(image):
    #if radius is not None or feature is not None:
    # Define spatial aggregation function
    if stats == 'mean':
        fun = ee.Reducer.mean()
    elif stats == 'median':
        fun = ee.Reducer.median()
    elif stats == 'max':
        fun = ee.Reducer.max()
    elif stats == 'min':
        fun = ee.Reducer.min()
    elif stats == 'minMax':
        fun = ee.Reducer.minMax()
    else:
        raise ValueError('Unknown spatial aggregation function. Must be one of mean, median, max, or min')

    """Spatial aggregation function for a single image and a polygon feature"""
    stat_dict = image.reduceRegion(fun, geometry, 30);
    # FEature needs to be rebuilt because the backend doesn't accept to map
    # functions that return dictionaries
    return ee.Feature(None, stat_dict)

stats="mean"
fc = sentinel1.filterBounds(polygon).map(_reduce_region).getInfo()
out1 = geextract.simplify(fc)

df1=get_S1_date(out1)
df1.sort_index()

plt.style.use('ggplot')
fig, (ax, ax1,) = plt.subplots(2, 1, sharex=True,figsize=(19,9)) 
df['VH'].plot(ax=ax,marker='o', markersize=4, linestyle='--', linewidth=1,color="blue",label="Point")
df1['VH'].plot(ax=ax1,marker='o', markersize=4, linestyle='--', linewidth=1,color="red",label="polygon")
ax.legend()
ax1.legend()
plt.tight_layout()

#print(polygon.coordinates())

M06T01L237 = ee.Geometry.Polygon([[[-78.81792738406443, -8.491412672228353],
                                   [-78.81786032883906, -8.492288097425805],
                                   [-78.81688668696665, -8.492158110174252],
                                   [-78.8170020219543, -8.491298601525253]]])

fc1 = sentinel1.filterBounds(M06T01L237).map(_reduce_region).getInfo()
out2 = geextract.simplify(fc1)

df2=get_S1_date(out2)
df2.sort_index

fig, ax = plt.subplots(1, 1, sharex=True,figsize=(19,9)) 
#df1['VH'].plot(ax=ax,marker='o', markersize=4, linestyle='--', linewidth=1,color="blue",label="M14T04L490")
df2['VH'].plot(ax=ax,marker='o', markersize=4, linestyle='--', linewidth=1,color="red",label="M06T01L237")
ax.legend()
#ax1.legend()
plt.tight_layout()

geo_df = geopandas.read_file("D:\\PhD Info\\EO4 Cultivar\\Ground data\\20180504_Danper_Ground_data\Shapefiles\\Danper_Campositan_parcels_boundaries_lotes_EPSG 4326 - WGS 84.shp")
a=geo_df['geometry'][0]

b=str(a)
c=b.split(sep='POLYGON Z ((')
print(c[1])

