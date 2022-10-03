# -*- coding: utf-8 -*- adding some text here
"""
Created on Sun Aug  4 17:00:53 2019

@author: crisj
"""
import ee
ee.Initialize()
#import time
from datetime import datetime
#from pprint import pprint
#import geextract
#from geextract import ts_extract, get_date
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
#import geopandas
#import folium
import matplotlib.dates as mdates
import math


def simplify(fc):
    """Take a feature collection, as returned by mapping a reducer to a ImageCollection,
        and reshape it into a simpler list of dictionaries

    Args:
        fc (dict): Dictionary representation of a feature collection, as returned
            by mapping a reducer to an ImageCollection

    Returns:
        list: A list of dictionaries.

    Examples:
        >>> fc = {u'columns': {},
        ...       u'features': [{u'geometry': None,
        ...                      u'id': u'LC81970292013106',
        ...                      u'properties': {u'B1': 651.8054424353023,
        ...                                      u'B2': 676.6018246419446},
        ...                      u'type': u'Feature'},
        ...                     {u'geometry': None,
        ...                      u'id': u'LC81970292013122',
        ...                      u'properties': {u'B1': 176.99323997958842,
        ...                                      u'B2': 235.83196553144882},
        ...                      u'type': u'Feature'}]}
        >>> simplify(fc)
    """
    def feature2dict(f):
        id = f['id']
        out = f['properties']
        out.update(id=id)
        return out
    out = [feature2dict(x) for x in fc['features']]
    return out

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

def set_plot_ylims(automatic_ax_lim,df,polarisation):
#-------------------------- if user selects automatic axis, get the max y min for each polarisation as y limits in the plot ------------
    a=df.filter(like=polarisation)
    
    if automatic_ax_lim=="Yes":                                 # user selection
        bb=a.max()                                              # Max of VH
        if math.isfinite(bb.max()) == True:                     # Check that is neither NaN nor inf
            max_y_lim=bb.max()
        else: 
            max_y_lim=0 
                                              # Default value in case is NaN or inf
        dd=a.min()    
        if math.isfinite(dd.min())==True:                       # Repeat for minimun y lim
            min_y_lim=dd.min()
        else: 
            min_y_lim=-30  
    else:
        if polarisation == "VH":
            max_y_lim=-13
            min_y_lim=-29
        elif polarisation == "VV":
            max_y_lim=-9
            min_y_lim=-15
        elif polarisation == "Ratio":
            max_y_lim=-2.5
            min_y_lim=-13
    return(min_y_lim,max_y_lim)    
    
def Time_series_of_a_point(Img_collection,point,title):
    l = Img_collection.filterBounds(point).getRegion(point, 30).getInfo()
    out = [dict(zip(l[0], values)) for values in l[1:]]
    df=get_S1_date(out)
    df=df.sort_index(axis = 0)
    
#    fig, (ax, ax1,ax2) = plt.subplots(3, 1, sharex=True,figsize=(19,9)) 
#    ax.plot(df['VH'],marker='o', markersize=4, linestyle='--', linewidth=1,color="blue",label='VH '+title)
#    ax1.plot(df['VV'],marker='o', markersize=4, linestyle='--', linewidth=1,color="red",label='VV '+title)
#    ax2.plot((df['VH']-df['VV']),marker='o', markersize=4, linestyle='--', linewidth=1,color="green",label='Ratio '+title)
#    ax.legend()
#    ax1.legend()
#    ax2.legend()
#    ax.set_ylim(-29,-2)
#    ax1.set_ylim(-15,-2)
#    ax2.set_ylim(-13,-2)
#    ejes = fig.axes #which is used to extract the axes
#    # For every axis, set the x and y major locator
#    for axi in ejes:
#        axi.xaxis.set_major_locator(mdates.WeekdayLocator(byweekday=1, interval=2, tz=None))
#    plt.setp(ax2.get_xticklabels(), rotation=30, horizontalalignment='right',)
#    ax2.tick_params(axis='x', which='major', labelsize=9)
#    plt.tight_layout()
    
    return(df)
    
def Time_series_of_a_region(Img_collection,geometry,stats,title,automatic_ax_lim):
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

    def _reduce_region(image):
                    """Spatial aggregation function for a single image and a polygon feature"""
                    #the reduction is specified by providing the reducer (ee.Reducer.mean()), the geometry  (a polygon coord), at the scale (30 meters)
                    stat_dict = image.reduceRegion(fun, geometry, 30);
                    # FEature needs to be rebuilt because the backend doesn't accept to map
                    # functions that return dictionaries
                    return ee.Feature(None, stat_dict)
    fc = Img_collection.filterBounds(geometry).map(_reduce_region).getInfo()
    out = simplify(fc)
    df=get_S1_date(out)
    df=df.sort_index(axis = 0)
    
#    fig, (ax, ax1,ax2) = plt.subplots(3, 1, sharex=True,figsize=(19,9)) 
#    ax.plot(df['VH'],marker='o', markersize=4, linestyle='--', linewidth=1,color="blue",label='VH '+title)
#    ax1.plot(df['VV'],marker='o', markersize=4, linestyle='--', linewidth=1,color="red",label='VV '+title)
#    ax2.plot((df['VH']-df['VV']),marker='o', markersize=4, linestyle='--', linewidth=1,color="green",label='Ratio '+title)
#    ax.legend()
#    ax1.legend()
#    ax2.legend()
#
#    #-------------------------- if user selects automatic axis, get the max y min for each polarisation as y limits in the plot ------------
#    polarisation="VH"
#    min_y_lim,max_y_lim=set_plot_ylims(automatic_ax_lim,df,polarisation)
#    ax.set(ylim=(min_y_lim, max_y_lim))                                  
#    polarisation="VV"
#    min_y_lim,max_y_lim=set_plot_ylims(automatic_ax_lim,df,polarisation)
#    ax1.set(ylim=(min_y_lim, max_y_lim))
#    polarisation="Ratio"
#    min_y_lim,max_y_lim=set_plot_ylims(automatic_ax_lim,df,polarisation)
#    ax2.set(ylim=(min_y_lim, max_y_lim))
#
#    ejes = fig.axes #which is used to extract the axes
#    # For every axis, set the x and y major locator
#    for axi in ejes:
#        axi.xaxis.set_major_locator(mdates.WeekdayLocator(byweekday=1, interval=2, tz=None))
#    plt.setp(ax2.get_xticklabels(), rotation=30, horizontalalignment='right',)
#    ax2.tick_params(axis='x', which='major', labelsize=9)
#    plt.tight_layout()
    
    return(df)    
    
def time_series_from_shp(Img_collection,geo_df_IDMT,ID_MT,stats,save_ts,out_path,automatic_ax_lim):
    # Fancy plot style
    plt.style.use('ggplot')
    df1=pd.DataFrame()
    for k in range(geo_df_IDMT.shape[0]):
        x, y = geo_df_IDMT['geometry'][geo_df_IDMT.index[k]].exterior.coords.xy
        
        e_full=[]
        e_full_1=[]
        for h in range(len(geo_df_IDMT['geometry'][geo_df_IDMT.index[k]].exterior.coords)):
            e=[]
            e.append(x[h])
            e.append(y[h])
            e_full.append(e)        
        #print(e_full)
        e_full_1.append(e_full)
        geom = ee.Geometry.Polygon(e_full_1)
        geometry= geom   
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
        
        def _reduce_region(image):
                    """Spatial aggregation function for a single image and a polygon feature"""
                    #the reduction is specified by providing the reducer (ee.Reducer.mean()), the geometry  (ee.Geometry.Polygon[]), at the scale (30 meters)
                    stat_dict = image.reduceRegion(fun, geometry, 30);
                    # FEature needs to be rebuilt because the backend doesn't accept to map
                    # functions that return dictionaries
                    return ee.Feature(None, stat_dict)
        fc = Img_collection.filterBounds(geometry).map(_reduce_region).getInfo()
        out = simplify(fc)
        df=get_S1_date(out)
        df=df.sort_index(axis = 0)
        df['Ratio']=df['VH']-df['VV']
        new_col_name_VH='VH_'+str(geo_df_IDMT['ID_MTL'][geo_df_IDMT.index[k]])
        new_col_name_VV='VV_'+str(geo_df_IDMT['ID_MTL'][geo_df_IDMT.index[k]])
        new_col_name_ratio='Ratio_'+str(geo_df_IDMT['ID_MTL'][geo_df_IDMT.index[k]])
        df.rename(columns={'VH':new_col_name_VH,'VV':new_col_name_VV,'Ratio':new_col_name_ratio},inplace=True)
        df1=pd.concat([df1, df], axis=1) 
        
    # --------------------------------------- plot each lote individually --------------------------------------    
    #    fig, (ax, ax1,ax2) = plt.subplots(3, 1, sharex=True,figsize=(19,9)) 
    #    ax.plot(df[new_col_name_VH],marker='o', markersize=4, linestyle='--', linewidth=1,color="blue",label=new_col_name_VH)
    #    ax1.plot(df[new_col_name_VV],marker='o', markersize=4, linestyle='--', linewidth=1,color="red",label=new_col_name_VV)
    #    ax2.plot((df[new_col_name_VH]-df[new_col_name_VV]),marker='o', markersize=4, linestyle='--', linewidth=1,color="green",label="Ratio "+str(geo_df_IDMT['ID_MTL'][geo_df_IDMT.index[k]]))
    #    ax.legend()
    #    ax1.legend()
    #    ax2.legend()
    #    ax.set_ylim(-28.5,-13)
    #    ax1.set_ylim(-15,-9)
    #    ax2.set_ylim(-13,-2.5)
    #    ejes = fig.axes #which is used to extract the axes
    #    # For every axis, set the x and y major locator
    #    for axi in ejes:
    #        axi.xaxis.set_major_locator(mdates.WeekdayLocator(byweekday=1, interval=2, tz=None))
    #    plt.setp(ax2.get_xticklabels(), rotation=30, horizontalalignment='right',)
    #    ax2.tick_params(axis='x', which='major', labelsize=9)
    #    plt.tight_layout()
    a=df1.filter(like="VH")
    b=df1.filter(like="VV")
    c=df1.filter(like="Ratio")
#    fig, (ax, ax1, ax2) = plt.subplots(3, 1, sharex=True,figsize=(19,9))                        
#    # Plot three dataframes in the same chart, assigning the same colors to the columns
#    a.plot(ax=ax,  marker='o', figsize=(19, 9), markersize=4, linestyle='--', linewidth=1,color=["blue","red","green","black","magenta","yellow","orange","brown","pink"] )
#    a.mean(axis=1).plot(ax=ax, marker='s',figsize=(19, 9), markersize=5, linestyle='-', linewidth=3,color="black",label=ID_MT+"Avg")
#    b.plot(ax=ax1,marker='s',                  markersize=4, linestyle='--', linewidth=1,color=["blue","red","green","black","magenta","yellow","orange","brown","pink"] )
#    b.mean(axis=1).plot(ax=ax1, marker='s',figsize=(19, 9), markersize=5, linestyle='-', linewidth=3,color="black",label=ID_MT+"Avg")
#    c.plot(ax=ax2,marker='*',                  markersize=4, linestyle='--', linewidth=1,color=["blue","red","green","black","magenta","yellow","orange","brown","pink"] )
#    c.mean(axis=1).plot(ax=ax2,marker='s',figsize=(19, 9), markersize=5, linestyle='-', linewidth=3,color="black",label=ID_MT+"Avg")#               
#
#    ax.set_ylabel('VH (dB)')
#    ax1.set_ylabel('VV (dB)')
#    ax2.set_ylabel('VH-VV')
#    
##-------------------------- if user selects automatic axis, get the max y min for each polarisation as y limits in the plot ------------
#    if automatic_ax_lim=="Yes":                                 # user selection
#        bb=a.max()                                              # Max of VH
#        if math.isfinite(bb.max()) == True:                     # Check that is neither NaN nor inf
#            cc=bb.max()
#        else: 
#            cc=-5                                               # Default value in case is NaN or inf
#        dd=a.min()    
#        if math.isfinite(dd.min())==True:                       # Repeat for minimun y lim
#            eee=dd.min()
#        else: 
#            eee=-30  
#        
#        ax.set(ylim=(eee, cc))                                  
#        
#        bb=b.max()                                              # Repeat for VH channel
#        if math.isfinite(bb.max()) == True:
#            cc=bb.max()
#        else: 
#            cc=-6
#        dd=b.min()    
#        if math.isfinite(dd.min())==True:
#            eee=dd.min()
#        else: 
#            eee=-20          
#
#        ax1.set(ylim=(eee, cc))
#        
#        bb=c.max()
#        if math.isfinite(bb.max()) == True:
#            cc=bb.max()
#        else: 
#            cc=0
#        dd=b.min()    
#        if math.isfinite(dd.min())==True:
#            eee=dd.min()
#        else: 
#            eee=-17
#        ax2.set(ylim=(eee, cc))
#    else:
#        ax.set_ylim(-29,-13)
#        ax1.set_ylim(-17,-10)
#        ax2.set_ylim(-13,-2.5)        
#    ax.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
#    ax1.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
#    ax2.legend(loc='center left', bbox_to_anchor=(1.0, 0.5)) 
#    
#    ejes = fig.axes #which is used to extract the axes
#    # For every axis, set the x and y major locator
#    for axi in ejes:
#        axi.xaxis.set_major_locator(mdates.WeekdayLocator(byweekday=1, interval=2, tz=None))
#    plt.setp(ax2.get_xticklabels(), rotation=30, horizontalalignment='right',)
#    ax2.tick_params(axis='x', which='major', labelsize=9)
#     
#    plt.suptitle(ID_MT,size=16)
#    plt.tight_layout()
#    fig.subplots_adjust(top=0.95)
#    
#    if save_ts =="Yes":
#        print("saving "+ID_MT)
#        Name=ID_MT
#        plt.savefig((out_path+Name),bbox_inches='tight')
        
    return(df1)