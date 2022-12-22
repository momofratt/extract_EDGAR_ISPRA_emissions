#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 12 13:07:11 2021

@author: cosimo
"""
from netCDF4 import Dataset
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import pandas as pd
import matplotlib.pyplot as plt 
import matplotlib.colors as colors
from shapely.geometry import Point
import numpy as np
import select_emissions_param as par
import select_emissions_netcdf as net
import sys


def plot_emissions(specie, years_edgar, emis, emis_err, anni_ISPRA='', emis_ISPRA='', emis_ISPRA_err=[]):

    IPR_suffix = '' # suffisso per nome plot
    
    fig, ax = plt.subplots(1,1, figsize = (8,4))
    

    plt.rcParams.update({'font.size':15, 'axes.labelsize' : 15, 'xtick.labelsize' : 15, 'ytick.labelsize':15})
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)


    ax.errorbar(years_edgar,emis/1e6,[emis_err[0]/1e6,emis_err[1]/1e6],fmt='o', elinewidth=1, capsize=3, c='grey',label = 'EDGAR', barsabove=False)
    
    
    ################################################
    if (anni_ISPRA!='') & (emis_ISPRA!='') & (emis_ISPRA_err!=[]): 
        ax.errorbar(anni_ISPRA, [e*1e-6 for e in emis_ISPRA], [[e*1e-6 for e in emis_ISPRA_err[0]], [e*1e-6 for e in emis_ISPRA_err[1]]], 
                    fmt='^', elinewidth=1, capsize=3, c='black', label = 'ISPRA', barsabove=True)
        IPR_suffix='ISPRA_'
       
      
    ax.set_xlabel('years', fontsize=16)
    
    ax.set_ylim(bottom =0)
    ax.set_ylabel(specie+' emissions [Mt]', fontsize=16)
    ax.grid()
    ax.legend(loc='best')
    

    fig.savefig(par.res_folder+'total_emissions_EDGAR_'+IPR_suffix+par.region+'_'+specie+'.png', format = 'png', bbox_inches='tight', dpi=200)
    plt.close(fig)
    plt.rcdefaults()

def plot_2D_emis(spec, year,plot_included_bins, edgar_version = '', sector_suff=''): 
    # plot 2d emission and selected area for a given specie in a giver year
    dlat = 0.5 + abs((max(par.max_lat,par.min_lat)-par.lat_ref))
    dlon =   1 + abs((max(par.max_lon,par.min_lon)-par.lon_ref))

    latinf, loninf = net.lat_to_index(par.lat_ref-dlat, par.lon_ref-dlon)
    latsup, lonsup = net.lat_to_index(par.lat_ref+dlat, par.lon_ref+dlon)
    
    if sector_suff =='':  # define filenames for input netcdf EDGAR data. 
        file_path = './data_'+spec.upper()+'/'+edgar_version+'_'+spec.upper()+'_'+str(year)+'_TOTALS.0.1x0.1/'
        file_name = edgar_version+'_'+spec.upper()+'_'+str(year)+'_TOTALS.0.1x0.1.nc'
    else: # If the sector suffix is specified, reads the emissions for the given sector
        file_path = './data_'+spec.upper()+'/sectors/'+edgar_version+'_'+spec.upper()+'_'+str(year)+'_'+sector_suff+'.0.1x0.1/'
        file_name = edgar_version+'_'+spec.upper()+'_'+str(year)+'_'+sector_suff+'.0.1x0.1.nc'
    

    data = Dataset(file_path+file_name, 'r')

    ### Plot selected region ###
    lats = data['lat'][latinf:latsup]
    lons = data['lon'][loninf:lonsup]
    emis = data['emi_'+spec.lower()][latinf:latsup, loninf:lonsup] 

    # set up a map
    ax1 = plt.axes(projection=ccrs.Stereographic(central_latitude=par.lat_stat, central_longitude=par.lon_stat))

    # define the coordinate system that the grid lons and grid lats are on
    data_transf = ccrs.PlateCarree()
    plt.pcolormesh(lons, lats, emis, transform=data_transf, norm=colors.LogNorm())

    # Add Colorbar
    plt.colorbar()

    # Add Title
    plt.title(spec+' emission [kg m$^{-2}$ s$^{-1}$] for '+str(year))

    # Add borders
    border_lats = [b[1] for b in par.bound_coords]
    border_lons = [b[0] for b in par.bound_coords]
    ax1.plot(border_lons,border_lats, transform=data_transf,color='m', lw = .5)

    # Add station mark
    ax1.plot(par.lon_stat, par.lat_stat,transform=data_transf, marker='^',color='red', label=par.stat_nm, markersize=3.5)

    if plot_included_bins == True:
        for j in range(0, par.nbin): 
                for i in range(0, par.nbin):  
                      lat_ind = par.j_ref - par.nbin_2 + j
                      lon_ind = par.i_ref - par.nbin_2 + i
                      point_coords = net.index_to_lat( lat_ind, lon_ind)
                      point = Point(point_coords[1], point_coords[0])

                      if point.within(par.regione):
                          ax1.plot(point_coords[1], point_coords[0], transform=data_transf, marker='o', color='y', markersize=2)

                      if net.within_internal_boundary(point, par.regione):
                          ax1.plot(point_coords[1], point_coords[0], transform=data_transf, marker='x', color='g', markersize=2)

                      if net.within_external_boundary(point, par.regione):
                          ax1.plot(point_coords[1], point_coords[0], transform=data_transf, marker='x', color='g', markersize=1)


    # Add coastlines and rivers
    ax1.add_feature(cfeature.COASTLINE)
    ax1.add_feature(cfeature.RIVERS)

    gl = ax1.gridlines(draw_labels=True, x_inline=False, y_inline=False)
    gl.top_labels = False
    gl.right_labels = False
    
    if sector_suff!='':
        sector_suffix='_'+sector_suff
    else:
        sector_suffix=''
        
    plt.savefig('./results/'+spec+'_2D_emission_'+par.region+'_'+str(year)+ sector_suffix + '.png', format='png', dpi=200)
    plt.close('all')
