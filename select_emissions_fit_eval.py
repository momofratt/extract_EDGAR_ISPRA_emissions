#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 12 12:36:55 2021

@author: cosimo
"""
from netCDF4 import Dataset
import sys
import numpy as np 
import math 
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score
from shapely.geometry import Point
import select_emissions_param as par
import select_emissions_netcdf as net

def eval_emission(spec, edgar_version, start_end_years, relative_emi_err_asymm, sector_suff=''):
    """
    evalutate total emissions for a given gas specie spec from EDGAR netcdf files.

    Parameters
    ----------
    spec : str
        specie to be analysed. It is used to define input and output filenames.
    edgar_version : str
        version of the edgar files. It is used to read input filenames
    start_end_years : list
        list with two elements [start_year,end_year].
    relative_emi_err_asymm: list
        list with lower and upper relative errors
    sector_suff : str, optional
        suffix to analyse sector specific netcdf files from edgar inventory. The default is ''.
        
    Returns
    -------
    emissions: numpy.array:
        total emissions for each year
    tot_error_inf: numpy.array: 
        lower error for each year
    tot_error_sup: numpy.array: 
         upper error for each year
    list:
        years corresponding to the emissions
            

    """
    
    if sector_suff == '': # suffix for outfile name
        sector_suffix = ''
    else:
        sector_suffix = '_'+sector_suff
    
    emission = []
    emission_err = []
    
    start_year, end_year = start_end_years
    years = range(start_year, end_year+1)

    for year in years: # loop over years
        print('evaluate emissions for', year)
        if sector_suff =='':  # define filenames for input netcdf EDGAR data. 
            file_path = './data_'+spec.upper()+'/'+edgar_version+'_'+spec.upper()+'_'+str(year)+'_TOTALS.0.1x0.1/'
            file_name = edgar_version+'_'+spec.upper()+'_'+str(year)+'_TOTALS.0.1x0.1.nc'
        else: # If the sector suffix is specified, reads the emissions for the given sector
            file_path = './data_'+spec.upper()+'/sectors/'+edgar_version+'_'+spec.upper()+'_'+str(year)+'_'+sector_suff+'.0.1x0.1/'
            file_name = edgar_version+'_'+spec.upper()+'_'+str(year)+'_'+sector_suff+'.0.1x0.1.nc'
                
        data = Dataset(file_path+file_name, 'r') # read inpput netcdf
        
        # define matrices to store selected data
        sel_data     = np.zeros((par.nbin,par.nbin)) # array to store data inside the boundary 
        sel_data_int = np.zeros((par.nbin,par.nbin)) # array to store data inside the internal boundary 
        sel_data_ext = np.zeros((par.nbin,par.nbin)) # array to store data inside the external boundary 
       
        # output file with three columns (lat, lon, emi_ch4)
        outfile = open('./res_files/selected_emissions/selected_emissions'+sector_suff+'_'+par.region+'_'+spec+'_'+str(year)+'.txt', 'w')
        outfile.write('Lat Lon emi_'+spec+'[kg/m2/s]\n')
        
        for j in range(0, par.nbin): # loop over data matrix

            for i in range(0,par.nbin):
                lat_ind = par.j_ref - par.nbin_2 + j
                lon_ind = par.i_ref - par.nbin_2 + i
                point_coords = net.index_to_lat( lat_ind, lon_ind)
                point = Point(point_coords[1], point_coords[0])
               
                if point.within(par.regione): # if the point is included inside the sudy area, then it is added to the sel_data array
                    sel_data[j][i] = (( data['emi_'+spec.lower()][ lat_ind ][ lon_ind ] ))
                    outfile.write(str(point_coords[0]) + " " + str(point_coords[1]) + " " + str(sel_data[j][i]) + "\n") # write on file
                
                if net.within_internal_boundary(point, par.regione): # if the point is included inside the inner boundary area, then it is added to the sel_data_int array 
                    sel_data_int[j][i] = (( data['emi_'+spec.lower()][ lat_ind ][ lon_ind ] ))
                
                if net.within_external_boundary(point, par.regione): # same as before but for external boundary
                    sel_data_ext[j][i] = (( data['emi_'+spec.lower()][ lat_ind ][ lon_ind ] ))

        outfile.close()
        
        # ###############################################################################
        ###                      Evaluate total emission                            ###
        ###############################################################################
        # initalise to 0 total emissions for the study area, inner and external boundaries
        tot_emi = 0
        tot_emi_int = 0
        tot_emi_ext = 0
        earth_rad2 = pow(par.earth_rad*1000,2) # square earth radius [m^2]
        d_lat_rad = par.d_lat *math.pi/180     # convert grid point separation lat and lon to radiants
        d_lon_rad = par.d_lon *math.pi/180
        pi = math.pi
       
        for j in range(0, par.nbin):    # loop over latitudes to evaluate total emission
            temp_lat, _ = net.index_to_lat(par.j_ref-par.nbin_2 + j , 1) # get latitude of jth iteration
            bin_surf = earth_rad2 * d_lat_rad * math.sin( 0.5*pi - abs(temp_lat) *pi/180) * d_lon_rad # [m^2] evaluate the surface of the grid point. dS = r*Dtheta * r*sin(theta)*Dphi
            for i in range(1, par.nbin-1):      # loop over longitudes
                tot_emi     = tot_emi     + sel_data[j][i]*bin_surf  # add the grid point emission to total emission
                tot_emi_int = tot_emi_int + sel_data_int[j][i]*bin_surf
                tot_emi_ext = tot_emi_ext + sel_data_ext[j][i]*bin_surf

        emission.append(tot_emi*60*60*24*365/1000)
        emission_err.append( max(abs(tot_emi-tot_emi_int), abs(tot_emi-tot_emi_ext))*60*60*24*365/1000 )
   
    # evaluate total error by adding the spatial error and the relative error from AD and EFs
    tot_error_inf, tot_error_sup = eval_tot_error_asym(emission, emission_err, relative_emi_err_asymm) 
    
    #open total emission file and write header
    outfile_emi = open('./res_files/total_emission'+sector_suffix+'_'+par.region+'_'+spec+'_'+str(start_year)+'-'+str(end_year)+'.txt', 'w')
    outfile_emi.write('year tot_emi[t] inf_err[t] sup_err[t]\n')

    for i in range(len(emission)):
        outfile_emi.write(str(years[i]) + ' ' + str(round(emission[i])) + ' ' + str(round(tot_error_inf[i])) + ' ' +str(round(tot_error_sup[i]))  +'\n')

    outfile_emi.close() # close total emission file

    return np.array(emission), np.array(tot_error_inf), np.array(tot_error_sup), [y for y in range(start_year, end_year+1)]


def eval_tot_error_asym(emi, emi_err_spatial, relative_emi_err_asymm):
    """
    evaluate the total error by adding the error due to the domain selection with the error relative to Emission Factors (EFs) and Activity Data (AD)

    Parameters
    ----------
    emi : list
        list containing emission values.
    emi_err_spatial : list
        list containing emission error due to domain selection.
    relative_emi_err_asymm : list
        2 elements list containing the lower and the upper relative error due to EFs and AD.

    Returns
    -------
    new_error_inf : list
        list conatining the lower absolute error for each element of the input emi list.
    new_error_sup : list
        list conatining the upper absolute error for each element of the input emi list.

    """

    relative_err_inf, relative_err_sup = relative_emi_err_asymm

    new_error_sup = []
    new_error_inf = []
    
    if len(emi) != len(emi_err_spatial):
        print('ERROR in eval_tot_error_asym function: emi list and emi_err_spatial list length do not correspond')
        sys.exit()
   
    else:  
        for i in range(len(emi)): # add error due to spatial boundary with the relative error from the Emission Factors and Activity Data
            new_error_sup.append( relative_err_sup*emi[i] + emi_err_spatial[i] )
            new_error_inf.append( relative_err_inf*emi[i] + emi_err_spatial[i] )

    return new_error_inf, new_error_sup



