#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: Cosimo Fratticioli
@contact: c.fratticioli@isac.cnr.it
"""
import read_ISPRA
import select_emissions_fit_eval as fitev
import select_emissions_param as par
import select_emissions_plot as plot

world_provinces = par.states_geoframe # used to see  and get the names of the provinces 

###############################################################################
###   Select emissions inside emission region and write output on files     ###
###############################################################################
# evaluate CH4 emissions for 2015 from sector AWB
emi_CH4, err_inf_CH4,  err_sup_CH4, years_edgar_CH4 = fitev.eval_emission('CH4', edgar_version = 'v6.0',        
                                                                          start_end_years=[2015,2015], 
                                                                          relative_emi_err_asymm=[0.17,0.23], 
                                                                          sector_suff='AWB')

# eval CH4 total emissions for years 2015-2018
emi_CH4, err_inf_CH4,  err_sup_CH4, years_edgar_CH4 = fitev.eval_emission('CH4', edgar_version = 'v6.0',        
                                                                          start_end_years=[2015,2018], 
                                                                          relative_emi_err_asymm=[0.17,0.23], 
                                                                          sector_suff='')
# eval CO total emissions for years 2015-2018
emi_CO,  err_inf_CO,   err_sup_CO,  years_edgar_CO  = fitev.eval_emission('CO',  edgar_version  = 'EDGARv6.1',  
                                                                          start_end_years=[2015,2018], 
                                                                          relative_emi_err_asymm=[0.2,0.2],   
                                                                          sector_suff='')
# eval N2O total emissions for years 2015-2018
emi_N2O, err_N2O_inf,  err_N2O_sup, years_edgar_N2O = fitev.eval_emission('N2O', edgar_version = 'v7.0_FT2021', 
                                                                          start_end_years=[2015,2018], 
                                                                          relative_emi_err_asymm=[0.2,0.2],   
                                                                          sector_suff='')


if par.IPR:
    # get ISPRA emissions and write on file
    anni_ISPRA = [1990, 1995, 2000, 2005, 2010, 2015, 2019]
    # eval CH4 emissionsfrom ISPRA
    emi_CH4_ISPRA, err_CH4_ISPRA_inf, err_CH4_ISPRA_sup = read_ISPRA.select_ISPRA(par.prov_num, anni_ISPRA, 'Metano', [0.13,0.18], sector='')
    # eval CH4 emissionsfrom ISPRA from 
    emi_CH4_ISPRA, err_CH4_ISPRA_inf, err_CH4_ISPRA_sup = read_ISPRA.select_ISPRA(par.prov_num, anni_ISPRA, 'Metano', [0.13,0.18], sector='')
    # eval CO emissions from ISPRA
    emi_CO_ISPRA, err_CO_ISPRA_inf, err_CO_ISPRA_sup    = read_ISPRA.select_ISPRA(par.prov_num, anni_ISPRA, 'Monossido di carbonio', [0.2,0.2], sector='')


###############################################################################
###                               PLOT                                      ###
###############################################################################
plot.plot_emissions('CH4', 
                    years_edgar    = years_edgar_CH4,
                    emis           = emi_CH4, 
                    emis_err       = [err_inf_CH4, err_sup_CH4] , 
                    anni_ISPRA     = anni_ISPRA,
                    emis_ISPRA     =  emi_CH4_ISPRA, 
                    emis_ISPRA_err = [err_CH4_ISPRA_inf ,err_CH4_ISPRA_sup] )

plot.plot_emissions('CO', 
                    years_edgar    = years_edgar_CO,
                    emis           = emi_CO, 
                    emis_err       = [err_inf_CO, err_sup_CO] , 
                    anni_ISPRA     = anni_ISPRA,
                    emis_ISPRA     =  emi_CO_ISPRA, 
                    emis_ISPRA_err = [err_CO_ISPRA_inf ,err_CO_ISPRA_sup] )

plot.plot_emissions('N2O', 
                    years_edgar = years_edgar_N2O,
                    emis        = emi_N2O, 
                    emis_err    = [err_N2O_inf, err_N2O_sup] )
                  
# plot 2D for CH4 emissions from sector AWB
plot.plot_2D_emis('CH4', 2015, False, edgar_version = 'v6.0',      sector_suff='AWB')

# plot 2D for CO emissions
plot.plot_2D_emis('CO' , 2014, False, edgar_version = 'EDGARv6.1', sector_suff='')

