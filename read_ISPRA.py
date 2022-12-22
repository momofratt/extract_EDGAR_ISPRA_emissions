#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 28 14:41:43 2021

@author: Cosimo Fratticioli 
@contact: c.fratticioli@isac.cnr.it
"""
# la selezione può essere effettuata anche utilizzando le province ed i rispettivi codici:
# province = [('Piacenza',33), ('Parma',34),('Reggio Emilia',35),('Modena',36),('Bologna',37),('Ferrara',38),('Ravenna',39),('Forli',40),('Rimini',99)]
import pandas as pd
import matplotlib.pyplot as plt
import select_emissions_param as par
import numpy as np
import select_emissions_fit_eval as fitev
from os import sys


def select_prov_data(prov_number):
    # seleziona dati riferiti alle province di interesse
    sel_data_frame = pd.DataFrame()
    frame = pd.read_csv( './data_ISPRA/DB_ON_LINE_CampiIncrociati.csv', sep = ' ', thousands=',')    # legge frame 

    for prov in prov_number:
        prov_frame = frame[frame['COD_PROV'] == prov]  # seleziona dati regionali per il macrosettore
        sel_data_frame = pd.concat([sel_data_frame,prov_frame], ignore_index=True)

    return sel_data_frame


def select_ISPRA(prov_num, anni, specie, inf_sup_err, sector=''):
    # returns ISPRA total emissions for the given specie
    inf_err, sup_err = inf_sup_err
    regional_frame = select_prov_data(prov_num)

    regional_frame = regional_frame[regional_frame['NOMPOL']==specie] # select emissions for the specie
    if regional_frame.empty:
        print('ERRORE: specie non conosciuta')
        sys.exit()
        
    if sector!='':
        regional_frame = regional_frame[regional_frame['SNAP'].astype(str).apply(lambda x: x.startswith(sector))]  # leggi solo emissioni dal settore selezionato

    else:
        regional_frame = regional_frame[regional_frame['SNAP'].astype(str).apply(lambda x: not x.startswith('11'))]  # elimina settore LULUCF

    emi_annuale  = []

    for yr in anni:
        emi_annuale.append( regional_frame[ ' '+str(yr)+' '].sum())

    emi_annuale_err_sup = [emi_annuale[i]*inf_err for i in range(len(emi_annuale))]
    emi_annuale_err_inf = [emi_annuale[i]*sup_err for i in range(len(emi_annuale))]



    if sector == '': # suffix for outfile name
        sector_suffix = ''
    else:
        sector_suffix = '_'+sector

    file_emi_ISPRA = open('./res_files/total_emissions'+sector_suffix+'_'+par.region+'_CH4_ISPRA.txt', 'w')
    file_emi_ISPRA.write('year tot_emi[t] inf_err[t] sup_err[t]\n')
    
    for i in range(len(emi_annuale)):
        file_emi_ISPRA.write(str(anni[i]) +' '+ str(round(emi_annuale[i])) +' '+str(round(emi_annuale_err_inf[i]))+' '+str(round(emi_annuale_err_sup[i]))+'\n')
    
    file_emi_ISPRA.close()

    return  emi_annuale, emi_annuale_err_inf, emi_annuale_err_sup

