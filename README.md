Python code for the evaluation of total emissions of chemica species from the EDGAR netcdf emissions database (https://edgar.jrc.ec.europa.eu/emissions_data_and_maps). 
The spatial selection is performed over provinces that are selected from the naturalhearth admin1.geojson file (https://www.naturalearthdata.com/downloads/10m-cultural-vectors/).

Only for Italy, also the evaluation of total emissions from the ISPRA database (http://emissioni.sina.isprambiente.it/serie-storiche-emissioni/ , http://emissioni.sina.isprambiente.it/wp-content/uploads/2021/11/DB_ON_LINE_CampiIncrociati.xlsx) can be performed.

The code need to access to the admin1.geojson file from the naturalhearth database as well as to the netcdf emission database from EDGAR. Both the geojson file and the EDGAR databases must be stored locally.

The code is divided in a main program (select_emissions.py) and other modules that contain specific functions for the data elaboration.

The file **select_emissions_param.py** is used to specify the region over wich the evaluation of total emissions must be performed, as well as other parameters.
The provinces over which to perform the evaluation are defined using a list:
    prov_numname = [(3, 'Novara'),(12, 'Varese'),(13, 'Como')]
    For each province a touple with a number and a string is required. The number corresponds to the province number in the ISPRA database (thus the province number is required only for italian provinces in the case of the evaluation of ISPRA emissions), while the string corresponds to the province name on the admin1.geojson database. The case reported above will then select emissions from the provinces of Novara, Varese and Como. For non-italian provinces a random number can be used as province number. NB1 in case a province is divided in more parts (e.g. some islands are included in the province boundary) only the biggest part of the province is retained. NB2 provinces must be adjacent (e.g the evaluation of the emissions from the provinces of Paris and Beijing is not possible).
The evaluation of ISPRA emission is performed by using the IPR variable. If True, the evaluation of ISPRA emissions is performed, if False not.
Other parameters in the select_emissions_param.py file are used only to define filenames and plot titles.

The script can perform different tasks: 
1) evaluate yearly total emissions over some provinces with their relative error. The error is defined as the sum of the error related to the spatial domain selection and a relative error that is provided by the user. The error related to the spatial domain selection is due to the fact that the 0.1x0.1 bins from the EDGAR netcdf file can not follow precisely the provinces boundaries.
An algorithm for the spatial domain error evaluation was defined: the error is evaluated by selecting an inner and an outer boundaries and by evaluating the differences in the emission estimation for the "standard" case and the inner/outer ones. The error that is obtained with this routine is than compared with the EDGAR relative error and the larger error is retained. 
The error provided by the user is passed to the function select_emissions_fit_eval.eval_emission() using the argument relative_emi_err_asymm. This is a list on two elements containing the lower and upper relative errors, expressed in decimal (i.e. [0.15, 0.23] for a lower and upper error of 15% and 23% respectively).
The selection can be performed also on single sectors from EDGAR or ISPRA. For references to EDGAR and ISPRA (https://www.isprambiente.gov.it/files2022/pubblicazioni/rapporti/dis_inv_naz_12luglio_2022_rev-1.pdf) sectors refer to their websites.                                                                                                                                                                                                                            2) results from emission selection are saved on files. Each file contains mean falue, lower and upper error for each year. In the res_files/selecte_emissions are stored files that contains the flus value from EDGAR for each selected bin with its relative latitude and longitude
3) plot total emissions from EDGAR and ISPRA
4) plot 2D map showing 2D emissions from EDGAR and the used boundary for emission selection.
                                                                                                                                                                                                                                        
                
NOTE:
- input netcdf file must have 0.1x0.1 degree resolution
- data for ISPRA are stored in the data_ISPRA folder
- data for EDGAR are stored in the data_SPECIE folder, where specie must be replaced with the proper specie
- data for EDGAR sectors are stored in the data_SPECIE/sectors folder
- EDGAR data must have yearly resolution (i.e. one file for each year)
