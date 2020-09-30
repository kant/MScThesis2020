# TripsCountPy2 Readme 

This file contains the scripts for transcript quantification using the exisiting data on Trips-Viz. 

The tripsSplicepy2.py file is the same as that found in the TripsSplice directory but written in python 2. This was done for consistencey with the exisitng Trips-Viz codebase. This is a utility file containing various functions used in the processing of annotations and read alignments from the sqlite file. 

The tripsCountpy2.py file contains functions for further processing of the read alignment data and the implementation of the feature counts-like algorithm described in section 2.6 in the thesis. 

The orfQuant.py script contains the implementation of the ORF quant-like algorithm described in section 2.7 of the thesis. 


To run these scripts copy the data from the data directory into this directory (or alter the paths to the files in the code). This script will require the use of python 2 and the installation of sqlite and sqlitedict modules
