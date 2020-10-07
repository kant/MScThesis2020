# TripsSplice Readme

Language: Python 
Version: 3.7.6

Matplotlib version: 3.1


This directory contains the scripts for processing the data stored on Trips-Viz and the production of my implementations of the splice graph and the super transcript concepts described in the thesis. 

The TripsSplice.py script carries out the processing of the annotation and read alignment data stored in sqlite files and sqlitedicts on Trips-Viz. 

The plotSpliceGraph_genomicPositioning.py script will produce a splice graph where the layout of the nodes in the visualisation is determined by the genomic coordinates of the start bases of each exon. Nodes are coloured based on support by aligned reads and edges are numbered similarly. This file requires matplotlib and networkx packages. 

The plotSpliceGraph_kamadalayout.py script will produce a splice graph similar to described above but using a built in layout method. This is included to demonstrate the need for the genomic positioning layout that I devised. 

The plotSuperTranscript.py script will implement (as best possible) the concept of a supertranscript. Although not useful as a visualisation, this script was informative regarding viable visualisations in the future. From this I concluded that it is more difficult to visualise detail in a linear plot such as this. This require the matplotlib package only. 
