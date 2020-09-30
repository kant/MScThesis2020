# MScThesis2020
 Towards the Development of Accurate Representations of Protein Coding Loci
 
By Jack Tierney for the degree of Masters in Bioinformatics and Computational Biology
Under the supervision of Prof. Pasha Baranov

School of Biochemistry and Cell Biology, University College Cork

The python scripts relating to the three main branches of this project are available here.

# Repository Structure

TripsSplice directory relates to the implementation of splice graphs and super transcripts using python 3

TripsCountPy2 contains the scripts used to quantify the reads across the locus using both approaches described in the thesis.

TripsViz is a clone of the official github.com/skiniry/Trips-Viz in which I have implemented the transcript quantification.


# Note

Each directory has their own README.md files with specifics on those scripts.

Please find a copy of the thesis in thesis.pdf. 

# Project Abstract

The development of accurate representations of protein coding loci has the potential to enable researchers, who are not familiar with the fallacies associated with existing annotations, to reliably interpret the regulation of gene expression at a transcript level. This project explores  alternative representations of protein coding loci from existing literature through implementation for potential use on trips.ucc.ie. The transcriptome browser, TripsViz, allows users to visualise individual transcripts alongside aligned reads (Kiniry et al., 2019). However, knowing which transcript is of most interest to the user is not clear cut. Primarily, the focus of this project is on the visual representation of protein coding loci using splice graphs. Splice graphs are commonly used in methods for transcript abundance estimation as both computational objects and visualisations (Rogers et al., 2012; Ryan et al., 2012). Here, splice graphs are produced and visualised from the available annotation and read alignment data stored on TripsViz in a manner that is consistent with existing operations carried out on the web application.

Following the implementation of these visualisations, quantification methods were also implemented as an alternative way to inform the TripsViz user of the transcript of interest. Transcript quantification is a key aspect of most RNA-Seq and Ribo-Seq analysis pipelines and as a result a well informed decision could be made on which algorithm to adapt in this instance. Due to the nature of the data stored on TripsViz, the FeatureCounts (Liao et al., 2014) and ORFquant (Calviello et al., 2020) algorithms were adapted in this case. 
