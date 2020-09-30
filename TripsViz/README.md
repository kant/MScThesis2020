# Trips-Viz Readme

This is essentially a clone of github.com/skiniry/Trips-Viz.

There are however a number of changes. The scripts required for the quantification of the transcripts of a selected gene are included along with a number of changes:

Recorded sqlite_path_organism as a string in single_transcript_routes.py and comparison.py at the if owner == 1 statement. 

Introduced the counts function from TripsCount to get the scores per transcript. 
 
Return a string in the single_transcript_routes.py scripts Query() function now including the vcalculated coverage scores. 

Made change in JavaScript common.js to replace type with coverage in table column headers of the transcript table. 

# --------------------------------------------- FROM TRIPS-VIZ README GITHUB -------------------------------------------------------------------



A transcriptome browser for Ribo-Seq and RNA-Seq data

Trips-viz is a transcriptome browser designed to visualize Ribosome profiling and RNA-seq data at the level of a single gene/transcript isoform as opposed to at the genome level. Trips-viz also provides you the ability to vizualize data from a gene under different conditions, get meta-information at an individual dataset level such as read length distribution or triplet periodicity, and provides the ability to find differentially expressed or translated genes.

To run Trips-Viz locally rename trips_empty_db.sqlite to trips.sqlite, change the value of SCRIPT_LOC in config_template.py to the full path of the Trips-Viz directory on your computer and rename config_template.py to config.py.

Run the init.py script and navigate to 0.0.0.0:5000/ in a browser (if running on a server replace 0.0.0.0 with the server IP).
