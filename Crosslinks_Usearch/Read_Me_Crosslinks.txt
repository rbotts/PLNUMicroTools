Crosslinks.py

Tool used to visually compare plasmids and draw links between genes that share similar sequences
Implemented usearch so that the script can be run locally

Reads in two gen bank format files to create the graphical output from
	all the features of the genbank first are put into separate fasta files to be compared by usearch
	then the fasta files are compared and an output tab file is created and converted into csv
	then crosslinks generates the image
	
	** When running select yes for highlighting a specific gene and name "C" to highlight all concatenation sites red **
	If no gene is selected to be highlighted the genes are colored randomly and not based on the gene type
	
	Currently there are errors when forming the links based from the usearch output
	## I have not identified the error or why they are not forming ##
	Links formed in lines 219-234 for debugging
	
The old script requiring manual BLAST on ncbi is included in the directory Old_Blast_Crosslinks
	This requires the user to manually create the fasta files from the genbank through the toFasta.py script
	then the created files must be run through BLAST and the user must download the output as a .csv file
	the user is prompted through this script and is asked for the names of the different input as the script runs.
	## This has no issues with forming links but requires multiple manual steps ##