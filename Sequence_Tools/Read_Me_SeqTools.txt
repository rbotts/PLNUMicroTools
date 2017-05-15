SequenceTools.py

Workflow:

1. GeneMarkS: the assembled contigs must be run through GeneMarkS and saved as a .fnn
file which will contain all of the open reading frames found in the sequence.

2. Usearch_Features(seq,database): this function will create annotations based on a
local database. Seq is the name of the sample you are annotation with no extensions, 
and database is the name of the usearch database you want to create annotations from. 
After annotating locally all ORF’s with no hits will be aligned to BLAST on NCBI. The annotations from NCBI must be hand edited after running the script. In order for the 
script to run you need a .fnn and .fasta file of your sample, as well you much save a 
copy of the fast file as temp.fasta to run the no hits against NCBI or else you will
 run an error. The Usearch_Features will output a .tab file with the same name as seq
 which will contain the annotation information.

3. extract_features_tab(seq): this function will create a GenBank annotation based on
 the .tab file from Usearch_Features. First the .tab file must be saved as a .csv
 after any manual curation then can be run to create the GenBank file.

4. Create_Usearch_Database(): this function creates a Usearch database with the only
 input being the file name you want to be made into the database.

5. composite_fasta(folder_directory,outfile): this function will take the output from
 the backbone gene database and place all of the sequences in one fasta file.

6. Translate_Fasta('ResistanceGenes.fasta','ResistanceGenesAA.fasta’): this tool will
 translate a NA fast file into AA.