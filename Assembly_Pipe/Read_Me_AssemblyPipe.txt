AssemblyPipe.py

Used to clean and assemble Illumina MiSeq reads into larger contigs

Dependencies:
	Bowtie2 http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
	Spades http://cab.spbu.ru/software/spades/
	BioPython http://biopython.org/wiki/Documentation

Functions:

Bowtie2_Build_reference(reference,outname)
	Used to create the database for Bowtie2 to align too.
	reference = name of Fasta file to create database from
	outname = the desired name of the database

Bowtie2_align(reference,path,cleanreadfname,outreadfname)
	Removes all Illumina reads that match to a desired database, used to remove genomic reads.
	reference = name of the reference database for the alignment
	path = path where the folder with reads can be found
	cleanreadfname = name of the folder with the cleaned reads within the specified path
	outreadfname = name of the new directory that will contain the reads minus those that aligned to the database
	
AssembleEach_Plasmid_Spades(path, cleanreadfname ,outreadfname)
	Performes a de-novo assembly of Ilumina MiSeq reads
	path = path were the directory containing all the reads is found
	cleanreadfname = directory that contains all reads you want to be assembled
	outreadfname = name of new directory that will contain all the assembled reads