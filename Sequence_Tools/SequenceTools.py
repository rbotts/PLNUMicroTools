from Bio.Seq import Seq, MutableSeq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature,FeatureLocation
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO, Entrez, Alphabet
import os
import csv
import sys
import numpy as np
from Bio.Blast.Applications import NcbiblastxCommandline

#from user_plasmids import *
#from expanding import *
#from copy import deepcopy

###########################################


# routines for sequence annotation
def composite_fasta(folder_dir,new_fasta_name):
	outfile = open(new_fasta_name, "a")
	
	for filename in os.listdir(folder_dir):
		filedir = folder_dir+'/'+filename
		
		## open each fasta file in the directory
		record = SeqIO.parse(filedir, "fasta")
		
		for feature in record:
			## rename each feature based on the filename
			feature.id = (filename[:-3])
			feature.name = (filename[:-3])
			feature.description = (filename[:-3])
			SeqIO.write(feature,outfile,"fasta")
			#print feature
					
	outfile.close()

def Usearch_Features(seqname,database):
	import datetime
	# open contents of seq file, add features from file of orfs
	# compares with local database
	usearch_out=seqname+'_usearch.tab'
	usearch_out_no_hits=seqname+'_no_hits.tab'
	outFileName=seqname+'.tab'
	os.system('./usearch -ublast '+seqname+'.fnn -db '+database+' -evalue 0.04 -blast6out '+usearch_out+' -notmatched temp.fnn -strand plus -top_hit_only -id .70')
	#os.system('./usearch -ublast '+seqname+'.fnn -db '+database+' -evalue 0.04 -blast6out '+usearch_out+' -strand plus -top_hit_only -output_no_hits -id .70')
	
	# Reads in the results and makes them into a list of list object 
	results_table = list(csv.reader(open(usearch_out, 'rb'), delimiter='\t'))
	
	# Create File with readable format for the annotations
	outfile = open(outFileName, "a") # use a to avoid overwriting old data
	outfile.write(seqname+'\t'+datetime.datetime.now().strftime("%m/%d/%y")+'\n')
	outfile.write('ORF\t'+'start\t'+'stop\t'+'direction\t'+'size (bp)\t'+'size (aa)\t'+'Gene\t'+\
		'definition\t'+'GB \t'+'Type \t'+'query cov\t'+'query id\t'+'comment\n')
	
	#make data into previous format	
	for row in results_table:
		aa_num= str(int(row[3])/3)
		start= str((int(row[0].split('|')[4])+int(row[6]))-1)
		stop= str((int(row[0].split('|')[4])+int(row[7]))-1)
		#print aa_num
		outfile.write(row[0].split('|')[0]+'\t'+start+'\t'+stop+'\t'+row[0].split('|')[3]+'1\t'+row[3]+'\t'+aa_num+'\t'+row[1]+ \
		'\t\t\tCDS\t\t'+row[2]+'\t')
		outfile.write('\n')
	
	# searches all ORF's that did not have a hit and runs blast on them
	# identical to Multi_Brooke, except uses features on a seqRecord
	outcsv=True #True for outputting to csv file, false produces tab delimited output
	blast_type="blastx" # set type of blast: n, p or x, translate seq record for blastp
	record=add_features_fn('temp')
	#record = SeqIO.parse('temp.fasta', "fasta")
	geneList=[]

	for feature in record.features:
		#call blast with nr database
		# print feature.id
		result_handle = 3(blast_type, "nr", feature.extract(record).seq)
		giList = []
	
		#open file for writing and then save the blast results
		#xml format is chosen because parser is more stable
		save_file = open ("temp_blast.xml", "w")
		save_file.write(result_handle.read())
		save_file.close()
		result_handle.close()
		
		# open contents of Blast output and write into nice tabular format
		# works the same as Brook_Blast, but output format is different
		# outfield = open(outFileName, "w") #open for appending so don't write over previous results

		#bring file contents back into handle for reading
		result_handle = open("temp_blast.xml")
		blast_record = NCBIXML.read(result_handle)
		E_VALUE_THRESH = 0.04
		found_best = False
		count = 0

		for alignment in blast_record.alignments:
			for hsp in alignment.hsps:
				if hsp.expect < E_VALUE_THRESH:
					#Extract accession number and description from title
					title_info = alignment.title.split('|')
					accession=title_info[3]
					descript=title_info[4]
					gi_num = title_info[1]

					#find the best match that is not hypothetical
					#save all others with same max_iden and query cover as that best
					if descript.find('unknown')==-1 and descript.find('hypothetical')==-1 and not (gi_num in giList):
						giList.append(gi_num)
						#calc max iden and query cover
						cover = (len(hsp.query)*1.0/len(feature))
						max_iden = (hsp.identities * 1.0)/len(hsp.query)

						#save the information about the best match
						if not found_best:
							best_cover=cover
							best_iden=max_iden
							found_best=True
							print ("query cover: %s	max_iden: %s" % (cover, max_iden))

						#collect all other alignments that have same "best" scores
						if cover==best_cover and max_iden==best_iden:
							#calc subject start and end based on strand directions
							if hsp.sbjct_start < hsp.sbjct_end:
								sbjct_end = hsp.sbjct_end;
								sbjct_start = hsp.sbjct_start
							else:
								sbjct_end = hsp.sbjct_start
								sbjct_start = hsp.sbjct_end
							geneList.append((gi_num,sbjct_start,sbjct_end, accession, descript, cover, max_iden))
							count=count+1
							outfile.write("%s\t %i\t %i\t %i\t %i\t %i\t \t %s\t %s\t \t %f\t %f\t \n" \
							% (feature.id, int(feature.location.start), int(feature.location.end),\
							feature.strand, len(feature)+1, (len(feature)+1)/3, descript, accession,\
							cover, max_iden))
							print("%s \t %s \t %s \t %s \t %s" % (gi_num,accession,descript,\
								str(sbjct_start),str(sbjct_end)))
		# if no suitable alignments are  found insert a blank line in the annotation  table
		if giList==[]:
			outfile.write("\n")
		else:
			if outcsv:
				outfile.write("%s\t %i\t %i\t %i\t %i\t %i\t \t \t \t \t \t \t \n" \
							 % (feature.id, int(feature.location.start), int(feature.location.end),\
							 feature.strand, len(feature)+1, (len(feature)+1)/3))
			else:
				outfile.write("%s\t %i\t %i\t %i\t %i\t %i\t \t \t \t \t \t \t \n" \
						% (feature.id, int(feature.location.start), int(feature.location.end),\
						feature.strand, len(feature)+1, (len(feature)+1)/3))
			
			outfile.write("\n") #print a space between orfs
	print("----------------------- ORF BLAST COMPLETE ------- %s records found----------" % (str(count)))
	outfile.write("\n\n")

	
	outfile.close()

def Create_Usearch_Database(databasefile):
	# create usearch database to be used in usearch features
	os.system('./usearch -makeudb_ublast '+databasefile+' -output '+databasefile[:-6]+'.udb')
	
def Translate_Fasta(filename,outfile):
	# Translate a fasta file of nucleic acids to amino acids
	aa_file = open(outfile, "a")
	na_file = SeqIO.parse(filename, "fasta")
	for feature in na_file:
		feature.seq = feature.seq.translate()
		SeqIO.write(feature,aa_file,"fasta")
	aa_file.close()

def Blast_Features(seqname):
	# open contents of seq file, add features from file of orfs uses nucleotide
	# blast each file and return a single output table
	# identical to Multi_Brooke, except uses features on a seqRecord
	# loop through ORFs and print results in outFile
	import datetime
	outcsv=True #True for outputting to csv file, false produces tab delimited output
	blast_type="blastx" # set type of blast: n, p or x, translate seq record for blastp
	record=add_features_fn(seqname)
	outFileName=seqname+'.tab'
	geneList=[]
	outfile = open(outFileName, "a") # use a to avoid overwriting old data
	outfile.write(seqname+'\t'+datetime.datetime.now().strftime("%m/%d/%y")+'\n')
	outfile.write('ORF\t'+'start\t'+'stop\t'+'direction\t'+'size (bp)\t'+'size (aa)\t'+'Gene\t'+\
		'definition\t'+'GB \t'+'Type \t'+'query cov\t'+'query id\t'+'comment\n')

	for feature in record.features:
		#call blast with nr database
		# print feature.id
		result_handle = NCBIWWW.qblast(blast_type, "nr", feature.extract(record).seq)
		giList = []
	
		#open file for writing and then save the blast results
		#xml format is chosen because parser is more stable
		save_file = open ("temp_blast.xml", "w")
		save_file.write(result_handle.read())
		save_file.close()
		result_handle.close()
		
		# open contents of Blast output and write into nice tabular format
		# works the same as Brook_Blast, but output format is different
		# outfield = open(outFileName, "w") #open for appending so don't write over previous results

		#bring file contents back into handle for reading
		result_handle = open("temp_blast.xml")
		blast_record = NCBIXML.read(result_handle)
		E_VALUE_THRESH = 0.04
		found_best = False
		count = 0

		for alignment in blast_record.alignments:
			for hsp in alignment.hsps:
				if hsp.expect < E_VALUE_THRESH:
					#Extract accession number and description from title
					title_info = alignment.title.split('|')
					accession=title_info[3]
					descript=title_info[4]
					gi_num = title_info[1]

					#find the best match that is not hypothetical
					#save all others with same max_iden and query cover as that best
					if descript.find('unknown')==-1 and descript.find('hypothetical')==-1 and not (gi_num in giList):
						giList.append(gi_num)
						#calc max iden and query cover
						cover = (len(hsp.query)*1.0/len(feature))
						max_iden = (hsp.identities * 1.0)/len(hsp.query)

						#save the information about the best match
						if not found_best:
							best_cover=cover
							best_iden=max_iden
							found_best=True
							print ("query cover: %s	max_iden: %s" % (cover, max_iden))

						#collect all other alignments that have same "best" scores
						if cover==best_cover and max_iden==best_iden:
							#calc subject start and end based on strand directions
							if hsp.sbjct_start < hsp.sbjct_end:
								sbjct_end = hsp.sbjct_end;
								sbjct_start = hsp.sbjct_start
							else:
								sbjct_end = hsp.sbjct_start
								sbjct_start = hsp.sbjct_end
							geneList.append((gi_num,sbjct_start,sbjct_end, accession, descript, cover, max_iden))
							count=count+1
							outfile.write("%s\t %i\t %i\t %i\t %i\t %i\t \t %s\t %s\t \t %f\t %f\t \n" \
							% (feature.id, int(feature.location.start), int(feature.location.end),\
							feature.strand, len(feature)+1, (len(feature)+1)/3, descript, accession,\
							cover, max_iden))
							print("%s \t %s \t %s \t %s \t %s" % (gi_num,accession,descript,\
								str(sbjct_start),str(sbjct_end)))
		# if no suitable alignments are  found insert a blank line in the annotation  table
		if giList==[]:
			outfile.write("\n")
		else:
			if outcsv:
				outfile.write("%s\t %i\t %i\t %i\t %i\t %i\t \t \t \t \t \t \t \n" \
							 % (feature.id, int(feature.location.start), int(feature.location.end),\
							 feature.strand, len(feature)+1, (len(feature)+1)/3))
			else:
				outfile.write("%s\t %i\t %i\t %i\t %i\t %i\t \t \t \t \t \t \t \n" \
						% (feature.id, int(feature.location.start), int(feature.location.end),\
						feature.strand, len(feature)+1, (len(feature)+1)/3))
			
			outfile.write("\n") #print a space between orfs
	print("----------------------- ORF BLAST COMPLETE ------- %s records found----------" % (str(count)))
	outfile.write("\n\n")
	outfile.close()

def elim_Genomic_Contigs_BLAST(seqfname,minlength=25000,maxlength=200000, minlenalign= 6000, minperalign=40):
	# routine takes the result of blasting genomic sequence against contig file finds contigs between
	# the min and max length that do not align well with genomic DNA
	#
	# Removes contigs that have large regions that align with genomic contigs.
	# Should add criteria that then looks at ORF's on likely contigs and blasts each one,
	# keeping contigs that have orfs that have plasmid genes.
	#
	# inputs:
	#	falign - alignment text file, output from Blast
	#	fseqs - file of seqs
	#	fpathname -  optional pathname used for both files, adds to defaults
	#   minperaling - minimum percent of alignment to call genomic content
	#	minlength - minimum contig length worth looking at
	#	maxlength - longest reasonable contig length that could be plasmid
	#	minlenalign - minimum length alignment with genome that seems significant.
	# 	minperalign - minimum  percentage align to count
	#	minnumalign - minimum number of quality alignments to exclude as genomic DNA
	import csv
	path="//Users/lucasustick/SR/2015Contigs/"+seqfname+"/"
	fseqs=seqfname+"LargeContigs.fna"
	outname = seqfname+"-LikelyContigs.fasta"
	outhandle=open(outname,'w')
	outcontigs=[]
	data=[]
	
	E_VALUE_THRESH = 0.04
	
	# open multi fasta file
	inseqs = SeqIO.parse(open(path+fseqs,"rU"),"fasta")
	for record in inseqs:
		if len(record.seq)>minlength and len(record.seq) < maxlength:
			count=0
			print 'BLAST-ing record '+record.id
			# BLAST record
			result_handle = NCBIWWW.qblast('blastn', 'nr', record.seq)
			
			#open file for writing and then save the blast results
			#xml format is chosen because parser is more stable
			save_file = open("temp_blast.xml", "w")
			save_file.write(result_handle.read())
			save_file.close()
			result_handle.close()
			
			
			# read BLAST output
			result_handle = open("temp_blast.xml", "r")
			blast_record = NCBIXML.read(result_handle)
			# run through each hit and find hits that show high alignment to regions that are not plasmid
			# also check for multiple unique regions that are genomic
			notgenome = True
			for alignment in blast_record.alignments:
				# check if the alignment is genomic and if it is greater than minlealign
				if alignment.hit_def.find('complete genome')>-1:
					alignlen=0
					for hsp in alignment.hsps:
						if hsp.expect < E_VALUE_THRESH:
							alignlen+=hsp.align_length
					if alignlen>minlenalign or alignlen/len(record.seq)>minperalign:
						notgenome = False
						print 'Removing genomic contig, aligned with ' + alignment.hit_def
						print 'Length of alignment= %i' % alignlen
						break # break for loop because we can throw the seq out
			if notgenome:
				print 'Writing contig.  Length %i' % len(record.seq)
				#os.rename('temp_blast.xml',seq.id+'xml')
				SeqIO.write(record,outhandle,'fasta')
				outcontigs.append(record)

	print "You have found "+str(len(outcontigs))+" matching your criteria"
	outhandle.close()

def elim_Genomic_Contigs_BLAST2Seqs(seqfname,minlength=5000,maxlength=200000, minnumalign=5,minlenalign=300, minperalign=50):
	# routine takes the result of blasting genomic sequence against contig file finds contigs between
	# the min and max length that do not align well with genomic DNA
	# inputs:
	#	falign - alignment text file, output from Blast
	#	fseqs - file of seqs
	#	fpathname -  optional pathname used for both files, adds to defaults
	#   minperaling - minimum percent of alignment to call genomic content
	#	minlength - minimum contig length worth looking at
	#	maxlength - longest reasonable contig length that could be plasmid
	#	minlenalign - minimum length alignment with genome that seems significant.
	# 	minperalign - minimum  percentage align to count
	#	minnumalign - minimum number of quality alignments to exclude as genomic DNA

	path="//Users/rbotts/Documents/DATA/Projects/BioPythonTools/SequenceData/SU2013/"+seqfname+"/"
	fseqs=seqfname+"LargeContigs.fna"
	falign=seqfname+"-Alignment.txt"
	outcontigs=[]
	data=[]
	
	# read in data from blast alignment
	datareader =csv.reader(open(path+falign,"r"),delimiter="\t")
	for row in datareader:
		# skip blank rows
		if any(row):
			#skip header rows
			if row[0] != "#":
				#print row
				data.append(row)
	# open seq file with contigs
	inseqs = SeqIO.parse(open(path+fseqs,"rU"),"fasta")
	for record in inseqs:
		if len(record.seq)>minlength and len(record.seq) < maxlength:
			#print record.id, len(record.seq)
			count=0
			for align in data:
				if len(align)>1:
					if align[1] == record.id and int(align[3]) > minlenalign and float(align[2]) > minperalign:
						count+=1
			if count<= minnumalign:
				outcontigs.append(record)
				print count, record.id, len(record.seq)
			#print record.id, len(record.seq)
	print "You have found "+str(len(outcontigs))+" matching your criteria"
	# write a single output file
	output_handle = open(path+seqfname+"LikelyContigs.fasta", "w")
	SeqIO.write(outcontigs, output_handle, "fasta")
	output_handle.close()

def add_features_fn(seqname):
	# open contents of seq file (.fa) in fasta format and associated fnn file with
    # multiple features as separate fasta sequences (output from GenemarkS
	#
	# inserts contents of fnn file as seq features on the seq object
	# generates an output file of the same name with features to be used later
	
    
    sqname = seqname
	
    # read in seq to add features to
    handle = open(sqname+".fasta","rU")
    record = SeqIO.read(handle,"fasta")
    handle.close()
    
    record.id = sqname
    
    #read contents of gmm file.  gmm output is not standard tab or comma delimited file
    gmm_handle = open(sqname+".fnn", "rU")
    
    features_rec = list(SeqIO.parse(gmm_handle, "fasta"))
    for fture in features_rec:
        feat_info = fture.id.split('|')
        if feat_info[-3] == '+':
            feat_dir = 1
        else:
            feat_dir = -1
        record.features.append(SeqFeature(id=feat_info[0],location=FeatureLocation(int(feat_info[-2])-1,\
                                                                                   int(feat_info[-1]),strand=feat_dir)))
	
    print "Extracted %i features" % len(record.features)
    return record

def extract_features_tab(seqname,writegb=True):
	# Extract the annotation features from a csv delimited file of the form annotation file
	#
	# writegb switches writing output on
	# returns a record with the annotations stored in seqname.tab
	


	# open seq record file
	handle = open(seqname+".fasta","rU") 
	record = SeqIO.read(handle,"fasta") 
	handle.close()
	
	record.seq.alphabet=Alphabet.generic_dna #set alphabet for use later
	
	# open annotation file, first two lines are headers 
	ann_handle = open(seqname+".csv","rU") 

	# read in data
	csv_file = csv.reader(ann_handle,delimiter=',',dialect=csv.excel)
	
	# extract info in each row and add create temp seq objects to write to image
	record.name=seqname[0:5] # seq record prints the name as locus_id, and there is a maximum length
	record.id = csv_file.next()[0] # skip header row
	print record.id
	csv_file.next() # skip date row
	orfcount=1  # count unnamed orfs for naming
	i = 0
	for row in csv_file:
		i = i+1
		# catch cases where row contains invalid info 
		print row
		try:
			if row[6]=='':
				feat_name='orf'+str(orfcount)
				orfcount=orfcount+1
			else:
				feat_name = row[6]
			print feat_name
			feat_start = int(row[1])
			feat_end = int(row[2])
			feat_dir = int(row[3])
			#feat_acc = row[8] Not workingright
			feat_note=row[7]
			if row[9]== '':
				feat_type='CDS'
			else:
				feat_type = row[9]
			record.features.append(SeqFeature(id=feat_name,type=feat_type,location=FeatureLocation(feat_start-1,feat_end,strand = feat_dir),\
						qualifiers={"product":feat_name,"protein_id":str(i)}))

		except ValueError:
			print "A row contains invalid entries, row skipped"

	ann_handle.close()
	for f in record.features:
		print f
	if writegb:
		outhandle=open(seqname+".gb","w")
		SeqIO.write(record,outhandle,"genbank")
		outhandle.close()
		print "Successfully wrote " + seqname + ".gb"
	return record


def Brooke_blast(seqStr,anOutFile, geneList,blasttype="blastn"):
	#change string into a Seq object if necessary
	if type(seqStr) == str:
		my_seq = Seq(seqStr)
	else:
		my_seq=seqStr
	#call blast using blastn algorithm with nr database
	result_handle = NCBIWWW.qblast(blasttype, "nr", my_seq)
	giList = []
	
	#open file for writing and then save the blast results
	#xml format is chosen because parser is more stable
	save_file = open ("my_blast.xml", "w")
	save_file.write(result_handle.read())
	save_file.close()
	result_handle.close()
	outfile = open(anOutFile, "a") #open for appending so don't write over previous results

	#bring file contents back into handle for reading
	safe_file = open("my_blast.xml")
	blast_record = NCBIXML.read(result_handle)
	E_VALUE_THRESH = 0.04
	found_best = False
	count =0
	#this deals with a single query
  
	for alignment in blast_record.alignments:
		for hsp in alignment.hsps:
			if hsp.expect < E_VALUE_THRESH:
  
				#Extract accession number and description from title
				title_info = alignment.title.split('|')
				accession=title_info[3]
				descript=title_info[4]
				gi_num = title_info[1]

				#find the best match that is not hypothetical
				#save all others with same max_iden and query cover as that best
				if descript.find('unknown')==-1 and descript.find('hypothetical')==-1 and not (gi_num in giList):	
					giList.append(gi_num)
					
					#calc max iden and query cover
					cover = (len(hsp.query)*1.0/len(seqStr))
					max_iden = (hsp.identities * 1.0)/len(hsp.query)

					#save the information about the best match
					if not found_best:
						best_cover=cover
						best_iden=max_iden
						found_best=True
						print ("query cover: %s	max_iden: %s" % (cover, max_iden))
						outfile.write(str(cover)+"\t"+str(max_iden)+"\n")

					#collect all other alignments that have same "best" scores
					if cover==best_cover and max_iden==best_iden:
						
						#calc subject start and end based on strand directions 
						if hsp.sbjct_start < hsp.sbjct_end:
							sbjct_end = hsp.sbjct_end;
							sbjct_start = hsp.sbjct_start
						else:
							sbjct_end = hsp.sbjct_start
							sbjct_start = hsp.sbjct_end
						geneList.append((gi_num,sbjct_start,sbjct_end, accession, descript, cover, max_iden))
						count=count+1
						outfile.write(gi_num+"\t"+accession+"\t"+descript+"\t"+ str(sbjct_start)+"\t"+ str(sbjct_end)+"\n")
						print("%s \t %s \t %s \t %s \t %s" % (gi_num,accession,descript,str(sbjct_start),str(sbjct_end)))
	print("----------------------- ORF BLAST COMPLETE ------- %s records found----------" % (str(count)))
	outfile.write("\n\n")
	outfile.close()


def multi_Brooke(ORFfile,outFile):
	ORFsFile = open(ORFfile, "r")
	#loop through ORFs and print results in outFile
	geneList=[]
	for anORF in ORFsFile:
		Brooke_blast(anORF, outFile, geneList)

	#giList contains all gi_nums to check in nucleotide database
	Entrez.email = "rbotts@pointloma.edu" #Let NCBI know who I am

	#for each entry in geneList, get the tuple components
	for gi_num, sbj_start, sbj_end, accession, descript, cover, max_iden in geneList:
	   # get_genbank_info(gi_num,sbj_start,sbj_end, accession, descript,cover,max_iden, outFile)
		access_Genbank(gi_num,sbj_start,sbj_end, accession, descript,cover,max_iden, outFile)

def access_Genbank(gi_num,beginBLseq,endBLseq, accession, descript, cover, max_iden, anOutFile):
	#should be called by multi_Brooke after ORF is sent to BLAST
	#accesses the genbank record for the given gi number and finds annotation related to the sequence slice found by BLAST

	#for when ready to write output to file
   # outfile = open(anOutFile, "a") #open for appending so don't write over previous results

   #if call standalone, must give Entrez your email address so they can alert you of too  much usage
   #Entrez.email = "lcarter@pointloma.edu" #Let NCBI know who I am
   
	#access NCBI Genbank and turn result into SeqRecord for parsing
	handle=Entrez.efetch(db="nucleotide", id = str(gi_num), rettype = "gb", retmode="text")
	record = SeqIO.read(handle, "gb")
	handle.close()

	#Get basic feature information
	print "%s with %i features \n" %(record.id, len(record.features))

	#look through all of the features, searching for the one that matches the slice of sequence located by BLAST
	for feature in record.features:
		if feature.type=='CDS':
			# find the reference to the gene, and remove extraneous characters for printing
			gene = str(feature.qualifiers.get('gene'))
			if gene != 'None':
				gene = gene[2:(len(gene)-2)]

			# find the reference to the product, and remove extraneous characters for printing
			product = str(feature.qualifiers.get('product'))
			if product != 'None':
				product =product[2:(len(product)-2)]

			# find the location information, and extract start and end ints from the range given (appears as [start:end])
			location = str(feature.location)
			colon = location.find(':')
			location1= int(location[1:colon])
			location2=int(location[colon+1:location.find(']')])

			#we've found the right one, so print it
			if beginBLseq>=location1 and endBLseq<=location2:
				print("%s  %s   cover: %s   max_iden: %s" % (str(gi_num), str(accession), str(cover), str(max_iden)))
				print ("%s   gene: %s   product: %s   %s ... %s" % (descript, gene, product, str(location1), str(location2)))
				print("--------------------------------------------------------------")
				


def seq_map(record, DiagramLabel,outFile):
	# take all of the records in seqfile and include in a circular map
	# wrapper for biopython tools, from biopython cookbook
	# have not tested
	from reportlab.lib import colors
	from reportlab.lib.units import cm
	from Bio.Graphics import GenomeDiagram
	from Bio import SeqIO
	#record = SeqIO.read(seqfile, "fasta") #use if you are passing in a seqfile
	#Create the feature set and its feature objects,

	gd_feature_set = GenomeDiagram.FeatureSet()
	for feature in record.features:
		#if feature.type != "gene":
		#Exclude this feature
		#	continue
		if len(gd_feature_set) % 2 == 0:
			color = colors.blue
		else:
			color = colors.lightblue
		gd_feature_set.add_feature(feature, feature_id=feature.id,sigil="ARROW",color=color, label=True,label_size = 14, label_angle=0)
	#(this for loop is the same as in the previous example)
	#Create a track, and a diagram
	gd_track_for_features = GenomeDiagram.Track(name="Annotated Features",height=1)
	gd_diagram = GenomeDiagram.Diagram(DiagramLabel)
	#Now have to glue the bits together...
	gd_track_for_features.add_set(gd_feature_set)
	gd_diagram.add_track(gd_track_for_features, 1)

	gd_diagram.draw(format="linear", orientation="landscape", pagesize='A4',
			  fragments=4, start=0, end=len(record)) # format may be linear or circular
	gd_diagram.write(outFile+'.eps', "EPS")


###########################################
# routines for clone curation and cleanup, from CTX-M project
prim1=Seq('ATGGTTAAA')
prim2=Seq('CGGTTTGTAA').reverse_complement()
#prim1=Seq('TTCGTCTCTTCCAGAATAAGG') #hard code CTX-M group 1 primers
#prim2=Seq('CAGCACTTTTGCCGTCTAAG')

def retRecords(foldername):
	# retrieves all fasta sequences in foldername
	# handles I/O
	# inputs: foldername
	# !!!!! only reads seq files in fasta format, no errors if in the incorrect format
	# outputs: seqarray (array of sequences from files)
	
	seqarray=[] 	

	# change dir to foldername
	os.chdir(foldername)
	for f in os.listdir("."):
		# find all seq files
		if f.endswith(".seq"):
			# read in seq to add features to
			handle = open(f,"rU") 
			record = SeqIO.read(handle,"fasta") 
			handle.close()
			record.name = f[f.find('X'):f.find('.seq')-5] # beware not a hard copy
			seqarray.append(record) # add the record to the list of outputs
			# print f+" successfully read, sequence length= %i" % len(record.seq)
	os.chdir('..')
	print 'successfully read %i sequences' % len(seqarray)
	print seqarray[-1]
	return seqarray

			
def trimSeqs(forward='pfasta-Nov26-11-35-54',rev='pfasta-Nov26-11-39-0'):
	# trims all sequences in foldername
	# writes two files, one with seqs successfully trimmed and one that is not
	# inputs:  
	#	forward-- file with forward reads
	# 	rev-- file with backward reads
	# outputs:
	#	outseqs---array of all successfully trimmed sequences
	#	failedtofind---array of sequences that were not successfully trimmed

	# output file names for fasta format
	outname="TrimSeqs.fasta"	
	untrname="UntrimmedSeqs.fasta"

	outseqs =[]
	ftfind=[]
	
	# read in sequences from files, no trimming
	seqs=retRecords(forward)
	revseqs=retRecords(rev)
	

	# select a maximum length to avoid having junk reads
	maxlength=800

	#trim sequences
	for i in range(len(seqs)):
		foundinboth=0 # keep track of successful trimming

		start=seqs[i].seq.find(prim1)
		if start>0:
			start=start+len(prim1)
			seqs[i].seq=seqs[i].seq[start:start+maxlength]  
			# seqs[i].seq=seqs[i].seq[:seqs[i].seq.find('NNNNN')]  
			print 'trimmed with primer 1, start=%i, length=%i' % (start, len(seqs[i]))
			
			#trim reverse read with primer 2
			start2=revseqs[i].seq.find(prim2)
			if start2>0:
				start2=start2+len(prim2)
				print 'start site for reverse read: %i' % start2
				revseqs[i].seq=revseqs[i].seq[start2:start2+maxlength].reverse_complement() 
				# revseqs[i].seq=revseqs[i].seq[:revseqs[i].seq.find('NNNNN')].reverse_complement() #
				
				outseqs+=[seqs[i], revseqs[i]] # save the forward and reverse seqs for later

				foundinboth=True
		else:
			start=seqs[i].seq.find(prim2)
			if start>0:
				start=start+len(prim2)
				seqs[i].seq=seqs[i].seq[start:start+maxlength].reverse_complement()
				# seqs[i].seq=seqs[i].seq[:seqs[i].seq.find('NNNNN')] 
				print 'trimmed with primer 2, start=%i, length=%i' % (start,  len(seqs[i]))

				#trim reverse read with primer 2
				start2=revseqs[i].seq.find(prim1)
				if start2>0:
					start2=start2+len(prim1)
					print 'start site for reverse read: %i' % start2
					revseqs[i].seq=revseqs[i].seq[start2+len(prim1):start2+len(prim1)+maxlength]
					# revseqs[i].seq=revseqs[i].seq[:revseqs[i].seq.find('NNNNN')] 
					outseqs+=[revseqs[i],seqs[i]] 
					# all are oriented so that the first is trimmed by primer1 and the second sequence 
					# is the end trimmed by primer2
					foundinboth=True
		if foundinboth:
			print  'found both primers in %s' % seqs[i].name
		
		else:
			ftfind+=[seqs[i],revseqs[i]]
			print 'did not find both primers'
	
	print len(outseqs)
	print len(ftfind)
	outhandle=open(outname, "w")
	SeqIO.write(outseqs, outhandle, "fasta")
	outhandle.close()

	out_untrimmed=open(untrname, "w")
	SeqIO.write(ftfind, out_untrimmed, "fasta")
	out_untrimmed.close()
	
	return outseqs,ftfind

def FindIdenticalSeqs(seqfile,outfasta,outcsv,fixedlength=876):
	# find identical seqs in a fasta file and output a new fasta file with all duplicates removed along with a csv file
	# listing which sequences are identical
	# seqfile- name of the fasta file containing the seqs
	# outfasta - the  output fasta file
	# outcsv - name of the output csv file
	# fixed length specifies the length of all seqs in the file,
	# this is optional and not really necessary, but a good test that all sequences have the same length
	path = "//Users/rbotts/Documents/DATA/Projects/BioPythonTools/CTX-MStudy/CTX-MFinal/"
	isnucleotide = True

	templibrary={} 	# array of all seqs read in
	maxdist=0	# maximum distance found in the list

	# change dir to filename
	os.chdir(path)
	# open input and output files
	inhandle = open(seqfile,"rU")
	outseqs = open(outfasta,"w")
	
	firsttime = True

	for record in SeqIO.parse(inhandle,"fasta"):
		# check that the record has the write length and no odd characters
		if len(record) != fixedlength:
			print record.id + ' is not the correct length, its length is ',len(record)
		# check sequence for errors if it is nucleotide, and notify the user
		if isnucleotide:
			for i in range(len(record)):
				if not record[i] in ['A','T','C','G']:
					print record.id+' has a '+ record[i] + ' at location '+ str(i)

		# keep track of whether we have found a duplicate
		isfound = False
		for dups in templibrary:
			# save the reference to the known seq for ease of use
			knownseq = templibrary[dups][0]
			
			if NumDiffs(record,knownseq)==0:
				# find duplicate records and save the name for the csv
				print 'Duplicate sequence found: '+record.id+' is the same as '+ knownseq.id
				# append the duplicate to the end of the list to the one it is a duplicate of.
				templibrary[dups].append(record)
				isfound=True
		if not isfound:
			# add the new item to templibrary and create a key to it
			templibrary[record.id]=[record]
			SeqIO.write(record,outseqs,"fasta")

	print 'The number of distinct sequences is ',len(templibrary)
	inhandle.close()
	outseqs.close()
	outcsvhandle = open(outcsv,"w")
	# write each list of duplicate sequences to an output csv file.
	for dups in templibrary:
		if len(templibrary[dups])>0:
			for seq in templibrary[dups]:
				outcsvhandle.write(seq.id+',')
			outcsvhandle.write('\n')
	outcsvhandle.close()





def GroupSeqs(max_allow_dist=0,filename="//Users/rbotts/Documents/DATA/Projects/BioPythonTools/CTX-MClones/"):
	# open all seqs in filename
	# max_allow_dist maximum accepted distance, inclusive
	# change name to standard naming convention
	# identify duplicates at some threshold
	# output single file in fasta format
	# output table of sequence differences	
	# output file names for fasta format

	outname="//Users/rbotts/Documents/DATA/Projects/BioPythonTools/CTX-MLibTest"+str(max_allow_dist)+".fasta"
	   #"//Volumes/Users/rbotts/DATA/Projects/BioPythonTools/CTXM-Clones/CTX-MLib"+str(max_allow_dist)+".fasta"

	outdataname="//Users/rbotts/Documents/DATA/Projects/BioPythonTools/CTX-MLibTest"+str(max_allow_dist)+"Counts.txt"
	fixedlength=876 # check to make sure the sequences are this length before adding to the list

	templibrary=[] 	# array of all seqs read in
	maxdist=0	# maximum distance found in the list

	# change dir to filename
	os.chdir(filename)
	# print os.listdir(".")
	for f in os.listdir("."):
		# find all seq files
		if f.endswith(".gb"):
			# read in seq to add features to
			handle = open(f,"rU") 
			record = SeqIO.read(handle,"gb") 
			handle.close()

			record.id = 'X-1_'+f[f.find('X')+1:record.id.find('.')-2] # beware not a hard copy
			record.name = ''#record.id
			record.description=''
			record.seq=record.seq.ungap('-') # trims out any gaps
			# identify distinct sequences any of the inappropriate bases

			for i in range(len(record)):
				if not record[i] in ['A','T','C','G']:
					print record.id+' has a '+ record[i] + ' at location '+ str(i)

			templibrary.append(record) # add the record to the list of outputs
			# print f+" successfully read, sequence length= %i" % len(record.seq)
	os.chdir('..')
	print 'successfully read %i sequences' % len(templibrary)

	
	seqlibrary=[templibrary[0]] # stores distinct seqs
	wronglength=[] # store sequences of the wrong length separately


	# identify distinct sequences and put them in a new list, remove redundant sequences according to cutoff
	for f in templibrary:
		if len(f) != fixedlength:
			print f.id +' is not 876 bp, actual length is '+ str(len(f))
			wronglength.append(f)
		else:
			mindist=fixedlength #used to keep track of the smallest observed distance
			for g in seqlibrary:
				diff = NumDiffs(f,g)
				#print f.id + ' is '+str(diff)+' bp different from '+g.id
				if diff < mindist:	
				 	mindist=diff
				# keep track of the maximum observed diet
				if diff > maxdist:
					maxdist=diff # update maximum observed distance
				#print f.id + ' is ' +str(mindset)

			# if f is not a duplicate add it to the list of distinct sequences
			if mindist  > max_allow_dist:
				# should check here if the seqs have any uncalled bases
				seqlibrary.append(f) # add distinct sequence to list
				
			
	#output distinct sequences and sequences of incorrect length that cannot be easily compared		
	print 'Number of distinct sequences '+str(len(seqlibrary))
	print 'Number of incorrect length seqs '+str(len(wronglength))
	

	print 'Maximum observed distance ' + str(maxdist)

	if True:
		outhandle=open(outname, "w")
		SeqIO.write(seqlibrary, outhandle, "fasta")
		outhandle.close()
	
	# Count the number of replicates of each
	if True:
		count=[0 for f in seqlibrary]
		for i in range(len(count)):
			for g in templibrary:
				if NumDiffs(seqlibrary[i],g) <= max_allow_dist:
					count[i]+=1
					#print seqlibrary[i].id + ' is identical to ' + g.id
					templibrary.remove(g) #remove already counted sequence from the list
					#seems to be a problem when the max_dist is > 0, some sequences aren't counted
			print str(count[i]) + ' replicates of ' + seqlibrary[i].id
		
		print sum(count)

		#print count
		outhandle=open(outdataname, "w")
		for i in range(len(count)):
			outhandle.write(seqlibrary[i].id + "\t" + str(count[i]) + "\n")
		outhandle.close()
	return seqlibrary

def NumDiffs(seq1,seq2):
	count = 0

	# Ignore differences that are ambiguous bases ['N','R','Y','K','M','S','W','B','D','H','V','X']
	# could improve this logic
	for i in range(len(seq1)):  
		if (seq1[i] != seq2[i]): #and (seq1[i] in ['A','T','C','G']) and (seq2[i] in ['A','T','C','G']) :
			count+=1
	# print 'The difference between '+ seq1.id +' and ' +  seq2.id + ' is ' + str(count)
	return count

def ShannonEntropy(seqfile,outfile,lenseqs=291):
	# computes the shannon entropy across a list of amino acid sequences
	# SE=sum_x p(x)log(p(x))
	# saves output as a tab delimited file that can be read in Excel.
	# ignores blanks, or uncalled bases.  Requires seqs to be the same length
	# lenseqs, just used to check that they all have the same length
	from numpy import log
	print " ALL SEQUENCES MUST BE THE SAME LENGTH"
	listvals=[{"A":0.,"B":0.,"C":0.,"D":0.,"E":0.,"F":0.,"G":0.,"H":0.,"I":0.,"K":0.,"L":0.,"M":0.,"N":0.,\
			"O":0.,"P":0.,"Q":0.,"R":0.,"S":0.,"T":0.,"U":0.,"V":0.,"W":0.,"Y":0.,"Z":0.} for i in range(lenseqs)]
	# open seqlibrary file
	path="//Users/rbotts/Documents/DATA/Projects/BioPythonTools/"
	inseqs = SeqIO.parse(open(path+seqfile,"rU"),"fasta")
	for record in inseqs:
		# account for duplicates
		if record.id[-1]==")":
			numtimes=int(record.id[record.id.find("(")+1:-1])
		else:
			numtimes=1
		print numtimes
		for i in range(lenseqs):
			if record[i] in listvals[i]:
				listvals[i][record[i]]+=1*numtimes

	SE=[-1*sum([(x/sum(site.values()))*log(x/sum(site.values())) for x in site.values() if x !=0])for site in listvals]
	print SE

	# write tab delimited output:
	import csv
	with open(path+outfile+".csv","wb") as csvfile:
		out=csv.writer(csvfile)
		out.writerow(SE)
		csvfile.close()
	

def listResChanges(simfile='CTX-MSimilarities.csv',seqfile='CompleteLibraryAA.fasta',outfile='CTX-Msubs.csv'):
	#  simfile contains list of clone and most similar CTX-M variant
	#  opens seqs (AA) in seqfile and adds discrepencies to seqfile
   
	refnames=[]
	innames=[]
	refseqs=[]
	inseqs=[]
	
	handle=open(simfile,"rU")
	datareader =csv.reader(handle,delimiter=",")
	for row in datareader:
		refnames.append(row[1])
		innames.append(row[0])
	handle.close()

	seqhandle = open(seqfile,"rU")
	seqs = SeqIO.parse(seqhandle,"fasta")
	# open up seqs and place them into the reference list if they are a potential reference
	for seq in seqs:
		if seq.id.find("CTX-M")==0:
			if seq.id in refnames:
				refseqs.append(seq)
				print seq.id+' added as reference'
		else:
			inseqs.append(seq)
			print seq.id+' added to library'

	with open(outfile,"wb") as csvfile:
		out = csv.writer(csvfile)
		for i in range(len(innames)):
			inname=innames[i]
			refname=refnames[i]
			print inname + ' ' + refname
			# find input seq and reference seq
			for seq in inseqs:
				if seq.id==inname:#)>-1:
					inseq=seq
					break
			for seq in refseqs:
				if seq.id==refname:#)>-1:
					refseq=seq
					break
			r=[inname,refname]
			subs=''
			for j in range(len(refseq.seq)):
				if refseq.seq[j] != inseq.seq[j]:
					# add name of substitution to the list
					# +1 shifts j to 1 through 271 sites, -3 or -2 to count according to naming convention
					if j+1 < 242:
						tempcount = j+1-3
					if j+1>241:
						tempcount = j+1-2
					temp=refseq[j]+str(tempcount)+inseq.seq[j]
				
					#print 'Comparing '+ inname + ' to ' + refname
					
					subs=subs+temp+'_:_' # use concatenation marker to be replaced with , later
					print inseq.id+' has '+temp
			r.append(subs)
			# write output
			out.writerow(r)
		csvfile.close()


def retCTXMfromCSV(infile='CTX-MSeqs2015.csv',outfile='CTX-MSeqsN.fasta'):
	# retrieves CTX-M seqs from a file where the accession numnber is in the second column and the name is in the first
	# file retrieved from Lahey clinic
	# outputs all records in fasta format
	Entrez.email = "rbotts@pointloma.edu"
	seqs=[]
	with open(infile,'rbU') as f:
		reader = csv.reader(f,delimiter=',')
		for row in reader:
			# check that the accession is really given in the table
			if not (row[1] == "Withdrawn" or row[1] == "Assigned" or row[1] == ""):
				print "Fetching "+row[1]
				seqhandle=Entrez.efetch(db="nucleotide", id=row[1], rettype="gb", retmode="text")
				seq = SeqIO.read(seqhandle,"genbank")

				for f in seq.features:
					if f.id.find("CTX-M")>0:
						s = f.extract(seq)
						s.id = row[0]
						seqs.append(s)
	print len(seqs)
	SeqIO.write(seqs, outfile,"fasta")

#######################################
# Routine's from Kristen Peterson
def findPlasmidInfo(email,filename,path):
	
	Entrez.email = email
	list=[]
	with open(filename,'rbU') as f:
		next(f) # skip headings
		reader=csv.reader(f,delimiter='\t')
		for plasmid, num, inc in reader:
			if ((os.path.isfile(path+'plasmid_k_'+inc+'_k_'+plasmid+'.gb'))!=True):
				print plasmid+num+inc
				handle=Entrez.efetch(db='nucleotide',id=num,rettype='gb') # Accession id works, returns genbank format, looks in the 'nucleotide' database
				#store locally
				local_file=open(path+'plasmid_k_'+inc+'_k_'+plasmid+'.gb','wb')
				local_file.write(handle.read())
				handle.close()
				local_file.close()
			else:
				print "done"
				seq=SeqIO.parse(path+'plasmid_k_'+inc+'_k_'+plasmid+'.gb','genbank').next()
				print seq

def findAccs(filename,path,email="rbotts@pointloma.edu"):
	# download accessory genes from file
	# first name in file is the name of the gene, the second is the accession number
	Entrez.email = email
	list=[]
	
	with open(filename,'rbU') as f:
		next(f) # skip headings
		reader=csv.reader(f,delimiter=',')
		for gene, num in reader:
			# skip any genes that have already been downloaded
			if ((os.path.isfile(path+gene+'_k_'+num+'.fasta'))!=True):
				print 'Searching for ' + gene + ' on ' + num
				try:
					handle=Entrez.efetch(db='nucleotide',id=num,rettype='gb') # Accession id works, returns genbank format, looks in the 'nucleotide' database
					#extract feature
				
					# read sequence
					record = SeqIO.read(handle,"genbank")
				
					# extract feature
					if len(record.features)==0:
						outhandle=open(path+gene+'_k_'+num+'.fasta','w')
						print 'No features on this accession, using entire sequence'
						print gene + " found on "+ num
						SeqIO.write(record,outhandle,"fasta")
						outhandle.close()
					else:
						for feature in record.features:
							try:
								if feature.type == "gene":
									#print feature.qualifiers['gene'][0]
									if feature.qualifiers['gene'][0].find(gene)>-1:
										outhandle=open(path+gene+'_k_'+num+'.fasta','w')
										print gene + " found on "+ num
										seq=feature.extract(record)
										seq.id=gene+'_k_'+num
										SeqIO.write(seq,outhandle,"fasta")
										outhandle.close()
										break
							except KeyError:
								print num + ' has no gene features'
				except urllib2.HTTPError, error:
					print 'HTTP error'
					print 'skipping '+gene + 'on ' + num
				handle.close()





##############################################
# Homeless routines #			
def map_annotation_file(seq_file,ann_file):
	# routine reads in seq_file and the corresponding annotation file
	#  annotation file should be .csv
	# 6th column in file is location in the form 123-1840
	# 10th column is name
	# 18th column is the direction of the reading frame.
	# saves image as seq_file.pdf
	#from reportlab.lib import colors
	#from reportlab.lib.units import cm
	#from Bio.Graphics import GenomeDiagram
	#from Bio import SeqIO
	import csv

	outname = seq_file[0:-4]+'.pdf'
	ann_data = open(ann_file,"rU") # rU is for universal line end character


	# read in data
	csv_file = csv.reader(ann_data,delimiter=',',dialect=csv.excel_tab)
	
	# extract info in each row and add create temp seq objects to write to image
	csv_file.next() # skip header row
	for row in csv_file:
		# catch cases where row contains invalid info 
		try:
			seq_name = row[9]
			print seq_name
			seq_loc = row[5]
			tempindex = seq_loc.find('-') # find - and separate name into two addresses
			seq_start = int(seq_loc[:tempindex])
			seq_finish = int(seq_loc[tempindex+1:])
			print seq_start
			print seq_finish

		except ValueError:
			print "A row contains invalid entries, row skipped"

def find_orfs(contig,min_pro_length=200,table=11):
	# currently imports  contig and uses in house orf finder and blasts each separate ors
	# find_orfs_with_trans may not be best technique so we will not develop this further
	# import contig file and identify orf's
	# input: read in contig file, string of basepairs
	# 	min_pro_length - minimum protein sequence length
	#	table - identifies NCBI table to use, default is bacterial sequence
	# output: 
	#	addresses -  list of addresses defining orf's in original sequence
	# Currently only looks for stop codon in protein sequence and returns the 
	# protein sequences with the amino acid sequence addresses
	# Need to modify to work properly, but I am not sure how to best find orfs
	tab=table
	min_pro=min_pro_length

	record = SeqIO.read(contig,"fasta") 
	# read() reads the entire file at once, then parse it if it has multiple records, use open(filename,"ru") to open and 	
	# then SeqIO.parse(filename,"fasta") to read one record at a time
	addresses=find_orfs_with_trans(record.seq,tab,min_pro)
	for start, end, strand, sq in addresses:
		 print "%s...%s - length %i, strand %i, %i:%i" \
			 % (sq[:30], sq[-3:], len(sq), strand, start, end)
	print "Found %i orf's" % len(addresses)

	return addresses
	

def find_orfs_with_trans(seq, trans_table, min_protein_length): 
	# find orf's using translated sequence, returns original sequence
	# too lazy right now to efficiently do this with the original sequence
	# It appears GeneMark does more than identify start and stop codons, so we will use it for now
	# Original Code from BioPython Cookbook Sec 16.1.13 returns nucleotide sequence
	answer = []
	seq_len = len(seq) 
	for strand, nuc in [(+1, seq), (-1, seq.reverse_complement())]:
		for frame in range(3): 
			trans = str(nuc[frame:].translate(trans_table)) 
			trans_len = len(trans) 
			aa_start = 0 
			aa_end = 0 
			while aa_start < trans_len:
				aa_end = trans.find("*", aa_start) 
				if aa_end == -1:
					aa_end = trans_len 
				if aa_end-aa_start >= min_protein_length:
					if strand == 1: 
						start = frame+aa_start*3 
						end = min(seq_len,frame+aa_end*3+3)
					else: 
						start = seq_len-frame-aa_end*3-3 
						end = seq_len-frame-aa_start*3
					answer.append((start, end, strand, str(nuc)[start:end]))
				aa_start = aa_end+1 
	answer.sort()

	return answer
###########################################################
### Tools for ISCR1 analysis

###########################################################
### Tools for ISCR1 analysis
def Find_Sequence(target,outname="ISCR1hits.gb",outextendedname = "ISCR1alignments.fasta",E_VALUE_THRESH=0.0001, MINCOVER = 0.70):
	# open a fasta target sequence, BLAST against Genbank
	# download up and downstream regions of hit and save in output file
	import datetime
	Entrez.email = "rbotts@pointloma.edu"
	blast_type="tblastx" # set type of blast: n, p or x, translate seq record for blastp
	tseqhandle=open(target,"rU")
	targetseq = SeqIO.read(tseqhandle,"fasta")
	print "Opening sequence " + targetseq.id +" length nt %i /AA %i " % (len(targetseq),len(targetseq)/3)
	outhandle = open(outname, "w") # open outfile to store genbank format target sequences.
	outextendedhandle = open(outextendedname,"w")
	
	result_handle = NCBIWWW.qblast(blast_type, "nr", targetseq.seq,hitlist_size=100)

	
	#open file for writing and then save the blast results
	#xml format is chosen because parser is more stable
	save_file = open ("temp_blast.xml", "w")
	save_file.write(result_handle.read())
	save_file.close()
	result_handle.close()

	save_file = open("temp_blast.xml")
	blast_record = NCBIXML.read(save_file)
	
	found_best = False
	count = 0 # count how many hits we are downloading

	hitExtseqs = []# for storing the extended alignments
	aligns = [] # for storing only the aligned region
	#  should also sort for chromosomal vs plasmid, remove integrons

	for alignment in blast_record.alignments:
		for hsp in alignment.hsps:
			#calc max iden and query cover
			#print "query len %s , align length %s" % (len(hsp.query),alignment.length)
			cover = (hsp.align_length*3.0/len(targetseq))
			max_iden = (hsp.identities * 1.0)/len(hsp.query)
			#print ("query cover: %s	max_iden: %s" % (cover, max_iden))
			#save the information about the best match
			if not found_best:
				best_cover=cover
				best_iden=max_iden
				found_best=True

			#collect all other alignments for a particular target that have same "best" scores
			if cover==best_cover and max_iden==best_iden:
				if hsp.sbjct_start < hsp.sbjct_end:
					sbjct_end = hsp.sbjct_end;
					sbjct_start = hsp.sbjct_start
				else:
					sbjct_end = hsp.sbjct_start
					sbjct_start = hsp.sbjct_end
				print "start: %s , end: %s " % (sbjct_start,sbjct_end)
				# check that alignment is high quality and high coverage
				if hsp.expect < E_VALUE_THRESH and cover >  MINCOVER:
					print "e value %s" % hsp.expect
					print "length %s" % alignment.length
					#Extract accession number and description from title
					title_info = alignment.title.split('|')
					accession=title_info[3]
					descript=title_info[4]
					gi_num = title_info[1]
				
					count+=1
					print "Alignment %i" % count
					print title_info
					print ("query cover: %s	max_iden: %s" % (cover, max_iden))
					handle=Entrez.efetch(db='nucleotide',id=gi_num,rettype='gb') # Accession id works, returns genbank format, looks in the 'nucleotide' database
				
					# read sequence
					record = SeqIO.read(handle,"genbank")
					aligns.append(record[sbjct_start:sbjct_end])
					hitExtseqs.append(record[sbjct_start-5000:sbjct_end+5000])
					print record

	print "%s Alignments have been found" % count
	SeqIO.write(hitExtseqs, outhandle,"genbank")
	SeqIO.write(aligns,outextendedhandle,"fasta")
	save_file.close()
	outhandle.close()
	outextendedhandle.close()
	print("----------------------- ORF BLAST COMPLETE ------- %s records found----------" % (str(count)))

def Add_aligned_features(infile = "ISCR1Align.csv",seqfile = "ISCR1extAlign.gb", updatedseqfile = "ISCR1extAlignWAlignFeature.gb"):
	# Ensures that each sequence has the aligned region from the alignment file
	feature = SeqFeature(FeatureLocation(25, 125), strand=+1)
	#gds_features.add_feature(feature, name="Forward", label=True)

def StripCDS(infile="ISCR1hits.gb", outfile = "ISCR1FeaturesAA.fasta"):
	# open all genbank records in infile and write each of the reads to outfile, while maintaining their name
	in_handle = open(infile,"rU")
	inseqs = SeqIO.parse(in_handle,"genbank")
	outhandle = open(outfile,"w")
	k=0 #count the number of records
	m=0 # count the number of features
	for record in inseqs:
		k+=1
		j = 0 # count the number of orfs on each plasmid
		# check that the sequence is a plasmid
		#print record.features
		if True: #record.description.find('plasmid')>-1 and record.description.find('integron')==-1:
			for feature in record.features:
				
				try:
					if feature.type == "CDS":
						m+=1
						j+=1
						print feature.qualifiers
						if 'gene' in feature.qualifiers:
							feat_name=feature.qualifiers['gene'][0].replace(" ","")
						elif 'product' in feature.qualifiers:
							feat_name=feature.qualifiers['product'][0].replace(" ","")
						else:
							feat_name=str(['orfX'+str(j)]).replace(" ","")
						print feat_name + " found on "+ record.id
						sq=feature.extract(record)
						sq.seq=sq.seq.translate(table="Bacterial", to_stop=True)
						org=record.description.split(' ')[0]+'_'+record.description.split(' ')[1] # try to find the organism and add it to the name
						sq.id=feat_name+"-"+str(j)+"_:_"+"ISCR1"+"_:_"+record.id#+'_:_'+org #record.organism.replace(' ','_')
						SeqIO.write(sq,outhandle,"fasta")
			
				except KeyError:
					print record.id + ' has no gene features'
	outhandle.close()
	print "%s sequences in analysis" % k
	print "%s features for analysis" % m

def get_Target_Names(infile="ISCR1extAlign.gb"):
	# get record ids
	in_handle = open(infile,"rU")
	inseqs = SeqIO.parse(in_handle,"genbank")
	return [record.id for record in inseqs]


def open_blast_table(infile = "ISCR1Align.csv", outfile = "ISCRAlign.fasta", extoutfile = "ISCR1extAlign.gb", targetlen = 513):
	# Blast table is available from web based alignment, export .csv
	# opens a blast table and downloads the best alignment for each hit
	# writes all alignments to outfile
	# writes extended regions to extoutfile
	Entrez.email = "ryanbotts@pointloma.edu"
	regionext = 10000 # amount of sequence to gather in both directions from target
	inhandle = open(infile,"rU")
	csv_file = csv.reader(inhandle,delimiter=',',dialect=csv.excel)
	print "opening input file"
	out_handle = open(outfile,"w")
	out_ext_handle = open(extoutfile,"w")
	Accnums=[] # store all accession numbers
	Best_aligns={}
	for row in csv_file:
		try:
			gi_num = row[1].split("|")[3]
			sim = float(row[3])/100
			print row[6]
			align_length = int(row[6])
			print row[11]
			sbjct_start = int(row[11])
			print row[12]
			sbjct_end = int(row[12])
			
			print "Acc Number: %s, Similarity: %s, Align Len: %s " % (gi_num, sim, align_length)
			
			
			# check the alignment is of high enough quality
			if sim >.80 and align_length > targetlen*.90:
				# check if we have found the alignment before, only uses the first, however likely needs improvement.
				if not gi_num in Accnums:
					Accnums.append(gi_num)
					if sbjct_start < sbjct_end:
						true_end = sbjct_end;
						true_start = sbjct_start
					else:
						true_end = sbjct_start
						true_start = sbjct_end
				
				
					handle=Entrez.efetch(db='nucleotide',id=gi_num,rettype='gb')
					record = SeqIO.read(handle,"genbank")
					record.id = gi_num+"_l_%s-%s" % (true_start,true_end)
					shortrec=record
					
					
					shortrec = shortrec[true_start : true_end]
				 

					print "Original record length is: %s" % len(record)
					
					
					# define the sites of the actual start and stop of the target region
					start = true_start - regionext
					end = true_end + regionext
					print "Expected length: %s" % str(end - start)
					if end - start > len(record):
						longrec = record
					else:
						longrec = record[start : end]
						if start < 0:
							longrec= record[start:]+longrec
							print "length of first region %s" % len(longrec)
						if end > len(record):
							longrec = longrec + record[:end-len(record)]
							print "length with end added %s" % len(longrec)
						print "The final lenth is %s" % len(longrec)

					print "Number of features: %s" % len(longrec.features)
					
					if len(longrec.features)>0:
						SeqIO.write(shortrec, out_handle, "fasta")
						SeqIO.write(longrec, out_ext_handle, "genbank")
					print "wrote " + gi_num
	
		except IndexError:
			break
	inhandle.close()
	out_handle.close()
	out_ext_handle.close()
	print "Number of sequences: %s" % len(Accnums)


###------GET CLUSTERS -------####
def get_clusters(clusterfilename): #Dr. B.
	### finds clusters from default tab delimited output from USEARCH
	# creates multiple lists each containing all items in a cluster
	# input:
	#	clusterfilename - tab delimited default output from USEARCH
	#		center of cluster begins with a S in the first column and name in the 9th column
	# output:
	# 	groups - list containing all elements of each cluster
	print "May not correctly extract clusters of singletons"
	groups = []  # empty dictionary to store each cluster
	
	with open(clusterfilename, 'rU') as f:
		reader = csv.reader(f, delimiter='\t')
		notfirst = False # identifies if we are starting at the first group or not
		group=[]
		for row in reader:
			if row[0].strip() == "S" and notfirst:
				# add previous group to dictionary
				groups.append(group)
				# start a new group
				group=[]
			# the rows marked with a C indicate the center of a cluster and may be omitted.
			elif row[0].strip() == "C":
				groups.append(group)
				break
			notfirst = True
			name = row[8].strip() # extract the name from the row
			group.append(name)
		# append the last cluster to the list
		groups.append(group)
	return groups
###----- END GET CLUSTERS ----- ###




###----- GET CLUSTER NAME ------###
def get_clust_name(list): #Dr. B.
	### searches the list of plasmid feature names to find a feature that is not an orf
	# if none is found, returns the name orf
	group_name = list[0]
	
	for name in list:
		if 'orf' not in name.split('_:_')[0]:
			group_name = name
			break
	return group_name


def write_matrix_file(filename='matrix_noDupNames.csv',targetnamefile = 'ISCR1extAlign.gb',clusterfile = 'ISCR1Clusters.tab',writefullnames=False):
	##-----This function creates a matrix of what clusters all the plasmids belong to... so plasmid1 will have a 1 in every column that it is found with the cluster for that column and a 0 in every coumn otherwise
	##--- input: filename: the name of what you want the file to be called
	##		  names: the plasmid names
	##		  groups: the list of the groups from the get_clusters function
	##		  writefullnames = Used to edit how pretty the feature matrix is, this will help to make the heat map pretty
	## *****Be sure to change the cluster file name to what it is suppose to be if you change it else where******
	##--- output: A csv file that contains a matrix of 1s and 0s with the plasmid names down the side and the cluster names at the top.
	row=""
	notEmp=0
	names = get_Target_Names()
	groups =  get_clusters(clusterfile)
	print "Check that the feature matrix has features named according to the convention if plamsid maps are going to be created later."
	with open(filename,'w') as f:
		plasmid_matrix=np.zeros(shape=(len(names),len(groups))) #creates a matrix with a lot of zeros that has the shape of length of plasmid names list by length of cluster groups
		f.write(' ,')
		#writing out the cluster group name at the top of every column
		for c in range(len(groups)):
			if writefullnames:
				f.write(get_clust_name(groups[c])+',')
			else:
				f.write(get_clust_name(groups[c]).split('_:_')[0]+',')
		f.write('\n')
		k=0
		for i in range(len(names)):
			row+=str(names[i])+','
			#f.write(names[i]+',') #writing the plasmid name on the side of every row
			for j in range(len(groups)):
				#plasmid_matrix[i][j]=0
				for gene in groups[j]:
					if gene.split('_:_')[-1].find(names[i])>-1:
						plasmid_matrix[k][j]=1
						notEmp=1
						break
				row+=str(plasmid_matrix[k][j])+','
			if notEmp==1:
				f.write(row)
				k+=1
				f.write('\n')
			
			#f.write(str(plasmid_matrix[i][j])+',') # writing the value into the file
				#if row.find(str(1))>-1:
				#	f.write(row)
			
			notEmp=0
			row=""
		f.close()

def gather_row_seqs(inseqs = "ISCR1Align1.fasta", plasmidmatrix = "ISCR1featurematrix.csv", outseqs = "ISCR1AlignClean.fasta"):
	# make sure that the sequence file only contains sequences which are included in the final feature matrix.
	# tool is used in case feature matrix has been edited by hand
	seq_handle =  open(inseqs,"r")
	
	seqs = SeqIO.parse(seq_handle, "fasta")
	seqlist = {s.id : s for s in seqs}
	
	seq_handle.close()
	
	
	in_handle =  open(feattable,"rU")
	table = csv.reader(in_handle,delimiter=',')
	table.next()
	
	seqnames =[row[0] for row in table]
	
	out = [seqlist[k] for s in seqnames for k in seqlist.keys() if k.find(s)>-1]
	print len(out)
	out_seq_handle = open(outseqs,"w")
	SeqIO.write(out,out_seq_handle, "fasta")
	out_seq_handle.close()

def map_seq_groups(outfilename='ISCR1MapGroup', plasmidftfile = 'ISCR1featurematrix.csv', plasmidclustfile = 'ISCR1Groups.csv', plasmidseqs='ISCR1extAlign.gb',\
	seqclusters = 'ISCR1Clusters.tab', alignfile = 'ISCR1Align.csv', pathtoplasmids='./'):
	# modification of previous version, accounts for having all seqs in one file
	# creates map of all plasmids within a cluster and colors the genes
	# according to the percent of the time the gene occurs in the all plasmids within cluster.
	# inputs:
	#	plasftfile - feature matrix name for each group  All such files will be opened and read in
	#	plasmidgroup - csv file with plasmid name in col1 and group number in col2
	#	plasmidseqs - file with all plasmid genbank files
	#	plasmidmatrix - matrix of plasmids with consensus groups
	#	seqclusters - csv file containing each of the clusters of genes
	#	outfilename - name of file, -i.csv will be appended to the name
	#	pathtoplasmids - path extension to the file with all of the plasmid sequences
	# output:
	#	map for each group of plasmids, written in the folder containing the plasmids
	import csv

	# Read in the clusters of plasmids, create a list of all plasmids in each cluster
	#  the keys for groups are the plasmid cluster numbers
	
	groups={}
	
	print 'Feature names in '+plasmidclustfile+' must follow the convention feature_:_something_:_SeqIdnumber'
	print 'Sequences retrieved by Id number'
	handle = open(plasmidclustfile,'rU')
	tablegroups = csv.reader(handle, delimiter=',',dialect=csv.excel_tab)
	tablegroups.next() # skip header line
	
	
	# col 1 (row[0]) has the plasmid name, col2 (row[1]) is the group number
	for row in tablegroups:
		if row[1] in groups:
			groups[row[1]].append(row[0]) # if the group has been started add the name of the plasmid to the row
		else:
			groups[row[1]]=[row[0]]
		
	handle.close()
	print 'Found %s sequence clusters' % len(groups)

	# identify the location of the ISCR region to mark on the annotation diagram.
	handle = open(alignfile,'rU')
	csvmatrix = csv.reader(handle, delimiter=',',dialect=csv.excel_tab)
	alignloc = {} # stores alignment locations and directions.  key = gi_num, item = feature
	for row in csvmatrix:
		gi_num = row[1].split("|")[3].split(".")[0]
		sim = float(row[3])/100
		print row[6]
		align_length = int(row[6])
		print row[11]
		sbjct_start = int(row[11])
		print row[12]
		sbjct_end = int(row[12])
		if sbjct_start < sbjct_end:
			true_end = sbjct_end;
			true_start = sbjct_start
			feature = SeqFeature(FeatureLocation(true_start, true_end), strand=1)
			feature.id = "ISCR"
			alignloc[gi_num] = feature
		
		else:
			true_end = sbjct_start
			true_start = sbjct_end
			feature = SeqFeature(FeatureLocation(true_start, true_end), strand=-1)
			feature.id = "ISCR"
			alignloc[gi_num] = feature
	handle.close()


	# First identify the percent of sequences in a cluster with particular genes:
	handle = open(plasmidftfile,'rU')
	csvmatrix = csv.reader(handle, delimiter=',',dialect=csv.excel_tab)
	
	# identify the names of the features used in the plasmid matrix
	# these are the first row of the matrx file with the first entry removed.
	names=csvmatrix.next()

	names=names[1:-1] #there is an extra blank entry in the matrix
	print 'The total number of protein families used in analysis '+str(len(names))


	# read the entries in the feature matrix
	matrix={}
	for row in csvmatrix:
		matrix[row[0]]=[float(x) for x in row[1:-1]] # read in row and convert strings to ints, the -1 is used as there is an empty column at the end of each row.
	handle.close()
	#print matrix.keys()

	# get all of the groups of features and place them in a dictionary with the feature names as keys
	# watch that naming of clusters is unique and that the names can be found
	# genelist is used to search for a feature on a sequence and identify the gene cluster name used in the plasmid matrix
	gene_groups = get_clusters(seqclusters)

	#create dictionary of clusters so that we can cross reference the name of the feature on the initial plasmid to the name of the cluster
	# names should correspond to names used in the feature matrix
	genelist={get_clust_name(clust).split('_:_')[0]:clust for clust in gene_groups}


	# find all plasmids in each group, draw the plasmid,
	#  write all features that are in the list, and color them according to the percent used by the group

	# first obtain list of plasmid files

	#filenames=os.listdir(pathtoplasmids)
	# read in all seqs for plotting, not efficient, but works for small batches of seqs
	handle = open(pathtoplasmids+plasmidseqs,'rU')
	seqs = {s.id : s for s in SeqIO.parse(handle,"genbank")}

	for grp in groups.keys():
		print "Starting group " + grp
		name = "Groups"
		gd_diagram = GenomeDiagram.Diagram(name)
		max_len=0 # store the maximum length of all tracks in a file
		k = 0 # keep track of the the track number
		
	
		rowsfromgroups=[]
		
		for plasmid in groups[grp]:
			try:
				rowsfromgroups.append(matrix[plasmid])
			except KeyError:
				print 'Sequence '+plasmid+' not found'

		genepercents = np.mean(rowsfromgroups,axis=0)
		#print genepercents
		# print len(genepercents)
		# print len(names)

		nonzeroGroups = {names[i]: genepercents[i] for i in range(len(names)) if genepercents[i]>0}
		# print nonzeroGroups
		#print nonzeroGroups
		# open the genbank record file for each plasmid in the group
		for plasmid in groups[grp]:

			# find all nonzero features on plasmid, record the name and the percent in group with gene
			grpfeatplasmid = {}
			for feat in nonzeroGroups.keys():
				# keys for genelist are short names and not full names, so we must first find the group containing the feature name
				# Find the key for the group containing the feature
				genekeys = [genekey for genekey in genelist if any([feat in gene for gene in genelist[genekey]])]
				#print feat
				#print genekeys
				#print genelist
				try:
					# use the first key because the second seems to be some strange artifact in the clustering
					for gene in genelist[genekeys[0]]:
						# check if the last portion of the name of the gene is our plasmid name
						if gene.split('_:_')[-1] == plasmid:
							#print gene
							# note that this adds the name of the feat to the key
							print gene.split('_:_')[0]+ " found on "+plasmid
							# If the feature is found on the plasmid assign it a value for how often it is found.
							grpfeatplasmid[gene.split('_:_')[0]] = nonzeroGroups[feat]
							#grpfeatplasmid[gene.split('_:_')[0].split('-')[0]+'_:_'+feat.split('-')[0]] = nonzeroGroups[feat]
				except IndexError:
					print feat + " not found in seqs, must have been removed from the cluster table"
			record = seqs[plasmid]
			# reverse the sequence if the alignement is on the reverse complement

			if alignloc[plasmid].strand == -1:
				print record.id + " has %s features" % len(record.features)
				record = record.reverse_complement()
				print record.id+" Reversed, has %s features" % len(record.features)
				# set the feature to now be in the correct orientation
				alignloc[plasmid].strand = 1
				
			max_len = max(max_len, len(record))
			gd_track_for_features = gd_diagram.new_track(k,
												 name=plasmid,
												 greytrack=True,
												 start=0, end=len(record),
												 height=1.0,smallticklabels=1)
			gd_feature_set = gd_track_for_features.new_set()
			k+=1 # keep track of the feature count
			j = 0
			print plasmid
			# get the plasmid features in the group that are nonzero
			for feature in record.features:
				if feature.type=="CDS":
					j+=1
					# name feature appropriately
					if 'gene' in feature.qualifiers:
						feat_name=feature.qualifiers['gene'][0].replace(" ","")
					elif 'product' in feature.qualifiers:
						feat_name=feature.qualifiers['product'][0].replace(" ","")
					else:
						feat_name=str(['orfX'+str(j)]).replace(" ","")
					feature.id = feat_name+"-"+str(j)
					
					# only draw features that are common to the group
					for grpfeat in grpfeatplasmid.keys():
						try:
							if feature.id in grpfeat:
								col=colors.linearlyInterpolatedColor(colors.white, colors.blue, 0, 1, grpfeatplasmid[grpfeat])
								gd_feature_set.add_feature(feature, sigil="BIGARROW",
								   arrowshaft_height=.8,
								   arrowhead_length=.25,
								   color=col, border=col,label=True,
								   name = feature.id,
								   label_position="center",
								   label_size = 7, label_angle=45)
							else:
								print feature
								print '%s not found on %s, skipping' % (feature.id, plasmid)
						except KeyError:
							print 'CDS does not have gene name on '+plasmid
			# add ISCR region to the map and
			alignft = alignloc[plasmid]
			gd_feature_set.add_feature(alignft, sigil="BIGARROW",
								   arrowshaft_height=.8,
								   arrowhead_length=.25,
								   color=colors.green, border=colors.green,label=True,
								   name = "ISCR1",
								   label_position="center",
								   label_size = 7, label_angle=45)

		#width = max(1500,max_length/100)
		height= len(groups[grp])*100
		gd_diagram.draw(format="linear", pagesize=(2000,height), fragments=1,
				 track_size =.2,start=0)
		gd_diagram.write(pathtoplasmids + outfilename +'-'+ grp+ ".pdf", "PDF")

##########################################################


###########################################################

if __name__=="__main__":

	## SR Tools
	#seq = "Temp"
	#Blast_Features(seq)
	#extract_features_tab(seq)
	
	## Honors tools
	seq = 'PC1concat'
	database = 'SeqsByCluster.udb'
	Usearch_Features(seq,database)
	#extract_features_tab(seq)
	##For database creation place the desired database in fasta format into Create_Usearch_Database
	#Create_Usearch_Database('ResistanceGenes.fasta')
	#Create_Usearch_Database('SeqsByCluster.fasta')
	#Create_Usearch_Database('SmallTest1.fasta')
	folder_directory = '/Users/lucasustick/SR/CrossLink_Honors/SequenceToolsH/SeqsByCluster'
	outfile='SeqsByCluster.fasta'
	#composite_fasta(folder_directory,outfile)
	#Translate_Fasta('ResistanceGenes.fasta','ResistanceGenesAA.fasta')
	

	
	#Blast_Features("T11concat")
	#GroupSeqs("//Volumes/Users/rbotts/DATA/Projects/BioPythonTools/CTXM-Clones/")
	#elim_Genomic_Contigs_BLAST('8PlasmidContigs2015.fasta')
	#elim_Genomic_Contigs_BLAST('pTRE.J53-9')
	#elim_Genomic_Contigs_BLAST('pTRE.J53-27')
	#elim_Genomic_Contigs_BLAST('pTRE.J53-7')
	#elim_Genomic_Contigs_BLAST('pTRE.PP5')
	#elim_Genomic_Contigs_BLAST('pTRE.P11')
	#ShannonEntropy('TRE-1LibraryAA.fasta','TRESE')
	#FindIdenticalSeqs("CompleteLibrary.fasta","CompleteLibraryUnique.fasta","DuplicateTable.csv")
	#listResChanges()
	#findAccs('ARGAccNumbers.csv','//Users/rbotts/Documents/DATA/Projects/BioPythonTools/Accessories/')
	#retCTXMfromCSV()
	
	
	
	### tools for downloading ISCR sequences
	#Find_Sequence("M1-ISCR1.fasta")
	#open_blast_table(infile = "ISCR1Align.csv", outfile = "ISCRAlignTest.fasta", extoutfile = "ISCR1extAlign.gb", targetlen = 514)
	# note that if usearch fails for empty file it is likely due to having written a duplicate sequence in the AA file, remove sequence
	#StripCDS(infile="TestAlign.gb", outfile = "TestAlign.fasta")
	#StripCDS(infile="ISCR1extAlign.gb", outfile = "ISCR1FeaturesAA.fasta")
	#os.system('./usearch -cluster_fast ISCR1FeaturesAA.fasta -id 0.7 -target_cov 0.7 -centroids clust.fasta -uc ISCR1Clusters.tab -userfields query+target+id+ql+tl+alnlen+qcov+tcov')
	#write_matrix_file('ISCR1featurematrix.csv',targetnamefile = 'ISCR1extAlign.gb',clusterfile = 'ISCR1Clusters.tab')
	#map_seq_groups(outfilename='ISCR1MapGroup', plasmidftfile = 'ISCR1featurematrix.csv', plasmidclustfile = 'ISCR1Groups.csv', plasmidseqs='ISCR1extAlign.gb',\
	#seqclusters = 'ISCR1Clusters.tab', pathtoplasmids='./')
	# prior to running matrices through heat map, removed all proteins found in 2 or fewer sequences.

	#gather_row_seqs()