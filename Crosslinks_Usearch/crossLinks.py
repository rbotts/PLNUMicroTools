from Bio.Blast.Applications import NcbiblastpCommandline
from StringIO import StringIO
from Bio.Blast import NCBIXML
from Bio.Blast import NCBIWWW
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Graphics import GenomeDiagram
from Bio.Graphics.GenomeDiagram import CrossLink
from reportlab.lib import colors
from Bio.SeqFeature import FeatureLocation
from Bio import SeqFeature
from Bio.SeqFeature import SeqFeature
from reportlab.lib.colors import red, grey, orange, green, brown, blue, lightblue, purple, yellow
import re
import csv
from Bio.Blast.Applications import NcbiblastpCommandline
from StringIO import StringIO
import sys
import os

def ublastfeatures(features1,features2):
	## Make database
	#outname = raw_input("Enter the name you would lke the Ublast output to be under, do not include extension\n")
	outname = features1[:-6]+'_'+features2[:-6]
	outfile = outname+'.tab'
	print 'creating usearch database based on the first sequence'
	os.system('./usearch -makeudb_ublast '+features1+' -output '+features1[:-6]+'.udb')
	## Do the allignment
	print 'running usearch comparison'
	os.system('./usearch -ublast '+features2+' -db '+features1[:-6]+'.udb -evalue 0.04 -blast6out '+outfile+' -strand plus -top_hit_only') 
	## make the .tab file into a csv
	print 'Converting Usearch output to .csv from tab delimited file'
	csv_w = csv.writer(open(outfile[:-4]+'.csv', 'w'))
	for row in csv.reader(open(outfile, 'r'), delimiter = '\t'):
		csv_w.writerow(row)
	#csv_w.close()
	os.remove(outfile)
	return outfile[:-4]+'.csv'
	return outfile

def writeFasta(Filename, FastaFileName):
	#creates a fasta file with each feature as a its own sequence to be BLAST aganst
	try:
		handle1=open(Filename,"rU")
		record1 = SeqIO.read(handle1, "gb")
		outhandle1 = open(FastaFileName,'w')
		for feature in record1.features:
			if feature.type=='CDS':
				proteinseq = feature.extract(record1)
				try:
					proteinseq.id=feature.qualifiers['protein_id'][0]
				except KeyError:#if no protein id to use
					try:
						print 'protein_id qualifier not found. Using gene qualifier'
						print feature.qualifiers['gene'][0]
						proteinseq.id=feature.qualifiers['gene'][0]
					except KeyError:#if no gene name
						try:
							print 'gene qualifier not found using product qualifier'
							print feature.qualifiers['product'][0]
							proteinseq.id=feature.qualifiers['product'][0]
						except KeyError:#if no product name
							proteinseq.id='No_id'
							print "qualifier not found. Using No_id as identifier"
				SeqIO.write(proteinseq,outhandle1,"fasta")
		outhandle1.close()
	except IOError:
		print "file not found please try again"
	
def get_feature(features, id, tags=["protein_id","product"]):

	"""Search list of SeqFeature objects for an identifier under the given tags."""
	for f in features:
		for key in tags:
			#tag may not be present in this feature
			for x in f.qualifiers.get(key, []):
				#if x == id and f.qualifiers.get("gene") != ['parA']:
				if x == id:
					return f
	raise KeyError(id)

def getCrossLinks(filename):
	A_v_B = [(1,"gg", "gg" )]#structure crosslink to start off (percentmatch, "gene1", "gene2")
	with open(filename) as f:
	# read the file line by line
		for line in f:
			try:
				words = line.split(',')
				A_v_B.append((float(words[2]), words[0] ,words[1]))
			except IndexError:
				print "empty line"
	return A_v_B

def ask_for_plasmid():
	records = 0
	while records==0:
		plasmid=raw_input("Please enter a file to be added to the diagram, include extension\n")
		try:
			if plasmid[-3:]==".gb":
				records = SeqIO.read(plasmid,"gb")
			elif plasmid[-6:]==".fasta": 
				records = SeqIO.read(plasmid,"fasta")
			else:
				print "extension needs to be .gb or .fasta"
		except IOError:
			print "file not found try again\n---------------"

	return records
	
def crosslinks(fileName,GenBank_1,GenBank_2):
	gd_diagram = GenomeDiagram.Diagram(fileName)
	max_len = 0
	
	#Open Files and create fasta files to be compared by Ublast
	A_rec = SeqIO.read(GenBank_1,'gb')
	GB_file_name = GenBank_1
	fasta_file_name_A=GB_file_name[:-3]+".fasta"
	writeFasta(GB_file_name,fasta_file_name_A)
	B_rec = SeqIO.read(GenBank_2,'gb')
	GB_file_name = GenBank_2
	fasta_file_name_B=GB_file_name[:-3]+".fasta"
	writeFasta(GB_file_name,fasta_file_name_B)
	
	
	#create the tab file with the Ublast output
	blastfile = ublastfeatures(fasta_file_name_A,fasta_file_name_B)
	
	Gname='nn'#name of gene to add
	
	#First section gets the crosslinks from the blast files
	A_vs_B = getCrossLinks(blastfile)
	
	#print ('(percent, Gene Query, Gene result)')#This prints the list of Blast results for reference
	for item in  A_vs_B:
		print item
	# asks user for a gene name to highlight 
	gene_search = raw_input("would you like to highlight a specific gene name?\n \t1) Yes\n \t2) No\n")
	if gene_search=="1" or gene_search.lower() =="yes":
		gene_highlight = raw_input("What is the name of the gene you would like to highlight?\n")
	
		print gene_highlight +" will be printed in red on the genome diagram, all other genes will be grey"
		C_colors =  [yellow]*1+[orange]*1+[brown]*1+[lightblue]*1+[purple]*1+[green]*1+[grey]*1
	else:
		gene_highlight = "NONE"
		C_colors = [yellow]*1+[orange]*1+[brown]*1+[lightblue]*1+[purple]*1+[green]*1+[grey]*1#this creates an array of color for the arrows in the GUI
	i = 0 #index of random color to add
	
	geneColor = grey#color of gene. Grey= no name
	
	# Create new features for concatenations
	recs = ("A","B")
	for rec in recs:
		if rec == "A":
			for loc_a in re.finditer('NNNNNCACACACTTAATTAATTAAGTGTGTGNNNNN', str(A_rec.seq)):
				concat_feature = SeqFeature(FeatureLocation(loc_a.start(),loc_a.start()+35,strand=-1), id="Concat",type="CDS",qualifiers={'product':'Concat'})
				A_rec.features.append(concat_feature)
		else:
			for loc_b in re.finditer('NNNNNCACACACTTAATTAATTAAGTGTGTGNNNNN', str(B_rec.seq)):
				concat_feature = SeqFeature(FeatureLocation(loc_b.start(),loc_b.start()+35,strand=-1), id="Concat",type="CDS",qualifiers={'product':'Concat'})
				B_rec.features.append(concat_feature)
	
	#Read in lists of gene names and types
	with open('Backbones_2_Clean.csv','r') as hand1:
		back_b = csv.reader(hand1)
		backbone = list(back_b)
	with open('AntibioticResistanceGenesClean.csv','r') as hand2:
		AnRe = csv.reader(hand2)
		An_Re = list(AnRe)
	
	
	#this loop adds each gene feature to the record with a color and name
	for record, gene_colors in zip([A_rec, B_rec], [C_colors, C_colors]):
	
		max_len = max(max_len, len(record))
		gd_track_for_features = gd_diagram.new_track(1,
								name=record.name,
								greytrack=True,
								start=0, end=len(record))
		gd_feature_set = gd_track_for_features.new_set()
	
		for feature in record.features:
			if feature.type != "CDS":
				#Exclude this feature
				continue
			## Chose Colors of annotations based on gene name
			try:
				Gname=feature.qualifiers['product'][0]
				if Gname == gene_highlight:
					geneColor= red
				elif gene_highlight == "NONE":
					geneColor=gene_colors[i%6]
				# Backbone genes
				elif Gname in backbone[0]:
					geneColor= blue
				#Antibiotic Resistance Genes
				elif Gname in An_Re[0]:
					geneColor= green
				#Transposease & Intergrase
				elif Gname == 'Tnp' or Gname == 'Int':
					geneColor= orange
				else:
					geneColor=grey
			except KeyError:#if no gene name make it grey
				Gname ='No Name'
				geneColor = grey;
			gd_feature_set.add_feature(feature, sigil="BIGARROW",#this adds gene features to gd_feature_set
									   arrowhead_length = .25,
									   color = geneColor, label = True,
									   name = Gname,
									   label_position = "start",
									   label_size = 6, label_angle = 45)
			i+=1#increment i so that arrows will have a random color
	
	track_X = gd_diagram.tracks[2]
	track_Y = gd_diagram.tracks[1]
	
	#this loop adds the cross links so they point to their feature in the diagram
	for score, id_X, id_Y in A_vs_B:
		try:
			feature_X = get_feature(A_rec.features, id_X)
			feature_Y = get_feature(B_rec.features, id_Y)
			color = colors.linearlyInterpolatedColor(colors.white, colors.firebrick, 0, 100, score)
			link_xy = CrossLink((track_X, feature_X.location.start, feature_X.location.end),
							 (track_Y, feature_Y.location.start, feature_Y.location.end),
							 color, colors.lightgrey)
			print "Link made	"
			gd_diagram.cross_track_links.append(link_xy)
		except KeyError:
			print "Feature qualifier for crosslink not found"# for those pesky nameless genes
	
	gd_diagram.draw(format="linear", pagesize=(1200,2400), fragments=1,start=0, end=max_len)
	print max_len
	gd_diagram.write(fileName + ".pdf", "PDF")


	
if __name__=="__main__":

	# XXXXXXXXXXXXXXXXXXXXXXXXXX #
	#        Main Methods        #
	# XXXXXXXXXXXXXXXXXXXXXXXXXX #
	
	## Name of the image file that will be produced at the end without any extension
	filename = 'graph1_10'
	## First sequence in GenBank Format with .gb extension
	Seq1= 'PC1Concat.gb'
	## Second sequence in GenBank Format with .gb extension
	Seq2= 'PC10Concat.gb'
	
	## Remove comments to be prompted to put in input
	# fileName = raw_input("enter the name of the output file\n")
	# Seq1 = ask_for_plasmid()
	# Seq2 = ask_for_plasmid()
	
	crosslinks(filename,Seq1,Seq2)
	