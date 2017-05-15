# Assembly tools to automate assemblies of clean reads in a given folder.
# Uses SPAdes
# Assumes reads have not been merged prior to running and uses BayesHammer, which is biult into SPAdes
import os
import subprocess
import re
import Bio
from Bio import SeqIO
import time


def FolderStructure(path,cleanreadfname, outreadfname):
	#  Check structure of current folder and create new folder for storing reads
	#  Looks for all folders in 02-Cleaned
	#  Creates new folder for assemblies called '03-Complete'
	#  Creates corresponding folder for each sample in '03-Complete'
	if not os.path.exists(path + '/' + outreadfname):
		os.makedirs(path + '/' + outreadfname)
		names = os.listdir(path + '/'+cleanreadfname)
		for f in names:
			os.makedirs(path + '/' + outreadfname + '/' + f)
	else:
		print 'Check that '+ outreadfname +' does not contain assembled reads and move prior to assembly'

def AssembleEach(path, cleanreadfname = '02-CleanedReady',outreadfname = '03-Assembled'):
	# assembles paired and unpaired reads in each folder in cleanreadfname
	# assumes that there is only one file with forward/reverse/unpaired reads
	FolderStructure(path,cleanreadfname,outreadfname)
	names = [ d for d in os.listdir(path + '/'+cleanreadfname) if os.path.isdir(path + '/'+cleanreadfname+'/'+d)]
	pathtosample = path + '/'+cleanreadfname + '/'
	print 'Forward reads named **R1.fastq, Reverse **R2.fastq, and unpaired **SE.fastq'
	for f in names:
		outdirname = path + '/' + outreadfname + '/' + f
		if not os.path.exists(outdirname+'/contigs.fasta'): 
			# get the forward, reverse and unpaired read file names
			for fi in os.listdir(pathtosample + f):
				if 'R1.fastq' in fi:
					forwardseqfname = pathtosample + f +'/' +fi
				elif 'R2.fastq' in fi: # SHOULD THIS BE AN ELIF?
					reverseseqfname = pathtosample + f +'/' +fi
				elif 'SE.fastq' in fi:
					unseqfname = pathtosample + f +'/' +fi 
		
			print 'Assembling sample ' + f
			subprocess.call(['spades.py','-k', '21,33,55,77', '--careful', '-1', forwardseqfname, '-2', reverseseqfname,\
				'-s', unseqfname, '-o', outdirname])
			time.sleep(600)
		else:
			print f+'sample already assembled'

def AssembleEach_Plasmid_Spades(path, cleanreadfname = '02-CleanedReady',outreadfname = '03-Assembled'):
	# assembles paired and unpaired reads in each folder in cleanreadfname
	# assumes that there is only one file with forward/reverse/unpaired reads
	FolderStructure(path,cleanreadfname,outreadfname)
	names = [ d for d in os.listdir(path + '/'+cleanreadfname) if os.path.isdir(path + '/'+cleanreadfname+'/'+d)]
	pathtosample = path + '/'+cleanreadfname + '/'
	print 'Forward reads named **R1.fastq, Reverse **R2.fastq, and unpaired **SE.fastq'
	for f in names:
		outdirname = path + '/' + outreadfname + '/' + f
		if f != '.DS_Store':
			if not os.path.exists(outdirname+'/contigs.fasta'): 
				# get the forward, reverse and unpaired read file names
				for fi in os.listdir(pathtosample + f):
					if 'R1.fastq' in fi:
						forwardseqfname = pathtosample + f +'/' +fi
					elif 'R2.fastq' in fi: # SHOULD THIS BE AN ELIF?
						reverseseqfname = pathtosample + f +'/' +fi
					elif 'SE.fastq' in fi:
						unseqfname = pathtosample + f +'/' +fi 
		
				print 'Assembling sample ' + f
				subprocess.call(['./SPAdes-3.8.0-Darwin/bin/plasmidspades.py','-k', '21,33,55,77', '--careful', '-1', forwardseqfname, '-2', reverseseqfname,\
						 	'-s', unseqfname, '-o', outdirname])
				time.sleep(600)
			else:
				print f+'sample already assembled'

def ReformatFastQ(path):
	# Find all sample folders in path, within these folders reformat the F and R fastq files to have the @ character at the
	# start of the line

	# Assumes a subfolder for each sample, finds all subfolders
	fnames = [ d for d in os.listdir(path) if os.path.isdir(path + '/' + d)]
	
	# within each subfolder find the forward and reverse files
	for f in fnames:
		#os.makedirs(path+'/'+f+'out')
		for fi in [ d for d in os.listdir(path+'/'+f) if (d.find("_R2.fastq")>0 or d.find("_R1.fastq")>0)]:
			pathtoin = path+'/'+f+'/'+fi
			pathtoout = path+'/'+f+'/'+fi+'2'
			print pathtoin
			print pathtoout
			file = open(pathtoin,'r')
			fout = open (pathtoout,'w')
			for line in file:
				replaced = re.sub(r'(^WI)|(^HWI)','@NHWI',line)
				fout.write(replaced);
			file.close()
			fout.close()
			# remove the old file and replace it with the new one
			os.remove(pathtoin)
			os.rename(pathtoout, pathtoin)
		#os.remove(path+'/'+f+'out')

def RestrictContigsByLength(minlen = 5000, maxlen = 150000, contigfname = 'contigs.fasta', \
							infolder = '03-Assembled/',outfolder = '04-RestrictedContigs/'):
	# Finds all contigs in the file
	try:
		os.makedirs('./' + outfolder)
	except OSError:
		os.makedirs('./' + outfolder + '-2')
	fnames = [ d for d in os.listdir('./'+infolder) if os.path.isdir('./' + infolder)]
	for f in fnames:
		# Keep track of the number of contigs that meet our criteria for each sample
		initialcounts = 0
		qualifyingcounts = 0
		try:
			incontigf = open('./'+infolder+f+'/'+contigfname,'rU')
			outcontigf = open('./'+outfolder+'/'+f+'-contigs.fasta', 'w')
			for s in SeqIO.parse(incontigf,"fasta"):
				initialcounts+=1
				if len(s)<maxlen and len(s) > minlen:
					s.id=f+'-contig-'+str(initialcounts)+'.fasta'
					SeqIO.write(s,outcontigf,'fasta')
					qualifyingcounts+= 1
			incontigf.close()
			outcontigf.close()
			print 'Sample %s had %i contigs, of which %i were between %i and %i bp long', (f,initialcounts,qualifyingcounts, minlen,maxlen)

				
		except IOError:
			print "Contig file not found in sample %s", f
			
def Bowtie2_align(reference,path,cleanreadfname,outreadfname):
	#Call Bowtie2
	#--un <path> Write unpaired reads that fail to align to file at <path>.
	#--un-conc <path> Write paired-end reads that fail to align concordantly to file(s) at <path>.
	# -x reference index
	# -1 -2 comma-separated list of files containing paired end reads to be aligned
	# -U Comma-separated list of files containing unpaired reads to be aligned
	# -S File to write SAM alignments to.
	FolderStructure(path,cleanreadfname,outreadfname)
	names = [ d for d in os.listdir(path + '/'+cleanreadfname) if os.path.isdir(path + '/'+cleanreadfname+'/'+d)]
	pathtosample = path + '/'+cleanreadfname + '/'
	print 'Forward reads named **R1.fastq, Reverse **R2.fastq, and unpaired **SE.fastq'
	for f in names:
		outdirname = path + '/' + outreadfname + '/' + f
		if not (os.path.exists(outdirname+'/'+f+'_R1.fastq') and os.path.exists(outdirname+'/'+f+'_R2.fastq') and os.path.exists(outdirname+'/'+f+'_SE.fastq')):
			#get the forward, reverse and unpaired read file names
			for fi in os.listdir(pathtosample + f):
				if 'R1.fastq' in fi:
					forwardseqfname = pathtosample + f +'/' +fi
				elif 'R2.fastq' in fi: # SHOULD THIS BE AN ELIF?
					reverseseqfname = pathtosample + f +'/' +fi
				elif 'SE.fastq' in fi:
					unseqfname = pathtosample + f +'/' +fi
			outdirname = path + '/' + outreadfname + '/' + f
			print 'Aligning Sample ' + f
			#Run Bowtie2 on paired ends
			os.system('bowtie2 -x '+reference+' -1 '+forwardseqfname+' -2 '+reverseseqfname+' --un-conc '+outdirname+'/'+f+'_R%.fastq')
			#Run Bowtie2 on unpaired ends
			os.system('bowtie2 -x '+reference+' -U '+unseqfname+' --un '+outdirname+'/'+f+'_SE.fastq')
		else:
			print 'Sample '+f+' Already Aligned'
	
def Bowtie2_Build_reference(reference,outname):
	# Place file in directory you want the reference built in, with the fasta file
	os.system('bowtie2-build '+reference+' '+outname)

#build reference
#Bowtie2_Build_reference('Combined_Reference.fasta','Reference')

#ReformatFastQ('./02-Cleaned') # this line adds the @ character to the fastq file because some were intially missing it
Bowtie2_align('Reference',path='./',cleanreadfname='01-Test',outreadfname='02-Bowtie2_Cleaned')
AssembleEach_Plasmid_Spades('./',cleanreadfname='02-Bowtie2_Cleaned',outreadfname='03-Assembled_All_Bowtie2')

#AssembleEach('./',cleanreadfname='02-Cleaned',outreadfname='03-Test')
#RestrictContigsByLength(maxlen = 250000, infolder='03-Assembled_All_Bowtie2/',outfolder='04-RestrictedContigs_All_Bowtie2')

#Test alignment

