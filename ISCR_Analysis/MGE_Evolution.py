
from Bio.Seq import Seq, MutableSeq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature,FeatureLocation
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.SeqUtils import GC
from Bio import SeqIO, Entrez, Alphabet
from Bio.Graphics import GenomeDiagram
from Bio.Graphics.GenomeDiagram import CrossLink
from reportlab.lib.colors import *
from reportlab.lib import colors
import urllib2
import os
import csv
import numpy as np
import glob
from pandas.io.parsers import read_csv

#from user_plasmids import *
#from expanding import *
#from copy import deepcopy

###########################################################


### Tools for downloading and retrieving sequences

def retrieve_from_blast_table(infile = "ISCRHitTable.csv", compseqout = "ISCRCompleteAlignTest.fasta", outfile = "ISCRAlignTest.fasta", extoutfile = "ISCRTestExtAlign.gb", targetlen = 513, featname = "ISCR", down_only = False):
    # Blast table is available from web based alignment, export .csv
    # opens a blast table and downloads the best alignment for each hit
    # Inputs:
    #   infile - name of file for BLAST output in csv format
    #   compseqout - name of file for complete sequences of the alignment, useful for building phyogenetic trees, etc.
    #   outfile - name of file for sequences of the alignment region in fasta format, useful for building phyogenetic trees, etc.
    #   extoutfile - name of file for outputting the alignment regions
    #   targetlength - nucleotide length of the the target sequence, used to compute max percent coverage
    #   featname - name of the alignment region, used when adding a feature to the alignment region to ensure it is annotated on the sequence
    #
    # Retrieved sequences must meet the following criteria:
    #   Minlength = length of the extension region or twice its length if extending in both directions
    #   min_per_cov = minimum allowable percent coverage (hard coded)
    #   min_per_id =  minimum allowable percient identity in alignment (hard coded)
    #   region must containg annotated features
    #
    # The following are output:
    #   genbank file with extended window around alignment
    #   fasta file of alignment regions
    #   fasta file of complete records of aligned sequences
    Entrez.email = "ryanbotts@pointloma.edu"
    regionext = 10000 # length of sequence to retrieve beyond alignment
    min_per_cov = .70
    min_per_id = .30
    print 'Only using alignments greater than %.2f id and %.2f coverage.' % (min_per_id, min_per_cov)
    if down_only:
        print "Only extending in one direction, assumes protein used in search is pointing forward."
        total_len = regionext
        len_up = regionext
        len_down = 0 #regionext
    else:
        print "extracting up and downstream region"
        total_len = 2*regionext
        len_up = regionext
        len_down = regionext

    inhandle = open(infile,"rU")
    csv_file = csv.reader(inhandle,delimiter=',',dialect=csv.excel)
    
    save_complete_seqs = True # save the complete record for the alignment seq
    out_handle = open(outfile,"w")
    if save_complete_seqs:
        comp_out_handle = open(compseqout,"w")
    out_ext_handle = open(extoutfile,"w")
    Accnums=[] # store all accession numbers
    Best_aligns={}
    countseqs = 0 # keep track of how many sequences meet criteria

    for row in csv_file:
        try:
            gi_num = row[1].split("|")[3]
            id = float(row[3])/100
            print row[6]
            align_length = int(row[6])
            print row[11]
            sbjct_start = int(row[11])
            print row[12]
            sbjct_end = int(row[12])
            e_value = float(row[13])
            
            print "Acc Number: %s, Identity: %s, Align Len: %s " % (gi_num, id, align_length)
            
            #determine orientiation of alignment, for uniform ordering and determine the appropriate start and stop coordinates
            if sbjct_start > sbjct_end:
                rev_align = True
                true_end = sbjct_start
                true_start = sbjct_end-1
                len_left = len_down
                len_right = len_up
                feature = SeqFeature(FeatureLocation(true_start, true_end), strand=-1, type ='CDS', id = featname,qualifiers={"locus_tag":featname})
            else:
                rev_align = False
                true_end = sbjct_end;
                true_start = sbjct_start-1
                len_left = len_up
                len_right = len_down
                feature = SeqFeature(FeatureLocation(true_start, true_end), strand=1, type ='CDS', id = featname,qualifiers={"locus_tag":featname})

 
            # compute end of extended region
            start = true_start - len_left
            end = true_end + len_right
            
            
            # check the alignment is of high enough quality
            if e_value < .0001 and align_length > targetlen * min_per_cov:
                # link to record
                handle=Entrez.efetch(db='nucleotide',id=gi_num,rettype='gb')
                record = SeqIO.read(handle,"genbank")
                #print record.description
                if record.description.find('plasmid')>0 or record.description.find('Plasmid')>0:
                    print "Found plasmid"
                    record.id = gi_num+'_l_%s-%s_pl' % (true_start,true_end)
                    #record.name = record.id # should be the record name, but this is too long
                    record.description = record.id + ' ' + record.description
                else:
                    record.id = gi_num+'_l_%s-%s' % (true_start,true_end)
                    #record.name = record.id
                    record.description = record.id + ' ' + record.description
                #print record.id
                #record.id = gi_num+"_l_%s-%s" % (true_start,true_end)
                
                record.features = record.features + [feature]
                
                # check that the record is long enough
                if len(record) > total_len: #and (record.description.find("complete")> 0):
                    # order all of the aligned regions so that they are in the appropriate orientation
                    if rev_align:
                        shortrec = record[true_start : true_end].reverse_complement()
                    else:
                        shortrec = record[true_start : true_end]
                    
                    print "Original record length is: %s Window length: %s" % (str(len(record)),str(end - start))
                    
                    # extract the appropriate region if the extended window wraps around the end of the sequence or extends past the end, assumes circular DNA
                    #if start < 0:
                    #    longrec= record[start:]+record[:end]
                    #if end > len(record):
                    #    longrec = record[start:] + record[:end-len(record)]
                    #else:
                    #    longrec = record[start : end]
                    longrec = record[ max(0,start): min(len(record),end)]
                    longrec.description = record.description
                    shortrec.id = record.id
                    shortrec.description = ''
                    print "The final lenth is %s" % len(longrec)
                    print "Number of features in align region: %s" % len(longrec.features)
                    
                    # check that there aren't a bunch of NNNN's inalignment region, not sure why this would happen
                    if len(longrec.features)>5 and longrec.seq.find("NNNNNN")<0 and shortrec.seq.find("NNNNNN")< 0:
                            countseqs+=1
                            print len(record.features)
                            SeqIO.write(shortrec, out_handle, "fasta")
                            SeqIO.write(longrec, out_ext_handle, "genbank")
                            # check if we have found the record before, however likely needs improvement.
                            if not gi_num in Accnums:  # check to see if we have pulled the record yet
                                Accnums.append(gi_num)
                                if save_complete_seqs:
                                    SeqIO.write(record, comp_out_handle, "genbank")
                            print "wrote " + gi_num
    
        except IndexError:
            break
    inhandle.close()
    if save_complete_seqs:
        comp_out_handle.close()
    out_handle.close()
    out_ext_handle.close()
    print "Number of hits: %s" % len(Accnums)
    print "Number of qualifying records: %s" % countseqs
    
    
def retrieve_from_new_blast_table(infile = "ISCRHitTable.csv", compseqout = "ISCRCompleteAlignTest.fasta", outfile = "ISCRAlignTest.fasta", extoutfile = "ISCRTestExtAlign.gb", targetlen = 513, featname = "ISCR", down_only = False):
    # Blast table is available from web based alignment, export .csv
    # opens a blast table and downloads the best alignment for each hit
    # Inputs:
    #   infile - name of file for BLAST output in csv format
    #   compseqout - name of file for complete sequences of the alignment, useful for building phyogenetic trees, etc.
    #   outfile - name of file for sequences of the alignment region in fasta format, useful for building phyogenetic trees, etc.
    #   extoutfile - name of file for outputting the alignment regions
    #   targetlength - nucleotide length of the the target sequence, used to compute max percent coverage
    #   featname - name of the alignment region, used when adding a feature to the alignment region to ensure it is annotated on the sequence
    #
    # Retrieved sequences must meet the following criteria:
    #   Minlength = length of the extension region or twice its length if extending in both directions
    #   min_per_cov = minimum allowable percent coverage (hard coded)
    #   min_per_id =  minimum allowable percient identity in alignment (hard coded)
    #   region must containg annotated features
    #
    # The following are output:
    #   genbank file with extended window around alignment
    #   fasta file of alignment regions
    #   fasta file of complete records of aligned sequences
    Entrez.email = "ryanbotts@pointloma.edu"
    regionext = 10000 # length of sequence to retrieve beyond alignment
    min_per_cov = .70
    min_per_id = .30
    print 'Only using alignments greater than %.2f id and %.2f coverage.' % (min_per_id, min_per_cov)
    if down_only:
        print "Only extending in one direction, assumes protein used in search is pointing forward."
        total_len = regionext
        len_up = regionext
        len_down = 0 #regionext
    else:
        print "extracting up and downstream region"
        total_len = 2*regionext
        len_up = regionext
        len_down = regionext

    inhandle = open(infile,"rU")
    csv_file = csv.reader(inhandle,delimiter=',',dialect=csv.excel)
    
    save_complete_seqs = True # save the complete record for the alignment seq
    out_handle = open(outfile,"w")
    if save_complete_seqs:
        comp_out_handle = open(compseqout,"w")
    out_ext_handle = open(extoutfile,"w")
    Accnums=[] # store all accession numbers
    Best_aligns={}
    countseqs = 0 # keep track of how many sequences meet criteria

    for row in csv_file:
        try:
            gi_num = row[3]
            print gi_num
            id = float(row[4])/100
            print row[5]
            align_length = int(row[5])
            print row[10]
            sbjct_start = int(row[10])
            print row[11]
            sbjct_end = int(row[11])
            e_value = float(row[12])
            
            print "Acc Number: %s, Identity: %s, Align Len: %s " % (gi_num, id, align_length)
            
            #determine orientiation of alignment, for uniform ordering and determine the appropriate start and stop coordinates
            if sbjct_start > sbjct_end:
                rev_align = True
                true_end = sbjct_start
                true_start = sbjct_end-1
                len_left = len_down
                len_right = len_up
                feature = SeqFeature(FeatureLocation(true_start, true_end), strand=-1, type ='CDS', id = featname,qualifiers={"locus_tag":featname})
            else:
                rev_align = False
                true_end = sbjct_end;
                true_start = sbjct_start-1
                len_left = len_up
                len_right = len_down
                feature = SeqFeature(FeatureLocation(true_start, true_end), strand=1, type ='CDS', id = featname,qualifiers={"locus_tag":featname})

 
            # compute end of extended region
            start = true_start - len_left
            end = true_end + len_right
            
            
            # check the alignment is of high enough quality
            if e_value < .0001 and align_length > targetlen * min_per_cov:
                # link to record
                handle=Entrez.efetch(db='nucleotide',id=gi_num,rettype='gb')
                record = SeqIO.read(handle,"genbank")
                #print record.description
                if record.description.find('plasmid')>0 or record.description.find('Plasmid')>0:
                    print "Found plasmid"
                    record.id = gi_num+'_l_%s-%s_pl' % (true_start,true_end)
                    #record.name = record.id # should be the record name, but this is too long
                    record.description = record.id + ' ' + record.description
                else:
                    record.id = gi_num+'_l_%s-%s' % (true_start,true_end)
                    #record.name = record.id
                    record.description = record.id + ' ' + record.description
                #print record.id
                #record.id = gi_num+"_l_%s-%s" % (true_start,true_end)
                
                record.features = record.features + [feature]
                
                # check that the record is long enough
                if len(record) > total_len: #and (record.description.find("complete")> 0):
                    # order all of the aligned regions so that they are in the appropriate orientation
                    if rev_align:
                        shortrec = record[true_start : true_end].reverse_complement()
                    else:
                        shortrec = record[true_start : true_end]
                    
                    print "Original record length is: %s Window length: %s" % (str(len(record)),str(end - start))
                    
                    # extract the appropriate region if the extended window wraps around the end of the sequence or extends past the end, assumes circular DNA
                    #if start < 0:
                    #    longrec= record[start:]+record[:end]
                    #if end > len(record):
                    #    longrec = record[start:] + record[:end-len(record)]
                    #else:
                    #    longrec = record[start : end]
                    longrec = record[ max(0,start): min(len(record),end)]
                    longrec.description = record.description
                    shortrec.id = record.id
                    shortrec.description = ''
                    print "The final lenth is %s" % len(longrec)
                    print "Number of features in align region: %s" % len(longrec.features)
                    
                    # check that there aren't a bunch of NNNN's inalignment region, not sure why this would happen
                    if len(longrec.features)>5 and longrec.seq.find("NNNNNN")<0 and shortrec.seq.find("NNNNNN")< 0:
                            countseqs+=1
                            print len(record.features)
                            SeqIO.write(shortrec, out_handle, "fasta")
                            SeqIO.write(longrec, out_ext_handle, "genbank")
                            # check if we have found the record before, however likely needs improvement.
                            if not gi_num in Accnums:  # check to see if we have pulled the record yet
                                Accnums.append(gi_num)
                                if save_complete_seqs:
                                    SeqIO.write(record, comp_out_handle, "fasta")
                            print "wrote " + gi_num
    
        except IndexError:
            break
    inhandle.close()
    if save_complete_seqs:
        comp_out_handle.close()
    out_handle.close()
    out_ext_handle.close()
    print "Number of hits: %s" % len(Accnums)
    print "Number of qualifying records: %s" % countseqs

def retrieve_seqs_using_blast_table(infile = "ISCRHitTable.csv", ingb = "ISCRComplete.gb", outfile = "ISCRAlign.fasta", extoutfile = "ISCRExtAlign.gb", targetlen = 513, featname = "ISCR", down_only = False):
    # Same as retrieve_from_blast_table, with additional input
    # ingb is a file already containing the complete genbank records of the alignments
    # useful for modifying the criteria for the window around the sequences
    
    regionext = 10000 # length of sequence to retrieve beyond alignment
    min_per_cov = .70
    min_per_id = .30
    print 'Only using alignments greater than %.2f id and %.2f coverage.' % (min_per_id, min_per_cov)
    if down_only:
        print "Only extending in one direction, assumes protein used in search is pointing forward."
        total_len = regionext
        len_up = regionext
        len_down = 0 #regionext
    else:
        print "extracting up and downstream region"
        total_len = 2*regionext
        len_up = regionext
        len_down = regionext

    inhandle = open(infile,"rU")
    csv_file = csv.reader(inhandle,delimiter=',',dialect=csv.excel)
    out_handle = open(outfile,"w")
    out_ext_handle = open(extoutfile,"w")

    Accnums=[] # store all accession numbers
    Best_aligns={}
    countseqs = 0 # keep track of how many sequences meet criteria

    for row in csv_file:
        try:
            gi_num = row[1].split("|")[3]
            id = float(row[3])/100
            print row[6]
            align_length = int(row[6])
            print row[11]
            sbjct_start = int(row[11])
            print row[12]
            sbjct_end = int(row[12])
            e_value = float(row[13])
            
            print "Acc Number: %s, Identity: %s, Align Len: %s " % (gi_num, id, align_length)
            
            #determine orientiation of alignment, for uniform ordering and determine the appropriate start and stop coordinates
            if sbjct_start > sbjct_end:
                rev_align = True
                true_end = sbjct_start
                true_start = sbjct_end-1
                len_left = len_down
                len_right = len_up
                feature = SeqFeature(FeatureLocation(true_start, true_end), strand=-1, type ='CDS', id = featname,qualifiers={"locus_tag":featname})
            else:
                rev_align = False
                true_end = sbjct_end;
                true_start = sbjct_start-1
                len_left = len_up
                len_right = len_down
                feature = SeqFeature(FeatureLocation(true_start, true_end), strand=1, type ='CDS', id = featname,qualifiers={"locus_tag":featname})

 
            # compute end of extended region
            start = true_start - len_left
            end = true_end + len_right
            
            
            # check the alignment is of high enough quality
            if e_value < .0001 and align_length > targetlen * min_per_cov:
                # find and load record from gb file
                print gi_num
                record = next(seq for seq in SeqIO.parse(open(ingb,'rU'),'genbank') if gi_num in seq.id )
                print record
                
                #print record.description
                if record.description.find('plasmid')>0 or record.description.find('Plasmid')>0:
                    print "Found plasmid"
                    record.id = gi_num+'_l_%s-%s_pl' % (true_start,true_end)
                    #record.name = record.id # should be the record name, but this is too long
                    record.description = record.id + ' ' + record.description
                else:
                    record.id = gi_num+'_l_%s-%s' % (true_start,true_end)
                    #record.name = record.id
                    record.description = record.id + ' ' + record.description
                #print record.id
                #record.id = gi_num+"_l_%s-%s" % (true_start,true_end)
                
                record.features = record.features + [feature]
                
                # check that the record is long enough
                if len(record) > total_len: #and (record.description.find("complete")> 0):
                    # order all of the aligned regions so that they are in the appropriate orientation
                    if rev_align:
                        shortrec = record[true_start : true_end].reverse_complement()
                    else:
                        shortrec = record[true_start : true_end]
                    
                    print "Original record length is: %s Window length: %s" % (str(len(record)),str(end - start))
                    
                    # extract the appropriate region if the extended window wraps around the end of the sequence or extends past the end, assumes circular DNA
                    #if start < 0:
                    #    longrec= record[start:]+record[:end]
                    #if end > len(record):
                    #    longrec = record[start:] + record[:end-len(record)]
                    #else:
                    #    longrec = record[start : end]
                    longrec = record[ max(0,start): min(len(record),end)]
                    longrec.description = record.description
                    shortrec.id = record.id
                    shortrec.description = ''
                    print "The final lenth is %s" % len(longrec)
                    print "Number of features in align region: %s" % len(longrec.features)
                    
                    # check that there aren't a bunch of NNNN's inalignment region, not sure why this would happen
                    if len(longrec.features)>5 and longrec.seq.find("NNNNNN")<0 and shortrec.seq.find("NNNNNN")< 0:
                            countseqs+=1
                            print len(record.features)
                            SeqIO.write(shortrec, out_handle, "fasta")
                            SeqIO.write(longrec, out_ext_handle, "genbank")
                            # check if we have found the record before, however likely needs improvement.
                            if not gi_num in Accnums:  # check to see if we have pulled the record yet
                                Accnums.append(gi_num)
                            print "wrote " + gi_num
    
        except IndexError:
            break
    inhandle.close()
    out_handle.close()
    out_ext_handle.close()
    print "Number of hits: %s" % len(Accnums)
    print "Number of qualifying records: %s" % countseqs


def Find_Sequence(target,outname="ISCR1hits.gb",outextendedname = "ISCR1alignments.fasta",E_VALUE_THRESH=0.0001, MINCOVER = 0.70):
	# open a fasta target sequence, BLAST against Genbank
    # download up and downstream regions of hit and save in output file
    # designed to automate the BLAST search for downloading alignment regions, in practice this was too slow and abandonded
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
            #print ("query cover: %s    max_iden: %s" % (cover, max_iden))
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
                    print ("query cover: %s    max_iden: %s" % (cover, max_iden))
                    handle=Entrez.efetch(db='nucleotide',id=gi_num,rettype='gb') # Accession id works, returns genbank format, looks in the 'nucleotide' database
                
                    # read sequence
                    record = SeqIO.read(handle,"genbank")
                    aligns.append(record[sbjct_start:sbjct_end])
                    hitExtseqs.append(record[sbjct_start-5000:sbjct_end+5000])
                    print record.id

    print "%s Alignments have been found" % count
    SeqIO.write(hitExtseqs, outhandle,"genbank")
    SeqIO.write(aligns,outextendedhandle,"fasta")
    save_file.close()
    outhandle.close()
    outextendedhandle.close()
    print("----------------------- ORF BLAST COMPLETE ------- %s records found----------" % (str(count)))


###  Tools for extracting the proteins on a list of records and naming them in a convenient manner

def StripCDS(infile="ISCR1hits.gb", outfile = "ISCR1FeaturesAA.fasta", nucoutfile = "ISCR1ReaturesNN.fasta"):
    # Utility for extracting all CDS features from a file of genbank records.
    # Exports both nucleotide and amino acid sequences of all records
    #
    # Naming convention follows "gene name"-"numerical position on genome"_:_"CDS"_:_"Accession number or first word in description"
    # Numerical position on genome was added due to multiple genes of the same name on a genome.
    #
    in_handle = open(infile,"rU")
    inseqs = SeqIO.parse(in_handle,"genbank")
    outhandle = open(outfile,"w")
    nucouthandle = open(nucoutfile,'w')
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
                        #print feature.qualifiers
                        # name the feature based on the gene or product name.  If neither is available, set to orf something
                        if 'gene' in feature.qualifiers:
                            feat_name=feature.qualifiers['gene'][0].replace(" ","")
                        elif 'product' in feature.qualifiers:
                            feat_name=feature.qualifiers['product'][0].replace(" ","")
                        else:
                            feat_name=str(['orfX'+str(j)]).replace(" ","")
                        print feat_name + " found on "+ record.id
                        sq=feature.extract(record)
                        
                        #org = record.description.split(' ')[0]+'_'+record.description.split(' ')[1] # try to find the organism and add it to the name
                        #sq.id = feat_name+"-"+str(j)+"_:_"+"CDS"+"_:_"+record.id#+'_:_'+org #record.organism.replace(' ','_')
                        # Using description to store name because the naems are too long with the positions in them
                        sq.id = feat_name+"-"+str(j)+"_:_"+"CDS"+"_:_"+record.description.split(' ')[0]#+'_:_'
                        # only right features that are nonempty (not sure why some are)
                        if len(sq)>0:
                            SeqIO.write(sq,nucouthandle,"fasta")
                            sq.seq = sq.seq.translate(table="Bacterial", to_stop=True)
                            SeqIO.write(sq,outhandle,"fasta")
                
            
                except KeyError:
                    print record.id + ' has no gene features'
    outhandle.close()
    nucouthandle.close()
    print "%s sequences in analysis" % k
    print "%s features for analysis" % m
    print "m is ", m

### General utilities for handling lists of sequences, retrieving records from a list, extracting plasmid names, downloading records from a list, etc.
def get_Target_Names(infile="ISCR1extAlign.gb"):
    # get all of the record id's for sequences in a multi-genbank file
    # built to allow a user to manually remove records from the analysis
    in_handle = open(infile,"rU")
    inseqs = SeqIO.parse(in_handle,"genbank")
    print "Using name stored in description."
    # Full name stored in the record because it was too long to store as the id or name.
    return [record.description.split(' ')[0] for record in inseqs]


def extract_plnames_from_file(infile = "Plasmids.gb"):
    # Given a file of genbank records, try to strip the plasmid names
    # outputs into a csv that has the name, accession number and can be modified to have the Inc group
    # Automates finding the plasmid name in the genbank record
    # In places where the name cannot be called well, the csv output can be modified to manually call the correct name
    # Output:
    #   csv file with the following columns:  AccNum, Plasmid name, genus, species
    
    outfname = infile.strip(".gb")+".csv"
    print "Writing names to %s." % outfname
    csv_file = csv.writer(open(outfname,"w"),delimiter=',',dialect=csv.excel)

    count = 0  # a counter for the number of plasmids in the file
    total = 0  # count the total number of sequences in the file
    print 'Extracting names from file, requires that the name of the plasmid follows the word plasmid in the genbank description'

    with open(infile,"rU") as inhandle:
        for sq in SeqIO.parse(inhandle,"genbank"):
            print sq.description
            total+=1
            Genus=sq.description.split(" ")[0]
            Species=sq.description.split(" ")[1]
            AccNum=sq.name
            print AccNum
            if sq.description.find('plasmid')>0:
                # logic takes the name after the word plasmid if it exists, otherwise takes the name before it.
                #prior = sq.description.split("plasmid")[0].split(" ")
                
                name = sq.description.split("plasmid")[1].split(" ")[1].strip(",").strip(".")
                print name

                #if len(sq.description.split('plasmid ')[1].split(' ')) > 0:
                #    name = sq.description.split('plasmid ')[1].split(' ')[0] # take the word after plasmid as the name
                #elif len(sq.description.split('plasmid ')[1].split(' ')) == 0:
                #    name = sq.description.split('plasmid ')[0].split(' ')[0] # take the word after plasmid as the name
                count+=1
                print count
            else:
                continue
            #if count > 30:
            #    break
            row = [AccNum,name,Genus,Species]
            csv_file.writerow(row)
            total+=1
    print "Found %i seqs in expected list of %i." % (count, total)

def StripCDSAddInc(ingb,incsv):
    # routine opens a file with plasmid names and inc groups, if known, along with gb file with seqs
    # strips all features from each genbank record found in the csv file and saves the amino acid sequence, named with inc group
    # open all genbank records in ingb and reads all names and groups in the csv
    # naming convention of genes is:
    #   "gene_name"-"number"_:_"IncGroup"_:_"PlasmidName"
    # Inc Group is set as IncUnk if it is not known
    # Inputs:
    #   ingb-input file of all of the genbank records for analysis
    #   incsv-input file with the plasmid name, inc group, and accession number, used to name the features appropriately
    #
    # Outputs:
    #   outfa-fasta output of amino acid sequences of all of the genes, named accordingly, file is named the same as the input
    #           ingb with the appendix changed to AA.fa.
    inseqs = SeqIO.parse(open(ingb,"rU"),"genbank")
    
    seqdict = {x[0]:x for x in GetSeqsAndGroupsFromCSV(incsv)}
    
    outfa = ingb.strip(".gb")+"AA.fa"
    outhandle = open(outfa,"w")

    print "Stripping all CDS's from seqs"
    
    k=0 #count the number of records
    m=0 # count the number of features
    for record in inseqs:
        if record.name in seqdict:
            print seqdict[record.name]
            plname=seqdict[record.name][1] # Save the plasmid name
            group = seqdict[record.name][4] # Save the group name
            k+=1
            j = 0 # count the number of orfs on each plasmid
            for feature in record.features:
                try:
                    if feature.type == "CDS":
                        m+=1
                        j+=1
                        #print feature.qualifiers
                        # name the feature based on the gene or product name.  If neither is available, set to orf something
                        if 'gene' in feature.qualifiers:
                            feat_name=feature.qualifiers['gene'][0].replace(" ","")
                        elif 'product' in feature.qualifiers:
                            feat_name=feature.qualifiers['product'][0].replace(" ","")
                        else:
                            feat_name=str(['orfX'+str(j)]).replace(" ","")
                        print feat_name + " found on "+ record.id
                        sq=feature.extract(record)
                        # Using description to store name because the naems are too long with the positions in them
                        sq.id = feat_name+"-"+str(j)+"_:_"+group+"_:_"+plname
                        # only right features that are nonempty (not sure why some are)
                        if len(sq)>0:
                            sq.seq = sq.seq.translate(table="Bacterial", to_stop=True)
                            SeqIO.write(sq,outhandle,"fasta")
                
            
                except KeyError:
                    print record.id + ' has no gene features'
    outhandle.close()

    print "%s sequences in analysis" % k
    print "%s features for analysis" % m


def extract_seqs_from_file(infile = "ISCR1Keeps.csv",targetfile = "ISCR1Align.gb", outseqfile = "ISCR1AlignSubset"):
    # Utility for selecting a subset of sequences from a file with multiple records
    # Given a file of seqs and a csv file with names, extract all sequences in the file with those names
    # csv file must have the seqs listed in the first column
    # Inputs:
    #   infile-csv file with all of the names of the desired sequences, searches by accession number
    #   targetfile-fasta or genbank file with all records
    #   outseqfile-name of the file for outputting only the extracted sequences
    # Outputs:
    #   outseqfile-all of the records retrieved from the larger set of records
    inhandle = open(infile,"rU")
    csv_file = csv.reader(inhandle,delimiter=',',dialect=csv.excel)
    seqnames =  [row[0] for row in csv_file]
    inhandle.close()

    seqhandle = open(targetfile,"rU")
    
    # identify the type of reads we are using
    if targetfile.find("gb") > 0:
        type = "genbank"
        outname = outseqfile + '.gb'
    elif targetfile.find("fasta") > 0:
        type = "fasta"
        outname =  outseqfile + '.fasta'

    count = 0
    print 'Extracting features from file.  Sequence name must be a subset of names in'
    outseqhandle = open(outname, 'w')
    
    for sqname in seqnames:
        print "Searching for %s" % sqname
        for sq in SeqIO.parse(open(targetfile,"rU"),type):
        
            if type == "genbank":
                sqnm = sq.description.split(' ')[0]
            
            else:
                sqnm = sq.id
            
            if sqnm in sqname:
                print "record id %s" % sq.id
    
    #for sq in SeqIO.parse(seqhandle,type):
    #    print sq.id
    #    if any(sq.id in sqname for sqname in seqnames):
                SeqIO.write(sq, outseqhandle, type)
                count+=1
    outseqhandle.close()

    print "Found %i seqs in expected list of %i." % (count, len(seqnames))

def GetSeqsAndGroupsFromCSV(infile):
    ### Reads csv file with Acc number, name and Inc Group
    # Assumes table format AccNum, Name, Genus, Species, Group(optional)
    # Keeps track of Inc Group if known, otherwise assigned IncUKN
    # Returns list, where each entry has this info
    print "Reading sequences and groups from CSV file, check that file contains updated Inc groups"
    inhandle = open(infile,"rU")
    csv_file = csv.reader(inhandle,delimiter=',',dialect=csv.excel)
    
    #Reads csv file and adds unknown group if that column is unknown
    out = [row if row[4]!= "" else row[0:4] +["IncUNK"] for row in csv_file]
    return out


# extracting gc content needs to be modified
def compare_gc_content(infile = "ISCR1Align.csv", outfile = "ISCR1GC.csv", outseqfile = "ISCR1AlignsComplete.gb",targetlen = 513):
    # Blast table is available from web based alignment, export .csv
    # opens a blast table and retrieves the GC content of the entire record
    # and of a smaller insertion sequence
    # also retrieves the complete genbank record of the ISCR element
    Entrez.email = "ryanbotts@pointloma.edu"
    regionext = 10000 # amount of sequence to gather in both directions from target
    inhandle = open(infile,"rU")
    out_handle = open(outseqfile,"w")
    csv_file = csv.reader(inhandle,delimiter=',',dialect=csv.excel)
    print "opening input file"
    Accnums=[] # store all accession numbers
    Best_aligns={}
    GCvals = {} # key will be accession = [len of record, GC overall, GC element]
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
            if sim >.50 and align_length > targetlen*.50:
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
                    
                    totGC= GC(record.seq)
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
                    localGC = GC(longrec.seq) # store the GC content of the IS sequence
                    print "Number of features: %s" % len(longrec.features)
    
                    # save GC contents
                    GCvals[gi_num]=[len(record),totGC,localGC]
                    print "%s Overall GC content: %.2f, Insertion content: %.2f" % (gi_num, totGC, localGC)
                    
                    SeqIO.write(record,out_handle,"genbank")
    
        except IndexError:
            break
    inhandle.close()
    out_handle.close()
    writer = csv.writer(open(outfile, 'wb'))
    writer.writerow(['AccNum','Record Length','Overall GC', 'Region GC'])
    for key, value in GCvals.items():
        writer.writerow([key, value[0],value[1]])
    print GCvals
    print "Number of sequences: %s" % len(Accnums)


### Tools for handling the output of USEARCH clustering
def get_clusters(clusterfilename):
    ### finds clusters from default tab delimited output from USEARCH
    # creates multiple lists each containing all items in a cluster
    # input:
	#	clusterfilename - tab delimited default output from USEARCH
	#		center of cluster begins with a S in the first column and name in the 9th column
    # output:
	# 	groups - list containing all elements of each cluster
    
    groups = []  # empty dictionary to store each cluster
    
    #with open(clusterfilename, 'rU') as f:
     #   reader = csv.reader(f, delimiter='\t')
      #  found_clusters = []
       # for row in reader:
        #    if row[0].strip() == "C":
                #do nothing
        #        trash = row
        #    else:
        #        num = int(row[1].strip())
        #        print num
        #        if num not in found_clusters:
        #            found_clusters.append(num)
        #print found_clusters
        #Num_of_Clusters = max(n for n in found_clusters)
        #Total_groups = []
        #group = []
        #for x in range(0, Num_of_Clusters+1):
        #    f.seek(0)
        #    for row in reader:
        #        if row[0].strip() == "C":
        #            #do nothing
        #            trash = row
        #        else:
        #            num2 = int(row[1].strip())
        #            if num2 == x:
        #                print num2, " equals ", int(row[1].strip())
        #                name = row[8].strip()
        #                group.append(name)
        #    print group
        #    Total_groups.append(group)
        #    print "group added"
        #    #print Total_groups
        #    del group[:]
            
        #print Total_groups
            
    #return Total_groups
    
    with open(clusterfilename, 'rU') as f:
        reader = csv.reader(f, delimiter='\t')
        first = True # identifies if we are on the first group or not
        group=[]
        for row in reader:
            name = row[8].strip()
            #Fullname = row[8].strip()
            #name = Fullname.split("_:_",2)[0]
            #print name
            if row[0].strip() == "S" and first:
                group=[name]
                first = False
                groups.append(group)
            elif row[0].strip() == "S" and not first:
            	# add previous group to dictionary
                # start a new group
                group=[name]
                groups.append(group)
            elif row[0].strip() == "H":
                num = int(row[1].strip())
                temp = groups[num]
                temp.append(name)                
                groups[num] = temp
            # the rows marked with a C indicate the center of a cluster and may be omitted
            elif row[0].strip() == "C":
                #groups.append(group)
                break
    return groups

def get_clust_name(pro_list):
    ### Attempts to determine an appropriate name for all proteins in a cluster
    # logic needs improvement,
    # names the group with the following priorities
    #   1.  a protein in the list has one of the preferred names, then that is used first
    #   2.  the first 4 letter protein name in the list
    #   3.  the first 6 (for 4 letter names with a number) letter named protein in the list
    #   4.  the first protein name in the list that is not hypothetical or orf
   
    #print 'Creating protein family names, attempting to name the cluster'
    #print 'Uses the following names if possible:  sul1, qacEdelta, blaCTX-M,aac, aad4'
    
    # keep track of how long the name is and whether it is the best we have or not.
    min_len_name = 0
    found_best = False
    
    pref_names = ['sul1','qacEdelta','blaCTX-M','aac','aad4','ISCR', 'tnpR','tnpA','ins1','mer', 'TEM']
    
    # default name is the first name in the group of nothing better is found
    group_name = pro_list[0]

    # take the next item satsfying the following criteria, if one is not found, procede to the next statement
    # if none of these criteria works, use default name

    group_name = next((x for x in pro_list for y in pref_names if y in x.split('_:_')[0]),\
        next((x for x in pro_list if len(x.split('_:_')[0].split('-')[0]) == 4),\
        next((x for x in pro_list if len(x.split('_:_')[0].split('-')[0]) == 6),\
            group_name)
        ))

    return group_name

def remove_clusters(groups):
    # Attempt to remove clusters that aren't interesting
    # current criteria for removal
    # >50% of the cluster is named using one of the unwanted names
    # removes clusters with fewer than 5 in the cluster
    # remove clusters that are transposases, hypothetical proteins, or a collection of antibiotic resistance proteins
    avoid_names = ['ins','int','tnp','hypotheticalprotein','transposase','integrase','sul1','qacedelta','bla','aac','aad4','iscr', 'tnpr','tnpa','ins1','mer', 'tem','ndm']
    return  [ group for group in groups if sum([any([a in g.split("_:_")[0].lower() for a in avoid_names]) for g in group ])/len(group)<.5 and len(group)>5]

def save_qualifying_clusters(fname,clusters):
    # removes clusters that don't satisfy the criteria from clusters using the function remove_clusters
    # writes output to a csv.
    # each row begins with the number of proteins in the cluster followed by all of the items in the cluster
    # inputs:
    #   fname-name of output csv file
    #   clusters -  list of lists of all members of each gene cluster
    # ouputs:
    #   csv file with all remaining clusters
    cs = csv.writer(open(fname,"w"),delimiter=',',dialect=csv.excel)
    
    print "Initially %i clusters" % len(clusters)
    clusters = remove_clusters(clusters)
    print "Remaining %i clusters" % len(clusters)
    for group in clusters:
        cs.writerow([len(group)]+ group)


def write_matrix_file(filename='matrix_noDupNames.csv',targetnamefile = 'ISCR1extAlign.gb',clusterfile = 'ISCR1Clusters.tab',writefullnames=False):
    # This function creates a matrix that indicates the presence or absence of a particular family of proteins on each sequence in a list
    # Idea:
    #   1. List of records of interest, complete records or partial records
    #   2. Use Usearch clustering of proteins into families, so that we know which proteins are the same regardless of name
    #   3. (this function) Using the distinct clusters of proteins found by USEARCH, create a matrix indicating the presence or absence of each protein family on every record in the list
    #   4. Analyze these features using heirarchical clustering or summary statistics.
    # inputs:
    #   filename:
    # output: A csv file that contains a matrix of 1s and 0s with the plasmid names down the side and the cluster names at the top.
    row=""
    notEmp=0
    names = get_Target_Names(targetnamefile)
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
                f.write(get_clust_name(groups[c]).split('_:_')[0].split('-')[0]+',')
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
                #    f.write(row)
            
            notEmp=0
            row=""
        f.close()

def gather_row_seqs(inseqs = "ISCR1Align.fasta", plasmidmatrix = "ISCR1FeatureMatrixShortNames.csv", outseqs = "ISCR1AlignClean.fasta"):
    # make sure that the sequence file only contains sequences which are included in the final feature matrix.
    # tool is used in case feature matrix has been edited by hand
    seq_handle =  open(inseqs,"r")
    
    seqs = SeqIO.parse(seq_handle, "fasta")
    seqlist = {s.id : s for s in seqs}
    
    seq_handle.close()
    
    
    in_handle =  open(plasmidmatrix,"rU")
    table = csv.reader(in_handle,delimiter=',')
    table.next()
    
    seqnames =[row[0] for row in table]
    
    out = [seqlist[k] for s in seqnames for k in seqlist.keys() if k.find(s)>-1]
    print len(out)
    out_seq_handle = open(outseqs,"w")
    SeqIO.write(out,out_seq_handle, "fasta")
    out_seq_handle.close()

def map_seq_groups(outfilename='ISCR1MapGroup', plasmidftfile = 'ISCR1featurematrix.csv', plasmidclustfile = 'ISCR1Groups.csv', plasmidseqs='ISCR1extAlign.gb',\
    seqclusters = 'ISCR1Clusters.tab', pathtoplasmids='./', featname = 'ISCR1'):
	# modification of previous version, accounts for having all seqs in one file
    # creates map of all plasmids within a cluster and colors the genes
    # according to the percent of the time the gene occurs in the all plasmids within cluster.
    # inputs:
    #    plasftfile - feature matrix name for each group  All such files will be opened and read in
    #    plasmidgroup - csv file with plasmid name in col1 and group number in col2
    #    plasmidseqs - file with all plasmid genbank files
    #    plasmidmatrix - matrix of plasmids with consensus groups
    #    seqclusters - csv file containing each of the clusters of genes
    #    outfilename - name of file, -i.csv will be appended to the name
    #    pathtoplasmids - path extension to the file with all of the plasmid sequences
    # output:
    #    map for each group of plasmids, written in the folder containing the plasmids
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
    # name had to be stored in definition instead of in id
    seqs = {s.description.split(' ')[0] : s for s in SeqIO.parse(handle,"genbank")}

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
        
        print len(genepercents)

        nonzeroGroups = {names[i]: genepercents[i] for i in range(len(genepercents)) if genepercents[i]>0}
        print nonzeroGroups
        
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
                            #grpfeatplasmid[gene.split('_:_')[0]] = nonzeroGroups[feat]
                            grpfeatplasmid[gene.split('_:_')[0].split('-')[0]+'_:_'+feat.split('-')[0]] = nonzeroGroups[feat]
                except IndexError:
                    print feat + " not found in seqs, must have been removed from the cluster table"
            record = seqs[plasmid]
            
            
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
            print record.features[-1].qualifiers['locus_tag'][0]
            for feature in record.features:
                # identify the alignment feature
                if 'locus_tag' in feature.qualifiers:
                    if feature.qualifiers['locus_tag'][0] == featname:
                        print 'Found %s region on sequence.' % featname
                        gd_feature_set.add_feature(feature, sigil="BIGARROW",
                                   arrowshaft_height=.8,
                                   arrowhead_length=.25,
                                   color=colors.green, border=colors.green,label=True,
                                   name = featname,
                                   label_position="center",
                                   label_size = 7, label_angle=45)
                
                elif feature.type=="CDS":
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
                                x = 1
                                #print feature
                                #print '%s not found on %s, skipping' % (feature.id, plasmid)
                        except KeyError:
                            print 'CDS does not have gene name on '+plasmid



        #width = max(1500,max_length/100)
        height= len(groups[grp])*100
        gd_diagram.draw(format="linear", pagesize=(2000,height), fragments=1,
                 track_size =.2,start=0)
        gd_diagram.write(pathtoplasmids + outfilename +'-'+ grp+ ".pdf", "PDF")

def map_seq_groups_old(outfilename='ISCR1MapGroup', plasmidftfile = 'ISCR1featurematrix.csv', plasmidclustfile = 'ISCR1Groups.csv', plasmidseqs='ISCR1extAlign.gb',\
    seqclusters = 'ISCR1Clusters.tab', alignfile = 'ISCR1Align.csv', pathtoplasmids='./', featname = 'ISCR'):
	# modification of previous version, accounts for having all seqs in one file
    # creates map of all plasmids within a cluster and colors the genes
    # according to the percent of the time the gene occurs in the all plasmids within cluster.
    # inputs:
    #    plasftfile - feature matrix name for each group  All such files will be opened and read in
    #    plasmidgroup - csv file with plasmid name in col1 and group number in col2
    #    plasmidseqs - file with all plasmid genbank files
    #    plasmidmatrix - matrix of plasmids with consensus groups
    #    seqclusters - csv file containing each of the clusters of genes
    #    outfilename - name of file, -i.csv will be appended to the name
    #    pathtoplasmids - path extension to the file with all of the plasmid sequences
    # output:
    #    map for each group of plasmids, written in the folder containing the plasmids
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
            feature.id = featname
            alignloc[gi_num] = feature
        
        else:
            true_end = sbjct_start
            true_start = sbjct_end
            feature = SeqFeature(FeatureLocation(true_start, true_end), strand=-1)
            feature.id = featname
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
                            #grpfeatplasmid[gene.split('_:_')[0]] = nonzeroGroups[feat]
                            grpfeatplasmid[gene.split('_:_')[0].split('-')[0]+'_:_'+feat.split('-')[0]] = nonzeroGroups[feat]
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
                                   name = featname,
                                   label_position="center",
                                   label_size = 7, label_angle=45)

        #width = max(1500,max_length/100)
        height= len(groups[grp])*100
        gd_diagram.draw(format="linear", pagesize=(2000,height), fragments=1,
                 track_size =.2,start=0)
        gd_diagram.write(pathtoplasmids + outfilename +'-'+ grp+ ".pdf", "PDF")

def convert_feature_matrix_for_mba(matrixfilename, outfilename):
    # reads in a file of binary presence absence data and returns rows with list of features in each basket
    handle = open(matrixfilename,'rU')
    tablegroups = csv.reader(handle, delimiter=',',dialect=csv.excel_tab)
    outhandle = open(outfilename, 'w')
    basketwriter = csv.writer(outhandle, delimiter=',',dialect=csv.excel_tab)
    
    print "Obtaining item names for basket, must be unique names."
    itemnames = tablegroups.next() # obtain row with headers
    #print itemnames
    #print itemnames
    for row in tablegroups:
        #print row
        basket = [itemnames[i] for i in range(1,len(row)) if row[i] == '1.0' or row[i] == '1']
        basketwriter.writerow(basket)
        #print basket
    

    handle.close()
    outhandle.close()

def extract_short_seqs(seqfile, outfilename, minlen):
    # extracts all of the records in a file that are greater than some min length
    outhandle = open(outfilename,'w')
    inhandle = open(seqfile,'rU')
    if seqfile.find('.gb')>0:
        filetype = 'genbank'
    else:
        filetype = 'fasta'
    for sq in SeqIO.parse(inhandle, filetype):
        if len(sq) >  minlen-1:
            SeqIO.write(sq, outhandle, filetype)
    outhandle.close()
    inhandle.close()

def get_seq_properties(seqfile):
    #  Input a gb file and get relevant features in one place
    #  Outputs a csv file with the relevant features from the genbank metadata
    #       outfilename = seqfile + .csv

    writer = csv.writer(open(seqfile.strip('.gb')+'.csv', 'wb'))
    writer.writerow(['AccNum','PLorCHR','Species','GC','Description', 'Source'])
    for sq in SeqIO.parse(open(seqfile,'rU'), 'genbank'):
        if sq.description.find('plasmid')>0 or sq.description.find('Plasmid')>0:
            tmp_source = 'Plasmid'
        else:
            tmp_source = 'Chromosome'
        writer.writerow([sq.id,tmp_source,sq.description.split(' ')[1]+' '+sq.description.split(' ')[2], \
            GC(sq.seq),sq.description.strip(sq.description.split(' ')[0]), sq.annotations['source']])
        print 'Writing metadata for %s.' % sq.id
        

#Top Genes
#   Input: The matrix csv file, the alignment fasta file, the name of an output csv file, and the minimum number of genes to return
#   Output: A fasta a file of the top n genes that can be used to create a phologenetic tree
def TopGenes(infile, fastaFile, outCSVfile, num_of_genes):
    #open the input file (matrix file)
    with open(infile, 'r') as inFile, open(outCSVfile, 'wb') as outFile:
        in_handle = csv.reader(inFile, delimiter = ',')
        out_handle = csv.writer(outFile, delimiter = ',')
        #declare all variables 
        #List of genes in the file and the total number of 1s in every column (each column represents a gene)
        GeneNames = []
        TotalOnes = []
        GenesReturned = []
        Sequences = []
        FullGenes = []        
        
        GeneNames = in_handle.next()
        #print len(GeneNames)
        
        for x in range(1, len(GeneNames)):
            total = 0
            inFile.seek(0)
            in_handle.next()
            for row in in_handle:
                total = total + int(row[x])
            TotalOnes.append(total)
        
        MaxList = []
        templist = []
        while (len(GenesReturned) <= num_of_genes):
            MaxNum = max(n for n in TotalOnes if n not in MaxList)
            print MaxNum
            MaxList.append(MaxNum)
            for i, j in enumerate(TotalOnes):
                if j == MaxNum:
                    templist.append(i)
            for item in templist:
                #print item
                gene = GeneNames[item+1]
                #print gene
                inFile.seek(0)
                for row in in_handle:
                    if row[item+1] == "1":
                        FullGenes.append(row[0])                
                geneName = gene.split("_:_")[0]
                name = gene.split("_:_")[2]
                NewName = geneName + "_:_" + name
                if "tnp" not in geneName:
                    GenesReturned.append(NewName)
                    input_row = [geneName, MaxNum]
                    out_handle.writerow(input_row)
                    for record in SeqIO.parse(fastaFile, 'fasta'):
                        if record.id in FullGenes:
                            #print record.id
                            #RecordName = record.id.split("_l_")[0]
                            #record.id = RecordName
                            Sequences.append(record)
                            print record.id, " added"
                
                    SeqIO.write(Sequences, geneName+"Fasta.fasta", 'fasta')
                    del Sequences[:]
                    print "record written"
                    #print FullGenes
                    del FullGenes[:]
                
            del templist[:]
                
        
        print GenesReturned
        
        
        
        
        
def TopGenesAll(infile, fastaFile, outCSVfile, min_num):
    #open the input file (matrix file) and the output csv file
    with open(infile, 'r') as inFile, open(outCSVfile, 'wb') as outFile:
        in_handle = csv.reader(inFile, delimiter = ',')
        out_handle = csv.writer(outFile, delimiter = ',')
        #declare all variables 
        #List of genes in the file and the total number of 1s in every column in parallel lists (each column represents a gene)
        GeneNames = []
        TotalOnes = []
        
        GenesReturned = []
        Sequences = []
        FullGenes = []        
        
        GeneNames = in_handle.next()
        
        #Loops through the csv file to find the told number of 1s in each column        
        for x in range(1, len(GeneNames)):
            total = 0
            inFile.seek(0)
            in_handle.next()
            for row in in_handle:
                if(row[x] != ''):
                    if row[x] == "1.0" or row[x] == '1':
                        num = '1'
                    elif row[x] == "0.0" or row[x] == '0':
                        num = '0'
                    else:
                        num = row[x]
                    adding = int(num)
                    total = total + adding
            TotalOnes.append(total)
        
        MaxList = []
        templist = []
        MaxNum = max(n for n in TotalOnes if n not in MaxList)
        print MaxNum
        MaxList.append(MaxNum)
        while(MaxNum > min_num):
            for i, j in enumerate(TotalOnes):
                if j == MaxNum:
                    templist.append(i)
            for item in templist:
                #print item
                gene = GeneNames[item+1]
                #print gene
                inFile.seek(0)
                for row in in_handle:
                    if row[item+1] == "1":
                        FullGenes.append(row[0])
                        #print row[0]
                geneName = gene.split("_:_")[0]
                name = gene.split("_:_")[2]
                NewName = geneName + "_:_" + name
                if ("tnp" not in geneName and "transposase" not in geneName) :
                    GenesReturned.append(NewName)
                    input_row = [geneName, MaxNum]
                    out_handle.writerow(input_row)
                    for record in SeqIO.parse(fastaFile, 'fasta'):
                        if record.id in FullGenes:
                            #print record.id
                            #RecordName = record.id.split("_l_")[0]
                            #record.id = RecordName
                            Sequences.append(record)
                            #print record.id, " added"
                
                    SeqIO.write(Sequences, geneName+ "Fasta.fasta", 'fasta')
                    del Sequences[:]
                    #print geneName, " ", name, " ", NewName, " ", "record written"
                    #print FullGenes
                    del FullGenes[:]
                
            del templist[:]
            MaxNum = max(n for n in TotalOnes if n not in MaxList)
            print MaxNum
            MaxList.append(MaxNum)           
        
        print GenesReturned        
     
     
#Total Genes
#   Input: Matrix name csv file
#   Output: The total number of genes found in the file is outputed to the screen 
def total_genes(infile):
    with open(infile, 'r') as inFile:
        in_handle = csv.reader(inFile, delimiter = ',')
        count = -1
        for row in in_handle:
            for field in row:
                count += 1
            return count
                
  

##Gene Summary: 
#   Inputs: the Matrix name csv.file, name of the output csv. file to be created
#   Output: A csv file is created with all of the names of the genes found and how many times
#           it was found 
def GeneSummary(infile, outfile):
    with open(infile, 'r') as inFile, open(outfile, 'wb') as outFile:
        in_handle = csv.reader(inFile, delimiter = ',')
        out_handle = csv.writer(outFile, delimiter = ',')
        total = -1
        
        for row in in_handle:
            for field in row:
                total += 1
            break
        inFile.seek(0)
        GeneNames = in_handle.next()
        Totals = []
        Totals.append(0)
        print total
        
        for x in range(1,total):
            ones = 0
            for row in in_handle:
                ones = ones + int(row[x])
                #print ones
            Totals.append(ones)
            inFile.seek(0)
            in_handle.next()
            
            
        for y in range(1, len(GeneNames)-1):
            gene = GeneNames[y]
            geneName = gene.split("_:_")[0]
            input_row = [geneName, Totals[y]]
            out_handle.writerow(input_row)

            
##GetTrees
def GetTrees():
    for fasta in glob.glob("Genes\*.fasta"):
        geneName = fasta.split('.')[0]
        geneName = geneName.split('\\')[1]
        os.system("megacc -a muscle_align_nucleotide_Non_Coding.mao -d "+fasta+" -o MegaFiles\\"+geneName+".meg")
        print "mega file done"
        os.system("megacc -a Neighbor_Joining_nucleotide_Non_Coding.mao -d MegaFiles\\"+geneName+".meg -o Trees\\"+geneName+".nwk")
        print "Tree done"              
         
##Delete empty columns in matrix
##https://gist.github.com/aculich/fb2769414850d20911eb - panda code
def DeleteEmptyColums(input_file):
    data = read_csv(input_file)
    filtered_data = data.dropna(axis='columns', how='all')
    filtered_data.to_csv(input_file, mode = 'w', index=False)
    
    
def ClusterAnalysis(clusterfilename, outputfile):

    groups = []  # empty dictionary to store each cluster

    
    with open(clusterfilename, 'rU') as f, open(outputfile, 'wb') as outFile:
        reader = csv.reader(f, delimiter='\t')
        out_handle = csv.writer(outFile, delimiter = ',')
        first = True # identifies if we are on the first group or not
        group=[]
        for row in reader:
            name = row[8].strip()
            #Fullname = row[8].strip()
            #name = Fullname.split("_:_",2)[0]
            #print name
            if row[0].strip() == "S" and first:
                group=[name]
                first = False
                groups.append(group)
            elif row[0].strip() == "S" and not first:
            	# add previous group to dictionary
                # start a new group
                group=[name]
                groups.append(group)
            elif row[0].strip() == "H":
                num = int(row[1].strip())
                temp = groups[num]
                temp.append(name)                
                groups[num] = temp
            # the rows marked with a C indicate the center of a cluster and may be omitted
            elif row[0].strip() == "C":
                #groups.append(group)
                break
        for x in range(0, len(groups)):
            temp = groups[x]
            InputRow = [temp[0], len(temp)]
            out_handle.writerow(InputRow)
                
        

    
    

##########################################################


###########################################################

if __name__=="__main__":

    ### tools for downloading ISCR sequences
    #Find_Sequence("M1-ISCR1.fasta")
    proj_name = 'ISCR'
    # 1.  Done through BLAST, download all sequences in BLAST output
    retrieve_from_new_blast_table(infile = proj_name+"HitTable.csv", compseqout = proj_name+"CompleteAlign.fasta", outfile = proj_name+"Align.fasta", extoutfile = proj_name+"ExtAlign.gb", targetlen = 513, featname = "ISCR", down_only = False)
    print "Step 1 retrieve from blast table done"
    
    # note that if usearch fails for empty file it is likely due to having written a duplicate sequence in the AA file, remove sequence
    # 2. extract all proteins from all records
    StripCDS(infile=proj_name+'ExtAlign.gb', outfile = proj_name+'FeaturesAA.fasta',nucoutfile = proj_name+'FeaturesNN.fasta')
    print "Step 2 StripCDS done"
    
    # 2b. Cluster all proteins into protein families 
    # Need to download Microsoft Visual C++ Library
    os.system('usearch8.1.1861_win32.exe -cluster_fast '+proj_name+'FeaturesAA.fasta -id 0.5 -target_cov 0.5 -centroids clust.fasta -uc '+proj_name+'Clusters.tab -userfields query+target+id+ql+tl+alnlen+qcov+tcov')
    print "Step 2b Usearch  done"
    
    # 3.  Identify protein families on each sequence and export into matrix
    write_matrix_file(proj_name+'FeatureMatrixShortNames.csv',targetnamefile = proj_name+'ExtAlign.gb',clusterfile = proj_name+'Clusters.tab',writefullnames=True)
    DeleteEmptyColums(proj_name+'FeatureMatrixShortNames.csv')
    print "Step 3 create matrix done"

    # 4. Identify the top genes inside the families 
    os.mkdir("Genes\\")  
    TopGenesAll(proj_name+"FeatureMatrixShortNames.csv","ISCRAlign.fasta", proj_name+"TopGenes.csv", 10)
    print "Step 4 finding the top genes done"

    #5. Use mega to find alignments
    GetTrees()
    
    
    #genes = total_genes("ISCRFeatureMatrixShortNames.csv")
    #print genes    
   
  