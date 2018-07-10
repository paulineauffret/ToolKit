#!/usr/bin/env python
#-*- coding: utf-8 -*-
import re
import sys
import codecs

###########################################################################################################
#	A little python script to get nucleotide sequence from a list of platax transcript IDs
#	Tested with Python 3.4
#	pauline.auffret@ifremer.fr
#	last update July, 5th 2017
###########################################################################################################
#Input and output files (required)
fileRef=sys.argv[1]		#Platax reference transcriptome (fasta file)
fileInList=sys.argv[2]		#List of transcript IDS to get the sequence (txt file, one ID per line)
fileOut=sys.argv[3]		#List of selected transcript sequences (fasta file)


#Open files
fref = codecs.open(fileRef, "r")
finList = codecs.open(fileInList, "r")
fout = codecs.open(fileOut, "w")



#Initialize a dictionary 
#platax_transcriptome[transcript_id]=nucleotide sequence
platax_transcriptome=dict()

#Read fasta file
i=0
seq=""
line=fref.readline()
while line :
	if line[0] == ">" :
		if i!=0 :
			platax_transcriptome[str(tr_id)]=seq
		seq=""
		tr_id=line.replace(">","").split(" ",1)[0]
		i+=1
	else :
		seq=seq+line.replace("\n","")
	line=fref.readline()
platax_transcriptome[str(tr_id)]=seq
		
line=finList.readline()
while line :
	line=line.replace("\n","")
	fout.write(">"+line+"\n"+str(platax_transcriptome[str(line)])+"\n")
	line=finList.readline()




