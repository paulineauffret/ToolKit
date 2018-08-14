#!/usr/bin/env python
# -*- coding: utf-8 -*-

#Import libraries
import sys

###################################################################
#Set input files
fasta=sys.argv[1] 			#transcriptome or genome fasta file
idlist=sys.argv[2]			#IDs list to get fasta sequence (text file; one id per line, no header)
output=sys.argv[3]			#output file

#Open files
try:
	fileIN=open(fasta, "rU")
except IOError, e:      		
	print "Error : cannot open file ", fasta

try:
	fileID=open(idlist, "rU")
except IOError, e:      		
	print "Error : cannot open file ", idlist

try:
	fileOUT=open(output, "w")
except IOError, e:      		
	print "Error : cannot open file ", output

###################################################################
#Initialize dictionary
ref=dict()

#Read fasta file
line=fileIN.readline()
seq=""
id=""
nbSeq=0
while line :
	if line.startswith(">") :
		if nbSeq!=0 :
			ref[str(id)]=str(seq)
		id=line.split(" ")[0].replace(">","")
		seq=""
		nbSeq=nbSeq+1
	else :
		seq=seq+line
	line=fileIN.readline()
ref[str(id)]=str(seq)
print("There are "+str(nbSeq)+" sequences in input transcriptome/genome.")

#print(ref)
#print(len(ref.keys()))

###################################################################
#Read id file
dup=dict()
line=fileID.readline()
nb_id=0
while(line) :
	line=line.replace("\n","").replace(">","")
	if str(line) != "" :
		nb_id=nb_id+1
		if not ref.has_key(line) :
			print("Error : sequence "+str(line)+" not present in input transcriptome/genome.")	
		elif not dup.has_key(line) :
			fileOUT.write(">"+line+"\n")
			fileOUT.write(ref[line])
	dup[str(line)]=1
	line=fileID.readline()
print(str(nb_id)+" sequences printed in "+str(output)+" file.")


