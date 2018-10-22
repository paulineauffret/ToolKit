#!/usr/bin/env python
# -*- coding: utf-8 -*-

#Input : vcf file from freebayes /!\ update parameters in this script from line 18 to line 27 if needed 
#Output : vcf-like format file, with "CHROM POS REF ALT" field and then one colum per sample containing 3 values separated with a colon (:). 
#First value is the DP (unchanged from the original VCF file), the second value is the frequencie of the reference allele (AFr).
#Third value is the frequencie of the alternative allele (AFa).
#Example : if the original vcf indicates 0/1:14:12,2:12:488:2:82:-3.56526,0,-40.056 for sample 1, and if 14 is the DP, 12,2 the frequencies of alleles,
#Then AFr=12/14=0.85714 and AFa=2/14=0.14286 rounded to 5 digits. So the field for sample 1 in the custom vcf will be 14:0.85714:0.14286
#/!\ fields with NA:NA:NA indicates missing values in the original vcf file.


#Import libraries
import sys

###################################################################
#Set parameters
column_names=["CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","pool9KV","pool6GV","pool4GR","pool8KJ","pool2TJ","pool10ER","pool5GJ","pool1TR","pool7KR","pool3TV","pool11EJ","pool12EV"]

comment_char="#"	#comment character for header in vcf file
split_char="\t"		#sample field separator in vcf file
nb_samples=12		#number of samples in vcf file
indx_columns_start=10	#index of colum for the first sample if colum number 1 = CHROM
field_sep=":"		#for exemple in GT:DP:AD:RO:QR:AO:QA:GL field separator = ":"
DP_field=2		#for exemple in GT:DP:AD:RO:QR:AO:QA:GL DP field = 2
AD_field=3		#for exemple in GT:DP:AD:RO:QR:AO:QA:GL AD field = 3
AD_field_sep=","	#separator in AD field


###################################################################
#Set input files
vcf=sys.argv[1] 			#vcf file, freebayes format
output=sys.argv[2]			#output file

#Open files
try:
	fileIN=open(vcf, "rU")
except IOError, e:      		
	print "Error : cannot open file ", vcf

try:
	fileOUT=open(output, "w")
except IOError, e:      		
	print "Error : cannot open file ", output

###################################################################
#Write header on new file
for elem in range(len(column_names)) :
	fileOUT.write(str(column_names[elem])+split_char)
fileOUT.write("\n")

#Get number of columns
nbcol=int(indx_columns_start)+int(nb_samples)

#Loop on vcf file to get allele frequencies
line=fileIN.readline()
while line :
	if not line.startswith(comment_char) and not len(line) <= 1:
		line=line.split(split_char)
		newline=line[0]+split_char+line[1]+split_char+line[3]+split_char+line[4]
		for index in range(int(indx_columns_start),nbcol) :
			#print(index)
   			section=line[index-1]
			DP=section.split(field_sep)[int(DP_field)-1].replace("\n","")
			if str(DP) == "." :
				newfield="NA"+field_sep+"NA"+field_sep+"NA"
			else :
				ADr=section.split(field_sep)[int(AD_field)-1].split(AD_field_sep)[0].replace("\n","")
				ADa=section.split(field_sep)[int(AD_field)-1].split(AD_field_sep)[1].replace("\n","")
				AFr=round(float(ADr)/float(DP),5)
				AFa=round(float(ADa)/float(DP),5)
				newfield=str(DP)+field_sep+str(AFr)+field_sep+str(AFa)
			newline=newline+split_char+newfield
			newfield=""
		fileOUT.write(str(newline)+"\n")
	line=fileIN.readline()

#Close files
fileIN.close()
fileOUT.close()

