#!/usr/bin/env python
# -*- coding: utf-8 -*-

#Input : 1: vcf file from freebayes /!\ update parameters in this script from line 21 to line 38 if needed. 
#	2: alpha (pvalue cutoff)
#	3: pvalue correction method
#	4: temporary output file name
#	5: final output file name
	
#Output : vcf-like format file, with "CHROM POS REF ALT raw p-values (KW test)	adjusted p-values" fields. Raw p-value from Kruskal Wallis test. 

#Usage : python get_AF_from_vcf_v2.py input_test.vcf 0.05 fdr_bh out out2 

###################################################################
#Import libraries
import sys
import numpy
import scipy.stats
import statsmodels.stats.multitest as multi

###################################################################
#Set parameters
original_column_names=["CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","pool9KV","pool6GV","pool4GR","pool8KJ","pool2TJ","pool10ER","pool5GJ","pool1TR","pool7KR","pool3TV","pool11EJ","pool12EV"]
column_names=["CHROM","POS","REF","ALT","raw p-values (KW test)","adjusted p-values"]
sample_info={"pool9KV":"V","pool6GV":"V","pool3TV":"V","pool12EV":"V","pool4GR":"R","pool10ER":"R","pool1TR":"R","pool7KR":"R","pool8KJ":"J","pool2TJ":"J","pool5GJ":"J","pool11EJ":"J"}


comment_char="#"	#comment character for header in vcf file
split_char="\t"		#sample field separator in vcf file
nb_samples=12		#number of samples in vcf file
indx_columns_start=10	#index of colum for the first sample if colum number 1 = CHROM
field_sep=":"		#for exemple in GT:DP:AD:RO:QR:AO:QA:GL field separator = ":"
DP_field=2		#for exemple in GT:DP:AD:RO:QR:AO:QA:GL DP field = 2
AD_field=3		#for exemple in GT:DP:AD:RO:QR:AO:QA:GL AD field = 3
AD_field_sep=","	#separator in AD field
#alpha=0.05		#adjusted p-value cutoff
#adj_method="fdr_bh"	#p-value correction method 

###################################################################
#Set input files and command-line parameters
vcf=sys.argv[1] 			#vcf file, freebayes format
alpha=float(sys.argv[2])		#adjusted p-value cutoff
adj_method=sys.argv[3]			#p-value correction method {'b': 'Bonferroni',
                         		#'s': 'Sidak',
                           		#'h': 'Holm',
                           		#'hs': 'Holm-Sidak',
                           		#'sh': 'Simes-Hochberg',
                           		#'ho': 'Hommel',
                           		#'fdr_bh': 'FDR Benjamini-Hochberg',
                          		#'fdr_by': 'FDR Benjamini-Yekutieli',
                           		#'fdr_tsbh': 'FDR 2-stage Benjamini-Hochberg',
                           		#'fdr_tsbky': 'FDR 2-stage Benjamini-Krieger-Yekutieli',
                           		#'fdr_gbs': 'FDR adaptive Gavrilov-Benjamini-Sarkar'}

output=sys.argv[4]                      #temprary output file, will be deleted
output2=sys.argv[5]                     #final output file

#Open files
try:
	fileIN=open(vcf, "rU")
except IOError, e:      		
	print "Error : cannot open file ", vcf

try:
	fileOUT=open(output, "w")
except IOError, e:      		
	print "Error : cannot open file ", output

try:
        fileOUT2=open(output2, "w")
except IOError, e:
        print "Error : cannot open file ", output2

###################################################################
#Write header on new file
for elem in range(len(column_names)) :
	fileOUT2.write(str(column_names[elem])+split_char)
fileOUT2.write("\n")

#initialize raw-pvalues list 
rawp=[]

#Get number of columns
nbcol=int(indx_columns_start)+int(nb_samples)

#Loop on vcf file to get allele frequencies
line=fileIN.readline()
while line :
	if not line.startswith(comment_char) and not len(line) <= 1:
		line=line.split(split_char)
		#Get chrom, pos, ref and alt fields
		newline=line[0]+split_char+line[1]+split_char+line[3]+split_char+line[4]
		#Initialize 3 lists for allele frequencies
		vert=[]
		rouge=[]
		jaune=[]
		for index in range(int(indx_columns_start),nbcol) :
			missing=0
			sample_name=original_column_names[index-1]
			color=sample_info[sample_name]
   			section=line[index-1]
			DP=section.split(field_sep)[int(DP_field)-1].replace("\n","")
			if str(DP) == "." :
				missing=1	#to do : missing case
			#Compute Allele frequencies
			else :
				ADr=section.split(field_sep)[int(AD_field)-1].split(AD_field_sep)[0].replace("\n","")
				ADa=section.split(field_sep)[int(AD_field)-1].split(AD_field_sep)[1].replace("\n","")
				AFr=round(float(ADr)/float(DP),5)
				AFa=round(float(ADa)/float(DP),5)
				if color == "V" :
					vert.append(AFr)
				elif color == "R" :
					rouge.append(AFr)
				else :
					jaune.append(AFr)
		#Kruskal Wallis test
		kw=scipy.stats.kruskal(vert,jaune,rouge)
		rawp.append(kw[1])
		newline=newline+split_char+str(kw[1])
		#Write in temporary file
		fileOUT.write(str(newline)+"\n")
	line=fileIN.readline()

#Close files
fileIN.close()
fileOUT.close()

#Correct p-values
adjp=multi.multipletests(rawp, alpha, adj_method, is_sorted=False, returnsorted=False)[1]

#Reopen temporary output file
try:
        fileOUT=open(output, "r")
except IOError, e:                      
        print "Error : cannot open file ", output

#Adding adjusted pvalue and printing into final output file
i=0
line=fileOUT.readline()
while line :
	line=line.replace("\n","")+split_char+str(adjp[i])+"\n"
	fileOUT2.write(line)
	i+=1
	line=fileOUT.readline()

#Close files
fileOUT.close()
fileOUT2.close()



