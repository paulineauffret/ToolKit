#!/bin/bash
			
###################################################################
#Pipeline to evaluate if a SNP is leading to a synonymous mutation.
#######Input : 
#[1] 	: SNP info text file (subset from vcf file, ie 4 columns : CHROM	POS	REF	ALT)
#[2]	: orf fasta file ('longest_orfs.cds' obtained from Transdecoder.LongOrfs)

				
#######Output : a text file with 15 columns and one line per ORF. Each SNP is evaluated on all the ORF existing on the same transcript/gene.
#CHROM 		: transcript / gene id ;
#ORF_start 	: start position of the current orf ;
#ORF_end 	: end position of the current orf ;
#syn 		: either "SYN" if the mutation is synonymous ; "NON SYN" if the mutation is non synonymous, and "NOT IN ORF" if the mutation is not in the current orf ;
#strand		: + or - strand ;
#orf_nt_seq	: nucleotide orf sequence with REF nucleotide ('NA' if syn='NOT IN ORF') ;
#orf_prot_seq	: translated nucleotide orf sequence with REF nucleotide ('NA' if syn='NOT IN ORF') ;
#mut_nt_seq	: nucleotide orf sequence with ALT nucleotide ('NA' if syn='NOT IN ORF') ;
#mut_prot_seq	: translated nucleotide orf sequence with ALT nucleotide ('NA' if syn='NOT IN ORF') ;
#aa_ref		: translated codon in current orf with REF nucleotide ('NA' if syn='NOT IN ORF') ;
#aa_mut 	: translated codon in current orf with ALT nucleotide ('NA' if syn='NOT IN ORF') ;
#CHROM		: transcript / gene id ;
#POS		: evaluated snp position ;
#REF		: evaluated reference nucleotide ;
#ALT		: evaluated alt nucleotide.
				
###################################################################
#Set input/output files and directories
working_directory=									#working/output directory
get_prot_seq_script=									#python script 'get_prot_seq_custom.py'
snp_file=${working_directory}/info_file_snp.txt 					#SNP info text file (subset from vcf file, ie 4 columns : CHROM	POS	REF	ALT)
tmp_file=${working_directory}/tmp							#a temporary file
orf_file=${working_directory}/longest_orfs.cds						#orf fasta file ('longest_orfs.cds' obtained from Transdecoder.LongOrfs)
out_file=${working_directory}/output.txt						#output file

#Print header line in output file
echo "CHROM	ORF_start	ORF_end	syn	strand	orf_nt_seq	orf_prot_seq	mut_nt_seq	mut_prot_seq	aa_ref	aa_mut	CHROM	POS	REF	ALT" > $out_file

#Get the number of snp
nb_line=$(awk -F, 'END{print NR}' $snp_file) 

#Run the python script to evaluate each snp 
for (( i=1; i<=$nb_line; i++)) ; 
do 
	head -$i $snp_file | tail -1 > $tmp_file ; 
	python2 $get_prot_seq_script $orf_file $tmp_file $out_file ;
done

#Remove temporary file
rm $tmp_file


