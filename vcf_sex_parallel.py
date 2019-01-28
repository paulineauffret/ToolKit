#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Input: vcf file from freebayes

# Output: vcf-like format file, with "CHROM POS REF ALT raw p-values (Chi2 test)" fields.

# Help: python get_AF_from_vcf_v4.py --help

#AUTHOR : LEO MILHADE (milhadel)

###################################################################
#Import libraries
import sys
import os
import numpy as np

import argparse

#from FisherExact import fisher_exact
from scipy.stats import chisquare

import multiprocessing as mp

import warnings
warnings.filterwarnings("error")

###################################################################
# Define parser
parser = argparse.ArgumentParser(description='Compute significance scores from a vcf file')

parser.add_argument('input_file', type=str,
                   help='vcf file, freebayes format')

parser.add_argument('-o', '--output', dest='output_file', default="parse_vcf.out",
                    help='output file name')
parser.add_argument('-t', '--threads', dest='processes', default=2, type = int,
                    help='number of worker processes to use')
parser.add_argument('--chunk_size', dest='chunk_size', default=1024*1024, type = int,
                    help='size of the chunk of lines read in the input file')

args = parser.parse_args()

###################################################################
#Set parameters
###################################################################
#Set parameters
#original_column_names = ["CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO",
                         #"FORMAT","pool9KV","pool6GV","pool4GR","pool8KJ","pool2TJ",
                         #"pool10ER","pool5GJ","pool1TR","pool7KR","pool3TV",
                         #"pool11EJ","pool12EV"]
original_column_names = ["CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","yellow","green","red"]
original_column_names = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "male1", "male2", "poolFemales", "male3", "male4","male5","male6","male7", "male8","male9","male10","male11","male12"]
column_names = ["CHROM","POS","REF","ALT","QUAL","Female missing", "Male missing", "raw p-values (Chi2 contingency test)"]

comment_char       = "#"        #comment character for header in vcf file
split_char         = "\t"       #sample field separator in vcf file
nb_samples         = 13		#number of samples in vcf file
indx_columns_start = 9         #index of colum for the first sample if colum number 1 = CHROM
nbcol              = int(indx_columns_start)+int(nb_samples)
field_sep          = ":"        #for exemple in GT:DP:AD:RO:QR:AO:QA:GL field separator = ":"
DP_field           = 2          #for exemple in GT:DP:AD:RO:QR:AO:QA:GL DP field = 2
AD_field           = 3          #for exemple in GT:DP:AD:RO:QR:AO:QA:GL AD field = 3
GT_field           = 0          #for exemple in GT:DP:AD:RO:QR:AO:QA:GL AD field = 3
AD_field_sep       = ","        #separator in AD field
#alpha              = 0.05           #adjusted p-value cutoff
#adj_method         = "fdr_bh"   #p-value correction method 

end_signal = "PROCESSING DONE"

###################################################################

def process(line):

    line = line.split(split_char)

    #Get chrom, pos, ref and alt fields

    newline = split_char.join(line[i] for i in [0, 1, 3, 4, 5])

    male_het_count = 0
    male_count = 0
    female_dp = 0
    female_ad = 0
    female_missing = 0
    male_missing = 0

    for index in range(int(indx_columns_start),nbcol) :
        sample_name = original_column_names[index]
        section = line[index].split(field_sep)
        if sample_name == "poolFemales":
            if section[0] == ".":
                female_missing = 1
            else:
                female_dp = int(section[DP_field])
                female_ad = max(int(ad) for ad in section[AD_field].split(AD_field_sep))
        else : #sample_name == "male."
            if section[0] == ".":
                male_missing += 1
            else:
                male_het_count += section[GT_field] == "0/1"
                male_count += 1

    try:
        raw_pval = chisquare([female_ad, male_het_count], f_exp=[female_dp, male_count])[1]
    except ValueError, e:
        raw_pval = np.nan 
    except RuntimeWarning, e:
        raw_pval = np.nan 
    
    #newline = split_char.join((newline, str(female_missing), str(male_missing), str([female_ad, male_het_count, female_dp, male_count]), str(raw_pval)))
    newline = split_char.join((newline, str(female_missing), str(male_missing), str(raw_pval)))

    return newline

###################################################################

def write_header():
    with open(args.output_file, 'w') as out_file:
        header = split_char.join(column_names)
        out_file.write(header + "\n")
        out_file.flush()

def write_line(queue):
    with open(args.output_file, 'a+') as f:
        while True:
            line = queue.get()
            if line == end_signal:
                break
            else:
                f.write(line + "\n")
                f.flush()

def process_wrapper(chunk_start, chunk_size, queue):
    with open(args.input_file, 'rU') as in_file:
        in_file.seek(chunk_start)
        lines = in_file.read(chunk_size).splitlines()
        for line in lines:
            queue.put(process(line))

def skip_comments(f):
    line = f.readline()
    while line.startswith(comment_char):
        chunk_start = f.tell()
        line = f.readline()
    return f, chunk_start

def chunkify(fname, size=args.chunk_size):
    file_end = os.path.getsize(fname)
    with open(args.input_file,'r') as f:
        chunk_end = f.tell()
        
        f, chunk_end = skip_comments(f)
        

        while True:
            chunk_start = chunk_end
            f.seek(size,1)
            f.readline()
            chunk_end = f.tell()
            yield chunk_start, chunk_end - chunk_start
            if chunk_end > file_end:
                break

def main():
    
    write_header()

    manager = mp.Manager()
    queue = manager.Queue()
    pool = mp.Pool(args.processes)

    watcher = pool.apply_async(write_line, (queue, ))

    jobs = []

    for chunk_start, chunk_size in chunkify(args.input_file):
        job = pool.apply_async(process_wrapper, (chunk_start, chunk_size, queue))
        jobs.append(job)

    for job in jobs:
        job.get()

    queue.put(end_signal)

    pool.close()

if __name__ == "__main__":
    main()
