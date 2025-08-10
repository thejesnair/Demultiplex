#!/usr/bin/env python

import bioinfo
import argparse
import gzip

'''
demux.py demultiplexes data from 4 FASTQ files (2 biological reads and 2 indexes)
It will output a FASTQ file corresponding to every matched index pair for each biological read (totalling 48 FASTQ files)
as well as 2 FASTQ files per bio read for hopped and unknown indexes
A summarized report of findings (including totals and percentages) will be generated and outputted as a .md file as well as .txt
files for hopped and matched index pairs
'''

'''
CHECKLIST:
    Fxn: Argparse
    Argparse variables
    Counters
    index set
    FASTQ dict
    Fxn: read FASTQ record, assign to each variable
    Fxn: assign indexes to a set for membership checking
    Fxn: open and prepare FASTQ output files for demux pipeline (52 files total)
    Fxn: close files

    output:
    number of read pairs with unknown index pairs
    number of read pairs hopped index pairs
    number of read pairs that are matched
    total number of reads
    hopped/mapped txt files with stats SORTED
    results.md with summary of stats SORTED
    
'''

def get_args():
    '''Takes in arguments for 4 FASTQ files, known indexes: bio reads 1,2 and index reads 1,2 for demultiplexing'''
    parser = argparse.ArgumentParser(description="Program to demultiplex sequencing data")
    parser.add_argument("-r1", "--read1", help="bio read 1 (R1) pathway/filename", type=str, required=True)
    parser.add_argument("-r2", "--read2", help="index read 1 (R2) pathway/filename", type=str, required=True)
    parser.add_argument("-r3", "--read3", help="index read 2 (R3) pathway/filename", type=str, required=True)
    parser.add_argument("-r4", "--read4", help="bio read 2 (R4) pathway/filename", type=str, required=True)
    parser.add_argument("-i", "--indexfile", help="pathway/filename of known indexes file", type=str, required=True)
    parser.add_argument("-o", "--outputdir", help="pathway/filename output dir for generated txt files", type=str, required=True)

    return parser.parse_args()

args = get_args()
r1 = args.read1
r2 = args.read2
r3 = args.read3
r4 = args.read4
indexfile = args.indexfile
outputdir = args.outputdir

#Counters
matched_count = 0
hopped_count = 0
unknown_count = 0

#Dict for txt files/results.md
matched_pairs = {}
hopped_pairs = {}

#Read FASTQ records
def FASTQ_record(fq: str) -> str:
    ''' Read FASTQ file and assign each line of record to appropriate variable
requires file to already be opened '''
    header = fq.readline().strip()  #read one line at a time, strip newline
    seq = fq.readline().strip()
    plus = fq.readline().strip()
    QS = fq.readline().strip()

    return header,seq,plus,QS

#Get 24 known indexes, assign to set
def getIndexSet(indexfile: str) -> set:
    ''' Reads index file, grabs fifth column indexes and assigns to a set '''
    with open(indexfile, "r") as f:
        next(f)
        index_set = {(line.strip().split()[4]) for line in f} #strips 5th col, splits into list of strs
    return(index_set)

index_set = getIndexSet(indexfile)

#Create files to write the 48 FASTQ files for matched index pairs
FASTQ_files = {}

def open_files(index_set: set, outputdir: str) -> dict:
    ''' Opens and prepares files to be written for matched index pairs (48 total), hopped (2), and unknown (2)
    key: type of file, value: file names
    R1, R2 within loop because file for each index and each read (24 indexes for R1,R2 = 48)
    Hopped and unknown OUTSIDE of loop because only 2 files total for each'''
    for i in index_set:
        R1 = open(outputdir + i + "_R1.fq", "w") 
        R2 = open(outputdir + i + "_R2.fq", "w")
        FASTQ_files[i] = (R1, R2) #where i is for index in index_set #0:R1, 1:R2, careful with indenting
    FASTQ_files["hopped"] = open(outputdir + "hopped_R1.fq", "w"), open(outputdir + "hopped_R2.fq", "w") #2:R1 3:R2
    FASTQ_files["unknown"] = open(outputdir + "unknown_R1.fq", "w"), open(outputdir + "unknown_R2.fq", "w") #4:R1, 5:R2
    return FASTQ_files


FASTQ_files = open_files(index_set, outputdir)

#Close files!
def close_files(FASTQ_files: dict) -> None:
    '''Closes all files in the FAST_files dict'''
    for key, (f1, f2) in FASTQ_files.items(): #every value is a tuple of two open file handles, r1 and r2
        f1.close()
        f2.close()

#Open files to start parsing through records and determining if indexes match
with gzip.open(r1, "rt", encoding='utf-8') as r1, gzip.open(r2, "rt", encoding='utf-8') as r2, \
gzip.open(r3, "rt", encoding='utf-8') as r3, gzip.open(r4, "rt", encoding='utf-8') as r4:

    while True:
        header1, seq1, p1, QS1 = FASTQ_record(r1) #assigns record lines to variables
        header2, seq2, p2, QS2 = FASTQ_record(r2) #only need seq2
        header3, seq3, p3, QS3 = FASTQ_record(r3)
        header4, seq4, p4, QS4 = FASTQ_record(r4)
        
        if not header1: #EOR/EOF: end of reading records for all reads in all files, break loop
            break

        seq3_RC = bioinfo.reverse_complement(seq3)  #need RC of index 2 to check against index 1
        index_pair = f"{seq2}-{seq3_RC}"
        reverse_indexpair = f"{seq3_RC}-{seq2}"

        #Dual matched
        if seq2 == seq3_RC and seq2 in index_set:   #if index matched AND known 
            matched_count += 1
            FASTQ_files[seq2][0].write(f"{header1} {index_pair}\n{seq1}\n{p1}\n{QS1}\n") #index1,r1
            FASTQ_files[seq2][1].write(f"{header4} {reverse_indexpair}\n{seq4}\n{p4}\n{QS4}\n") #index2,r2
            matched_pairs[index_pair] = matched_pairs.get(index_pair, 0) +1 #if no value assoc, increments by 1 for index pair

        #Hopped
        elif seq2 != seq3_RC and seq2 in index_set and seq3_RC in index_set:  #if NOT index matched AND known index
            hopped_count += 1
            FASTQ_files['hopped'][0].write(f"{header1} {index_pair}\n{seq1}\n{p1}\n{QS1}\n")
            FASTQ_files['hopped'][1].write(f"{header4} {reverse_indexpair}\n{seq4}\n{p4}\n{QS4}\n")
            hopped_pairs[index_pair] = hopped_pairs.get(index_pair, 0) +1 #if no value assoc, increments by 1 for index pair

        #Unknown
        else:
            unknown_count += 1
            FASTQ_files['unknown'][0].write(f"{header1} {seq2}-{seq3_RC}\n{seq1}\n{p1}\n{QS1}\n")
            FASTQ_files['unknown'][1].write(f"{header4} {seq3_RC}-{seq2}\n{seq4}\n{p4}\n{QS4}\n")

close_files(FASTQ_files)

#Percent variables
total_counts = matched_count + hopped_count + unknown_count
matched_percent = (matched_count/total_counts) * 100
hopped_percent = (hopped_count/total_counts) * 100
unknown_percent = (unknown_count/total_counts) * 100

print(f"Total reads: {total_counts}")
print(f"Matched index read pairs: {matched_count}")
print(f"Hopped index read pairs: {hopped_count}")
print(f"Unknown index read pairs: {unknown_count}")

#Write hopped, matched indexpairs to txt file
with open("matched_index_pairs.txt", "w") as matched_output:
    matched_output.write(f'Index Pair\tCount\tPercentage of Reads\n')
    for indexpair, count in sorted(matched_pairs.items(), key=lambda x:x[1], reverse=True): 
        #remember lambda is just way to look at pattern, x is for (indexpair, count) x[1] is saying to look at the count value
        #reverse=True, descending orders
        matched_output.write(f"{indexpair}\t{count}\t{matched_percent:.2f}\n")

with open("hopped_index_pairs.txt", "w") as hopped_output:
    hopped_output.write(f'Index Pair\tCount\tPercentage of Reads\n')
    for indexpair, count in sorted(hopped_pairs.items(), key=lambda x:x[1], reverse=True):
        percent_reads = (count/total_counts) * 100
        hopped_output.write(f"{indexpair}\t{count}\t{hopped_percent:.2f}\n")

#Write summary report to .md for github
with open("results.md", "w") as results:
    results.write('# Demultiplexing Report #\n')
    results.write(f'**Total Reads**: {total_counts}\n\n')
    results.write(f'**Matched indexes** : {matched_count} ({matched_percent:.2f}%)\n\n')
    results.write(f'**Hopped indexes**: {hopped_count} ({hopped_percent:.2f}%)\n\n')
    results.write(f'**Unknown indexes**: {unknown_count} ({unknown_percent:.2f}%)\n')

    results.write("## Matched Indexes ##\n\n")
    results.write("| Matched Index Pair | Count | Percentage Reads (Total Matched) |\n")
    results.write("|------------|-------|---------------------|\n")
    for indexpair, count in sorted(matched_pairs.items(), key=lambda x:x[1], reverse=True):
        percent_reads = (count/matched_count) * 100
        results.write(f'{indexpair} | {count} | {percent_reads:.2f}% |\n')

    results.write("## Hopped Indexes ##\n\n")
    results.write("| Hopped Index Pair | Count | Percentage Reads (Total Hopped) |\n")
    results.write("|------------|-------|---------------------|\n")
    for indexpair, count in sorted(hopped_pairs.items(), key=lambda x:x[1], reverse=True):
        percent_reads = (count/hopped_count) * 100
        results.write(f'{indexpair} | {count} | {percent_reads:.2f}% |\n')
