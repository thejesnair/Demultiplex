#!/usr/bin/env python

import bioinfo
import gzip
import argparse
import matplotlib.pyplot as plt


def get_args():
    '''Takes in arguments for filename, output png file name, and list length to produce
    QS distribution plots for zipped fastq files'''
    parser = argparse.ArgumentParser(description="Program to map QS distribution per NT.")
    parser.add_argument("-f", "--filename", help="specify input fastq file/pathway", type=str)
    parser.add_argument("-o", "--outputpng", help="name for output png", type=str)
    parser.add_argument("-l", "--listlength", help="lenght of list", type=int)

    return parser.parse_args()

args = get_args()
filename = args.filename
outputpng = args.outputpng
listlength = args.listlength

def init_list(lst: list, listlength: int, value: float=0.0) -> list:
    '''This function takes an empty list and will populate it with
    the value passed in "value". If no value is passed, initializes list
    with 101 values of 0.0.'''
    return [value] * listlength

my_list: list = []
my_list = init_list(my_list, listlength)

def populate_list(file: str) -> tuple[list, int]:
    """Creates an empty list, reads fastq file, strips line of \n, 
    converts quality score and sums all QS at each position, adds to list"""
    my_list = init_list([], listlength)
    line_count = 0
    char_count = 0
    with gzip.open(filename, 'rt', encoding='utf-8') as fq:
        for line in fq:
            line = line.strip()
            line_count += 1
            if line_count % 4 == 0:
                for i, QS in enumerate(line):
                    my_list[i] += bioinfo.convert_phred(QS)
    return my_list, line_count

my_list, num_lines = populate_list(filename)

num_records = num_lines/4

#Calculate mean qscore at each base
for i,QS in enumerate(my_list):
    my_list[i] = (QS/num_records) 

x = range(listlength)
y = my_list

plt.bar(x, y)
plt.xlabel("Position on Sequence")
plt.title("QS Distribution per Base Position")
plt.ylabel("Avg Quality Score")
plt.savefig(outputpng)

