# Demultiplexing and Index Swapping

**Goals**: Determine the level of index swapping and undetermined index-pairs, before and after quality filtering, from sequencing done by the 2017 BGMP cohort. 

The data must first be demultiplexed as demultiplexing is necessary for downstream analyses.

The goal of the `demux.py` script is to demultiplex samples to create 48 FASTQ files that contain the acceptable index pairs (read1 and read2 for 24 different index pairs), and 2 FASTQ files with index-hopped read-pairs, and 2 FASTQ files containing undetermined (non-matching or low quality_ index-pairs)


There are 24 indexed (dual matched) libraries for sequencing. The indexes are:

```
B1	GTAGCGTA    A5	CGATCGAT    C1	GATCAAGG
B9	AACAGCGA    C9	TAGCCATG    C3	CGGTAATC
B3	CTCTGGAT    C4	TACCGGAT    A11	CTAGCTCA
C7	CACTTCAC    B2	GCTACTCT    A1	ACGATCAG
B7	TATGGCAC    A3	TGTTCCGT    B4	GTCCTAAG
A12	TCGACAAG    C10	TCTTCGAC    A2	ATCATGCG
C2	ATCGTGGT    A10	TCGAGAGT    B8	TCGGATTC
A7	GATCTTGC    B10	AGAGTCCA    A8	AGGATAGC
```

The 4 FASTQ files are: 
```bash
1294_S1_L008_R1_001.fastq.gz
1294_S1_L008_R2_001.fastq.gz
1294_S1_L008_R3_001.fastq.gz
1294_S1_L008_R4_001.fastq.gz
```

## Part 1 – Quality Score Distribution per-nucleotide


## Part 2 – Demultiplexing algorithm

- one R1 FASTQ file and one R2 FASTQ file **per** matching index-pair, 
- another two FASTQ files for non-matching index-pairs (index-hopping)
- two additional FASTQ files when one or both index reads are unknown or low quality (do not match the 24 known indexes [this includes indexes with 'N's in them] or do not meet a quality score cutoff)
    
Add the sequence of the index-pair to the header of BOTH reads in all of your FASTQ files for all categories (e.g. add “AAAAAAAA-CCCCCCCC” to the end of headers of every read pair that had an index1 of AAAAAAAA and an index2 of CCCCCCCC; this pair of reads would be in the unknown category as one or both of these indexes do not match the 24 known indexes).

Additionally, your algorithm should report: 
- the number of read-pairs with properly matched indexes (per index-pair), 
- the number of read pairs with index-hopping observed, and
- the number of read-pairs with unknown index(es).

