## Part 2 – Demultiplexing algorithm

- `bioinfo.py` commonly used functions
- `demux.py` demultiplexing script
    - one R1 FASTQ file and one R2 FASTQ file **per** matching index-pair, 
    - another two FASTQ files for non-matching index-pairs (index-hopping)
    - two additional FASTQ files when one or both index reads are unknown or low quality (do not match the 24 known indexes [this includes indexes with 'N's in them] or do not meet a quality score cutoff)
- `demux.sh` SLURM script for use on Talapas HPC
- `results.md` summarized output report containing:
    - the number of read-pairs with properly matched indexes (per index-pair), 
    - the number of read pairs with index-hopping observed
    - the number of read-pairs with unknown index(es)