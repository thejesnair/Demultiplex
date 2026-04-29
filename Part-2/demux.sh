#!/bin/bash
#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --cpus-per-task=8
#SBATCH --mem=16GB
#SBATCH --job-name=Demux
#SBATCH --output=Demux%j

index="/projects/bgmp/shared/2017_sequencing/indexes.txt"
R1="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz"
R2="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz"
R3="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz"
R4="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz"
outputdir="/projects/bgmp/tnair/bioinfo/Bi622/Demultiplex/DemuxP3_FASTQ/"


#Run demux.py
/usr/bin/time -v ./demux.py -r1 $R1 -r2 $R2 -r3 $R3 -r4 $R4 -i $index -o $outputdir

