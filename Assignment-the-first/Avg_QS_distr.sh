#!/bin/bash
#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --cpus-per-task=8
#SBATCH --job-name=Avg_QS_distr
#SBATCH --output=Avg_QS_distr%j


#Run test.fastq first!
#/usr/bin/time -v ./QS_distr.py -f "/projects/bgmp/tnair/bioinfo/Bi622/Demultiplex/test.fastq" -o "test2.png" -l 101

#Running R1
#/usr/bin/time -v ./QS_distr.py -f "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz" -o "1294_S1_L008_R1_001.png" -l 101

#Running R2
#/usr/bin/time -v ./QS_distr.py -f "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz" -o "1294_S1_L008_R2_001.png" -l 8

#Running R3
#/usr/bin/time -v ./QS_distr.py -f "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz" -o "1294_S1_L008_R3_001.png" -l 8

#Running R4
/usr/bin/time -v ./QS_distr.py -f "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz" -o "1294_S1_L008_R4_001.png" -l 101
