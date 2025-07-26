### GOAL:

This pseudocode details the algorithmic outline that will be used to demultiplex data generated from a previous cohort's library preps

- Create 48 FASTQ files contaiing acceptable index pairs
    - Biological read 1 and biological read 2
- 2 FASTQ files with index-hopped pairs
    - meaning indexes do not match
- 2 FASTQ files with undetermined index-pairs
    - generally indexes containing an N


What files are we working with?

- 24 indexed (dual matched) libraries, indexes.txt
- R1_001.fastq.gz
    - biological Read 1
- R2_001.fastq.gz
    - index 1
- R3_001.fastq.gz
    - index 2
- R4_001.fastq.gz
    - biological Read 2

*Note:* R3 will need to be reverse complemented

Output files:

48 FASTQ matched index pairs<br>
2 FASTQ for index-hopped pairs<br>
2 FASTQ for unknown<br>


**PSEUDOCODE**

Function that will take in a string and return its reverse complement (RC) (This will be used for I2 [R3])<br>

Function that will automatically write contents to a file with appropriate file name<br>
&nbsp;&nbsp;&nbsp;&nbsp;Input filename, filecontents<br>
&nbsp;&nbsp;&nbsp;&nbsp;Extract line1-4 of record<br>
&nbsp;&nbsp;&nbsp;&nbsp;Output FASTQ file

Parse through all files simultaneously(R1, R2, I1, I2)<br>
&nbsp;&nbsp;&nbsp;&nbsp;If I1 is the same as I2_RC AND I1 and I2_RC are in our group of defined indexes<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Then write sequence header and I1-I2_RC (line1), sequence (line2, coming from either R1 or R2), identifier (line3), and QS (line4) to appropriate FASTQ file (based on the known index pair) for R1 and a separate file for R2, using write file fxn<br>
&nbsp;&nbsp;&nbsp;&nbsp;Else if I1 and I2_RC do not match AND are still in the list of known indexes AND do not contain 'N'<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Write to index hopped file for R1 and R2, using write file fxn<br>
&nbsp;&nbsp;&nbsp;&nbsp;Else<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Write to unknown file for R1 and R2, using write file fxn



