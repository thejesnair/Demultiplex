## Demultiplex Lab Notebook ##


*GOAL:*<br>
Perform initial data exploration on four fastq files. Write py script to create QS distribution plot of average QS per base position. Write pseudocode for demultiplexing reads & indexes and assign to appropriate fastq files (dual-matched, hopped, unknown). 

Files found in :
```/projects/bgmp/shared/2017_sequencing/```

**FILES**<br>
- R1_001.fastq.gz : Biological R1<br>
- R2_001.fastq.gz : Index1 (for bioR1)<br>
- R3_001.fastq.gz : Index2 (for bioR2)<br>
- R4_001.fastq.gz : Biological R2<br>

For initial data exploration

```
srun -A bgmp -p bgmp --time=1:00:00 --pty bash
zcat 1294_S1_L008_R1_001.fastq.gz | head    #view data
zcat 1294_S1_L008_R1_001.fastq.gz | head -2 | tail -1 | wc #how long seq is
```
seq is 101 NTs long (102 - \n)<br>
indexes are 8 NTs long (9 - \n)

Determine number of Ns in index files
```
zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz| awk 'NR % 4 == 2 && /N/'| wc -l
```
R2_001.fastq.gz : 3976613<br>
R3_001.fastq.gz : 3328051<br>
Total : 7304664


Determine Phred encoding, can do through bash commands!<br>
- Initially was trying to use sed but decided to use awk for practice<br>
- Grabbing JUST ONE quality score from each read, want to get a general idea.<br>
- od -An -t u1 will converts ASCII to number values! *** very handy
    - looked for command that would allow me to calcualte phred scores and found this
    - quick way to get a feel for the range of scores
        - if Phred+33, scores will be in range 33-73~
        - if Phred+64, scores will be in range 64-106~

```
zcat 1294_S1_L008_R1_001.fastq.gz | awk 'NR % 4 ==0' | head -1 | od -An -t u1
```
Copy of bioinfo.py to Demultiplex dir on Talapas

Remember this?<br>
Done in terminal on comp<br>

scp'd bioinfo.py AND test.fastq used in PS4 so I wouldn't have to make new test file
```
scp bioinfo.py tnair@login.talapas.uoregon.edu:/projects/bgmp/tnair/bioinfo/Bi622/Demultiplex

scp test.fastq tnair@login.talapas.uoregon.edu:/projects/bgmp/tnair/bioinfo/Bi622/Demultiplex
```

**Slurm script outputs**
-

```
    Command being timed: "./QS_distr.py -f /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz -o 1294_S1_L008_R2_001.fastq.gz.png -l 8"
	User time (seconds): 647.15
	System time (seconds): 1.18
	Percent of CPU this job got: 99%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 10:50.28
    Maximum resident set size (kbytes): 66888
```

```
    Command being timed: "./QS_distr.py -f /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz -o 1294_S1_L008_R3_001.fastq.gz.png -l 8"
	User time (seconds): 641.44
	System time (seconds): 1.29
	Percent of CPU this job got: 99%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 10:44.47
    Maximum resident set size (kbytes): 63096

```

```
    Command being timed: "./QS_distr.py -f /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz -o 1294_S1_L008_R1_001.fastq.gz.png -l 101"
	User time (seconds): 4017.38
	System time (seconds): 8.75
	Percent of CPU this job got: 99%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 1:07:16
    Maximum resident set size (kbytes): 67868
```

```
    Command being timed: "./QS_distr.py -f /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz -o 1294_S1_L008_R4_001.fastq.gz.png -l 101"
	User time (seconds): 4124.66
	System time (seconds): 8.64
	Percent of CPU this job got: 99%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 1:09:07
    Maximum resident set size (kbytes): 71248

```

Some troubleshooting with git:

Had issues git pushing some changes, probably because I was cleaning up and removing things from github directly (was reorganizing and wanted to avoid issue I had with PS7 where couldn't push at all)

This created diverged branches and wasn't able to push new files...great

used:<br>
```git stash```
to basically "hide" the changes that had been made/converged branches?

then able to do<br>
```git pull --rebase origin master```

then could git add, commit, push files!

But take more caution next time to avoid this issue again, similar to PS7 but different?




### PART2 ###
**Pseudocode**

Strategy for algorithm to demultiplex files and report index hopping.
Given 4 input FASTQ files (2 with biological reads, 2 with index) and 24 known indexes, demultiplex reads by index-pair and output these FASTQ files:


- per matching index-pair (24 indexes)
- non-matching index-pairs (index-hopping)
- unknown/low quality (do not match 24 known indexes [includes indexes with 'N'] or do not meet QS cutoff)

one R1 FASTQ file and one R2 FASTQ file per matching index-pair,
another two FASTQ files for non-matching index-pairs (index-hopping), and
two additional FASTQ files when one or both index reads are unknown or low quality (do not match the 24 known indexes [this includes indexes with 'N's in them] or do not meet a quality score cutoff)

Notes from peer review:
- need to define quality score check, define threshold (30)
    - feedback from Anna and Emma, Emma explained that barcodes with too low of an avg QS should go into the unknown output
- Emma: switch order of if, else statements so that checking unknown comes first, that way won't accidentally put low quality barcodes and reads as dual paired output
    - i had done this prev and then switched the order for some reason, so dont overthink it
- Annika mentioned same thing with if/else order
    - after checking if indexes are valid, THEN sort unknown into unknown file before checking for hopped, matched. 
- E mentioned it's unclear if using boolean statements for if statements or nested if statements (either way could probably work)
- Annika: With large files, defining standalone fxn to handle writing FASTQ won't work. Can only parse one line at a time since files are so large, so need to be written concurrently while parsing. 
- Only need to check if indexes are valid once so...<br><br>
   Instead of:<br> ```Else if I1 and I2_RC do not match AND are still in the list of known indexes AND do not contain 'N' ```<br>
    
   Can do: <br> ```Else if I1 and I2_RC do not match  AND do not contain 'N'```

### Feedback from P1 and P2 ###

- Need to adjust unit tests, however many records there are need to have same number of outputs
    - So minimum records needed is 3 to accomplish dual matched, index hopped, and unknown
- Can fix pseudocode if it helps me but not necessary
- Additionally, needed to put in p2 helpful outputs
    - so printing total number of records of dual matched, hopped, unknown
    - any figures
    - whatever else may be helpful

### PART3 ###
**Writing script**

- Did not implement QS cutoff/quality filter data
- demux.py to demultiplex data
    - created randomtests.py to test smaller bits of code while creating algorithim
    - was able to test getIndexSet fxn using indexes.txt
    - testdemuxscript.py edited script *slightly* to work with test files
    - test output dir "demuxTEST" to see if correct output files
- tried to use try: except: but didn't work, endless loop
    - used if not header1 to break loop

### Talapas Issues ###
8/8/25:
- Talapas has been incredibly slow, not sure what's going on, both login4 and login.

8/9/25:
- Having issues with Talapas again, even slower than last night. Tried adding new login node to the remote SSH, but then couldn't log into Talapas at all. 
    - Tried adjusting the .ssh/config file (?) but that made it worse
    - Had to comment out the commands with ____ and uncomment top
        - I know Jason helped me with this earlier in 621 when I was having issues, so I switched it back to what it was before.
        - No idea why all of a sudden this became an issue, but it seems to be working for now.

    ```bash
    Host *
    ServerAliveInterval 60 # Send a "keep-alive" package every 60 secs
    ServerAliveCountMax 480 # Keep connection alive for 8 hours (60s*480=28,800s)
    TCPKeepAlive yes

    # Host *
    # ServerAliveInterval 240# Send a "keep-alive" package every ___ secs
    # ServerAliveCountMax 240# Keep connection alive for ___ hours
    # TCPKeepAlive yes

    ```
# Part3 Script #

**Scripts:**
- demux.py
- demux.sh

- testdemuxscript.py (for running testfiles/troubleshooting)

### Troubleshooting ###  
- Rerun test script an error pops up:

```bash
(base) [tnair@login1 Demultiplex]$ ./testdemuxscript.py 
Traceback (most recent call last):
  File "/gpfs/projects/bgmp/tnair/bioinfo/Bi622/Demultiplex/./testdemuxscript.py", line 146, in <module>
    FASTQ_files[seq2][0].write(f"{header1} {seq2}-{seq3_RC}\n{seq1}\n{p1}\n{QS1}\n") #index1,r1
    ~~~~~~~~~~~^^^^^^
KeyError: 'ATCATCAG'
```

This error only happens when I rerun script, but if I run it initially there's no issues--files are outputted correctly and contents are correct.

When running actual data, didn't output anything into the files

Identified the issue:
```python
    for i in index_set:
        R1 = open(outputdir + i + "_R1.fq", "w") 
        R2 = open(outputdir + i + "_R2.fq", "w")
    *FASTQ_files[i] = (R1, R2) #where i is for index in index_set #0:R1, 1:R2
    FASTQ_files["hopped"] = open(outputdir + "hopped_R1.fq", "w"), open(outputdir + "hopped_R2.fq", "w") #2:R1 3:R2
    FASTQ_files["unknown"] = open(outputdir + "unknown_R1.fq", "w"), open(outputdir + "unknown_R2.fq", "w") #4:R1, 5:R2
    return FASTQ_files
```
```
*needs to be indented
with the way its written now, is rewriting file--everything was going into unknown output
```

Another issue:
Script didn't write index pairs, count, percentage to results.md
- didn't populate list in if/else loop


### Running script ###

Demux37135739
```
Command being timed: "./demux.py -r1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz -r2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz -r3 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz -r4 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz -i /projects/bgmp/shared/2017_sequencing/indexes.txt -o /projects/bgmp/tnair/bioinfo/Bi622/Demultiplex/DemuxP3_FASTQ/"
User time (seconds): 2449.43
System time (seconds): 90.19
Percent of CPU this job got: 96%
Elapsed (wall clock) time (h:mm:ss or m:ss): 43:56.75
Maximum resident set size (kbytes): 244840
Exit status: 0
```

Quick summary stats of demux
```
Total reads: 363246735
Matched index read pairs: 331755033
Hopped index read pairs: 707740
Unknown index read pairs: 30783962

```

*Notes after completing script:*

- Wanted to write results.md file in same .py script, but bc of troubleshooting had to rerun code several times
    - Could have separated into two scripts, or done a more thorough run through with test files. Or created bigger test files.
- Having a separate test .py was helpful but have to be careful in making sure changes in test are copied over to actual .py

**demux.py summary in results.md**