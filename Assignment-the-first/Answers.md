# Assignment the First

## Part 1
1. Be sure to upload your Python script. Provide a link to it here:<br>
[QS_distr.py](https://github.com/thejesnair/Demultiplex/blob/ce617b2011b562f83f5a5d6ecece875de927f9dd/Assignment-the-first/QS_distr.py)<br>

| File name | label | Read length | Phred encoding |
|---|---|---|---|
| 1294_S1_L008_R1_001.fastq.gz | Biological Read 1 | 101 | 33 |
| 1294_S1_L008_R2_001.fastq.gz | Index 1 | 8 | 33 |
| 1294_S1_L008_R3_001.fastq.gz | Index 2 | 8 | 33 |
| 1294_S1_L008_R4_001.fastq.gz | Biological Read 2 | 101 | 33 |

2. Per-base NT distribution
    1. Use markdown to insert your 4 histograms here.

- 1294_S1_L008_R1_001.fastq.gz.png<br>
![R1 histogram](https://github.com/thejesnair/Demultiplex/blob/master/Assignment-the-first/1294_S1_L008_R1_001.png)

- 1294_S1_L008_R2_001.fastq.gz.png<br>
![I1 histogram](https://github.com/thejesnair/Demultiplex/blob/master/Assignment-the-first/1294_S1_L008_R2_001.png)

- 1294_S1_L008_R4_001.fastq.gz.png<br>
![R2 histogram](https://github.com/thejesnair/Demultiplex/blob/master/Assignment-the-first/1294_S1_L008_R4_001.png)

- 1294_S1_L008_R3_001.fastq.gz.png<br>
![I2 histogram](https://github.com/thejesnair/Demultiplex/blob/master/Assignment-the-first/1294_S1_L008_R3_001.png)


    2. Generally for index reads a quality cutoff of 30 is preferred because it equates to 99.9% accuracy and we require high confidence since we are demultiplexing indexes and assigning them to reads. A single error may cause reads to be incorrectly assigned or tossed which we want to avoid. For biological read pairs, a cutoff of Q20 _is_ the minimum acceptable threshold, and generally acceptable, as it equates to a error rate of 1%. However, a Q30 is preferred if we want a higher sensitivity for analyses. Higher sensitivity would be beneficial for situations like detecting low-frequeny variants, mutation detection, limiting false positives, or variant calling.<br><br>

        Keeping that in mind, and given the results above for the biological reads in these sequencing results we could increase the cutoff to ~35 for the biological reads and indexes to account for the initial lower-quality reads, however I believe Q30 should still be sufficient. 

    3. Determine how many indexes have undetermined (N) base calls: <br>
    ```
    zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz| awk 'NR % 4 == 2 && /N/'| wc -l

    and

    zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz| awk 'NR % 4 == 2 && /N/'| wc -l
    ```
    ```
    R2_001.fastq.gz : 3976613
    R3_001.fastq.gz : 3328051
    Total : 7304664
    ```
    
## Part 2
1. Define the problem

    When samples are run on Illumina they are assigned barcodes (indexes) as a way of identifying them. So biological read 1 (forward) and biological read 2 (reverse) are assigned dual matched barcodes index 1 and index 2 (which is the reverse complement of index 1). When sequencing is complete the biological reads and their indexes are outputted to four fastq files. Biological Reads 1 and 2 are labeled as R1 and R4 and indexes 1 and 2 are labelled as R2 and R3, respectively. The issue that arises at this point is that these reads have been multiplexed, they have been pooled together and sequenced simultaneously. What we want to do is demultiplex them--which means to sort the samples according to their assigned barcodes. During this sorting, we will need to account for low quality reads, unknown reads, and index hopped reads, with the ultiamte goal being to separate and assigned these reads, plus dual-matched, to appropriate files.

2. Describe output

Multiple fastq files will be outputted by our script. There will be two fastqs for each dual-matched indexes. In the case of dual-matched indexes, biological read 1 and 2 will have the same dual-matched barcode. The outputted files will look like R1 with index1-index2 and another file for R2 with index2-index1. All unknown reads in R1 will be assigned to an unknown R1 file; all index-hopped reads (mean indexes that are not appropriately matched with their counterpart) will be assigned to R1 index hopped, the same formatting goes for R2.

3. Upload your [4 input FASTQ files](../TEST-input_FASTQ) and your [>=6 expected output FASTQ files](../TEST-output_FASTQ).
4. [Pseudocode here!](./Pseudocode.md)

5. High level functions. For each function, be sure to include:
    1. Description/doc string
    2. Function headers (name and parameters)
    3. Test examples for individual functions
    4. Return statement

```python
def convert_phred(letter: str) -> int:
    '''Converts a single character into a phred score'''
    return QS

Input: I
Expected output: 40


def QS_check(QS_str: str, threshold: int = 30) -> bool:
    '''Takes average of QS and compares quality to defined threshold of 30. 
    If low quality returns False, if >=30 returns True. Uses convert_phred function'''
    return bool

Input: @AB
Expected output: True


def qual_score(phred_score: str) -> float:
    '''Calculates the quality score of string, phred_score, by using function convert_phred and computing the avg'''
    return avg QS

Input: #I
Expected output: 21


def reverse_complement(seq: str) -> str:
    '''This function takes in a string of DNA/RNA and returns the reverse complement. Will return A for U and N for N'''
    return rc_seq

Input: AGTCA
Expected output: TGACT
```