# Assignment the First

## Part 1
1. Be sure to upload your Python script. Provide a link to it here:

| File name | label | Read length | Phred encoding |
|---|---|---|---|
| 1294_S1_L008_R1_001.fastq.gz | Biological Read 1 | 101 | 33 |
| 1294_S1_L008_R2_001.fastq.gz | Index 1 | 8 | 33 |
| 1294_S1_L008_R3_001.fastq.gz | Index 2 | 8 | 33 |
| 1294_S1_L008_R4_001.fastq.gz | Biological Read 2 | 101 | 33 |

2. Per-base NT distribution
    1. Use markdown to insert your 4 histograms here.
    2. For index reads a quality cutoff of 30 is preferred because it equates to 99.9% accuracy and we require high confidence since we are demultiplexing indexes and assigning them to reads. A single error may cause reads to be incorrectly assigned or tossed which we want to avoid. For biological read pairs, a cutoff of Q20 _is_ the minimum acceptable threshold, and generally acceptable, as it equates to a error rate of 1%. However, a Q30 is preferred if we want a higher sensitivity for analyses. Higher sensitivity would be beneficial for situations like detecting low-frequeny variants, mutation detection, limiting false positives, or variant calling.
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
2. Describe output
3. Upload your [4 input FASTQ files](../TEST-input_FASTQ) and your [>=6 expected output FASTQ files](../TEST-output_FASTQ).
4. Pseudocode
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