#!/usr/bin/env python

# Author: tnair@uoregon.edu

# Check out some Python module resources:
#   - https://docs.python.org/3/tutorial/modules.html
#   - https://python101.pythonlibrary.org/chapter36_creating_modules_and_packages.html
#   - and many more: https://www.google.com/search?q=how+to+write+a+python+module

'''This module is a collection of useful bioinformatics functions
written during the Bioinformatics and Genomics Program coursework.'''

__version__ = "0.5"         # Read way more about versioning here:
                            # https://en.wikipedia.org/wiki/Software_versioning

DNA_bases = set("ATCGN") #use set for membership checking, easier than using a list like you had in prev versions
RNA_bases = set("AUCGN")

complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N':'N', 'U':'A'}

def convert_phred(letter: str) -> int:
    '''Converts a single character into a phred score'''
    return ord(letter) - 33

def qual_score(phred_score: str) -> float:
    '''Calculates the quality score of string, phred_score, by using function convert_phred and computing the avg'''
    total = 0
    for letter in phred_score:
        total += convert_phred(letter)
    return total/len(phred_score)

def QS_check(QS_str: str, threshold: int = 30) -> bool:
    '''Takes average of QS and compares quality to defined threshold of 30. 
    If low quality returns False, if >=30 returns True. Uses convert_phred function'''
    return (int(qual_score(QS_str))) >= threshold

def validate_base_seq(seq: str, RNAflag: bool=False) -> bool:
    '''This function takes a string. Returns True if string is composed
    of only As, Ts (or Us if RNAflag), Gs, Cs. False otherwise. Case insensitive.'''
    seq = set(seq.upper())
    return seq <= (RNA_bases if RNAflag else DNA_bases) #is seq a subset of RNA_bases or DNA_bases 

def gc_content(DNA):
    '''Returns GC content of a DNA sequence as a decimal between 0 and 1.'''
    assert validate_base_seq(DNA), "String contains invalid characters - are you sure you used a DNA sequence?"
    
    DNA = DNA.upper()
    Gs = DNA.count("G")
    Cs = DNA.count("C")
    return (Gs+Cs)/len(DNA)

def calc_median(lst: list) -> float:
    '''Fxn calculates the median of a list and returns a float. If list length is 
    even will take avg of middle two values. If list len is odd will use floor division'''
    list_length = len(lst)
    if list_length % 2 == 0: #even
        x = lst[int(list_length/2)]
        y = lst[int(list_length/2)-1]
        return (x+y)/2
    else:
        return lst[list_length//2]


def oneline_fasta(input_fa: str, output_fa:str) -> None:
    '''Takes a multiline fasta file and returns single line file: where seq on one line only'''
    with open(input_fa, 'r') as in_fa, open(output_fa, 'w') as out_fa:
        for i, line in enumerate(in_fa):
            if line[0] == '>': #if first position in line is equal to >, True
                if i > 0: #if position of line is greater than 0
                    out_fa.write('\n') #write a new line at the end of header
                out_fa.write(line) #write the whole line to the output file
            else:
                out_fa.write(line.strip('\n')) #if not header line then strip line of newline and write

def reverse_complement(seq: str) -> str:
    '''This function takes in a string of DNA/RNA and returns the reverse complement. Will return A for U and N for N'''
    seq = seq.upper()
    rc_seq = ''.join(complement_dict[bp] for bp in seq)[::-1] #.join will cat chars together into str 
    return rc_seq


if __name__ == "__main__":
    # Write tests for functions above, Leslie has already populated some tests for convert_phred
    # These tests are run when you execute this file directly (instead of importing it)

    #Phred Score Conversion Test
    assert convert_phred("I") == 40, "wrong phred score for 'I'"
    assert convert_phred("C") == 34, "wrong phred score for 'C'"
    assert convert_phred("2") == 17, "wrong phred score for '2'"
    assert convert_phred("@") == 31, "wrong phred score for '@'"
    assert convert_phred("$") == 3, "wrong phred score for '$'"
    print("Your convert_phred function is working! Nice job")

    #Quality Score Conversion Test
    phred_score: str = "FFHHHHHJJJJIJIJJJIJJJJJJIIIJJJEHJJJJJJJIJIDGEHIJJFIGGGHFGHGFFF@EEDE@C??DDDDDDD@CDDDDBBDDDBDBDD@"
    assert qual_score("EEE") == 36
    assert qual_score("#I") == 21
    assert qual_score("EJ") == 38.5
    assert qual_score(phred_score) == 37.62105263157895, "wrong average phred score"
    print("You calculated the correct average phred score")

    #Quality Score Threshold Check
    
    assert QS_check("IIIIIIII") == True, "Wrong avg calculated" #40
    assert QS_check("@AB") == True, "Wrong avg calculated" #32
    assert QS_check(")*+,-") == False, "Wrong avg calculated" #10
    assert QS_check("??????") == True, "Wrong avg calculated" #30
    print("Function is properly calculating avg QS and comparing to threshold of 30")

    #Median Calculation Test
    assert calc_median([1,2,100]) == 2, "calc_median function does not work for odd length list"
    assert calc_median([1,2]) == 1.5, "calc_median function does not work for even length list"
    assert calc_median([1,1,1,1,1,1,1,1,1,5000]) == 1
    assert calc_median([1,2,3,4,5,6,7,13]) == 4.5
    print("Median calculated successfully")

    #Validate Seq is DNA/RNA
    assert validate_base_seq("AATAGAT") == True, "Validate base seq does not work on DNA"
    assert validate_base_seq("AAUAGAU", True) == True, "Validate base seq does not work on RNA"
    assert validate_base_seq("GGCANT") == True, "Validate base seq does not work with N in DNA"
    assert validate_base_seq("UAAGCN", True) == True, "Validate base seq does not work with N in RNA"
    assert validate_base_seq("Hello hello hello") == False, "Goodbye goodbye goodbye, validate base seq doesn't recognize nonDNA/RNA"
    assert validate_base_seq("Hi hi hi", True) == False, "Bye bye bye, validate base seq doesn't recognize nonDNA/RNA"
    print("Passed DNA, RNA, and non-nucleic tests")

    #Determine GC Content
    assert gc_content("GCGCGC") == 1
    assert gc_content("AATTATA") == 0
    assert gc_content("GCATCGAT") == 0.5
    print("Correctly calculated GC content")

    #Oneline Fasta test in separate script

    #Reverse Complement Test
    assert reverse_complement("AGTCA") == True, "Did not return rc seq"
    assert reverse_complement("UUUUUA") == True, "Did not return rc for RNA"
    assert reverse_complement("NNNTAC") == True, "Did not return rc for ambiguous bp, N"
    print("Correctly returned reverse complement sequences")
