# SeqMate
This readme describes the user-friendly program SeqMate for performing various operations with fastq files, nucleic and amino acid sequences.

SeqMate can perform different operations:
* Transcribe, reverse, complement, reverse_complement DNA or RNA sequences
* Find for a particular amino acid(s) in the entire amino acid sequence(s)
* Calculate amino acid's occurrence in amino acid sequences
* Calculate GC content, sequences length, and quality for fastq sequences
* Filter fastq sequences according to counted values of GC content, length, and quality

## Usage
1. Clone this repo using SSH or HTTPS with 'modules' directory:
```bash
git clone https://github.com/uzunmasha/Bioinformatics_utilities.git
``` 
**or**
```bash
git clone git@github.com:uzunmasha/Bioinformatics_utilities.git
``` 
2. Import common_utility.py with the required function (listed below) in a code interpreter like Jupyter Notebook.
3. Enjoy your results!

## List of functions
For all functions, nucleic acids in the sequences should be indicated as one-letter symbols. For operations with DNA, RNA, and amino acids letters can be uppercase or lowercase.
For operations with fastq files, sequence letters should be uppercase.
### Operations with DNA and RNA 
All sequences in the input should be comma-separated. Any number of nucleic acid sequences is possible. The operation should be one and it should be pointed last.  
#### transcribe
This function takes DNA sequence(s) as input and converts it from DNA nucleotides to RNA nucleotides based on complementarity.

Usage:

```python
run_dna_rna_tools('ATG', 'transcribe') # 'AUG'
```
#### reverse
This function takes a DNA or RNA sequence as input and reverses it.

Usage:

```python
run_dna_rna_tools('ATG', 'aT', 'reverse') # ['GTA', 'Ta']
```
#### complement
This function takes a DNA or RNA sequence as input and generates its complementary sequence

Usage:

```python
run_dna_rna_tools('AtG', 'complement') # 'TaC'
```

#### reverse_complement
This function takes a DNA or RNA sequence as input, generates its complementary sequence, and then reverses it

Usage:
```python
run_dna_rna_tools('ATg', 'reverse_complement') #'cAT'
```
### Operations with amino acids
#### aa_pattern_position
This function searches for the presence of particular amino acid(s) in the entire amino acid sequence. As input, it takes a string of amino acids and a substring (pattern) that needs to be found. All sequences and pattern should be comma-separated. Any number of amino acid sequences is possible. The searched pattern should be one and it should be pointed last.  As an output, the function returns the position in the original sequence where the searched element was found for the first time.
Usage example:
```python
run_aa_tools('RNwDeACEQEZ', 'E','aa_pattern_position') #4
run_aa_tools('RNwDeACEQEZ', 'DFKAaaE','A','aa_pattern_position') #[5, 3]
```

#### count_pattern_in_aa_sequences
This function finds how many times a particular amino acid or sequence of several amino acids (pattern) occurs in the original sequence. As input, it takes a string of amino acids and a pattern that needs to be counted. All sequences and pattern should be comma-separated. Any number of amino acid sequences is possible. The searched substring should be one and it should be pointed last. As an output, the function returns the count of searched amino acid(s).
Usage example:
```python
run_aa_tools('GHcLfKF','f','count_pattern_in_aa_sequences') #2
run_aa_tools('HILAKMaF', 'GDaKFAAE','A','count_pattern_in_aa_sequences') #[2, 3]
```
### Fastq files filtering function
#### run_fastq_filter
This function calculates sequences GC content, length, and quality and filters fastq sequences according to counted values. As input, it takes four parameters:
* **seqs** - input dictionary with fastq files. An example can be found [here](https://github.com/Python-BI-2023/HW5_Modules/blob/main/example_data.py)
* **gc_bounds** - boundaries for filtering according to sequence GC content. Two numbers indicate the lower and upper limits (inclusively) of sequence filtering. If the user specifies only one value to filter, it will be considered the upper bound.
* **length_bounds** - boundaries for filtering according to sequence length. Default - (0, 2**32). Two numbers indicate the lower and upper limits (inclusively) of sequence filtering. If the user specifies only one value to filter, it will be considered the upper bound.
* **quality_threshold** - boundaries for filtering according to sequence quality. Default - 0.
As output, it returns a dictionary with filtered fastq files.

Usage:
```python
run_fastq_filter(seqs, gc_bounds = (0,100), length_bounds = (0, 2**32), quality_threshold = 0) # list of filtered fastq dictionary sequences
```
  
## Troubleshooting
* In function `'aa_pattern_position'` the position counting starts at 0, so don't be confused if the second element in the sequence has the output [1]. 
* In functions `'aa_pattern_position'` and `'count_pattern_in_aa_sequences'` output [-1] means that there is no such element in the sequence.
* In functions `'aa_pattern_position'` and `'count_pattern_in_aa_sequences'` the error message "name '..' is not defined" means that the given argument is not quoted in the input string.

**In case of non-working code:**

* Report any problems directly to the GitHub issue tracker

or

* Send your feedback to uzunmasha@gmail.com
