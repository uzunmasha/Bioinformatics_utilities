# bio_files_processor
This readme describes the program `bio_files_processor` converts a multi-line FASTA file into a one-line FASTA format.

## Usage
1. Clone this repo using SSH or HTTPS with the modules directory:
```bash
git clone git@github.com:uzunmasha/Bioinformatics_utilities.git
``` 
**or**
```bash
git clone https://github.com/uzunmasha/Bioinformatics_utilities.git
``` 
2. Launch the program with the required function (listed below) in a code interpreter like Jupyter Notebook.
3. Enjoy your results!

## List of functions
### convert_multiline_fasta_to_oneline
This function calculates sequences GC content, length, and quality and filters fastq sequences according to counted values. As input, it takes four parameters:
* **input_fasta** - eads the input FASTA file, in which the sequence (DNA/RNA/protein/...) can be split into multiple lines. An example can be found [here](https://github.com/Python-BI-2023/HW6_Files/blob/main/example_data/example_multiline_fasta.fasta)
* **output_fasta** - saves sequence data to a new FASTA file in which each sequence fits into a single line.

Usage:
```python
convert_multiline_fasta_to_oneline(input_fasta = "example_multiline_fasta.fasta", output_fasta = 'myfasta') # fasta file of converted sequences
```
  
## Troubleshooting
* Argument `input_fasta` takes a full path to the input file (in quotes) if you are running the code from another directory. Or you can give only the file's name if you are in the same directory.
* Argument `output_fasta` should have a name in quotes or None.

**In case of non-working code:**

* Report any problems directly to the GitHub issue tracker

or

* Send your feedback to uzunmasha@gmail.com
