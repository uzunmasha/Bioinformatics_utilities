# Bioinformatics_utilities
This repository contains homework scripts developed in the Python course at the [Bioinformatics Institute](https://bioinf.me/education/program) during the 2023-2024 academic year. Code with examples in the showcases.py notebook is marked with *.

### bioseq.py
This script provides various functionalities for working with biological data.
* `RNASequence/DNASequence/AminoAcidSequence` classes *

Assists in working with DNA, RNA, and amino acid sequencing data. 
* `filter_fastq` function

Filters FASTQ files based on GC content, sequence length, and quality threshold.
* `run_genscan` function *

Uses the [Genscan](http://hollywood.mit.edu/GENSCAN.html) prediction tool for DNA sequences, and extracts predicted peptide sequences, intron, and exon information.
* `telegram_logger` decorator

Sends messages and log files of run scripts to a Telegram chat for notification purposes. Implementation of this function was based on [Telegram bot API](https://core.telegram.org/bots/api).

### bio_files_processor.py
* `convert_multiline_fasta_to_oneline` function

Converts any number of DNA/RNA/protein sequences in FASTA file from multi-line FASTA files into one-line FASTA format.
* `OpenFasta` context manager *

Opens FASTA files, like the `open` built-in function. Returns separate FASTA records including ID, description, and sequence.

### custom_random_forest.py
* `RandomForestClassifierCustom` class *

Allows to apply parallelization for custom random forest class for faster usage

### test_modules.py
Contains tests to verify the functionality of the code in `bioseq.py` and `bio_files_processor.py`

### showcases.py
Demonstrates examples of using functions and classes from other files in the repository.
