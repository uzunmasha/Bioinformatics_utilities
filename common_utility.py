from modules.dna_rna import check_nucleotides, transcribe, reverse, complement, reverse_complement
from modules.amino_acids import aa_pattern_position, count_pattern_in_aa_sequences
from modules.fastq import read_fastq_file, write_filtered_fastq, calculate_gc, calculate_length, calculate_quality
import os

def run_dna_rna_tools(*args):
    """
    Main function for nucleic acid sequences processing.
    Parameters: *args - nucleic acid sequences and operation.
    Returns: List of results according to called operations.
    """
    sequences = args[:-1]
    string = [value.strip().strip("'") for value in sequences]
    operation = args[-1]
    dna, rna = check_nucleotides(string)

    if operation == "transcribe":
        transcribe_result = transcribe(dna)
        if len(transcribe_result) > 1:
            return list(transcribe_result)
        else:
            return transcribe_result[0]

    if operation == "reverse":
        reverse_result = reverse(dna, rna)
        if len(reverse_result) > 1:
            return list(reverse_result)
        else:
            return reverse_result[0]

    if operation == "complement":
        complement_result = complement(dna, rna)
        if len(complement_result) > 1:
            return list(complement_result)
        else:
            return complement_result[0]

    if operation == "reverse_complement":
        reverse_complement_result = reverse_complement(dna, rna)
        if len(reverse_complement_result) > 1:
            return list(reverse_complement_result)
        else:
            return reverse_complement_result[0]


def run_aa_tools(*args):
    """
    Main function for amino acid sequences processing.
    Parameters: *args - amino acid sequences and operation.
    Returns: List of results according to called operations.
    """
    if len(args) < 3:
        return "Invalid input format: At least three arguments are required."

    amino_acids = [arg.upper() for arg in args[:-2]]
    pattern = args[-2].upper()
    operation = args[-1]

    if operation == "aa_pattern_position":
        result = aa_pattern_position(amino_acids, pattern)
    elif operation == "count_pattern_in_aa_sequences":
        result = count_pattern_in_aa_sequences(amino_acids, pattern)
    else:
        return "Invalid operation"

    return result


def run_fastq_filter(input_path: str, gc_bounds: tuple, length_bounds: tuple, quality_threshold: float, output_filename: str):
    """
    Main function for fastq files processing.
    Parameters: seqs - input dictionary with fastq files.
                gc_bounds - boundaries for filtering according to GC content. Default - (0,100).
                length_bounds - boundaries for filtering according to sequences length. Default - (0, 2**32).
                quality_threshold -boundaries for filtering according to sequences quality. Default - 0.
    Returns: fastq file with filtered data.
    """
    if not isinstance(gc_bounds, tuple):
        gc_bounds = (0, gc_bounds)

    if not isinstance(length_bounds, tuple):
        length_bounds = (0, length_bounds)

    seqs = read_fastq_file(input_path)
    filtered_seqs = {}

    gc_contents = calculate_gc(seqs)
    seq_lengths = calculate_length(seqs)
    quality_scores = calculate_quality(seqs)

    for key, value in seqs.items():
        gc_content = gc_contents.pop(0)
        seq_length = seq_lengths.pop(0)
        quality_score = quality_scores.pop(0)

        if gc_bounds[0] <= gc_content <= gc_bounds[1] and \
           length_bounds[0] <= seq_length <= length_bounds[1] and \
           quality_score >= quality_threshold:
            filtered_seqs[key] = value

    if output_filename is None:
        output_filename = f"filtered_{input_path}"
    elif not output_filename.endswith('.fastq'):
        output_filename += '.fastq'

    write_filtered_fastq(filtered_seqs, input_path, output_filename)

    return "Filtered data was saved into output file"
