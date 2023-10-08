from modules.dna_rna import check_nucleotides, transcribe, reverse, complement, reverse_complement
from modules.amino_acids import aa_pattern_position, count_pattern_in_aa_sequences
from modules.fastq import calculate_gc, calculate_length, calculate_quality


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

