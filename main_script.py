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

    amino_acids = args[:-2]
    pattern = args[-2]
    operation = args[-1]

    if operation == "aa_pattern_position":
        result = aa_pattern_position(amino_acids, pattern)
    elif operation == "count_pattern_in_aa_sequences":
        result = count_pattern_in_aa_sequences(amino_acids, pattern)
    else:
        return "Invalid operation"

    return result


def run_fastq_filter(seqs, gc_bounds=(0,100), length_bounds=(0, 2**32), quality_threshold=0):
    """
    Main function for fastq files processing.
    Parameters: seqs - input dictionary with fastq files.
                gc_bounds - boundaries for filtering according to GC content. Default - (0,100).
                length_bounds - boundaries for filtering according to sequences length. Default - (0, 2**32).
                quality_threshold -boundaries for filtering according to sequences quality. Default - 0.
    Returns: Dictionary with filtered fastq files.
    """
    if not isinstance(gc_bounds, tuple):
        gc_bounds = (0, gc_bounds)

    if not isinstance(length_bounds, tuple):
        length_bounds = (0, length_bounds)

    gc_contents = calculate_gc(seqs)
    seq_lenghts = calculate_length(seqs)
    quality_contents = calculate_quality(seqs)

    filtered_seqs = {}

    for key, gc_content, seq_length, quality_content in zip(seqs.keys(), gc_contents, seq_lenghts, quality_contents):
        if gc_bounds[0] <= gc_content <= gc_bounds[1] and \
           length_bounds[0] <= seq_length <= length_bounds[1] and \
           quality_content > quality_threshold:
                filtered_seqs[key] = seqs[key]
    result = "\n".join([f'{key}: {value[0]}, {value[1]}' for key, value in filtered_seqs.items()])

    return result
