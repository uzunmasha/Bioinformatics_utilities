def calculate_gc(seqs: dict) -> list:
    """
    Calculate the GC content for given sequences.
    As input takes a list of strings representing nucleic acid sequences.
    Input should be obtained from a FASTQ dictionary.
    Returns a list containing the GC content (%) for each sequence.

    """
    gc_contents = []
    for key, value in seqs.items():
        sequence_value = value[0]
        gc_count = (sequence_value.count('G') + sequence_value.count('C')) / len(sequence_value) * 100
        gc_contents.append(gc_count)
    return gc_contents


def calculate_length(seqs: dict) -> list:
    """
    Calculate length of given sequences.
    As input takes a list of strings representing nucleic acid sequences.
    Input should be obtained from a FASTQ dictionary.
    Returns a list with lengths for each sequence.

    """
    seq_lenghts = []
    for key, value in seqs.items():
        sequence_value = value[0]
        length_count = len(sequence_value)
        seq_lenghts.append(length_count)

    return seq_lenghts