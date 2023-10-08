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


def calculate_quality(seqs: dict) -> list:
    """
    Calculate quality of given sequences.
    As input takes a list of strings representing nucleic acid's quality symbols according to Phred33 scale.
    Input should be obtained from a FASTQ dictionary.
    Returns a list with mean values of quality for each sequence.
    """
    QUALITY_VOCAB = {'!': 0, '"': 1, '#': 2, '$': 3, '%': 4, '&': 5, "'": 6, '(': 7, ')': 8, '*': 9, '+': 10,
                     ',': 11, '-': 12, '.': 13, '/': 14, '0': 15, '1': 16, '2': 17, '3': 18, '4': 19, '5': 20,
                     '6': 21, '7': 22, '8': 23, '9': 24, ':': 25, ';': 26, '<': 27, '=': 28, '>': 29, '?': 30,
                     '@': 31, 'A': 32, 'B': 33, 'C': 34, 'D': 35, 'E': 36, 'F': 37, 'G': 38, 'H': 39, 'I': 40}

    quality_contents = []
    for key, value in seqs.items():
        quality_string = value[1]
        quality_sum = 0

        for sign in quality_string:
            quality_value = QUALITY_VOCAB.get(sign)
            quality_sum += int(quality_value)

        quality_count = quality_sum / len(quality_string)
        quality_contents.append(quality_count)

    return quality_contents
