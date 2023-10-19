def read_fastq_file(input_path: str) -> dict:
    """
    Reads a FASTQ file and returns the data as a dictionary.
    """
    seqs = {}
    with open(input_path, mode='r') as fastq_file:
        lines = fastq_file.readlines()
        i = 0
        while i < len(lines):
            if lines[i].startswith("@"):
                name = lines[i].strip()
                sequence = lines[i + 1].strip()
                comment = lines[i + 2].strip()
                quality = lines[i + 3].strip()
                seqs[name] = (sequence, comment, quality)
                i += 4
    return seqs


def write_filtered_fastq(filtered_seqs: dict, input_path: str, output_filename: str):
    """
    Writes the filtered FASTQ data to a file in the 'fastq_filtered_results' folder
    """
    if output_filename is None:
        output_filename = os.path.splitext(os.path.basename(input_path))[0]

    output_directory = 'fastq_filtered_results'
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    output_path = os.path.join(output_directory, output_filename)

    with open(output_path, mode='w') as output_file:
        for key, value in filtered_seqs.items():
            output_file.write(f'{key}\n{value[0]}\n{value[1]}\n{value[2]}\n')


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
    quality_contents = []

    for key, value in seqs.items():
        quality_string = value[2]
        quality_sum = 0

        for sign in quality_string:
            quality_value = ord(sign) - 33
            quality_sum += quality_value

        quality_count = quality_sum / len(quality_string)
        quality_contents.append(quality_count)

    return quality_contents
