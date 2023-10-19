import os
def convert_multiline_fasta_to_oneline(input_fasta: str, output_fasta: str) -> str:
    """
    This function reads a multi-line FASTA file and converts it into a one-line FASTA format.
    The input FASTA file is read, for each sequence name, its multi-line parts are merged into a single line.
    The resulting sequences are written to the output FASTA file in one-line format.

    Args:
        input_fasta: The path to the input multi-line FASTA file.
        output_fasta (optional): The path to the output one-line FASTA file. If not provided,
        a default filename based on the input filename will be generated. 
        To provide, a name should be given in parentheses.

    Returns a message indicating the status of the operation, including the name of the output file.
    """
        with open(input_fasta, mode='r') as fasta_file:
        sequence_name = ""
        sequence = ""
        results = []
        for line in fasta_file:
            line = line.strip()
            if line.startswith(">"):
                if sequence_name:
                    results.append((sequence_name, sequence))
                sequence_name = line
                sequence = ""
            else:
                sequence += line

        if sequence_name:
            results.append((sequence_name, sequence))

    if output_fasta is None:
        input_filename = os.path.splitext(os.path.basename(input_fasta))[0]
        output_fasta = f"{input_filename}_long_seq.fasta"
    elif not output_fasta.endswith('.fasta'):
        output_fasta += '.fasta'

    with open(output_fasta, mode='w') as output_file:
        for sequence_name, sequence in results:
            output_file.write(f'{sequence_name}\n{sequence}\n')
    return "Your data was saved into output file"
