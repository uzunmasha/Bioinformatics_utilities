import os

from dataclasses import dataclass


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


@dataclass
class FastaRecord:
    """
    Represents a single record in a FASTA file.

    Attributes:
        id (str): The identifier of the record.
        description (str): The description of the record.
        sequence (str): The sequence data of the record.
    """
    id: str
    description: str
    sequence: str

    def __repr__(self) -> str:
        """
        Return a string representation of the FastaRecord object.

        Returns:
            str: A string representation of the FastaRecord object.
        """
        truncated_seq = self.sequence[:100] + "..." if len(self.sequence) > 100 else self.sequence
        return f"{self.id} {self.description}\n{truncated_seq}\n"


class OpenFasta:
    """
    A context manager for reading records from a FASTA file.

    Attributes:
        filename (str): The path to the FASTA file.
    """
    def __init__(self, filename: str, mode: str = 'r'):
        self.filename = filename
        self.file = None
        self.name = None
        self.seq = []

    def __enter__(self):
        if self.filename:
            self.file = open(self.filename, 'r')
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        if self.file:
            self.file.close()

    def __iter__(self):
        return self

    def __next__(self):
        for line in self.file:
            line = line.strip()
            if line.startswith(">"):
                if self.name:
                    record_id = self.name[1:].strip()
                    parts = record_id.split(' ', 1)
                    record_id = parts[0]
                    description = parts[1] if len(parts) > 1 else ''
                    sequence = ''.join(self.seq)
                    self.name = line
                    self.seq = []
                    return FastaRecord(record_id, description, sequence)
                else:
                    self.name = line
            else:
                self.seq.append(line)
        if self.name:
            record_id = self.name[1:].strip()
            parts = record_id.split(' ', 1)
            record_id = parts[0]
            description = parts[1] if len(parts) > 1 else ''
            sequence = ''.join(self.seq)
            self.name = None
            self.seq = []
            return FastaRecord(record_id, description, sequence)
        raise StopIteration

    def read_record(self):
        record = next(self)
        return record

    def read_records(self):
        return list(self)
