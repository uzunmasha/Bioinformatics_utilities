from abc import ABC, abstractmethod
from dataclasses import dataclass
import os
import sys
import time
from io import StringIO

from Bio import SeqIO
from Bio.SeqUtils import gc_fraction 
from dotenv import load_dotenv
import pandas as pd
import requests
from bs4 import BeautifulSoup

load_dotenv()


class BiologicalSequence(ABC):
    """Abstract base class for biological sequences."""
    @abstractmethod
    def __init__(self, sequence: str):
        """Initialize a BiologicalSequence object with a given sequence."""
        self.sequence = sequence
    @abstractmethod
    def __len__(self) -> int:
        """Return the length of the sequence."""
        return len(self.sequence)
    @abstractmethod
    def __getitem__(self, index: int) -> str:
        """Return the item at the specified index."""
        return self.sequence[index]
    @abstractmethod
    def __str__(self) -> str:
        """Return the string representation of the sequence."""
        return self.sequence

    @abstractmethod
    def check_alphabet(self) -> bool:
        """Check if the sequence contains valid alphabet characters."""


class NucleicAcidSequence(BiologicalSequence):
    """Abstract base class"""

    def __init__(self, sequence: str):
        """Initialize a NucleicAcidSequence object with a given sequence."""
        super().__init__(sequence)

    def __len__(self) -> int:
        """Return the length of the sequence."""
        return len(self.sequence)

    def __getitem__(self, index: int) -> str:
        """Return the item at the specified index."""
        return self.sequence[index]

    def __str__(self) -> str:
        """Return the string representation of the sequence."""
        return self.sequence

    def check_alphabet(self) -> bool:
        """Check if the sequence contains valid nucleic acid alphabet characters."""
        raise NotImplementedError("Method 'check_alphabet' must be implemented in subclasses.")

    def complement(self) -> str:
        """Return the complement sequence."""
        comp_map_dna = {"A": "T", "G": "C", "T": "A", "C": "G", "a": "t", "t": "a", "g": "c", "c": "g"}
        return ''.join([comp_map_dna[base] for base in self.sequence])

    def gc_content(self) -> float:
        """Return the GC content of the sequence."""
        gc_count = (self.sequence.upper().count('G') + self.sequence.upper().count('C')) / len(self.sequence) * 100
        return gc_count


class DNASequence(NucleicAcidSequence):
    """Class representing a DNA sequence."""
    TRANSCRIBE_DICT = {
        'T': 'U',
        't': 'u'
    }
    alphabet = set("ATGCatgc")

    def __init__(self, sequence: str):
        """Initialize a DNASequence object with a given sequence."""
        super().__init__(sequence)
        if not self.check_alphabet():
            raise ValueError("Invalid DNA sequence")

    def check_alphabet(self) -> bool:
        """Check if the sequence contains valid DNA alphabet characters."""
        return set(self.sequence).issubset(self.alphabet)

    def transcribe(self) -> 'RNASequence':
        """Transcribe the DNA sequence into an RNA sequence."""
        transcribed_seq = ''.join(self.TRANSCRIBE_DICT[base] if base in self.TRANSCRIBE_DICT else base for base in self.sequence)
        return RNASequence(transcribed_seq)


class RNASequence(NucleicAcidSequence):
    """Class representing an RNA sequence."""
    alphabet = set("AUGCaugc")

    def __init__(self, sequence: str):
        """Initialize an RNASequence object with a given sequence."""
        super().__init__(sequence)
        if not self.check_alphabet():
            raise ValueError("Invalid RNA sequence")

    def check_alphabet(self) -> bool:
        """Check if the sequence contains valid RNA alphabet characters."""
        return set(self.sequence).issubset(self.alphabet)

    def reverse(self) -> 'RNASequence':
        """Return the reverse of the RNA sequence."""
        return type(self)(self.sequence[::-1])


class AminoAcidSequence(BiologicalSequence):
    """Class representing an amino acid sequence."""
    alphabet = set("ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy")

    def __init__(self, sequence: str):
        """Initialize an AminoAcidSequence object with a given sequence."""
        super().__init__(sequence)
        if not self.check_alphabet():
            raise ValueError("Invalid amino acid sequence")

    def __len__(self) -> int:
        """Return the length of the sequence."""
        return len(self.sequence)

    def __getitem__(self, index: int) -> str:
        """Return the item at the specified index."""
        return self.sequence[index]

    def __str__(self) -> str:
        """Return the string representation of the sequence."""
        return self.sequence

    def check_alphabet(self) -> bool:
        """Check if the sequence contains valid amino acid alphabet characters."""
        return set(self.sequence).issubset(self.alphabet)

    def amino_acid_profile(self):
        """Return the profile of the amino acid sequence."""
        self.check_alphabet()
    
        aa_biochemistry = {
            'hydrophobic': ['G', 'A', 'V', 'L', 'I', 'P', 'F', 'M', 'W'],
            'polar': ['S', 'T', 'C', 'N', 'Q', 'Y'],
            '- charged': ['E', 'D'],
            '+ charged': ['K', 'H', 'R']
        }
        profile = {}

        for group in aa_biochemistry:
            profile[group] = 0.0

        for amino_acid in self.sequence:
            for group_name, group_list in aa_biochemistry.items():
                if amino_acid.upper() in group_list:
                    profile[group_name] += 1

        total_length = len(self.sequence)
        for group, count in profile.items():
            profile[group] = round((count / total_length), 2)
        return profile


@dataclass
class GenscanOutput:
    """
    Dataclass represents the output of the run_genscan function.

    Attributes:
        status (int): The HTTP status code of the response.
        cds_list (list): A list of predicted peptide sequences from the analysed site.
        intron_list (list): A list of predicted intron information.
        exon_list (list): A list of predicted exon information.
    """
    status: int
    cds_list: str
    intron_list: str
    exon_list: str


def run_genscan(sequence: str = "",
                sequence_file: str = "",
                organism: str = "Vertebrate",
                exon_cutoff: float = 1.00,
                sequence_name: str = "") -> GenscanOutput:
    """
    Runs the Genscan prediction for the given DNA sequence or sequence file.

    Args:
        sequence (str): The DNA sequence.
        sequence_file (str): The path to the file containing the DNA sequence.
        organism (str, optional): The organism for which to perform the prediction. Defaults to "Vertebrate".
        exon_cutoff (float, optional): The exon probability cutoff. Defaults to 1.00.
        sequence_name (str, optional): The name of the sequence. Defaults to "".

    Returns:
        GenscanOutput: An object containing the prediction results from the site.
    """

    url = "http://argonaute.mit.edu/cgi-bin/genscanw_py.cgi"

    if sequence_file:
        with open(sequence_file, "r") as file:
            sequence = file.read().strip()

    payload = {
        "-o": organism,
        "-e": exon_cutoff,
        "-n": sequence_name,
        "-p": "Predicted peptides only",
        "-s": sequence,
    }

    response = requests.post(url, data=payload)
    output_html = response.text

    cds_list = []
    intron_list = []
    exon_list = []

    soup = BeautifulSoup(output_html, 'html.parser')
    result_data = soup.find('pre').string
    lines = result_data.split('\n')
    lines = list(filter(lambda x: x.strip(), lines))

    #  cds_list
    protein_indices = [i for i, line in enumerate(lines) if line.startswith('>')]
    for i in range(len(protein_indices)):
        start_index = protein_indices[i]
        end_index = protein_indices[i + 1] if i < len(protein_indices) - 1 else len(lines)
        protein_sequence = ''.join(lines[start_index+1:end_index])
        cds_list.append(protein_sequence)

    #  exon_list
    exon_table_start_index = lines.index('Gn.Ex Type S .Begin ...End .Len Fr Ph I/Ac Do/T CodRg P.... Tscr..')
    exon_table_end_index = None
    if exon_cutoff == 1.000:
        exon_table_end_index = lines.index('Suboptimal exons with probability > 1.000')
    else:
        exon_table_end_index = lines.index('Suboptimal exons with probability > 0.010')

    exon_table_data = lines[exon_table_start_index:exon_table_end_index]

    for line in exon_table_data[2:]:
        parts = line.split()
        exon_info = [parts[0], parts[1], int(parts[3]), int(parts[4])]
        exon_list.append(exon_info)

    exon_df = pd.DataFrame(exon_list, columns=['Index number', 'Type', 'Start', 'End'])

    #  intron_list
    #  Обрабатываем первый экзон отдельно
    for i in range(len(exon_df) - 1):
        current_exon = exon_df.iloc[i]
        next_exon = exon_df.iloc[i + 1]

        #  Вычисляем границы интрона для + и - цепей
        if current_exon['End'] < current_exon['Start'] and next_exon['End'] < next_exon['Start']:
            intron_start = current_exon['Start'] + 1
            intron_end = next_exon['End'] - 1
        elif current_exon['End'] < current_exon['Start'] and next_exon['End'] > next_exon['Start']:
            intron_start = current_exon['Start'] + 1
            intron_end = next_exon['Start'] - 1
        elif current_exon['End'] > current_exon['Start'] and next_exon['End'] < next_exon['Start']:
            intron_start = current_exon['End'] + 1
            intron_end = next_exon['End'] - 1
        else:
            intron_start = current_exon['End'] + 1
            intron_end = next_exon['Start'] - 1

        intron_list.append([current_exon['Index number'], 'Intron', intron_start, intron_end])

    return GenscanOutput(response.status_code, cds_list, intron_list, exon_list)


def send_telegram_message(chat_id: str, message: str, log_content: str, filename: str) -> None:
    """
    Sends a message along with a document to a Telegram chat.

    Args:
        chat_id (str): The ID of the Telegram chat.
        message (str): The message to be sent.
        log_content (str): The content of the log file.
        filename (str): The name of the file to be sent.

    Returns:
        None
    """
    if log_content.strip():
        bot_token = os.getenv("TG_API_TOKEN")
        files = {'document': (filename, log_content)}
        data = {'chat_id': chat_id, 'caption': message, 'parse_mode': 'Markdown'}
        requests.post(f'https://api.telegram.org/bot{bot_token}/sendDocument', files=files, data=data)


def format_time(seconds: float) -> str:
    """
    Formats the given time in seconds into a human-readable format.

    Args:
        seconds (float): The time in seconds.

    Returns:
        str: The formatted time string.
    """
    if seconds < 86400:
        formatted_time = time.strftime('%H:%M:%S', time.gmtime(seconds))
        milliseconds = int((seconds - int(seconds)) * 1000)
        formatted_time += f".{milliseconds:03}"
        return formatted_time
    else:
        days = seconds // 86400
        remaining_seconds = seconds % 86400
        return f'{days} days, {time.strftime("%H:%M:%S", time.gmtime(remaining_seconds))}'


def telegram_logger(chat_id: str):
    """
    Decorator function to log function execution time and exceptions to Telegram.

    Args:
        chat_id (str): The ID of the Telegram chat.

    Returns:
        callable: Decorator function.
    """
    def decorator(func: callable) -> callable:
        def wrapper(*args, **kwargs):
            start_time = time.time()
            output_buffer = StringIO()
            sys.stdout = output_buffer
            sys.stderr = output_buffer
            try:
                func(*args, **kwargs)
                execution_time = time.time() - start_time
                message = f"🥳 Function `{func.__name__}` successfully finished in `{format_time(execution_time)}`"
            except Exception as e:
                execution_time = time.time() - start_time
                message = f"😢 Function `{func.__name__}` failed with an exception: `{type(e).__name__}: {str(e)}`"
                raise
            finally:
                sys.stdout = sys.__stdout__
                sys.stderr = sys.__stderr__
                log_content = output_buffer.getvalue()
                filename = f"{func.__name__}.log"
                send_telegram_message(chat_id, message, log_content, filename)
        return wrapper
    return decorator


def filter_fastq(input_path: str, gc_lower_bound: float = 0, gc_upper_bound: float = 100,
                 length_lower_bound: float = 0, length_upper_bound: float = float('inf'),
                 quality_threshold: int = 0, output_filename: str = None) -> None:
    """
    Filters a FASTQ file based on GC content, sequence length, and quality threshold using Biopython.

    Args:
    - input_path (str): Path to the input FASTQ file.
    - gc_bounds (tuple): Tuple specifying the minimum and maximum GC content for filtering. Default is (0, 100).
    - length_bounds (tuple): Tuple specifying the minimum and maximum sequence length for filtering. Default is (0, infinity).
    - quality_threshold (float): Minimum quality score for filtering. Default is 0.
    - output_filename (str): Name of the output file. If None, the default filename will be used.

    Returns:
    - str: Message indicating the success of the filtering process.
    """
    filtered_seqs = []

    with open(input_path, 'r') as fastq_file:
        for record in SeqIO.parse(fastq_file, 'fastq'):
            gc_content = gc_fraction(record.seq)
            seq_length = len(record.seq)
            quality_score = sum(record.letter_annotations["phred_quality"]) / seq_length

            if (
                gc_bounds[0] <= gc_content <= gc_bounds[1] and
                length_bounds[0] <= seq_length <= length_bounds[1] and
                quality_score >= quality_threshold
            ):
                filtered_seqs.append(record)

    if output_filename is None:
        output_filename = f"filtered_{input_path}"
    elif not output_filename.endswith('.fastq'):
        output_filename += '.fastq'

    with open(output_filename, 'w') as output_file:
        SeqIO.write(filtered_seqs, output_file, 'fastq')

    return "Filtered data was saved into output file"