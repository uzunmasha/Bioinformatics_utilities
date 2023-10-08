def check_nucleotides(string):
    """
    Cheks if the input values are DNA, RNA or something else.
    Takes a string of amino acids as input.
    Appends sequences to dna or rna variables.

    """
    dna_chars = set('CAGTcagt')
    rna_chars = set('CAGUcagu')

    is_dna = True
    is_rna = True

    for sequence in string:
        for char in sequence:
            if char not in dna_chars:
                is_dna = False
                break
        if not is_dna:
            break

    for sequence in string:
        for char in sequence:
            if char not in rna_chars:
                is_rna = False
                break
        if not is_rna:
            break

    dna = None
    rna = None

    if is_dna and not is_rna:
        dna = string
    elif is_rna and not is_dna:
        rna = string
    elif is_rna and is_dna:
        additional_input = input("Is this DNA or RNA?: ")
        if additional_input == 'DNA':
            dna = string
        elif additional_input == 'RNA':
            rna = string
        else:
            print("Input 'DNA' or 'RNA'")
    else:
        print("Neither DNA nor RNA")
    return dna, rna


def transcribe(dna):
    """
    Takes a DNA sequence as input.
    Converts it from nucleotides to RNA nucleotides based on complementarity.

    """
    transcription_map = {"A": "A", "T": "U", "G": "G", "C": "C", "a": "a", "t": "u", "g": "g", "c": "c"}
    transcr_rna = []
    for sequence in dna:
        transcr_sequence = ""
        for nucleotide in sequence:
            if nucleotide in transcription_map:
                transcr_sequence += transcription_map[nucleotide]
            else:
                transcr_sequence += nucleotide
        transcr_rna.append(transcr_sequence)
    return transcr_rna


def reverse(dna, rna):
    """
    This function takes a DNA or RNA sequence as input and reverses it.

    """
    if dna is not None:
        rev_dna = []
        for sequence in dna:
            reversed_sequence = sequence[::-1]
            rev_dna.append(reversed_sequence)
        return rev_dna
    elif rna is not None:
        rev_rna = []
        for sequence in rna:
            reversed_sequence = sequence[::-1]
            rev_rna.append(reversed_sequence)
        return rev_rna
    else:
        return "Empty string"


def complement(dna, rna):
    """
    This function takes a DNA or RNA sequence as input and generates its complementary sequence.

    """
    comp_map_dna = {"A": "T", "G": "C", "T": "A", "C": "G", "a": "t", "t": "a", "g": "c", "c": "g"}
    comp_map_rna = {"A": "U", "G": "C", "U": "A", "C": "G", "a": "u", "u": "a", "g": "c", "c": "g"}
    compl_seq = []

    if dna is not None:
        for sequence in dna:
            complemented_sequence = ''
            for nucl in sequence:
                if nucl in comp_map_dna:
                    complemented_sequence += comp_map_dna[nucl]
                else:
                    complemented_sequence += nucl
            compl_seq.append(complemented_sequence)
        return compl_seq
    elif rna is not None:
        for sequence in rna:
            complemented_sequence = ''
            for nucl in sequence:
                if nucl in comp_map_rna:
                    complemented_sequence += comp_map_rna[nucl]
                else:
                    complemented_sequence += nucl
            compl_seq.append(complemented_sequence)
        return compl_seq
    else:
        return "Empty string"


def reverse_complement(dna, rna):
    """
    This function takes a DNA or RNA sequence as input, generates its complementary sequence, and then reverses it.

    """
    comp_map_dna = {"A": "T", "G": "C", "T": "A", "C": "G", "a": "t", "t": "a", "g": "c", "c": "g"}
    comp_map_rna = {"A": "U", "G": "C", "U": "A", "C": "G", "a": "u", "u": "a", "g": "c", "c": "g"}
    compl_seq_rev = []

    if dna is not None:
        for sequence in dna:
            complemented_sequence = ''
            for nucl in sequence:
                if nucl in comp_map_dna:
                    complemented_sequence += comp_map_dna[nucl]
                else:
                    complemented_sequence += nucl
            compl_seq_rev.append(complemented_sequence[::-1])
        return compl_seq_rev
    elif rna is not None:
        for sequence in rna:
            complemented_sequence = ''
            for nucl in sequence:
                if nucl in comp_map_dna:
                    complemented_sequence += comp_map_rna[nucl]
                else:
                    complemented_sequence += nucl
            compl_seq_rev.append(complemented_sequence[::-1])
        return compl_seq_rev
    else:
        return "Empty string"
