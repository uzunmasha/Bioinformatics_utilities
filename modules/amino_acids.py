def aa_pattern_position(amino_acids: str, substring: str) -> list:
    """
    Searches for a substring of amino acids in the entire amino acid sequence.
    Takes a string of amino acids and a substring, which should be found.
    Returns the position where the searched one was found for the first time.
    Amino acids in the string should be indicated as one-letter symbols.

    """
    results = []
    for sequence in amino_acids:
        position = sequence.find(substring)
        results.append(position)
    return results


def count_pattern_in_aa_sequences(amino_acids: str, pattern: str) -> list:
    """
    Finds how many times a particular sequence(s) occurs in the original one.
    Takes a string of amino acids and a substring, which should be counted.
    Returns the count of searched amino acids.
    Amino acids in the string should be indicated as one-letter symbols.

    """
    results = []
    for sequence in amino_acids:
        aa_count = sequence.count(pattern)
        results.append(aa_count)
    return results
