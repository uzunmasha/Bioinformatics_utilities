import inspect
import os
import re
import tempfile
from typing import List, Tuple

import pytest

from bio_files_processor import OpenFasta, FastaRecord, convert_multiline_fasta_to_oneline
from bioseq import DNASequence, RNASequence, AminoAcidSequence, run_genscan


class TestDNASequence:
    """
    Test methods of the DNASequence class.
    """
    def test_complement(self):
        """
        Test the complement method of the DNASequence class.
        """
        dna_seq = DNASequence("ATGC")
        assert dna_seq.complement() == "TACG"

    def test_gc_content(self):
        """
        Test the gc_content method of the DNASequence class.
        """
        dna = DNASequence("AGGC")
        gc_content = dna.gc_content()
        assert gc_content == 75.0


class TestRNASequence:
    """
    Test methods of the RNASequence class.
    """
    def test_check_alphabet_with_mixed_case(self):
        """
        Test the check_alphabet method of the RNASequence class with mixed case input.
        """
        rna_seq = RNASequence("AUGcu")
        assert rna_seq.check_alphabet() == True


class TestAminoAcidSequence:
    """
    Test methods of the AminoAcidSequence class.
    """
    def test_invalid_sequence(self):
        """
        Test raising ValueError for an invalid amino acid sequence.
        """
        invalid_sequence = "GAVLIz"

        with pytest.raises(ValueError):
            AminoAcidSequence(invalid_sequence)


class TestFastaProcessing:
    """
    Test functions related to processing FASTA files.
    """

    @pytest.fixture
    def multiline_fasta_content(self) -> str:
        """
        Fixture providing multiline FASTA content.
        """

        return ">Sequence1 Species1\nATGC\nCAATCG\nGAT\n>Sequence2 Species2\nTTAA\nCCGG\n>Sequence3 Species3\nGATTACA\n"

    @pytest.fixture
    def multiline_fasta_file(self, tmp_path: str, multiline_fasta_content: str) -> str:
        """
        Fixture creating a temporary multiline FASTA file.
        """
        filepath = tmp_path / "multiline.fasta"
        with open(filepath, "w") as f:
            f.write(multiline_fasta_content)
        return filepath

    def test_convert_multiline_fasta_to_oneline_output(self, multiline_fasta_file: str, tmp_path: str) -> None:
        """
        Test conversion of multiline FASTA to oneline format and check output file existence.
        """
        output_file = tmp_path / "one_line.fasta"
        convert_multiline_fasta_to_oneline(multiline_fasta_file, str(output_file))
        assert output_file.exists()

    def test_convert_multiline_fasta_to_oneline_content(self, multiline_fasta_file: str, tmp_path: str) -> None:
        """
        Test conversion of multiline FASTA to oneline format and check output content.
        """
        output_file = tmp_path / "one_line.fasta"
        convert_multiline_fasta_to_oneline(multiline_fasta_file, str(output_file))
        with open(output_file, "r") as f:
            output_content = f.read()
            assert output_content == ">Sequence1 Species1\nATGCCAATCGGAT\n>Sequence2 Species2\nTTAACCGG\n>Sequence3 Species3\nGATTACA\n"

    def test_open_fasta_read_records(self, tmp_path: str) -> None:
        """
        Test reading records from a FASTA file and comparing with expected records.
        """
        one_line_content = ">Sequence1 Species1\nATGCCAATCGGAT\n>Sequence2 Species2\nTTAACCGG\n>Sequence3 Species3\nGATTACA\n"
        one_line_file = tmp_path / "one_line.fasta"
        with open(one_line_file, "w") as f:
            f.write(one_line_content)

        expected_records = [
            FastaRecord("Sequence1", "Species1", "ATGCCAATCGGAT"),
            FastaRecord("Sequence2", "Species2", "TTAACCGG"),
            FastaRecord("Sequence3", "Species3", "GATTACA")
        ]
        with OpenFasta(one_line_file) as fasta_file:
            for expected_record in expected_records:
                record = fasta_file.read_record()
                assert record == expected_record


class TestGenscan:
    """
    Test the Genscan functionality.
    """
    def test_genscan_url(self) -> None:
        """
        Test that the correct URL for Genscan is present in the source code.
        """
        source_code = inspect.getsource(run_genscan)
        pattern = r'url\s*=\s*"http://argonaute\.mit\.edu/cgi-bin/genscanw_py\.cgi"'
        assert re.search(pattern, source_code) is not None
