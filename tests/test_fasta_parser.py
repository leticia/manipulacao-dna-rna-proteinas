import pytest
from src.core.fasta_parser import FastaParser

class TestFastaParser:
    """
    Testes unitários para a classe FastaParser.
    """

    def test_parse_valid_fasta(self):
        """Testa parsing de FASTA válido."""
        fasta = ">gene1 description\nATGCGATCG\nATCGATCG"
        result = FastaParser.parse(fasta)

        assert result['header'] == "gene1 description"
        assert result['sequence'] == "ATGCGATCGATCGATCG"
        assert result['id'] == "gene1"
        assert result['description'] == "description"

    def test_parse_multiline_sequence(self):
        """Testa parsing com sequência em múltiplas linhas."""
        fasta = ">gene1\nATGC\nGATC\nTAGC"
        result = FastaParser.parse(fasta)

        assert result['sequence'] == "ATGCGATCTAGC"
        assert len(result['sequence']) == 12

    def test_parse_lowercase_converted_to_uppercase(self):
        """Testa que sequência minúscula é convertida para maiúscula."""
        fasta = ">gene1\natgcgatcg"
        result = FastaParser.parse(fasta)

        assert result['sequence'] == "ATGCGATCG"
        assert result['sequence'].isupper()

    def test_parse_with_spaces_in_sequence(self):
        """Testa que espaços na sequência são removidos."""
        fasta = ">gene1\nATGC GATC TAGC"
        result = FastaParser.parse(fasta)

        assert result['sequence'] == "ATGCGATCTAGC"
        assert ' ' not in result['sequence']

    def test_parse_header_without_description(self):
        """Testa parsing de header sem descrição."""
        fasta = ">gene1\nATGC"
        result = FastaParser.parse(fasta)

        assert result['id'] == "gene1"
        assert result['description'] == ""
