"""
Testes unitários para FastaReader.
"""

import pytest
from pathlib import Path
from src.fasta_reader import FastaReader


class TestFastaReader:
    """Testes para a classe FastaReader."""

    def test_read_valid_fasta(self, tmp_path):
        """Testa leitura de arquivo FASTA válido."""
        # Cria arquivo FASTA temporário
        fasta_file = tmp_path / "test.fasta"
        fasta_content = """>test_sequence
ATGCGATCGATCG
ATCGATCGATCGA"""
        fasta_file.write_text(fasta_content)

        # Lê o arquivo
        reader = FastaReader()
        result = reader.read(str(fasta_file))

        # Validações
        assert result is not None
        assert result['header'] == 'test_sequence'
        assert result['sequence'] == 'ATGCGATCGATCGATCGATCGATCGA'

    def test_read_multiple_lines(self, tmp_path):
        """Testa leitura de sequência em múltiplas linhas."""
        fasta_file = tmp_path / "multi.fasta"
        fasta_content = """>multi_line
ATGC
GATC
GATA"""
        fasta_file.write_text(fasta_content)

        reader = FastaReader()
        result = reader.read(str(fasta_file))

        assert result['sequence'] == 'ATGCGATCGATA'

    def test_read_with_spaces(self, tmp_path):
        """Testa que espaços são removidos."""
        fasta_file = tmp_path / "spaces.fasta"
        fasta_content = """>with_spaces
ATG C GAT C"""
        fasta_file.write_text(fasta_content)

        reader = FastaReader()
        result = reader.read(str(fasta_file))

        assert result['sequence'] == 'ATGCGATC'

    def test_read_lowercase_converted(self, tmp_path):
        """Testa que minúsculas são convertidas."""
        fasta_file = tmp_path / "lower.fasta"
        fasta_content = """>lowercase
atgcgatc"""
        fasta_file.write_text(fasta_content)

        reader = FastaReader()
        result = reader.read(str(fasta_file))

        assert result['sequence'] == 'ATGCGATC'

    def test_read_nonexistent_file(self):
        """Testa leitura de arquivo inexistente."""
        reader = FastaReader()
        result = reader.read("nonexistent.fasta")

        assert result is None

