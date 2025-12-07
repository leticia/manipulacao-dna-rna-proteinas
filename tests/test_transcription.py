"""
Testes unitários para Transcription.
"""

import pytest
from src.transcription import Transcription

class TestTranscription:
    """Testes para a classe Transcription."""

    def test_transcribe_basic(self):
        """Testa transcrição básica DNA -> RNA."""
        transcription = Transcription()
        dna = "ATGC"

        rna = transcription.transcribe(dna)

        assert rna == "AUGC"

    def test_transcribe_all_t_to_u(self):
        """Testa que todos os T são convertidos em U."""
        transcription = Transcription()
        dna = "TTTAAA"

        rna = transcription.transcribe(dna)

        assert rna == "UUUAAA"
        assert 'T' not in rna

    def test_transcribe_preserves_other_bases(self):
        """Testa que outras bases são preservadas."""
        transcription = Transcription()
        dna = "ACGACG"

        rna = transcription.transcribe(dna)

        assert rna == "ACGACG"

    def test_transcribe_empty_sequence(self):
        """Testa transcrição de sequência vazia."""
        transcription = Transcription()

        rna = transcription.transcribe("")

        assert rna == ""

    def test_transcribe_long_sequence(self):
        """Testa transcrição de sequência longa."""
        transcription = Transcription()
        dna = "ATGC" * 1000

        rna = transcription.transcribe(dna)

        assert len(rna) == 4000
        assert 'T' not in rna
        assert rna.count('U') == 1000

