"""
Testes unitários para Translation.
"""

import pytest
from src.translation import Translation

class TestTranslation:
    """Testes para a classe Translation."""

    def test_translate_start_codon(self):
        """Testa tradução do códon de início (AUG -> M)."""
        translation = Translation()
        rna = "AUG"

        protein = translation.translate(rna)

        assert protein == "M"

    def test_translate_multiple_codons(self):
        """Testa tradução de múltiplos códons."""
        translation = Translation()
        # AUG (M) + GCU (A) + UGG (W)
        rna = "AUGGCUUGG"

        protein = translation.translate(rna)

        assert protein == "MAW"

    def test_translate_stops_at_stop_codon(self):
        """Testa que tradução para em códon de parada."""
        translation = Translation()
        # AUG (M) + UAA (stop) + GCU (A)
        rna = "AUGUAAGCU"

        protein = translation.translate(rna)

        # Deve parar no UAA e não traduzir GCU
        assert protein == "M"

    def test_translate_unknown_codon(self):
        """Testa tradução de códon desconhecido/inválido."""
        translation = Translation()
        rna = "NNN"  # Códon inválido

        protein = translation.translate(rna)

        assert protein == "X"

    def test_translate_empty_sequence(self):
        """Testa tradução de sequência vazia."""
        translation = Translation()

        protein = translation.translate("")

        assert protein == ""

    def test_translate_incomplete_codon(self):
        """Testa que códons incompletos não são traduzidos."""
        translation = Translation()
        # AUG (M) + GC (incompleto)
        rna = "AUGGC"

        protein = translation.translate(rna)

        # Apenas AUG deve ser traduzido
        assert protein == "M"
