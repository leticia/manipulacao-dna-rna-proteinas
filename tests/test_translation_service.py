"""
Testes para TranslationService.
"""
import pytest
from src.core.translation_service import TranslationService, CodonTable


class TestCodonTable:
    """Testes para tabela de códons."""

    def test_singleton(self):
        """Testa que CodonTable é singleton."""
        table1 = CodonTable()
        table2 = CodonTable()
        assert table1 is table2

    def test_start_codon(self):
        """Testa identificação de códon de início."""
        table = CodonTable()
        assert table.is_start_codon("AUG")
        assert not table.is_start_codon("UUU")

    def test_stop_codons(self):
        """Testa identificação de códons de parada."""
        table = CodonTable()
        assert table.is_stop_codon("UAA")
        assert table.is_stop_codon("UAG")
        assert table.is_stop_codon("UGA")
        assert not table.is_stop_codon("AUG")

    def test_get_amino_acid(self):
        """Testa obtenção de aminoácidos."""
        table = CodonTable()
        assert table.get_amino_acid("AUG") == "M"
        assert table.get_amino_acid("UUU") == "F"
        assert table.get_amino_acid("GGG") == "G"


class TestTranslationService:
    """Testes de tradução."""

    def test_translate_simple(self):
        """Testa tradução simples."""
        service = TranslationService()
        result = service.translate("AUGUUUUAA", find_orf=False)
        assert result == "MF"
