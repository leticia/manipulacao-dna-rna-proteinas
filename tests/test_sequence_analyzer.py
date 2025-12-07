"""
Testes unitários para SequenceAnalyzer.
"""

import pytest
from src.sequence_analyzer import SequenceAnalyzer


class TestSequenceAnalyzer:
    """Testes para a classe SequenceAnalyzer."""

    def test_analyze_dna_sequence(self):
        """Testa análise de sequência de DNA."""
        analyzer = SequenceAnalyzer()
        sequence = "ATGCGATCG"

        result = analyzer.analyze(sequence)

        assert result['length'] == 9
        assert result['counts']['A'] == 2
        assert result['counts']['T'] == 2
        assert result['counts']['G'] == 3
        assert result['counts']['C'] == 2
        assert 'U' not in result['counts']

    def test_gc_percent_calculation(self):
        """Testa cálculo de GC%."""
        analyzer = SequenceAnalyzer()

        # 50% GC
        sequence = "ATGC"
        result = analyzer.analyze(sequence)
        assert result['gc_percent'] == 50.0

        # 100% GC
        sequence = "GGCC"
        result = analyzer.analyze(sequence)
        assert result['gc_percent'] == 100.0

        # 0% GC
        sequence = "AATT"
        result = analyzer.analyze(sequence)
        assert result['gc_percent'] == 0.0

    def test_analyze_rna_sequence(self):
        """Testa análise de sequência de RNA."""
        analyzer = SequenceAnalyzer()
        sequence = "AUGCGAUCG"

        result = analyzer.analyze(sequence)

        assert result['counts']['U'] == 2
        assert 'T' not in result['counts']

    def test_empty_sequence(self):
        """Testa análise de sequência vazia."""
        analyzer = SequenceAnalyzer()
        result = analyzer.analyze("")

        assert result['length'] == 0
        assert result['gc_percent'] == 0.0
        assert result['counts'] == {}
