"""
Testes unitários para StatisticsCalculator.
Cobertura completa de ambas as estratégias e do contexto.
"""
import pytest
from src.core.statistics_calculator import (
    NucleotideStatistics,
    ProteinStatistics,
    StatisticsCalculator,
    StatisticsStrategy
)


class TestNucleotideStatistics:
    """Testes para estatísticas de nucleotídeos."""

    def test_calculate_simple_sequence(self):
        """Testa cálculo com sequência simples."""
        strategy = NucleotideStatistics()
        result = strategy.calculate("ATGC")

        assert result['length'] == 4
        assert result['counts'] == {'A': 1, 'T': 1, 'G': 1, 'C': 1}
        assert result['gc_content'] == 50.0
        assert result['type'] == 'nucleotide'

    def test_calculate_gc_content_100(self):
        """Testa cálculo de GC% = 100%."""
        strategy = NucleotideStatistics()
        result = strategy.calculate("GCGCGCGC")

        assert result['gc_content'] == 100.0
        assert result['at_content'] == 0.0

    def test_calculate_gc_content_0(self):
        """Testa cálculo de GC% = 0%."""
        strategy = NucleotideStatistics()
        result = strategy.calculate("ATATATAT")

        assert result['gc_content'] == 0.0
        assert result['at_content'] == 100.0

    def test_calculate_frequencies(self):
        """Testa cálculo de frequências relativas."""
        strategy = NucleotideStatistics()
        result = strategy.calculate("AAATGC")  # 3A, 1T, 1G, 1C

        assert result['frequencies']['A'] == 50.0
        assert result['frequencies']['T'] == 16.67
        assert result['frequencies']['G'] == 16.67
        assert result['frequencies']['C'] == 16.67

    def test_purine_pyrimidine_ratio(self):
        """Testa razão purina/pirimidina."""
        strategy = NucleotideStatistics()

        # AG (purinas) vs CT (pirimidinas)
        result = strategy.calculate("AGAGCT")  # 4 purinas, 2 pirimidinas
        assert result['purine_pyrimidine_ratio'] == 2.0

    def test_shannon_entropy_homogeneous(self):
        """Testa entropia de Shannon para sequência homogênea."""
        strategy = NucleotideStatistics()
        result = strategy.calculate("AAAA")

        # Entropia = 0 para sequência homogênea
        assert result['shannon_entropy'] == 0.0

    def test_shannon_entropy_maximum(self):
        """Testa entropia máxima (sequência equilibrada)."""
        strategy = NucleotideStatistics()
        result = strategy.calculate("ATGC")

        # Entropia máxima = 2 bits para 4 bases equiprováveis
        assert result['shannon_entropy'] == pytest.approx(2.0, abs=0.01)

    def test_empty_sequence(self):
        """Testa sequência vazia."""
        strategy = NucleotideStatistics()
        result = strategy.calculate("")

        assert result['length'] == 0
        assert result['gc_content'] == 0.0

    def test_rna_sequence(self):
        """Testa sequência de RNA (com U ao invés de T)."""
        strategy = NucleotideStatistics()
        result = strategy.calculate("AUGC")

        assert 'U' in result['counts']
        assert result['counts']['U'] == 1


class TestProteinStatistics:
    """Testes para estatísticas de proteínas."""

    def test_calculate_simple_protein(self):
        """Testa cálculo com proteína simples."""
        strategy = ProteinStatistics()
        result = strategy.calculate("ACDEFGHIKLMNPQRSTVWY")

        assert result['length'] == 20
        assert 'molecular_weight' in result
        assert result['type'] == 'protein'

    def test_molecular_weight_alanine(self):
        """Testa peso molecular de poli-alanina."""
        strategy = ProteinStatistics()
        result = strategy.calculate("AAA")  # 3 alaninas

        # 3 * 89.09 - 2 * 18.015 = 231.24 Da
        assert result['molecular_weight'] == pytest.approx(231.24, abs=0.1)

    def test_isoelectric_point_neutral(self):
        """Testa pI para proteína sem AAs ionizáveis."""
        strategy = ProteinStatistics()
        result = strategy.calculate("AAA")

        # Sem AAs ionizáveis na cadeia lateral
        assert result['isoelectric_point'] > 5.0

    def test_gravy_hydrophobic(self):
        """Testa GRAVY para proteína hidrofóbica."""
        strategy = ProteinStatistics()
        # I, V, L são muito hidrofóbicos
        result = strategy.calculate("IVLIVL")

        assert result['gravy_index'] > 2.0  # Positivo = hidrofóbico

    def test_gravy_hydrophilic(self):
        """Testa GRAVY para proteína hidrofílica."""
        strategy = ProteinStatistics()
        # D, E, K, R são hidrofílicos
        result = strategy.calculate("DEKRDE")

        assert result['gravy_index'] < -2.0  # Negativo = hidrofílico

    def test_composition_polar(self):
        """Testa composição de aminoácidos polares."""
        strategy = ProteinStatistics()
        # S, T, N, Q são polares
        result = strategy.calculate("STNQSTNQ")

        assert result['composition']['polar'] == 100.0

    def test_composition_nonpolar(self):
        """Testa composição de aminoácidos apolares."""
        strategy = ProteinStatistics()
        # A, V, I, L são apolares
        result = strategy.calculate("AVILGP")

        assert result['composition']['nonpolar'] > 90.0

    def test_composition_aromatic(self):
        """Testa composição de aminoácidos aromáticos."""
        strategy = ProteinStatistics()
        result = strategy.calculate("FWYFWY")  # F, W, Y são aromáticos

        assert result['composition']['aromatic'] == 100.0

    def test_remove_stop_codons(self):
        """Testa que códons de parada (*) são removidos."""
        strategy = ProteinStatistics()
        result = strategy.calculate("MFG*KWI*")

        assert '*' not in result['counts']
        assert result['length'] == 6  # 8 - 2 stops

    def test_empty_protein(self):
        """Testa proteína vazia."""
        strategy = ProteinStatistics()
        result = strategy.calculate("")

        assert result['length'] == 0
        assert result['molecular_weight'] == 0.0

class TestStatisticsCalculator:
    """Testes para o contexto StatisticsCalculator."""
    def test_default_strategy_nucleotide(self):
        """Testa que estratégia padrão é NucleotideStatistics."""
        calculator = StatisticsCalculator()

        strategy = calculator.get_strategy()
        assert isinstance(strategy, NucleotideStatistics)

    def test_change_strategy(self):
        """Testa mudança de estratégia em tempo de execução."""
        calculator = StatisticsCalculator(NucleotideStatistics())
        result1 = calculator.calculate_statistics("ATGC")

        assert 'gc_content' in result1

        # Mudar para proteína
        calculator.set_strategy(ProteinStatistics())
        result2 = calculator.calculate_statistics("ACDE")

        assert 'molecular_weight' in result2
        assert 'gc_content' not in result2

    def test_invalid_strategy_raises_error(self):
        """Testa que estratégia inválida lança erro."""
        calculator = StatisticsCalculator()

        with pytest.raises(TypeError, match="deve ser instância de StatisticsStrategy"):
            calculator.set_strategy("invalid")

    def test_none_sequence_raises_error(self):
        """Testa que sequência None lança erro."""
        calculator = StatisticsCalculator()

        with pytest.raises(ValueError, match="não pode ser None"):
            calculator.calculate_statistics(None)

    def test_non_string_sequence_raises_error(self):
        """Testa que sequência não-string lança erro."""
        calculator = StatisticsCalculator()

        with pytest.raises(TypeError, match="deve ser string"):
            calculator.calculate_statistics(12345)

    def test_auto_detect_dna(self):
        """Testa detecção automática de DNA."""
        calculator = StatisticsCalculator()
        result = calculator.auto_detect_and_calculate("ATGCGATCG")

        assert result['type'] == 'nucleotide'
        assert 'gc_content' in result

    def test_auto_detect_rna(self):
        """Testa detecção automática de RNA."""
        calculator = StatisticsCalculator()
        result = calculator.auto_detect_and_calculate("AUGCGAUCG")

        assert result['type'] == 'nucleotide'

    def test_auto_detect_protein(self):
        """Testa detecção automática de proteína."""
        calculator = StatisticsCalculator()
        result = calculator.auto_detect_and_calculate("MFKGDWIV")

        assert result['type'] == 'protein'
        assert 'molecular_weight' in result

@pytest.fixture
def dna_sequence():
    """Fixture com sequência DNA."""
    return "ATGCGATCGATCGATCGATCG"

@pytest.fixture
def rna_sequence():
    """Fixture com sequência RNA."""
    return "AUGCGAUCGAUCGAUCGAUCG"

@pytest.fixture
def protein_sequence():
    """Fixture com sequência proteína."""
    return "MFKGDWIVACDEFGHIKLMNPQRSTVWY"

class TestWithFixtures:
    """Testes usando fixtures."""
    def test_dna_statistics(self, dna_sequence):
        """Testa estatísticas de DNA."""
        calc = StatisticsCalculator(NucleotideStatistics())
        stats = calc.calculate_statistics(dna_sequence)

        assert stats['length'] == len(dna_sequence)
        assert 0 <= stats['gc_content'] <= 100

    def test_protein_statistics(self, protein_sequence):
        """Testa estatísticas de proteína."""
        calc = StatisticsCalculator(ProteinStatistics())
        stats = calc.calculate_statistics(protein_sequence)

        assert stats['length'] == len(protein_sequence)
        assert stats['molecular_weight'] > 0
