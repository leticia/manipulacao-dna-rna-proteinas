"""
Módulo para cálculo de estatísticas bioinformáticas.
Implementa Strategy Pattern para diferentes tipos de análise.

Design Pattern: Strategy
- Permite trocar algoritmo de cálculo em tempo de execução
- Facilita adição de novas estratégias sem modificar código existente
- Cada estratégia é testável independentemente
"""
from abc import ABC, abstractmethod
from typing import Dict, Optional
from src.core.fasta_parser import FastaParser
import math

class StatisticsStrategy(ABC):
    """
    Interface base para estratégias de cálculo estatístico.
    Pattern: Strategy - define contrato que todas as estratégias devem seguir.
    """

    @abstractmethod
    def calculate(self, sequence: str) -> Dict:
        """
        Calcula estatísticas para uma sequência.

        Args:
            sequence: Sequência de DNA/RNA/Proteína

        Returns:
            dict: Dicionário com estatísticas calculadas
        """
        pass

class NucleotideStatistics(StatisticsStrategy):
    """
    Estratégia de cálculo para sequências de nucleotídeos (DNA/RNA).

    Estatísticas calculadas:
    - Tamanho da sequência
    - Contagem absoluta de cada base
    - Frequências relativas (%)
    - Conteúdo GC%
    - Razão purina/pirimidina
    - Conteúdo AT%
    - Entropia de Shannon (bits)
    """

    # Classificação de bases
    PURINES = {'A', 'G'}      # Adenina, Guanina
    PYRIMIDINES = {'C', 'T', 'U'}  # Citosina, Timina, Uracila

    def calculate(self, sequence: str) -> Dict:
        """
        Calcula estatísticas para sequências de nucleotídeos.

        Args:
            sequence: Sequência de nucleotídeos (DNA ou RNA)

        Returns:
            dict: Estatísticas calculadas
        """
        sequence = sequence.upper()
        length = len(sequence)

        if length == 0:
            return self._empty_stats()

        # Contagem de cada nucleotídeo
        counts = {}
        for base in sequence:
            counts[base] = counts.get(base, 0) + 1

        # Cálculo de frequências relativas
        frequencies = {
            base: round((count / length) * 100, 2)
            for base, count in counts.items()
        }

        # Cálculo de GC%
        g_count = counts.get('G', 0)
        c_count = counts.get('C', 0)
        gc_content = round(((g_count + c_count) / length) * 100, 2)

        # Cálculo de AT% (ou AU% para RNA)
        a_count = counts.get('A', 0)
        t_count = counts.get('T', 0)
        u_count = counts.get('U', 0)
        at_content = round(((a_count + t_count + u_count) / length) * 100, 2)

        # Razão purina/pirimidina
        purine_count = sum(counts.get(base, 0) for base in self.PURINES)
        pyrimidine_count = sum(counts.get(base, 0) for base in self.PYRIMIDINES)

        if pyrimidine_count > 0:
            purine_pyrimidine_ratio = round(purine_count / pyrimidine_count, 2)
        else:
            purine_pyrimidine_ratio = float('inf') if purine_count > 0 else 0.0

        # Entropia de Shannon (medida de complexidade)
        shannon_entropy = self._calculate_shannon_entropy(frequencies)

        return {
            'length': length,
            'counts': counts,
            'frequencies': frequencies,
            'gc_content': gc_content,
            'at_content': at_content,
            'purine_count': purine_count,
            'pyrimidine_count': pyrimidine_count,
            'purine_pyrimidine_ratio': purine_pyrimidine_ratio,
            'shannon_entropy': shannon_entropy,
            'type': 'nucleotide'
        }

    def _calculate_shannon_entropy(self, frequencies: Dict[str, float]) -> float:
        """
        Calcula entropia de Shannon (bits).

        Fórmula: H = -Σ(p_i * log2(p_i))
        Onde p_i é a probabilidade (frequência) de cada símbolo.

        Entropia indica complexidade:
        - 0 bits: sequência homogênea (ex: "AAAA")
        - ~2 bits: máxima complexidade para DNA

        Args:
            frequencies: Dicionário de frequências relativas

        Returns:
            float: Entropia em bits
        """
        entropy = 0.0
        for freq_percent in frequencies.values():
            if freq_percent > 0:
                # Converter % para probabilidade (0-1)
                p = freq_percent / 100.0
                entropy -= p * math.log2(p)

        return round(entropy, 3)

    def _empty_stats(self) -> Dict:
        """Retorna estatísticas vazias para sequência vazia."""
        return {
            'length': 0,
            'counts': {},
            'frequencies': {},
            'gc_content': 0.0,
            'at_content': 0.0,
            'purine_count': 0,
            'pyrimidine_count': 0,
            'purine_pyrimidine_ratio': 0.0,
            'shannon_entropy': 0.0,
            'type': 'nucleotide'
        }


class ProteinStatistics(StatisticsStrategy):
    """
    Estratégia de cálculo para sequências de proteínas.

    Estatísticas calculadas:
    - Tamanho (número de aminoácidos)
    - Contagem de cada aminoácido
    - Frequências relativas
    - Peso molecular (Da)
    - Ponto isoelétrico (pI) aproximado
    - Índice de instabilidade
    - Índice GRAVY (hidrofobicidade)
    - Composição por classe de aminoácido
    """

    # Pesos moleculares médios (Da) - resíduo em cadeia peptídica
    AA_WEIGHTS = {
        'A': 89.09,  'C': 121.15, 'D': 133.10, 'E': 147.13,
        'F': 165.19, 'G': 75.07,  'H': 155.16, 'I': 131.17,
        'K': 146.19, 'L': 131.17, 'M': 149.21, 'N': 132.12,
        'P': 115.13, 'Q': 146.15, 'R': 174.20, 'S': 105.09,
        'T': 119.12, 'V': 117.15, 'W': 204.23, 'Y': 181.19
    }

    # Valores pK para cálculo de pI
    PK_VALUES = {
        'N_TERM': 9.69,  # N-terminal
        'C_TERM': 2.34,  # C-terminal
        'C': 8.33,   # Cisteína
        'D': 3.86,   # Ácido aspártico
        'E': 4.25,   # Ácido glutâmico
        'H': 6.00,   # Histidina
        'K': 10.53,  # Lisina
        'R': 12.48,  # Arginina
        'Y': 10.07   # Tirosina
    }

    # Índice de hidrofobicidade (Kyte-Doolittle)
    HYDROPHOBICITY = {
        'A': 1.8,  'C': 2.5,  'D': -3.5, 'E': -3.5,
        'F': 2.8,  'G': -0.4, 'H': -3.2, 'I': 4.5,
        'K': -3.9, 'L': 3.8,  'M': 1.9,  'N': -3.5,
        'P': -1.6, 'Q': -3.5, 'R': -4.5, 'S': -0.8,
        'T': -0.7, 'V': 4.2,  'W': -0.9, 'Y': -1.3
    }

    # Classificações de aminoácidos
    POLAR = {'S', 'T', 'N', 'Q', 'C', 'Y'}
    NONPOLAR = {'A', 'V', 'I', 'L', 'M', 'F', 'W', 'P', 'G'}
    ACIDIC = {'D', 'E'}
    BASIC = {'K', 'R', 'H'}
    AROMATIC = {'F', 'W', 'Y'}
    ALIPHATIC = {'A', 'V', 'I', 'L', 'M'}

    def calculate(self, sequence: str) -> Dict:
        """
        Calcula estatísticas para sequências de aminoácidos.

        Args:
            sequence: Sequência de aminoácidos (código de 1 letra)

        Returns:
            dict: Estatísticas calculadas

        Examples:
            >>> strategy = ProteinStatistics()
            >>> stats = strategy.calculate("MFKGDWIV")
            >>> stats['length']
            8
            >>> 'molecular_weight' in stats
            True
        """
        sequence = sequence.upper().replace('*', '')  # Remove stops
        length = len(sequence)

        if length == 0:
            return self._empty_stats()

        # Contagem de aminoácidos
        counts = {}
        for aa in sequence:
            counts[aa] = counts.get(aa, 0) + 1

        # Frequências relativas
        frequencies = {
            aa: round((count / length) * 100, 2)
            for aa, count in counts.items()
        }

        # Peso molecular
        molecular_weight = self._calculate_molecular_weight(sequence)

        # Ponto isoelétrico
        isoelectric_point = self._calculate_pi(sequence)

        # Índice GRAVY
        gravy_index = self._calculate_gravy(sequence)

        # Composição por classe
        composition = self._calculate_composition(sequence)

        return {
            'length': length,
            'counts': counts,
            'frequencies': frequencies,
            'molecular_weight': molecular_weight,
            'isoelectric_point': isoelectric_point,
            'gravy_index': gravy_index,
            'composition': composition,
            'type': 'protein'
        }

    def _calculate_molecular_weight(self, sequence: str) -> float:
        """
        Calcula peso molecular aproximado (Da).

        Peso = Soma(peso_aa) - (n-1) * peso_água
        onde n = número de aminoácidos
        """
        weight = sum(
            self.AA_WEIGHTS.get(aa, 0)
            for aa in sequence
        )

        # Subtrair peso de moléculas de água das ligações peptídicas
        if len(sequence) > 1:
            weight -= (len(sequence) - 1) * 18.015  # H2O

        return round(weight, 2)

    def _calculate_pi(self, sequence: str) -> float:
        """
        Calcula ponto isoelétrico aproximado.

        Método simplificado baseado em pK values.
        """
        # Contar aminoácidos ionizáveis
        counts = {}
        for aa in sequence:
            if aa in self.PK_VALUES:
                counts[aa] = counts.get(aa, 0) + 1

        # Cálculo simplificado: média ponderada dos pKs
        if not counts:
            return 7.0  # Neutro se não houver AAs ionizáveis

        total_pk = 0.0
        total_count = 0

        # N-terminal e C-terminal
        total_pk += self.PK_VALUES['N_TERM']
        total_pk += self.PK_VALUES['C_TERM']
        total_count += 2

        # Cadeias laterais
        for aa, count in counts.items():
            total_pk += self.PK_VALUES[aa] * count
            total_count += count

        pi = total_pk / total_count
        return round(pi, 2)

    def _calculate_gravy(self, sequence: str) -> float:
        """
        Calcula índice GRAVY (Grand Average of Hydropathy).

        GRAVY = Soma(hidrofobicidade) / n

        Interpretação:
        - GRAVY > 0: Proteína hidrofóbica
        - GRAVY < 0: Proteína hidrofílica
        """
        if not sequence:
            return 0.0

        total_hydro = sum(
            self.HYDROPHOBICITY.get(aa, 0)
            for aa in sequence
        )

        gravy = total_hydro / len(sequence)
        return round(gravy, 3)

    def _calculate_composition(self, sequence: str) -> Dict:
        """
        Calcula composição por classe de aminoácido.

        Returns:
            dict: Percentual de cada classe
        """
        length = len(sequence)

        polar = sum(1 for aa in sequence if aa in self.POLAR)
        nonpolar = sum(1 for aa in sequence if aa in self.NONPOLAR)
        acidic = sum(1 for aa in sequence if aa in self.ACIDIC)
        basic = sum(1 for aa in sequence if aa in self.BASIC)
        aromatic = sum(1 for aa in sequence if aa in self.AROMATIC)
        aliphatic = sum(1 for aa in sequence if aa in self.ALIPHATIC)

        return {
            'polar': round((polar / length) * 100, 2),
            'nonpolar': round((nonpolar / length) * 100, 2),
            'acidic': round((acidic / length) * 100, 2),
            'basic': round((basic / length) * 100, 2),
            'aromatic': round((aromatic / length) * 100, 2),
            'aliphatic': round((aliphatic / length) * 100, 2)
        }

    def _empty_stats(self) -> Dict:
        """Retorna estatísticas vazias."""
        return {
            'length': 0,
            'counts': {},
            'frequencies': {},
            'molecular_weight': 0.0,
            'isoelectric_point': 7.0,
            'gravy_index': 0.0,
            'composition': {},
            'type': 'protein'
        }


class StatisticsCalculator:
    """
    Context class que utiliza Strategy pattern.

    Permite trocar a estratégia de cálculo conforme o tipo de sequência,
    sem necessidade de modificar código cliente.

    Examples:
        >>> # Para DNA
        >>> calc = StatisticsCalculator(NucleotideStatistics())
        >>> stats = calc.calculate_statistics("ATGC")
        >>> 
        >>> # Mudar para proteína
        >>> calc.set_strategy(ProteinStatistics())
        >>> stats = calc.calculate_statistics("MFKG")
    """

    def __init__(self, strategy: Optional[StatisticsStrategy] = None):
        """
        Inicializa calculadora com uma estratégia específica.

        Args:
            strategy: Estratégia de cálculo. Se None, usa NucleotideStatistics
        """
        self._strategy = strategy or NucleotideStatistics()

    def set_strategy(self, strategy: StatisticsStrategy) -> None:
        """
        Altera a estratégia de cálculo em tempo de execução.

        Args:
            strategy: Nova estratégia a ser utilizada
        """
        if not isinstance(strategy, StatisticsStrategy):
            raise TypeError(
                f"Estratégia deve ser instância de StatisticsStrategy, "
                f"recebido {type(strategy).__name__}"
            )
        self._strategy = strategy

    def get_strategy(self) -> StatisticsStrategy:
        """Retorna estratégia atual."""
        return self._strategy

    def calculate_statistics(self, sequence: str) -> Dict:
        """
        Calcula estatísticas usando a estratégia configurada.

        Args:
            sequence: Sequência a ser analisada

        Returns:
            dict: Estatísticas calculadas pela estratégia atual

        Raises:
            ValueError: Se sequência for None ou vazia
        """
        if sequence is None:
            raise ValueError("Sequência não pode ser None")

        if not isinstance(sequence, str):
            raise TypeError(
                f"Sequência deve ser string, recebido {type(sequence).__name__}"
            )

        return self._strategy.calculate(sequence)

    def auto_detect_and_calculate(self, sequence: str) -> Dict:
        """
        Detecta automaticamente o tipo de sequência e calcula estatísticas.

        Args:
            sequence: Sequência a ser analisada

        Returns:
            dict: Estatísticas calculadas com estratégia apropriada
        """
        seq_type = FastaParser.detect_sequence_type(sequence)

        if seq_type in ['dna', 'rna']:
            self.set_strategy(NucleotideStatistics())
        elif seq_type == 'protein':
            self.set_strategy(ProteinStatistics())
        else:
            # Default para nucleotídeo
            self.set_strategy(NucleotideStatistics())

        return self.calculate_statistics(sequence)
