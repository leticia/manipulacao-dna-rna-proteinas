
"""
Módulo para análise estatística de sequências.
Single Responsibility Principle: apenas cálculos estatísticos.
"""


class SequenceAnalyzer:
    """
    Analisador de estatísticas de sequências de DNA/RNA.

    Single Responsibility: Responsável apenas por cálculos estatísticos,
    sem conhecer detalhes de I/O ou formato de dados.
    """

    def analyze(self, sequence: str) -> dict:
        """
        Calcula estatísticas da sequência.

        Args:
            sequence: Sequência de DNA ou RNA

        Returns:
            Dicionário com estatísticas:
            - length: tamanho da sequência
            - counts: contagem de cada nucleotídeo
            - gc_percent: percentual de GC
        """
        # Tamanho da sequência
        length = len(sequence)

        # Contagem de nucleotídeos
        counts = self._count_nucleotides(sequence)

        # Cálculo de GC%
        gc_percent = self._calculate_gc_percent(counts, length)

        return {
            'length': length,
            'counts': counts,
            'gc_percent': gc_percent
        }

    def _count_nucleotides(self, sequence: str) -> dict:
        """
        Conta a ocorrência de cada nucleotídeo.

        Args:
            sequence: Sequência de DNA/RNA

        Returns:
            Dicionário com contagens
        """
        # Nucleotídeos possíveis
        nucleotides = ['A', 'C', 'G', 'T', 'U']

        counts = {nt: sequence.count(nt) for nt in nucleotides}

        # Remove contagens zero
        return {k: v for k, v in counts.items() if v > 0}

    def _calculate_gc_percent(self, counts: dict, length: int) -> float:
        """
        Calcula o percentual de GC.

        Args:
            counts: Dicionário com contagens de nucleotídeos
            length: Tamanho total da sequência

        Returns:
            Percentual de GC (0-100)
        """
        if length == 0:
            return 0.0

        gc_count = counts.get('G', 0) + counts.get('C', 0)
        return (gc_count / length) * 100
