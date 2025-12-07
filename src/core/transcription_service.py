"""
Serviço de transcrição de DNA para RNA.
Implementa Single Responsibility Principle (SRP).

Processo biológico:
    DNA (template) → RNA (transcrito)
    T (Timina) → U (Uracila)
    Outras bases permanecem iguais (A, G, C)
"""
from typing import Optional


class TranscriptionError(Exception):
    """Exceção customizada para erros de transcrição."""
    pass


class TranscriptionService:
    """
    Serviço responsável pela transcrição de DNA em RNA.

    Responsabilidade única: converter sequências de DNA para RNA
    seguindo as regras da biologia molecular.

    Regra principal:
        T (Timina) → U (Uracila)
        A, G, C permanecem iguais
    """

    # Bases válidas para validação
    VALID_DNA_BASES = set('ATGCN')
    VALID_RNA_BASES = set('AUGCN')

    @staticmethod
    def transcribe(dna_sequence: str, validate: bool = True) -> str:
        """
        Transcreve sequência de DNA para RNA.

        Regra biológica: T → U

        Args:
            dna_sequence: Sequência de DNA
            validate: Se True, valida bases antes de transcrever

        Returns:
            str: Sequência de RNA correspondente

        Raises:
            TranscriptionError: Se validação falhar
            ValueError: Se entrada for inválida
        """
        # Validação de entrada
        if dna_sequence is None:
            raise ValueError("Sequência de DNA não pode ser None")

        if not isinstance(dna_sequence, str):
            raise ValueError(
                f"Sequência deve ser string, recebido {type(dna_sequence).__name__}"
            )

        if not dna_sequence.strip():
            raise TranscriptionError("Sequência de DNA está vazia")

        # Converter para uppercase
        dna_sequence = dna_sequence.upper()

        # Validação de bases (se solicitado)
        if validate:
            # Remove espaços para validação
            clean_seq = dna_sequence.replace(' ', '').replace('\n', '')
            invalid_bases = set(clean_seq) - TranscriptionService.VALID_DNA_BASES

            if invalid_bases:
                raise TranscriptionError(
                    f"Sequência de DNA contém bases inválidas: {invalid_bases}. "
                    f"Bases válidas: {TranscriptionService.VALID_DNA_BASES}"
                )

        # Transcrição: T → U
        rna_sequence = dna_sequence.replace('T', 'U')

        return rna_sequence

    @staticmethod
    def reverse_transcribe(rna_sequence: str, validate: bool = True) -> str:
        """
        Realiza transcrição reversa (RNA → DNA).

        Útil para análise de sequências de RNA ou vírus de RNA.

        Regra: U → T

        Args:
            rna_sequence: Sequência de RNA
            validate: Se True, valida bases antes de transcrever

        Returns:
            str: Sequência de DNA correspondente

        Raises:
            TranscriptionError: Se validação falhar

        Examples:
            >>> TranscriptionService.reverse_transcribe("AUGC")
            'ATGC'
        """
        if rna_sequence is None:
            raise ValueError("Sequência de RNA não pode ser None")

        if not isinstance(rna_sequence, str):
            raise ValueError(
                f"Sequência deve ser string, recebido {type(rna_sequence).__name__}"
            )

        if not rna_sequence.strip():
            raise TranscriptionError("Sequência de RNA está vazia")

        rna_sequence = rna_sequence.upper()

        # Validação
        if validate:
            clean_seq = rna_sequence.replace(' ', '').replace('\n', '')
            invalid_bases = set(clean_seq) - TranscriptionService.VALID_RNA_BASES

            if invalid_bases:
                raise TranscriptionError(
                    f"Sequência de RNA contém bases inválidas: {invalid_bases}"
                )

        # Transcrição reversa: U → T
        dna_sequence = rna_sequence.replace('U', 'T')

        return dna_sequence

    @staticmethod
    def get_complement(sequence: str, sequence_type: str = 'dna') -> str:
        """
        Obtém sequência complementar.

        DNA complementar:
            A ↔ T
            G ↔ C

        RNA complementar:
            A ↔ U
            G ↔ C

        Args:
            sequence: Sequência de DNA ou RNA
            sequence_type: 'dna' ou 'rna'

        Returns:
            str: Sequência complementar
        """
        sequence = sequence.upper()
        sequence_type = sequence_type.lower()

        if sequence_type == 'dna':
            complement_map = {
                'A': 'T', 'T': 'A',
                'G': 'C', 'C': 'G',
                'N': 'N'
            }
        elif sequence_type == 'rna':
            complement_map = {
                'A': 'U', 'U': 'A',
                'G': 'C', 'C': 'G',
                'N': 'N'
            }
        else:
            raise ValueError(f"Tipo de sequência inválido: {sequence_type}")

        complement = ''.join(complement_map.get(base, base) for base in sequence)
        return complement

    @staticmethod
    def get_reverse_complement(sequence: str, sequence_type: str = 'dna') -> str:
        """
        Obtém sequência reverso-complementar.
        Útil para análise de fita oposta do DNA.

        Args:
            sequence: Sequência de DNA ou RNA
            sequence_type: 'dna' ou 'rna'

        Returns:
            str: Sequência reverso-complementar
        """
        complement = TranscriptionService.get_complement(sequence, sequence_type)
        return complement[::-1]  # Reverse
