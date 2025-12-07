"""
Módulo para transcrição de DNA para RNA.
Single Responsibility: apenas transcrição.
"""


class Transcription:
    """
    Realiza a transcrição de DNA para RNA.

    Single Responsibility: Responsável apenas pela transcrição,
    seguindo as regras biológicas (T -> U).
    """

    def transcribe(self, dna_sequence: str) -> str:
        """
        Transcreve DNA para RNA.

        Regra: Substitui todas as timinas (T) por uracilas (U).

        Args:
            dna_sequence: Sequência de DNA

        Returns:
            Sequência de RNA
        """
        # Validação básica
        if not dna_sequence:
            return ""

        # Transcrição: T -> U
        rna_sequence = dna_sequence.replace('T', 'U')

        return rna_sequence
