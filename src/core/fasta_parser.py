"""
Módulo responsável por ler e processar arquivos FASTA.
Implementa o padrão Single Responsibility Principle (SRP).
"""

# from typing import Dict, List

class FastaParseError(Exception):
    """Exceção personalizada para erros de parsing FASTA."""

class FastaParser:
    """
    Parser para arquivos no formato FASTA.

    Responsabilidade única: converter texto FASTA em estrutura de dados
    utilizável pelo sistema, validando o formato e conteúdo.

    Attributes:
        VALID_DNA_BASES: Conjunto de bases válidas para DNA
        VALID_RNA_BASES: Conjunto de bases válidas para RNA
        VALID_PROTEIN_AA: Conjunto de aminoácidos válidos
    """

    # Constantes para validação
    VALID_DNA_BASES = set('ACGTN')
    VALID_RNA_BASES = set('ACGUN')
    VALID_PROTEIN_AA = set('')

    @staticmethod
    def parse(fasta_text: str) -> dict:
        """
        Converte uma string no formato FASTA em um dicionário com
        informações estruturadas.

        Args:
            fasta_text (str): Conteúdo do arquivo FASTA como string.

        Returns:
            dict: {
                'header': str,
                'sequence': str,
                'id': str,
                'description': str
            }

        Raises:
            FastaParseError: Se o formato do FASTA for inválido.
        """

        if not fasta_text or not isinstance(fasta_text, str):
            raise FastaParseError("Entrada FASTA vazia ou conteúdo inválido.")

        lines = fasta_text.strip().split('\n')

        # Validação do formato FASTA
        if not lines or not lines[0].startswith('>'):
            raise FastaParseError("Formato FASTA inválido: cabeçalho malformatado ou ausente.")

        header = lines[0][1:].strip() # Remove '>'
        if not header:
            raise FastaParseError("Header FASTA está vazio após '>'")

        # Separa ID e descrição do cabeçalho
        header_parts = header.split(maxsplit=1)
        sequence_id = header.split()[0] if header else "ID desconhecido"
        description = header_parts[1] if len(header_parts) > 1 else ""

        # Concatena as linhas de sequência, removendo espaços
        # em branco e convertendo para maiúsculas
        sequence_lines = lines[1:]
        if not sequence_lines:
            raise FastaParseError("Formato FASTA inválido: sequência ausente.")
        sequence = ''.join(sequence_lines).upper().replace(' ', '').replace('\t', '').replace('\n', '')

        return {
            'header': header,
            'sequence': sequence,
            'id': sequence_id,
            'description': description
        }

    @staticmethod
    def validate_sequence(sequence: str, sequence_type: str = 'dna') -> bool:
        """
        Valida se a sequência contém apenas caracteres válidos para o tipo especificado.

        Args:
            sequence: Sequência a ser validada
            sequence_type: Tipo de sequência ('dna', 'rna', ou 'protein')

        Returns:
            bool: True se válido, False caso contrário

        Raises:
            ValueError: Se sequence_type for inválido
        """

        sequence_type = sequence_type.lower()

        if sequence_type == 'dna':
            valid_chars = FastaParser.VALID_DNA_BASES
        elif sequence_type == 'rna':
            valid_chars = FastaParser.VALID_RNA_BASES
        elif sequence_type == 'protein':
            valid_chars = FastaParser.VALID_PROTEIN_AA
        else:
            raise ValueError(
                f"Tipo de sequência inválido: {sequence_type}. "
                "Use 'dna', 'rna', ou 'protein'"
            )

        sequence_upper = sequence.upper()
        return all(char in valid_chars for char in sequence_upper)
