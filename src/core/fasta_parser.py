"""
Módulo responsável por ler e processar arquivos FASTA.
Implementa o padrão Single Responsibility Principle (SRP).
"""

from typing import Dict, List

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
    VALID_PROTEIN_AA = set('ACDEFGHIKLMNPQRSTVWY*X') # * = stop, X = desconhecido

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

    @staticmethod
    def detect_sequence_type(sequence: str) -> str:
        """
        Detecta o tipo de sequência (DNA, RNA, ou Proteína) com base nos caracteres presentes.

        Regras:
        - Se contém 'U' mas não 'T' → RNA
        - Se contém 'T' mas não 'U' → DNA
        - Se contém aminoácidos além de ACGTU → Proteína

        Args:
            sequence: Sequência a ser analisada

        Returns:
            str: 'dna', 'rna', 'protein', ou 'unknown'
        """

        sequence_upper = sequence.upper()
        unique_chars = set(sequence_upper)

        has_t = 'T' in sequence_upper
        has_u = 'U' in sequence_upper

        if has_t and not has_u:
            if FastaParser.validate_sequence(sequence, 'dna'):
                return 'dna'
        if has_u and not has_t:
            if FastaParser.validate_sequence(sequence, 'rna'):
                return 'rna'

        # Se não é DNA nem RNA, assume proteína
        nucleotide_chars = FastaParser.VALID_DNA_BASES | FastaParser.VALID_RNA_BASES
        has_protein_chars = any(char not in nucleotide_chars for char in unique_chars)

        if has_protein_chars:
            if FastaParser.validate_sequence(sequence, 'protein'):
                return 'protein'

        return 'unknown'

    @staticmethod
    def parse_multi_fasta(fasta_text: str) -> List[Dict[str, str]]:
        """
        Parse arquivo FASTA com múltiplas sequências.

        Args:
            fasta_text: Texto contendo múltiplas sequências FASTA

        Returns:
            list: Lista de dicionários, cada um representando uma sequência

        Raises:
            FastaParseError: Se formato for inválido
        """
        if not fasta_text or not fasta_text.strip():
            raise FastaParseError("Texto FASTA está vazio")

        sequences = []
        current_header = None
        current_sequence = []

        for line in fasta_text.strip().split('\n'):
            line = line.strip()

            if not line:
                continue

            if line.startswith('>'):
                # Salvar sequência anterior se existir
                if current_header is not None:
                    fasta_block = f"{current_header}\n{''.join(current_sequence)}"
                    sequences.append(FastaParser.parse(fasta_block))

                # Iniciar nova sequência
                current_header = line
                current_sequence = []
            else:
                current_sequence.append(line)

        # Salvar última sequência
        if current_header is not None:
            fasta_block = f"{current_header}\n{''.join(current_sequence)}"
            sequences.append(FastaParser.parse(fasta_block))

        if not sequences:
            raise FastaParseError("Nenhuma sequência válida encontrada")

        return sequences
