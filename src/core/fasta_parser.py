"""
Módulo responsável por ler e processar arquivos FASTA.
Implementa o padrão Single Responsibility Principle (SRP).
"""

# from typing import Dict, List

class FastaParserError(Exception):
    """Exceção personalizada para erros de parsing FASTA."""

class FastaParser:
    """
    Parser para arquivos no formato FASTA.

    Responsabilidade única: converter texto FASTA em estrutura de dados
    utilizável pelo sistema, validando o formato e conteúdo.
    """

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
            ValueError: Se o formato do FASTA for inválido.
        """

        if not fasta_text or not isinstance(fasta_text, str):
            raise ValueError("Entrada FASTA vazia ou conteúdo inválido.")

        lines = fasta_text.strip().split('\n')

        # Validação do formato FASTA
        if not lines or not lines[0].startswith('>'):
            raise FastaParserError("Formato FASTA inválido: cabeçalho malformatado ou ausente.")

        header = lines[0][1:].strip() # Remove '>'
        if not header:
            raise FastaParserError("Header FASTA está vazio após '>'")

        # Separa ID e descrição do cabeçalho
        header_parts = header.split(maxsplit=1)
        sequence_id = header.split()[0] if header else "ID desconhecido"
        description = header_parts[1] if len(header_parts) > 1 else ""

        # Concatena as linhas de sequência, removendo espaços
        # em branco e convertendo para maiúsculas
        sequence_lines = lines[1:]
        if not sequence_lines:
            raise FastaParserError("Formato FASTA inválido: sequência ausente.")
        sequence = ''.join(sequence_lines).upper().replace(' ', '').replace('\t', '')

        if not sequence:
            raise FastaParserError("Formato FASTA inválido: sequência vazia.")

        return {
            'header': header,
            'sequence': sequence,
            'id': sequence_id,
            'description': description
        }
