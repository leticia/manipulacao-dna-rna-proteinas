"""
Módulo responsável pela leitura de arquivos FASTA.
Implementa Strategy Pattern para diferentes formatos.
"""

class FastaReader:
    """
    Leitor de arquivos FASTA.

    Strategy Pattern: Permite extensão futura para outros formatos
    (FASTQ, GenBank, etc.) sem modificar o código existente.
    """

    def read(self, filepath: str) -> dict:
        """
        Lê um arquivo FASTA e retorna header e sequência.

        Args:
            filepath: Caminho para o arquivo FASTA

        Returns:
            Dicionário com 'header' e 'sequence'
        """
        try:
            with open(filepath, 'r', encoding='utf-8') as file:
                lines = file.readlines()

            # Interpretação do formato FASTA
            header = ""
            sequence_lines = []

            for line in lines:
                line = line.strip()

                if line.startswith('>'):
                    # Linha de cabeçalho
                    header = line[1:]  # Remove o '>'
                elif line:
                    # Linha de sequência (remove espaços e converte para maiúscula)
                    sequence_lines.append(line.upper().replace(' ', ''))

            sequence = ''.join(sequence_lines)

            return {
                'header': header,
                'sequence': sequence
            }

        except Exception as e:
            print(f"Erro ao ler arquivo: {e}")
            return None
