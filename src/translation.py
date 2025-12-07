"""
Módulo para tradução de RNA para proteína.
Single Responsibility: apenas tradução usando código genético.
"""


class Translation:
    """
    Realiza a tradução de RNA para proteína.

    Single Responsibility: Responsável apenas pela tradução,
    usando o código genético padrão.
    """

    # Código genético padrão (tabela de códons)
    CODON_TABLE = {
        'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
        'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
        'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*',
        'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W',
        'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
        'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',
        'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
        'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    }

    def translate(self, rna_sequence: str) -> str:
        """
        Traduz RNA para proteína usando o código genético.

        Args:
            rna_sequence: Sequência de RNA

        Returns:
            Sequência de aminoácidos (proteína)
        """
        if not rna_sequence:
            return ""

        protein = []

        # Processa a sequência em códons (grupos de 3 nucleotídeos)
        for i in range(0, len(rna_sequence) - 2, 3):
            codon = rna_sequence[i:i+3]

            # Traduz o códon
            amino_acid = self.CODON_TABLE.get(codon, 'X')  # 'X' para códons desconhecidos

            # Para no códon de parada
            if amino_acid == '*':
                break

            protein.append(amino_acid)

        return ''.join(protein)
