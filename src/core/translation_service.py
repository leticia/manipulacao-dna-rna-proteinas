"""
Serviço de tradução de RNA para proteína.
Implementa padrão Singleton para tabela de códons.

Design Pattern: Singleton
- Garante única instância da tabela de códons em memória
- Economiza recursos (64 códons são constantes)
- Evita reinicializações desnecessárias
"""
from typing import List, Dict
from src.core.transcription_service import TranscriptionService

class CodonTable:
    """
    Tabela de códons genéticos (código genético padrão).

    Pattern: Singleton - garante única instância da tabela de códons.

    O código genético padrão (universal) possui 64 códons:
    - 61 códons que codificam aminoácidos
    - 3 códons de parada (UAA, UAG, UGA)
    - 1 códon de início (AUG) que também codifica Metionina
    """

    _instance = None
    _codon_table = None

    def __new__(cls):
        """
        Implementação do padrão Singleton.

        Garante que apenas uma instância da tabela exista em memória,
        economizando recursos e evitando recriações desnecessárias.

        Returns:
            CodonTable: Instância única da tabela de códons
        """
        if cls._instance is None:
            cls._instance = super().__new__(cls)
            cls._instance._initialize_table()
        return cls._instance

    def _initialize_table(self):
        """
        Inicializa tabela de códons (código genético padrão).

        Baseado no código genético universal (NCBI genetic code table 1).
        Cada códon (triplete de nucleotídeos) mapeia para um aminoácido
        representado por código de 1 letra.
        """
        self._codon_table = {
            # Fenilalanina (F)
            'UUU': 'F', 'UUC': 'F',

            # Leucina (L)
            'UUA': 'L', 'UUG': 'L', 'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',

            # Isoleucina (I)
            'AUU': 'I', 'AUC': 'I', 'AUA': 'I',

            # Metionina (M) - também códon de início
            'AUG': 'M',

            # Valina (V)
            'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',

            # Serina (S)
            'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S', 'AGU': 'S', 'AGC': 'S',

            # Prolina (P)
            'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',

            # Treonina (T)
            'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',

            # Alanina (A)
            'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',

            # Tirosina (Y)
            'UAU': 'Y', 'UAC': 'Y',

            # Códons de parada (stop codons)
            'UAA': '*', 'UAG': '*', 'UGA': '*',

            # Histidina (H)
            'CAU': 'H', 'CAC': 'H',

            # Glutamina (Q)
            'CAA': 'Q', 'CAG': 'Q',

            # Asparagina (N)
            'AAU': 'N', 'AAC': 'N',

            # Lisina (K)
            'AAA': 'K', 'AAG': 'K',

            # Ácido aspártico (D)
            'GAU': 'D', 'GAC': 'D',

            # Ácido glutâmico (E)
            'GAA': 'E', 'GAG': 'E',

            # Cisteína (C)
            'UGU': 'C', 'UGC': 'C',

            # Triptofano (W)
            'UGG': 'W',

            # Arginina (R)
            'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AGA': 'R', 'AGG': 'R',

            # Glicina (G)
            'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
        }

    def get_amino_acid(self, codon: str) -> str:
        """
        Retorna aminoácido correspondente a um códon.

        Args:
            codon: Códon de 3 nucleotídeos (RNA)

        Returns:
            str: Código de 1 letra do aminoácido
                 '*' para stop codon
                 'X' para códon desconhecido/inválido
        """
        return self._codon_table.get(codon.upper(), 'X')

    def is_start_codon(self, codon: str) -> bool:
        """
        Verifica se códon é de início (AUG).

        Args:
            codon: Códon de 3 nucleotídeos

        Returns:
            bool: True se for códon de início
        """
        return codon.upper() == 'AUG'

    def is_stop_codon(self, codon: str) -> bool:
        """
        Verifica se códon é de parada (UAA, UAG, UGA).

        Args:
            codon: Códon de 3 nucleotídeos

        Returns:
            bool: True se for códon de parada
        """
        return self.get_amino_acid(codon) == '*'

    def get_all_start_codons(self) -> List[str]:
        """
        Retorna todos os códons de início.

        No código genético padrão, apenas AUG é códon de início.

        Returns:
            list: Lista de códons de início
        """
        return ['AUG']

    def get_all_stop_codons(self) -> List[str]:
        """
        Retorna todos os códons de parada.

        Returns:
            list: Lista de códons de parada
        """
        return ['UAA', 'UAG', 'UGA']


class TranslationError(Exception):
    """Exceção customizada para erros de tradução."""
    pass


class TranslationService:
    """
    Serviço responsável pela tradução de RNA em proteína.

    Processo: RNA → Proteína (códons → aminoácidos)

    Funcionalidades:
    - Tradução direta de RNA
    - Busca de ORF (Open Reading Frame)
    - Tradução em múltiplas fases de leitura
    - Detecção de múltiplos ORFs
    """

    def __init__(self):
        """Inicializa serviço com tabela de códons (Singleton)."""
        self.codon_table = CodonTable()

    def translate(self, rna_sequence: str, find_orf: bool = True) -> str:
        """
        Traduz sequência de RNA em proteína.

        Args:
            rna_sequence: Sequência de RNA
            find_orf: Se True, busca por ORF (Open Reading Frame)

        Returns:
            str: Sequência de aminoácidos

        Raises:
            ValueError: Se sequência for inválida
        """
        if rna_sequence is None or not rna_sequence.strip():
            raise TranslationError("Não há sequência para traduzir")
        if not isinstance(rna_sequence, str):
            raise TranslationError(f"Sequência deve ser um texto")

        rna_sequence = rna_sequence.upper()

        # Validação
        valid_bases = set('AUGCN')
        if not all(base in valid_bases for base in rna_sequence):
            raise TranslationError("Sequência de RNA contém bases inválidas")

        if find_orf:
            return self._translate_with_orf(rna_sequence)
        else:
            return self._translate_direct(rna_sequence)

    def _translate_direct(self, rna_sequence: str) -> str:
        """
        Tradução direta sem busca de ORF.

        Traduz desde o início da sequência até o final ou códon de parada.
        """
        protein = []

        # Itera em passos de 3 (códons)
        for i in range(0, len(rna_sequence) - 2, 3):
            codon = rna_sequence[i:i+3]

            if len(codon) != 3:
                break

            amino_acid = self.codon_table.get_amino_acid(codon)

            # Para se encontrar códon de parada
            if amino_acid == '*':
                break

            protein.append(amino_acid)

        return ''.join(protein) if protein else "Nenhuma tradução possível"

    def _translate_with_orf(self, rna_sequence: str) -> str:
        """
        Tradução com busca de ORF (Open Reading Frame).

        Procura por códon de início (AUG) e traduz até códon de parada.
        """
        protein = []
        start_found = False

        # Busca códon de início
        for i in range(0, len(rna_sequence) - 2):
            codon = rna_sequence[i:i+3]

            if len(codon) != 3:
                break

            # Encontrou códon de início
            if not start_found and self.codon_table.is_start_codon(codon):
                start_found = True
                protein.append('M')  # Metionina

                # Continua traduzindo a partir daqui
                for j in range(i + 3, len(rna_sequence) - 2, 3):
                    next_codon = rna_sequence[j:j+3]

                    if len(next_codon) != 3:
                        break

                    amino_acid = self.codon_table.get_amino_acid(next_codon)

                    # Para se encontrar códon de parada
                    if amino_acid == '*':
                        break

                    protein.append(amino_acid)

                break  # ORF encontrado, não precisa continuar

        if not protein:
            return "Nenhum ORF encontrado"

        return ''.join(protein)

    def find_all_orfs(self, rna_sequence: str, min_length: int = 30) -> List[Dict[str, any]]:
        """
        Encontra todos os ORFs possíveis na sequência.

        Args:
            rna_sequence: Sequência de RNA
            min_length: Tamanho mínimo do ORF (aminoácidos)

        Returns:
            list: Lista de dicionários com informações de cada ORF:
                {
                    'frame': int,          # Fase de leitura (0, 1, 2)
                    'start': int,          # Posição inicial (nucleotídeo)
                    'end': int,            # Posição final (nucleotídeo)
                    'length': int,         # Comprimento (aminoácidos)
                    'protein': str         # Sequência de aminoácidos
                }
        """
        orfs = []
        rna_sequence = rna_sequence.upper()

        # Verifica as 3 fases de leitura
        for frame in range(3):
            i = frame
            while i < len(rna_sequence) - 2:
                codon = rna_sequence[i:i+3]

                if len(codon) != 3:
                    break

                # Encontrou códon de início
                if self.codon_table.is_start_codon(codon):
                    start_pos = i
                    protein = ['M']

                    # Traduz até encontrar códon de parada
                    for j in range(i + 3, len(rna_sequence) - 2, 3):
                        next_codon = rna_sequence[j:j+3]

                        if len(next_codon) != 3:
                            break

                        amino_acid = self.codon_table.get_amino_acid(next_codon)

                        if amino_acid == '*':
                            end_pos = j + 3
                            protein_seq = ''.join(protein)

                            # Adiciona se atender tamanho mínimo
                            if len(protein_seq) >= min_length:
                                orfs.append({
                                    'frame': frame,
                                    'start': start_pos,
                                    'end': end_pos,
                                    'length': len(protein_seq),
                                    'protein': protein_seq
                                })

                            i = j + 3  # Pula para depois do códon de parada
                            break

                        protein.append(amino_acid)
                    else:
                        # Chegou ao fim sem encontrar stop codon
                        i += 3
                else:
                    i += 3

        # Ordenar por tamanho (maior primeiro)
        orfs.sort(key=lambda x: x['length'], reverse=True)

        return orfs

    def translate_six_frames(self, dna_sequence: str) -> Dict[str, List[str]]:
        """
        Traduz sequência de DNA nas 6 fases de leitura possíveis.

        3 frames na fita direta (forward)
        3 frames na fita reverso-complementar (reverse)

        Args:
            dna_sequence: Sequência de DNA

        Returns:
            dict: {
                'forward': [frame0, frame1, frame2],
                'reverse': [frame0, frame1, frame2]
            }
        """

        # Transcrever DNA → RNA
        rna_forward = TranscriptionService.transcribe(dna_sequence)

        # Obter reverso-complemento e transcrever
        dna_reverse = TranscriptionService.get_reverse_complement(
            dna_sequence, 'dna'
        )
        rna_reverse = TranscriptionService.transcribe(dna_reverse)

        # Traduzir 3 frames da fita direta
        forward_frames = []
        for frame in range(3):
            protein = self._translate_direct(rna_forward[frame:])
            forward_frames.append(protein)

        # Traduzir 3 frames da fita reversa
        reverse_frames = []
        for frame in range(3):
            protein = self._translate_direct(rna_reverse[frame:])
            reverse_frames.append(protein)

        return {
            'forward': forward_frames,
            'reverse': reverse_frames
        }

    def get_codon_usage(self, rna_sequence: str) -> Dict[str, int]:
        """
        Calcula uso de códons na sequência.

        Args:
            rna_sequence: Sequência de RNA

        Returns:
            dict: Dicionário com contagem de cada códon

        Examples:
            >>> service = TranslationService()
            >>> usage = service.get_codon_usage("AUGAAAUAA")
            >>> usage['AUG']
            1
        """
        rna_sequence = rna_sequence.upper()
        codon_counts = {}

        for i in range(0, len(rna_sequence) - 2, 3):
            codon = rna_sequence[i:i+3]

            if len(codon) == 3:
                codon_counts[codon] = codon_counts.get(codon, 0) + 1

        return codon_counts
