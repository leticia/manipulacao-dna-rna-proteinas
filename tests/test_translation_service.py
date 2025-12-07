"""
Testes para TranslationService.
"""
import pytest
from src.core.translation_service import TranslationService, CodonTable, TranslationError


class TestCodonTable:
    """Testes para tabela de códons."""

    def test_singleton(self):
        """Testa que CodonTable é singleton."""
        table1 = CodonTable()
        table2 = CodonTable()
        assert table1 is table2

    def test_start_codon(self):
        """Testa identificação de códon de início."""
        table = CodonTable()
        assert table.is_start_codon("AUG")
        assert not table.is_start_codon("UUU")

    def test_stop_codons(self):
        """Testa identificação de códons de parada."""
        table = CodonTable()
        assert table.is_stop_codon("UAA")
        assert table.is_stop_codon("UAG")
        assert table.is_stop_codon("UGA")
        assert not table.is_stop_codon("AUG")

    def test_get_amino_acid(self):
        """Testa obtenção de aminoácidos."""
        table = CodonTable()
        assert table.get_amino_acid("AUG") == "M"
        assert table.get_amino_acid("UUU") == "F"
        assert table.get_amino_acid("GGG") == "G"


class TestTranslationService:
    """Testes de tradução."""
    def test_translate_simple_no_orf(self):
        """Testa tradução simples sem busca de ORF."""
        service = TranslationService()
        result = service.translate("AUGUUUUAA", find_orf=False)

        # AUG=M, UUU=F, UAA=*
        assert result == "MF"

    def test_translate_simple_with_orf(self):
        """Testa tradução com busca de ORF."""
        service = TranslationService()
        # Garbage antes do AUG
        result = service.translate("GCGCAUGUUUUAA", find_orf=True)

        assert result == "MF"
        assert result[0] == "M"  # Começa com Metionina

    def test_translate_no_start_codon(self):
        """Testa tradução sem códon de início."""
        service = TranslationService()
        result = service.translate("UUUGGGCCC", find_orf=True)

        assert result == "Nenhum ORF encontrado"

    def test_translate_only_start_codon(self):
        """Testa tradução com apenas códon de início."""
        service = TranslationService()
        result = service.translate("AUG", find_orf=False)

        assert result == "M"

    def test_translate_multiple_stop_codons(self):
        """Testa que tradução para no primeiro stop."""
        service = TranslationService()
        result = service.translate("AUGUUUUAAGGG", find_orf=False)

        # Para em UAA, não traduz GGG
        assert result == "MF"

    def test_translate_long_sequence(self):
        """Testa tradução de sequência longa."""
        service = TranslationService()
        # 10 códons válidos
        rna = "AUGUUUCCCAAAGGGCCCUUUAAAUAA"
        result = service.translate(rna, find_orf=False)

        # AUG UUU CCC AAA GGG CCC UUU AAA UAA
        # M   F   P   K   G   P   F   K   *
        assert len(result) == 8

    def test_translate_empty_raises_error(self):
        """Testa que sequência vazia lança erro."""
        service = TranslationService()

        with pytest.raises(TranslationError, match="Não há sequência para traduzir"):
            service.translate("")

    def test_translate_none_raises_error(self):
        """Testa que None lança erro."""
        service = TranslationService()

        with pytest.raises(TranslationError, match="Não há sequência para traduzir"):
            service.translate(None)

    def test_translate_invalid_bases_raises_error(self):
        """Testa que bases inválidas lançam erro."""
        service = TranslationService()

        with pytest.raises(TranslationError, match="bases inválidas"):
            service.translate("AUGXYZ")


class TestORFDetection:
    """Testes de detecção de ORFs."""

    def test_find_all_orfs_single(self):
        """Testa detecção de um único ORF."""
        service = TranslationService()
        # ORF: AUG...UAA (34 aminoácidos = 102 nt + 3 stop)
        rna = "AUG" + "AAA" * 33 + "UAA"

        orfs = service.find_all_orfs(rna, min_length=30)

        assert len(orfs) >= 1
        assert orfs[0]['length'] >= 30
        assert orfs[0]['protein'][0] == 'M'

    def test_find_all_orfs_multiple_frames(self):
        """Testa detecção de ORFs em múltiplas fases."""
        service = TranslationService()
        # Frame 0: AUG no início
        # Frame 1: AUG no meio
        rna = "AUGAAA" + "AUG" + "CCC" * 30 + "UAA"

        orfs = service.find_all_orfs(rna, min_length=5)

        # Deve encontrar pelo menos 1 ORF
        assert len(orfs) >= 1

    def test_find_all_orfs_min_length_filter(self):
        """Testa que filtro de tamanho mínimo funciona."""
        service = TranslationService()
        # ORF curto (5 aminoácidos)
        rna = "AUGAAACCCGGGUUUUAA"

        # Com min_length=10, não deve encontrar
        orfs = service.find_all_orfs(rna, min_length=10)
        assert len(orfs) == 0

        # Com min_length=3, deve encontrar
        orfs = service.find_all_orfs(rna, min_length=3)
        assert len(orfs) >= 1

    def test_find_all_orfs_sorted_by_length(self):
        """Testa que ORFs são ordenados por tamanho."""
        service = TranslationService()
        # Criar sequência com múltiplos ORFs de tamanhos diferentes
        # ORF longo seguido de ORF curto
        rna = "AUG" + "AAA" * 40 + "UAA" + "GCGC" + "AUG" + "CCC" * 35 + "UAA"

        orfs = service.find_all_orfs(rna, min_length=30)

        if len(orfs) >= 2:
            # Primeiro deve ser o maior
            assert orfs[0]['length'] >= orfs[1]['length']

    def test_find_all_orfs_no_orfs(self):
        """Testa que retorna lista vazia sem ORFs."""
        service = TranslationService()
        rna = "UUUGGGCCC"  # Sem AUG

        orfs = service.find_all_orfs(rna)

        assert orfs == []


class TestSixFrameTranslation:
    """Testes de tradução em 6 frames."""

    def test_translate_six_frames_structure(self):
        """Testa estrutura de retorno."""
        service = TranslationService()
        dna = "ATGAAACCCGGGTTT"

        result = service.translate_six_frames(dna)

        assert 'forward' in result
        assert 'reverse' in result
        assert len(result['forward']) == 3
        assert len(result['reverse']) == 3

    def test_translate_six_frames_forward(self):
        """Testa que frames forward funcionam."""
        service = TranslationService()
        dna = "ATGAAACCCGGGTTT"

        result = service.translate_six_frames(dna)

        # Todos os 3 frames devem retornar strings
        for frame in result['forward']:
            assert isinstance(frame, str)

    def test_translate_six_frames_reverse(self):
        """Testa que frames reverse funcionam."""
        service = TranslationService()
        dna = "ATGAAACCCGGGTTT"

        result = service.translate_six_frames(dna)

        # Todos os 3 frames devem retornar strings
        for frame in result['reverse']:
            assert isinstance(frame, str)


class TestCodonUsage:
    """Testes de análise de uso de códons."""

    def test_get_codon_usage_simple(self):
        """Testa contagem de uso de códons."""
        service = TranslationService()
        rna = "AUGAAACCCAUG"

        usage = service.get_codon_usage(rna)

        assert usage['AUG'] == 2
        assert usage['AAA'] == 1
        assert usage['CCC'] == 1

    def test_get_codon_usage_empty(self):
        """Testa uso de códons em sequência vazia."""
        service = TranslationService()

        usage = service.get_codon_usage("")

        assert usage == {}

    def test_get_codon_usage_incomplete_codon(self):
        """Testa que códons incompletos são ignorados."""
        service = TranslationService()
        rna = "AUGAA"  # Último códon incompleto

        usage = service.get_codon_usage(rna)

        assert usage == {'AUG': 1}  # AA é ignorado


# Fixtures
@pytest.fixture
def simple_rna():
    """Fixture com RNA simples."""
    return "AUGUUUCCCAAAGGG"

@pytest.fixture
def rna_with_orf():
    """Fixture com RNA contendo ORF claro."""
    return "GCGCGCAUGUUUCCCAAAGGGCCCUAA"

@pytest.fixture
def long_rna():
    """Fixture com RNA longo (300 nt)."""
    return "AUG" + "AAA" * 97 + "UAA"


class TestTranslationIntegration:
    """Testes de integração."""

    def test_complete_workflow(self, rna_with_orf):
        """Testa workflow completo de tradução."""
        service = TranslationService()

        # 1. Traduzir com ORF
        protein = service.translate(rna_with_orf, find_orf=True)
        assert protein[0] == 'M'

        # 2. Encontrar todos os ORFs
        orfs = service.find_all_orfs(rna_with_orf, min_length=5)
        assert len(orfs) >= 1

        # 3. Uso de códons
        usage = service.get_codon_usage(rna_with_orf)
        assert 'AUG' in usage

    def test_long_sequence_performance(self, long_rna):
        """Testa performance com sequência longa."""
        service = TranslationService()

        # Deve completar rapidamente mesmo com 300 nt
        protein = service.translate(long_rna, find_orf=False)

        assert len(protein) >= 90  # ~97 aminoácidos

    def test_real_world_example(self):
        """Testa com exemplo do mundo real."""
        service = TranslationService()

        # Início do gene da insulina humana (simplificado)
        # Preproinsulina
        rna = "AUGGGCUCUGUGGUUCCACGCCUCCAAGCUCCUUCUGGGAGCCCCUG"

        protein = service.translate(rna, find_orf=True)

        # Deve começar com Met
        assert protein[0] == 'M'
        assert len(protein) > 5
