"""
Testes unitários para TranscriptionService.
"""
import pytest
from src.core.transcription_service import TranscriptionService, TranscriptionError


class TestTranscription:
    """Testes de transcrição DNA → RNA."""

    def test_transcribe_simple(self):
        """Testa transcrição simples."""
        result = TranscriptionService.transcribe("ATGC")
        assert result == "AUGC"

    def test_transcribe_only_t_to_u(self):
        """Testa que apenas T é convertido em U."""
        result = TranscriptionService.transcribe("ATGCATGC")
        assert result == "AUGCAUGC"
        assert 'T' not in result
        assert 'U' in result

    def test_transcribe_preserves_other_bases(self):
        """Testa que outras bases (A, G, C) são preservadas."""
        result = TranscriptionService.transcribe("AGCAGC")
        assert result == "AGCAGC"  # Sem T, nada muda

    def test_transcribe_lowercase_converted(self):
        """Testa que minúsculas são convertidas."""
        result = TranscriptionService.transcribe("atgc")
        assert result == "AUGC"

    def test_transcribe_with_whitespace(self):
        """Testa transcrição com espaços."""
        result = TranscriptionService.transcribe("ATG TAC", validate=False)
        assert result == "AUG UAC"

    def test_transcribe_invalid_bases_raises_error(self):
        """Testa que bases inválidas lançam erro."""
        with pytest.raises(TranscriptionError, match="bases inválidas"):
            TranscriptionService.transcribe("ATGXYZ")

    def test_transcribe_with_n(self):
        """Testa transcrição com N (nucleotídeo ambíguo)."""
        result = TranscriptionService.transcribe("ATGCNNN")
        assert result == "AUGCNNN"

    def test_transcribe_empty_raises_error(self):
        """Testa que sequência vazia lança erro."""
        with pytest.raises(TranscriptionError, match="vazia"):
            TranscriptionService.transcribe("")

    def test_transcribe_none_raises_error(self):
        """Testa que None lança erro."""
        with pytest.raises(ValueError, match="não pode ser None"):
            TranscriptionService.transcribe(None)

    def test_transcribe_non_string_raises_error(self):
        """Testa que não-string lança erro."""
        with pytest.raises(ValueError, match="deve ser string"):
            TranscriptionService.transcribe(12345)

    def test_transcribe_without_validation(self):
        """Testa transcrição sem validação."""
        # Permite bases inválidas se validate=False
        result = TranscriptionService.transcribe("ATXYZ", validate=False)
        assert result == "AUXYZ"


class TestReverseTranscription:
    """Testes de transcrição reversa RNA → DNA."""

    def test_reverse_transcribe_simple(self):
        """Testa transcrição reversa simples."""
        result = TranscriptionService.reverse_transcribe("AUGC")
        assert result == "ATGC"

    def test_reverse_transcribe_only_u_to_t(self):
        """Testa que apenas U é convertido em T."""
        result = TranscriptionService.reverse_transcribe("AUGCAUGC")
        assert result == "ATGCATGC"
        assert 'U' not in result
        assert 'T' in result

    def test_reverse_transcribe_invalid_bases_raises_error(self):
        """Testa que bases inválidas lançam erro."""
        with pytest.raises(TranscriptionError):
            TranscriptionService.reverse_transcribe("AUGXYZ")

    def test_reverse_transcribe_empty_raises_error(self):
        """Testa que sequência vazia lança erro."""
        with pytest.raises(TranscriptionError):
            TranscriptionService.reverse_transcribe("")


class TestComplement:
    """Testes de sequência complementar."""

    def test_dna_complement(self):
        """Testa complemento de DNA."""
        result = TranscriptionService.get_complement("ATGC", "dna")
        assert result == "TACG"

    def test_rna_complement(self):
        """Testa complemento de RNA."""
        result = TranscriptionService.get_complement("AUGC", "rna")
        assert result == "UACG"

    def test_complement_self_complementary(self):
        """Testa que GC é auto-complementar."""
        result = TranscriptionService.get_complement("GCGC", "dna")
        assert result == "CGCG"

    def test_complement_invalid_type_raises_error(self):
        """Testa que tipo inválido lança erro."""
        with pytest.raises(ValueError, match="Tipo de sequência inválido"):
            TranscriptionService.get_complement("ATGC", "invalid")


class TestReverseComplement:
    """Testes de sequência reverso-complementar."""

    def test_reverse_complement_dna(self):
        """Testa reverso-complemento de DNA."""
        result = TranscriptionService.get_reverse_complement("ATGC", "dna")
        # Complemento: TACG, Reverso: GCAT
        assert result == "GCAT"

    def test_reverse_complement_rna(self):
        """Testa reverso-complemento de RNA."""
        result = TranscriptionService.get_reverse_complement("AUGC", "rna")
        # Complemento: UACG, Reverso: GCAU
        assert result == "GCAU"

    def test_reverse_complement_palindrome(self):
        """Testa reverso-complemento de palíndromo."""
        # GAATTC é palíndromo (sítio de restrição EcoRI)
        result = TranscriptionService.get_reverse_complement("GAATTC", "dna")
        assert result == "GAATTC"


# Fixtures
@pytest.fixture
def dna_sample():
    """Fixture com sequência DNA de exemplo."""
    return "ATGCGATCGATCGATCG"

@pytest.fixture
def rna_sample():
    """Fixture com sequência RNA de exemplo."""
    return "AUGCGAUCGAUCGAUCG"


class TestTranscriptionIntegration:
    """Testes de integração."""

    def test_transcribe_and_reverse(self, dna_sample):
        """Testa transcrição ida e volta."""
        rna = TranscriptionService.transcribe(dna_sample)
        dna_back = TranscriptionService.reverse_transcribe(rna)

        assert dna_back == dna_sample

    def test_complement_and_reverse_complement(self, dna_sample):
        """Testa complemento e reverso-complemento."""
        complement = TranscriptionService.get_complement(dna_sample, "dna")
        rev_comp = TranscriptionService.get_reverse_complement(dna_sample, "dna")

        # Reverso do complemento deve ser igual ao reverso-complemento
        assert complement[::-1] == rev_comp

    def test_transcribe_long_sequence(self):
        """Testa transcrição de sequência longa."""
        # Sequência de 1000 bases
        dna = "ATGC" * 250
        rna = TranscriptionService.transcribe(dna)

        assert len(rna) == 1000
        assert rna.count('U') == 250
        assert 'T' not in rna
