import pytest
from src.core.fasta_parser import FastaParser, FastaParseError

class TestFastaParser:
    def test_parse_valid_fasta(self):
        """Testa parsing de FASTA válido."""
        fasta = ">gene1 description\nATGCGATCG\nATCGATCG"
        result = FastaParser.parse(fasta)

        assert result['header'] == "gene1 description"
        assert result['sequence'] == "ATGCGATCGATCGATCG"
        assert result['id'] == "gene1"
        assert result['description'] == "description"

    def test_parse_multiline_sequence(self):
        """Testa parsing com sequência em múltiplas linhas."""
        fasta = ">gene1\nATGC\nGATC\nTAGC"
        result = FastaParser.parse(fasta)

        assert result['sequence'] == "ATGCGATCTAGC"
        assert len(result['sequence']) == 12

    def test_parse_lowercase_converted_to_uppercase(self):
        """Testa que sequência minúscula é convertida para maiúscula."""
        fasta = ">gene1\natgcgatcg"
        result = FastaParser.parse(fasta)

        assert result['sequence'] == "ATGCGATCG"
        assert result['sequence'].isupper()

    def test_parse_with_spaces_in_sequence(self):
        """Testa que espaços na sequência são removidos."""
        fasta = ">gene1\nATGC GATC TAGC"
        result = FastaParser.parse(fasta)

        assert result['sequence'] == "ATGCGATCTAGC"
        assert ' ' not in result['sequence']

    def test_parse_header_without_description(self):
        """Testa parsing de header sem descrição."""
        fasta = ">gene1\nATGC"
        result = FastaParser.parse(fasta)

        assert result['id'] == "gene1"
        assert result['description'] == ""

    def test_parse_complex_header(self):
        """Testa parsing de header complexo com múltiplas palavras."""
        fasta = ">NM_004006.3 Homo sapiens dystrophin (DMD), transcript variant 1"
        fasta += "\nATGC"
        result = FastaParser.parse(fasta)

        assert result['id'] == "NM_004006.3"
        assert "Homo sapiens" in result['description']

class TestFastaParseErrors:
    def test_parse_empty_string_raises_error(self):
        """Testa que string vazia lança exceção."""
        with pytest.raises(FastaParseError, match="Entrada FASTA vazia ou conteúdo inválido."):
            FastaParser.parse("")

    def test_parse_none_raises_error(self):
        """Testa que None lança exceção."""
        with pytest.raises(FastaParseError, match="Entrada FASTA vazia ou conteúdo inválido."):
            FastaParser.parse(None)

    def test_parse_non_string_raises_error(self):
        """Testa que tipo não-string lança exceção."""
        with pytest.raises(FastaParseError, match="Entrada FASTA vazia ou conteúdo inválido."):
            FastaParser.parse(12345)

    def test_parse_no_header_raises_error(self):
        """Testa que FASTA sem '>' lança exceção."""
        with pytest.raises(FastaParseError, match="Formato FASTA inválido: cabeçalho malformatado ou ausente."):
            FastaParser.parse("ATGCGATCG")

    def test_parse_empty_header_raises_error(self):
        """Testa que header vazio lança exceção."""
        with pytest.raises(FastaParseError, match="Header FASTA está vazio após '>'"):
            FastaParser.parse(">\nATGC")

    def test_parse_no_sequence_raises_error(self):
        """Testa que FASTA sem sequência lança exceção."""
        with pytest.raises(FastaParseError, match="Formato FASTA inválido: sequência ausente."):
            FastaParser.parse(">gene1")

    def test_parse_only_whitespace_sequence_raises_error(self):
        """Testa que sequência apenas com espaços lança exceção."""
        with pytest.raises(FastaParseError, match="Formato FASTA inválido: sequência ausente."):
            FastaParser.parse(">gene1\n   \n  ")

class TestSequenceValidation:
    def test_validate_dna_valid(self):
        """Testa validação de sequência DNA válida."""
        assert FastaParser.validate_sequence("ATGCGATCG", "dna")
        assert FastaParser.validate_sequence("NNNNATGC", "dna")

    def test_validate_dna_invalid(self):
        """Testa validação de DNA inválido."""
        assert not FastaParser.validate_sequence("ATGCUXYZ", "dna")
        assert not FastaParser.validate_sequence("ATGCMKWS", "dna")

    def test_validate_invalid_type_raises_error(self):
        """Testa que tipo inválido lança exceção."""
        with pytest.raises(ValueError, match="Tipo de sequência inválido"):
            FastaParser.validate_sequence("ATGC", "invalid_type")

    def test_validate_rna_valid(self):
        """Testa validação de RNA válido."""
        assert FastaParser.validate_sequence("AUGCGAUCG", "rna")
        assert not FastaParser.validate_sequence("ATGC", "rna")

    def test_validate_protein_valid(self):
        """Testa validação de proteína válida."""
        assert FastaParser.validate_sequence("MFKGDWIV", "protein")
        assert FastaParser.validate_sequence("ACDEFGHIKLMNPQRSTVWY", "protein")

    def test_validate_protein_with_stop(self):
        """Testa validação de proteína com códon de parada."""
        assert FastaParser.validate_sequence("MFKG*", "protein")

class TestSequenceTypeDetection:
    """Testes de detecção automática de tipo."""

    def test_detect_dna(self):
        """Testa detecção de DNA."""
        assert FastaParser.detect_sequence_type("ATGCGATCG") == "dna"
        assert FastaParser.detect_sequence_type("NNNATGC") == "dna"

    def test_detect_rna(self):
        """Testa detecção de RNA."""
        assert FastaParser.detect_sequence_type("AUGCGAUCG") == "rna"

    def test_detect_protein(self):
        """Testa detecção de proteína."""
        assert FastaParser.detect_sequence_type("MFKGDWIV") == "protein"
        assert FastaParser.detect_sequence_type("ACDEFGHIKLMNPQRSTVWY") == "protein"

    def test_detect_unknown(self):
        """Testa detecção de tipo desconhecido."""
        result = FastaParser.detect_sequence_type("XYZ123")
        assert result == "unknown"
