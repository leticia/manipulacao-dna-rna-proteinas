"""
Ponto de entrada principal do analisador de arquivos FASTA.
Coordena a leitura, análise e geração de relatório HTML.
"""

import sys
from pathlib import Path
from src.fasta_reader import FastaReader
from src.sequence_analyzer import SequenceAnalyzer
from src.transcription import Transcription
from src.translation import Translation
#from src.html_generator import HTMLGenerator

def main():
    """
    Função principal que chama o fluxo de análise.
    """
    # Validação de argumentos
    if len(sys.argv) != 2:
        print("Uso: python main.py <arquivo_fasta>")
        sys.exit(1)

    fasta_path = Path(sys.argv[1])

    # Validação de arquivo
    if not fasta_path.exists():
        print(f"Erro: Arquivo {fasta_path} não encontrado.")
        sys.exit(1)

    # 1. Leitura do arquivo FASTA (Strategy Pattern)
    reader = FastaReader()
    sequence_data = reader.read(str(fasta_path))

    if not sequence_data:
        print("Erro: Não foi possível ler o arquivo FASTA.")
        sys.exit(1)

    header = sequence_data['header']
    sequence = sequence_data['sequence']

    print(f"Sequência carregada: {header}")
    print(f"Tamanho: {len(sequence)} bases\n")

    # 2. Análise estatística (Single Responsibility)
    analyzer = SequenceAnalyzer()
    stats = analyzer.analyze(sequence)

    print("Estatísticas:")
    print(f"  Tamanho: {stats['length']}")
    print(f"  GC%: {stats['gc_percent']:.2f}%")
    print(f"  Contagens: {stats['counts']}\n")

    # 3. Transcrição DNA -> RNA (Single Responsibility)
    transcription = Transcription()
    rna_sequence = transcription.transcribe(sequence)

    print(f"RNA (primeiros 60 nt): {rna_sequence[:60]}...\n")

    # 4. Tradução RNA -> Proteína (Single Responsibility)
    translation = Translation()
    protein_sequence = translation.translate(rna_sequence)

    print(f"Proteína (primeiros 20 aa): {protein_sequence[:20]}...\n")

    # 5. Geração do relatório HTML (Template Method Pattern)
    generator = HTMLGenerator()
    html_content = generator.generate(
        header=header,
        dna_sequence=sequence,
        rna_sequence=rna_sequence,
        protein_sequence=protein_sequence,
        stats=stats
    )

    # Salvamento do relatório
    output_path = Path("resultado_analise.html")
    output_path.write_text(html_content, encoding='utf-8')

    print(f"Relatório gerado: {output_path.absolute()}")


if __name__ == "__main__":
    main()
