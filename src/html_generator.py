"""
Módulo para geração de relatórios HTML.
Template Method Pattern: estrutura fixa com conteúdo variável.
"""


class HTMLGenerator:
    """
    Gerador de relatórios HTML.

    Template Method Pattern: Define a estrutura do HTML (template),
    permitindo que o conteúdo seja preenchido dinamicamente.
    """

    def generate(self, header: str, dna_sequence: str, rna_sequence: str,
                 protein_sequence: str, stats: dict) -> str:
        """
        Gera relatório HTML completo.

        Args:
            header: Cabeçalho do FASTA
            dna_sequence: Sequência de DNA
            rna_sequence: Sequência de RNA
            protein_sequence: Sequência de proteína
            stats: Dicionário com estatísticas

        Returns:
            String com HTML completo
        """
        # Template HTML
        html_parts = []
        html_parts.append('<!DOCTYPE html>')
        html_parts.append('<html lang="pt-BR">')
        html_parts.append('<head>')
        html_parts.append('    <meta charset="UTF-8">')
        html_parts.append('    <meta name="viewport" content="width=device-width, initial-scale=1.0">')
        html_parts.append('    <title>Relatório estatístico</title>')
        html_parts.append('    <style>')
        html_parts.append('        body { font-family: Arial, sans-serif; max-width: 1200px; margin: 0 auto; padding: 20px; background: #f5f5f5; }')
        html_parts.append('        .container { background: white; padding: 30px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }')
        html_parts.append('        h1 { color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 10px; }')
        html_parts.append('        h2 { color: #34495e; margin-top: 30px; }')
        html_parts.append('        .stats { background: #ecf0f1; padding: 15px; border-radius: 5px; margin: 20px 0; }')
        html_parts.append('        .stats-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 15px; margin-top: 10px; }')
        html_parts.append('        .stat-item { background: white; padding: 10px; border-radius: 4px; border-left: 4px solid #3498db; }')
        html_parts.append('        .stat-label { font-weight: bold; color: #7f8c8d; font-size: 0.9em; }')
        html_parts.append('        .stat-value { font-size: 1.5em; color: #2c3e50; margin-top: 5px; }')
        html_parts.append('        .sequence-box { background: #2c3e50; color: #2ecc71; padding: 15px; border-radius: 5px; font-family: monospace; overflow-x: auto; word-break: break-all; max-height: 200px; overflow-y: auto; }')
        html_parts.append('        .info { background: #e8f4f8; border-left: 4px solid #3498db; padding: 10px 15px; margin: 15px 0; }')
        html_parts.append('    </style>')
        html_parts.append('</head>')
        html_parts.append('<body>')
        html_parts.append('    <div class="container">')
        html_parts.append('        <h1>Relatório</h1>')
        html_parts.append('        <div class="info">')
        html_parts.append(f'            <strong>Sequência:</strong> {self._escape_html(header)}')
        html_parts.append('        </div>')
        html_parts.append('        <h2>Estatísticas da Sequência</h2>')
        html_parts.append('        <div class="stats">')
        html_parts.append('            <div class="stats-grid">')
        html_parts.append('                <div class="stat-item">')
        html_parts.append('                    <div class="stat-label">Tamanho</div>')
        html_parts.append(f'                    <div class="stat-value">{stats["length"]:,} bp</div>')
        html_parts.append('                </div>')
        html_parts.append('                <div class="stat-item">')
        html_parts.append('                    <div class="stat-label">GC%</div>')
        html_parts.append(f'                    <div class="stat-value">{stats["gc_percent"]:.2f}%</div>')
        html_parts.append('                </div>')
        html_parts.append(self._generate_nucleotide_counts(stats['counts']))
        html_parts.append('            </div>')
        html_parts.append('        </div>')
        html_parts.append('        <h2>Sequência de DNA</h2>')
        html_parts.append('        <div class="sequence-box">')
        html_parts.append(self._format_sequence(dna_sequence, 60))
        html_parts.append('        </div>')
        html_parts.append('        <h2>Sequência de RNA (Transcrita)</h2>')
        html_parts.append('        <div class="sequence-box">')
        html_parts.append(self._format_sequence(rna_sequence, 60))
        html_parts.append('        </div>')
        html_parts.append('        <h2>Sequência de Proteína (Traduzida)</h2>')
        html_parts.append('        <div class="sequence-box">')
        html_parts.append(self._format_sequence(protein_sequence, 60))
        html_parts.append('        </div>')
        html_parts.append('        <div class="info" style="margin-top: 30px;">')
        html_parts.append('            <strong>Nota:</strong> A tradução usa o código genético padrão.')
        html_parts.append('        </div>')
        html_parts.append('    </div>')
        html_parts.append('</body>')
        html_parts.append('</html>')

        return '\n'.join(html_parts)

    def _escape_html(self, text: str) -> str:
        """Escapa caracteres HTML especiais."""
        return (text.replace('&', '&amp;')
                   .replace('<', '&lt;')
                   .replace('>', '&gt;'))

    def _format_sequence(self, sequence: str, width: int) -> str:
        """
        Formata sequência em linhas de largura fixa.
        """
        lines = []
        for i in range(0, len(sequence), width):
            lines.append(sequence[i:i+width])
        return '\n'.join(lines)

    def _generate_nucleotide_counts(self, counts: dict) -> str:
        """
        Gera HTML para contagens de nucleotídeos.
        """
        html_parts = []
        for nt, count in sorted(counts.items()):
            html_parts.append('                <div class="stat-item">')
            html_parts.append(f'                    <div class="stat-label">{nt}</div>')
            html_parts.append(f'                    <div class="stat-value">{count:,}</div>')
            html_parts.append('                </div>')
        return '\n'.join(html_parts)
