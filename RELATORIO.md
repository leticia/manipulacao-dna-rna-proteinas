# Analisador de Arquivos FASTA

## 1. Processo de Desenvolvimento

O desenvolvimento do analisador foi feito em módulos e seguindo uma abordagem orientada a objetos, dividindo o problema em componentes independentes e testáveis.

### Etapas do desenvolvimento:

1. **Planejamento inicial**: Definição do problema a ser resolvido (simulação
   do Dogma Central da Biologia Molecular).

2. **Análise de requisitos**: Identificação das funcionalidades necessárias (leitura FASTA, estatísticas, transcrição, tradução, geração de relatório HTML).

2. **Design da arquitetura**: Divisão em módulos separados, cada um com uma responsabilidade única:
   - `fasta_reader.py`: Leitura de arquivos FASTA
   - `sequence_analyzer.py`: Cálculos estatísticos
   - `transcription.py`: Transcrição DNA → RNA
   - `translation.py`: Tradução RNA → Proteína
   - `html_generator.py`: Geração de relatórios

3. **Implementação**: Codificação dos módulos seguindo boas práticas (nomes
   descritivos, legibilidade).

4. **Testes unitários**: Criação de testes para cada módulo, cobrindo casos normais e extremos.

5. **Integração**: Criação do script principal (main.py) que coordena todos os módulos.

6. **Documentação**: Comentários no código e README com instruções de uso.

---

## 2. Design Patterns Utilizados

### 2.1 Single Responsibility Principle (SRP)
O SRP estabelece que uma classe deve ter apenas uma razão para mudar. Cada módulo tem uma função específica e não mistura responsabilidades.

Exemplos: `sequence_analyzer.py, transcription.py, translation.py`.

Cada classe tem uma única responsabilidade bem definida:
- `SequenceAnalyzer`: apenas calcula estatísticas
- `Transcription`: apenas faz transcrição
- `Translation`: apenas faz tradução

Isso facilita manutenção, testes e extensão. Se mudar a regra de tradução, apenas `translation.py` precisa ser modificado.

---

### 2.2 Strategy Pattern
Define uma família de algoritmos, encapsula cada um deles e os torna intercambiáveis. Permite que o algoritmo varie independentemente dos clientes que o utilizam.

Exemplo: `fasta_reader.py`

O algoritmo de leitura de FASTA na classe `FastaReader` fica encapsulado. No futuro, poderíamos ter `FastqReader`, `GenbankReader`, etc., sem modificar o código existente (Open/Closed Principle).

---

### 2.3 Template Method Pattern
Define o esqueleto de um algoritmo em uma operação, delegando alguns passos para as subclasses. Neste caso, o template HTML é fixo, mas os dados são variáveis.

Exemplo: `html_generator.py`

A classe `HTMLGenerator` define a estrutura fixa do HTML (template), mas permite que o conteúdo específico seja preenchido dinamicamente.

---

## 3. Limitações do Projeto

### 3.1 Limitações biológicas

- **Fase de leitura única**: O código traduz a partir da primeira base, sem buscar ORFs (Open Reading Frames) alternativos.
- **Códon de início fixo**: Assume que a tradução sempre começa no início da sequência, sem buscar AUG.
- **Código genético padrão**: Não suporta códigos genéticos alternativos (mitocondrial, bacteriano).
- **Splicing não implementado**: Se o FASTA contiver sequência genômica (com íntrons), a tradução será incorreta.

### 3.2 Limitações técnicas

- **Performance**: Não otimizado para sequências muito longas (>100MB).
- **Validação limitada**: Não valida completamente se a sequência contém apenas nucleotídeos válidos.
- **Formato único**: Suporta apenas FASTA, não outros formatos bioinformáticos.
- **Estatísticas básicas**: Calcula apenas tamanho, contagens e GC%. Não calcula complexidade, repetições, etc.

### 3.3 Limitações de interface

- **CLI apenas**: Interface apenas por linha de comando.
- **HTML estático**: Relatório gerado é estático, sem interatividade.
- **Sem visualização gráfica**: Não gera gráficos de composição, distribuição de GC ao longo da sequência, etc.

---

## 4. Sugestões para o Futuro

### 4.1 Melhorias bioinformáticas

1. **Detecção de ORFs**: Buscar todos os possíveis frames de leitura e identificar ORFs válidos (começam com AUG, terminam com stop codon).

2. **Análise de múltiplos frames**: Traduzir nos 6 frames de leitura possíveis (3 forward + 3 reverse complement).

3. **Suporte a splicing**: Processar arquivos GFF/GTF para remover íntrons antes da tradução.

4. **Códigos genéticos alternativos**: Permitir seleção de diferentes tabelas de códons.

5. **Alinhamento de sequências**: Comparar a proteína traduzida com bancos de dados (BLAST local).

### 4.2 Melhorias técnicas

1. **Validação**: Verificar nucleotídeos inválidos, sequências ambíguas (N, R, Y, etc.).

2. **Performance**: Usar generators/iterators para processar sequências muito grandes sem carregar tudo na memória.

3. **Suporte a múltiplos formatos**: FASTQ, GenBank, EMBL, etc.

4. **Logging**: Sistema de logs para debugging e monitoramento.

5. **Configuração externa**: Arquivo de configuração para parâmetros (tabela de códons, cutoffs, etc.).

### 4.3 Melhorias de interface

1. **GUI**: Interface gráfica?

2. **Web app**: Versão web com Flask/FastAPI + frontend React.

3. **Visualizações**: Gráficos interativos com Plotly ou Bokeh:
   - Distribuição de GC ao longo da sequência (sliding window)
   - Composição de nucleotídeos/aminoácidos (pie charts)

4. **Batch processing**: Processar múltiplos arquivos FASTA de uma vez (em um
   futuro mais distante).
