# Manipulação de arquivos e sequências biológicas

Projeto prático desenvolvido como parte dos estudos no PyLadies Bioinfo. Dupla: Ilária ([@ilariacs](https://github.com/ilariacs)) e Letícia ([@leticia](https://github.com/leticia)).

> Criar um pacote que permita a manipulação de arquivos e sequências de DNA, RNA e proteínas.

# Objetivo do Projeto
O objetivo deste projeto é criar um _pipeline_ que faça os passos do Dogma Central da
Biologia Molecular:

![Dogma Central da Biologia](dogma_central.png "Dogma Central da Biologia")

# Estrutura do projeto
```
manipulaçao-dna-rna-proteina/
├── .gitignore
├── main.py
├── requirements.txt
├── LICENSE
├── README.md
├── pytest.ini
├── src/
│   ├── __init__.py
│   ├── fasta_reader.py
│   ├── sequence_analyzer.py
│   ├── transcription.py
│   ├── translation.py
│   └── html_generator.py
└── tests/
    ├── __init__.py
    ├── test_fasta_reader.py
    ├── test_sequence_analyzer.py
    ├── test_transcription.py
    └── test_translation.py
```

## Requerimentos
- Python 3.10+
- Bibliotecas: pytest, pycov (para testes)
- Ambiente virtual: venv ou conda (opcional, mas recomendado)

## Instalação
1. Clone o repositório;
2. Crie um ambiente virtual (se usar venv ou conda) e ative-o:
   ```bash
   python -m venv venv
   source venv/bin/activate  # No Windows use `venv\Scripts\activate`
   ```
3. Instale as dependências:
    ```bash
    pip install -r requirements.txt
    ```

## Uso

```bash
python main.py caminho/para/arquivo.fasta
```

## Testes

```bash
pytest
```

# Cronograma de desenvolvimento
- Semana 1: Planejamento e definição dos requisitos do pacote ✅
    - Escolher projeto: Ilária e Lele
    - Tirar dúvidas do escopo biológico do projeto: Ilária e Lele
    - Criação do GitHub + [projeto](https://github.com/users/leticia/projects/8): Lele
    - Terminar as aulas + criar a conta no GitHub: Ilária
- Semana 2: Implementação das funções básicas para leitura (e escrita?) de arquivos/sequências ⚠️
    - Dividir tarefas no GitHub Projects
    - Começar a escrever código (e primeiros testes)
- Semana 3: Adição de funcionalidades para manipulação de sequências (transcrição, tradução, etc. ?). ❌
    - 50% do código pronto
    - Início da documentação
- Semana 4: Testes e documentação do pacote ⚠️
    - Finalização da documentação
    - Testes automatizados
    - Testes manuais: casos não cobertos pelos testes automatizados
    - Últimos ajustes
