# Manipulação de arquivos e sequências biológicas

Projeto prático desenvolvido como parte dos estudos no PyLadies Bioinfo. Dupla: Ilária ([@ilariacs](https://github.com/ilariacs)) e Letícia ([@leticia](https://github.com/leticia)).

> Criar um pacote que permita a manipulação de arquivos e sequências de DNA, RNA e proteínas.

# Objetivo do Projeto
O objetivo deste projeto é criar um _pipeline_ que faça os passos do Dogma Central da
Biologia Molecular:

> DNA -> transcrição -> RNA -> tradução -> Proteína

# Requerimentos
- Python 3.10+
- Bibliotecas: Biopython, Flask/FastAPI (para a interface web), pytest (para
  testes)
- Ambiente virtual: venv ou conda (opcional, mas recomendado)

# Instalação
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

# Estrutura do projeto
```
manipulaçao-dna-rna-proteina/
├── src/
│   ├── __init__.py
│   ├── core/
│   │   ├── __init__.py
│   │   └── fasta_parser.py
│   ├── services/
│   │   └── __init__.py
│   ├── utils/
│   │   └── __init__.py
│   └── web/
│       └── __init__.py
├── tests/
│   ├── __init__.py
│   └── test_fasta_parser.py
├── .gitignore
├── requirements.txt
├── pytest.ini
├── setup.py
├── LICENSE
└── README.md
```

# Cronograma de desenvolvimento
- Semana 1: Planejamento e definição dos requisitos do pacote
    - Escolher projeto: Ilária e Lele
    - Tirar dúvidas do escopo biológico do projeto: Ilária e Lele
    - Criação do GitHub + [projeto](https://github.com/users/leticia/projects/8): Lele
    - Terminar as aulas + criar a conta no GitHub: Ilária
- Semana 2: Implementação das funções básicas para leitura (e escrita?) de arquivos/sequências
    - Dividir tarefas no GitHub Projects
    - Começar a escrever código (e primeiros testes)
- Semana 3: Adição de funcionalidades para manipulação de sequências (transcrição, tradução, etc. ?).
    - 50% do código pronto
    - Início da documentação
- Semana 4: Testes e documentação do pacote
    - Finalização da documentação
    - Testes automatizados
    - Testes manuais: casos não cobertos pelos testes automatizados

# Apresentação do projeto
- Data: 07/12/2025
