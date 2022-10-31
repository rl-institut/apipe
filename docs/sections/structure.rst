.. _structure_label:

Struktur des Repos
==================

Das Repository ist wie folgt aufgebaut:

.. code-block::

    .
    ├── digipipe
    │   ├── config
    │   │   └── global.yml              # Global config
    │   ├── scenarios                   # Scenario definition
    │   │   ├── .TEMPLATE
    │   │   └── SCENARIOS.md
    │   ├── scripts                     # Main scripts
    │   │   ├── esm                     # - for energy system modelling
    │   ├── store                       # Data store
    │   │   ├── 0_raw                   # - Raw datasets
    │   │   │   ├── .TEMPLATE
    │   │   ├── 1_preprocessed          # - Preprocessed datasets
    │   │   │   └── .TEMPLATE
    │   │   ├── 2_datasets              # - Processed datasets
    │   │   │   └── .TEMPLATE
    │   │   ├── 3_appdata               # - App data
    │   │   │   ├── data
    │   │   │   ├── metadata
    │   │   │   └── scenarios
    │   │   ├── temp                    # - Store temp files here
    │   │   └── DATASETS.md
    │   ├── workflow
    │   │   ├── helpers.py
    │   │   ├── Snakefile               # Main snakefile
    │   │   └── WORKFLOW.md
    ├── docs                            # ReadTheDocs content
    ├── .github
    ├── tests
    ├── CHANGELOG.md
    ├── CONTRIBUTING.md
    ├── environment.yml
    ├── .gitignore
    ├── LICENSE
    ├── README.md
    ├── requirements.txt
    └── setup.py


Erstellt mit:
.. code-block::

    tree --dirsfirst -L 4 -a -I '__*|*log|.gitkeep|PKG-INFO|*egg-info*|img|.git|.idea|venv|.snakemake' . > dirtree.txt