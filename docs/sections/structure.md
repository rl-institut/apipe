# Struktur des Repos

Das Repository ist wie folgt aufgebaut:
```
.
├── digipipe
│   ├── config
│   │   └── global.yml                          # Global config
│   ├── logs                                    # Place for log files
│   ├── scenarios                               # Scenario definition
│   │   ├── .TEMPLATE                           # - Template
│   │   └── SCENARIOS.md
│   ├── scripts                                 # Main scripts
│   │   ├── datasets                            # - for data processing which is shared by different datasets
│   │   ├── esm                                 # - for energy system modelling
│   │   ├── config.py                           # - config-related functions
│   │   └── geo.py                              # - spatial functions
│   ├── store                                   # Data store
│   │   ├── appdata                             # - App-ready data
│   │   │   ├── data
│   │   │   ├── metadata
│   │   │   └── scenarios
│   │   ├── datasets                            # - Processed datasets
│   │   │   ├── bkg_vg250_districts_abw
│   │   │   ├── bkg_vg250_muns_abw
│   │   │   ├── bkg_vg250_region_abw
│   │   │   ├── osm_forest
│   │   │   ├── .TEMPLATE
│   │   │   └── module.smk
│   │   ├── preprocessed                        # - Preprocessed datasets
│   │   │   ├── bkg_vg250
│   │   │   ├── osm_filtered
│   │   │   ├── .TEMPLATE
│   │   │   └── module.smk
│   │   ├── raw                                 # - Raw datasets
│   │   │   ├── bkg_vg250
│   │   │   ├── destatis_gv
│   │   │   ├── osm_sachsen-anhalt
│   │   │   └── .TEMPLATE
│   │   ├── temp                                # - Temporary files
│   │   ├── DATASETS.md
│   │   └── utils.py
│   └── workflow
│       ├── Snakefile                           # Main snakefile
│       ├── utils.py
│       └── WORKFLOW.md
├── docs                                        # Documentation
│   ├── sections
│   │   ├── data.rst
│   │   ├── installation.rst
│   │   ├── scenarios.rst
│   │   ├── structure.rst
│   │   └── workflow.rst
│   ├── conf.py
│   ├── index.rst
│   └── Makefile
├── tests
├── CHANGELOG.md
├── CONTRIBUTING.md
├── environment.yml
├── LICENSE
├── README.md
├── requirements.txt
└── setup.py
```

Erstellt mit:
```
tree --dirsfirst -L 4 -a -I '__*|*log|.gitkeep|PKG-INFO|*egg-info*|img|.git|.idea|venv|.snakemake' . > dirtree.txt
```