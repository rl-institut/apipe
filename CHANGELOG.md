# Changelog
All notable changes to this project will be documented in this file.

The format is inspired from [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and the versioning aim to respect [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- Setup initial structure and files
- Add first bunch of datasets for testing the data flow
- Modularize datasets
- Add first draft of RTD docs
- Add dataset: BNetzA MaStR
- Add Nominatim geocoder
- Add dataset: population
- Clean rule
- Datasets attribute captions
- Create list of region-specific datasets in the docs
- pyproject.toml and poetry.lock file have been added with the conversion to poetry
- Add pre-commit in order to check for errors and linting bugs before commits
- Add types-pyyaml package
- Integrate building of energy system for appdata in pipeline via dir *esys*
- Update store with dir structure for *esys* data

### Changed

- Move dataset docs from config to md files
- Retain mastr_id in MaStR datasets
- Fix loading of empty yml files
- Fix loading of global config when in workflow dir
- Integrate esys Snakefile in workflow Snakefile and update clean rule
- Fix shapely deprecation warning
- Fix ogr2ogr conversion with recent GDAL version (v3.6.2)
- Fix conda installation by removing python-gdal from environment.yml
- The package management in digipipe has been changed to poetry.
- The installation of a virtual environment is done only from the environment.yml file and via conda.
- Apply linters on repo among others: black, isort, check-json and end-of-file-fixer
- Fix C419 flake8 error
- Exchange the Test_scenario with the scenarios status_quo_scenario and
  2045_scenario of *esys* in digipipe/esys/scenarios directory

### Removed

- setup.py and requirements.txt files are omitted with the conversion to poetry
