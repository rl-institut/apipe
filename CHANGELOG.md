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

### Changed

- Move dataset docs from config to md files
- Retain mastr_id in MaStR datasets
- Fix loading of empty yml files
- Fix loading of global config when in workflow dir
- The package management in digipipe has been changed to poetry.
- The installation of a virtual environment is done only from the environment.yml file and via conda.

### Removed

- setup.py and requirements.txt files are omitted with the conversion to poetry
