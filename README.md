<p align="left">
    <img alt="Digiplan logo" align="left" height="150" src="docs/img/logos/digiplan-logo.png">
    <img alt="RLI logo" align="right" height="140" src="docs/img/logos/rli_logo.png">
</p>
<br/><br/><br/><br/><br/><br/><br/>

----------

# Digipipe

Pipeline for data and energy system in the Digiplan project.

## Installation 

**Note: Linux only, Windows is currently not supported.**

First, clone via SSH using

    git clone git@github.com:rl-institut-private/digipipe.git /local/path/to/digipipe/

### Install using pip:

Make sure you have Python >= 3.6 installed, let's create a virtual env:

    virtualenv --python=python3.8 venv
    source venv/bin/activate

Some additional system packages are required, install them by

    apt install gdal-bin python-gdal libspatialindex-dev imagemagick

Notes:
- Make sure you have GDAL>=3.0 as older versions will not work
- `imagemagick` is optional and only required for report creation

Install package with

    pip install -e /local/path/to/djagora_data/

### Install using conda

Make sure you have conda installed, e.g. miniconda. Then create the env:
    
    conda create -n digipipe /local/path/to/digipipe/environment.yml
    conda activate digipipe

## Further reading on structure, pipeline and conventions

- Datasets/data flow: [DATASETS.md](digipipe/store/DATASETS.md)
- Workflow: [WORKFLOW.md](digipipe/workflow/WORKFLOW.md)
- Scenarios: [SCENARIOS.md](digipipe/scenarios/SCENARIOS.md)

## Runtime and resources

**TO BE UPDATED**

**Warning:** A full pipeline run takes 10 hours on a Xeon E5-2690 using 14
cores and requires about 600 GB of disk space.
