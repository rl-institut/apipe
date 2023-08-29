# Installation

**Note: Linux only, Windows is currently not supported.**

First, clone via SSH using

    git clone git@github.com:rl-institut-private/digipipe.git /local/path/to/digipipe/

Prerequisite:

- _conda_ ([Anaconda](https://docs.anaconda.com/anaconda/install/)
/[Miniconda](https://conda.io/en/latest/miniconda.html))
- `python <https://www.python.org/downloads/>`_

Enter repo folder. Set up a conda environment and activate it with:

```
conda env create -f environment.yml
conda activate digipipe
```

Install [poetry](https://python-poetry.org/) (python dependency manager used
in this project) and dependencies for the project (Note: Installing poetry via
pip into same environment is not recommended and can cause trouble! Instead, it
should be installed system-wide via command below or
[pipx](https://python-poetry.org/docs/#installing-with-pipx)):

```
curl -sSL https://install.python-poetry.org | python3 -
poetry install
```

Some additional system packages are required, install them by

    sudo apt install gdal-bin python3-gdal libspatialindex-dev imagemagick osmium-tool graphviz graphviz-dev

Notes:

- Make sure you have GDAL>=3.0 as older versions will not work.
- `imagemagick` is optional and only required for report creation

## Contributing to digipipe

You can write [issues](https://github.com/rl-institut-private/digipipe/issues>)
to announce bugs or to propose enhancements.

If you want to participate in the development of digipipe, please make sure you
use pre-commit.

You activate it with:

    pre-commit install

To trigger a check manually, execute:

    pre-commit run -a

## Runtime and resources

**Warning:** Conversion and extraction process needs ~50 GB disk space and may
take a lot of time!
