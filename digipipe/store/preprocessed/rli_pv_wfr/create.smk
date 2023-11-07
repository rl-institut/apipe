"""
Snakefile for this dataset

Note: To include the file in the main workflow, it must be added to the respective module.smk .
"""

import json
import re
import geopandas as gpd
from pathlib import Path
from digipipe.store.utils import get_abs_dataset_path

DATASET_PATH = get_abs_dataset_path(
    "preprocessed", "rli_pv_wfr", data_dir=True)

rule extract:
    """
    Extract files: geodata, datapackage and metadata
    """
    input:
        get_abs_dataset_path(
            "raw", "rli_pv_wfr") / "data" /
            "rli_pv_windflaechenrechner_geodaten_v1.0.zip"
    output:
        [DATASET_PATH / f
         for files in config["files_extract"].values() for f in files]
    params:
        outpath=DATASET_PATH,
        files_extract=" ".join(
            [f for files in config["files_extract"].values() for f in files]
        )
    shell:
        """
        unzip {input} {params.files_extract} -d {params.outpath}
        """
