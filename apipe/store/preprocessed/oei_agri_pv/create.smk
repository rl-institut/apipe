"""
Snakefile for this dataset

Note: To include the file in the main workflow, it must be added to the respective module.smk .
"""

from apipe.store.utils import get_abs_dataset_path
from pathlib import Path

DATASET_PATH = get_abs_dataset_path(
    "preprocessed", "oei_agri_pv", data_dir=True
)
BASE_INPUT_PATH = get_abs_dataset_path("raw", "oei_agri_pv") / "data"

DATASET_PATH = get_abs_dataset_path(
    "preprocessed", "oei_agri_pv", data_dir=True
)
BASE_INPUT_PATH = get_abs_dataset_path("raw", "oei_agri_pv") / "data"

rule copy_datasets:
    input:
        expand(BASE_INPUT_PATH / "{file}.tif", file=config["datasets"])
    output:
        expand(DATASET_PATH / "{file}.tif", file=config["datasets"])
    params:
        outpath=DATASET_PATH,
        inpath=BASE_INPUT_PATH,
        files_convert=config["datasets"]
    shell:
        """
        {{
            for file in {params.files_convert}; do
                infile="{params.inpath}/$file.tif"
                outfile="{params.outpath}/$file.tif"
                gdal_translate "$infile" "$outfile" -ot Float32
            done
        }}
        """
