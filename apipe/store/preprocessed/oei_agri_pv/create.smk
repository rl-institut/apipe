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

rule convert_to_gpkg:
    input:
        expand(BASE_INPUT_PATH / "{file}.tif", file=config["files_convert"])
    output:
        expand(DATASET_PATH / "{file}.gpkg", file=config["files_convert"])
    params:
        outpath=DATASET_PATH,
        inpath=BASE_INPUT_PATH,
        files_convert=config["files_convert"]
    shell:
        """
        {{
            for file in {params.files_convert}; do
                infile="{params.inpath}/$file.tif"
                outfile="{params.outpath}/$file.gpkg"
                gdal_translate "$infile" "$outfile" -ot Float32 -of GPKG
            done
        }}
        """
