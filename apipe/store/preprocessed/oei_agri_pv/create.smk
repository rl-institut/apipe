"""
Snakefile for this dataset

Note: To include the file in the main workflow, it must be added to the respective module.smk .
"""

from apipe.store.utils import get_abs_dataset_path

DATASET_PATH = get_abs_dataset_path(
    "preprocessed", "oei_agri_pv", data_dir=True
)

rule convert_to_gpkg:
    input:
        get_abs_dataset_path("raw", "oei_agri_pv")
        / "data"
        / "Agri-PV-Potenziale_SQR_50-70_100x100_EPSG3035_pos.tif",
    output:
        DATASET_PATH / "data" / "oei_agri_pv.gpkg",
    shell:
        "gdal_translate {input} {output} -of GPKG"
