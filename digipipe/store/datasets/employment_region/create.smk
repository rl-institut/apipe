"""
Snakefile for this dataset

Note: To include the file in the main workflow, it must be added to the respective module.smk .
"""

from digipipe.store.utils import (
    get_abs_dataset_path,
    PATH_TO_REGION_MUNICIPALITIES_GPKG
)

DATASET_PATH = get_abs_dataset_path("datasets", "employment_region")

rule create:
    """
    Create employment dataset for region
    """
    input:
        employment=get_abs_dataset_path("preprocessed", "ba_employment") / "data" / "employees.csv",
        region_muns=PATH_TO_REGION_MUNICIPALITIES_GPKG
    output:
        DATASET_PATH / "data" / "employees.csv"
    script:
        DATASET_PATH / "scripts" / "create.py"
