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
        employment_total=get_abs_dataset_path("preprocessed", "ba_employment") / "data" / "employment_muns.csv",
        employment_ind=get_abs_dataset_path("preprocessed", "regiostat") / "data" / "employment_industry_muns.csv",
        region_muns=PATH_TO_REGION_MUNICIPALITIES_GPKG
    output:
        DATASET_PATH / "data" / "employment.csv"
    script:
        DATASET_PATH / "scripts" / "create.py"
