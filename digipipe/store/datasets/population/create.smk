"""
Snakefile for this dataset

Note: To include the file in the main workflow, it must be added to the respective module.smk .
"""

from digipipe.store.utils import (
    get_abs_dataset_path,
    get_abs_store_root_path,
    PATH_TO_REGION_MUNICIPALITIES_GPKG
)

DATASET_PATH = get_abs_dataset_path("datasets", "population")
STORE_PATH = get_abs_store_root_path()

rule create:
    """
    Create full population dataset for region
    """
    input:
        pop_history=expand(
            get_abs_dataset_path("preprocessed", "destatis_gv") / "data" / "3112{year}_Auszug_GV.csv",
            year=[2010, 2015, 2020, 2021, 2022]
        ),
        pop_prognosis=get_abs_dataset_path("preprocessed", "stala_st_pop_prog") / "data" / "population_prognosis.csv",
        region_muns=PATH_TO_REGION_MUNICIPALITIES_GPKG
    output:
        DATASET_PATH / "data" / "population.csv"
    log:
        STORE_PATH / "datasets" / ".log" / "population.log"
    script:
        DATASET_PATH / "scripts" / "create.py"
