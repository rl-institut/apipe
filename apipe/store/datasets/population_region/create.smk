"""
Snakefile for this dataset

Note: To include the file in the main workflow, it must be added to the respective module.smk .
"""

from apipe.store.utils import (
    get_abs_dataset_path,
    PATH_TO_REGION_MUNICIPALITIES_GPKG,
    PATH_TO_REGION_DISTRICTS_GPKG
)

DATASET_PATH = get_abs_dataset_path("datasets", "population_region")

rule create:
    """
    Create full population dataset for region
    """
    input:
        pop_history=expand(
            get_abs_dataset_path("preprocessed", "destatis_gv") / "data" / "3112{year}_Auszug_GV.csv",
            year=[2010, 2015, 2020, 2021, 2022, 2023]
        ),
        # prognosis_fstate_munlevel=get_abs_dataset_path(
        #     "preprocessed", "stala_st_pop_prog") / "data" /
        #     "population_prognosis_st_muns.csv",
        # prognosis_germany_districtlevel=get_abs_dataset_path(
        #     "preprocessed", "demandregio") / "data" / "dr_hh_population.csv",
        prognosis_germany_fstatelevel=get_abs_dataset_path(
            "preprocessed", "destatis_pop_prog") / "data" / "population_prognosis_federal_states.csv",
        region_muns=PATH_TO_REGION_MUNICIPALITIES_GPKG,
        region_districts=PATH_TO_REGION_DISTRICTS_GPKG,
        federal_states=get_abs_dataset_path(
            "datasets", "bkg_vg250_federal_states") / "data" / "bkg_vg250_federal_states.gpkg"
    output:
        DATASET_PATH / "data" / "population.csv"
    script:
        DATASET_PATH / "scripts" / "create.py"
