"""
Snakefile for this dataset

Note: To include the file in the main workflow, it must be added to the respective module.smk .
"""

from digipipe.store.utils import (
    get_abs_dataset_path,
    PATH_TO_ABW_MUNICIPALITIES_GPKG,
    PATH_TO_ABW_DISTRICTS_GPKG
)

DATASET_PATH = get_abs_dataset_path("datasets", "bnetza_mastr_biomass_abw")
SOURCE_DATASET_PATH = get_abs_dataset_path("preprocessed", "bnetza_mastr", data_dir=True)

rule create:
    """
    Extract biomass units for ABW region
    """
    input:
        units=SOURCE_DATASET_PATH / "bnetza_mastr_biomass_raw.csv",
        locations=SOURCE_DATASET_PATH / "bnetza_mastr_locations_extended_raw.csv",
        gridconn=SOURCE_DATASET_PATH / "bnetza_mastr_grid_connections_raw.csv",
        abw_muns=PATH_TO_ABW_MUNICIPALITIES_GPKG,
        abw_districts=PATH_TO_ABW_DISTRICTS_GPKG
    output:
        outfile=DATASET_PATH / "data" / "bnetza_mastr_biomass_abw.gpkg",
        outfile_agg=DATASET_PATH / "data" / "bnetza_mastr_biomass_agg_abw.gpkg"
    params:
        config_file=DATASET_PATH / "config.yml"
    script:
        DATASET_PATH / "scripts" / "create.py"
