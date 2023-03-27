"""
Snakefile for this dataset

Note: To include the file in the main workflow, it must be added to the respective module.smk .
"""

from digipipe.store.utils import (
    get_abs_dataset_path,
    PATH_TO_REGION_MUNICIPALITIES_GPKG,
    PATH_TO_REGION_DISTRICTS_GPKG
)

DATASET_PATH = get_abs_dataset_path("datasets", "bnetza_mastr_hydro_region")
SOURCE_DATASET_PATH = get_abs_dataset_path("preprocessed", "bnetza_mastr", data_dir=True)

rule create:
    """
    Extract hydro units for region
    """
    input:
        units=SOURCE_DATASET_PATH / "bnetza_mastr_hydro_raw.csv",
        locations=SOURCE_DATASET_PATH / "bnetza_mastr_locations_extended_raw.csv",
        gridconn=SOURCE_DATASET_PATH / "bnetza_mastr_grid_connections_raw.csv",
        region_muns=PATH_TO_REGION_MUNICIPALITIES_GPKG,
        region_districts=PATH_TO_REGION_DISTRICTS_GPKG
    output:
        outfile=DATASET_PATH / "data" / "bnetza_mastr_hydro_region.gpkg",
        outfile_agg=DATASET_PATH / "data" / "bnetza_mastr_hydro_agg_region.gpkg"
    params:
        config_file=DATASET_PATH / "config.yml"
    log: DATASET_PATH / "data" / "bnetza_mastr_hydro_region.log"
    script:
        DATASET_PATH / "scripts" / "create.py"
