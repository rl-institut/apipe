"""
Snakefile for this dataset

Note: To include the file in the main workflow, it must be added to the respective module.smk .
"""

from digipipe.store.utils import get_abs_dataset_path
#from digipipe.store.datasets.bnetza_mastr_wind_abw.scripts.create import process

DATASET_PATH = get_abs_dataset_path("datasets", "bnetza_mastr_wind_abw")
SOURCE_DATASET_PATH = get_abs_dataset_path("preprocessed", "bnetza_mastr", data_dir=True)

rule create:
    """
    Extract districts of ABW region
    """
    input:
        units=SOURCE_DATASET_PATH / "bnetza_mastr_wind_raw.csv",
        locations=SOURCE_DATASET_PATH / "bnetza_mastr_locations_extended_raw.csv",
        gridconn=SOURCE_DATASET_PATH / "bnetza_mastr_grid_connections_raw.csv",
        abw_districts=get_abs_dataset_path("datasets", "bkg_vg250_districts_abw", data_dir=True) / "bkg_vg250_districts_abw.gpkg",
        abw_muns=get_abs_dataset_path("datasets", "bkg_vg250_muns_abw", data_dir=True) / "bkg_vg250_muns_abw.gpkg"
    output:
        outfile=DATASET_PATH / "data" / "bnetza_mastr_wind_abw.gpkg"
    params:
        config_file=DATASET_PATH / "config.yml"
    script:
        DATASET_PATH / "scripts" / "create.py"
