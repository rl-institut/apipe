"""
Snakefile for this dataset

Note: To include the file in the main workflow, it must be added to the respective module.smk .
"""

import geopandas as gpd

from digipipe.scripts.datasets.mastr import create_stats_per_municipality
from digipipe.store.utils import (
    get_abs_dataset_path,
    PATH_TO_REGION_MUNICIPALITIES_GPKG,
    PATH_TO_REGION_DISTRICTS_GPKG
)

DATASET_PATH = get_abs_dataset_path("datasets", "bnetza_mastr_storage_region")
SOURCE_DATASET_PATH = get_abs_dataset_path("preprocessed", "bnetza_mastr", data_dir=True)

rule create:
    """
    Extract storages for region
    """
    input:
        units=SOURCE_DATASET_PATH / "bnetza_mastr_storage_raw.csv",
        units_capacity=SOURCE_DATASET_PATH / "bnetza_mastr_storage_unit_raw.csv",
        locations=SOURCE_DATASET_PATH / "bnetza_mastr_locations_extended_raw.csv",
        gridconn=SOURCE_DATASET_PATH / "bnetza_mastr_grid_connections_raw.csv",
        region_muns=PATH_TO_REGION_MUNICIPALITIES_GPKG,
        region_districts=PATH_TO_REGION_DISTRICTS_GPKG
    output:
        outfile=DATASET_PATH / "data" / "bnetza_mastr_storage_region.gpkg",
        outfile_agg=DATASET_PATH / "data" / "bnetza_mastr_storage_agg_region.gpkg"
    params:
        config_file=DATASET_PATH / "config.yml"
    script:
        DATASET_PATH / "scripts" / "create.py"

rule create_power_stats_muns:
    """
    Create stats on installed count of units and power per mun
    """
    input:
        units=DATASET_PATH / "data" / "bnetza_mastr_storage_region.gpkg",
        region_muns=PATH_TO_REGION_MUNICIPALITIES_GPKG,
    output:
        total=DATASET_PATH / "data" / "bnetza_mastr_storage_stats_muns.csv",
        large=DATASET_PATH/ "data" / "bnetza_mastr_storage_large_stats_muns.csv",
        small=DATASET_PATH/ "data" / "bnetza_mastr_storage_small_stats_muns.csv",
    run:
        units = gpd.read_file(input.units)
        units["storage_capacity"] = units["storage_capacity"].div(1e3)  # kW to MW
        muns = gpd.read_file(input.region_muns)
        print("Battery storages:")

        # All storage units
        units_total = create_stats_per_municipality(
            units_df=units,
            muns=muns,
            column="storage_capacity"
        )
        print(
            f"  Storage capacity (total): "
            f"{units_total.storage_capacity.sum().round(3)} MWh"
        )
        units_total.to_csv(output.total)

        # Large storage units
        units_large = create_stats_per_municipality(
            units_df=units.loc[
                units.storage_capacity >= config.get("battery_size_threshold")
            ],
            muns=muns,
            column="storage_capacity"
        )
        print(
            f"  Storage capacity (large): "
            f"{units_large.storage_capacity.sum().round(3)} MWh"
        )
        units_large.to_csv(output.large)

        # Small storage units
        units_small = create_stats_per_municipality(
            units_df=units.loc[
                units.storage_capacity < config.get("battery_size_threshold")
            ],
            muns=muns,
            column="storage_capacity"
        )
        print(
            f"  Storage capacity (small): "
            f"{units_small.storage_capacity.sum().round(3)} MWh"
        )
        units_small.to_csv(output.small)
