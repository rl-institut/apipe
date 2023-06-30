"""
Snakefile for this dataset

Note: To include the file in the main workflow, it must be added to the respective module.smk .
"""

import json
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
        pv_roof_units=rules.datasets_bnetza_mastr_pv_roof_region_create.output.outfile,
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
    Create stats on installed count of units and capacity per mun
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

rule create_storage_pv_roof_stats:
    input:
        units=DATASET_PATH / "data" / "bnetza_mastr_storage_region.gpkg",
        pv_roof_units=rules.datasets_bnetza_mastr_pv_roof_region_create.output.outfile,
    output: DATASET_PATH / "data" / "bnetza_mastr_storage_pv_roof.json"
    run:
        # PV home units: guess PV home systems
        pv_roof_units = gpd.read_file(input.pv_roof_units)[
            ["mastr_id", "capacity_net"]
        ]
        hs_cfg = config.get("home_storages")
        mask = (
            (pv_roof_units.capacity_net >=
             hs_cfg.get("pv_roof_capacity_thres_min")) &
            (pv_roof_units.capacity_net <=
             hs_cfg.get("pv_roof_capacity_thres_max"))
        )
        pv_roof_units_small = pv_roof_units.loc[mask]

        # Storage units: calc specific values
        units = gpd.read_file(input.units)[[
            "capacity_net",
            "storage_capacity",
            "mastr_location_id",
            "pv_roof_unit_count",
            "pv_roof_unit_capacity_sum"
        ]]
        units_unique_loc = units.groupby("mastr_location_id").agg(
            storage_count=("storage_capacity", "count"),
            storage_capacity=("storage_capacity", "sum"),
            capacity_net=("capacity_net", "mean"),
            pv_roof_unit_count=("pv_roof_unit_count", "first"),
            pv_roof_unit_capacity_sum=("pv_roof_unit_capacity_sum", "first"),
        )

        # Filter storages
        mask = (
            (units_unique_loc.storage_count == 1
             if hs_cfg.get("only_single_storages")
             else units_unique_loc.storage_count > 0
             ) &
            (units_unique_loc.storage_capacity <=
             hs_cfg.get("storage_capacity_thres_max")) &
            (units_unique_loc.capacity_net <=
             hs_cfg.get("storage_power_thres_max")) &
            (units_unique_loc.pv_roof_unit_count == 1
             if hs_cfg.get("only_single_pv_roof_units")
             else units_unique_loc.pv_roof_unit_count > 0
             ) &
            (units_unique_loc.pv_roof_unit_capacity_sum >=
             hs_cfg.get("pv_roof_capacity_thres_min")) &
            (units_unique_loc.pv_roof_unit_capacity_sum <=
             hs_cfg.get("pv_roof_capacity_thres_max"))
        )
        units_unique_loc_small = units_unique_loc.loc[mask]

        with open(output[0], "w", encoding="utf8") as f:
            json.dump(
                {
                    "pv_roof_share": {
                        "all_storages": round(
                            len(units_unique_loc) /
                            len(pv_roof_units),
                            3,
                        ),
                        "home_storages": round(
                            len(units_unique_loc_small) /
                            len(pv_roof_units_small),
                            3,
                        ),
                    },
                    "specific_capacity": {
                        "all_storages": round(
                            units_unique_loc.storage_capacity.sum()
                            / units_unique_loc.pv_roof_unit_capacity_sum.sum(),
                            2,
                        ),
                        "home_storages": round(
                            units_unique_loc_small.storage_capacity.sum()
                            / units_unique_loc_small.pv_roof_unit_capacity_sum.sum(),
                            2,
                        ),
                    },
                    "specific_power": {
                        "all_storages": round(
                            units_unique_loc.capacity_net.sum()
                            / units_unique_loc.pv_roof_unit_capacity_sum.sum(),
                            2,
                        ),
                        "home_storages": round(
                            units_unique_loc_small.capacity_net.sum()
                            / units_unique_loc_small.pv_roof_unit_capacity_sum.sum(),
                            2,
                        ),
                    },
                },
                f,
                indent=4,
            )
