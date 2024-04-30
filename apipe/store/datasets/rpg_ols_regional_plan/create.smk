"""
Snakefile for this dataset

Note: To include the file in the main workflow, it must be added to the respective module.smk .
"""

import geopandas as gpd
import pandas as pd

from apipe.scripts.datasets.mastr import create_stats_per_municipality
from apipe.scripts.geo import (
    convert_to_multipolygon,
    overlay,
    write_geofile,
)
from apipe.store.utils import (
    get_abs_dataset_path,
    PATH_TO_REGION_MUNICIPALITIES_GPKG,
)

DATASET_PATH = get_abs_dataset_path("datasets", "rpg_ols_regional_plan")


rule create_pv_ground_criteria:
    """
    Freiflächen-Photovoltaikanlagen Negativkriterien
    """
    input:
        get_abs_dataset_path(
            "preprocessed", "rpg_ols_regional_plan"
        ) / "data" / "{file}.gpkg"
    output:
        DATASET_PATH / "data" / "{file}.gpkg"
    shell: "cp -p {input} {output}"


rule create_pv_ground_units_filtered:
    """
    Filter PV units for different status and add municipality ids
    """
    input:
        units=DATASET_PATH / "data" / "rpg_ols_pv_ground.gpkg",
        region_muns=PATH_TO_REGION_MUNICIPALITIES_GPKG,
    output:
        units_all=DATASET_PATH / "data" / "rpg_ols_pv_ground_all.gpkg",
        units_filtered=expand(
            DATASET_PATH / "data" / "rpg_ols_pv_ground_{status}.gpkg",
            status=["operating", "approved", "planned"]
        )
    run:
        units = gpd.read_file(input.units)

        # assign mun id
        units = overlay(
            gdf=units,
            gdf_overlay=gpd.read_file(input.region_muns),
            retain_rename_overlay_columns={"id": "municipality_id"},
        )

        # Write all
        write_geofile(
            gdf=convert_to_multipolygon(units),
            file=DATASET_PATH / "data" / "rpg_ols_pv_ground_all.gpkg"
        )

        # Write filtered
        for status, file_suffix in {
            "realisiert": "operating",
            "Planung": "planned",
            "genehmigt": "approved",
        }.items():
            write_geofile(
                gdf=convert_to_multipolygon(units.loc[units.status == status].copy()),
                file=DATASET_PATH / "data" / f"rpg_ols_pv_ground_{file_suffix}.gpkg"
            )


rule create_pv_ground_power_stats_muns:
    """
    Create stats on installed count of units and power per mun
    """
    input:
        units=expand(
            DATASET_PATH / "data" / "rpg_ols_pv_ground_{status}.gpkg",
            status=["all", "operating", "approved", "planned"]
        ),
        region_muns=PATH_TO_REGION_MUNICIPALITIES_GPKG,
    output:
        units_filtered=expand(
            DATASET_PATH/ "data" / "rpg_ols_pv_ground_stats_muns_{status}.csv",
            status=["all", "operating", "approved", "planned"]
        )
    run:
        for status in ["all", "operating", "approved", "planned"]:
            units = gpd.read_file(
                DATASET_PATH / "data" / f"rpg_ols_pv_ground_{status}.gpkg"
            )
            units = create_stats_per_municipality(
                units_df=units,
                muns=gpd.read_file(input.region_muns),
                column="capacity_net",
                only_operating_units=False  # Disable MaStR-specific setting
            )
            print(f"Total capacity for {status} units: {units.capacity_net.sum()}")
            units.to_csv(
                DATASET_PATH / "data" /
                f"rpg_ols_pv_ground_stats_muns_{status}.csv"
            )
