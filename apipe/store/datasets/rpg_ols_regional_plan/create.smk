"""
Snakefile for this dataset

Note: To include the file in the main workflow, it must be added to the respective module.smk .
"""

import geopandas as gpd
import numpy as np
import pandas as pd

from apipe.scripts.data_io import load_json
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

DATASET_PATH = get_abs_dataset_path(
    "datasets", "rpg_ols_regional_plan", data_dir=True)


rule create_pv_ground_criteria_single:
    """
    Freiflächen-Photovoltaikanlagen Negativkriterien einzeln
    """
    input:
        get_abs_dataset_path(
            "preprocessed", "rpg_ols_regional_plan"
        ) / "data" / "{file}.gpkg"
    output:
        DATASET_PATH / "{file}.gpkg"
    shell: "cp -p {input} {output}"


rule create_pv_ground_units_filtered:
    """
    Filter PV units for different status and add municipality ids.
    For approved and planned units add missing power (=0) using power density
    from tech data.
    """
    input:
        units=DATASET_PATH / "rpg_ols_pv_ground.gpkg",
        region_muns=PATH_TO_REGION_MUNICIPALITIES_GPKG,
        tech_data=rules.datasets_technology_data_copy_files.output,
    output:
        units_all=DATASET_PATH / "rpg_ols_pv_ground_all.gpkg",
        units_filtered=expand(
            DATASET_PATH / "rpg_ols_pv_ground_{status}.gpkg",
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

        # Add capacity to units where capacity is 0
        # (approximation for planned units)
        tech_data = load_json(input.tech_data[0])
        print(
            "Capacity in original data: ",
            units[["status", "capacity_net"]].groupby(
                "status").capacity_net.sum()
        )
        units = units.assign(capacity_net_inferred=0)
        mask = (
            (units.status.isin(["Planung", "genehmigt"])) &
            (units.capacity_net == 0)
        )
        units["capacity_net"].update(
            units[mask].area / 1e6 *
            tech_data["power_density"]["pv_ground"]
        )
        units.loc[mask, "capacity_net_inferred"] = 1
        units["capacity_net"] = units["capacity_net"].mul(1e3)  # MW to kW
        print(
            "Capacity inferred: ",
            units[["status", "capacity_net"]].groupby(
                "status").capacity_net.sum()
        )

        # Write all
        write_geofile(
            gdf=convert_to_multipolygon(units),
            file=output.units_all
        )

        # Write filtered
        for status, file_suffix in {
            "realisiert": "operating",
            "Planung": "planned",
            "genehmigt": "approved",
        }.items():
            write_geofile(
                gdf=convert_to_multipolygon(units.loc[units.status == status].copy()),
                file=DATASET_PATH / f"rpg_ols_pv_ground_{file_suffix}.gpkg"
            )


rule create_pv_ground_power_stats_muns:
    """
    Create stats on installed count of units and power per mun
    """
    input:
        units=expand(
            DATASET_PATH / "rpg_ols_pv_ground_{status}.gpkg",
            status=["all", "operating", "approved", "planned"]
        ),
        region_muns=PATH_TO_REGION_MUNICIPALITIES_GPKG,
    output:
        stats=expand(
            DATASET_PATH / "rpg_ols_pv_ground_stats_muns_{status}.csv",
            status=["all", "operating", "approved", "planned"]
        )
    run:
        for status in ["all", "operating", "approved", "planned"]:
            units = gpd.read_file(
                DATASET_PATH / f"rpg_ols_pv_ground_{status}.gpkg"
            )
            units = create_stats_per_municipality(
                units_df=units,
                muns=gpd.read_file(input.region_muns),
                column="capacity_net",
                only_operating_units=False  # Disable MaStR-specific setting
            )
            units["capacity_net"] = units["capacity_net"].div(1e3)  # kW to MW
            print(f"Total capacity for {status} units: {units.capacity_net.sum()}")
            units.to_csv(
                DATASET_PATH /
                f"rpg_ols_pv_ground_stats_muns_{status}.csv"
            )


rule create_wind_units:
    """
    Ddd municipality ids to wind units
    """
    input:
        units=expand(
            get_abs_dataset_path("preprocessed", "rpg_ols_regional_plan") /
            "data" / "rpg_ols_wind_{status}.gpkg",
            status=["approved", "planned", "operating"]
        ),
        region_muns=PATH_TO_REGION_MUNICIPALITIES_GPKG
    output:
        units=expand(
            DATASET_PATH / "rpg_ols_wind_{status}.gpkg",
            status=["approved", "planned", "operating"]
        ),
        units_all=DATASET_PATH / "rpg_ols_wind_all.gpkg"
    run:
        all_units = list()
        for file_in, file_out, status in zip(
                input.units, output.units, ["approved", "planned", "operating"]
        ):
            units = gpd.read_file(file_in)
            # assign mun id
            units = overlay(gdf=units,
                gdf_overlay=gpd.read_file(input.region_muns),
                retain_rename_overlay_columns={"id": "municipality_id"},
            )
            all_units.append(units.copy().assign(status=status))
            write_geofile(
                gdf=units,
                file=file_out
            )
        write_geofile(
            gdf=pd.concat(all_units, axis=0),
            file=output.units_all
        )


rule create_wind_power_stats_muns:
    """
    Create stats on installed count of units and power per mun
    """
    input:
        units=expand(
            DATASET_PATH / "rpg_ols_wind_{status}.gpkg",
            status=["all", "approved", "planned", "operating"]
        ),
        region_muns=PATH_TO_REGION_MUNICIPALITIES_GPKG,
    output:
        stats=expand(
            DATASET_PATH / "rpg_ols_wind_stats_muns_{status}.csv",
            status=["all", "operating", "approved", "planned"]
        )
    run:
        for status in ["all", "operating", "approved", "planned"]:
            units = gpd.read_file(
                DATASET_PATH / f"rpg_ols_wind_{status}.gpkg"
            )
            units = create_stats_per_municipality(
                units_df=units,
                muns=gpd.read_file(input.region_muns),
                column="capacity_net",
                only_operating_units=False  # Disable MaStR-specific setting
            )
            print(f"Total capacity for {status} units: {units.capacity_net.sum()}")
            units.to_csv(
                DATASET_PATH /
                f"rpg_ols_wind_stats_muns_{status}.csv"
            )
