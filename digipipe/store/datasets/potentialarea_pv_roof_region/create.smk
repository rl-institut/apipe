"""
Snakefile for this dataset

Note: To include the file in the main workflow, it must be added to the respective module.smk .
"""
import json
import re
import geopandas as gpd
import pandas as pd
from pathlib import Path
from digipipe.scripts.geo import (
    overlay,
    convert_to_multipolygon,
    write_geofile
)
from digipipe.store.utils import (
    get_abs_dataset_path,
    PATH_TO_REGION_MUNICIPALITIES_GPKG
)

DATASET_PATH = get_abs_dataset_path(
    "datasets", "potentialarea_pv_roof_region", data_dir=True
)

rule overlay_muns:
    """
    Overlay potential area with municipalities
    """
    input:
        area=get_abs_dataset_path(
            "preprocessed", "rpg_abw_pv_roof_potential", data_dir=True
        ) / "rpg_abw_pv_roof_potential.gpkg",
        region_muns=PATH_TO_REGION_MUNICIPALITIES_GPKG
    output:
        area=DATASET_PATH / "potentialarea_pv_roof_region.gpkg"
    run:
        data = gpd.read_file(input.area)
        data = data.assign(
            historic_preservation=data.historic_preservation.fillna(
                False).replace({"true": True})
        )
        data = overlay(
            gdf=data,
            gdf_overlay=gpd.read_file(input.region_muns),
            retain_rename_overlay_columns={"id": "municipality_id"},
            gdf_use_centroid=True
        )
        write_geofile(
            gdf=convert_to_multipolygon(data),
            file=output.area,
        )

rule create_area_stats_muns:
    """
    Create stats on PV roof potentials per mun for 1) all and 2) all but
    historic buildings.
    """
    input:
        area=DATASET_PATH / "potentialarea_pv_roof_region.gpkg",
        region_muns=PATH_TO_REGION_MUNICIPALITIES_GPKG,
    output:
        total=DATASET_PATH / "potentialarea_pv_roof_area_stats_muns.csv",
        wo_historic=(
            DATASET_PATH /
            "potentialarea_pv_roof_wo_historic_area_stats_muns.csv"
        )
    run:
        print("PV roof potential area stats:")
        muns = gpd.read_file(input.region_muns)
        potential_all = gpd.read_file(input.area).fillna(0)

        # define columns
        cols_base = [
            "municipality_id", "historic_preservation", "building_area_sqm"
        ]
        orientation_suffix = ["south", "north", "east", "west", "flat"]
        cols_power = [
            f"installable_power_kw_{orient}"
            for orient in orientation_suffix
        ]
        cols_power_new = {
            c.replace("_kw_", "_"): (c, "sum") for c in cols_power
        }
        cols_energy = [
            f"energy_annual_mwh_{orient}"
            for orient in orientation_suffix
        ]
        cols_energy_new = {
            c.replace("_mwh_","_"): (c, "sum") for c in cols_energy
        }

        # aggregate per mun
        agg_cols = dict(
            roof_count=("building_area_sqm", "count"),
            building_area_sqm=("building_area_sqm", "sum"),
            historic_preservation_count=("historic_preservation", "sum"),
            **cols_power_new,
            **cols_energy_new
        )
        potential_wo_historic = potential_all.copy()
        potential_wo_historic = potential_wo_historic.loc[
            ~potential_wo_historic.historic_preservation
        ]

        for df, file in zip(
            [potential_all, potential_wo_historic],
            [output.total, output.wo_historic]
        ):
            df = df[
                cols_base + cols_power + cols_energy
            ].groupby("municipality_id").agg(**agg_cols)
            # kW -> MW
            df[[_ for _ in cols_power_new.keys()]] = df[
                cols_power_new.keys()].div(1e3)
            # Produce totals
            df = df.assign(
                installable_power_total=df[
                    [_ for _ in cols_power_new.keys()]].sum(axis=1),
                energy_annual_total=df[
                    [_ for _ in cols_energy_new.keys()]].sum(axis=1),
            )
            # Dump
            df.to_csv(file)

rule create_relative_deployment_stats_muns:
    """
    Create stats on how much of the theoretically installable PV rooftop
    potential is used per mun for 1) all and 2) all but historic buildings.
    """
    input:
        area_stats_total=(
            rules.datasets_potentialarea_pv_roof_region_create_area_stats_muns.output.total
        ),
        area_stats_wo_historic=(
            rules.datasets_potentialarea_pv_roof_region_create_area_stats_muns.output.wo_historic
        ),
        unit_stats=(
            rules.datasets_bnetza_mastr_pv_roof_region_create_power_stats_muns.output[0]
        )
    output:
        total=DATASET_PATH / "potentialarea_pv_roof_deployment_stats_muns.csv",
        wo_historic=(
            DATASET_PATH/
            "potentialarea_pv_roof_wo_historic_deployment_stats_muns.csv"
        )
    run:
        orientation_suffix = ["south", "north", "east", "west", "flat"]
        cols_power = [
            f"installable_power_{orient}"
            for orient in orientation_suffix
        ]
        pv_installed = pd.read_csv(
            input.unit_stats,
            usecols=["municipality_id", "capacity_net"],
            index_col="municipality_id",
        )

        for file_in, file_out in zip(
            [input.area_stats_total, input.area_stats_wo_historic],
            [output.total, output.wo_historic]
        ):
            pv_potential = pd.read_csv(
                file_in,
                usecols=["municipality_id"] + cols_power,
                index_col="municipality_id",
            )
            pv_potential = pd.DataFrame(
                pv_potential.sum(axis=1),columns=["installable_power"]
            )
            pv_deployed = pd.concat(
                [pv_installed.capacity_net, pv_potential.installable_power],
                axis=1,
            )
            pv_deployed = pv_deployed.assign(
                relative_deployment=(
                    pv_installed.capacity_net.div(
                        pv_potential.installable_power
                    )
                )
            )
            pv_deployed.to_csv(file_out)
