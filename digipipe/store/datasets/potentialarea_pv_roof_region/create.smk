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
    Create stats on PV potential areas per mun
    """
    input:
        area=DATASET_PATH / "potentialarea_pv_roof_region.gpkg",
        region_muns=PATH_TO_REGION_MUNICIPALITIES_GPKG,
    output: DATASET_PATH / "potentialarea_pv_roof_area_stats_muns.csv"
    run:
        print("PV roof potential area stats:")
        muns = gpd.read_file(input.region_muns)
        data = gpd.read_file(input.area).fillna(0)

        # define columns
        cols_base = ["municipality_id", "historic_preservation", "building_area_sqm"]
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
        agg_cols=dict(
            roof_count=("building_area_sqm", "count"),
            building_area_sqm=("building_area_sqm", "sum"),
            historic_preservation_count=("historic_preservation", "sum"),
            **cols_power_new,
            **cols_energy_new
        )
        data = data[cols_base + cols_power + cols_energy].groupby(
            "municipality_id").agg(**agg_cols)

        # kW -> MW
        data[[_ for _ in cols_power_new.keys()]] = data[
            cols_power_new.keys()].div(1e3)

        # Produce totals
        data = data.assign(
            installable_power_total=data[
                [_ for _ in cols_power_new.keys()]].sum(axis=1),
            energy_annual_total=data[
                [_ for _ in cols_energy_new.keys()]].sum(axis=1),
        )

        # Dump
        data.to_csv(output[0])

rule create_relative_deployment_stats_muns:
    """
    Create stats on how much of the theoretically installable PV rooftop
    potential is used, per mun.
    """
    input:
        area_stats=(
            rules.datasets_potentialarea_pv_roof_region_create_area_stats_muns.output[0]
        ),
        pv_roof_stats=(
            rules.datasets_bnetza_mastr_pv_roof_region_create_power_stats_muns.output[0]
        )
    output: DATASET_PATH / "potentialarea_pv_roof_deployment_stats_muns.csv"
    run:
        orientation_suffix = ["south", "north", "east", "west", "flat"]
        cols_power = [
            f"installable_power_{orient}"
            for orient in orientation_suffix
        ]
        pv_potential = pd.read_csv(
            input.area_stats,
            usecols=["municipality_id"]+cols_power,
            index_col="municipality_id",
        )
        pv_potential = pd.DataFrame(
            pv_potential.sum(axis=1), columns=["installable_power"]
        )

        pv_installed = pd.read_csv(
            input.pv_roof_stats,
            usecols=["municipality_id", "capacity_net"],
            index_col="municipality_id",
        )

        pv_deployed = pd.DataFrame(
            pv_installed.capacity_net.div(pv_potential.installable_power),
            columns=["relative_deployment"],
        )

        # Dump
        pv_deployed.to_csv(output[0])
