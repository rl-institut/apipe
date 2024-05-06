"""
Snakefile for this dataset

Note: To include the file in the main workflow, it must be added to the respective module.smk .
"""
import json
import re
import geopandas as gpd
import pandas as pd
from apipe.config import GLOBAL_CONFIG
from apipe.scripts.data_io import load_json
from apipe.scripts.geo import (
    overlay,
    convert_to_multipolygon,
    write_geofile
)
from apipe.store.utils import (
    get_abs_dataset_path,
    PATH_TO_REGION_MUNICIPALITIES_GPKG
)

DATASET_PATH = get_abs_dataset_path(
    "datasets", "potentialarea_pv_roof_region2", data_dir=True
)

rule clip_to_region_and_filter:
    """
    Clip to region and filter roof suitability
    """
    input:
        geodata=get_abs_dataset_path(
            "preprocessed", "wfbb_pv_roof_potential", data_dir=True
        ) / "solaratlas_eignung_dachflaechen_pv.gpkg",
        region=rules.datasets_bkg_vg250_region_create.output
    output:
        geodata=DATASET_PATH / "solaratlas_eignung_dachflaechen_pv_region.gpkg"
    run:
        print("Reading file...")
        geodata = gpd.read_file(input.geodata).to_crs(
            GLOBAL_CONFIG["global"]["geodata"]["crs"]
        )

        print("Overlaying...")
        geodata = overlay(
            gdf=geodata,
            gdf_overlay=gpd.read_file(input.region[0]),
        )

        # Delete non-suitable roofs
        print("Filtering...")
        geodata = geodata.loc[
            geodata["roof_suitability_perc"] >=
            config.get("roof_suitability_threshold") * 100
        ]

        print("Writing result file...")
        write_geofile(
            gdf=convert_to_multipolygon(geodata),
            file=output.geodata,
            layer_name="solaratlas_eignung_dachflaechen_pv_region",
        )

rule overlay_muns:
    """
    Overlay potential area with municipalities
    """
    input:
        area=DATASET_PATH / "solaratlas_eignung_dachflaechen_pv_region.gpkg",
        region_muns=PATH_TO_REGION_MUNICIPALITIES_GPKG
    output:
        area=DATASET_PATH / "potentialarea_pv_roof_region.gpkg"
    run:
        print("Reading file...")
        data = gpd.read_file(input.area)

        print("Overlaying...")
        data = overlay(
            gdf=data,
            gdf_overlay=gpd.read_file(input.region_muns),
            retain_rename_overlay_columns={"id": "municipality_id"},
            gdf_use_centroid=True
        )
        print("Writing result file...")
        write_geofile(
            gdf=data,
            file=output.area,
        )

rule create_area_stats_muns:
    """
    Create stats on PV roof potentials per mun.
    """
    input:
        area=DATASET_PATH / "potentialarea_pv_roof_region.gpkg",
        region_muns=PATH_TO_REGION_MUNICIPALITIES_GPKG
    output:
        stats=DATASET_PATH / "potentialarea_pv_roof_area_stats_muns.csv"
    run:
        muns = gpd.read_file(input.region_muns)
        potential = gpd.read_file(input.area).fillna(0)

        # define columns
        cols_base = [
            "municipality_id", "roof_area_sqm", "roof_area_pv_potential_sqm",
            "roof_suitability_perc", "installable_power_kw", "energy_annual_mwh"
        ]

        # aggregate per mun
        agg_cols = dict(
            roof_count=("roof_area_sqm", "count"),
            roof_area_sqkm=("roof_area_sqm", "sum"),
            roof_area_pv_potential_sqkm=("roof_area_pv_potential_sqm", "sum"),
            roof_suitability_perc=("roof_suitability_perc", "mean"),
            installable_power=("installable_power_kw", "sum"),
            energy_annual=("energy_annual_mwh", "sum"),
        )

        potential = potential[
            cols_base].groupby("municipality_id").agg(**agg_cols)

        # kW -> MW, sqm -> sqkm
        potential = potential.assign(
            installable_power=potential["installable_power"].div(1e3),
            roof_area_sqkm=potential["roof_area_sqkm"].div(1e6),
            roof_area_pv_potential_sqkm=potential["roof_area_pv_potential_sqkm"].div(1e6),
        )
        print(f"Total installable power: {potential.installable_power.sum()} MW")
        print(f"Total energy annual: {potential.energy_annual.sum()} MWh")

        # Dump
        potential.to_csv(output.stats)

rule create_relative_deployment_stats_muns:
    """
    Create stats on how much of the theoretically installable PV rooftop
    potential is used per mun.
    """
    input:
        area_stats=(
            rules.datasets_potentialarea_pv_roof_region2_create_area_stats_muns.output.stats
        ),
        unit_stats=(
            rules.datasets_bnetza_mastr_pv_roof_region_create_power_stats_muns.output[0]
        )
    output:
        total=DATASET_PATH / "potentialarea_pv_roof_deployment_stats_muns.csv"
    run:
        pv_installed = pd.read_csv(
            input.unit_stats,
            usecols=["municipality_id", "capacity_net"],
            index_col="municipality_id",
        )

        pv_potential = pd.read_csv(
            input.area_stats,
            usecols=["municipality_id", "installable_power"],
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
        pv_deployed.to_csv(output.total)

rule regionalize_state_targets:
    """
    Calculate PV roof targets of region
    """
    input:
        osm_buildings_stats=get_abs_dataset_path(
            "datasets", "osm_buildings", data_dir=True
        ) / "osm_buildings_ground_area.json",
        el_capacity_targets_bmwk_de=get_abs_dataset_path(
            "preprocessed", "bmwk_long_term_scenarios"
        ) / "data" / "T45-Strom_electricity_installed_power_reformatted.csv",
        el_capacity_targets_mwae_bb=get_abs_dataset_path(
            "datasets","mwae_bb_energy_strategy_region"
        ) / "data" / "mwae_bb_energy_strategy_region.json"
    output:
        DATASET_PATH / "potentialarea_pv_roof_regionalized_targets.json"
    run:
        osm_buildings_stats = load_json(input.osm_buildings_stats)

        # Power target from scenarios
        targets_bmwk_de = pd.read_csv(input.el_capacity_targets_bmwk_de, index_col="year")
        with open(input.el_capacity_targets_mwae_bb, "r") as f:
            targets_mwae_bb = json.load(f)

        target_cap_bmwk_de = targets_bmwk_de.loc[
                         targets_bmwk_de.technology == "pv"
                     ].loc[2045].capacity * 1e3 * config.get("pv_roof_share")

        targets_output = {
            "bmwk_de": {
                "target_power_total":  {
                    "2045": round(
                        target_cap_bmwk_de *
                        osm_buildings_stats["building_ground_area_share_region"]
                    ),
                },
            },
            "mwae_bb": {
                "target_power_total": {
                    year: cap * config.get("pv_roof_share")
                    for year, cap in targets_mwae_bb["re_installed_capacity_pv_mw"].items()
                },
            },
        }

        with open(output[0], "w", encoding="utf8") as f:
            json.dump(targets_output, f, indent=4)
