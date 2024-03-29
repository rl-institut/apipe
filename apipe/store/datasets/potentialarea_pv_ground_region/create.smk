"""
Snakefile for this dataset

Note: To include the file in the main workflow, it must be added to the respective module.smk .
"""
import json
import re
import geopandas as gpd
import pandas as pd
from pathlib import Path

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
    "datasets", "potentialarea_pv_ground_region", data_dir=True
)

rule overlay_muns:
    """
    Overlay potential area with municipalities
    """
    input:
        area=get_abs_dataset_path(
            "datasets", "rli_pv_wfr_region", data_dir=True) / "{file}",
        region_muns=PATH_TO_REGION_MUNICIPALITIES_GPKG
    output:
        area=DATASET_PATH / "{file}"
    run:
        data = gpd.read_file(input.area)
        data = overlay(
            gdf=data,
            gdf_overlay=gpd.read_file(input.region_muns),
            retain_rename_overlay_columns={"id": "municipality_id"},
        )
        write_geofile(
            gdf=convert_to_multipolygon(data),
            file=output.area,
        )

rule create_area_stats_muns:
    """
    Create stats on pv potential areas per mun
    """
    input:
        area=expand(
            DATASET_PATH / "potentialarea_pv_{area}.gpkg",
            area=config["areas"],
        ),
        region_muns=PATH_TO_REGION_MUNICIPALITIES_GPKG,
    output: DATASET_PATH / "potentialarea_pv_ground_area_stats_muns.csv"
    run:
        print("PV ground potential area stats:")
        muns = gpd.read_file(input.region_muns)
        area_dict = {}

        # Calc areas per area type file
        for file in input.area:
            area_name = re.findall(
                "potentialarea_pv_(.*).gpkg",
                Path(file).name,
            )[0]
            data = gpd.read_file(file)
            data["area_km2"] = data.area / 1e6
            area_km2 = data[
                ["municipality_id", "area_km2"]
            ].groupby("municipality_id").sum()

            # Set area of non-occurring muns to 0
            area_km2 = area_km2.reindex(muns.id, fill_value=0)
            area_dict[area_name] = area_km2.to_dict()["area_km2"]
            print(
                f"  Total area for {area_name}: "
                f"{round(float(area_km2.sum()), 1)} sqkm"
            )

        area_df = pd.DataFrame(area_dict)
        area_df.index.name="municipality_id"
        area_df.to_csv(output[0])

rule create_potarea_shares:
    """
    Calc shares of actual potential areas in total potential areas (per type)
    """
    input:
        potarea_pv_road_railway=get_abs_dataset_path(
            "datasets", "rli_pv_wfr_region", data_dir=True
        ) / "potentialarea_pv_road_railway_region.gpkg",
        road_railway_500m=get_abs_dataset_path(
            "datasets", "rli_pv_wfr_region", data_dir=True
        ) / "road_railway-500m_region.gpkg",
        potarea_pv_agri=get_abs_dataset_path(
            "datasets", "rli_pv_wfr_region", data_dir=True
        ) / "potentialarea_pv_agriculture_lfa-off_region.gpkg",
        soil_quality_low=get_abs_dataset_path(
            "datasets", "rli_pv_wfr_region", data_dir=True
        ) / "soil_quality_low_region.gpkg",

    output: DATASET_PATH / "potentialarea_pv_ground_area_shares.json"
    run:
        area_dict = {
            "road_railway": round(
                gpd.read_file(input.potarea_pv_road_railway).area.sum() /
                gpd.read_file(input.road_railway_500m).area.sum(),
                3
            ),
            "agri": round(
                gpd.read_file(input.potarea_pv_agri).area.sum() /
                gpd.read_file(input.soil_quality_low).area.sum(),
                3
            ),
        }

        # Dump
        with open(output[0], "w", encoding="utf8") as f:
            json.dump(area_dict, f, indent=4)

rule regionalize_state_targets:
    """
    Calculate PV ground targets of region
    """
    input:
        potarea_pv_road_railway=get_abs_dataset_path(
            "preprocessed", "rli_pv_wfr", data_dir=True
        ) / "potentialarea_pv_road_railway.gpkg",
        potarea_pv_agri=get_abs_dataset_path(
            "preprocessed", "rli_pv_wfr", data_dir=True
        ) / "potentialarea_pv_agriculture_lfa-off.gpkg",
        potarea_pv_road_railway_region=get_abs_dataset_path(
            "datasets", "rli_pv_wfr_region", data_dir=True
        ) / "potentialarea_pv_road_railway_region.gpkg",
        potarea_pv_agri_region=get_abs_dataset_path(
            "datasets", "rli_pv_wfr_region", data_dir=True
        ) / "potentialarea_pv_agriculture_lfa-off_region.gpkg",
        el_capacity_targets=get_abs_dataset_path(
            "preprocessed", "bmwk_long_term_scenarios"
        ) / "data" / "T45-Strom_electricity_installed_power_reformatted.csv",
        tech_data=get_abs_dataset_path(
            "datasets", "technology_data", data_dir=True
        ) / "technology_data.json"
    output:
        DATASET_PATH / "potentialarea_pv_ground_regionalized_targets.json"
    run:
        potarea_pv_rr = gpd.read_file(input.potarea_pv_road_railway).to_crs(
            GLOBAL_CONFIG["global"]["geodata"]["crs"]
        ).area.sum()
        potarea_pv_agri = gpd.read_file(input.potarea_pv_agri).to_crs(
            GLOBAL_CONFIG["global"]["geodata"]["crs"]
        ).area.sum()
        potarea_pv_rr_region = gpd.read_file(
            input.potarea_pv_road_railway_region).area.sum()
        potarea_pv_agri_region = gpd.read_file(
            input.potarea_pv_agri_region).area.sum()
        tech_data = load_json(input.tech_data)

        # Power target from longterm scenario
        targets = pd.read_csv(input.el_capacity_targets, index_col="year")
        target_cap = targets.loc[
                         targets.technology == "pv"
                     ].loc[2045].capacity * 1e3 * config.get("pv_ground_share")
        targets_region = (
            target_cap * (
                (potarea_pv_rr_region + potarea_pv_agri_region) /
                (potarea_pv_rr + potarea_pv_agri)
            )
        )

        with open(output[0], "w", encoding="utf8") as f:
            json.dump(
                {
                    # Power targets (disaggregated)
                    "target_power_total": round(
                        target_cap * (
                            (potarea_pv_rr_region + potarea_pv_agri_region) /
                            (potarea_pv_rr + potarea_pv_agri)
                        )
                    ),
                    "target_power_road_railway": round(
                        target_cap * potarea_pv_rr_region / (
                            potarea_pv_rr + potarea_pv_agri)
                    ),
                    "target_power_agri": round(
                        target_cap * potarea_pv_agri_region / (
                            potarea_pv_rr + potarea_pv_agri)
                    ),
                    # Areas targets (from power targets)
                    "target_area_total": round(
                        target_cap * (
                            (potarea_pv_rr_region + potarea_pv_agri_region) /
                            (potarea_pv_rr + potarea_pv_agri)
                        ) / tech_data["power_density"]["pv_ground"]
                        ,2
                    ),
                    "target_area_road_railway": round(
                        target_cap * potarea_pv_rr_region / (
                            potarea_pv_rr + potarea_pv_agri
                        ) / tech_data["power_density"]["pv_ground"]
                        ,2
                    ),
                    "target_area_agri": round(
                        target_cap * potarea_pv_agri_region / (
                            potarea_pv_rr + potarea_pv_agri
                        ) / tech_data["power_density"]["pv_ground"]
                        ,2
                    ),
                },
                f,
                indent=4
            )
