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
    Create JSON file with stats on pv potential areas per mun
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
                f"{round(float(area_km2.sum()), 1)} sqm"
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
