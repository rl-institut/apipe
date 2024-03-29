"""
Snakefile for this dataset

Note: To include the file in the main workflow, it must be added to the respective module.smk .
"""
import json
import re
import geopandas as gpd
import pandas as pd
from pathlib import Path
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
    "datasets", "potentialarea_wind_region", data_dir=True
)

rule overlay_muns:
    """
    Overlay potential area with municipalities
    """
    input:
        area=get_abs_dataset_path(
            "preprocessed", "rpg_ols_regional_plan") / "data" / "{file}.gpkg",
        region_muns=PATH_TO_REGION_MUNICIPALITIES_GPKG
    output:
        area=DATASET_PATH / "potentialarea_wind_{file}.gpkg"
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
    Create stats on wind potential areas per mun
    """
    input:
        area=expand(
            DATASET_PATH / "potentialarea_wind_{area}.gpkg",
            area=config["areas"],
        ),
        region_muns=PATH_TO_REGION_MUNICIPALITIES_GPKG,
    output: DATASET_PATH / "potentialarea_wind_area_stats_muns.csv"
    run:
        print("Wind potential area stats:")
        muns = gpd.read_file(input.region_muns)
        area_dict = {}

        # Calc areas per area type file
        for file in input.area:
            area_name = re.findall(
                "potentialarea_wind_(.*).gpkg",
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

rule create_captions:
    """
    Create attribute captions for app
    """
    input: rules.datasets_potentialarea_wind_region_create_area_stats_muns.input.area
    output: DATASET_PATH / "potentialarea_wind_attribute_captions.json"
    run:
        captions = {
            "datasets_caption_map": {
                Path(f).stem: "potentialarea_wind" for f in input
            },
            "captions": {
                "potentialarea_wind": config["captions"]
            }
        }
        with open(output[0], "w", encoding="utf8") as f:
            json.dump(captions, f, indent=4)
