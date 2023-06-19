"""
Snakefile for this dataset

Note: To include the file in the main workflow, it must be added to the respective module.smk .
"""
import json
import re
import geopandas as gpd
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
    "datasets", "potentialarea_wind_region", data_dir=True
)

rule overlay_municipalities:
    """
    Overlay potential area with municipalities
    """
    input:
        area=get_abs_dataset_path(
            "preprocessed", "rpg_abw_regional_plan") / "data" / "{file}.gpkg",
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

rule create_area_stats:
    """
    Create JSON file with stats on wind potential areas per mun
    """
    input:
        area=expand(
            DATASET_PATH / "potentialarea_wind_{area}.gpkg",
            area=config["areas"],
        ),
        region_muns=PATH_TO_REGION_MUNICIPALITIES_GPKG,
    output: DATASET_PATH / "potentialarea_wind_area_stats_muns.json"
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
            data["area_ha"] = data.area / 1e4
            area_ha = data[
                ["municipality_id", "area_ha"]
            ].groupby("municipality_id").sum()

            # Set area of non-occurring muns to 0
            area_ha = area_ha.reindex(muns.id, fill_value=0)
            area_dict[area_name] = area_ha.to_dict()["area_ha"]
            print(
                f"  Total area for {area_name}: "
                f"{round(float(area_ha.sum()), 1)} ha"
            )

        # Dump
        with open(output[0], "w", encoding="utf8") as f:
            json.dump(area_dict, f, indent=4)
