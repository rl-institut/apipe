"""
Snakefile for this dataset

Note: To include the file in the main workflow, it must be added to the respective module.smk .
"""
import geopandas as gpd
import pandas as pd
import rasterio
import json
import re
from pathlib import Path

from apipe.config import GLOBAL_CONFIG
from apipe.scripts.data_io import load_json
from apipe.store.utils import (
    get_abs_dataset_path,
    PATH_TO_REGION_MUNICIPALITIES_GPKG,
)


DATASET_PATH = get_abs_dataset_path(
    "datasets", "potentialarea_pv_ground_region2"
)

potentialarea_pv_ground_path = get_abs_dataset_path(
    "datasets", "potentialarea_pv_ground", data_dir=True
)


rule clip_raster_to_region_muns:
    """
    Clip given raster to boundaries of specified municipalities,
    creating new raster dataset that includes only the area within boundaries.
    """
    input:
        potentialarea_pv_ground=potentialarea_pv_ground_path / "{file}",
        region_muns=PATH_TO_REGION_MUNICIPALITIES_GPKG,
    output:
        clipped_potentialarea_pv_ground=DATASET_PATH / "data" / "clipped_{file}",
    shell:
        """
        gdalwarp -cutline {input.region_muns} -crop_to_cutline -dstalpha \
        {input.potentialarea_pv_ground} {output.clipped_potentialarea_pv_ground}
        """


rule vectorize_and_add_zonal_stats:
    """
    Convert clipped raster to vector and calculate zonal statistics for each
    polygon, adding these statistics as new attributes to resulting GeoPackage.
    """
    input:
        clipped_raster=DATASET_PATH / "data" / "clipped_{file}.tif",
        region_muns=PATH_TO_REGION_MUNICIPALITIES_GPKG,
    output:
        vector_overlay_gpkg=DATASET_PATH / "data" / "{file}_region.gpkg",
    params:
        script=DATASET_PATH / "scripts" / "create.py",
        area_threshold=config["area_threshold"],
    script:
        DATASET_PATH / "scripts" / "create.py"


rule create_area_stats_muns:
    """
    Create stats on pv potential areas per mun.
    """
    input:
        area=expand(
            DATASET_PATH
            / "data"
            / "potentialarea_pv_ground_{area}_region.gpkg",
            area=config["areas"],
        ),
        region_muns=PATH_TO_REGION_MUNICIPALITIES_GPKG,
    output:
        DATASET_PATH / "data" / "potentialarea_pv_ground_area_stats_muns.csv",
    run:
        print("PV ground potential area stats:")
        muns = gpd.read_file(input.region_muns)
        area_dict = {}

        # Calc areas per area type file
        for file in input.area:
            area_name = re.findall(
                "potentialarea_pv_ground_(.*).gpkg",
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

rule regionalize_state_targets:
    """
    Calculate PV ground targets of region
    """
    input:
        potarea_pv_ground_soil_quality_low=get_abs_dataset_path(
            "datasets", "potentialarea_pv_ground", data_dir=True
        )
            / "potentialarea_pv_ground_soil_quality_low.tif",
        potarea_pv_ground_soil_quality_medium=get_abs_dataset_path(
            "datasets", "potentialarea_pv_ground", data_dir=True
        )
            / "potentialarea_pv_ground_soil_quality_medium.tif",
        potarea_pv_ground_permanent_crops=get_abs_dataset_path(
            "datasets", "potentialarea_pv_ground", data_dir=True
        )
            / "potentialarea_pv_ground_permanent_crops.tif",
        potarea_pv_ground_soil_quality_low_region=get_abs_dataset_path(
            "datasets", "potentialarea_pv_ground_region2", data_dir=True
        )
            / "clipped_potentialarea_pv_ground_soil_quality_low.tif",
        potarea_pv_ground_soil_quality_medium_region=get_abs_dataset_path(
            "datasets", "potentialarea_pv_ground_region2", data_dir=True
        )
            / "clipped_potentialarea_pv_ground_soil_quality_medium.tif",
        potarea_pv_ground_permanent_crops_region=get_abs_dataset_path(
            "datasets", "potentialarea_pv_ground_region2", data_dir=True
        )
            / "clipped_potentialarea_pv_ground_permanent_crops.tif",
        el_capacity_targets=get_abs_dataset_path(
            "preprocessed", "bmwk_long_term_scenarios"
        )
        / "data"
        / "T45-Strom_electricity_installed_power_reformatted.csv",
        tech_data=get_abs_dataset_path(
            "datasets", "technology_data", data_dir=True
        )
            / "technology_data.json",
    output:
        DATASET_PATH
        / "data"
        / "potentialarea_pv_ground_regionalized_targets.json",
    run:
        def calculate_nonzero_pixel_sum(raster_path):
            with rasterio.open(raster_path) as src:
                return (src.read(1) != 0).sum()

        area_low = calculate_nonzero_pixel_sum(input.potarea_pv_ground_soil_quality_low)
        area_medium = calculate_nonzero_pixel_sum(input.potarea_pv_ground_soil_quality_medium)
        area_crops = calculate_nonzero_pixel_sum(input.potarea_pv_ground_permanent_crops)

        area_low_region = calculate_nonzero_pixel_sum(input.potarea_pv_ground_soil_quality_low_region)
        area_medium_region = calculate_nonzero_pixel_sum(input.potarea_pv_ground_soil_quality_medium_region)
        area_crops_region = calculate_nonzero_pixel_sum(input.potarea_pv_ground_permanent_crops_region)

        # Total targets in scenario
        tech_data = json.load(open(input.tech_data))
        targets = pd.read_csv(input.el_capacity_targets, index_col="year")
        target_cap = targets.loc[targets.technology == "pv"].loc[2045].capacity * 1e3

        target_power_soil_quality_low = target_cap * config.get("pv_ground_share") * (area_low_region / area_low)
        target_power_soil_quality_medium = (
            target_cap *  # total capacity Germany
            config.get("pv_ground_agri_share") *  # share of agri PV
            area_medium_region / (area_medium_region + area_crops_region) *  # region's area share medium of total agri PV
            (area_medium_region / area_medium)  #  region's medium area share of Germany's
        )
        target_power_permanent_crops = (
            target_cap *
            config.get("pv_ground_agri_share") *
            area_crops_region / (area_medium_region + area_crops_region) *
            (area_crops_region / area_crops)
        )

        target_area_soil_quality_low = target_power_soil_quality_low / tech_data["power_density"]["pv_ground"]
        target_area_soil_quality_medium = target_power_soil_quality_medium / tech_data["power_density"]["pv_ground_vertical_bifacial"]
        target_area_permanent_crops = target_power_permanent_crops / tech_data["power_density"]["pv_ground_high"]

        target_power_total = target_power_soil_quality_low + target_power_soil_quality_medium + target_power_permanent_crops
        target_area_total = target_area_soil_quality_low + target_area_soil_quality_medium + target_area_permanent_crops

        output_data = {
            "target_power_total": round(target_power_total,2),
            "target_power_agri_soil_quality_low": round(target_power_soil_quality_low,2),
            "target_power_agri_soil_quality_medium": round(target_power_soil_quality_medium,2),
            "target_power_agri_permanent_crops": round(target_power_permanent_crops,2),
            "target_area_total": round(target_area_total,2),
            "target_area_agri_soil_quality_low": round(target_area_soil_quality_low,2),
            "target_area_agri_soil_quality_medium": round(target_area_soil_quality_medium,2),
            "target_area_agri_permanent_crops": round(target_area_permanent_crops,2),
        }

        with open(output[0], "w", encoding="utf8") as f:
            json.dump(output_data, f, indent=4)
