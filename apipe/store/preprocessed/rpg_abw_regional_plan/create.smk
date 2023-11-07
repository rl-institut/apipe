"""
Snakefile for this dataset

Note: To include the file in the main workflow, it must be added to the respective module.smk .
"""
import geopandas as gpd
from apipe.scripts.geo import (
    rename_filter_attributes,
    reproject_simplify,
    write_geofile
)
from apipe.store.utils import get_abs_dataset_path

DATASET_PATH = get_abs_dataset_path("preprocessed", "rpg_abw_regional_plan")

rule create_stp_2018:
    """
    Sachlicher Teilplan Wind 2018: Preprocess VR/EG
    """
    input:
        get_abs_dataset_path(
            "raw", "rpg_abw_regional_plan") / "data" / "stp_2018_vreg.gpkg"
    output:
        DATASET_PATH / "data" / "stp_2018_vreg.gpkg"
    run:
        data = gpd.read_file(input[0])
        data = reproject_simplify(
            rename_filter_attributes(
                gdf=data,
                attrs_mapping=config["stp_2018"]["attributes"],
            )
        )
        write_geofile(
            gdf=data,
            file=output[0],
            layer_name=config["stp_2018"]["layer"],
        )

rule create_stp_2027_plan_intent:
    """
    Sachlicher Teilplan Wind 2027: Preprocess plan intent
    """
    input:
        get_abs_dataset_path(
            "raw", "rpg_abw_regional_plan") / "data" /
            "stp_2027_ideen_{category}.gpkg"
    output:
        DATASET_PATH / "data" / "stp_2027_{category}.gpkg",
    run:
        data = gpd.read_file(input[0])
        data = reproject_simplify(
            rename_filter_attributes(
                gdf=data,
                attrs_mapping=config["stp_2027"]["vr"]["attributes"],
            )
        )
        write_geofile(
            gdf=data,
            file=output[0],
            layer_name=config["stp_2027"]["vr"]["layer"],
        )

rule create_stp_2027_search_area:
    """
    Sachlicher Teilplan Wind 2027: Preprocess search areas
    """
    input:
        get_abs_dataset_path(
            "raw", "rpg_abw_regional_plan") / "data" /
            "stp_2027_suchraum.gpkg"
    output:
        forest=DATASET_PATH / "data" / "stp_2027_search_area_forest_area.gpkg",
        open_land=DATASET_PATH / "data" / "stp_2027_search_area_open_area.gpkg"
    run:
        # Forest
        data = gpd.read_file(
            input[0],
            layer="suchraum_wald_03032023",
        )
        data = reproject_simplify(
            rename_filter_attributes(
                gdf=data,
                attrs_mapping=config[
                    "stp_2027"]["search_area_forest_area"]["attributes"],
            )
        )
        write_geofile(
            gdf=data,
            file=output.forest,
        )

        # Open land
        data = gpd.read_file(
            input[0],
            layer="suchraum_offenland_03032023",
        )
        data = data.explode(index_parts=False)
        data = reproject_simplify(
            rename_filter_attributes(
                gdf=data,
                attrs_mapping=config[
                    "stp_2027"]["search_area_open_area"]["attributes"],
            )
        )
        write_geofile(
            gdf=data,
            file=output.open_land,
        )
