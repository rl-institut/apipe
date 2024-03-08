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

DATASET_PATH = get_abs_dataset_path("preprocessed", "rpg_ols_regional_plan")

rule create_stp_2018:
    """
    Sachlicher Teilplan Wind 2018: Preprocess VR/EG
    """
    input:
        get_abs_dataset_path(
            "raw", "rpg_ols_regional_plan") / "data" /
        "Windeignungsgebiete_Satzung_2018_OLS.gpkg"
    output:
        DATASET_PATH / "data" / "stp_2018_eg.gpkg"
    run:
        data = reproject_simplify(
            rename_filter_attributes(
                gdf=gpd.read_file(input[0]),
                attrs_mapping=config["stp_2018"]["attributes"],
            )
        )
        write_geofile(
            gdf=data,
            file=output[0],
            layer_name=config["stp_2018"]["layer"],
        )

rule create_stp_2024:
    """
    Sachlicher Teilplan Wind 2024: Preprocess VR
    """
    input:
        get_abs_dataset_path(
            "raw", "rpg_ols_regional_plan") / "data" /
            "VR_Windenergienutzung_2024_RPG_OLS.gpkg"
    output:
        DATASET_PATH / "data" / "stp_2024_vr.gpkg"
    run:
        data = reproject_simplify(
            rename_filter_attributes(
                gdf=gpd.read_file(input[0]),
                attrs_mapping=config["stp_2024"]["attributes"],
            )
        )
        write_geofile(
            gdf=data,
            file=output[0],
            layer_name=config["stp_2024"]["layer"],
        )

rule create_pv_ground:
    """
    Freiflächen-Photovoltaikanlagen: Preprocess
    """
    input:
        get_abs_dataset_path(
            "raw", "rpg_ols_regional_plan") / "data" /
            "PV_FFA_OLS_Stand_Sommer2023.gpkg"
    output:
        DATASET_PATH / "data" / "pv_ground.gpkg"
    run:
        data = reproject_simplify(
            rename_filter_attributes(
                gdf=gpd.read_file(input[0]),
                attrs_mapping=config["pv_ground"]["attributes"],
            )
        )
        write_geofile(
            gdf=data,
            file=output[0],
            layer_name=config["pv_ground"]["layer"],
        )
