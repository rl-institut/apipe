"""
Snakefile for this dataset

Note: To include the file in the main workflow, it must be added to the respective module.smk .
"""
import fiona
import os
import geopandas as gpd
import pandas as pd
from apipe.scripts.geo import (
    rename_filter_attributes,
    reproject_simplify,
    write_geofile,
    convert_to_multipolygon
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
            "PV_FFA_OLS_Stand_November2023.gpkg"
    output:
        DATASET_PATH / "data" / "rpg_ols_pv_ground.gpkg"
    run:
        data = reproject_simplify(
            rename_filter_attributes(
                gdf=gpd.read_file(input[0]),
                attrs_filter_by_values=config["pv_ground"]["attributes_filter"],
                attrs_mapping=config["pv_ground"]["attributes"],
            )
        )
        write_geofile(
            gdf=data,
            file=output[0],
            layer_name=config["pv_ground"]["layer"],
        )

rule create_pv_ground_criteria:
    """
    Freiflächen-Photovoltaikanlagen Negativkriterien: Preprocess
    """
    input:
        get_abs_dataset_path(
            "raw", "rpg_ols_regional_plan") / "data" /
            "Negativkriterien_PV_RPG_OLS_07032024.gpkg"
    output:
        [DATASET_PATH / "data" / f"{fname}.gpkg"
         for fname in set(config["pv_ground_criteria"]["layers"].values())
         if fname != ""]
    run:
        target_layers = [
            layer for layer in fiona.listlayers(input[0])
            if config["pv_ground_criteria"]["layers"].get(layer) != ""
        ]
        for target_layer in target_layers:
            print(f"Processing layer: {target_layer}"),
            data = reproject_simplify(
                rename_filter_attributes(
                    gdf=gpd.read_file(input[0], layer=target_layer),
                    attrs_mapping=config["pv_ground"]["attributes"],
                ),
                fix_geom=True,
            )
            target_file = (
                DATASET_PATH / "data" /
                f'{config["pv_ground_criteria"]["layers"].get(target_layer)}.gpkg'
            )

            # If 2 layers have same target: Append and union data
            if os.path.exists(target_file):
                print(
                    f"File {target_file} exists, merging new data and overwrite "
                    f"exiting file..."
                )
                data_existing = gpd.read_file(target_file)
                data = pd.concat([data, data_existing])
                os.remove(target_file)

            data = gpd.GeoDataFrame(
                crs=data.crs.srs, geometry=[data.unary_union]
            )

            write_geofile(
                gdf=convert_to_multipolygon(data),
                file=target_file,
                layer_name=target_layer,
            )
