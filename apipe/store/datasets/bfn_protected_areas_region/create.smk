"""
Snakefile for this dataset

Note: To include the file in the main workflow, it must be added to the respective module.smk .
"""

import geopandas as gpd
from apipe.config import GLOBAL_CONFIG
from apipe.scripts.geo import (
    overlay,
    convert_to_multipolygon,
    write_geofile
)
from apipe.store.utils import get_abs_dataset_path

DATASET_PATH = get_abs_dataset_path(
    "datasets", "bfn_protected_areas_region", data_dir=True)


rule clip_to_region:
    """
    Clip to region and reproject
    """
    input:
        geodata=get_abs_dataset_path(
            "raw", "bfn_protected_areas", data_dir=True
        ) / "{file}.gpkg",
        region=rules.datasets_bkg_vg250_region_create.output
    output:
        geodata=DATASET_PATH / "{file}_region.gpkg"
    run:
        print(f"Clipping layer {wildcards.file} to region...")

        geodata = gpd.read_file(input.geodata).to_crs(
            GLOBAL_CONFIG["global"]["geodata"]["crs"]
        )
        geodata = overlay(
            gdf=geodata,
            gdf_overlay=gpd.read_file(input.region[0]),
        )

        if len(geodata) == 0:
            print("  Layer has no data in region!")

        write_geofile(
            gdf=convert_to_multipolygon(geodata),
            file=output.geodata,
            layer_name=wildcards.file,
        )
