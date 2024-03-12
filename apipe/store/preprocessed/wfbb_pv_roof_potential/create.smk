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

DATASET_PATH = get_abs_dataset_path("preprocessed", "wfbb_pv_roof_potential")

rule create:
    input:
        get_abs_dataset_path(
            "raw", "wfbb_pv_roof_potential"
        ) / "data" / "solaratlas_eignung_dachflaechen_pv.gpkg"
    output:
        DATASET_PATH / "data" / "solaratlas_eignung_dachflaechen_pv.gpkg"
    run:
        data = reproject_simplify(
            rename_filter_attributes(
                gdf=gpd.read_file(input[0]),
                attrs_mapping=config["attributes"]
            )
        )
        write_geofile(
            gdf=data,
            file=output[0],
            layer_name=config["layer"],
        )
