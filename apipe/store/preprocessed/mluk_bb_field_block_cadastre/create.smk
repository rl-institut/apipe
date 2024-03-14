"""
Snakefile for this dataset

Note: To include the file in the main workflow, it must be added to the respective module.smk .
"""
import shutil
import geopandas as gpd
from apipe.config import GLOBAL_CONFIG
from apipe.scripts.geo import write_geofile
from apipe.store.utils import get_abs_dataset_path

DATASET_PATH = get_abs_dataset_path("preprocessed", "mluk_bb_field_block_cadastre")

rule unzip:
    """
    Unzip file
    """
    input:
        get_abs_dataset_path("raw", "mluk_bb_field_block_cadastre") / "data" / "dfbk.zip"
    output:
        files = DATASET_PATH / "data" / "temp" / "DFBK_FB.shp"
    params:
        outpath=DATASET_PATH / "data" / "temp"
    shell:
        """
        unzip -j {input} -d {params.outpath}
        """

rule convert:
    """
    Convert to geopackage and reproject
    """
    input:
        DATASET_PATH / "data" / "temp" / "DFBK_FB.shp"
    output:
        DATASET_PATH / "data" / "DFBK_FB.gpkg"
    run:
        data = gpd.read_file(input[0], layer="DFBK_FB").to_crs(
            GLOBAL_CONFIG["global"]["geodata"]["crs"]
        )
        write_geofile(
            gdf=data,
            file=output[0],
            layer_name="DFBK_FB",
        )
        shutil.rmtree(DATASET_PATH / "data" / "temp")
