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

rule convert_filter:
    """
    Convert to geopackage, reproject and filter for permanent crops
    """
    input:
        DATASET_PATH / "data" / "temp" / "DFBK_FB.shp"
    output:
        DATASET_PATH / "data" / "DFBK_FB.gpkg"
    run:
        data = gpd.read_file(input[0], layer="DFBK_FB").to_crs(
            GLOBAL_CONFIG["global"]["geodata"]["crs"]
        )
        # filter for permanent crops
        data = data.loc[data["HBN_KAT"] == "DK"]

        write_geofile(
            gdf=data,
            file=output[0],
            layer_name="DFBK_FB",
        )
        shutil.rmtree(DATASET_PATH / "data" / "temp")

rule rasterize:
    """
    Rasterize vector data for further processing in
    potentialarea_pv_ground_region2.
    """
    input:
        DATASET_PATH / "data" / "DFBK_FB.gpkg"
    output:
        DATASET_PATH / "data" / "DFBK_FB.tif"
    params:
        res_m=100
    run:
        shell(
            f"gdal_rasterize -burn 1 -a_nodata 0 "
            f"-tr {params.res_m} {params.res_m} "
            f"{input} {output}"
        )
