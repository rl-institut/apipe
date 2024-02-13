"""
Snakefile for this dataset

Note: To include the file in the main workflow, it must be added to the respective module.smk .
"""

from apipe.store.utils import get_abs_dataset_path

DATASET_PATH = get_abs_dataset_path("preprocessed", "bgr_sqr", data_dir=True)

rule unzip:
    input:
        get_abs_dataset_path("raw", "bgr_sqr") / "data" / "sqr1000_250_v10.zip"
    output:
        DATASET_PATH / "sqr1000_250_v10.tif"
    params:
        outpath=DATASET_PATH,
        file_to_extract= "sqr1000_250_v10.tif"
    shell:
        """
        unzip -j {input} {params.file_to_extract} -d {params.outpath}
        """

rule fix_crs:
    input:
        tif_in=DATASET_PATH / "sqr1000_250_v10.tif"
    output:
        tif_out=DATASET_PATH / "sqr1000_250_v10_3035.tif"
    shell:
        """
        gdalwarp -t_srs EPSG:3035 -s_srs EPSG:3034 {input.tif_in} {output.tif_out}
        """


# rule convert_to_gpkg:
#     input:
#         tif_file=DATASET_PATH / "sqr1000_250_v10.tif"
#     output:
#         gpkg_file=DATASET_PATH / "sqr1000_250_v10.gpkg"
#     shell:
#         "gdal_translate {input.tif_file} {output.gpkg_file} -of GPKG"
