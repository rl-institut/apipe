"""
Snakefile for this dataset

Note: To include the file in the main workflow, it must be added to the respective module.smk .
"""

from apipe.store.utils import get_abs_dataset_path
import numpy as np
import rasterio
from rasterio.warp import reproject
from rasterio.enums import Resampling

DATASET_PATH = get_abs_dataset_path("datasets", "potentialarea_pv_ground")

OEI_AGRI_PATH = get_abs_dataset_path("preprocessed", "oei_agri_pv", data_dir=True)

# Path to 'BGR SQR' (3035)
BGR_SQR_PATH = get_abs_dataset_path("preprocessed", "bgr_sqr", data_dir=True) / "sqr1000_250_v10_3035.tif"

# Path to 'SQR Gesamt Agri'
AGRI_SQR_TOTAL_PATH = OEI_AGRI_PATH / "Agri-PV-Potenziale_Gesamt_100x100_EPSG3035.tif"

# Path to 'SQR 50-70'
AGRI_SQR_50_70_PATH = OEI_AGRI_PATH / "Agri-PV-Potenziale_SQR_50-70_100x100_EPSG3035_pos.tif"

# Path to 'Dauerkulturen'
PERMANENT_CROPS_PATH = OEI_AGRI_PATH / "Agri-PV-Potenziale_Dauerkulturen_100x100_EPSG3035.tif"


rule calculate_raster_intersection_sq_low:
    input:
        bgr_sqr=BGR_SQR_PATH,
        agri_sqr_total=AGRI_SQR_TOTAL_PATH,
        permanent_crops=PERMANENT_CROPS_PATH
    output:
        DATASET_PATH / "data" / "potentialarea_pv_ground_soil_quality_low.tif"
    run:
        with rasterio.open(input.bgr_sqr) as bgr_sqr_src, rasterio.open(input.agri_sqr_total) as agri_sqr_total_src, rasterio.open(input.permanent_crops) as permanent_crops_src:
            out_meta = agri_sqr_total_src.meta.copy()

            # Create empty data structures for the reprojected rasters
            bgr_sqr_reproj = np.empty(
                (agri_sqr_total_src.height, agri_sqr_total_src.width), dtype=rasterio.float32
            )
            permanent_crops_reproj = np.empty(
                (agri_sqr_total_src.height, agri_sqr_total_src.width), dtype=rasterio.float32
            )

            # Bring raster data for 'bgr_sqr' and 'permanent_crops' to the same resolution and extent as 'agri_sqr_total'
            reproject(
                source=rasterio.band(bgr_sqr_src, 1),
                destination=bgr_sqr_reproj,
                src_transform=bgr_sqr_src.transform,
                src_crs=bgr_sqr_src.crs,
                dst_transform=agri_sqr_total_src.transform,
                dst_crs=agri_sqr_total_src.crs,
                resampling=Resampling.nearest,
            )

            reproject(
                source=rasterio.band(permanent_crops_src, 1),
                destination=permanent_crops_reproj,
                src_transform=permanent_crops_src.transform,
                src_crs=permanent_crops_src.crs,
                dst_transform=agri_sqr_total_src.transform,
                dst_crs=agri_sqr_total_src.crs,
                resampling=Resampling.nearest,
            )

            agri_sqr_total_data = agri_sqr_total_src.read(1)

            mask = (bgr_sqr_reproj > 50) | (permanent_crops_reproj > 0) | (agri_sqr_total_data <= 0)
            agri_sqr_total_data[mask] = 0
            out_meta['nodata'] = 0

            with rasterio.open(output[0], "w", **out_meta) as out_raster:
                out_raster.write(agri_sqr_total_data, 1)

rule calculate_raster_intersection_sq_medium:
    input:
        agri_sqr_50_70=AGRI_SQR_50_70_PATH,
        permanent_crops=PERMANENT_CROPS_PATH
    output:
        DATASET_PATH / "data" / "potentialarea_pv_ground_soil_quality_medium.tif"
    run:
        with rasterio.open(input.agri_sqr_50_70) as agri_sqr_50_70_src, rasterio.open(input.permanent_crops) as permanent_crops_src:
            out_meta = agri_sqr_50_70_src.meta.copy()

            # Create an empty data structure for the reprojected raster 'permanent_crops'
            permanent_crops_reproj = np.empty(
                (agri_sqr_50_70_src.height, agri_sqr_50_70_src.width), dtype=rasterio.float32
            )

            # Bring raster data for 'permanent_crops' to the same resolution and extent as 'agri_sqr_50_70'
            reproject(
                source=rasterio.band(permanent_crops_src, 1),
                destination=permanent_crops_reproj,
                src_transform=permanent_crops_src.transform,
                src_crs=permanent_crops_src.crs,
                dst_transform=agri_sqr_50_70_src.transform,
                dst_crs=agri_sqr_50_70_src.crs,
                resampling=Resampling.nearest,
            )

            agri_sqr_50_70_data = agri_sqr_50_70_src.read(1)

            mask = (permanent_crops_reproj > 0) | (agri_sqr_50_70_data <= 0)
            agri_sqr_50_70_data[mask] = 0
            out_meta['nodata'] = 0

            with rasterio.open(output[0], "w", **out_meta) as out_raster:
                out_raster.write(agri_sqr_50_70_data, 1)

rule create_adjusted_permanent_crops_tif:
    input:
        permanent_crops=PERMANENT_CROPS_PATH
    output:
        DATASET_PATH / "data" / "potentialarea_pv_ground_permanent_crops.tif"
    shell:
        """
        gdal_translate {input.permanent_crops} {output[0]} -of GTiff -ot Float32 -a_nodata 0
        """
