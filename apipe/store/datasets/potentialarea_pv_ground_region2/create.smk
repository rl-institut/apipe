"""
Snakefile for this dataset

Note: To include the file in the main workflow, it must be added to the respective module.smk .
"""
import geopandas as gpd
import rasterio
from rasterio.features import shapes
import pandas as pd
from shapely.geometry import shape

from apipe.store.utils import (
    get_abs_dataset_path,
    PATH_TO_REGION_MUNICIPALITIES_GPKG,
)
from apipe.scripts.geo import (
    overlay,
    convert_to_multipolygon,
    write_geofile,
    raster_zonal_stats

)

### ToDo
### Expected Result files (statistics, see dataset "potentialarea_pv_ground_region" in digipipe)
# - [ ] per mun: potentialarea_pv_ground_area_stats_muns.csv
# - [ ] total shares: potentialarea_pv_ground_area_shares.json
# - [ ] targets: potentialarea_pv_ground_regionalized_targets.json
#             ~~~~ sehr an rule in "potentialarea_pv_ground_region" orientieren, oder hier was grundlegend anders?


DATASET_PATH = get_abs_dataset_path(
    "datasets", "potentialarea_pv_ground_region2"
)

potentialarea_pv_ground_path = get_abs_dataset_path(
    "datasets", "potentialarea_pv_ground", data_dir=True
)

rule clip_raster_to_region_muns:
    input:
        potentialarea_pv_ground =
        potentialarea_pv_ground_path / "{file}",
        region_muns=PATH_TO_REGION_MUNICIPALITIES_GPKG
    output:
        clipped_potentialarea_pv_ground = DATASET_PATH / "data" / "clipped_{file}",
    shell:
        """
        gdalwarp -cutline {input.region_muns} -crop_to_cutline -dstalpha \
        {input.potentialarea_pv_ground} {output.clipped_potentialarea_pv_ground}
        """

rule vectorize_and_add_zonal_stats:
    input:
        clipped_raster=DATASET_PATH / "data" / "clipped_{file}.tif",
        region_muns=PATH_TO_REGION_MUNICIPALITIES_GPKG
    output:
        vector_overlay_gpkg=DATASET_PATH / "data" / "{file}_region.gpkg"
    params:
        script=DATASET_PATH / "scripts" / "create.py",
        area_threshold=config["area_threshold"]
    script:
        DATASET_PATH / "scripts" / "create.py"

rule create_area_stats_muns:
    input:
        area=expand(
            DATASET_PATH / "data" / "overlay_potentialarea_pv_ground_{area}_region.gpkg",
            area=config["areas"],
        ),
        region_muns=PATH_TO_REGION_MUNICIPALITIES_GPKG
    output:
        DATASET_PATH / "data" / "potentialarea_pv_ground_area_stats_muns.csv"
    run:
        print("PV ground potential area stats:")
        muns = gpd.read_file(input.region_muns)
        area_dict = {}

        # Calc areas per area type file
        for file in input.area:
            area_name = re.findall(
                "potentialarea_pv_(.*).gpkg",
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
