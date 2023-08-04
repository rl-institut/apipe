"""
Snakefile for this dataset

Note: To include the file in the main workflow, it must be added to the respective module.smk .
"""

import geopandas as gpd
import pandas as pd

from digipipe.scripts.datasets.mastr import create_stats_per_municipality
from digipipe.store.utils import (
    get_abs_dataset_path,
    PATH_TO_REGION_MUNICIPALITIES_GPKG,
    PATH_TO_REGION_DISTRICTS_GPKG,
)

DATASET_PATH = get_abs_dataset_path("datasets", "bnetza_mastr_pv_roof_region")
SOURCE_DATASET_PATH = get_abs_dataset_path(
    "preprocessed", "bnetza_mastr", data_dir=True
)


rule create:
    """
    Extract roof-mounted PV plants for region
    """
    input:
        units=SOURCE_DATASET_PATH / "bnetza_mastr_solar_raw.csv",
        locations=SOURCE_DATASET_PATH
        / "bnetza_mastr_locations_extended_raw.csv",
        gridconn=SOURCE_DATASET_PATH / "bnetza_mastr_grid_connections_raw.csv",
        unit_correction=get_abs_dataset_path(
            "raw", "bnetza_mastr_correction_region"
        )
        / "data"
        / "bnetza_mastr_pv_roof_region_correction.csv",
        region_muns=PATH_TO_REGION_MUNICIPALITIES_GPKG,
        region_districts=PATH_TO_REGION_DISTRICTS_GPKG,
    output:
        outfile=DATASET_PATH / "data" / "bnetza_mastr_pv_roof_region.gpkg",
        outfile_agg=DATASET_PATH
        / "data"
        / "bnetza_mastr_pv_roof_agg_region.gpkg",
    params:
        config_file=DATASET_PATH / "config.yml",
    script:
        DATASET_PATH / "scripts" / "create.py"


rule create_power_stats_muns:
    """
    Create stats on installed count of units and power per mun
    """
    input:
        units=DATASET_PATH / "data" / "bnetza_mastr_pv_roof_region.gpkg",
        region_muns=PATH_TO_REGION_MUNICIPALITIES_GPKG,
    output:
        DATASET_PATH / "data" / "bnetza_mastr_pv_roof_stats_muns.csv",
    run:
        units = create_stats_per_municipality(
            units_df=gpd.read_file(input.units),
            muns=gpd.read_file(input.region_muns),
            column="capacity_net",
        )
        units["capacity_net"] = units["capacity_net"].div(1e3)  # kW to MW
        units.to_csv(output[0])


rule create_capacity_change_per_year:
    """
    Create stats on temporal change (per year) of total installed capacity and number of units
    """
    input:
        agg_region=DATASET_PATH
        / "data"
        / "bnetza_mastr_pv_roof_agg_region.gpkg",
    output:
        DATASET_PATH
        / "data"
        / "bnetza_mastr_pv_roof_capacity_change_per_year.csv",
    run:
        df = gpd.read_file(input.agg_region)
        df["commissioning_date"] = pd.to_datetime(df["commissioning_date"])
        df_capacity_over_time = (
            df.groupby(df["commissioning_date"].dt.year)["capacity_net"]
            .sum()
            .reset_index()
        )
        df_capacity_over_time.columns = ["year", "capacity_net"]
        df_capacity_over_time["year"] = df_capacity_over_time["year"].astype(
            int
        )
        df_capacity_over_time.to_csv(output[0])
