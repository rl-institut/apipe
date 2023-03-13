"""
Snakefile for this dataset

Note: To include the file in the main workflow, it must be added to the respective module.smk .
"""
import geopandas as gpd
import pandas as pd

from digipipe.store.utils import (
    get_abs_dataset_path,
    PATH_TO_REGION_MUNICIPALITIES_GPKG,
    PATH_TO_REGION_DISTRICTS_GPKG
)
from digipipe.store.datasets.demand_electricity_region.scripts import create

DATASET_PATH = get_abs_dataset_path("datasets", "demand_electricity_region")

rule hh_normalize_timeseries:
    """
    Extract household demand timeseries for districts, merge and normalize them
    """
    input:
        timeseries=get_abs_dataset_path("preprocessed", "demandregio") /
                   "data" / "dr_hh_power_timeseries_2022.csv",
        region_districts=PATH_TO_REGION_DISTRICTS_GPKG
    output:
        timeseries=DATASET_PATH / "data" / "demand_hh_power_timeseries.csv"
    run:
        create.normalize_filter_timeseries(
            infile=input.timeseries,
            outfile=output.timeseries,
            region_nuts=gpd.read_file(input.region_districts).nuts.to_list(),
        )

rule hh_disaggregate_consumption:
    """
    Disaggregate household consumption from districts to municipalities for one year
    """
    input:
        consumption=get_abs_dataset_path("preprocessed", "demandregio") /
                    "data" / "dr_hh_power_consumption_{year}.csv",
        population=get_abs_dataset_path("datasets", "population_region") /
                   "data" / "population.csv",
        region_muns=PATH_TO_REGION_MUNICIPALITIES_GPKG,
        region_districts=PATH_TO_REGION_DISTRICTS_GPKG
    output:
        consumption=DATASET_PATH / "data" / "demand_hh_power_consumption_{year}.csv"
    run:
        create.disaggregate_consumption_by_pop(
            infile=input.consumption,
            outfile=output.consumption,
            muns=gpd.read_file(input.region_muns),
            districts=gpd.read_file(input.region_districts),
            population=pd.read_csv(
                input.population,
                header=[0, 1],
                index_col=0
            ),
            year=wildcards.year
        )

rule hh_merge_consumption_years:
    """
    Merge the consumptions from different years into one
    """
    input:
        consumption=expand(
            DATASET_PATH / "data" / "demand_hh_power_consumption_{year}.csv",
            year=config["hh_electricity_demand"]["years"]
        )
    output:
        consumption = DATASET_PATH / "data" / "demand_hh_power_consumption.csv"
    run:
        create.merge_consumption_multiple_years(
            infiles=input.consumption,
            outfile=output.consumption,
        )
