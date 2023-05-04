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
    Extract household electricity demand timeseries for districts, merge and
    normalize them
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
    Disaggregate household electricity consumption from districts to
    municipalities for one year
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
        population = pd.read_csv(input.population, header=[0, 1], index_col=0)
        population.columns = population.columns.droplevel(1)
        consumption_district = pd.read_csv(
            input.consumption
        ).set_index("nuts3").sum(axis=1).to_frame(name="consumption_district") * 1e3
        consumption = create.disaggregate_consumption_to_municipality(
            consumption_district=consumption_district,
            muns=gpd.read_file(input.region_muns),
            districts=gpd.read_file(input.region_districts),
            disagg_data=population,
            disagg_data_col=str(wildcards.year),
        )
        consumption.to_csv(output.consumption)

rule hh_merge_consumption_years:
    """
    Merge the electricity consumptions from different years into one
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

rule cts_normalize_timeseries:
    """
    Extract CTS electricity demand timeseries for districts, merge and
    normalize them
    """
    input:
        timeseries=get_abs_dataset_path("preprocessed", "demandregio") /
                   "data" / "dr_cts_power_timeseries_2022.csv",
        region_districts=PATH_TO_REGION_DISTRICTS_GPKG
    output:
        timeseries=DATASET_PATH / "data" / "demand_cts_power_timeseries.csv"
    run:
        create.normalize_filter_timeseries(
            infile=input.timeseries,
            outfile=output.timeseries,
            region_nuts=gpd.read_file(input.region_districts).nuts.to_list(),
        )

rule cts_disaggregate_consumption:
    """
    Disaggregate CTS electricity consumption from districts to municipalities
    for one year
    """
    input:
        consumption=get_abs_dataset_path("preprocessed", "demandregio") /
                    "data" / "dr_cts_power_consumption_{year}.csv",
        employment=get_abs_dataset_path("datasets", "employment_region") /
                   "data" / "employees.csv",
        region_muns=PATH_TO_REGION_MUNICIPALITIES_GPKG,
        region_districts=PATH_TO_REGION_DISTRICTS_GPKG
    output:
        consumption=DATASET_PATH / "data" / "demand_cts_power_consumption_{year}.csv"
    run:
        consumption_district = pd.read_csv(
            input.consumption,
            index_col=0
        ).sum(axis=0).T.to_frame(name="consumption_district")
        consumption = create.disaggregate_consumption_to_municipality(
            consumption_district=consumption_district,
            muns=gpd.read_file(input.region_muns),
            districts=gpd.read_file(input.region_districts),
            disagg_data=pd.read_csv(
                input.employment,
                index_col=0,
            ),
            disagg_data_col="employees"
        )
        consumption.rename(columns={"employees": wildcards.year}).to_csv(output.consumption)

rule cts_merge_consumption_years:
    """
    Merge the electricity consumptions from different years into one
    """
    input:
        consumption=expand(
            DATASET_PATH / "data" / "demand_cts_power_consumption_{year}.csv",
            year=config["cts_electricity_demand"]["years"]
        )
    output:
        consumption = DATASET_PATH / "data" / "demand_cts_power_consumption.csv"
    run:
        create.merge_consumption_multiple_years(
            infiles=input.consumption,
            outfile=output.consumption,
        )

rule ind_normalize_timeseries:
    """
    Extract industry electricity demand timeseries for districts, merge and
    normalize them
    """
    input:
        timeseries=get_abs_dataset_path("preprocessed", "demandregio") /
                   "data" / "dr_ind_power_timeseries_2022.csv",
        region_districts=PATH_TO_REGION_DISTRICTS_GPKG
    output:
        timeseries=DATASET_PATH / "data" / "demand_ind_power_timeseries.csv"
    run:
        create.normalize_filter_timeseries(
            infile=input.timeseries,
            outfile=output.timeseries,
            region_nuts=gpd.read_file(input.region_districts).nuts.to_list(),
        )

rule ind_disaggregate_consumption:
    """
    Disaggregate industry electricity consumption from districts to
    municipalities for one year
    """
    input:
        consumption=get_abs_dataset_path("preprocessed", "demandregio") /
                    "data" / "dr_ind_power_consumption_{year}.csv",
        employment=get_abs_dataset_path("datasets", "employment_region") /
                   "data" / "employees.csv",
        region_muns=PATH_TO_REGION_MUNICIPALITIES_GPKG,
        region_districts=PATH_TO_REGION_DISTRICTS_GPKG
    output:
        consumption=DATASET_PATH / "data" / "demand_ind_power_consumption_{year}.csv"
    run:
        consumption_district = pd.read_csv(
            input.consumption,
            index_col=0
        ).sum(axis=0).T.to_frame(name="consumption_district")
        consumption = create.disaggregate_consumption_to_municipality(
            consumption_district=consumption_district,
            muns=gpd.read_file(input.region_muns),
            districts=gpd.read_file(input.region_districts),
            disagg_data=pd.read_csv(
                input.employment,
                index_col=0,
            ),
            disagg_data_col="employees"
        )
        consumption.rename(columns={"employees": wildcards.year}).to_csv(output.consumption)

rule ind_merge_consumption_years:
    """
    Merge the electricity consumptions from different years into one
    """
    input:
        consumption=expand(
            DATASET_PATH / "data" / "demand_ind_power_consumption_{year}.csv",
            year=config["ind_electricity_demand"]["years"]
        )
    output:
        consumption = DATASET_PATH / "data" / "demand_ind_power_consumption.csv"
    run:
        create.merge_consumption_multiple_years(
            infiles=input.consumption,
            outfile=output.consumption,
        )
