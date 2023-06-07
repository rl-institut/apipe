"""
Snakefile for this dataset

Note: To include the file in the main workflow, it must be added to the respective module.smk .
"""
import difflib
import geopandas as gpd
import pandas as pd

from digipipe.store.utils import (
    get_abs_dataset_path,
    PATH_TO_REGION_MUNICIPALITIES_GPKG,
    PATH_TO_REGION_DISTRICTS_GPKG
)
from digipipe.scripts.datasets.demand import (
    demand_prognosis,
    disaggregate_demand_to_municipality,
    merge_demand_multiple_years,
    normalize_filter_timeseries
)

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
        normalize_filter_timeseries(
            infile=input.timeseries,
            outfile=output.timeseries,
            region_nuts=gpd.read_file(input.region_districts).nuts.to_list(),
        )

rule hh_disaggregate_demand:
    """
    Disaggregate household electricity demand from districts to
    municipalities for one year and create prognosis
    """
    input:
        demand_today_region=get_abs_dataset_path("preprocessed", "demandregio") /
                    "data" / "dr_hh_power_demand_2022.csv",
        demand_future_TN=get_abs_dataset_path(
            "preprocessed","bmwk_long_term_scenarios"
            ) / "data" / "TN-Strom_hh_demand.csv",
        demand_future_T45=get_abs_dataset_path(
            "preprocessed","bmwk_long_term_scenarios"
            ) / "data" / "T45-Strom_hh_demand.csv",
        population=get_abs_dataset_path("datasets", "population_region") /
                   "data" / "population.csv",
        region_muns=PATH_TO_REGION_MUNICIPALITIES_GPKG,
        region_districts=PATH_TO_REGION_DISTRICTS_GPKG
    output:
        demand=DATASET_PATH / "data" / "demand_hh_power_demand_{year}.csv"
    run:
        population = pd.read_csv(input.population, header=[0, 1], index_col=0)
        population.columns = population.columns.droplevel(1)
        # Today's demand
        demand_districts = pd.read_csv(
            input.demand_today_region
        ).set_index("nuts3").sum(axis=1).to_frame(name="demand_districts") * 1e3
        demand = disaggregate_demand_to_municipality(
            demand_districts=demand_districts,
            muns=gpd.read_file(input.region_muns),
            districts=gpd.read_file(input.region_districts),
            disagg_data=population,
            disagg_data_col=str(wildcards.year),
        )
        # Future demand
        if int(wildcards.year) > 2022:
            demand = demand_prognosis(
                demand_future_T45=input.demand_future_T45,
                demand_future_TN=input.demand_future_TN,
                demand_region=demand,
                year_base=2022,
                year_target=int(wildcards.year),
                scale_by="carrier",
                carrier="Strom",
            )
        demand.to_csv(output.demand)

rule hh_merge_demand_years:
    """
    Merge the electricity demands from different years into one
    """
    input:
        demand=expand(
            DATASET_PATH / "data" / "demand_hh_power_demand_{year}.csv",
            year=config["hh_electricity_demand"]["years"]
        )
    output:
        demand = DATASET_PATH / "data" / "demand_hh_power_demand.csv"
    run:
        merge_demand_multiple_years(
            infiles=input.demand,
            outfile=output.demand,
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
        normalize_filter_timeseries(
            infile=input.timeseries,
            outfile=output.timeseries,
            region_nuts=gpd.read_file(input.region_districts).nuts.to_list(),
        )

rule cts_disaggregate_demand:
    """
    Disaggregate CTS electricity demand from districts to municipalities
    for one year
    """
    input:
        demand_today_region=get_abs_dataset_path(
            "preprocessed", "demandregio") / "data" /
            "dr_cts_power_demand_2022.csv",
        demand_future_TN=get_abs_dataset_path(
            "preprocessed","bmwk_long_term_scenarios"
            ) / "data" / "TN-Strom_cts_demand.csv",
        demand_future_T45=get_abs_dataset_path(
        "preprocessed","bmwk_long_term_scenarios"
            ) / "data" / "T45-Strom_cts_demand.csv",
        employment=get_abs_dataset_path("datasets", "employment_region") /
                   "data" / "employment.csv",
        region_muns=PATH_TO_REGION_MUNICIPALITIES_GPKG,
        region_districts=PATH_TO_REGION_DISTRICTS_GPKG
    output:
        demand=DATASET_PATH / "data" / "demand_cts_power_demand_{year}.csv"
    run:
        # Today's demand
        demand_districts = pd.read_csv(
            input.demand_today_region,
            index_col=0
        ).sum(axis=0).T.to_frame(name="demand_districts")
        demand = disaggregate_demand_to_municipality(
            demand_districts=demand_districts,
            muns=gpd.read_file(input.region_muns),
            districts=gpd.read_file(input.region_districts),
            disagg_data=pd.read_csv(
                input.employment,
                index_col=0,
            ),
            disagg_data_col="employees_total"
        )
        # Future demand
        if int(wildcards.year) > 2022:
            demand = demand_prognosis(
                demand_future_T45=input.demand_future_T45,
                demand_future_TN=input.demand_future_TN,
                demand_region=demand,
                year_base=2022,
                year_target=int(wildcards.year),
                scale_by="carrier",
                carrier="Strom",
            )
        demand.rename(columns={"employees_total": wildcards.year}).to_csv(output.demand)

rule cts_merge_demand_years:
    """
    Merge the electricity demands from different years into one
    """
    input:
        demand=expand(
            DATASET_PATH / "data" / "demand_cts_power_demand_{year}.csv",
            year=config["cts_electricity_demand"]["years"]
        )
    output:
        demand = DATASET_PATH / "data" / "demand_cts_power_demand.csv"
    run:
        merge_demand_multiple_years(
            infiles=input.demand,
            outfile=output.demand,
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
        normalize_filter_timeseries(
            infile=input.timeseries,
            outfile=output.timeseries,
            region_nuts=gpd.read_file(input.region_districts).nuts.to_list(),
        )

rule ind_disaggregate_demand:
    """
    Disaggregate industry electricity demand from districts to
    municipalities for one year
    """
    input:
        demand_today_region_dr=get_abs_dataset_path(
            "preprocessed", "demandregio") /
            "data" / "dr_ind_power_demand_2022.csv",
        demand_today_region_stala=get_abs_dataset_path(
            "preprocessed", "stala_st_energy") /
            "data" / "power_demand_industry_st_districts.csv",
        demand_future_TN=get_abs_dataset_path(
            "preprocessed","bmwk_long_term_scenarios"
            ) / "data" / "TN-Strom_ind_demand.csv",
        demand_future_T45=get_abs_dataset_path(
            "preprocessed","bmwk_long_term_scenarios"
            ) / "data" / "T45-Strom_ind_demand.csv",
        employment=get_abs_dataset_path("datasets", "employment_region") /
                   "data" / "employment.csv",
        region_muns=PATH_TO_REGION_MUNICIPALITIES_GPKG,
        region_districts=PATH_TO_REGION_DISTRICTS_GPKG
    output:
        demand=DATASET_PATH / "data" / "demand_ind_power_demand_{year}.csv"
    run:
        districts = gpd.read_file(input.region_districts)
        # Today's demand: use dataset defined in config
        if config["ind_electricity_demand_source"] == "demandregio":
            demand_districts = pd.read_csv(
                input.demand_today_region_dr,
                index_col=0
            ).sum(axis=0).T.to_frame(name="demand_districts")
        elif config["ind_electricity_demand_source"] == "stala":
            demand_districts = pd.read_csv(
                input.demand_today_region_stala,
                usecols=["name", "2021"],
                index_col="name"
            )
            districts["name"] = districts["name"].map(
                lambda _: difflib.get_close_matches(
                    _, demand_districts.index)[0]
            )
            demand_districts = demand_districts.merge(districts,on="name")[
                ["nuts", "2021"]].rename(columns={
                "2021": "demand_districts"}).set_index("nuts")
        else:
            raise ValueError(
                "ind_electricity_demand_source must be one of "
                "['demandregio', 'stala']"
            )

        # Disaggregate
        demand = disaggregate_demand_to_municipality(
            demand_districts=demand_districts,
            muns=gpd.read_file(input.region_muns),
            districts=districts,
            disagg_data=pd.read_csv(
                input.employment,
                index_col=0,
            ),
            disagg_data_col="employees_ind"
        )

        # Future demand
        if int(wildcards.year) > 2022:
            demand = demand_prognosis(
                demand_future_T45=input.demand_future_T45,
                demand_future_TN=input.demand_future_TN,
                demand_region=demand,
                year_base=2022,
                year_target=int(wildcards.year),
                scale_by="total",
                carrier="Strom",
            )
        demand.rename(
            columns={"employees_ind": wildcards.year}
        ).to_csv(output.demand)

rule ind_merge_demand_years:
    """
    Merge the electricity demands from different years into one
    """
    input:
        demand=expand(
            DATASET_PATH / "data" / "demand_ind_power_demand_{year}.csv",
            year=config["ind_electricity_demand"]["years"]
        )
    output:
        demand = DATASET_PATH / "data" / "demand_ind_power_demand.csv"
    run:
        merge_demand_multiple_years(
            infiles=input.demand,
            outfile=output.demand,
        )
