"""
Snakefile for this dataset

Note: To include the file in the main workflow, it must be added to the respective module.smk .
"""

import json
import geopandas as gpd
import pandas as pd
from digipipe.scripts.geo import clip_raster, raster_zonal_stats
from digipipe.store.utils import (
    get_abs_dataset_path,
    get_abs_store_root_path,
    PATH_TO_REGION_DISTRICTS_GPKG
)
from digipipe.store.datasets.demand_electricity_region.scripts.create import (
    normalize_filter_timeseries
)

DATASET_PATH = get_abs_dataset_path("datasets", "demand_heat_region", data_dir=True)

rule raster_clip:
    """
    Clip raster file to boundary
    """
    input:
        heat_demand=get_abs_dataset_path(
            "preprocessed", "seenergies_peta5") / "data" / "HD_2015_{sector}_Peta5_0_1_GJ.tif",
        boundary=get_abs_store_root_path() / "datasets" / "{boundary_dataset}" / "data" / "{boundary_dataset}.gpkg"
    output:
        heat_demand = DATASET_PATH / "HD_2015_{sector}_Peta5_0_1_GJ_{boundary_dataset}.tif"
    run:
        clip_raster(
            raster_file_in=input.heat_demand,
            clip_file=input.boundary,
            raster_file_out=output.heat_demand,
        )

rule raster_convert_to_mwh:
    """
    Convert TIFF raster to geopackage and change unit from GJ to MWh
    Note: Unit conversion takes place after clipping as it takes too long to
          apply to entire dataset for Europe.
    """
    input: DATASET_PATH / "HD_2015_{sector}_Peta5_0_1_GJ_{boundary_dataset}.tif"
    output: DATASET_PATH / "HD_2015_{sector}_Peta5_0_1_MWh_{boundary_dataset}.tif"
    shell:
        """
        gdal_calc.py -A {input} --A_band=1 --outfile={output} --calc="A/3.6"
        rm {input}
        """

rule create_raster_zonal_stats:
    """
    Create zonal heat statistics (heat demand per geographical unit)
    """
    input:
        heat_demand=DATASET_PATH / "HD_2015_{sector}_Peta5_0_1_MWh_{boundary_dataset}.tif",
        clip_file=get_abs_store_root_path() / "datasets" / "{boundary_dataset}" / "data" / "{boundary_dataset}.gpkg"
    output:
        heat_demand=DATASET_PATH / "demand_heat_zonal_stats-{sector}-{boundary_dataset}.gpkg"
    run:
        raster_zonal_stats(
            raster_file_in=input.heat_demand,
            clip_file=input.clip_file,
            zonal_file_out=output.heat_demand,
            var_name="heat_demand",
            stats="sum"
        )

SECTOR_MAPPING = {"hh": "res", "cts": "ser"}

rule heat_demand_shares:
    """
    Calculate buildings' heat demand shares of federal state, region and
    municipalities by merging zonal statistics (for HH and CTS)
    """
    input:
        state=lambda wildcards: DATASET_PATH / f"demand_heat_zonal_stats-{SECTOR_MAPPING[wildcards.sector]}-bkg_vg250_state.gpkg",
        federal_states=lambda wildcards: DATASET_PATH / f"demand_heat_zonal_stats-{SECTOR_MAPPING[wildcards.sector]}-bkg_vg250_federal_states.gpkg",
        region_muns=lambda wildcards: DATASET_PATH / f"demand_heat_zonal_stats-{SECTOR_MAPPING[wildcards.sector]}-bkg_vg250_muns_region.gpkg"
    output:
        #demand_shares=lambda wildcards: SECTOR_MAPPING[wildcards.sector]
        demand_shares=DATASET_PATH / "demand_heat_shares_{sector}.json"
    run:
        # Read demands
        demand_state = gpd.read_file(input.state)
        demand_federal_state = gpd.read_file(input.federal_states)
        demand_region_muns = gpd.read_file(input.region_muns)

        demand_state = float(demand_state.heat_demand)
        demand_federal_state = float(
            demand_federal_state.loc[
                demand_federal_state.nuts == "DEE"].heat_demand)
        demand_region_muns = demand_region_muns.heat_demand

        # Calculate shares
        demand_shares = {
            "federal_state_of_state": demand_federal_state / demand_state,
            "region_of_federal_state": (
                demand_region_muns.sum() / demand_federal_state
            ),
            "muns_of_region": (
                    demand_region_muns / demand_region_muns.sum()
            ).to_dict()
        }

        # Dump
        with open(output.demand_shares, "w", encoding="utf8") as f:
            json.dump(demand_shares, f, indent=4)

rule heat_demand:
    """
    Calculate absolute heat demands for municipalities for HH and CTS
    """
    input:
        demand_germany=get_abs_dataset_path(
            "preprocessed", "ageb_energy_balance") / "data" /
            "ageb_energy_balance_germany_{sector}_twh_2021.csv",
        demand_shares=DATASET_PATH / "demand_heat_shares_{sector}.json",
        demand_future_germany_TNStrom=get_abs_dataset_path(
            "preprocessed", "bmwk_long_term_scenarios") / "data" /
            "TN-Strom_buildings_heating_demand_by_carrier.csv",
        demand_future_germany_T45Strom=get_abs_dataset_path(
            "preprocessed", "bmwk_long_term_scenarios") / "data" /
            "T45-Strom_buildings_heating_demand_by_carrier.csv",
    output:
        DATASET_PATH / "demand_heat_{sector}_2021.csv"
    run:
        ### Demand 2021 ###
        demand_germany=pd.read_csv(input.demand_germany, index_col="carrier")
        with open(input.demand_shares, "r") as f:
            demand_shares=json.load(f)

        # Calc demand share for each municipality
        demand_muns = (
                pd.DataFrame.from_dict(
                    demand_shares.get("muns_of_region"),
                    orient="index"
                )
                * demand_shares.get("federal_state_of_state")
                * demand_shares.get("region_of_federal_state")
                * demand_germany.drop("Strom", axis=0)[
                    ["space_heating", "hot_water", "process_heat"]
                ].sum().sum()
                * 1e6 # TWh to MWh
        )
        demand_muns.index.name = "municipality_id"
        demand_muns.rename(columns={0: 2022}, inplace=True)
        print(f"Demand {wildcards.sector}: ", demand_muns.sum())

        ### Demand 2045: Calc reduction factor ###
        demand_future_germany_TNStrom = pd.read_csv(
            input.demand_future_germany_TNStrom,
            usecols=["Jahr", "Energiebedarf in TWh"],
        ).rename(columns={
            "Jahr": "year",
            "Energiebedarf in TWh": "demand",
        })
        demand_future_germany_T45Strom = pd.read_csv(
            input.demand_future_germany_T45Strom,
            usecols=[" Jahr / Year", "Energiebedarf in TWh / Energy Demand in TWh"],
        ).rename(columns={
            " Jahr / Year": "year",
            "Energiebedarf in TWh / Energy Demand in TWh": "demand",
        })

        reduction_factor = (
            demand_future_germany_T45Strom.loc[
                demand_future_germany_T45Strom.year == 2045].demand.sum() /
            demand_future_germany_TNStrom.loc[
                demand_future_germany_TNStrom.year == 2020].demand.sum()
        )
        print(
            f"Heat demand for sector {wildcards.sector} in 2045: "
            f"{round(reduction_factor, 2)} of 2022's demand."
        )
        demand_muns[2045] = demand_muns[2022] * reduction_factor

        # Dump as CSV
        demand_muns.to_csv(output[0])

rule normalize_timeseries:
    """
    Extract heat demand timeseries for districts, merge and normalize them
    """
    input:
        timeseries=get_abs_dataset_path("preprocessed", "demandregio") /
                   "data" / "dr_{sector}_gas_timeseries_2011.csv",
        region_districts=PATH_TO_REGION_DISTRICTS_GPKG
    output:
        timeseries=DATASET_PATH / "demand_{sector}_heat_timeseries.csv"
    run:
        normalize_filter_timeseries(
            infile=input.timeseries,
            outfile=output.timeseries,
            region_nuts=gpd.read_file(input.region_districts).nuts.to_list(),
        )
