"""Shared functions for processing electricity and heat demand datasets"""

from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd


def normalize_filter_timeseries(
    infile: Path, outfile: Path, region_nuts: list = None
) -> None:
    """Extract timeseries for specific districts, merge and normalize them to
    1 (MWh)

    Parameters
    ----------
    infile : pathlib.Path
        Path to timeseries CSV with NUTS 3 codes in columns and
        timesteps in rows (expects 15min ot 1h resolution)
    outfile : pathlib.Path
        Path to CSV outfile - aggregated and normalized timeseries
        (1h resolution)
    region_nuts : list of str
        List of NUTS 3 codes whose timeseries are to be used

    Returns
    -------
    None
    """
    # Get timeseries, create timeindex, filter region, resample to 1h
    timeseries = pd.read_csv(infile)
    timeseries = timeseries[region_nuts]

    if len(timeseries) == 8760:
        timeseries.index = pd.DatetimeIndex(
            pd.date_range(
                start="2017-01-01 00:00:00",
                end="2017-12-31 23:45:00",
                freq="1H",
            )
        )
    elif len(timeseries) == 35040:
        timeseries.index = pd.DatetimeIndex(
            pd.date_range(
                start="2017-01-01 00:00:00",
                end="2017-12-31 23:45:00",
                freq="15Min",
            )
        )
        timeseries = timeseries.resample("H").sum()
    else:
        raise ValueError("Invalid number of rows in timeseries!")

    # Average SLP timeseries and normalize to 1 MWh
    timeseries = timeseries.sum(axis=1)
    timeseries = timeseries.div(timeseries.sum()).reset_index(drop=True)
    timeseries.name = "demand_norm"

    # Check result
    np.testing.assert_almost_equal(
        timeseries.sum(),
        1,
        err_msg="Sum of normalized timeseries is not 1!",
    )

    # Write
    timeseries.to_csv(outfile)


def merge_demand_multiple_years(
    infiles: list,
    outfile: Path,
) -> None:
    """Merge demand from different years into one

    Parameters
    ----------
    infiles : list of Path
        CSV with each holding demand per municipality for one year
    outfile : pathlib.Path
        CSV with demand per municipality for all years

    Returns
    -------
    None
    """
    demand = pd.concat(
        [pd.read_csv(f, index_col="municipality_id") for f in infiles],
        axis=1,
    )
    demand.to_csv(outfile)


def disaggregate_demand_to_municipality(
    demand_districts: pd.DataFrame,
    muns: gpd.GeoDataFrame,
    districts: gpd.GeoDataFrame,
    disagg_data: pd.DataFrame,
    disagg_data_col: str,
) -> pd.DataFrame:
    """Disaggregates NUTS 3 consumption to municipalities by linear scaling
    using factors from `disagg_data`.

    Parameters
    ----------
    demand_districts : pd.DataFrame
        Annual demand per district (NUTS 3 level)
    muns : gpd.GeoDataFrame
        Municipalities
    districts : gpd.GeoDataFrame
        Federal districts
    disagg_data : pd.DataFrame
        DF that contains municipal data to scale demand by
    disagg_data_col : str
        Name of DF column used for scaling

    Returns
    -------
    pd.DataFrame
        Annual demand per municipality
    """
    # Get muns and their NUTS3 code
    muns = muns.merge(
        districts[["id", "nuts"]].rename(columns={"id": "district_id"}),
        left_on="district_id",
        right_on="district_id",
    )
    disagg_data = (
        disagg_data.reset_index()
        .merge(muns[["id", "nuts"]], left_on="municipality_id", right_on="id")
        .drop(columns=["id"])
    )

    # Calc relative share within each district
    disagg_data = disagg_data.assign(
        pop_share=(
            disagg_data[disagg_data_col]
            / disagg_data.groupby("nuts")[disagg_data_col].transform("sum")
        )
    ).set_index("nuts")

    # Filter demand dataset
    demand_districts = demand_districts.loc[districts.nuts.to_list()]

    # Merge demand and calc value per municipality
    demand = disagg_data.merge(
        demand_districts,
        left_index=True,
        right_index=True,
    )
    demand[disagg_data_col] = demand_districts.demand_districts.mul(
        demand.pop_share
    )
    demand = demand[["municipality_id", disagg_data_col]].set_index(
        "municipality_id"
    )

    # Check result
    np.testing.assert_almost_equal(
        demand.sum().sum(),
        demand_districts.sum().sum(),
        err_msg="Sum of disaggregated values does not match total value!",
    )

    return demand


def demand_prognosis(
    demand_future_T45: Path,
    demand_future_TN: Path,
    demand_region: pd.DataFrame,
    year_base: int,
    year_target: int,
    scale_by: str,
    carrier: str = None,
) -> pd.DataFrame:
    """Create demand prognosis for target year per municipality by scaling with
    country data (by total demand or carrier)

    Parameters
    ----------
    demand_future_T45 : pathlib.Path
        Path to total demand per carrier in Germany (BMWK T45 scenario)
    demand_future_TN : pathlib.Path
        Path to total demand per carrier in Germany (BMWK TN scenario)
    demand_region : pd.DataFrame
        Annual demand per municipality (today)
    year_base : int
        Base year to use for calculation of reduction
    year_target : int
        Target year for demand
    scale_by : str
        One of ["total", "carrier"], scale future demand by total energy or
        specific carrier given by `carrier`
    carrier : str
        Carrier


    Returns
    -------
    pd.DataFrame
        Annual demand per municipality (prognosis)
    """
    if scale_by not in ["total", "carrier"]:
        raise ValueError("scale_by must be one of ['total', 'carrier']")
    if scale_by == "carrier":
        if carrier is None:
            raise ValueError("carrier must be provided when scale_by='carrier'")

    # Get data from future scenarios
    demand_future_T45 = pd.read_csv(demand_future_T45).set_index(
        ["year", "carrier"]
    )
    demand_future_TN = pd.read_csv(demand_future_TN)

    # Interpolate for base year
    demand_future_TN = demand_future_TN.set_index(["year", "carrier"]).append(
        pd.DataFrame(
            index=pd.MultiIndex.from_product(
                [[year_base], demand_future_TN.carrier.unique(), [np.nan]],
                names=demand_future_TN.columns,
            )
        )
    )
    demand_future_TN.sort_index(inplace=True)
    demand_future_TN = demand_future_TN.unstack(level=1).interpolate().stack()

    # Scale by total or carrier demand
    if scale_by == "total":
        # Total reduction factor
        reduction_factor = (
            demand_future_T45.loc[year_target].demand.sum()
            / demand_future_TN.loc[year_base].demand.sum()
        )
    elif scale_by == "carrier":
        # Reduction factor for carrier
        reduction_factor = (
            demand_future_T45.loc[year_target, carrier].demand
            / demand_future_TN.loc[year_base, carrier].demand
        )
    print(
        f"Total demand in {year_target}: {round(reduction_factor, 2)} of "
        f"{year_base}'s demand."
    )

    return demand_region * reduction_factor
