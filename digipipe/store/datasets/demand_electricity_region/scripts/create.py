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
        timesteps in rows (expects 15min resolution)
    outfile : pathlib.Path
        Path to CSV outfile - aggregated and normalized timeseries
        (1h resolution)
    region_nuts : list of str
        List of NUTS 3 codes whose timeseries are to be used

    Returns
    -------
    None
    """
    # Create timeindex
    timeindex = pd.DatetimeIndex(
        pd.date_range(
            start="2017-01-01 00:00:00",
            end="2017-12-31 23:45:00",
            freq="15Min",
        )
    )
    # Get timeseries, filter region, resample to 1h
    timeseries = pd.read_csv(infile)
    timeseries.index = timeindex
    timeseries = timeseries[region_nuts].resample("H").sum()
    # Average SLP timeseries and normalize to 1 MWh
    timeseries = timeseries.sum(axis=1)
    timeseries = timeseries.div(timeseries.sum()).reset_index(drop=True)
    timeseries.name = "demand_norm"

    # Check result
    assert timeseries.sum() == 1, "Sum of normalized timeseries is not 1"

    # Write
    timeseries.to_csv(outfile)


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


def merge_demand_multiple_years(
    infiles: list,
    outfile: Path,
) -> None:
    """Merge demand from different years into one

    Parameters
    ----------
    infiles : list of Path
        CSV with each holding demand per municipality for one year
    outfile : Path
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


def demand_prognosis(
    demand_future_germany: Path,
    demand_districts: pd.DataFrame,
    demand_region: pd.DataFrame,
    year: int,
) -> pd.DataFrame:
    """Create demand prognosis for target year per municipality by scaling with
    country values but taking the municipal values into account.

    Parameters
    ----------
    demand_future_germany : pd.DataFrame
        Total demand per carrier in Germany
    demand_districts : pd.DataFrame
        Annual demand per district (NUTS 3 level)
    demand_region : pd.DataFrame
        Annual demand per municipality (today)
    year

    Returns
    -------
    pd.DataFrame
        Annual demand per municipality (prognosis)
    """

    demand_future_germany = pd.read_csv(
        demand_future_germany, index_col=0
    ).rename(
        columns={
            "Energiebedarf in TWh / Energy Demand in TWh": "demand",
            "Energieträger / Energy Carrier": "carrier",
        }
    )
    demand_future_germany = demand_future_germany.loc[
        demand_future_germany.carrier == "Strom"
    ].demand.to_frame()
    demand = (
        float(demand_future_germany.loc[year])
        * float(demand_region.sum())
        / demand_districts.sum().sum()
        * 1e6
        * (demand_region / demand_region.sum())
    ).rename(columns={"2022": year})
    return demand
