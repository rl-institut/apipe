import geopandas as gpd
import pandas as pd
import numpy as np
from pathlib import Path


def normalize_filter_timeseries(
        infile: Path, outfile: Path, region_nuts: list = []
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
    timeseries = (
        timeseries[region_nuts].resample("H").sum()
    )
    # Average SLP timeseries and normalize to 1 MWh
    timeseries = timeseries.sum(axis=1)
    timeseries = timeseries.div(
        timeseries.sum()).reset_index(drop=True)
    timeseries.name = "consumption_norm"

    # Check result
    assert timeseries.sum() == 1, "Sum of normalized timeseries is not 1"

    # Write
    timeseries.to_csv(outfile)


def disaggregate_consumption_by_pop(
        infile: Path,
        outfile: Path,
        muns: gpd.GeoDataFrame,
        districts: gpd.GeoDataFrame,
        population: pd.DataFrame,
        year: int,
) -> None:
    """Disaggregates NUTS 3 consumption to municipalities using linear scaling
    by population for one year.

    Parameters
    ----------
    infile : Path
        CSV with consumption per NUTS 3 for specific year
    outfile : Path
        CSV with consumption per municipality for specific year
    muns : gpd.GeoDataFrame
        Municipalities
    districts : gpd.GeoDataFrame
        Federal districts
    population : pd.DataFrame
        Population for year
    year : int
        Target year for which the population data is used

    Returns
    -------

    """
    # Get muns and their NUTS3 code
    muns = muns.merge(
        districts[["id", "nuts"]].rename(columns={"id": "district_id"}),
        left_on="district_id",
        right_on="district_id"
    )
    population.columns = population.columns.droplevel(1)
    population = population[str(year)].reset_index().merge(
        muns[["id", "nuts"]],
        left_on="municipality_id",
        right_on="id"
    ).drop(columns=["id"])

    # Calc relative population share within each district
    population = population.assign(
        pop_share=(
                population[str(year)] /
                population.groupby("nuts")[str(year)].transform("sum")
        )
    ).set_index("nuts")

    consumption_nuts = pd.read_csv(infile).set_index("nuts3").loc[
        districts.nuts.to_list()].sum(axis=1).to_frame(name="consumption_nuts")

    # Merge consumption and calc value per municipality
    consumption = population.merge(
        consumption_nuts,
        left_index=True,
        right_index=True,
    )
    consumption[str(year)] = consumption_nuts.consumption_nuts.mul(
        consumption.pop_share)
    consumption = consumption[
        ["municipality_id", str(year)]].set_index("municipality_id")

    # Check result
    np.testing.assert_almost_equal(
        consumption.sum().sum(),
        consumption_nuts.sum().sum(),
        err_msg="Sum of disaggregated values does not match total value!"
    )

    # Write
    consumption.to_csv(outfile)


def merge_consumption_multiple_years(
        infiles: list,
        outfile: Path,
) -> None:
    """Merge consumption from different years into one

    Parameters
    ----------
    infiles : list of Path
        CSV with each holding consumption per municipality for one year
    outfile : Path
        CSV with consumption per municipality for all years

    Returns
    -------
    None
    """
    consumption = pd.concat(
        [
            pd.read_csv(f, index_col="municipality_id")
            for f in infiles
        ],
        axis=1,
    )
    consumption.to_csv(outfile)
