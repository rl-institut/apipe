"""Shared functions for processing data from MaStR dataset"""

import geopandas as gpd
import pandas as pd
from typing import Union


def cleanse(
    units: Union[pd.DataFrame, gpd.GeoDataFrame]
) -> Union[pd.DataFrame, gpd.GeoDataFrame]:
    """Do some basic cleansing of MaStR unit data.

    This involves:
    * country

    Parameters
    ----------
    units : pd.DataFrame or gpd.GeoDataFrame
        Units from MaStR

    Returns
    -------
    pd.DataFrame or gpd.GeoDataFrame
        Filtered units
    """
    units = units.loc[(units.country == "Deutschland")]

    # Drop unnecessary columns
    units.drop(
        columns=[
            "country",
        ],
        inplace=True,
    )

    return units


def add_voltage_level(
    units: pd.DataFrame,
    locations_path: str,
    gridconn_path: str,
) -> pd.DataFrame:
    """Add voltage level to units from MaStR using locations and grid
    connection points.

    Parameters
    ----------
    units : pd.DataFrame
        Units from MaStR
    locations_path : str
        Path to MaStR locations file
    gridconn_path : str
        Path to MaStR grid connection points file

    Returns
    -------
    pd.DataFrame
        Units with column `voltage_level`
    """
    # Get locations and grid connection point and merge both
    locations = pd.read_csv(
        locations_path,
        usecols=["MastrNummer", "Netzanschlusspunkte"],
    ).rename(columns={"MastrNummer": "mastr_location_id2"})
    gridconn = pd.read_csv(
        gridconn_path,
        usecols=["NetzanschlusspunktMastrNummer", "Spannungsebene"]
    )
    locations = locations.merge(
        gridconn,
        left_on="Netzanschlusspunkte",
        right_on="NetzanschlusspunktMastrNummer",
        how="left",
    ).drop_duplicates().rename(columns={"Spannungsebene": "voltage_level"})

    # Add voltage level to units
    units = units.merge(
        locations[["mastr_location_id2", "voltage_level"]],
        left_on="mastr_location_id",
        right_on="mastr_location_id2",
        how="left",
    )

    # Drop unnecessary columns
    units.drop(
        columns=[
            "mastr_location_id",
            "mastr_location_id2",
        ],
        inplace=True,
    )

    return units


def add_geometry(
        units: pd.DataFrame,
        drop_units_wo_coords: bool = True,
) -> gpd.GeoDataFrame:
    """
    Add `geometry` column to MaStR unit data using `lat` and `lon` values.
    Columns `lat` and `lon` are dropped

    Parameters
    ----------
    units : pd.DataFrame
        Units with columns `lat` and `lon` in CRS WGS84 (EPSG:4326)
    drop_units_wo_coords : bool
        Drop units which do not have valid lat and lon values.
        Note: coordinates in the MaStR are only provided for plants>30 kW.

    Returns
    -------
    gpd.GeoDataFrame
        Units with geometry in CRS LAEA Europe (EPSG:3035)
    """
    # Drop units without coords idf requested
    if drop_units_wo_coords:
        units_count_orig = len(units)
        units = units.loc[(~units.lon.isna() & ~units.lat.isna())]
        print(
            f"{units_count_orig-len(units)} units have been dropped as they "
            f"have no or invalid coordinates."
        )

    units = gpd.GeoDataFrame(
        units,
        geometry=gpd.points_from_xy(
            units["lon"], units["lat"], crs=4326
        ),
        crs=4326,
    ).to_crs(3035)

    # Drop unnecessaryw columns
    units.drop(columns=["lon", "lat"], inplace=True)

    return units
