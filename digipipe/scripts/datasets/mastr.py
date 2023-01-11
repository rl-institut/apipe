"""Shared functions for processing data from MaStR dataset"""

import geopandas as gpd
import pandas as pd
from geopy.extra.rate_limiter import RateLimiter
from geopy.geocoders import Nominatim
from typing import Union

from digipipe.config import GLOBAL_CONFIG


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

    # Drop unnecessary columns
    units.drop(columns=["lon", "lat"], inplace=True)

    return units


def geocode(mastr_df: pd.DataFrame) -> gpd.GeoDataFrame:
    """
    Geocode locations from MaStR unit table using zip code and city.

    Parameters
    ----------
    mastr_df : pd.DataFrame
        Units from MaStR

    Returns
    -------
    gpd.GeoDataFrame
        Units with geometry
    """

    def geocoder(
            user_agent: str,
            min_delay_seconds: int,
    ) -> RateLimiter:
        """Setup Nominatim geocoding class.

        Parameters
        -----------
        user_agent : str
            The app name.
        min_delay_seconds : int
            Delay in seconds to use between requests to Nominatim.
            A minimum of 1 is advised.
        Returns
        -------
        geopy.extra.rate_limiter.RateLimiter
            Nominatim RateLimiter geocoding class to use for geocoding.
        """
        locator = Nominatim(user_agent=user_agent)
        return RateLimiter(
            locator.geocode,
            min_delay_seconds=min_delay_seconds,
        )

    # Define geocoder
    ratelimiter = geocoder(
        GLOBAL_CONFIG["global"]["geodata"]["geocoder"]["user_agent"],
        GLOBAL_CONFIG["global"]["geodata"]["geocoder"]["interval_sec"],
    )

    # Merge zip code and city and get unique values
    mastr_df = mastr_df.assign(
        zip_and_city=mastr_df.zip_code.astype(str) + " " + mastr_df.city,
    )
    unique_locations = pd.DataFrame(
        data=mastr_df.zip_and_city.unique(),
        columns=["zip_and_city"],
    )
    # Geocode unique locations!
    print(f"Geocoding {len(unique_locations)} locations...")
    unique_locations = unique_locations.assign(
        location=unique_locations.zip_and_city.apply(ratelimiter)
    )
    unique_locations = unique_locations.assign(
        point=unique_locations.location.apply(
            lambda loc: tuple(loc.point) if loc else None
        )
    )
    unique_locations[["latitude", "longitude", "altitude"]] = pd.DataFrame(
        unique_locations.point.tolist(), index=unique_locations.index
    )
    unique_locations = gpd.GeoDataFrame(
        unique_locations,
        geometry=gpd.points_from_xy(unique_locations.longitude,
                                    unique_locations.latitude),
        crs="EPSG:4326",
    )
    # Merge locations back in units
    mastr_df = gpd.GeoDataFrame(
        mastr_df.merge(
            unique_locations[["zip_and_city", "geometry"]],
            on="zip_and_city",
        ).drop(columns=["zip_and_city"])
    ).to_crs(GLOBAL_CONFIG["global"]["geodata"]["crs"])

    # Drop unnecessary columns
    mastr_df.drop(columns=["lon", "lat"], inplace=True)

    return mastr_df
