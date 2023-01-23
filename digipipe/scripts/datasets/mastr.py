"""Shared functions for processing data from MaStR dataset"""

import geopandas as gpd
import pandas as pd
from geopy.extra.rate_limiter import RateLimiter
from geopy.geocoders import Nominatim
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
    units_df: pd.DataFrame,
    locations_path: str,
    gridconn_path: str,
) -> pd.DataFrame:
    """Add voltage level to units from MaStR using locations and grid
    connection points.

    Parameters
    ----------
    units_df : pd.DataFrame
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
    units_df = units_df.merge(
        locations[["mastr_location_id2", "voltage_level"]],
        left_on="mastr_location_id",
        right_on="mastr_location_id2",
        how="left",
    )

    # Drop unnecessary columns
    units_df.drop(
        columns=[
            "mastr_location_id",
            "mastr_location_id2",
        ],
        inplace=True,
    )

    return units_df


def add_geometry(
        units_df: pd.DataFrame,
        drop_units_wo_coords: bool = True,
) -> gpd.GeoDataFrame:
    """
    Add `geometry` column to MaStR unit data using `lat` and `lon` values.
    Columns `lat` and `lon` are dropped

    Parameters
    ----------
    units_df : pd.DataFrame
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
        units_count_orig = len(units_df)
        units_df = units_df.loc[(~units_df.lon.isna() & ~units_df.lat.isna())]
        print(
            f"{units_count_orig-len(units_df)} units have no or invalid "
            f"coordinates."
        )

    units_gdf = gpd.GeoDataFrame(
        units_df,
        geometry=gpd.points_from_xy(
            units_df["lon"], units_df["lat"], crs=4326
        ),
        crs=4326,
    ).to_crs(3035)

    # Drop unnecessary columns
    units_gdf.drop(columns=["lon", "lat"], inplace=True)

    return units_gdf


def geocode(
        units_df: pd.DataFrame,
        user_agent: str = "geocoder",
        interval: int = 1,
        target_crs: str = "EPSG:3035",
) -> gpd.GeoDataFrame:
    """
    Geocode locations from MaStR unit table using zip code and city.

    Parameters
    ----------
    units_df : pd.DataFrame
        Units from MaStR. Must contain the following columns:
        * zip_code (str)
        * city (str)
        *
    user_agent : str
        Some app name. Defaults to "geocoder"
    interval : int
        Delay in seconds to use between requests to Nominatim.
        A minimum of 1 is advised (default)
    target_crs : str
        CRS the data should be reprojected to. Defaults to EPSG:3035.

    Returns
    -------
    gpd.GeoDataFrame
        Units with geometry
    """

    def geocoder(
            user_agent: str,
            interval: int,
    ) -> RateLimiter:
        """Setup Nominatim geocoding class.

        Parameters
        -----------
        user_agent : str
            Some app name.
        interval : int
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
            min_delay_seconds=interval,
        )

    # Define geocoder
    ratelimiter = geocoder(
        user_agent,
        interval,
    )

    # Merge zip code and city and get unique values
    units_df = units_df.assign(
        zip_and_city=units_df.zip_code.astype(str) + " " + units_df.city,
    )
    unique_locations = pd.DataFrame(
        data=units_df.zip_and_city.unique(),
        columns=["zip_and_city"],
    )
    # Geocode unique locations!
    print(
        f"Geocoding {len(unique_locations)} unique locations, this will take "
        f"about {round(len(unique_locations) * interval / 60, 1)} min..."
    )
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
    units_gdf = gpd.GeoDataFrame(
        units_df.merge(
            unique_locations[["zip_and_city", "geometry"]],
            on="zip_and_city",
        ).drop(columns=["zip_and_city"])
    ).to_crs(target_crs)

    return units_gdf
