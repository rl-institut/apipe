import pandas as pd
import geopandas as gpd
from digipipe.scripts.geo import (
    write_geofile,
    overlay
)


def process() -> None:
    columns_dict = snakemake.config["attributes"]

    units = pd.read_csv(
        snakemake.input.units,
        usecols=columns_dict.keys(),
        dtype={"Postleitzahl": str},
    ).rename(columns=columns_dict).set_index("mastr_id")

    # Drop units outside of Germany
    units = units.loc[(units.country == "Deutschland")]

    # Get locations and grid connection point and merge both
    locations = pd.read_csv(
        snakemake.input.locations,
        usecols=["MastrNummer", "Netzanschlusspunkte"],
    ).rename(columns={"MastrNummer": "mastr_location_id2"})
    gridconn = pd.read_csv(
        snakemake.input.gridconn,
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

    # Add geometry and drop units without coords
    units = units.loc[(~units.lon.isna() & ~units.lat.isna())]
    units = gpd.GeoDataFrame(
        units,
        geometry=gpd.points_from_xy(
            units["lon"], units["lat"], crs=4326
        ),
        crs=4326,
    ).to_crs(3035)
    units = units.rename(columns={"geometry": "geom"}).set_geometry("geom")

    # Drop unnecessary columns
    units.drop(
        columns=[
            "mastr_location_id",
            "mastr_location_id2",
            "lon",
            "lat",
            "country",
        ],
        inplace=True,
    )

    # Clip to ABW and add mun and district ids
    units = overlay(
        gdf=units,
        gdf_overlay=gpd.read_file(snakemake.input.abw_muns),
        retain_rename_overlay_columns={"id": "mun_id"}
    )
    units = overlay(
        gdf=units,
        gdf_overlay=gpd.read_file(snakemake.input.abw_districts),
        retain_rename_overlay_columns={"id": "district_id"}
    )

    write_geofile(
        gdf=units,
        file=snakemake.output.outfile,
        layer_name=snakemake.config["layer"],
    )


process()
