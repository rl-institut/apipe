import pandas as pd
import geopandas as gpd
from digipipe.scripts.datasets import mastr
from digipipe.scripts.geo import (
    write_geofile,
    overlay,
    rename_filter_attributes
)

from digipipe.config import GLOBAL_CONFIG


def process() -> None:
    attrs = snakemake.config["attributes"]
    attrs_filter = snakemake.config["attributes_filter"]

    units = pd.read_csv(
        snakemake.input.units,
        usecols=set(attrs.keys()) | set(attrs_filter.keys()),
        dtype={"Postleitzahl": str},
    )

    units = rename_filter_attributes(
        gdf=units,
        attrs_filter_by_values=attrs_filter,
        attrs_mapping=attrs,
    ).set_index("mastr_id")

    units = mastr.add_voltage_level(
        units=units,
        locations_path=snakemake.input.locations,
        gridconn_path=snakemake.input.gridconn
    )

    # Add geometry and drop units without coords and
    # add column to indicate that location from original data was used
    units_with_geom = mastr.add_geometry(units)
    units_with_geom = units_with_geom.assign(
        geometry_approximated=False,
    )

    # Add geometry for all units without coords (<=30 kW) and
    # add column to indicate that location was inferred by geocoding
    units_with_inferred_geom = mastr.geocode(
        units.loc[(units.lon.isna() | units.lat.isna())].drop(
            columns=["lon", "lat"]
        ),
        user_agent=GLOBAL_CONFIG["global"]["geodata"]["geocoder"]["user_agent"],
        interval=GLOBAL_CONFIG["global"]["geodata"]["geocoder"]["interval_sec"],
    )
    units_with_inferred_geom = units_with_inferred_geom.assign(
        geometry_approximated=True,
    )

    units = pd.concat([units_with_geom, units_with_inferred_geom])

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
