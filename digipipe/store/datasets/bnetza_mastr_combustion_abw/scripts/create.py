import pandas as pd
import geopandas as gpd
from digipipe.scripts.datasets import mastr
from digipipe.scripts.geo import (
    write_geofile,
    overlay,
    rename_filter_attributes
)
from digipipe.store.utils import df_merge_string_columns

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
        geometry_approximated=0,
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
        geometry_approximated=1,
    )

    # Merge both DFs
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

    # Aggregate units with approximated position
    units_with_inferred_geom["lon"] = units_with_inferred_geom.geometry.x
    units_with_inferred_geom["lat"] = units_with_inferred_geom.geometry.y

    units_with_inferred_geom["fuel_primary"].fillna("", inplace=True)
    units_with_inferred_geom["fuel_secondary"].fillna("", inplace=True)

    units_with_inferred_geom_agg = (
        units_with_inferred_geom[
            ["zip_code", "city", "capacity_net", "capacity_gross",
             "th_capacity", "fuel_primary", "fuel_secondary", "lat", "lon"]
        ].groupby(["lat", "lon", "zip_code", "city"], as_index=False).agg({
            "capacity_net": ["sum", "count"],
            "capacity_gross": "sum",
            "th_capacity": "sum",
            "fuel_primary": ";".join,
            "fuel_secondary": ";".join,
        })
    )
    units_with_inferred_geom_agg.columns = [
        '_'.join(tup).rstrip('_') for tup in
        units_with_inferred_geom_agg.columns
    ]
    units_with_inferred_geom_agg = units_with_inferred_geom_agg.rename(
        columns={
            "capacity_net_sum": "capacity_net",
            "capacity_gross_sum": "capacity_gross",
            "th_capacity_sum": "th_capacity",
            "capacity_net_count": "unit_count",
            "fuel_primary_join": "fuel_primary",
            "fuel_secondary_join": "fuel_secondary",
        }
    )

    # Merge fuel types into one and make unique
    units_with_inferred_geom_agg["fuels"] = df_merge_string_columns(
        units_with_inferred_geom_agg[["fuel_primary", "fuel_secondary"]]
    )

    units_with_inferred_geom_agg = gpd.GeoDataFrame(
        units_with_inferred_geom_agg,
        geometry=gpd.points_from_xy(units_with_inferred_geom_agg.lon,
                                    units_with_inferred_geom_agg.lat),
        crs="EPSG:3035",
    )[[
        "zip_code", "city", "capacity_net", "capacity_gross", "th_capacity",
        "fuels", "unit_count", "geometry"
    ]]
    units_with_inferred_geom_agg = units_with_inferred_geom_agg.assign(
        status="In Betrieb oder in Planung",
        geometry_approximated=1,
    )
    units_agg = pd.concat([
        units_with_geom.assign(unit_count=1),
        units_with_inferred_geom_agg
    ])

    write_geofile(
        gdf=units,
        file=snakemake.output.outfile,
        layer_name=snakemake.config["layer"],
    )
    write_geofile(
        gdf=units_agg,
        file=snakemake.output.outfile_agg,
        layer_name=snakemake.config["layer"],
    )


process()
