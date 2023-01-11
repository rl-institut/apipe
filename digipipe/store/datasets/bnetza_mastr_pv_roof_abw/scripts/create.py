import pandas as pd
import geopandas as gpd
from digipipe.scripts.datasets import mastr
from digipipe.scripts.geo import (
    write_geofile,
    overlay,
    rename_filter_attributes
)


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

    # Add geometry and drop units without coords
    units_with_geom = mastr.add_geometry(units)

    # Add geometry for all units without coords (<=30 kW)
    units_with_inferred_geom = mastr.geocode(
        units.loc[(units.lon.isna() | units.lat.isna())]
    )

    units = pd.concat([units_with_geom, units_with_inferred_geom])

    # Drop mun columns
    units.drop(
        columns=["municipality_name", "municipality_ags"], inplace=True
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
