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
    units = mastr.add_geometry(units)

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
