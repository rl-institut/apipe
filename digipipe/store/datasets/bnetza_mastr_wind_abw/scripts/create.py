import pandas as pd
import geopandas as gpd
from digipipe.scripts.datasets import mastr
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

    units = mastr.add_voltage_level(
        units=units,
        locations_path=snakemake.input.locations,
        gridconn_path=snakemake.input.gridconn
    )

    # Drop units outside of Germany
    units = units.loc[(units.country == "Deutschland")]

    # Add geometry and drop units without coords
    units = mastr.add_geometry(units)

    # Do some basic filtering
    units = mastr.cleanse(units)

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
