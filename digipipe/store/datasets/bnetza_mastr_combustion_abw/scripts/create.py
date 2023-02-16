import pandas as pd
import geopandas as gpd
from digipipe.scripts.datasets import mastr
from digipipe.scripts.geo import (
    write_geofile,
    overlay,
    rename_filter_attributes
)
from digipipe.store.utils import df_merge_string_columns


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
        units_df=units,
        locations_path=snakemake.input.locations,
        gridconn_path=snakemake.input.gridconn
    )

    # Add geometry and drop units without coords and
    # add column to indicate that location from original data was used
    units_with_geom = mastr.add_geometry(units)
    units_with_geom = units_with_geom.assign(
        geometry_approximated=0,
    )

    units_without_geom = (
        units.loc[(units.lon.isna() | units.lat.isna())].drop(
            columns=["lon", "lat"])
        )

    # Add geometry for all units without coords (<=30 kW) and
    # add column to indicate that location was inferred by geocoding
    if len(units_without_geom) > 0:
        units_without_geom["fuel_primary"].fillna("", inplace=True)
        units_without_geom["fuel_secondary"].fillna("", inplace=True)

        units_with_inferred_geom_gdf, units_with_inferred_geom_agg_gdf = (
            mastr.geocode_units_wo_geometry(
                units_without_geom,
                columns_agg_functions={
                    "capacity_net": ("capacity_net", "sum"),
                    "unit_count": ("capacity_net", "count"),
                    "capacity_gross": ("capacity_gross", "sum"),
                    "th_capacity": ("th_capacity", "sum"),
                    "fuel_primary": ("fuel_primary", ";".join),
                    "fuel_secondary": ("fuel_secondary", ";".join),
                }
            )
        )

        # Merge fuel types into one and make unique
        # import pdb
        # pdb.set_trace()
        units_with_inferred_geom_agg_gdf["fuels"] = df_merge_string_columns(
            units_with_inferred_geom_agg_gdf[
                ["fuel_primary", "fuel_secondary"]
            ]
        )
        # units_with_inferred_geom_agg_gdf[
        #     ["fuel_primary", "fuel_secondary"]
        # ] = None

        # Merge both GDFs
        units = pd.concat(
            [units_with_geom, units_with_inferred_geom_gdf]
        ).drop(columns=["fuel_primary", "fuel_secondary"])
        units_agg = pd.concat([
            units_with_geom.assign(
                unit_count=1,
                fuels=(
                    units_with_geom["fuel_primary"].fillna("") + "; " +
                    units_with_geom["fuel_secondary"].fillna("")
                ),
            ),
            units_with_inferred_geom_agg_gdf
        ]).drop(columns=["fuel_primary", "fuel_secondary"])
    else:
        units = units_with_geom.drop(
            columns=["fuel_primary", "fuel_secondary"]
        )
        units_agg = units_with_geom.assign(
            unit_count=1,
            fuels=(
                    units_with_geom["fuel_primary"].fillna("") + "; " +
                    units_with_geom["fuel_secondary"].fillna("")
            ),
        ).drop(columns=["fuel_primary", "fuel_secondary"])

    # Clip to ABW and add mun and district ids
    units = overlay(
        gdf=units,
        gdf_overlay=gpd.read_file(snakemake.input.abw_muns),
        retain_rename_overlay_columns={"id": "municipality_id"}
    )
    units = overlay(
        gdf=units,
        gdf_overlay=gpd.read_file(snakemake.input.abw_districts),
        retain_rename_overlay_columns={"id": "district_id"}
    )
    units_agg = overlay(
        gdf=units_agg,
        gdf_overlay=gpd.read_file(snakemake.input.abw_muns),
        retain_rename_overlay_columns={"id": "municipality_id"}
    )
    units_agg = overlay(
        gdf=units_agg,
        gdf_overlay=gpd.read_file(snakemake.input.abw_districts),
        retain_rename_overlay_columns={"id": "district_id"}
    )

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
