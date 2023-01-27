import geopandas as gpd
from digipipe.scripts.geo import (
    convert_to_multipolygon,
    write_geofile,
    rename_filter_attributes,
    reproject_simplify,
    overlay
)


def process():
    muns = gpd.read_file(snakemake.input.abw_muns[0], layer=config["layer"])
    muns = rename_filter_attributes(
        gdf=muns,
        attrs_filter_by_values=config["attributes_filter"],
        attrs_mapping=config["attributes"],
    )
    muns = reproject_simplify(
        gdf=muns,
        add_id_column=True,
    )

    muns = overlay(
        gdf=muns.rename(columns={"id": "mun_id"}),
        gdf_overlay=gpd.read_file(snakemake.input.abw_districts),
        retain_rename_overlay_columns={"id": "district_id"},
        gdf_use_centroid=True,
    ).rename(columns={"mun_id": "id"})

    muns = convert_to_multipolygon(muns)

    write_geofile(
        gdf=muns,
        file=snakemake.output[0],
        layer_name=config["layer"],
    )


if __name__ == "__main__":
    config = snakemake.config
    process()
