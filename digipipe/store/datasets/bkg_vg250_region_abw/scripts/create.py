import os
import sys
import geopandas as gpd
from digipipe.scripts.geo import (
    convert_to_multipolygon,
    write_geofile,
    reproject_simplify
)


def process():
    data = gpd.read_file(infile)
    data = gpd.GeoDataFrame(
        crs=data.crs.srs,
        geometry=[data.buffer(0.1).unary_union]
    )
    data = reproject_simplify(
        gdf=data,
        add_id_column=True
    )
    data = convert_to_multipolygon(data)
    write_geofile(
        gdf=data,
        file=outfile,
        layer_name=os.path.basename(outfile).split(".")[0],
    )


if __name__ == "__main__":
    infile = sys.argv[1]
    outfile = sys.argv[3]
    process()
