"""
Helper functions for geodata processing
"""

import fiona
import os
import geopandas as gpd
from collections import OrderedDict
from shapely.geometry.multipolygon import MultiPolygon
from shapely.ops import transform
from typing import Tuple

from digipipe.scripts.config import read_config

from digipipe.config import GLOBAL_CONFIG


def read_schema_from_file(file: str) -> Tuple[str, OrderedDict]:
    """Read a geo file and returns schema definition using fiona

    Parameters
    ----------
    file : str
        Full path to file to read schema from

    Returns
    -------
    str
        Schema of geometry
    OrderedDict
        Properties/Fields of dataset (str: str)
    """
    try:
        with fiona.open(file) as f:
            schema_in_geom = f.schema["geometry"]
            schema_in_props = f.schema["properties"]
    except:
        f.close()
        raise
    return schema_in_geom, schema_in_props


def convert_to_multipolygon(
        gdf: gpd.GeoDataFrame
) -> gpd.GeoDataFrame:
    """Convert geometry column to type MultiPolygon

    Parameters
    ----------
    gdf : :pandas:`geopandas.GeoDataFrame`
        Data to be aligned

    Returns
    -------
    :pandas:`geopandas.GeoDataFrame`
    """

    def remove_z(row):
        """
        Remove z coordinate from Geometry, e.g. MultiPolygon (MULTIPOLYGON Z)
        """
        return transform(lambda x, y, z=None: (x, y), row)

    gdf["geometry"] = [MultiPolygon([feature])
                       if feature.type == "Polygon"
                       else feature
                       for feature in gdf["geometry"]]

    gdf["geometry"] = gdf["geometry"].apply(remove_z)

    return gdf


def write_geofile(
        gdf: gpd.GeoDataFrame,
        file: str,
        layer_name: str = None,
        schema: dict = None,
        driver: str = "GPKG",
        encoding: str = "utf-8"
) -> None:
    """Write geodata to file

    Parameters
    ----------
    gdf : :pandas:`geopandas.GeoDataFrame`
    file : str
        Target file
    layer_name : str
        Name of layer, usually same as file basename
    schema : dict
        Output schema with keys "geometry" and "properties"
    driver : str
        Geofile driver, default is Geopackage
    encoding : str
        Encoding
    """
    if layer_name is None:
        layer_name = os.path.basename(file).split(".")[0]
    if schema is None:
        # TODO: Log warning
        pass

    # check if data contain multiple geometry types
    if len(gdf.geometry.type.unique()) > 1:
        types = gdf.geometry.type.unique()
        raise ValueError(f"Data contain multiple geometry types: {types} !")

    gdf.to_file(file,
                layer=layer_name,
                schema=schema,
                driver=driver,
                encoding=encoding)


def reproject_simplify_filter_rename(
        gdf: gpd.GeoDataFrame,
        attrs_filter_by_values: dict = None,
        attrs_mapping: dict = None,
        target_crs: str = GLOBAL_CONFIG["geodata"]["crs"].lower(),
        min_size: float = None,
        simplify_tol: float = None,
        fix_geom: bool = False,
        add_id_column: bool = False,
) -> gpd.GeoDataFrame:
    """General purpose function for processing of geodata

    Parameters
    ----------
    gdf : gpd.GeoDataFrame
        Geodata
    attrs_filter_by_values : dict
        Attributes whose values are to be filtered. Use attributes as dict
        keys and desired values as dict values (values can be of type str,
        int, float or list)
    attrs_mapping : dict
        Attributes to select and rename. Use original attributes' names as
        dict keys and new names as values.
    target_crs : str
        CRS the data should be reprojected to.
    min_size : float
        Min. size of area to select (in sqm). Use None for no filtering
        (default).
        Raises exception if `target_crs` is not LAEA Europe (EPSG:3035).
    simplify_tol : float
        Threshold for simplification of geometries. Use None for no
        simplification (default).
        Raises exception if `target_crs` is not LAEA Europe (EPSG:3035).
    fix_geom : bool
        If True, invalid geometries are fixed by buffering by the value
        specified in the global config (geodata -> fix_geom_buffer).
    add_id_column : bool
        If True, data is reindexed starting from 0 and a new column "id" is
        added with the same values.

    Returns
    -------
    gpd.GeoDataFrame
        Processed geodata
    """

    def check_crs(operation: str) -> None:
        """Check if requested CRS is LAEA Europe (EPSG:3035)"""
        if target_crs.lower() != "epsg:3035":
            raise ValueError(
                f"Cannot apply {operation} in non-equistant CRS "
                f"(requested CRS: {target_crs.lower()}) !"
            )

    # Transform to target CRS
    if gdf.crs is not None:
        if str(gdf.crs).lower() != target_crs.lower():
            gdf = gdf.to_crs(target_crs)
    else:
        raise ValueError("Geodata has not CRS assigned.")

    # Filter by min size
    if min_size is not None:
        check_crs("min size filtering")
        gdf = gdf[gdf.area > min_size]

    # Generalize
    if simplify_tol is not None:
        check_crs("simplification")
        gdf["geometry"] = gdf.simplify(
            simplify_tol,
            preserve_topology=True
        )

    # Fix invalid geometries
    if fix_geom is True:
        buffer = GLOBAL_CONFIG["geodata"]["fix_geom_buffer"]
        if buffer > 0:
            gdf["geometry"] = gdf.buffer(buffer)

    # Filter by attribute values, if defined
    if attrs_filter_by_values is not None:
        query = ""
        for k, v in attrs_filter_by_values.items():
            if isinstance(v, list):
                query += f" & {k} in @v"
            elif isinstance(v, (str, int, float)):
                query += f" & {k}=={v}"
            else:
                raise ValueError(
                    "Data type in attribute filter is not supported!"
                )
        query = query[2:]
        gdf = gdf.query(query)

    # Extract and rename fields
    if attrs_mapping is not None:
        gdf = gdf.filter(attrs_mapping.keys())
        gdf.rename(columns=attrs_mapping, inplace=True)

    # Reindex starting from 0 and add new column "id" with same values
    if add_id_column is True:
        gdf.reset_index(drop=True, inplace=True)
        gdf = gdf.assign(id=gdf.index)

    return gdf
