"""
Snakefile for this dataset

The file will be automatically detected and included in the Snakemake workflow.
"""

from digipipe.store.utils import get_abs_dataset_path, create_tag_string_ogr

configfile: get_abs_dataset_path("2_datasets", "osm_forest", data_dir=False) / "config.yml"

rule osm_forest_extract_tags:
    """
    Create layers from converted OSM file using requested tags per layer. Only requested attributes are retained.
    """
    input: get_abs_dataset_path("1_preprocessed", "osm_filtered") / "sachsen-anhalt-221003.osm.gpkg"
    output: get_abs_dataset_path("2_datasets", "osm_forest") / "osm_forest.gpkg"
    params:
        tags=create_tag_string_ogr(config["tags"]),
        geom_type=config["geom_type"]
    run:
        shell(
            f"ogr2ogr -f GPKG -select {params.tags['tags']} "
            f"{params.tags['conditions']} {output} {input} {params.geom_type}"
        )
