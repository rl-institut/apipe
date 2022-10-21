"""
Snakefile for this dataset

Note: To include the file in the main workflow, it must be added to the respective module.smk .
"""

from digipipe.store.utils import get_abs_dataset_path, create_tag_string_ogr

rule extract_tags:
    """
    Create layers from converted OSM file using requested tags per layer. Only requested attributes are retained.
    """
    input: rules.preprocessed_osm_filtered_convert.output
    output: get_abs_dataset_path("datasets", "osm_forest") / "data" / "osm_forest.gpkg"
    params:
        tags=create_tag_string_ogr(config["tags"]),
        geom_type=config["geom_type"]
    run:
        shell(
            f"ogr2ogr -f GPKG -select {params.tags['tags']} "
            f"{params.tags['conditions']} {output} {input} {params.geom_type}"
        )
