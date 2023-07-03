"""
Snakefile for this dataset

Note: To include the file in the main workflow, it must be added to the respective module.smk .
"""

from digipipe.store.utils import (
    get_abs_dataset_path,
    create_tag_string_osmium,
    create_tag_string_ogr
)

DATASET_PATH = get_abs_dataset_path(
    "preprocessed", "osm_filtered", data_dir=True)

rule convert:
    """
    Convert OSM pbf file while filtering for all requested tags. All attributes are retained.
    """
    input: get_abs_dataset_path("raw", "osm") / "data" / "germany-230630.osm.pbf"
    #input: get_abs_dataset_path("raw", "osm") / "data" / "sachsen-anhalt-230630.osm.pbf"
    output: DATASET_PATH / "germany-230630_filtered.osm.gpkg"
    #output: DATASET_PATH / "sachsen-anhalt-230630_filtered.osm.gpkg"
    params:
        temp_file=DATASET_PATH / "germany-230630_filtered.osm.pbf",
        #temp_file=DATASET_PATH / "sachsen-anhalt-230630_filtered.osm.pbf",
        tags=create_tag_string_osmium(config["tags"])
    shell:
        """
        osmium tags-filter --remove-tags -f osm {input} {params.tags} -o {params.temp_file}
        ogr2ogr -f GPKG -t_srs EPSG:3035 {output} {params.temp_file}
        rm {params.temp_file}
        """

rule extract_buildings:
    """
    Create layers from converted OSM file using requested tags per layer. Only requested attributes are retained.
    """
    input: DATASET_PATH / "germany-230630_filtered.osm.gpkg"
    #input: DATASET_PATH / "sachsen-anhalt-230630_filtered.osm.gpkg"
    output: DATASET_PATH / "osm_buildings.gpkg"
    params:
        tags=create_tag_string_ogr(config["osm_buildings"]["tags"]),
        geom_type=config["osm_buildings"]["geom_type"]
    run:
        conditions = (
            params.tags['conditions']
            if config["osm_buildings"]["use_conditions"]
            else ""
        )
        shell(
            f"ogr2ogr -f GPKG -select {params.tags['tags']} "
            f"{conditions} {output} {input} {params.geom_type}"
        )
