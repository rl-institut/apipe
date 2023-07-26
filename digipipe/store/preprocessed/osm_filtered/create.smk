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
    input: get_abs_dataset_path("raw", "osm") / "data" / "germany-230101.osm.pbf"
    output: DATASET_PATH / "germany-230101_filtered.osm.gpkg"
    params:
        temp_file=DATASET_PATH / "germany-230101_filtered.osm.pbf",
        tags=create_tag_string_osmium(config["tags"])
    shell:
        """
        osmium tags-filter --remove-tags -f osm {input} {params.tags} -o {params.temp_file}
        ogr2ogr -f GPKG -t_srs EPSG:3035 {output} {params.temp_file}
        rm {params.temp_file}
        """
