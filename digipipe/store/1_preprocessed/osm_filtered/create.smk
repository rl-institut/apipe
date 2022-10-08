"""
Snakefile for this dataset

The file will be automatically detected and included in the Snakemake workflow.
"""

from digipipe.store.utils import get_abs_dataset_path, create_tag_string_osmium

configfile: get_abs_dataset_path("1_preprocessed", "osm_filtered", data_dir=False) / "config.yml"

rule osm_filtered_convert:
    """
    Convert OSM pbf file while filtering for all requested tags. All attributes are retained.
    """
    input: get_abs_dataset_path("0_raw", "osm_openstreetmap") / "sachsen-anhalt-221003.osm.pbf"
    output: get_abs_dataset_path("1_preprocessed", "osm_filtered") / "sachsen-anhalt-221003.osm.gpkg"
    params: tags=create_tag_string_osmium(config["tags"])
    shell:
        "osmium tags-filter --remove-tags -f osm {input} {params.tags} | "
        "ogr2ogr -f GPKG -t_srs EPSG:3035 {output} /vsistdin/"
