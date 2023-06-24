"""
Snakefile for this dataset

Note: To include the file in the main workflow, it must be added to the respective module.smk .
"""

from digipipe.store.utils import get_abs_dataset_path

DATASET_PATH = get_abs_dataset_path("preprocessed", "rpg_abw_pv_roof_potential")

rule create:
    input:
        get_abs_dataset_path(
            "raw", "rpg_abw_pv_roof_potential"
        ) / "data" / "rpg_abw_pv_roof_potential.gpkg"
    output:
        DATASET_PATH / "data" / "rpg_abw_pv_roof_potential.gpkg"
    shell:
        """
        ogr2ogr -f GPKG -t_srs EPSG:3035 {output} {input}
        """
