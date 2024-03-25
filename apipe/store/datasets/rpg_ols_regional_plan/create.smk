"""
Snakefile for this dataset

Note: To include the file in the main workflow, it must be added to the respective module.smk .
"""
from apipe.store.utils import get_abs_dataset_path


rule create_pv_ground_criteria:
    """
    Freiflächen-Photovoltaikanlagen Negativkriterien
    """
    input:
        get_abs_dataset_path(
            "preprocessed", "rpg_ols_regional_plan"
        ) / "data" / "{file}.gpkg"
    output:
        get_abs_dataset_path(
            "datasets", "rpg_ols_regional_plan"
        ) / "data" / "{file}.gpkg"
    shell: "cp -p {input} {output}"
