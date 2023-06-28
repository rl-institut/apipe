"""
Snakefile for this dataset

Note: To include the file in the main workflow, it must be added to the respective module.smk .
"""

from digipipe.store.utils import get_abs_dataset_path

DATASET_PATH = get_abs_dataset_path(
    "datasets", "technology_data", data_dir=True)

rule copy_files:
    """
    Copy full load hours
    """
    input:
        flh=get_abs_dataset_path(
            "raw", "technology_data") / "data" / "full_load_hours.json"
    output:
        flh=DATASET_PATH / "full_load_hours.json"
    shell:
        """
        cp -p {input.flh} {output.flh}
        """
