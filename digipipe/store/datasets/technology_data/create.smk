"""
Snakefile for this dataset

Note: To include the file in the main workflow, it must be added to the respective module.smk .
"""

from digipipe.store.utils import get_abs_dataset_path

DATASET_PATH = get_abs_dataset_path(
    "datasets", "technology_data", data_dir=True)

rule copy_files:
    """
    Copy technology data from raw
    """
    input:
        get_abs_dataset_path(
            "raw", "technology_data") / "data" / "technology_data.json"
    output:
        DATASET_PATH / "technology_data.json"
    shell:
        """
        cp -p {input} {output}
        """
