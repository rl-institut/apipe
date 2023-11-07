"""
Snakefile for this dataset

Note: To include the file in the main workflow, it must be added to the respective module.smk .
"""

from apipe.store.utils import get_abs_dataset_path

DATASET_PATH = get_abs_dataset_path("preprocessed", "eurostat_lau")

rule create:
    """
    Extract Germany's LAUs from Excel files and save to CSV
    """
    input: get_abs_dataset_path("raw", "eurostat_lau") / "data" / "EU-27-LAU-2022-NUTS-2021.xlsx"
    output: DATASET_PATH / "data" / "germany_lau_codes.csv"
    script: DATASET_PATH / "scripts" / "create.py"
