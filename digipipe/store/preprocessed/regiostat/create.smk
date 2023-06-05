"""
Snakefile for this dataset

Note: To include the file in the main workflow, it must be added to the respective module.smk .
"""

from digipipe.store.utils import get_abs_dataset_path

DATASET_PATH = get_abs_dataset_path("preprocessed", "regiostat")

rule extract_demand_ind:
    """
    Extract energy demand for industry from Excel files and save to CSV
    """
    input: get_abs_dataset_path("raw", "regiostat") / "data" / "43531-01-02-4.xlsx"
    output:
        demand_states=DATASET_PATH / "data" / "demand_energy_industry_states.csv",
        demand_districts=DATASET_PATH/ "data" / "demand_energy_industry_districts.csv"
    script: DATASET_PATH / "scripts" / "create.py"
