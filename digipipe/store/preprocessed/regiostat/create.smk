"""
Snakefile for this dataset

Note: To include the file in the main workflow, it must be added to the respective module.smk .
"""

from digipipe.store.utils import get_abs_dataset_path
from digipipe.store.preprocessed.regiostat.scripts import create

DATASET_PATH = get_abs_dataset_path("preprocessed", "regiostat")

rule extract_demand_ind:
    """
    Extract energy demand for industry from Excel files and save to CSV
    """
    input: get_abs_dataset_path("raw", "regiostat") / "data" / "43531-01-02-4.xlsx"
    output:
        demand_states=DATASET_PATH / "data" / "demand_energy_industry_states.csv",
        demand_districts=DATASET_PATH/ "data" / "demand_energy_industry_districts.csv"
    run:
        create.extract_demand_ind(
            infile=input[0],
            outfile_states=output.demand_states,
            outfile_districts=output.demand_districts,
            cfg=config
        )

rule extract_employment_ind:
    """
    Extract employment for industry from Excel files and save to CSV
    """
    input: get_abs_dataset_path("raw", "regiostat") / "data" / "42111-01-04-5.xlsx"
    output:
        employment_states=DATASET_PATH / "data" / "employment_industry_states.csv",
        employment_districts=DATASET_PATH/ "data" / "employment_industry_districts.csv",
        employment_muns=DATASET_PATH/ "data" / "employment_industry_muns.csv"
    run:
        create.extract_employment_ind(
            infile=input[0],
            outfile_states=output.employment_states,
            outfile_districts=output.employment_districts,
            outfile_muns=output.employment_muns,
            cfg=config
        )
