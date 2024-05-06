"""
Snakefile for this dataset

Note: To include the file in the main workflow, it must be added to the respective module.smk .
"""

import pandas as pd

from apipe.store.utils import get_abs_dataset_path


DATASET_PATH = get_abs_dataset_path("preprocessed", "destatis_pop_prog")

rule create:
    """
    Extract municipality population data from Excel files and save to CSVs
    """
    input: get_abs_dataset_path("raw", "destatis_pop_prog") / "data" / "12421-0003-DLAND_$F.csv"
    output: DATASET_PATH / "data" / "population_prognosis_federal_states.csv"
    run:
        data = pd.read_csv(
            input[0], skiprows=6, encoding='iso-8859-1',
            delimiter=";", decimal=","
        )
        data = data.loc[data["Unnamed: 0"] == config["scenario"]]
        data = data.drop(
            columns=["Unnamed: 0", "Unnamed: 1"]
        ).rename(columns={"Unnamed: 2": "federal_state"}).set_index("federal_state")
        data.columns = [int(c[6:]) for c in data.columns]
        data = data * 1e3

        if config["years"]:
            data = data[config["years"]]

        data.to_csv(output[0])
