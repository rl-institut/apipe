"""
Snakefile for this dataset

Note: To include the file in the main workflow, it must be added to the respective module.smk .
"""

import pandas as pd
from digipipe.store.utils import get_abs_dataset_path

DATASET_PATH = get_abs_dataset_path("preprocessed", "dwd_temperature", data_dir=True)

rule create:
    """
    Calc mean temperatures for entire region
    """
    input:
        temperature=get_abs_dataset_path(
            "raw", "dwd_temperature") / "data" / "temperature_2011.csv"
    output:
        temperature=DATASET_PATH / "temperature_2011.csv"
    run:
        temp = pd.read_csv(
            input.temperature,
            index_col=["timestamp", "ags_id"],
        )
        # Calc regional mean from municipal values
        temp = temp.astype("float").reset_index().drop(
            columns=["ags_id"]).groupby("timestamp").agg("mean").round(2)
        # Remove timestamp
        temp = temp.reset_index().drop(columns=["timestamp"])
        # Dump
        temp.to_csv(output.temperature)
