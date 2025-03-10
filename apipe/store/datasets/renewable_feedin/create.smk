"""
Snakefile for this dataset

Note: To include the file in the main workflow, it must be added to the respective module.smk .
"""
import os
import pandas as pd
from apipe.store.utils import get_abs_dataset_path
from apipe.scripts.data_io import load_json

DATASET_PATH = get_abs_dataset_path("datasets", "renewable_feedin")

rule normalize_feedin_timeseries:
    """
    Normalize feedin timeseries and drop time index
    """
    input:
        get_abs_dataset_path(
            "raw", "renewables.ninja_feedin") / "data" /
            "{tech}_feedin_timeseries.csv",
        get_abs_dataset_path("raw", "technology_data") / "data" / "technology_data.json"
    output:
        DATASET_PATH / "data" / "{tech}_feedin_timeseries.csv",
    run:
        feedin_timeseries = pd.read_csv(input[0]).power
        feedin_timeseries = feedin_timeseries.div(feedin_timeseries.sum())

        # multiply with full load hours
        filename = output[0].split(os.sep)[-1]
        tech = "_".join(filename.split("_")[0:-2])
        tech_data = load_json(input[1])
        flh = tech_data["full_load_hours"][tech]["2045"]
        feedin_timeseries = feedin_timeseries.multiply(flh)

        feedin_timeseries.to_csv(output[0])
