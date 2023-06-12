"""
Snakefile for this dataset

Note: To include the file in the main workflow, it must be added to the respective module.smk .
"""
import json
import pandas as pd
from digipipe.store.utils import get_abs_dataset_path

DATASET_PATH = get_abs_dataset_path("datasets", "renewable_feedin")

rule rescale_timeseries:
    """
    Get raw feedin timeseries and scale to full load hours
    """
    input:
        raw_ts=get_abs_dataset_path(
            "raw", "renewables.ninja_feedin") / "data" /
            "{tech}_feedin_timeseries.csv",
        full_load_hours=get_abs_dataset_path(
            "raw", "renewables.ninja_feedin") / "data" / "full_load_hours.json"
    output:
        today_ts=DATASET_PATH / "data" / "{tech}_feedin_timeseries_today.csv",
        future_ts=DATASET_PATH / "data" / "{tech}_feedin_timeseries_future.csv"
    run:
        feedin_timeseries = pd.read_csv(input.raw_ts).power
        with open(input.full_load_hours, "r") as f:
            flh=json.load(f)["full_load_hours"].get(wildcards.tech)

        feedin_timeseries.mul(flh.get("today")).to_csv(output.today_ts)
        feedin_timeseries.mul(flh.get("future")).to_csv(output.future_ts)
