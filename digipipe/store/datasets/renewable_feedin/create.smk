"""
Snakefile for this dataset

Note: To include the file in the main workflow, it must be added to the respective module.smk .
"""
import pandas as pd
from digipipe.store.utils import get_abs_dataset_path

DATASET_PATH = get_abs_dataset_path("datasets", "renewable_feedin")

rule normalize_feedin_timeseries:
    """
    normalize feedin timeseries and drop time index
    """
    input:
        get_abs_dataset_path(
            "raw", "renewables.ninja_feedin") / "data" /
            "{tech}_feedin_timeseries.csv"
    output:
        DATASET_PATH / "data" / "{tech}_feedin_timeseries.csv",
    run:
        feedin_timeseries = pd.read_csv(input[0]).power
        feedin_timeseries = feedin_timeseries.div(feedin_timeseries.sum())
        feedin_timeseries.to_csv(output[0])

rule copy_full_load_hours:
    """
    copy full load hours
    """
    input:
        get_abs_dataset_path(
            "raw", "renewables.ninja_feedin") / "data" / "full_load_hours.json"
    output:
        DATASET_PATH / "data" / "full_load_hours.json"
    shell:
        """
        cp -p {input} {output}
        """
