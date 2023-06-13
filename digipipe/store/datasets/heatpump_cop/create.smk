"""
Snakefile for this dataset

Note: To include the file in the main workflow, it must be added to the respective module.smk .
"""
import pandas as pd
from digipipe.store.utils import get_abs_dataset_path

DATASET_PATH = get_abs_dataset_path("datasets", "heatpump_cop")

rule create:
    """
    Calculate COPs for air- and ground source heatpumps
    """
    input:
        temperature=rules.preprocessed_dwd_temperature_create.output
    output:
        cop_ashp=DATASET_PATH / "data" / "heatpump_cop_ashp.csv",
        cop_gshp=DATASET_PATH / "data" / "heatpump_cop_gshp.csv"
    script:
        DATASET_PATH / "scripts" / "create.py"

rule merge:
    """
    Create one heatpump COP timeseries by weighting both COP timeseries
    """
    input: rules.datasets_heatpump_cop_create.output
    output: DATASET_PATH / "data" / "heatpump_cop.csv"
    run:
        tech = config["heatpumps"].get("technology")
        cop = pd.concat(
            [
                pd.read_csv(f,  index_col=0)
                for f in input
            ],
            axis=1,
        )
        cop = pd.Series(
            cop.cop_ashp * tech.get("share_ASHP") +
            cop.cop_gshp * tech.get("share_GSHP"),
            name="cop"
        ).round(3)
        cop.to_csv(output[0])
