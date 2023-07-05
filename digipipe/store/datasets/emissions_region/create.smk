"""
Snakefile for this dataset

Note: To include the file in the main workflow, it must be added to the respective module.smk .
"""
import json
import pandas as pd
from digipipe.store.utils import get_abs_dataset_path

DATASET_PATH = get_abs_dataset_path(
    "datasets", "emissions_region", data_dir=True)

rule create:
    """
    Extract and aggregate emissions
    """
    input:
        get_abs_dataset_path(
            "raw", "emissions") / "data" / "emissions.csv"
    output:
        DATASET_PATH / "emissions.json"
    run:
        emissions = pd.read_csv(
            input[0],
            header=[0, 1],
            index_col=[0, 1, 2, 3],
        )

        emissions_dict = dict(
            # Sector: Energy industry (1.A.1)
            **emissions.loc[
                (1, "A", "1", "energy_industry")
            ].to_frame(name="energy_industry").to_dict(orient="list"),

            # Sector: Industry (1.A.2 + 2)
            **emissions.loc[[
                (1, "A", "2", "industry"),
                (2, "*", "*", "process_emissions"),
            ]].sum().to_frame(name="industry").to_dict(orient="list"),

            # Sector: Traffic (1.A.3)
            **emissions.loc[
                (1, "A", "3", "traffic"),
            ].to_frame(name="traffic").to_dict(orient="list"),

            # Sector: Buildings (1.A.4 and 1.A.5)
            **emissions.loc[
                (1, "A", "4-5", "buildings_firing"),
            ].to_frame(name="buildings_firing").to_dict(orient="list"),

            # Sector: Agriculture (3)
            **emissions.groupby(
                "sector"
            ).sum().loc[3].to_frame(name="agricultural").to_dict(orient="list"),

            # Sector: Waste and waste water
            **emissions.loc[
                (5, "*", "*", "waste_waste_water"),
            ].to_frame(name="waste_waste_water").to_dict(orient="list"),
        )

        with open(output[0], "w", encoding="utf8") as f:
            json.dump(emissions_dict, f, indent=4)
