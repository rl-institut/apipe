"""
Snakefile for this dataset

Note: To include the file in the main workflow, it must be added to the respective module.smk .
"""
import json
import os
import pandas as pd
from pathlib import Path

from apipe.store.utils import get_abs_dataset_path

DATASET_PATH = get_abs_dataset_path(
    "preprocessed", "bmwk_long_term_scenarios", data_dir=True)

rule create:
    """
    Extract files
    """
    input:
        get_abs_dataset_path(
            "raw", "bmwk_long_term_scenarios") / "data" /
            "bmwk_long_term_scenarios.zip"
    output:
        files=[DATASET_PATH / f for f in config["files_extract"]]
    params:
        outpath=DATASET_PATH,
        files_extract=" ".join(config["files_extract"])
    shell:
        """
        unzip -j {input} {params.files_extract} -d {params.outpath}
        """

rule rename_columns_carriers:
    """
    Rename columns, carriers and technologies
    """
    input: DATASET_PATH / "{file}.csv"
    output: DATASET_PATH / "{file}_reformatted.csv"
    run:
        data = pd.read_csv(input[0])

        # Rename columns
        data.rename(
            columns=config["rename_columns"].get(wildcards.file),
            inplace=True
        )
        if not all(
            c in config["rename_columns"].get(wildcards.file).values()
            for c in data.columns
        ):
            raise ValueError("At least one column was not renamed!")
        # Rename carriers and technologies
        data = data.replace(config["rename_entries"])

        data.to_csv(output[0], index=False)
        #os.remove(input[0])

rule create_captions:
    """
    Create attribute captions for app
    """
    input: [DATASET_PATH / f for f in config["files_extract"]]
    output: DATASET_PATH / "bmwk_long_term_scenarios_attribute_captions.json"
    run:
        captions = dict(
            datasets_caption_map={
                Path(f).stem: "bmwk_long_term_scenarios" for f in input},
            captions={"bmwk_long_term_scenarios": {
                v: k for k, v in config["rename_entries"].items()}
            },
        )
        with open(output[0], "w", encoding="utf8") as f:
            json.dump(captions, f, indent=4)
