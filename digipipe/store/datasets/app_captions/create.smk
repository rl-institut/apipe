"""
Snakefile for this dataset

Note: To include the file in the main workflow, it must be added to the respective module.smk .
"""
import json

from digipipe.scripts.data_io import load_json
from digipipe.store.utils import get_abs_dataset_path

DATASET_PATH = get_abs_dataset_path("datasets", "captions", data_dir=True)

rule create_captions:
    """
    Create captions for fields in app
    """
    input:
        [
            rules.datasets_bnetza_mastr_captions_create.output.outfile,
            rules.datasets_potentialarea_wind_region_create_captions.output[0],
            rules.datasets_demand_heat_region_create_captions.output[0],
        ]
    output:
        DATASET_PATH / "captions_fields.json"
    run:
        captions = dict(
            datasets_caption_map={},
            captions={},
        )

        for file in input:
            captions_file = load_json(file)
            captions["datasets_caption_map"].update(
                captions_file["datasets_caption_map"]
            )
            captions["captions"].update(captions_file["captions"])

        with open(output[0], "w", encoding="utf8") as f:
            json.dump(captions, f, indent=4)
