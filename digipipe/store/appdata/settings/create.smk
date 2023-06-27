"""
Snakefile for this dataset

Note: To include the file in the main workflow, it must be added to the respective module.smk .
"""
import json

from digipipe.store.utils import get_abs_dataset_path

DATASET_PATH = get_abs_dataset_path("appdata", "settings", data_dir=True)

rule create_map_panel_layer_list:
    """
    Create layer list for right map panel
    """
    input: rules.appdata_datapackage_create_datapackage.output
    output: DATASET_PATH / "map_panel_layer_list.json"
    run:
        print("Creating list of layers...")
        with open(output[0], "w", encoding="utf8") as f:
            json.dump(config["map_panel_layer_list"], f, indent=4)
