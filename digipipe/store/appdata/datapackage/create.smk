"""
Snakefile for this dataset

Note: To include the file in the main workflow, it must be added to the respective module.smk .
"""
import json
import shutil
from typing import Tuple

from digipipe.store.utils import get_abs_dataset_path

DATASET_PATH = get_abs_dataset_path("appdata", "datapackage", data_dir=True)

def collect_files() -> Tuple[list, list]:
    """Collect paths of source and target files for app datapackage

    Returns
    -------
    list
        Paths to original info layer files in store
    """
    source_files = []
    target_files = []

    for cat in config["resources"].keys():
        print(f"Processing {cat} ...")
        for subcat in config["resources"][cat].keys():
            print(f"  Processing {subcat} ...")
            for item, data in config["resources"][cat][subcat].items():
                print(f"  Processing {item} ...")
                source_file = (
                    get_abs_dataset_path(
                        "datasets",
                        data["_source_path"].get("dataset"),
                        data_dir=True
                    ) / data["_source_path"].get("file")
                )
                target_file = DATASET_PATH / data.get("path")
                if not target_file in target_files:
                    source_files.append(source_file)
                    target_files.append(target_file)
                else:
                    print("    Target file already collected, skipping...")

    return source_files, target_files

DATAPACKAGE_FILES = collect_files()

rule copy_files:
    """
    Copy required files
    """
    input: DATAPACKAGE_FILES[0]
    output: DATAPACKAGE_FILES[1]
    run:
        print("Copy required files for datapackage...")
        for src_file, dst_file in zip(*DATAPACKAGE_FILES):
            shutil.copy(src_file, dst_file)

rule create_datapackage:
    """
    Create datapackage for app
    """
    input: DATAPACKAGE_FILES[1]
    output: DATASET_PATH / "datapackage_app.json"
    run:
        print("Creating datapackage...")
        with open(output[0], "w", encoding="utf8") as f:
            json.dump(config, f, indent=4)
