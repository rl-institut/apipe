"""
Snakefile for this dataset

Note: To include the file in the main workflow, it must be added to the respective module.smk .
"""
import json
import shutil

from apipe.store.utils import get_abs_dataset_path
from apipe.store.appdata.datapackage.scripts.create import collect_files

DATASET_PATH = get_abs_dataset_path("appdata", "datapackage", data_dir=True)
DATAPACKAGE_FILES = collect_files(config, DATASET_PATH)

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
