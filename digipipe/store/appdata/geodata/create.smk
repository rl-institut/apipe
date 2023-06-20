"""
Snakefile for this dataset

Note: To include the file in the main workflow, it must be added to the respective module.smk .
"""
import shutil
from pathlib import Path

from digipipe.store.utils import get_abs_dataset_path

DATASET_PATH = get_abs_dataset_path(
    "appdata", "geodata", data_dir=True)


def collect_geodata() -> list:
    """Collect paths of files for info layers

    Returns
    -------
    list
        Paths to original info layer files in store
    """
    infolayers = []
    for dataset, layers in config["geodata"].items():
        infolayers.extend(
            [get_abs_dataset_path("datasets", dataset, data_dir=True) / l
             for l in layers]
        )
    return infolayers

rule copy_files:
    """
    Copy geodata listet in config
    """
    input: collect_geodata()
    output:
        expand(
            DATASET_PATH / "{file}",
            file=[Path(f).name for f in collect_geodata()]
        )
    run:
        for file_in, file_out in zip(input, output):
            print(f"Copy file {Path(file_in).name}...")
            shutil.copyfile(file_in, file_out)
