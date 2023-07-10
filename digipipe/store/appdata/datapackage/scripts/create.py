from pathlib import Path
from typing import Tuple

from digipipe.store.utils import get_abs_dataset_path


def collect_files(
    config: dict,
    dataset_path: Path,
) -> Tuple[list, list]:
    """Collect paths of source and target files for app datapackage

    Parameters
    ----------
    config : dict
        Dataset config
    dataset_path: pathlib.Path
        Path to datapackage

    Returns
    -------
    list
        Source files
    list
        Target files
    """
    source_files = []
    target_files = []

    for cat in config["resources"].keys():
        print(f"Processing {cat} ...")
        for subcat in config["resources"][cat].keys():
            print(f"  Processing {subcat} ...")
            for item, data in config["resources"][cat][subcat].items():
                print(f"  Processing {item} ...")
                source_file = get_abs_dataset_path(
                    "datasets",
                    data["_source_path"].get("dataset"),
                    data_dir=True,
                ) / data["_source_path"].get("file")
                target_file = dataset_path / data.get("path")
                if target_file not in target_files:
                    source_files.append(source_file)
                    target_files.append(target_file)
                else:
                    print("    Target file already collected, skipping...")

    return source_files, target_files
