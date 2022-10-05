"""
Helper functions for datasets
"""

import inspect
import os
from pathlib import Path
from digipipe import store


def get_abs_store_root_path():
    """Get absolute path to data store

    Returns
    -------
    PosixPath
        Path to data store
    """
    return Path(os.path.dirname(inspect.getfile(store)))


def get_abs_dataset_path(category, name, data_dir=True):
    """Get absolute path to a dataset

    Parameters
    ----------
    category : str
        Category in data store, one of
        ["0_raw", "1-preprocessed", "2_dataset"]
    name : str
        Name of dataset (subdir)
    data_dir : bool
        If True, the subdir "data" where data lives is added to the path

    Returns
    -------
    PosixPath
        Path to dataset
    """
    p = Path(get_abs_store_root_path()) / category / name
    if data_dir is True:
        p = p / "data"
    return p
