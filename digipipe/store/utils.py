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


def create_tag_string_osmium(taglist):
    """Create tag string required by osmium for extraction of OSM data.

    Parameters
    ----------
    taglist : list of lists
        List of tags, format: [[key_1, value_1], ... [key_n, value_n]]

    Returns
    -------
    str
        String with tags
    """
    tag_string = ""
    for k, v in taglist:
        tag_string += "=".join([k, str(v)]) + " "
    tag_string = tag_string[:-1]
    return tag_string


def create_tag_string_ogr(taglist):
    """Create tag string required by ogr2ogr for extraction of OSM data.

    Parameters
    ----------
    taglist : list of lists
        List of tags, format: [[key_1, value_1], ... [key_n, value_n]]

    Returns
    -------
    dict
        Tags (key: "tags") and filter conditions (key: "conditions")
    """
    tag_conditions = "-where \"" + " OR ".join(
        ["=".join(["\\\"" + tag + "\\\"", f"\\\"" + str(val) + "\\\""])
         for tag, val in taglist]) + "\""
    tags = ",".join([tag for tag, _ in taglist])
    return {"tags": tags, "conditions": tag_conditions}
