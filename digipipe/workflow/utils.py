"""
Helper functions for the workflow
"""

import os
from digipipe.store.utils import get_abs_store_root_path


def search_data_workflows():
    """Search for snakefiles (*.smk) in data store.

    Returns
    -------
    list
        Paths to Snakemake files
    """
    smk_files = []
    for root, dirs, files in os.walk(get_abs_store_root_path()):
        for file in files:
            if file.endswith(".smk") and ".TEMPLATE" not in str(root):
                smk_files.append(os.path.join(root, file))
    return smk_files
