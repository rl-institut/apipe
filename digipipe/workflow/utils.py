"""
Helper functions for the workflow
"""

import os
from digipipe.store.utils import get_abs_store_path


def search_data_workflows():
    """Search for snakefiles (*.smk) in data store"""
    smk_files = []
    for root, dirs, files in os.walk(get_abs_store_path()):
        for file in files:
            if file.endswith(".smk"):
                smk_files.append(os.path.join(root,file))
    return smk_files
