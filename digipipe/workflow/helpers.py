"""
Helper functions for the workflow
"""

import inspect
import os
from digipipe import store


def search_data_workflows():
    """Search for snakefiles (*.smk) in data store"""
    smk_files = []
    for root, dirs, files in os.walk(os.path.dirname(inspect.getfile(store))):
        for file in files:
            if file.endswith(".smk"):
                smk_files.append(os.path.join(root,file))
    return smk_files
