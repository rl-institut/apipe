"""
Helper functions for datasets
"""

import inspect
import os
from digipipe import store


def get_abs_store_path():
    return os.path.dirname(inspect.getfile(store))
