"""
Helper functions for configs
"""

import yaml


def read_config(file: str) -> dict:
    """Reads a yml config and returns as dict

    Parameters
    ----------
    file : str
        Full path to config file to read

    Returns
    -------
    dict
        Config dict
    """
    with open(file, 'r') as cfg_file:
        try:
            cfg = yaml.safe_load(cfg_file)
        except yaml.YAMLError as exc:
            print(exc)
    return cfg
