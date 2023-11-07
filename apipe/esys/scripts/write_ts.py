# coding: utf-8
r"""
This module provides functionality to map time series from the subdirectories
in "store/datasets/" to the empty time series in
"store/datasets/esys_raw/data/time_series/empty_ts_efficiencies.csv",
"store/datasets/esys_raw/data/time_series/empty_ts_feedin.csv" and
"store/datasets/esys_raw/data/time_series/empty_ts_load.csv".
The empty files are initially created using the "create_empty_ts" rule.

The mapped time series are then written to
"store/datasets/esys_raw/data/scalars/ts_efficiencies.csv",
"store/datasets/esys_raw/data/scalars/ts_feedin.csv" and
"store/datasets/esys_raw/data/scalars/ts_load.csv".

The original files with the empty time series are deleted.
"""

import os
import sys

import pandas as pd

from apipe.esys.esys.config.esys_conf import add_snake_logger, map_ts
from apipe.esys.esys.tools.data_processing import (
    HEADER_B3_TS,
    load_b3_timeseries,
    save_df,
    stack_timeseries,
    unstack_timeseries,
)
from apipe.store.utils import get_abs_dataset_path


def check_if_ts_name_in_maps_ts(ts_name, mapping_file_name):
    """
    Check if a time series name is present in the mapping of time series.

    Parameters
    ----------
    ts_name : str
        The name of the time series to check.
    mapping_file_name : str
        The name of the mapping file for reference in the error message.

    Raises
    ------
    ValueError
        If the time series name is missing in the mapping.

    """
    if (
        ts_name not in map_ts.efficiencies
        and ts_name not in map_ts.feedin
        and ts_name not in map_ts.load
    ):
        raise ValueError(
            f"Missing mapping of time series '{ts_name}'. "
            f"Please provide a key and value for each time "
            f"series with the file '{mapping_file_name}'."
        )


def get_datasets_data_name(key, which):
    """
    Retrieve the specific data value from the time series based on the provided
    key and description.

    Parameters
    ----------
    key : str
        The name of the data to retrieve.
    which : str
        Specifies the kind of time series data.
        Valid options are "efficiencies", "feedin", and "load".

    Returns
    -------
    value : object
        The data value corresponding to the given key and description.

    Raises
    ------
    ValueError
        If an invalid option is provided for the 'which' parameter.

    """
    if which == "efficiencies":
        value = map_ts.efficiencies[key]
    elif which == "feedin":
        value = map_ts.feedin[key]
    elif which == "load":
        value = map_ts.load[key]
    else:
        raise ValueError(
            "Please provide a valid description of the time series with "
            f"'which'. '{which}' is not a valid option."
        )

    return value


def check_file_exists(path, file_name, df, mapping_file_name):
    """
    Check if a file exists at the specified path and raise a FileNotFoundError
    if not.

    Parameters
    ----------
    path : str
        The path to the file.
    file_name : str
        The name of the file (without the file extension).
    df : pd.DataFrame
        The DataFrame containing the variable name to be mapped.
    mapping_file_name : str
        The name of the mapping file for reference in the error message.

    Raises
    ------
    FileNotFoundError
        If the file does not exist at the specified path.

    """
    if not path:
        raise FileNotFoundError(
            f"The file '{file_name}.csv' could not be found. Please provide a "
            f"valid file name with '{mapping_file_name}' to map "
            f"'{df['var_name']}'."
        )


def write_ts_data(stacked_ts, data):
    """
    Write time series data to the specified stacked time series DataFrame.

    Parameters
    ----------
    stacked_ts : pd.DataFrame
        The stacked time series DataFrame.
    data : pd.DataFrame
        The data to write into the stacked time series.

    Returns
    -------
    pd.DataFrame
        The updated stacked time series DataFrame.

    """
    unstacked_data = unstack_timeseries(stacked_ts)
    unstacked_data.iloc[:, 0] = data.iloc[:, 0].values

    stacked_data = stack_timeseries(unstacked_data)
    stacked_ts.loc[:, "series"] = stacked_data["series"].values
    return stacked_ts


def map_over_var_name_ts(path_ts, which, path_ts_updated):
    """
    Map and update time series data in directory 'path_ts' based on var_name.
    Save the updated time series to the directory 'path_ts_updated'.
    Delete the origin time series (empty_ts_...) in 'path_ts' directory.

    Parameters
    ----------
    path_ts : str
        The path to the original time series file.
    which : str
        Specifies the type of data value to map and update.
    path_ts_updated : str
        The path to save the updated time series file.

    Returns
    -------
    None

    """
    # Create empty Dataframe for the updated time series
    updated_ts_df = pd.DataFrame(columns=HEADER_B3_TS)

    # Load set of time series
    ts_set = load_b3_timeseries(path_ts)

    # Loop over rows (set) of ts
    for index, row in ts_set.iterrows():
        check_if_ts_name_in_maps_ts(row["var_name"], "map_ts.yml")
        datasets_file_name = get_datasets_data_name(row["var_name"], which)
        datasets_file_path = get_abs_dataset_path(
            "datasets", datasets_file_name.get("dataset"), data_dir=True
        ) / datasets_file_name.get("file")

        check_file_exists(
            datasets_file_path, datasets_file_name, row, "map_ts.yml"
        )

        df_new_data_ts = pd.read_csv(datasets_file_path, usecols=[1], sep=",")

        ts_updated = write_ts_data(
            ts_set.iloc[[index], :].copy(), df_new_data_ts
        )

        updated_ts_df = pd.concat(
            [updated_ts_df, ts_updated], ignore_index=False
        )

        updated_ts_df.index.name = ts_updated.index.name

    save_df(updated_ts_df, path_ts_updated)

    os.remove(path_ts)
    # Print user info
    logger.info(f"The file {path_ts} has been deleted.")


if __name__ == "__main__":
    path_empty_ts_efficiencies = sys.argv[1]
    path_empty_ts_feedin = sys.argv[2]
    path_empty_ts_load = sys.argv[3]

    path_ts_efficiencies = sys.argv[4]
    path_ts_feedin = sys.argv[5]
    path_ts_load = sys.argv[6]

    logger = add_snake_logger("data_processing")

    map_over_var_name_ts(
        path_empty_ts_efficiencies, "efficiencies", path_ts_efficiencies
    )
    map_over_var_name_ts(path_empty_ts_feedin, "feedin", path_ts_feedin)
    map_over_var_name_ts(path_empty_ts_load, "load", path_ts_load)
