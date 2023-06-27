# coding: utf-8
r"""
This module contains functionality to update empty scalars (generated from
scenarios with create_empty_scalars module) with default data and saves it to a
new csv file.
"""

import sys

import numpy as np

from digipipe.esys.esys.config.esys_conf import write_default_scalars
from digipipe.esys.esys.tools.data_processing import (
    filter_df,
    load_b3_scalars,
    save_df,
)


def get_var_value_and_comment(which):
    """
    Returns the variable value and a comment based on the input argument.

    Parameters
    ----------
    which : str
        string indicating what var_value and comment are chosen as default.
        Valid options are: zeros, empty, high_costs, false and empty_dict.

    Returns
    -------
    var_value : object
        The updated variable value
    comment : str
        A comment for the var_value update
    """
    if which == "zeros":
        var_value = 0
        comment = "Zero"
    elif which == "empty":
        var_value = np.nan
        comment = "Empty"
    elif which == "high_costs":
        var_value = 1000000000
        comment = "High slack cost on shortage"
    elif which == "false":
        var_value = False
        comment = "Empty"
    elif which == "empty_dict":
        var_value = "{}"
        comment = "Empty"
    elif which == "variable_costs":
        var_value = str({"variable_costs": 1e-7})
        comment = "Own assumption to prevent hidden curtailment"
    else:
        raise ValueError(
            f"'{which}' is not a valid option. Please provide a valid options. "
            f"Valid options are: zeros, empty, high_costs, false and "
            f"empty_dict."
        )

    return var_value, comment


def get_filter_condition(_df, _var_name, _type, _tech):
    """
    Returns a filter condition based on the input arguments.

    Parameters
    ----------
    _df : pd.DataFrame
        DataFrame containing the empty scalar data
    _var_name : str
        Name of the variable to be updated
    _type : str
        Type to be updated
    _tech : str
        Technology to be updated

    Returns
    -------
    condition : pd.Series
        A Boolean series used to filter the DataFrame to update the specified
        variable
    """
    precondition = _df["var_value"].isna()
    if (_type != "None") and (_tech != "None"):
        condition = (
            precondition
            & (_df["var_name"] == _var_name)
            & (_df["type"] == _type)
            & (_df["tech"] == _tech)
        )
    elif (_type != "None") and (_tech == "None"):
        condition = (
            precondition
            & (_df["var_name"] == _var_name)
            & (_df["type"] == _type)
        )
    elif (_type == "None") and (_tech != "None"):
        condition = (
            precondition
            & (_df["var_name"] == _var_name)
            & (_df["tech"] == _tech)
        )
    else:
        condition = precondition & (_df["var_name"] == _var_name)

    return condition


def update_df(_df, which, condition, unit):
    """
    Update the input DataFrame by setting values for specified columns based on
    the provided condition.

    Parameters
    ----------
    _df : pandas.DataFrame
        The DataFrame to update.
    which : str
        Specifies the type of value to set for the "var_value" column.
        Valid options are "zeros", "empty", "high_costs", "false", and
        "empty_dict".
    condition : pandas.Series of bools
        A boolean mask indicating which rows of `_df` should be updated.
    unit : str
        The unit to write into the "var_unit" column.

    Returns
    -------
    pandas.DataFrame
        The updated DataFrame.
    """
    var_value, comment = get_var_value_and_comment(which)

    _df.loc[condition, "var_value"] = var_value
    _df.loc[condition, "comment"] = comment
    _df.loc[condition, "var_unit"] = unit

    return _df


if __name__ == "__main__":
    path_empty_sc = sys.argv[1]
    path_default_sc = sys.argv[2]
    path_default_costs_eff = sys.argv[3]

    empty_sc_df = load_b3_scalars(path_empty_sc)

    write_empty_scalars_dict = write_default_scalars.write_default_scalars

    for key, value in write_empty_scalars_dict.items():

        condition = get_filter_condition(
            empty_sc_df, value["var_name"], value["type"], value["tech"]
        )
        df_updated = update_df(
            empty_sc_df, value["which"], condition, value["var_unit"]
        )

    # Get all unique values for var_name
    var_names_all = list(df_updated["var_name"].unique())
    # Get all values for var_name in context of costs and efficiencies
    var_names_costs_efficiencies = list(
        write_default_scalars.costs_efficiencies
    )
    # Get all remaining var_values needed for default_scalars.csv
    var_names_scalars = [
        var_name
        for var_name in var_names_all
        if var_name not in var_names_costs_efficiencies
    ]

    df_costs_efficiencies = filter_df(
        df_updated, "var_name", var_names_costs_efficiencies
    )

    df_scalars = filter_df(df_updated, "var_name", var_names_scalars)

    # Write all attributes attached to costs and efficiencies in separate
    # default_cost_efficiencies.csv file
    save_df(df_costs_efficiencies, path_default_costs_eff)

    # Write all other scalars in default_scalars.csv
    save_df(df_scalars, path_default_sc)
