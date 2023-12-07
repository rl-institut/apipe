"""
Snakefile for this dataset

Note: To include the file in the main workflow, it must be added to the respective module.smk .
"""
import pandas as pd
from apipe.store.utils import get_abs_dataset_path

DATASET_PATH = get_abs_dataset_path(
    "datasets", "demand_heat_region2", data_dir=True)

rule heat_demand_build_type:
    """
    Extract values of ags and heat demands for building types 'residential'/'non-residential'/'industry' and dump as seperate CSV files
    """
    input:
        get_abs_dataset_path(
            "preprocessed", "wfbb_heat_atlas_bb") / "data" / "wfbb_heat_atlas.csv"
    output:
        demand_res_heat_demand=DATASET_PATH / "demand_residential_heat_demand.csv",
        demand_nonres_heat_demand=DATASET_PATH / "demand_non_residential_heat_demand.csv",
        demand_ind_heat_demand=DATASET_PATH / "demand_industry_heat_demand.csv",
    run:
        # Load the data
        data = pd.read_csv(input[0], dtype=str)
        column_selection = ["build_type_residential_heat_demand", "build_type_non_residential_heat_demand", "build_type_industry_heat_demand"]
        # Extract and save as separate CSVs
        for column_name in column_selection:
            extracted_data = data[["ags", column_name]].rename(columns={column_name: "2022"})
            output_file = DATASET_PATH / f"demand_{column_name.split('_')[2]}_heat_demand.csv"
            extracted_data.to_csv(output[column_selection.index(column_name)], index=False)
