from pathlib import Path

import pandas as pd

from apipe.esys.esys.model import model_structures
from apipe.esys.esys.tools import update_foreign_keys_for_heatpumps as hp

# 1. Update foreign keys for heatpumps
hp.update_foreign_keys_hp()

# 2. emob and output parameters postprocessing
# Set Paths
BASE_PATH = Path(__file__).parents[3] / "store"
DATA_PATH = (
    BASE_PATH / "appdata" / "esys" / "2045_scenario" / "preprocessed" / "data"
)

DATASET_EMOB_PATH = BASE_PATH / "datasets" / "demand_emobility_region" / "data"
DEMAND_FILE_PATH = (
    DATA_PATH / "sequences" / "electricity-demand_mob_profile.csv"
)

# Load demand profile
demand_df = pd.read_csv(
    DEMAND_FILE_PATH, delimiter=";", encoding="utf-8", index_col=0
)

# Get regions from model structure
regions = model_structures["model_structure_full"]["regions"]


def update_demand_mob_profile():
    """
    Updates the electricity demand mobility profile by replacing
    time series data for each region found in the dataset path.

    - Reads CSV files for each region if available.
    - Updates the corresponding column in `demand_df`.
    - Saves the updated DataFrame back to the original CSV file.
    """

    for region in regions:
        region_file = (
            DATASET_EMOB_PATH
            / f"emobility_charging_demand_{region}_normalized.csv"
        )

        if not region_file.exists():
            print(
                f"⚠️ Warning: File not found for region {region}: "
                f"{region_file}"
            )
            continue  # Skip to the next region

        # Read region-specific time series
        region_df = pd.read_csv(
            region_file, delimiter=",", encoding="utf-8", index_col=0
        )

        # Ensure the region file is not empty before updating
        if region_df.empty:
            print(
                f"⚠️ Warning: File for region {region} is empty: "
                f"{region_file}"
            )
            continue

        # Update the demand profile DataFrame
        column_name = f"{region}-electricity-demand_mob-profile"
        demand_df[column_name] = region_df.iloc[:, 0].values
        print(f"✅ Successfully updated demand_mob-profile for {region}")

    # Save the updated demand profile
    demand_df.to_csv(DEMAND_FILE_PATH, encoding="utf-8", sep=";", index=True)
    print(f"📁 Demand profile updated and saved: {DEMAND_FILE_PATH}")


def update_output_parameters():
    """
    Updates the 'output_parameters' column in component CSV files
    for specified commodities.

    - Iterates over `output_parameters` dictionary.
    - Matches components from `model_structure_full`.
    - Reads the corresponding CSV file.
    - Updates the 'output_parameters' column for the given key.
    - Saves the modified CSV file.
    """

    output_parameters = {
        "r120640428428-h2-commodity": '{"full_load_time_max":6500.0}',
        "r120640428428-biomass_gas-commodity": '{"full_load_time_max":1}',
        "r120640428428-biomass_solid-commodity": '{"full_load_time_max":1}',
        "r120640472472-biomass_gas-commodity": '{"full_load_time_max":1}',
        "r120640472472-biomass_solid-commodity": '{"full_load_time_max":1}',
        "r120670124124-biomass_gas-commodity": '{"full_load_time_max":1}',
        "r120670124124-biomass_solid-commodity": '{"full_load_time_max":1}',
        "r120670201201-biomass_gas-commodity": '{"full_load_time_max":1}',
        "r120670201201-biomass_solid-commodity": '{"full_load_time_max":1}',
    }


    # Regions from model structure
    valid_regions = model_structures["model_structure_full"]["regions"]

    # Iterate through the output_parameters
    for key, value in output_parameters.items():
        # Extract region ID from the key (assuming it's the first part of the key before the first "-")
        region_id = key.split("-")[0]

        # Check if the region is in the valid region list
        if region_id not in valid_regions:
            print(f"Skipping {key}: Region {region_id} not in model structure")
            continue

        matching_components = [
            comp
            for comp in model_structures["model_structure_full"]["components"]
            if comp in key
        ]

        for component in matching_components:
            commodity_file = DATA_PATH / "elements" / f"{component}.csv"

            if not commodity_file.exists():
                print(f"File not found: {commodity_file}")
                continue  # Skip to the next component

            # Load CSV file
            commodity_df = pd.read_csv(
                commodity_file, delimiter=";", encoding="utf-8", index_col=0
            )

            # Update the DataFrame and save it back
            try:
                commodity_df.loc[key, "output_parameters"] = value
                print(f"Successfully updated output parameters for {key}")
                commodity_df.to_csv(
                    commodity_file, encoding="utf-8", sep=";", index=True
                )
            except Exception as e:
                print(f"Error updating {key}: {e}")


update_output_parameters()
update_demand_mob_profile()
