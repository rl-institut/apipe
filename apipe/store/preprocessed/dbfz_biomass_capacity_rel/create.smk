"""
Snakefile for this dataset

Note: To include the file in the main workflow, it must be added to the respective module.smk .
"""

import pandas as pd
from apipe.store.utils import get_abs_dataset_path


DATASET_PATH = get_abs_dataset_path("preprocessed", "dbfz_biomass_capacity_rel")

rule create:
    input:
        get_abs_dataset_path(
            "raw", "dbfz_biomass_heat_capacities"
        ) / "data" / "dbfz_biomass_heat_capacities.csv"
    output:
        output_cen=DATASET_PATH / "data" / "dbfz_biomass_capacity_rel_central.csv",
        output_dec=DATASET_PATH/ "data" / "dbfz_biomass_capacity_rel_decentral.csv"
    run:
        for supply in config["supplies"]:
            all_shares_biomass_cap = pd.DataFrame()
            for year in config["years"]:
                # Read raw file
                biomass_cap = pd.read_csv(
                    input[0],
                    sep=";",
                    usecols=["carrier", "tech", "capacity_[MW]_" + year, supply],
                )

                # Drop nan values
                biomass_cap = biomass_cap.dropna()

                # Group by the columns "carrier" and "tech"
                biomass_cap = biomass_cap.groupby(by=["carrier", "tech"]).sum(
                    numeric_only=True
                )

                # Calculate the share of the installed capacity per grouped
                # technologies
                share_biomass_cap = biomass_cap.apply(
                    lambda biomass_cap: biomass_cap / biomass_cap.sum()
                )
                # Rename the column and drop the year from the name of the column with
                # the capacity
                share_biomass_cap.rename(
                    {"capacity_[MW]_" + year: "capacity_rel"}, axis=1, inplace=True
                )

                # Drop multiindex made by grouping
                share_biomass_cap.reset_index(inplace=True)

                # Write tech into carrier col with the convention:
                # "carrier_tech"
                share_biomass_cap["carrier"] = share_biomass_cap[
                    ["carrier", "tech"]
                ].agg("_".join, axis=1)

                # Delete redundant column "tech"
                del share_biomass_cap["tech"]

                # Add year according to year mapping
                share_biomass_cap["year"] = config["year_mapping"][year]

                # Set index on the year to move the column to the front
                share_biomass_cap.set_index("year", inplace=True)
                share_biomass_cap.reset_index(inplace=True)

                # Update final df with all shares
                all_shares_biomass_cap = pd.concat(
                    [all_shares_biomass_cap, share_biomass_cap]
                )

            # Save the file
            if supply == "central":
                all_shares_biomass_cap.to_csv(output.output_cen, index=False)
            elif supply == "decentral":
                all_shares_biomass_cap.to_csv(output.output_dec, index=False)
            else:
                raise ValueError("Please provide either 'central' or 'decentral'"
                                 "with 'supply'")
