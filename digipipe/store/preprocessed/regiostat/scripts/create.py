import numpy as np
import pandas as pd


def process() -> None:
    cfg = snakemake.config
    demand = pd.read_excel(
        snakemake.input[0],
        **cfg["excel_file"],
        engine="openpyxl",
    )
    # Rename columns and drop unnecessary headers
    demand.columns = cfg["column_names"]
    demand = demand.iloc[3:]
    # Use correct zero value and set unavailable values to NaN
    demand = demand.replace("-", 0).replace(".", np.nan)

    # Convert from GJ to MWh
    demand = demand.set_index(["lau_code", "name"]).div(3.6).reset_index()

    # Filter: federal states
    demand_states = demand.loc[demand.lau_code.str.len() == 2].set_index(
        "lau_code"
    )
    demand_states.name = demand_states.name.str[2:]
    demand_states.to_csv(snakemake.output.demand_states)

    # Filter: districts
    demand_districts = demand.loc[demand.lau_code.str.len() == 5].set_index(
        "lau_code"
    )
    demand_districts.name = demand_districts.name.str[6:]
    demand_districts.to_csv(snakemake.output.demand_districts)


if __name__ == "__main__":
    process()
