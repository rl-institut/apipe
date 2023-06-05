from pathlib import Path

import numpy as np
import pandas as pd


def extract_demand_ind(
    infile: Path, outfile_states: Path, outfile_districts: Path, cfg: dict
) -> None:
    data = pd.read_excel(
        infile,
        **cfg["demand_ind"]["excel_file"],
    )
    # Drop unnecessary headers
    data = data.iloc[3:]
    # Use correct zero value and set unavailable values to NaN
    data = data.replace("-", 0).replace(".", np.nan)

    # Convert from GJ to MWh
    data = data.set_index(["lau_code", "name"]).div(3.6).reset_index()

    # Filter: federal states
    data_states = data.loc[data.lau_code.str.len() == 2].set_index("lau_code")
    data_states["name"] = data_states["name"].str[2:]
    data_states.to_csv(outfile_states)

    # Filter: districts
    data_districts = data.loc[data.lau_code.str.len() == 5].set_index(
        "lau_code"
    )
    data_districts["name"] = data_districts["name"].str[6:]
    data_districts.to_csv(outfile_districts)


def extract_employment_ind(
    infile: Path,
    outfile_states: Path,
    outfile_districts: Path,
    outfile_muns: Path,
    cfg: dict,
) -> None:
    data = pd.read_excel(
        infile,
        **cfg["employment_ind"]["excel_file"],
    )

    # Drop unnecessary headers
    data = data.iloc[3:]
    # Use correct zero value and set unavailable values to NaN
    data = data.replace("-", 0).replace(".", np.nan)

    # Filter: federal states
    print("Extracting industrial employees and companies: Federal state level")
    data_states = data.loc[data.ags.str.len() == 2].set_index("ags")
    data_states["name"] = data_states["name"].str[2:]
    data_states.to_csv(outfile_states)
    print(f"  Employees: {data_states.employees_ind.sum()}")
    print(f"  Companies: {data_states.companies_ind.sum()}")

    # Filter: districts
    print("Extracting industrial employees and companies: District level")
    data_districts = data.loc[data.ags.str.len() == 5].set_index("ags")
    data_districts["name"] = data_districts["name"].str[6:]
    data_districts.to_csv(outfile_districts)
    print(f"  Employees: {data_districts.employees_ind.sum()}")
    print(f"  Companies: {data_districts.companies_ind.sum()}")
    if (
        data_districts.employees_ind.sum() != data_states.employees_ind.sum()
    ) or (
        data_districts.companies_ind.sum() != data_states.companies_ind.sum()
    ):
        print(
            "  WARNING: Numbers do not equal state values, "
            "probably due to missing values!"
        )

    # Filter: municipalities
    # Note: Kreisfreie Städte are not included in muns, so they're inserted
    #       manually
    print("Extracting industrial employees and companies: Municipal level")
    data_muns = data.loc[data.ags.str.len() == 8].set_index("ags")
    data_muns["name"] = data_muns["name"].str[8:]

    missing_districts = data_districts.loc[
        ~data_districts.index.isin(data_muns.index.str[:5])
    ]
    missing_districts.index += "000"
    data_muns = pd.concat([data_muns, missing_districts], axis=0).sort_index()
    print(
        f"  {len(missing_districts)} "
        f"municipalities missing, filled up with districts"
    )
    print(f"  Employees: {data_muns.employees_ind.sum()}")
    print(f"  Companies: {data_muns.companies_ind.sum()}")
    if (data_muns.employees_ind.sum() != data_states.employees_ind.sum()) or (
        data_muns.companies_ind.sum() != data_states.companies_ind.sum()
    ):
        print(
            "  WARNING: Numbers do not equal state values, "
            "probably due to missing values!"
        )

    data_muns.to_csv(outfile_muns)
