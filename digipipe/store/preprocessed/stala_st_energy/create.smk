"""
Snakefile for this dataset

Note: To include the file in the main workflow, it must be added to the respective module.smk .
"""
import pandas as pd

from digipipe.store.utils import get_abs_dataset_path

DATASET_PATH = get_abs_dataset_path(
    "preprocessed", "stala_st_energy", data_dir=True
)

rule extract_power_demand_ind:
    """
    Extract industrial electricity demand from Excel file and save to CSV
    """
    input:
        get_abs_dataset_path("raw", "stala_st_energy") / "data" /
            "Stromverbrauch_nach_Kreisen_ab_Jahr_2003.xlsx"
    output: DATASET_PATH / "power_demand_industry_st_districts.csv"
    run:
        data = pd.read_excel(
            input[0],
            **config["electricity_demand_ind"]["excel_file"],
            engine="openpyxl",
        ).rename(columns={
            "Kreisfreie Stadt\nLandkreis\n\nLand": "name"
        }).set_index("name")
        data.columns = [int(_[1]) for _ in data.columns.str.split("\n")]

        # Select desired years
        print(
            f"Available years for industrial demand: "
            f"{data.columns.to_list()}"
        )
        years = config["electricity_demand_ind"]["years"]
        if len(years) > 0:
            data = data[years]
            print(f"  Selected: {years}")

        data.to_csv(output[0])
