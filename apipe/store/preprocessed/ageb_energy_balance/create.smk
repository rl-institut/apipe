"""
Snakefile for this dataset

Note: To include the file in the main workflow, it must be added to the respective module.smk .
"""
import tabula
from apipe.store.utils import get_abs_dataset_path

DATASET_PATH = get_abs_dataset_path("preprocessed", "ageb_energy_balance", data_dir=True)

rule convert:
    """
    Convert PDF to CSVs
    """
    input:
        get_abs_dataset_path("raw", "ageb_energy_balance")
        / "data" / "AGEB_22p2_rev-1.pdf"
    output:
        DATASET_PATH / "ageb_energy_balance_germany_{sector}_twh_2021.csv",
    run:
        # Read data from PDF
        demand = tabula.read_pdf(
            input[0],
            pages=[config["sectors"][wildcards.sector]["page"]],
            area=config["area"],
            relative_area=True
        )[0]

        # Processing
        demand.drop(columns=["Unnamed: 0"], inplace=True)
        demand.columns=config["sectors"][wildcards.sector]["column_names"]
        demand = demand.set_index("carrier").replace({"-": 0})
        for col in demand.columns:
            demand[col] = demand[col].str.replace(
                ',', '.').str.replace(" ", "").astype(float).fillna(0)
        demand.drop("Insgesamt", axis=0, inplace=True)

        # PJ to TWh
        demand = demand.div(3.6)

        # Dump as CSV
        demand.to_csv(output[0])
