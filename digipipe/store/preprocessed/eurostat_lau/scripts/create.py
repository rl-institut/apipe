import pandas as pd


def process() -> None:
    cfg = snakemake.config
    laus = pd.read_excel(
        snakemake.input[0],
        **cfg["excel_file"],
        engine="openpyxl",
    )

    # Rename columns
    laus.columns = cfg["column_names"]
    laus.set_index("lau_code", inplace=True)
    laus.dropna(inplace=True)

    # Dump
    laus.to_csv(snakemake.output[0])


if __name__ == "__main__":
    process()
