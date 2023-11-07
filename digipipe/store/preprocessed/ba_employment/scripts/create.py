import pandas as pd


def process() -> None:
    # Get Excel config
    excel_cfg = snakemake.config["excel_file"]

    # Read file
    data = pd.read_excel(
        snakemake.input[0],
        **excel_cfg,
        engine="pyxlsb",
    )

    # Set empty values to 0 and convert dtypes
    data = data.replace(["*", " "], 0)
    data = data.astype({"employees_total": int, "companies_total": int})

    # Drop empoty and aggregated rows
    data = data.loc[data.ags.str.len() == 8].set_index("ags")

    # Write
    data.to_csv(snakemake.output[0])


if __name__ == "__main__":
    process()
