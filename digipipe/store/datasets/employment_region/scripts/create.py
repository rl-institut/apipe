import geopandas as gpd
import pandas as pd


def process() -> None:
    # Get muns and their NUTS3 code
    muns = gpd.read_file(snakemake.input.region_muns)

    employment = pd.read_csv(
        snakemake.input.employment,
        index_col=0,
        dtype={"ags": str}
    ).loc[muns.ags.to_list()][["employees"]]

    employment = muns.rename(
        columns={"id": "municipality_id"}
    ).set_index("ags").merge(
        employment,
        left_index=True,
        right_index=True
    ).set_index("municipality_id")[["employees"]]

    employment.to_csv(snakemake.output[0])


if __name__ == "__main__":
    process()
