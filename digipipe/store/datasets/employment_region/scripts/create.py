import geopandas as gpd
import pandas as pd


def process() -> None:
    # Get muns
    muns = gpd.read_file(snakemake.input.region_muns)

    # Get employment data
    employment = pd.read_csv(
        snakemake.input.employment, index_col=0, dtype={"ags": str}
    ).loc[muns.ags.to_list()][["employees_total"]]

    # Join employment data with municipality ids
    employment = (
        muns.rename(columns={"id": "municipality_id"})
        .set_index("ags")
        .merge(employment, left_index=True, right_index=True)
        .set_index("municipality_id")[["employees_total"]]
    )

    employment.to_csv(snakemake.output[0])


if __name__ == "__main__":
    process()
