import geopandas as gpd
import pandas as pd


def process() -> None:
    # Get muns
    muns = gpd.read_file(snakemake.input.region_muns)

    # Get total employment data
    employment_total = pd.read_csv(
        snakemake.input.employment_total, index_col=0, dtype={"ags": str}
    ).loc[muns.ags.to_list()]

    # Get industrial employment data
    employment_ind = pd.read_csv(
        snakemake.input.employment_ind, index_col=0, dtype={"ags": str}
    ).loc[muns.ags.to_list()][["employees_ind", "companies_ind"]]

    # Fill missing data by average employees per company
    employees_per_company = employment_ind.employees_ind.div(
        employment_ind.companies_ind
    )

    if employees_per_company.isna().sum() > 0:
        print(
            f"WARNING: Number of employees or companies missing in "
            f"{employees_per_company.isna().sum()} of {len(muns)} "
            f"municipalities! Using average for employees per company "
            f"({employees_per_company.mean().round(1)}) for those..."
        )
        employment_ind = employment_ind.assign(
            employees_ind=employment_ind.employees_ind.fillna(
                employment_ind.companies_ind * employees_per_company.mean()
            )
        )
    employment_ind = employment_ind.assign(
        employees_ind=employment_ind.employees_ind.round().astype(int)
    )

    # Join total data with industry data
    employment = employment_total.join(employment_ind)

    # Join employment data with municipality ids
    employment = (
        muns.rename(columns={"id": "municipality_id"})
        .set_index("ags")
        .merge(employment, left_index=True, right_index=True)
        .set_index("municipality_id")
    )[["employees_total", "companies_total", "employees_ind", "companies_ind"]]

    employment.to_csv(snakemake.output[0])


if __name__ == "__main__":
    process()
