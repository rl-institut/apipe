import geopandas as gpd
import pandas as pd


def process() -> None:
    # pylint: disable=R0914
    # Get muns and their NUTS3 code
    muns = gpd.read_file(snakemake.input.region_muns)
    districts = gpd.read_file(snakemake.input.region_districts)
    federal_states = gpd.read_file(snakemake.input.federal_states).set_index(
        "nuts"
    )
    muns = muns.merge(
        districts[["id", "nuts"]].rename(columns={"id": "district_id"}),
        left_on="district_id",
        right_on="district_id",
    )

    # Historical data: load, get available years and filter by region
    pop_history = pd.concat(
        [
            pd.read_csv(f, dtype={"ags": str}, index_col="ags")
            for f in snakemake.input.pop_history
        ],
        axis=1,
    )
    pop_history = pop_history.loc[pop_history.index.isin(muns.ags)]
    avail_years_history = pop_history.columns.astype(int).to_list()
    print(f"Historical years: {avail_years_history}")

    # Prognosis based on federal state numbers (linear scaling)
    pop_prognosis_years_fstates = snakemake.config[
        "prognosis_country_fstatelevel"
    ]["years"]
    if len(pop_prognosis_years_fstates) > 0:
        pop_prognosis_fstate = pd.read_csv(
            snakemake.input.prognosis_germany_fstatelevel,
            index_col="federal_state",
        )
        avail_years_prognosis = pop_prognosis_fstate.columns.astype(
            str
        ).to_list()
        print(
            f"Prognosis years (country data on federal state level) requested: "
            f"{pop_prognosis_years_fstates}"
        )
        print(
            f"Prognosis years (country data on federal state level) available: "
            f"{avail_years_prognosis}"
        )

        # Get federal state data
        fstate_nuts = muns.nuts.str[:3].unique()
        if len(fstate_nuts) != 1:
            raise ValueError(
                "Region is located in more than one federal state!"
            )
        fstate_name = federal_states.loc[fstate_nuts].name[0]
        pop_prognosis_fstate = pop_prognosis_fstate.loc[fstate_name]

        # Extrapolate
        latest_hist_year = pop_history.columns.astype(int).max().astype(str)
        if latest_hist_year not in pop_prognosis_fstate.index:
            raise ValueError(
                "Latest historical year not found in prognosis data, "
                "cannot scale!"
            )
        if not all(
            _ > int(latest_hist_year) for _ in pop_prognosis_years_fstates
        ):
            raise ValueError(
                "At least one requested prognosis year is not in the future, "
                "cannot scale!"
            )
        pop_prognosis_fstate = pop_prognosis_fstate[
            [str(_) for _ in [latest_hist_year] + pop_prognosis_years_fstates]
        ]
        pop_prognosis_factors = pop_prognosis_fstate.div(
            pop_prognosis_fstate.loc[latest_hist_year]
        ).drop(latest_hist_year)
        pop_prognosis_mun = pd.DataFrame(
            {
                y: (f * pop_history[latest_hist_year]).round()
                for y, f in pop_prognosis_factors.items()
            }
        )
        population = pd.concat([pop_history, pop_prognosis_mun], axis=1)
    else:
        population = pop_history

    # Merge historical and prognosis data

    # Add municipality_id and data origin
    population = (
        pd.concat(
            [muns.set_index("ags")["id"].rename("municipality_id"), population],
            axis=1,
        )
        .sort_index()
        .set_index("municipality_id", drop=True)
    )
    population.columns = pd.MultiIndex.from_arrays(
        [
            population.columns,
            len(avail_years_history) * ["historic"]
            + len(pop_prognosis_mun.columns) * ["prognosis"],
        ],
        names=("year", "type"),
    )

    population.to_csv(snakemake.output[0])


if __name__ == "__main__":
    process()
