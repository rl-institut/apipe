import pandas as pd


def calc_heat_pump_cops(
    t_high: list,
    t_low: list,
    quality_grade: float,
    consider_icing: bool = False,
    temp_icing: float = None,
    factor_icing: float = None,
    spf: float = None,
) -> list:
    """Calculate temperature-dependent COP of heat pumps including efficiency
    gain over time.

    COP-Code was adapted from oemof-thermal:
    https://github.com/oemof/oemof-thermal/blob/features/cmpr_heatpumps_and_chillers/src/oemof/thermal/compression_heatpumps_and_chillers.py

    Efficiency corrections are based upon increase of seasonal performance
    factor (SPF) for scenario year since today SQ.
    """

    # Expand length of lists with temperatures and convert unit to Kelvin
    length = max([len(t_high), len(t_low)])
    if len(t_high) == 1:
        list_t_high_k = [t_high[0] + 273.15] * length
    elif len(t_high) == length:
        list_t_high_k = [t + 273.15 for t in t_high]
    if len(t_low) == 1:
        list_t_low_k = [t_low[0] + 273.15] * length
    elif len(t_low) == length:
        list_t_low_k = [t + 273.15 for t in t_low]

    # Calculate COPs
    if not consider_icing:
        cops = [
            quality_grade * t_h / (t_h - t_l)
            for t_h, t_l in zip(list_t_high_k, list_t_low_k)
        ]

    # Temperatures below 2 degC lead to icing at evaporator in
    # heat pumps working with ambient air as heat source.
    elif consider_icing:
        cops = []
        for t_h, t_l in zip(list_t_high_k, list_t_low_k):
            if t_l < temp_icing + 273.15:
                cops = cops + [factor_icing * quality_grade * t_h / (t_h - t_l)]
            if t_l >= temp_icing + 273.15:
                cops = cops + [quality_grade * t_h / (t_h - t_l)]

    # Efficiency gain for scenario year
    if spf is not None:
        cops = [_ * spf for _ in cops]

    return cops


def process() -> None:
    # Get heatpump params and temperature timeseries
    hp_params = snakemake.config["heatpumps"].get("params")
    temp = pd.read_csv(
        snakemake.input.temperature[0],
        index_col=0,
    )

    # Calculate COPs
    cops_ashp = calc_heat_pump_cops(
        t_high=[hp_params.get("heating_temp")],
        t_low=temp.temp_amb.to_list(),
        quality_grade=hp_params.get("quality_grade_ASHP"),
        consider_icing=True,
        temp_icing=hp_params.get("icing_temp"),
        factor_icing=hp_params.get("icing_factor"),
        spf=hp_params.get("spf_ASHP"),
    )
    cops_gshp = calc_heat_pump_cops(
        t_high=[hp_params.get("heating_temp")],
        t_low=temp.temp_soil.to_list(),
        quality_grade=hp_params.get("quality_grade_GSHP"),
        spf=hp_params.get("spf_GSHP"),
    )

    # Round and dump
    pd.Series(cops_ashp, name="cop_ashp",).round(
        3
    ).to_csv(snakemake.output.cop_ashp)
    pd.Series(cops_gshp, name="cop_gshp",).round(
        3
    ).to_csv(snakemake.output.cop_gshp)


if __name__ == "__main__":
    process()
