import geopandas as gpd
import pandas as pd


class PanelSettings:
    """Value store for settings panel"""

    def __init__(self, **settings):
        """Create attributes from yaml file"""
        self.__dict__.update(settings)

    def update(self, **kwargs):
        """Updates control element's attributes"""

        def _keys_available(keys, keys_required):
            return sorted(keys) == sorted(keys_required)

        for control, values in kwargs.items():
            if hasattr(self, control):
                values_store = getattr(self, control)
                if not isinstance(values, dict):
                    raise ValueError("Attributes are no dict.")
                if _keys_available(values.keys(), values_store.keys()):
                    values_store.update(values)
                    setattr(self, control, values_store)
                else:
                    raise ValueError(
                        "Attributes in given value dict do not match the "
                        "value store schema! Check the config.yml!"
                    )
            else:
                raise ValueError(
                    f"{control} is no control element! Check the config.yml"
                )

    def is_complete(self):
        """Returns True if all values are set (no default value "none" left)"""
        for control, values_store in self.__dict__.items():
            for val in values_store.values():
                if val == "none":
                    print(f"Control {control} is missing at least one value.")
                    return False
        return True

    def make_dict(self):
        """Make dictionary from value store"""
        raise NotImplementedError


def generate_energy_panel_data(
    panel_settings: PanelSettings,
    region: gpd.GeoDataFrame,
    tech_data: dict,
    wind_stats: pd.DataFrame,
    wind_area_stats: pd.DataFrame,
    pv_ground_stats: pd.DataFrame,
    pv_ground_area_stats: pd.DataFrame,
    pv_ground_area_shares: dict,
    pv_roof_stats: pd.DataFrame,
    pv_roof_area_stats: pd.DataFrame,
    pv_roof_area_deploy_stats: pd.DataFrame,
    pv_ground_targets: dict,
    pv_roof_targets: dict,
    hydro_stats: pd.DataFrame,
    demand_hh_power: pd.DataFrame,
    demand_cts_power: pd.DataFrame,
    demand_ind_power: pd.DataFrame,
    storage_large_stats: pd.DataFrame,
    # storage_small_stats: pd.DataFrame,
    storage_pv_roof: dict,
):
    # Wind energy
    panel_settings.update(
        **dict(
            s_w_1=dict(
                max=round(
                    wind_area_stats.stp_2018_vreg.sum()
                    * tech_data["power_density"]["wind"]
                ),
                min=0,
                start=round(wind_stats.capacity_net.sum()),
                step=10,
                status_quo=round(wind_stats.capacity_net.sum()),
                future_scenario=round(
                    float(region.area_km2)
                    * 0.02
                    * tech_data["power_density"]["wind"]
                ),
            ),
            s_w_3=dict(start=True),
            s_w_4=dict(start=False),
            s_w_4_1=dict(start=False),
            s_w_4_2=dict(start=False),
            s_w_5=dict(start=False),
        )
    )

    # PV ground and roof
    panel_settings.update(
        **dict(
            s_pv_ff_1=dict(
                max=round(
                    pv_ground_area_stats.sum().sum()
                    * tech_data["power_density"]["pv_ground"]
                ),
                min=0,
                start=round(pv_ground_stats.capacity_net.sum()),
                step=10,
                status_quo=round(pv_ground_stats.capacity_net.sum()),
                future_scenario=round(pv_ground_targets["target_power_total"]),
            ),
            s_pv_ff_3=dict(
                max=pv_ground_area_shares["road_railway"] * 100,
                min=0,
                start=0,
                step=0.25,
            ),
            s_pv_ff_4=dict(
                max=pv_ground_area_shares["agri"] * 100,
                min=0,
                start=0,
                step=5,
            ),
            s_pv_d_1=dict(
                max=round(
                    # max. 30 % of all roofs
                    pv_roof_area_stats.installable_power_total.sum()
                    * 0.3
                ),
                min=0,
                start=round(pv_roof_stats.capacity_net.sum()),
                step=10,
                status_quo=round(pv_roof_stats.capacity_net.sum()),
                future_scenario=round(pv_roof_targets["target_power_total"]),
            ),
            s_pv_d_3=dict(
                max=50,
                min=0,
                start=round(
                    pv_roof_area_deploy_stats.capacity_net.sum()
                    / pv_roof_area_deploy_stats.installable_power.sum()
                    * 100
                ),
                step=5,
                status_quo=round(
                    pv_roof_area_deploy_stats.capacity_net.sum()
                    / pv_roof_area_deploy_stats.installable_power.sum()
                    * 100
                ),
            ),
            s_pv_d_4=dict(
                max=100,
                min=0,
                start=round(
                    storage_pv_roof["pv_roof_share"]["home_storages"] * 100
                ),
                step=5,
                status_quo=round(
                    storage_pv_roof["pv_roof_share"]["home_storages"] * 100
                ),
            ),
        )
    )

    # Hydro
    panel_settings.update(
        **dict(
            s_h_1=dict(
                max=round(hydro_stats.capacity_net.sum()),
                min=0,
                start=round(hydro_stats.capacity_net.sum()),
                step=1,
                status_quo=round(hydro_stats.capacity_net.sum()),
                disable=True,
            ),
        )
    )

    # Demand
    total_demand = (demand_hh_power + demand_cts_power + demand_ind_power).sum()
    feedin_wind_pv_daily_mean = (
        (
            pv_roof_stats.capacity_net.sum()
            * tech_data["full_load_hours"]["pv_roof"]["2022"]
        )
        + (
            pv_ground_stats.capacity_net.sum()
            * tech_data["full_load_hours"]["pv_ground"]["2022"]
        )
        + (
            wind_stats.capacity_net.sum()
            * tech_data["full_load_hours"]["wind"]["2022"]
        )
    ) / 365

    panel_settings.update(
        **dict(
            s_v_1=dict(
                max=100,
                min=50,
                start=100,
                step=10,
                status_quo=100,
                future_scenario=round(
                    total_demand["2045"] / total_demand["2022"] * 100
                ),
            ),
            s_v_3=dict(
                max=200,
                min=50,
                start=100,
                step=10,
                status_quo=100,
                future_scenario=round(
                    demand_hh_power.sum()["2045"]
                    / demand_hh_power.sum()["2022"]
                    * 100
                ),
            ),
            s_v_4=dict(
                max=200,
                min=50,
                start=100,
                step=10,
                status_quo=100,
                future_scenario=round(
                    demand_cts_power.sum()["2045"]
                    / demand_cts_power.sum()["2022"]
                    * 100
                ),
            ),
            s_v_5=dict(
                max=200,
                min=50,
                start=100,
                step=10,
                status_quo=100,
                future_scenario=round(
                    demand_ind_power.sum()["2045"]
                    / demand_ind_power.sum()["2022"]
                    * 100
                ),
            ),
            s_s_g_1=dict(
                max=round(feedin_wind_pv_daily_mean / 10),
                min=0,
                start=round(storage_large_stats.storage_capacity.sum()),
                step=0.1,
                status_quo=round(storage_large_stats.storage_capacity.sum()),
            ),
            s_s_g_3=dict(
                max=10,
                min=0,
                start=round(
                    storage_large_stats.storage_capacity.sum()
                    / feedin_wind_pv_daily_mean
                    * 100,
                    1,
                ),
                step=0.25,
                status_quo=round(
                    storage_large_stats.storage_capacity.sum()
                    / feedin_wind_pv_daily_mean
                    * 100,
                    1,
                ),
            ),
        )
    )

    return panel_settings
