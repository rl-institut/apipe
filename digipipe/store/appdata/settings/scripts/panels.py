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
):
    wind = dict(
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
    panel_settings.update(**wind)

    pv = dict(
        s_pv_ff_1=dict(
            max=round(
                pv_ground_area_stats.sum().sum()
                * tech_data["power_density"]["pv_ground"]
            ),
            min=0,
            start=round(pv_ground_stats.capacity_net.sum()),
            step=10,
            status_quo=round(pv_ground_stats.capacity_net.sum()),
            future_scenario=round(
                float(region.area_km2)
                * 0.007
                * tech_data["power_density"]["pv_ground"]
            ),
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
                pv_roof_area_stats.installable_power_total.sum() * 0.3  # TODO
            ),
            min=0,
            start=round(pv_roof_stats.capacity_net.sum()),
            step=10,
            status_quo=round(pv_roof_stats.capacity_net.sum()),
            future_scenario="none",
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
            start="none",  # TODO
            step=5,
        ),
    )
    panel_settings.update(**pv)

    # import pdb
    # pdb.set_trace()

    return panel_settings
