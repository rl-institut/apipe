"""
Snakefile for this dataset

Note: To include the file in the main workflow, it must be added to the respective module.smk .
"""
import json
import geopandas as gpd
import pandas as pd

from digipipe.scripts.data_io import load_json
from digipipe.store.utils import get_abs_dataset_path
from digipipe.store.appdata.settings.scripts.panels import (
    PanelSettings,
    generate_energy_panel_data,
)

DATASET_PATH = get_abs_dataset_path("appdata", "settings", data_dir=True)


rule create_map_panel_layer_list:
    """
    Create layer list for right map panel
    """
    input: rules.appdata_datapackage_create_datapackage.output
    output: DATASET_PATH / "map_panel_layer_list.json"
    run:
        print("Creating list of layers...")
        with open(output[0], "w", encoding="utf8") as f:
            json.dump(config["map_panel_layer_list"], f, indent=4)

rule create_panel_settings:
    """
    Create startup settings for settings panels
    """
    input:
        region=rules.datasets_bkg_vg250_region_create.output,
        tech_data=rules.datasets_technology_data_copy_files.output,
        wind_stats=rules.datasets_bnetza_mastr_wind_region_create_power_stats_muns.output,
        wind_area_stats=get_abs_dataset_path(
            "datasets", "potentialarea_wind_region") / "data" /
            "potentialarea_wind_area_stats_muns.csv",
        pv_ground_stats=rules.datasets_bnetza_mastr_pv_ground_region_create_power_stats_muns.output,
        pv_ground_area_stats=rules.datasets_potentialarea_pv_ground_region_create_area_stats_muns.output,
        pv_ground_area_shares=rules.datasets_potentialarea_pv_ground_region_create_potarea_shares.output,
        pv_ground_targets=rules.datasets_potentialarea_pv_ground_region_regionalize_state_targets.output,
        pv_roof_stats=rules.datasets_bnetza_mastr_pv_roof_region_create_power_stats_muns.output,
        pv_roof_area_stats=rules.datasets_potentialarea_pv_roof_region_create_area_stats_muns.output,
        pv_roof_area_deploy_stats=rules.datasets_potentialarea_pv_roof_region_create_relative_deployment_stats_muns.output,
        pv_roof_targets=rules.datasets_potentialarea_pv_roof_region_regionalize_state_targets.output,
        hydro_stats=rules.datasets_bnetza_mastr_hydro_region_create_power_stats_muns.output,
        demand_hh_power=rules.datasets_demand_electricity_region_hh_merge_demand_years.output.demand,
        demand_cts_power=rules.datasets_demand_electricity_region_cts_merge_demand_years.output.demand,
        demand_ind_power=rules.datasets_demand_electricity_region_ind_merge_demand_years.output.demand,
        storage_large_stats=rules.datasets_bnetza_mastr_storage_region_create_power_stats_muns.output.large,
        #storage_small_stats=rules.datasets_bnetza_mastr_storage_region_create_power_stats_muns.output.small,
        storage_pv_roof=rules.datasets_bnetza_mastr_storage_region_create_storage_pv_roof_stats.output
    output:
        expand(
            DATASET_PATH / "{panel}_settings_panel.json",
            panel=["energy"]#, "heat", "traffic"]
        )
    run:
        print("Creating panel settings...")
        for panel in ["energy"]:#, "heat", "traffic"]:
            panel_settings = PanelSettings(
                **config["panel_settings_templates"][f"energy_settings_panel"]
            )
            panel_settings = generate_energy_panel_data(
                panel_settings,
                region=gpd.read_file(input.region[0]),
                tech_data=load_json(input.tech_data[0]),
                wind_stats=pd.read_csv(input.wind_stats[0]),
                wind_area_stats=pd.read_csv(input.wind_area_stats),
                pv_ground_stats=pd.read_csv(input.pv_ground_stats[0]),
                pv_ground_area_stats=pd.read_csv(input.pv_ground_area_stats[0], index_col="municipality_id"),
                pv_ground_area_shares=load_json(input.pv_ground_area_shares[0]),
                pv_ground_targets=load_json(input.pv_ground_targets[0]),
                pv_roof_stats=pd.read_csv(input.pv_roof_stats[0]),
                pv_roof_area_stats=pd.read_csv(input.pv_roof_area_stats[0], index_col="municipality_id"),
                pv_roof_area_deploy_stats=pd.read_csv(input.pv_roof_area_deploy_stats[0]),
                pv_roof_targets=load_json(input.pv_roof_targets[0]),
                hydro_stats=pd.read_csv(input.hydro_stats[0]),
                demand_hh_power=pd.read_csv(input.demand_hh_power, index_col="municipality_id"),
                demand_cts_power=pd.read_csv(input.demand_cts_power, index_col="municipality_id"),
                demand_ind_power=pd.read_csv(input.demand_ind_power, index_col="municipality_id"),
                storage_large_stats=pd.read_csv(input.storage_large_stats),
                #storage_small_stats=pd.read_csv(input.storage_small_stats),
                storage_pv_roof=load_json(input.storage_pv_roof[0]),
            )

            import pdb
            pdb.set_trace()

            # with open(
            #         DATASET_PATH / f"{panel}_settings_panel.json",
            #         "w",
            #         encoding="utf8"
            # ) as f:
            #     json.dump(
            #         config["panel_settings_templates"][f"{panel}_settings_panel"],
            #         f,
            #         indent=4
            #     )
