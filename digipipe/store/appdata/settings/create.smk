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
    add_electricity_panel_settings,
    add_heat_panel_settings,
    add_traffic_panel_settings
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
        storage_pv_roof=rules.datasets_bnetza_mastr_storage_region_create_storage_pv_roof_stats.output,

        heating_structure_decentral=rules.datasets_demand_heat_region_heating_structure_hh_cts.output.heating_structure_esys_dec,
        demand_hh_heat=get_abs_dataset_path("datasets", "demand_heat_region") / "data" / "demand_hh_heat_demand.csv",
        demand_cts_heat=get_abs_dataset_path("datasets", "demand_heat_region") / "data" / "demand_cts_heat_demand.csv",
        demand_ind_heat=get_abs_dataset_path("datasets", "demand_heat_region") / "data" / "demand_ind_heat_demand.csv",

    output:
        panel_settings_electricity=DATASET_PATH / "energy_settings_panel.json",
        panel_settings_heat=DATASET_PATH / "heat_settings_panel.json",
        panel_settings_traffic=DATASET_PATH / "traffic_settings_panel.json"
    run:
        print("Creating electricity panel settings...")
        panel_settings_electricity = PanelSettings(
            name="panel_settings_electricity",
            **config["panel_settings_templates"]["energy_settings_panel"]
        )
        panel_settings_electricity = add_electricity_panel_settings(
            panel_settings_electricity,
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

        print("Creating heat panel settings...")
        panel_settings_heat = PanelSettings(
            name="panel_settings_heat",
            **config["panel_settings_templates"]["heat_settings_panel"]
        )
        panel_settings_heat = add_heat_panel_settings(
            panel_settings_heat,
            heating_structure_decentral=pd.read_csv(input.heating_structure_decentral, index_col="year"),
            demand_hh_heat=pd.read_csv(input.demand_hh_heat, index_col="municipality_id"),
            demand_cts_heat=pd.read_csv(input.demand_cts_heat, index_col="municipality_id"),
            demand_ind_heat=pd.read_csv(input.demand_ind_heat, index_col="municipality_id"),
        )

        print("Creating traffic panel settings...")
        panel_settings_traffic = PanelSettings(
            name="panel_settings_traffic",
            **config["panel_settings_templates"]["traffic_settings_panel"]
        )
        panel_settings_traffic = add_traffic_panel_settings(
            panel_settings_traffic,
        )

        # Check and dump
        for panel, file in zip(
            [panel_settings_electricity,
             panel_settings_heat,
             panel_settings_traffic],
            [output.panel_settings_electricity,
             output.panel_settings_heat,
             output.panel_settings_traffic]
        ):
            if panel.is_complete():
                with open(file, "w", encoding="utf8") as f:
                    json.dump(panel.settings, f, indent=4)
            else:
                raise ValueError(f"{panel.name} has missing values!")
