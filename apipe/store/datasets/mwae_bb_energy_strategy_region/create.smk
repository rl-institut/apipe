"""
Snakefile for this dataset

Note: To include the file in the main workflow, it must be added to the respective module.smk .
"""
import json
import geopandas as gpd

from apipe.store.utils import get_abs_dataset_path

DATASET_PATH = get_abs_dataset_path(
    "datasets", "mwae_bb_energy_strategy_region", data_dir=True)

rule copy_files:
    """
    Copy technology data from raw
    """
    input:
        get_abs_dataset_path(
            "raw", "mwae_bb_energy_strategy"
        ) / "data" / "mwae_bb_energy_strategy.json"
    output:
        DATASET_PATH / "mwae_bb_energy_strategy.json"
    shell:
        """
        cp -p {input} {output}
        """

rule regionalize_state_targets:
    """
    Calculate targets for region (by area share)
    """
    input:
        targets=DATASET_PATH / "mwae_bb_energy_strategy.json",
        region=rules.datasets_bkg_vg250_region_create.output,
        states=rules.preprocessed_bkg_vg250_create.output
    output:
        DATASET_PATH / "mwae_bb_energy_strategy_region.json"
    run:
        states = gpd.read_file(input.states[0], layer="vg250_lan")
        area_state = states.loc[states["NUTS"] == "DE4"].area.sum()
        area_region = gpd.read_file(input.region[0]).area.sum()
        area_share = area_region / area_state

        with open(input.targets, "r") as f:
            targets = json.load(f)

        targets_region = dict()
        skip_keys = [
            "re_share_in_demand_electricity_perc",
            "re_share_in_demand_heat_perc"
        ]
        # Rescale targets with area share
        for sec1, val1 in targets.items():
            if sec1 in skip_keys:
                targets_region.update({sec1: val1})
            else:
                targets_region[sec1]= dict()
                for sec2, val2 in val1.items():
                    targets_region[sec1].update(
                        {sec2: round(val2 * area_share)}
                    )
        with open(output[0], "w", encoding="utf8") as f:
            json.dump(targets_region, f, indent=4)
