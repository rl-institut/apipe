"""
Dataset registry for datasets module which is loaded by main snakemake file.
All datasets in the datasets category must be added to this file.

Template:
---------
module <DATASET_NAME>:
    snakefile: "<DATASET_NAME>/create.smk"
    config: config["store"]["datasets"]["<DATASET_NAME>"]
use rule * from <DATASET_NAME> as datasets_<DATASET_NAME>_*

"""

module bkg_vg250_districts_region:
    snakefile: "bkg_vg250_districts_region/create.smk"
    config: config["store"]["datasets"]["bkg_vg250_districts_region"]
use rule * from bkg_vg250_districts_region as datasets_bkg_vg250_districts_region_*

module bkg_vg250_muns_region:
    snakefile: "bkg_vg250_muns_region/create.smk"
    config: config["store"]["datasets"]["bkg_vg250_muns_region"]
use rule * from bkg_vg250_muns_region as datasets_bkg_vg250_muns_region_*

module bkg_vg250_state:
    snakefile: "bkg_vg250_state/create.smk"
    config: config["store"]["datasets"]["bkg_vg250_state"]
use rule * from bkg_vg250_state as datasets_bkg_vg250_state_*

module bkg_vg250_federal_states:
    snakefile: "bkg_vg250_federal_states/create.smk"
    config: config["store"]["datasets"]["bkg_vg250_federal_states"]
use rule * from bkg_vg250_federal_states as datasets_bkg_vg250_federal_states_*

module bkg_vg250_region:
    snakefile: "bkg_vg250_region/create.smk"
    config: config["store"]["datasets"]["bkg_vg250_region"]
use rule * from bkg_vg250_region as datasets_bkg_vg250_region_*

module osm_forest:
    snakefile: "osm_forest/create.smk"
    config: config["store"]["datasets"]["osm_forest"]
use rule * from osm_forest as datasets_osm_forest_*

module bnetza_mastr_wind_region:
    snakefile: "bnetza_mastr_wind_region/create.smk"
    config: config["store"]["datasets"]["bnetza_mastr_wind_region"]
use rule * from bnetza_mastr_wind_region as datasets_bnetza_mastr_wind_region_*

module bnetza_mastr_pv_ground_region:
    snakefile: "bnetza_mastr_pv_ground_region/create.smk"
    config: config["store"]["datasets"]["bnetza_mastr_pv_ground_region"]
use rule * from bnetza_mastr_pv_ground_region as datasets_bnetza_mastr_pv_ground_region_*

module bnetza_mastr_pv_roof_region:
    snakefile: "bnetza_mastr_pv_roof_region/create.smk"
    config: config["store"]["datasets"]["bnetza_mastr_pv_roof_region"]
use rule * from bnetza_mastr_pv_roof_region as datasets_bnetza_mastr_pv_roof_region_*

module bnetza_mastr_biomass_region:
    snakefile: "bnetza_mastr_biomass_region/create.smk"
    config: config["store"]["datasets"]["bnetza_mastr_biomass_region"]
use rule * from bnetza_mastr_biomass_region as datasets_bnetza_mastr_biomass_region_*

module bnetza_mastr_hydro_region:
    snakefile: "bnetza_mastr_hydro_region/create.smk"
    config: config["store"]["datasets"]["bnetza_mastr_hydro_region"]
use rule * from bnetza_mastr_hydro_region as datasets_bnetza_mastr_hydro_region_*

module bnetza_mastr_combustion_region:
    snakefile: "bnetza_mastr_combustion_region/create.smk"
    config: config["store"]["datasets"]["bnetza_mastr_combustion_region"]
use rule * from bnetza_mastr_combustion_region as datasets_bnetza_mastr_combustion_region_*

module bnetza_mastr_gsgk_region:
    snakefile: "bnetza_mastr_gsgk_region/create.smk"
    config: config["store"]["datasets"]["bnetza_mastr_gsgk_region"]
use rule * from bnetza_mastr_gsgk_region as datasets_bnetza_mastr_gsgk_region_*

module bnetza_mastr_storage_region:
    snakefile: "bnetza_mastr_storage_region/create.smk"
    config: config["store"]["datasets"]["bnetza_mastr_storage_region"]
use rule * from bnetza_mastr_storage_region as datasets_bnetza_mastr_storage_region_*

module bnetza_mastr_captions:
    snakefile: "bnetza_mastr_captions/create.smk"
    config: config["store"]["datasets"]["bnetza_mastr_captions"]
use rule * from bnetza_mastr_captions as datasets_bnetza_mastr_captions_*

module population_region:
    snakefile: "population_region/create.smk"
    config: config["store"]["datasets"]["population_region"]
use rule * from population_region as datasets_population_region_*

module employment_region:
    snakefile: "employment_region/create.smk"
    config: config["store"]["datasets"]["employment_region"]
use rule * from employment_region as datasets_employment_region_*

module demand_electricity_region:
    snakefile: "demand_electricity_region/create.smk"
    config: config["store"]["datasets"]["demand_electricity_region"]
use rule * from demand_electricity_region as datasets_demand_electricity_region_*
