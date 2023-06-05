"""
Dataset registry for preprocessed datasets module which is loaded by main
snakemake file. All datasets in the preprocessed category must be added to
this file.

Template:
---------
module <DATASET_NAME>:
    snakefile: "<DATASET_NAME>/create.smk"
    config: config["store"]["preprocessed"]["<DATASET_NAME>"]
use rule * from <DATASET_NAME> as preprocessed_<DATASET_NAME>_*

"""

module bkg_vg250:
    snakefile: "bkg_vg250/create.smk"
    config: config["store"]["preprocessed"]["bkg_vg250"]
use rule * from bkg_vg250 as preprocessed_bkg_vg250_*

module osm_filtered:
    snakefile: "osm_filtered/create.smk"
    config: config["store"]["preprocessed"]["osm_filtered"]
use rule * from osm_filtered as preprocessed_osm_filtered_*

module bnetza_mastr:
    snakefile: "bnetza_mastr/create.smk"
    config: config["store"]["preprocessed"]["bnetza_mastr"]
use rule * from bnetza_mastr as preprocessed_bnetza_mastr_*

module destatis_gv:
    snakefile: "destatis_gv/create.smk"
    config: config["store"]["preprocessed"]["destatis_gv"]
use rule * from destatis_gv as preprocessed_destatis_gv_*

module stala_st_pop_prog:
    snakefile: "stala_st_pop_prog/create.smk"
    config: config["store"]["preprocessed"]["stala_st_pop_prog"]
use rule * from stala_st_pop_prog as preprocessed_stala_st_pop_prog_*

module demandregio:
    snakefile: "demandregio/create.smk"
    config: config["store"]["preprocessed"]["demandregio"]
use rule * from demandregio as preprocessed_demandregio_*

module ba_employment:
    snakefile: "ba_employment/create.smk"
    config: config["store"]["preprocessed"]["ba_employment"]
use rule * from ba_employment as preprocessed_ba_employment_ *

module bmwk_long_term_scenarios:
    snakefile: "bmwk_long_term_scenarios/create.smk"
    config: config["store"]["preprocessed"]["bmwk_long_term_scenarios"]
use rule * from bmwk_long_term_scenarios as preprocessed_bmwk_long_term_scenarios_ *

module seenergies_peta5:
    snakefile: "seenergies_peta5/create.smk"
    config: config["store"]["preprocessed"]["seenergies_peta5"]
use rule * from seenergies_peta5 as preprocessed_seenergies_peta5_ *

module ageb_energy_balance:
    snakefile: "ageb_energy_balance/create.smk"
    config: config["store"]["preprocessed"]["ageb_energy_balance"]
use rule * from ageb_energy_balance as preprocessed_ageb_energy_balance_ *

module dwd_temperature:
    snakefile: "dwd_temperature/create.smk"
    config: config["store"]["preprocessed"]["dwd_temperature"]
use rule * from dwd_temperature as preprocessed_dwd_temperature_ *

module regiostat:
    snakefile: "regiostat/create.smk"
    config: config["store"]["preprocessed"]["regiostat"]
use rule * from regiostat as preprocessed_regiostat_ *

module eurostat_lau:
    snakefile: "eurostat_lau/create.smk"
    config: config["store"]["preprocessed"]["eurostat_lau"]
use rule * from eurostat_lau as preprocessed_eurostat_lau_ *
