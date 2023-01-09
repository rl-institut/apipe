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

module bkg_vg250_districts_abw:
    snakefile: "bkg_vg250_districts_abw/create.smk"
    config: config["store"]["datasets"]["bkg_vg250_districts_abw"]
use rule * from bkg_vg250_districts_abw as datasets_bkg_vg250_districts_abw_*

module bkg_vg250_muns_abw:
    snakefile: "bkg_vg250_muns_abw/create.smk"
    config: config["store"]["datasets"]["bkg_vg250_muns_abw"]
use rule * from bkg_vg250_muns_abw as datasets_bkg_vg250_muns_abw_*

module bkg_vg250_region_abw:
    snakefile: "bkg_vg250_region_abw/create.smk"
    config: config["store"]["datasets"]["bkg_vg250_region_abw"]
use rule * from bkg_vg250_region_abw as datasets_bkg_vg250_region_abw_*

module osm_forest:
    snakefile: "osm_forest/create.smk"
    config: config["store"]["datasets"]["osm_forest"]
use rule * from osm_forest as datasets_osm_forest_*

module bnetza_mastr_wind_abw:
    snakefile: "bnetza_mastr_wind_abw/create.smk"
    config: config["store"]["datasets"]["bnetza_mastr_wind_abw"]
use rule * from bnetza_mastr_wind_abw as datasets_bnetza_mastr_wind_abw_*
