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
    replace_prefix: {"osm_": "hoho_"}
use rule * from osm_filtered as preprocessed_osm_filtered_*
