"""
Dataset registry for appdata module which is loaded by main snakemake file.
All datasets in the datasets category must be added to this file.

Template:
---------
module <DATASET_NAME>:
    snakefile: "<DATASET_NAME>/create.smk"
    config: config["store"]["appdata"]["<DATASET_NAME>"]
use rule * from <DATASET_NAME> as appdata_<DATASET_NAME>_*

"""

module datapackage:
    snakefile: "datapackage/create.smk"
    config: config["store"]["appdata"]["datapackage"]
use rule * from datapackage as appdata_datapackage_*

module settings:
    snakefile: "settings/create.smk"
    config: config["store"]["appdata"]["settings"]
use rule * from settings as appdata_settings_*
