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
