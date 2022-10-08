"""
Snakefile for this dataset

The file will be automatically detected and included in the Snakemake workflow.
"""

from digipipe.store.utils import get_abs_dataset_path

DATASET_PATH = get_abs_dataset_path("2_datasets", "bkg_vg250_muns_abw", data_dir=False)

rule bkg_vg250_muns_abw:
    """
    Extract municipalities of ABW region
    """
    input: get_abs_dataset_path("1_preprocessed", "bkg_vg250") / "bkg_vg250.gpkg"
    output: get_abs_dataset_path("2_datasets", "bkg_vg250_muns_abw") / "bkg_vg250_muns_abw.gpkg"
    params:
        script=DATASET_PATH / "scripts" / "create.py",
        config_path= DATASET_PATH / "config.yml"
    shell:
        "python {params.script} {input} {params.config_path} {output}"
