"""
Snakefile for this dataset

Note: To include the file in the main workflow, it must be added to the respective module.smk .
"""

from digipipe.store.utils import get_abs_dataset_path

DATASET_PATH = get_abs_dataset_path("preprocessed", "ba_employment")

rule unzip:
    """
    Unzip file
    """
    input:
        get_abs_dataset_path("raw", "ba_employment") / "data" / "gemband-dlk-0-202206-zip.zip"
    output:
        files = DATASET_PATH / "data" / "gemband_dlk_0.xlsb"
    params:
        outpath=DATASET_PATH / "data"
    shell:
        """
        unzip -j {input} -d {params.outpath}
        """

rule create:
    """
    Extract municipality data from Excel file and save to CSV
    """
    input:
        DATASET_PATH / "data" / "gemband_dlk_0.xlsb"
    output:
        DATASET_PATH / "data" / "employees.csv"
    script:
        DATASET_PATH / "scripts" / "create.py"
