"""
Snakefile for this dataset

Note: To include the file in the main workflow, it must be added to the respective module.smk .
"""

from digipipe.store.utils import get_abs_dataset_path

DATASET_PATH = get_abs_dataset_path("preprocessed", "destatis_gv")

rule create:
    """
    Extract municipality population data from Excel files and save to CSVs 
    """
    input:
        #lambda wildcards: expand(get_abs_dataset_path("raw", "destatis_gv") / "data" / "3112{year}_Auszug_GV.xlsx", year=wildcards.year)
        #x=lambda wildcards: get_abs_dataset_path("raw", "destatis_gv") / "data" / f"3112{wildcards.year}_Auszug_GV.xlsx"
        #get_abs_dataset_path("raw", "destatis_gv") / "data" / "31122022_Auszug_GV.xlsx"
        #expand(get_abs_dataset_path("raw", "destatis_gv") / "data" / "3112{year}_Auszug_GV.xlsx", year=[2010, 2015, 2020, 2021, 2022])
        get_abs_dataset_path("raw", "destatis_gv") / "data" / "3112{year}_Auszug_GV.xlsx"
    output:
        #expand(DATASET_PATH / "data" / "3112{year}_Auszug_GV.csv", year=[2010, 2015, 2020, 2021, 2022])
        #x=[DATASET_PATH / "data" / f"3112{year}_Auszug_GV.csv" for year in [2010, 2015, 2020, 2021, 2022]]
        #DATASET_PATH / "data" / "31122022_Auszug_GV.csv"
        DATASET_PATH / "data" / "3112{year}_Auszug_GV.csv"
    script:
        DATASET_PATH / "scripts" / "create.py"
