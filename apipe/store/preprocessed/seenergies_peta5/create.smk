"""
Snakefile for this dataset

Note: To include the file in the main workflow, it must be added to the respective module.smk .
"""

import zipfile
from apipe.store.utils import get_abs_dataset_path

DATASET_PATH = get_abs_dataset_path("preprocessed", "seenergies_peta5", data_dir=True)

rule extract:
    """
    Extract files from archives
    """
    input:
        residential=get_abs_dataset_path(
            "raw", "seenergies_peta5") / "data" / "Peta5_0_1_HD_res.zip",
        cts=get_abs_dataset_path(
            "raw", "seenergies_peta5") / "data" / "Peta5_0_1_HD_ser.zip"
    output:
        residential=[DATASET_PATH / f for f in config["files_extract"]["residential"]],
        cts=[DATASET_PATH/f for f in config["files_extract"]["cts"]]
    params:
        outpath=DATASET_PATH,
        files_extract_residential=config["files_extract"]["residential"],
        files_extract_cts=config["files_extract"]["cts"]
    run:
        # note: buildin zipfile used as bash unzip doesn't work (pkzip required)
        with zipfile.ZipFile(input.residential, "r") as zf:
            for f in params.files_extract_residential:
                zf.extract(f, params.outpath)
        with zipfile.ZipFile(input.cts, "r") as zf:
            for f in params.files_extract_cts:
                zf.extract(f, params.outpath)
