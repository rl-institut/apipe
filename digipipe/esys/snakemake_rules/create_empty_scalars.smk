import os
from digipipe.store.utils import get_abs_dataset_path

DATASET_ESYS_RAW_PATH = get_abs_dataset_path("datasets", "esys_raw")
APPDATA_ESYS_PATH = get_abs_dataset_path("appdata", "esys")

rule create_empty_scalars:
    input: directory("scenarios/")
    output: os.path.join(DATASET_ESYS_RAW_PATH, "scalars", "empty_scalars.csv")
    shell: "python scripts/create_empty_scalars.py {input} {output}"