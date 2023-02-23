import os
from digipipe.store.utils import get_abs_dataset_path

DATASET_ESYS_RAW_PATH = get_abs_dataset_path("datasets", "esys_raw")
APPDATA_ESYS_PATH = get_abs_dataset_path("appdata", "esys")

rule prepare_scalars:
    input:
        raw_scalars=os.path.join(DATASET_ESYS_RAW_PATH, "scalars", "costs_efficiencies.csv"),
    output: os.path.join(APPDATA_ESYS_PATH, "_resources", "scal_costs_efficiencies.csv")
    shell: "python scripts/prepare_scalars.py {input.raw_scalars} {output}"
