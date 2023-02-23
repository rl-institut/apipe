import os
from digipipe.store.utils import get_abs_dataset_path

APPDATA_ESYS_PATH = get_abs_dataset_path("appdata", "esys")

rule optimize:
    input: os.path.join(APPDATA_ESYS_PATH, "{scenario}", "preprocessed")
    output: directory(os.path.join(APPDATA_ESYS_PATH, "{scenario}", "optimized"))
    params:
        logfile=os.path.join(APPDATA_ESYS_PATH, "{scenario}", "{scenario}.log")
    shell: "python scripts/optimize.py {input} {output} {params.logfile}"
