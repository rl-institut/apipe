from digipipe.esys.esys.config.esys_conf import load_yaml
from digipipe.store.utils import get_abs_dataset_path
import os

APPDATA_ESYS_PATH = get_abs_dataset_path("appdata", "esys")


def get_paths_scenario_input(wildcards):
    scenario_specs = load_yaml(f"scenarios/{wildcards.scenario}.yml")
    paths_scenario_inputs = list()
    for key in ["paths_scalars", "paths_timeseries"]:
        paths = scenario_specs[key]
        if isinstance(paths, list):
            paths_scenario_inputs.extend(paths)
        elif isinstance(paths, str):
            paths_scenario_inputs.append(paths)
    return paths_scenario_inputs

rule build_datapackage:
    input:
        scenario="scenarios/{scenario}.yml"
    output: directory(APPDATA_ESYS_PATH/"{scenario}/preprocessed")
    params:
        logfile=os.path.join(APPDATA_ESYS_PATH, "{scenario}", "{scenario}.log")
    wildcard_constraints:
        # Do not use this rule for the examples. Use prepare_example instead
        scenario=r"(?!example_).*"
    shell: "python scripts/build_datapackage.py {input.scenario} {output} {params.logfile}"

