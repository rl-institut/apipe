import os
from digipipe.store.utils import get_abs_dataset_path

APPDATA_ESYS_PATH = get_abs_dataset_path("appdata", "esys")

rule postprocess:
    input: os.path.join(APPDATA_ESYS_PATH, "{scenario}", "optimized")
    output: directory(os.path.join(APPDATA_ESYS_PATH, "{scenario}", "postprocessed/"))
    params:
        logfile= os.path.join(APPDATA_ESYS_PATH, "{scenario}", "{scenario}.log")
    shell: "python scripts/postprocess.py {input} {wildcards.scenario} {output} {params.logfile}"

def get_scenarios_in_group(wildcards):
    return [os.path.join(APPDATA_ESYS_PATH, scenario, "postprocessed") for scenario in scenario_groups[wildcards.scenario_group]]

rule join_scenario_results:
    input: get_scenarios_in_group
    output: directory(os.path.join(APPDATA_ESYS_PATH, "joined_scenarios", "{scenario_group}","joined/"))
    shell: "python scripts/join_scenarios.py {input} {output}"

rule map_results_to_b3_format:
    input:
        os.path.join(APPDATA_ESYS_PATH, "{scenario}", "postprocessed")
    output:
        directory(os.path.join(APPDATA_ESYS_PATH, "{scenario}", "b3_results/data"))
    params:
        logfile=os.path.join(APPDATA_ESYS_PATH, "{scenario}", "{scenario}.log")
    shell:
        "python scripts/map_results_to_b3_format.py {input} {output} {params.logfile}"
