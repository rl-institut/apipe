from apipe.esys.esys.config.esys_conf import load_yaml

def get_paths_scenario_input(wildcards):
    scenario_specs = load_yaml(f"esys/scenarios/{wildcards.scenario}.yml")
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
        get_paths_scenario_input,
        scenario="esys/scenarios/{scenario}.yml"
    output: directory("store/appdata/esys/{scenario}/preprocessed")
    params:
        logfile="store/appdata/esys/{scenario}/{scenario}.log"
    wildcard_constraints:
        # Do not use this rule for the examples. Use prepare_example instead
        scenario=r"(?!example_).*"
    shell: "python esys/scripts/build_datapackage.py {input.scenario} {output} {params.logfile}"
