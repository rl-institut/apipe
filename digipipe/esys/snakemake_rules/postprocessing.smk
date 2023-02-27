rule postprocess:
    input: "store/appdata/esys/{scenario}/optimized"
    output: directory("store/appdata/esys/{scenario}/postprocessed/")
    params:
        logfile="store/appdata/esys/{scenario}/{scenario}.log"
    shell: "python esys/scripts/postprocess.py {input} {wildcards.scenario} {output} {params.logfile}"
