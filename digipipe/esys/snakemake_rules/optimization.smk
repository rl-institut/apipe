rule optimize:
    input: "store/appdata/esys/{scenario}/preprocessed"
    output: directory("store/appdata/esys/{scenario}/optimized/")
    params:
        logfile="store/appdata/esys/{scenario}/{scenario}.log"
    shell: "python esys/scripts/optimize.py {input} {output} {params.logfile}"
