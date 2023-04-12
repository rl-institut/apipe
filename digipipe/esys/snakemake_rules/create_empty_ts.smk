rule create_empty_ts:
    input: "esys/scenarios/"
    output:
        "store/datasets/esys_raw/data/time_series/empty_ts_load.csv",
        "store/datasets/esys_raw/data/time_series/empty_ts_feedin.csv",
        "store/datasets/esys_raw/data/time_series/empty_ts_efficiencies.csv"

    shell: "python esys/scripts/create_empty_ts.py {input} {output[0]} {output[1]} {output[2]}"