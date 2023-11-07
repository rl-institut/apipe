rule write_ts:
    input:
        "store/datasets/esys_raw/data/time_series/empty_ts_efficiencies.csv",
        "store/datasets/esys_raw/data/time_series/empty_ts_feedin.csv",
        "store/datasets/esys_raw/data/time_series/empty_ts_load.csv",
    output:
        "store/datasets/esys_raw/data/time_series/ts_efficiencies.csv",
        "store/datasets/esys_raw/data/time_series/ts_feedin.csv",
        "store/datasets/esys_raw/data/time_series/ts_load.csv",
    shell: "python esys/scripts/write_ts.py {input} {output}"
