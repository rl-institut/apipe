from apipe.esys.esys.config.esys_conf import map_ts
from apipe.store.utils import get_abs_dataset_path


rule write_ts:
    input:
        eff="store/datasets/esys_raw/data/time_series/empty_ts_efficiencies.csv",
        feedin="store/datasets/esys_raw/data/time_series/empty_ts_feedin.csv",
        load="store/datasets/esys_raw/data/time_series/empty_ts_load.csv",
        file_mapping_dummy=[  # Workaround to include data from pipeline
            [get_abs_dataset_path("datasets", d["dataset"]) / "data"/ d["file"]
             for d in map_ts[c].values()]
            for c in ["load", "feedin", "efficiencies"]
        ]
    output:
        "store/datasets/esys_raw/data/time_series/ts_efficiencies.csv",
        "store/datasets/esys_raw/data/time_series/ts_feedin.csv",
        "store/datasets/esys_raw/data/time_series/ts_load.csv",
    shell:
        """
        python esys/scripts/write_ts.py {input.eff} {input.feedin} {input.load} {output}
        """
