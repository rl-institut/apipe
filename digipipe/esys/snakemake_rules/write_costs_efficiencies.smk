rule write_costs_efficiencies:
    """
    Write costs and efficiencies from raw data to esys datasets
    """
    input:
        "store/datasets/esys_raw/data/scalars/default_costs_efficiencies.csv",
        "store/raw/costs_efficiencies/data/raw_costs_efficiencies.csv",
    output: "store/datasets/esys_raw/data/scalars/costs_efficiencies.csv"
    shell: "python esys/scripts/write_costs_efficiencies.py {input} {output}"
