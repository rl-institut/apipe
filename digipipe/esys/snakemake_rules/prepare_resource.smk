rule prepare_scalars:
    input:
        raw_scalars="store/datasets/esys_raw/scalars/costs_efficiencies.csv",
    output: "store/appdata/esys/_resources/scal_costs_efficiencies.csv"
    shell: "python esys/scripts/prepare_scalars.py {input.raw_scalars} {output}"