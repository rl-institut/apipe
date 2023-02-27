rule create_empty_scalars:
    input: "esys/scenarios/"
    output: "store/datasets/esys_raw/scalars/empty_scalars.csv"
    shell: "python esys/scripts/create_empty_scalars.py {input} {output}"