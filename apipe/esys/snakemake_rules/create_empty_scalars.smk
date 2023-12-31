rule create_empty_scalars:
    """
    Create empty scalar files
    """
    input: "esys/scenarios/"
    output: "store/datasets/esys_raw/data/scalars/empty_scalars.csv"
    shell: "python esys/scripts/create_empty_scalars.py {input} {output}"
