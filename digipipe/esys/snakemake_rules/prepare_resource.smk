rule prepare_scalars:
    input:
        raw_scalars="raw/scalars/costs_efficiencies.csv",
    output: "results/_resources/scal_costs_efficiencies.csv"
    shell: "python scripts/prepare_scalars.py {input.raw_scalars} {output}"
