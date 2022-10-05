"""
Snakefile for this dataset

The file will be automatically detected and included in the Snakemake workflow.
Please add a docstring with a short description to each rule.
"""
from pathlib import Path
from digipipe.store.utils import get_abs_dataset_path

configfile: get_abs_dataset_path("1_preprocessed", ".TEMPLATE", data_dir=False) / "config.yml"

# Use function `get_abs_dataset_path()` to get path to dataset
rule template_simple_copy_file1:
    input: get_abs_dataset_path("0_raw", ".TEMPLATE") / "some_timeseries.csv"
    output: get_abs_dataset_path("1_preprocessed", ".TEMPLATE") / "some_timeseries.csv"
    shell: "cp -p {input} {output}"

# Or you can use relative paths which may break as snakemake can be invoked
# from directories digipipe/ or digipipe/workflow/.
rule template_simple_copy_file2:
    input: Path(".") / "store" / "0_raw" / ".TEMPLATE" / "data" /  "some_timeseries.csv"
    output: Path(".") / "store" / "1_preprocessed" / ".TEMPLATE" / "data" / "some_timeseries.csv"
    shell: "cp -p {input} {output}"
