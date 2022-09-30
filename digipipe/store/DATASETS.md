# Data pipeline

## Structure of the data store

**TBD**

## Data flow

This section describes the workflow of the data pipeline.

(RAW) -> (PREPROCESSED) -> (DATASETS) -> (APPDATA)

Overview:

| **Step** | **Directory**         | **Description**                  | **Rule(s) for this target** | **Cfg section** |
|:--------:|-----------------------|----------------------------------|-----------------------------|-----------------|
|    0     | store/0_raw/          | Raw data as downloaded           |                             |                 |
|    1     | store/1_preprocessed/ | Preprocessed data, 1:1 from (0)  |                             |                 |
|    2     | store/2_datasets/     | Datasets, n:1 from (1) and (2)   |                             |                 |
|    3     | store/3_appdata/      | Data ready to be used in the app |                             |                 |

In the following each step is shortly described along a common example use
case.

**Example data flow:**

![example data flow](../../docs/img/datasets/pipeline_dataflow_example.png)

### (0) Raw

Immutable raw data as downloaded with 2 additional files:
[description](0_raw/.TEMPLATE/dataset.md) (see this file for further
instructions) and [metadata](0_raw/.TEMPLATE/metadata.json).

Note: Assumptions are to be defined in the scenarios, not the raw data.
See the scenario readme in [SCENARIOS.md](../scenarios/SCENARIOS.md). 

> **Example:**
> - Dataset A: ERA5 weather dataset for Germany
> - Dataset B: MaStR dataset on renewable generators
> - Dataset C: Shapefile of region of interest

### (1) Preprocessed

Data from `(0) Raw`  that has undergone some preprocesing such as:
 - Archive extracted
 - CRS transformed (see below for CRS conventions)
 - Fields filtered
 - **But NO merging/combining/clipping of multiple (raw) datasets! This can be 
   done in (2)**

Note: The directory name MUST be the same as in `0_raw`.

The preprocessing rules can be defined in the dataset's
[snakemake file](1_preprocessed/.TEMPLATE/create.smk).

> **Example:**
> - Dataset D: Extracted ERA5 weather dataset for Germany (from dataset A)
> - Dataset E: Wind energy turbines extracted from MaStR dataset, filter for
>   columns power and geometry (from dataset B)
> - Dataset F: Region of interest converted to Geopackage file, CRS
>   transformed (from dataset C)

### (2) Datasets

Datasets, created from arbitrary combinations of datasets from
`(1) Preprocessed` and/or `(2) Datasets`.

Notes:
- The creation rules can be defined in the dataset's
[snakemake file](2_datasets/.TEMPLATE/create.smk).
- Custom, dataset-specific configuration can be put into the
[dataset config](2_datasets/.TEMPLATE/config.yml)
- Custom, dataset-specific scripts are located in `scripts`

> **Example:**
> 
> Using datasets from (1) and (2):
> - Dataset G: Wind energy turbines in the region of interest (from datasets E+F)
> - Dataset H: Normalized wind energy feedin timeseries for the region (from
>   datasets D+G)
> - Dataset I: Region of interest (from dataset F)

### (3) App data

**TBD**

Data ready to be used in the app / as expected by the app.

### Temporary files

**TODO: REVISE**

Temporary files are stored in `store/temp/` by default and, depending on your
configuration, can get quite large.  You can change the directory in
`config.yml` -> `path` -> `temp`.

## Further notes

### No data files in the repository!

Make sure **not to commit any data files** located in `store/` to the
repository (except for the descriptive readme and metadata files). They should
be omitted by the rules defined in the `.gitignore` file but you better don't
count on it. Instead, save them in the designated directory on the
[RLI Wolke](https://wolke.rl-institut.de/f/160572).

### Coordinate reference system

Please use LAEA Europe (EPSG:3035) as default CRS when writing geodata.

**TODO: REVISE**

- The files in `store/raw/` can have an arbitrary CRS.
- In the preprocessing (step 1) it is converted to the CRS specified in the global `config.yml` -> `preprocessing` -> 
  `crs`. It is important to use a equal-area CRS to make sure operations such as buffering work properly. By default,
  it is set to LAEA Europe (EPSG:3035).
- The final output is written in CRS specified in the global `config.yml` -> `output` -> `crs`. By default, it is set
  to WGS84 (EPSG:4326) used by the app.

## HowTos

### Add a new dataset

**TODO: REVISE**
