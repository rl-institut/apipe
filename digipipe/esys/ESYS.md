# Energy system in digiplan

The energy system in digipipe, the pipeline for digiplan, is created using
[oemof-B3](https://github.com/rl-institut/oemof-B3)`.

## Build the energy system

To test if everything works, you can run the test scenario with

```
snakemake -j1 make_esys_appdata
```
For this you have to copy the corresponding raw data into the raw directory.
In the future, empty raw data (scalars and time series) will be created
automatically. Then, assumptions on constant parameters such as plant costs,
lifetime and efficiencies are mapped and set as values of the corresponding
variables in the scalars.

To set up an empty energy system, the following prompt automatically writes
default scalar values (such as zero or NaN) to the empty scalars:

```
snakemake -j1 write_default_scalars
```
With this the file `empty_scalars.csv` is automatically updated and saved to:
`digipipe/digipipe/store/datasets/esys_raw/data/scalars/default_scalars.csv`
