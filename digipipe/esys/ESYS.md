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

Empty scalars and time series can be created from the energy model setup with

```
snakemake -j1 create_empty_scalars
snakemake -j1 create_empty_ts
```