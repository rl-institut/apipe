# Energy system in digiplan

The energy system in digipipe, the pipeline for digiplan, is created using 
[oemof-B3](https://github.com/rl-institut/oemof-B3)`.
 
## Build the energy system

To test if everything works, you can run the test scenario with 

```
snakemake -j1 results/Test_scenario/preprocessed
```
For this you have to copy the corresponding raw data into the raw directory.