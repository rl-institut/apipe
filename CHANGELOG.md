# Changelog
All notable changes to this project will be documented in this file.

The format is inspired from [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and the versioning aim to respect [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- Setup initial structure and files
- Add first bunch of datasets for testing the data flow
- Modularize datasets
- Add first draft of RTD docs
- Add dataset: BNetzA MaStR
- Add Nominatim geocoder
- Add dataset: population
- Clean rule
- Datasets attribute captions
- Create list of region-specific datasets in the docs
- pyproject.toml and poetry.lock file have been added with the conversion to poetry
- Add pre-commit in order to check for errors and linting bugs before commits
- Add types-pyyaml package
- Add dataset: employees and companies
- Add dataset: demandregio electricity demand
- Add dataset: BMWK long-term scenarios
- Add rules to download raw data (zipped) from cloud, extract and copy to 'store/raw'
- Add module 'data_io', containing relevant helper functions (downloading, extracting, copying, cleaning)
- Create metadata
- Add script to generate dataset md files for documentation
- Add dataset: demand_heat_region
- Add dataset: heatpump_cop
- Add dataset: stala_st_energy
- Add dataset: eurostat_lau
- Add dataset: regstat_energy
- Add dataset: dwd_temperature
- Add dataset: ageb_energy_balance
- Add dataset: seenergies_peta5
- Add dataset: renewables.ninja_feedin
- Add dataset: renewable_feedin
- Add dataset: bnetza_mastr_correction_region and correct wrong units
- Integrate building of energy system for appdata in pipeline via dir *esys*
- Update store with dir structure for *esys* data
- Add creation of empty time series for the *esys*
- Add dataset: rpg_abw_regional_plan
- Add dataset: potentialarea_wind_region
- Add writing of default values to *esys* raw scalar data
- Add datasets: rli_pv_wfr and rli_pv_wfr_region
- Add module appdata to workflow
- Add dataset: geodata_infolayers
- Add dataset: potentialarea_pv_ground_region
- Add dataset: app datapackage
- Add dataset: potentialarea_pv_roof_region
- Add dataset: technology_data
- Add dataset: settings
- Calc panel settings from datasets
- Add dataset: osm (Germany)
- Add dataset: osm_buildings and add stats on ground areas
- Add mapping of costs and efficiencies from store/raw to store/datasets
- Add dataset: emissions
- Add captions to app datapackage (here: MaStR, heat, potentialarea_wind)
- Add mapping of time series data in datasets to empty time series according to
  the mapping provided in map_ts.yml
- Add build configuration for readthedocs
- Add creation of stats of development over time for bnetza_mastr_wind_region,
  bnetza_mastr_pv_ground_region, bnetza_mastr_pv_roof_region
- Add dataset: dbfz_biomass_heat_capacities
- Add the calculation of relative demand of biomass conversion technologies via
  their relative capacities
- Add deletion of all data in store/datasets/esys_raw/data
- Add notes on OSM download and run resources
- Add nominal power per wind turbine for 2045
- Add technology data for batteries
- Add technology data for thermal storages

### Changed

- Move dataset docs from config to md files
- Retain mastr_id in MaStR datasets
- Fix loading of empty yml files
- Fix loading of global config when in workflow dir
- Integrate esys Snakefile in workflow Snakefile and update clean rule
- Fix shapely deprecation warning
- Fix ogr2ogr conversion with recent GDAL version (v3.6.2)
- Fix conda installation by removing python-gdal from environment.yml
- The package management in digipipe has been changed to poetry.
- The installation of a virtual environment is done only from the environment.yml file and via conda.
- Apply linters on repo among others: black, isort, check-json and end-of-file-fixer
- Update population with prognoses from demandregio dataset
- Fix C419 flake8 error
- Switch to mkdocs for documentation (Sphinx deleted)
- Normalize renewable feedin timeseries
- Fix instruction to obtain raw files
- Translate all dataset.md files to German
- Exchange *Test_scenario* with *2045_scenario* in *digipipe/esys/scenarios*
- Split each demand per sector in *esys*
- File .gitignore again includes ignoring of esys appdata
- pv_roof area stats: distinguish between all and non-historic buildings
- storage units: add region-wide values for spec. capacity and power for those
  connected to PV roof units
- Add data on installed el. power to bmwk_long_term_scenarios
- Disaggregate PV state targets to region in potentialarea_pv_ground_region
- Adapt osm_filtered to ose osm dataset and extract building data
- Disaggregate PV state targets to region in potentialarea_pv_roof_region
- Changes were applied to the energy system. Among others RoR, small batteries
  and biogas were added. A distinction was made between centralized and
  decentralized CHPs
- Scenario 2045_scenario needs default_scalars.csv instead of scalars.csv
- By default set costs and efficiencies of esys are written to
  default_scalars.csv instead of default_costs_efficiencies.csv
- Default variable_costs are passed with input_parameters for storages
- Pass time series instead of scalar with efficiency for central heat pump
- Fix wind+pv_ground default values in panel settings
- Set all default control values in panel settings
- Kick biogas shortage
- Rename dataset captions to app_captions
- Move app settings to datasets and include in app datapackage
- Adapt 2045_scenario.yml so that time series with values are used instead of
  empty ts
- Suppress warning of loosing data in source and comment columns while
  unstacking if they are empty
- Change max. installable PV roof capacity in panel settings
- Fix panel settings for large batteries
- Add additional captions to MaStR captions
- Use LTS version of OSM
- The unstacking of time series in esys was fixed so that warning is given if
  there is at least one value in columns 'source' or 'comment'
- Minor fix applied reformatting with black
- Only use operating units from mastr for municipality stats and temporal
  development
- Heat pump ASHP/GSHP split fixed
- Replace the relative demand of biomass with the relative demand of each
  biomass conversion technology
- Fix clean rule
- Update raw datapackage URL
- Restrict snakemake version to v7.32.0
- Add central heat pump targets to slider
- Restrict heat pump sliders to not move under 50%
- Fix pv ground slider values to prevent app to alter SQ value from panel
  settings
- Fix PV roof slider values
- Add HP share slider from-max values to prevent 100 % HP share
- Updated technology_data dataset.md and metadata

### Removed

- setup.py and requirements.txt files are omitted with the conversion to poetry
- sphinx from poetry environment
- Remove dataset: osm_forest
- Obsolete targets from rule all
- Merge dataset costs_efficiencies into technology_data
- Merge dataset costs_efficiencies into technology_data
- Remove values for redundant subsliders from app datapackage
