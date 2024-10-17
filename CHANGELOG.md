# Changelog
All notable changes to this project will be documented in this file.

The format is inspired from [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and the versioning aim to respect [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- Add PV ground potentials for Germany (dataset `potentialarea_pv_ground`) and
  regionalization (dataset `potentialarea_pv_ground_region2`)
- Add function to reduce DF size
- Add protected areas for Germany (dataset `bfn_protected_areas`) and
  regionalization
- Add RLG OLS data (dataset `rpg_ols_regional_plan`)
- Add PV roof potentials for Brandenburg (dataset `wfbb_pv_roof_potential`) and
  regionalization (dataset `potentialarea_pv_roof2`)
- Add energy strategy targets for Brandenburg (dataset
  `mwae_bb_energy_strategy`) and regionalization (dataset
  `mwae_bb_energy_strategy_region`)
- Add permanent crops Brandenburg (dataset `mluk_bb_field_block_cadastre`)
- Combine negative PV ground criteria layers to additional layer
- Buffer negative PV ground criteria layer for open spaces to prevent hatch
  style

### Changed

- Regionalize PV ground and roof targets from energy strategy BB
- Replace permanent crops dataset `oei_agri_pv` by
  `mluk_bb_field_block_cadastre` for PV ground potential calculation

### Fixed
