# Regionsspezifischer Datensatz Case Study: Erkner, Rüdersdorf (b Berlin), Strausberg, Grünheide (Mark)

Die meisten der Rohdatensätze können für eine beliebige Region in Deutschland
verwendet werden. Einige sind jedoch nur für eine Teilregion verfügbar oder
spezifisch für die Region  **r120640428428 (Rüdersdorf b Berlin)**, **r120640472472 (Strausberg)**, **r120670124124 (Erkner)**, **r120670201201(Grünheide)** :

## Daten für die Energiesystemmodellierung

### Strombedarf

#### Gesamtstrombedarf (2045), normiert

| Name                     | Raw-Datensatz | Dataset                   | Kommentar |
|--------------------------|---------------|---------------------------|-----------|
| *electricity-demand_hh*  | demandregio   | demand_electricity_region |           |
| *electricity-demand_cts* | demandregio   | demand_electricity_region |           |
| *electricity-demand_ind* | demandregio   | demand_electricity_region |           |
| *electricity-demand_mob* | egon_ev       | demand_emobility_region   |           |
 

*[Platzhalter für weiterführende Erklärungen und Annahmen]*



#### Stromlastprofile, normiert

| Name                             | Raw-Datensatz | Dataset                   | Kommentar |
|----------------------------------|---------------|---------------------------|-----------|
| *electricity-demand_hh_profile*  | demandregio   | demand_electricity_region |           |
| *electricity-demand_cts_profile* | demandregio   | demand_electricity_region |           |
| *electricity-demand_ind_profile* | demandregio   | demand_electricity_region |           |
| *electricity-demand_mob_profile* |               |                           |           |

*[Platzhalter für weiterführende Erklärungen und Annahmen]*

------------------------------------------------------------------------------------------------------------------------
### Wärmebedarf Niedrigtemperaturwärme

#### Gesamtwärmebedarf Niedrigtemperaturwärme zentral (2045) in MWh

| Name                          | Raw-Datensatz      | Dataset             | Kommentar                                                               |
|-------------------------------|--------------------|---------------------|-------------------------------------------------------------------------|
| *heat_low_central-demand_hh*  | wfbb_heat_atlas_bb | demand_heat_region2 | Manueller Checkout von `features/#3_add_dataset_heat_atlas_brandenburg` |
| *heat_low_central-demand_cts* | wfbb_heat_atlas_bb | demand_heat_region2 | Manueller Checkout von `features/#3_add_dataset_heat_atlas_brandenburg` |
| *heat_low_central-demand_ind* | wfbb_heat_atlas_bb | demand_heat_region2 | Manueller Checkout von `features/#3_add_dataset_heat_atlas_brandenburg` | 

BMWK lfs

*[Platzhalter für weiterführende Erklärungen]*

#### Gesamtwärmebedarf Niedrigtemperaturwärme dezentral (2045) in MWh

| Name                            | Raw-Datensatz      | Dataset             | Kommentar                                                               |
|---------------------------------|--------------------|---------------------|-------------------------------------------------------------------------|
| *heat_low_decentral-demand_hh*  | wfbb_heat_atlas_bb | demand_heat_region2 | Manueller Checkout von `features/#3_add_dataset_heat_atlas_brandenburg` |
| *heat_low_decentral-demand_cts* | wfbb_heat_atlas_bb | demand_heat_region2 | Manueller Checkout von `features/#3_add_dataset_heat_atlas_brandenburg` |
| *heat_low_decentral-demand_ind* | wfbb_heat_atlas_bb | demand_heat_region2 | Manueller Checkout von `features/#3_add_dataset_heat_atlas_brandenburg` |

BMWK lfs

*[Platzhalter für weiterführende Erklärungen]*

#### Wärmelastprofile Niedrigtemperaturwärme zentral und dezentral, normiert

| Name                                  | Raw-Datensatz | Dataset            | Kommentar |
|---------------------------------------|---------------|--------------------|-----------|
| *heat_low_central-demand_hh_profile*  | demandregio   | demand_heat_region |           |
| *heat_low_central-demand_cts_profile* | demandregio   | demand_heat_region |           |
| *heat_low_central-demand_ind_profile* | demandregio   | demand_heat_region |           |

*[Platzhalter für weiterführende Erklärungen]*

 
--------------------------------------------------------------------------------------------------------------------
### Wärmebedarf Prozesswärme (>80 °C)
Annahmen:
- Temperaturniveau: >100°C (high)
- Anwendungsgebiet: Industrieöfen und Hochtemperaturverfahren, Prozessdampf- und Warmwasser für Industrieprozesse [Quelle:Langfristszenarien]

#### Gesamtwärmebedarf Prozesswärme (2045) in MWh

| Name                   | Raw-Datensatz | Dataset                                                       | Kommentar       |
|------------------------|---------------|---------------------------------------------------------------|-----------------|
| *heat_high-demand_ind* |               |                                                               |                 |

-Quelle: https://github.com/asandhaa/ElectricalAndHeatProfiles/blob/main/IAEE%20Conference%20Paper%20Anna%20Sandhaas.pdf
basiert auf AGEB-Anwendungsbilanzen, durchschnittlicher Anteil für Prozesswärme 0.8 und Raumwärme und Warmwasser 0.2 über alle Industriezweige

#### Wärmelastprofile Prozesswärme 

| Name                           | Raw-Datensatz | Dataset                   | Kommentar                                                                |
|--------------------------------|---------------|---------------------------|--------------------------------------------------------------------------|
| *heat_high-demand_ind_profile* |  industry_heat_profiles              | demand_heat_high_ind  | Annahme: Profile von WZ08 Glas und Keramik, da durchschnittliches Profil | 

-----------------------------------------------------------------------------------------------------------------------------------------------------------------
### EE-Technologien

#### Einpeisezeitreihen EE-Technologien, normiert

| Name                                            | Raw-Datensatz | Dataset | Kommentar |
|-------------------------------------------------|---------------|---------|-----------|
| _electricity-wind-profile_                      |               |         |           | 
| _electricity-pv_ground-profile_                 |               |         |           | 
| _electricity-pv_agri_vertical-profile_          |               |         |           | 
| _electricity-pv_agri_horizontal-profile_        |               |         |           |
| _electricity-pv_rooftop-profile_                |               |         |           | 
| _heat_low_decentral-solarthermal_plant-profile_ |               |         |           |


#### Ausbaupotentiale EE-Technologie in MW

##### Wind

- In Rüdersdorf keine Windpotenzialflächen nach Regionalplan 2024 (1. Entwurf)
- Spalte `stp_2024_vr` gibt Flächenpotenziale in km², Umrechnung in MW mittels
  24 MW/km² (s. `technology_data.json`)

| Name                                                       | Raw-Datensatz          | Dataset                         | Kommentar                                                      |
|------------------------------------------------------------|------------------------|---------------------------------|----------------------------------------------------------------|
| _electricity-wind-capacity_potential_                      | rpg_ols_regional_plan  | potentialarea_wind_region       | Keine Windpotenzialflächen nach Regionalplan 2024 (1. Entwurf) |
| _electricity-pv_ground-capacity_potential_                 | oei_agri_pv            | potentialarea_pv_ground_region2 | Bezeichnung in Daten: *soil_quality_low*                       |
| _electricity-pv_agri_vertical-capacity_potential_          | oei_agri_pv            | potentialarea_pv_ground_region2 | Bezeichnung in Daten: *soil_quality_medium*                    |
| _electricity-pv_agri_horizontal-capacity_potential_        | oei_agri_pv            | potentialarea_pv_ground_region2 | Bezeichnung in Daten: *permanent_crops*                        |
| _electricity-pv_rooftop-capacity_potential_                | wfbb_pv_roof_potential | potentialarea_pv_roof_region2   | Annahme: 78 % von potentialarea_pv_roof_region2                |
| _heat_low_decentral-solarthermal_plant-capacity_potential_ |                        |                                 | Annahme: 22 % von potentialarea_pv_roof_region2                |

----------------------------------------------------------------------------------------------------------------------
### Wärmepumpen 

#### Ausbaupotentiale für Wärmepumpen in MW

| Name                                                | Raw-Datensatz | Dataset | Kommentar                   |
|-----------------------------------------------------|---------------|---------|-----------------------------|
| _electricity-heatpump_central-capacity_potential_   | -             | -       | Annahme: keine Ausbaugrenze | 
| _electricity-heatpump_decentral-capacity_potential_ | -             | -       | Annahme: keine Ausbaugrenze | 
| _electricity-heatpump_heat_high-capacity_potential_ | -             | -       | Annahme: keine Ausbaugrenze | 



#### COP-Zeitreihen für Wärmepumpen, normiert

| Name                                     | Raw-Datensatz          | Dataset                          | Kommentar                                                                 |
|------------------------------------------|------------------------|----------------------------------|---------------------------------------------------------------------------|
| _electricity-heatpump_central-profile_   |                        | `datasets/heatpump_cop`          | Zeitreihe wird auch für electricity-heatpump_decentral-profile angenommen | 
| _electricity-heatpump_decentral-profile_ |                        | `datasets/heatpump_cop`          |                                                                           | 
| _electricity-heatpump_heat_high-profile_ | heatpump_cop_heat_high | `datasets/heatpump_cop_heat_high`| Hochtemperaturwärmepumpe                                                  |

----------------------------------------------------------------------------------------------------------------------
## Kraft-Wärme-Kopplung

#### Ausbaupotentiale für KWK in MW

| Name                                                | Raw-Datensatz | Dataset | Kommentar |
|-----------------------------------------------------|---------------|---------|-----------|
| _biomass_solid-bpchp_heat_low_decentral-capacity_potential_   |               |         |           | 
| _biomass_solid-bpchp_heat_low_central-capacity_potential_ |               |         |           | 
| _biomass_solid-bpchp_heat_med-capacity_potential_  |               |         |           | 
| _biomass_gas-bpchp_heat_low_decentral-capacity_potential_  |               |         |           | 
| _biomass_gas-bpchp_heat_low_central-capacity_potential_  |               |         |           | 
| _biomass_gas-bpchp_heat_med-capacity_potential_  |               |         |           | 
| _biomass_solid-extchp_heat_med-capacity_potential_  |               |         |           | 
| _biomass_solid-extchp_heat_high-capacity_potential_  |               |         |           | 
| _biomass_gas-extchp_heat_med-capacity_potential_  |               |         |           | 
| _biomass_gas-extchp_heat_high-capacity_potential_  |               |         |           | 
| _h2-extchp_heat_high-capacity_potential_  |               |         |           | 
| _residual_waste-extchp_heat_high-capacity_potential_  |               |         |           | 
| _residual_waste-extchp_heat_low_central-capacity_potential_  |               |         |           | 


#### feste Kapazitäten für KWK in MW

| Name                                                | Raw-Datensatz | Dataset | Kommentar |
|-----------------------------------------------------|---------------|---------|-----------|
| _biomass_solid-bpchp_heat_low_decentral-capacity_   |               |         |           | 
| _biomass_solid-bpchp_heat_low_central-capacity_ |               |         |           | 
| _biomass_solid-bpchp_heat_med-capacity_  |               |         |           | 
| _biomass_gas-bpchp_heat_low_decentral-capacity_  |               |         |           | 
| _biomass_gas-bpchp_heat_low_central-capacity_  |               |         |           | 
| _biomass_gas-bpchp_heat_med-capacity_  |               |         |           | 
| _biomass_solid-extchp_heat_med-capacity_  |               |         |           | 
| _biomass_solid-extchp_heat_high-capacity_  |               |         |           | 
| _biomass_gas-extchp_heat_med-capacity_  |               |         |           | 
| _biomass_gas-extchp_heat_high-capacity_  |               |         |           | 
| _h2-extchp_heat_high-capacity_  |               |         |           | 
| _residual_waste-extchp_heat_high-capacity_  |               |         |           | 
| _residual_waste-extchp_heat_low_central-capacity_  |               |         |           | 

----------------------------------------------------------------------------------------------------------------------
## Brennkessel

#### Ausbaupotentiale für Brennkessel in MW

| Name                                               | Raw-Datensatz | Dataset | Kommentar                   |
|----------------------------------------------------|---------------|---------|-----------------------------|
| _biomass_solid-boiler_heat_med-capacity_potential_ |               |         | Annahme: keine Ausbaugrenze | 
| _biomass_solid-boiler_heat_high-capacity_potential_         |               |         |  Annahme: keine Ausbaugrenze                          | 
| _electricity-boiler_heat_med-capacity_potential_            |               |         |  Annahme: keine Ausbaugrenze                          | 
| _electricity-boiler_heat_high-capacity_potential_           |               |         |  Annahme: keine Ausbaugrenze                          | 
| _h2-boiler_heat_med-capacit_potential_                      |               |         |  Annahme: keine Ausbaugrenze                          | 
| _h2-boiler_heat_high-capacity_potential_                    |               |         |  Annahme: keine Ausbaugrenze                          | 


#### feste Kapazitäten für Brennkessel in MW

| Name                                          | Raw-Datensatz | Dataset | Kommentar |
|-----------------------------------------------|---------------|---------|-----------|
| _biomass_solid-boiler_heat_med-capacity_      |               |         | Annahme: keine Ausbaugrenze          | 
| _biomass_solid-boiler_heat_high-capacity_     |               |         | Annahme: keine Ausbaugrenze          | 
| _electricity-boiler_heat_med-capacity_        |               |         | Annahme: keine Ausbaugrenze          | 
| _electricity-boiler_heat_high-capacity_       |               |         | Annahme: keine Ausbaugrenze          | 
| _h2-boiler_heat_med-capacity_ |               |         | Annahme: keine Ausbaugrenze          | 
| _h2-boiler_heat_high-capacity_         |               |         | Annahme: keine Ausbaugrenze          | 

----------------------------------------------------------------------------------------------------------------------
## Speicher

#### Ausbaupotentiale für Speicher in MW

| Name                                                   | Raw-Datensatz | Dataset | Kommentar |
|--------------------------------------------------------|---------------|---------|-----------|
| _electricity-large_battery_storage-capacity_potential_ |               |         |   Annahme: keine Ausbaugrenze        |
| _heat_low_central-storage-capacity_potential_          |               |         |   Annahme: keine Ausbaugrenze        | 
| _heat_low_decentral-storage-capacity_  |               |         |                                            |


#### feste Kapazitäten für Speicher in MWh

| Name                                   | Raw-Datensatz | Dataset | Kommentar                                  |
|----------------------------------------|---------------|---------|--------------------------------------------|
| _electricity-small_battery_storage-capacity_ |               |         | Annahme: storage_capacity=pv_roof_capacity |

----------------------------------------------------------------------------------------------------------------------
## Elektrolyseur

#### Ausbaupotentiale für Elektrolyseur  in MW

| Name                                          | Raw-Datensatz | Dataset | Kommentar                   |
|-----------------------------------------------|---------------|---------|-----------------------------|
| _electricity-electrolyzer-capacity_potential_ |               |         | Annahme: keine Ausbaugrenze | 


#### feste Kapazitäten für Elektrolyseur in MW


| Name                                | Raw-Datensatz | Dataset | Kommentar |
|-------------------------------------|---------------|---------|-----------|
| _electricity-electrolyzer-capacity_ |               |         | -         |

----------------------------------------------------------------------------------------------------------------------
## Export

#### feste Kapazitäten für Export in MW

| Name                          | Raw-Datensatz | Dataset | Kommentar |
|-------------------------------|---------------|---------|-----------|
| _h2-export-capacity_          |               |         | -         |
| _electricity-export-capacity_ |               |         | -         |

----------------------------------------------------------------------------------------------------------------------
## Commodities

#### feste Kapazitäten für Commodities in MW

| Name                              | Raw-Datensatz | Dataset                                                                                                                                                      | Kommentar                                                                                                                                                                                                            |
|-----------------------------------|---------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| _biomass_gas-commodity-capacity_  |               |                                                                                                                                                              | 3% des Gesamtenergiebedarfs in https://www.agora-energiewende.de/fileadmin/Projekte/2023/2023-30_DE_KNDE_Update/A-EW_344_Klimaneutrales_Deutschland_WEB.pdf; summed_energy_demand*3%/0.5(Wirkungsgrad bpchp average) |
| _biomass_solid-commodity-capacity_ |               |  | 13% des Gesamtenergiebedarfs in https://www.agora-energiewende.de/fileadmin/Projekte/2023/2023-30_DE_KNDE_Update/A-EW_344_Klimaneutrales_Deutschland_WEB.pdf;summed_energy_demand*3%/0.5(Wirkungsgrad bpchp average) |
| _h2-commodity-capacity_ |               |                                                                                                                                                              | Rüdersdorf: 550 MW genehmigte Anschlussleistung an Wasserstoffkernnetz                                                                                                                                               |
| _residual_waste-commodity-capacity_ |               |                                                                                                                                                              | Rüdersdorf:  132MW;  installierte Leistung elektrisch der Müllverbrennungsanalge, die 2025 in Betrieb ist in  MW / Wirkungsgrad Müllverbrennungsanlage                                                               |
