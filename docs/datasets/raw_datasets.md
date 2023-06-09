# 'Raw' Datasets 

------------------------------
## Energiedaten Sachsen-Anhalt

Datensätze zur Energie- und Wasserversorgung des Statistischen Landesamtes
Sachsen-Anhalt.

### Daten

Stromverbrauch der Industriebetriebe nach Kreisen 2003-2021 in MWh
- [Quelle](https://statistik.sachsen-anhalt.de/themen/wirtschaftsbereiche/energie-und-wasserversorgung/tabellen-energieverwendung#c206986)

**Dataset: `raw/stala_st_energy`**

??? metadata "Metadata"
    ```json
    {
        "Datenquellen": {
            "Stromverbrauch der Industriebetriebe nach Kreisen": "https://statistik.sachsen-anhalt.de/themen/wirtschaftsbereiche/energie-und-wasserversorgung/tabellen-energieverwendung#c206986"
        }
    }
    ```

------------------------------
## OpenStreetMap

OpenStreetMap data extract for Sachsen-Anhalt.

**Dataset: `raw/osm_sachsen-anhalt`**

??? metadata "Metadata"
    ```json
    {
        "name": "openstreetmap",
        "title": "",
        "id": "openstreetmap",
        "description": "OpenStreetMap extract for federal state of Sachsen-Anhalt",
        "language": [
            "de-DE",
            "en-GB"
        ],
        "subject": [],
        "keywords": [
            "openstreetmap",
            "osm"
        ],
        "publicationDate": "2022-10-09",
        "context": {
            "homepage": "https://abw.rl-institut.de",
            "documentation": "https://digiplan.readthedocs.io",
            "sourceCode": "https://github.com/rl-institut/digipipe/",
            "contact": "https://reiner-lemoine-institut.de/ueber-uns/kontakt/",
            "grantNo": "None",
            "fundingAgency": "https://www.region-gestalten.bund.de",
            "fundingAgencyLogo": "https://www.region-gestalten.bund.de/Region/SiteGlobals/Frontend/Images/logo.svg",
            "publisherLogo": "https://reiner-lemoine-institut.de//wp-content/uploads/2015/09/rlilogo.png"
        },
        "spatial": {
            "location": "Sachsen-Anhalt",
            "extent": "Sachsen-Anhalt",
            "resolution": ""
        },
        "temporal": {
            "referenceDate": "2022-10-09",
            "timeseries": []
        },
        "sources": [
            {
                "title": "OpenStreetMap Data Extracts (Geofabrik)",
                "description": "Full data extract of OpenStreetMap data",
                "path": "https://download.geofabrik.de/europe/germany/sachsen-anhalt-221003.osm.pbf",
                "licenses": [
                    {
                        "name": "ODbL-1.0",
                        "title": "Open Data Commons Open Database License 1.0",
                        "path": "https://opendatacommons.org/licenses/odbl/1.0/",
                        "instruction": "You are free: To Share, To Create, To Adapt; As long as you: Attribute, Share-Alike, Keep open!",
                        "attribution": "\u00a9 OpenStreetMap contributors"
                    }
                ]
            }
        ],
        "licenses": [
            {
                "name": "ODbL-1.0",
                "title": "Open Data Commons Open Database License 1.0",
                "path": "https://opendatacommons.org/licenses/odbl/1.0/",
                "instruction": "You are free: To Share, To Create, To Adapt; As long as you: Attribute, Share-Alike, Keep open!",
                "attribution": "\u00a9 OpenStreetMap contributors"
            }
        ],
        "contributors": [
            {
                "title": "nesnoj",
                "email": "jonathan.amme@rl-institut.de",
                "date": "2022-10-09",
                "object": "metadata",
                "comment": "Create metadata"
            },
            {
                "title": "hedwig-lieselotte",
                "email": "hedwig.bartels@rl-institut.de",
                "date": "2023-03-28",
                "object": "metadata",
                "comment": "Check and update metadata"
            }
        ],
        "resources": [
            {
                "profile": "tabular-data-resource",
                "name": "model_draft.openstreetmap",
                "path": "",
                "format": "csv",
                "encoding": "UTF-8",
                "schema": {
                    "fields": [],
                    "primaryKey": [],
                    "foreignKeys": []
                },
                "dialect": {
                    "delimiter": "",
                    "decimalSeparator": ""
                }
            }
        ],
        "@id": null,
        "@context": "https://raw.githubusercontent.com/OpenEnergyPlatform/oemetadata/develop/metadata/latest/context.json",
        "review": {
            "path": "",
            "badge": ""
        },
        "metaMetadata": {
            "metadataVersion": "oemetadata_v1.5.1",
            "metadataLicense": {
                "name": "CC0-1.0",
                "title": "Creative Commons Zero v1.0 Universal",
                "path": "https://creativecommons.org/publicdomain/zero/1.0/"
            }
        },
        "_comment": {
            "metadata": "Metadata documentation and explanation (https://github.com/OpenEnergyPlatform/oemetadata)",
            "dates": "Dates and time must follow the ISO8601 including time zone (YYYY-MM-DD or YYYY-MM-DDThh:mm:ss\u00b1hh)",
            "units": "Use a space between numbers and units (100 m)",
            "languages": "Languages must follow the IETF (BCP47) format (en-GB, en-US, de-DE)",
            "licenses": "License name must follow the SPDX License List (https://spdx.org/licenses/)",
            "review": "Following the OEP Data Review (https://github.com/OpenEnergyPlatform/data-preprocessing/blob/master/data-review/manual/review_manual.md)",
            "null": "If not applicable use: null",
            "todo": "If a value is not yet available, use: todo"
        }
    }
    ```

------------------------------
## sEEnergies Pan-European Thermal Atlas 5.2 (Peta5)

Wärmebedarf für Europa 2015 in GJ (1ha Auflösung) für
- Haushalte: Raumwärme und Warmwasser
- GHD: Raumwärme, Warmwasser und Prozesswärme

Die Daten können auf der
[Projektseite](https://s-eenergies-open-data-euf.hub.arcgis.com)
eingesehen werden.

### Haushalte

Abgerufen mittels

```commandline
wget -O Peta5_0_1_HD_res.zip https://arcgis.com/sharing/rest/content/items/d7d18b63250240a49eb81db972aa573e/data
```

### GHD und Industrie

Abgerufen mittels

```commandline
wget -O Peta5_0_1_HD_ser.zip https://arcgis.com/sharing/rest/content/items/52ff5e02111142459ed5c2fe3d80b3a0/data
```

**Dataset: `raw/seenergies_peta5`**

??? metadata "Metadata"
    ```json
    {
        "Quellen": {
            "project": "https://www.seenergies.eu/peta5/",
            "data residential sector": "https://s-eenergies-open-data-euf.hub.arcgis.com/maps/d7d18b63250240a49eb81db972aa573e/about",
            "data service sector": "https://s-eenergies-open-data-euf.hub.arcgis.com/maps/52ff5e02111142459ed5c2fe3d80b3a0/about"
        }
    }
    ```

------------------------------
## BMWK Langfristszenarien

Langfristszenarien des Bundesministerium für Wirtschaft und Klimaschutz, Daten
auf Deutschlandebene.

Die Daten wurden über den
[Szenario Explorer](https://langfristszenarien.de/enertile-explorer-de/szenario-explorer/)
abgerufen.

### Verwendete Szenarien

- **T45-Strom:** Stromfokussiertes Szenario aus den T45-Szenarien aus 2023, die
  Wege zur Treibhausgasneutralität bis 2045 unter Einhaltung aktueller
  politischer Vorgaben erreichen. Die Daten dieses Szenarios werden als
  Grundlage für das Zielszenario in der Region verwendet.
- **TN-Strom:** Stromfokussiertes Szenario aus den TN-Szenarien aus 2021, die
  unterschiedliche Pfade für Deutschland mit dem Ziel treibhausgasneutral bis
  2050 zu werden. Die Daten dieses Szenarios werden als Grundlage für den
  Status quo verwendet.

### Daten

#### T45-Strom

| Datensatz                                      | Quelle                                                                                                    | Datei                                                     |
|------------------------------------------------|-----------------------------------------------------------------------------------------------------------|-----------------------------------------------------------|
| Gebäude: Haushalte und GHD Energiebedarf       | [Link](https://enertile-explorer.isi.fraunhofer.de:8443/open-view/51944/21559a9532131c061668bf0751e519e3) | `T45-Strom_buildings_heating_demand_by_carrier.csv`       |
| Gebäude: Anzahl der Heizungen nach Technologie | [Link](https://enertile-explorer.isi.fraunhofer.de:8443/open-view/51944/21559a9532131c061668bf0751e519e3) | `T45-Strom_buildings_heating_structure_by_technology.csv` |
| GHD Energieträger                              | [Link](https://enertile-explorer.isi.fraunhofer.de:8443/open-view/52700/c6980ea467bb26a922d34617b4fd4798) | `T45-Strom_cts_demand.csv`                                |
| Haushalte Energieträger                        | [Link](https://enertile-explorer.isi.fraunhofer.de:8443/open-view/52700/c6980ea467bb26a922d34617b4fd4798) | `T45-Strom_hh_demand.csv`                                 |
| Industrie Energiebedarf                        | [Link](https://enertile-explorer.isi.fraunhofer.de:8443/open-view/52612/9de48084ac2d54c418daaf02a6ee26e0) | `T45-Strom_ind_demand.csv`                                |

#### TN-Strom

| Datensatz                                      | Quelle                                                                                                    | Datei                                                    |
|------------------------------------------------|-----------------------------------------------------------------------------------------------------------|----------------------------------------------------------|
| Gebäude: Haushalte und GHD Energiebedarf       | [Link](https://enertile-explorer.isi.fraunhofer.de:8443/open-view/8198/698cee83d667a2f44fdea7e78ee799a2)  | `TN-Strom_buildings_heating_demand_by_carrier.csv`       |
| Gebäude: Anzahl der Heizungen nach Technologie | [Link](https://enertile-explorer.isi.fraunhofer.de:8443/open-view/8198/698cee83d667a2f44fdea7e78ee799a2)  | `TN-Strom_buildings_heating_structure_by_technology.csv` |
| GHD Energieträger                              | [Link](https://enertile-explorer.isi.fraunhofer.de:8443/open-view/8660/ae5a14ff0c320cbd31c5eeff2ede54ba)  | `TN-Strom_cts_demand.csv`                                |
| Haushalte Energieträger                        | [Link](https://enertile-explorer.isi.fraunhofer.de:8443/open-view/8660/ae5a14ff0c320cbd31c5eeff2ede54ba)  | `TN-Strom_hh_demand.csv`                                 |
| Industrie Energiebedarf                        | [Link](https://enertile-explorer.isi.fraunhofer.de:8443/open-view/29085/084bd7f45f40d31fd53341e6a94f532c) | `TN-Strom_ind_demand.csv`                                |

**Dataset: `raw/bmwk_long_term_scenarios`**

??? metadata "Metadata"
    ```json
    {
        "Datenquellen": {
            "Hauptseite": "https://langfristszenarien.de/enertile-explorer-de/szenario-explorer/"
        }
    }
    ```

------------------------------
## AGEB – Anwendungsbilanzen für die Endenergiesektoren 2011 bis 2021

Detaillierte Anwendungsbilanzen der Endenergiesektoren für 2020 und 2021 sowie
zusammenfassende Zeitreihen zum Endenergieverbrauch nach Energieträgern und
Anwendungszwecken für Jahre von 2011 bis 2021 der AG Energiebilanzen.

**Dataset: `raw/ageb_energy_balance`**

??? metadata "Metadata"
    ```json
    {
        "Quellen": {
            "Website": "https://ag-energiebilanzen.de/daten-und-fakten/anwendungsbilanzen/",
            "File": "https://ag-energiebilanzen.de/wp-content/uploads/2023/01/AGEB_21p2_V3_20221222.pdf"
        }
    }
    ```

------------------------------
## DemandRegio

Regionalisierte Bevölkerungsprognose, Haushalte sowie Strom- und Gasbedarfe
inkl. Zeitreihen auf Landkreisebene.

Die Daten wurden abgerufen mit einer
[modifizierten Version des DemandRegio disaggregators](https://github.com/nesnoj/disaggregator),
in der softwareseitige, jedoch keine methodischen Änderungen vorgenommen wurden.

Der disaggregator basiert auf Daten bis 2017, anschließende Jahre werden
fortgeschrieben.

Weitere Informationen zum Projekt DemandRegio
- [Abschlussbericht](https://www.ffe.de/wp-content/uploads/2020/10/DemandRegio_Abschlussbericht.pdf)
- [Abschlussworkshop](https://www.tu.berlin/er/forschung/projekte/demandregio-2)

Die erzeugten Rohdaten wie unten beschrieben wurden mittels
[API](http://opendata.ffe.de:4000/) abgerufen. Diese können alternativ direkt
vom [OpenData-Portal der FfE](https://opendata.ffe.de/project/demandregio/)
bezogen werden.

Verwendetes Wetterjahr für Gasbedarfszeitreihen: 2011

**Installation (in separater venv):**

```commandline
pip install disaggregator@git+https://github.com/nesnoj/disaggregator.git#egg=disaggregator
```

### Details zum Datenabruf

#### Bevölkerung

Bevölkerung (Summe) und Bevölkerung je Haushaltsgröße (1, 2, 3, 4, 5, >5) je
NUTS3.

Jahre
- Bevölkerung bis 2017 historische Werte
- Bevölkerung ab 2018 prognostizierte Werte basierend auf der 14. koordinierten
  Bevölkerungsvorausberechnung der Statistischen Ämter von Bund und Ländern.
- Haushalte nur 2011

```python
import pandas as pd
from disaggregator import data

## Population
dr_hh_population = pd.DataFrame()
for year in [2010, 2015, 2017, 2020, 2021, 2022, 2025, 2030, 2035, 2040, 2045]:
    dr_hh_population[year] = round(data.population(year=year)).astype(int)

dr_hh_population.to_csv("dr_hh_population.csv")

## Households
data.households_per_size().to_csv("dr_hh_households_2011.csv")
```

#### Haushalte: Strom

Bedarfe und SLP-Zeitreihen je NUTS3 mit Bottom-Up-Methode nach Haushaltsgröße.

Jahre
- 2017: Letzte verfügbare Daten
- 2022: Status quo, Fortschreibung mit Berücksichtigung Demografie und
  Wanderung
- 2035: Fortschreibungsjahr mit Berücksichtigung Demografie und Wanderung
- 2045: Fortschreibungsjahr

```python
from disaggregator import spatial, temporal

## Consumption
spatial.disagg_households_power(
    by="households",
    weight_by_income=True,
    year=2022,
    scale_by_pop=True,
).to_csv(f"dr_hh_power_demand_2022.csv")

## Timeseries
temporal.disagg_temporal_power_housholds_slp(
    use_nuts3code=True,
    by="households",
    weight_by_income=True,
    year=2022,
    scale_by_pop=True,
).to_csv(f"dr_hh_power_timeseries_2022.csv")
```

#### Haushalte: Gas

Zeitreihen je NUTS3

```python
from disaggregator import temporal

## Timeseries
temporal.disagg_temporal_gas_households(
    use_nuts3code=True,
    how='top-down',
    year=2011,
).to_csv(f"dr_hh_gas_timeseries_2011.csv")
```


#### GHD und Industrie: Strom

Bedarfe und Zeitreihen je NUTS3
- Bedarfe: Je Wirtschaftszweig (WZ), abzüglich Eigenerzeugung
- Zeitreihen: Für alle WZ bedarfsgewichtet aggregiert, Einzelprofile basieren
  je nach WZ auf gemessenen oder SLP inkl. Wanderung
- Letzte verfügbare Daten aus 2017, Fortschreibung für 2022 mit
  Berücksichtigung Beschäftigte und Effizienzgewinne

```python
from disaggregator import spatial, temporal

########
## CTS #
########

## Consumption
spatial.disagg_CTS_industry(
    sector='CTS',
    source='power',
    use_nuts3code=True,
    year=2022,
).to_csv("dr_cts_power_demand_2022.csv")
## Timeseries
temporal.disagg_temporal_power_CTS(
    detailed=False,
    use_nuts3code=True,
    year=2022,
).to_csv("dr_cts_power_timeseries_2022.csv")

#############
## Industry #
#############

## Consumption
spatial.disagg_CTS_industry(
    sector='industry',
    source='power',
    use_nuts3code=True,
    year=2022,
).to_csv("dr_ind_power_demand_2022.csv")
## Timeseries
temporal.disagg_temporal_industry(
    source="power",
    detailed=False,
    use_nuts3code=True,
    no_self_gen=False,
    year=2022,
).to_csv("dr_ind_power_timeseries_2022.csv")
```

#### GHD: Gas

Zeitreihen je NUTS3 für alle WZ bedarfsgewichtet aggregiert, Einzelprofile
basieren je nach WZ auf gemessenen oder SLP inkl. Wanderung. Letzte verfügbare
Daten aus 2017, Fortschreibung für 2022 mit Berücksichtigung Beschäftigte und
Effizienzgewinne.

```python
from disaggregator import spatial, temporal

## Timeseries
x=temporal.disagg_temporal_gas_CTS(
    detailed=False,
    use_nuts3code=True,
    year=2011,
).to_csv("dr_cts_gas_timeseries_2011.csv")
```

#### Industrie: Gas

Bedarfe und Zeitreihen je NUTS3
- Bedarfe: Je Wirtschaftszweig (WZ), abzüglich Eigenerzeugung
- Zeitreihen: Für alle WZ bedarfsgewichtet aggregiert, Einzelprofile basieren
  je nach WZ auf gemessenen oder SLP inkl. Wanderung
- Letzte verfügbare Daten aus 2017, Fortschreibung für 2022 mit
  Berücksichtigung Beschäftigte und Effizienzgewinne

```python
from disaggregator import spatial, temporal

## Consumption
spatial.disagg_CTS_industry(
    sector='industry',
    source='gas',
    use_nuts3code=True,
    year=2022,
).to_csv("dr_ind_gas_demand_2022.csv")
## Timeseries
x=temporal.disagg_temporal_industry(
    source="gas",
    detailed=False,
    use_nuts3code=True,
    no_self_gen=False,
    year=2011,
).to_csv("dr_ind_gas_timeseries_2011.csv")
```

**Dataset: `raw/demandregio`**

??? metadata "Metadata"
    ```json
    {
        "Quellen": {
            "Originales Tool": "https://github.com/DemandRegioTeam/disaggregator/",
            "Modifiziertes Tool": "https://github.com/nesnoj/disaggregator/"
        }
    }
    ```

------------------------------
## Sozialversicherungspflichtig Beschäftigte und Betriebe

Gemeindedaten der sozialversicherungspflichtig Beschäftigten am 30.06.2022 nach
Wohn- und Arbeitsort - Deutschland, Länder, Kreise und Gemeinden (Jahreszahlen)
der Bundesagentur für Arbeit.

**Dataset: `raw/ba_employment`**

??? metadata "Metadata"
    ```json
    {
        "Quellen": {
            "Website": "https://statistik.arbeitsagentur.de/SiteGlobals/Forms/Suche/Einzelheftsuche_Formular.html?nn=15024&topic_f=beschaeftigung-sozbe-gemband",
            "File": "https://statistik.arbeitsagentur.de/Statistikdaten/Detail/202206/iiia6/beschaeftigung-sozbe-gemband/gemband-dlk-0-202206-zip.zip?__blob=publicationFile&v=2}"
        }
    }
    ```

------------------------------
## Power units from Marktstammdatenregister

Power units from Marktstammdatenregister obtained and dumped with the
[open-mastr](https://github.com/OpenEnergyPlatform/open-MaStR) tool. The
dump was created with:

```
from open_mastr import Mastr
db = Mastr()
db.download("bulk")
db.to_csv(None)  # (None for all data)
```

The dumped CSV files (all tables) have been extended by a custom export of
storage units with
`sqlite3 -header -csv -separator "," open-mastr.db "select * from storage_units;" > bnetza_mastr_storage_unit_raw.csv`.
Subsequently, all files were zipped.

The Marktstammdatenregister (MaStR) is a German register provided by the
German Federal Network Agency (Bundesnetzagentur / BNetza) that keeps
track of all power and gas units located in Germany.

**Dataset: `raw/bnetza_mastr`**

??? metadata "Metadata"
    ```json
    {
        "name": "bnetza_mastr",
        "title": "Marktstammdatenregisterdaten",
        "id": "bnetza_mastr",
        "description": "Daten aus dem Marktstammdatenregister der Bundesnetzagentur",
        "language": [
            "en-GB",
            "de-DE"
        ],
        "subject": null,
        "keywords": [
            "Markstammdatenregister",
            "openmastr",
            "mastr"
        ],
        "publicationDate": "2022-12-19",
        "context": {
            "homepage": "https://abw.rl-institut.de",
            "documentation": "https://digiplan.readthedocs.io",
            "sourceCode": "https://github.com/rl-institut/digipipe/",
            "contact": "https://reiner-lemoine-institut.de/ueber-uns/kontakt/",
            "grantNo": "None",
            "fundingAgency": "https://www.region-gestalten.bund.de",
            "fundingAgencyLogo": "https://www.region-gestalten.bund.de/Region/SiteGlobals/Frontend/Images/logo.svg",
            "publisherLogo": "https://reiner-lemoine-institut.de//wp-content/uploads/2015/09/rlilogo.png"
        },
        "spatial": {
            "location": "Germany",
            "extent": "Germany",
            "resolution": null
        },
        "temporal": {
            "referenceDate": "2022-12-19",
            "timeseries": null
        },
        "sources": [
            {
                "title": "Marktstammdatenregister",
                "description": "Marktstammdatenregister der Bundesnetzagentur Deutschland",
                "path": "https://www.marktstammdatenregister.de/MaStR/Datendownload",
                "licenses": [
                    {
                        "name": "DL-DE-BY-2.0",
                        "title": "Open Data Datenlizenz Deutschland \u2013 Namensnennung \u2013 Version 2.0",
                        "path": "http://www.govdata.de/dl-de/by-2-0",
                        "instruction": "The data and meta-data provided may, for commercial and non-commercial use, in particular be copied, printed, presented, altered, processed and transmitted to third parties; be merged with own data and with the data of others and be combined to form new and independent datasets;be integrated in internal and external business processes, products and applications in public and non-public electronic networks.",
                        "attribution": "\u00a9 Marktstammdatenregister 2023"
                    }
                ]
            }
        ],
        "licenses": [
            {
                "name": "DL-DE-BY-2.0",
                "title": "Open Data Datenlizenz Deutschland \u2013 Namensnennung \u2013 Version 2.0?",
                "path": "http://www.govdata.de/dl-de/by-2-0",
                "instruction": "The data and meta-data provided may, for commercial and non-commercial use, in particular be copied, printed, presented, altered, processed and transmitted to third parties; be merged with own data and with the data of others and be combined to form new and independent datasets;be integrated in internal and external business processes, products and applications in public and non-public electronic networks.",
                "attribution": "\u00a9 Marktstammdatenregister 2023"
            }
        ],
        "contributors": [
            {
                "title": "hedwiglieselotte",
                "email": "hedwig.bartels@rl-institut.de",
                "date": "2023-03-28",
                "object": "metadata",
                "comment": "Create metadata"
            }
        ],
        "resources": [
            {
                "profile": null,
                "name": null,
                "path": null,
                "format": null,
                "encoding": null,
                "schema": {
                    "fields": [],
                    "primaryKey": [],
                    "foreignKeys": []
                },
                "dialect": {
                    "delimiter": "",
                    "decimalSeparator": "."
                }
            }
        ],
        "@id": null,
        "@context": "https://raw.githubusercontent.com/OpenEnergyPlatform/oemetadata/develop/metadata/latest/context.json",
        "review": {
            "path": "",
            "badge": ""
        },
        "metaMetadata": {
            "metadataVersion": "OEP-1.5.2",
            "metadataLicense": {
                "name": "CC0-1.0",
                "title": "Creative Commons Zero v1.0 Universal",
                "path": "https://creativecommons.org/publicdomain/zero/1.0/"
            }
        },
        "_comment": {
            "metadata": "Metadata documentation and explanation (https://github.com/OpenEnergyPlatform/oemetadata)",
            "dates": "Dates and time must follow the ISO8601 including time zone (YYYY-MM-DD or YYYY-MM-DDThh:mm:ss\u00b1hh)",
            "units": "Use a space between numbers and units (100 m)",
            "languages": "Languages must follow the IETF (BCP47) format (en-GB, en-US, de-DE)",
            "licenses": "License name must follow the SPDX License List (https://spdx.org/licenses/)",
            "review": "Following the OEP Data Review (https://github.com/OpenEnergyPlatform/data-preprocessing/blob/master/data-review/manual/review_manual.md)",
            "null": "If not applicable use: null",
            "todo": "If a value is not yet available, use: todo"
        }
    }
    ```

------------------------------
## Temperatur

Stündliche Mittelwerte der Luft- und Erdbodentemperatur des Deutschen
Wetterdienstes
([Climate Data Center](https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/hourly/))
für das Jahr 2011 je Gemeinde in der Region ABW, vorverarbeitet im Projekt
[WindNODE](https://windnode-abw.readthedocs.io/en/latest/energy_system_model.html#energy-demand-today).

Werte
- `temp_amb`: Lufttemperatur in 2 m Höhe
- `temp_soil`: Erdbodentemperatur in 1 m Tiefe

Verwendete Stationen
- Wittenberg
- Köthen
- Jessnitz
- Seehausen
- Holzdorf

Die Zuordnung der Stationsmesswerte zu Gemeinden erfolgte über die jeweils
nächstgelegene Wetterstation.

**Dataset: `raw/dwd_temperature`**

??? metadata "Metadata"
    ```json
    {
        "Datenquellen": {
            "Open Data Bereich des Climate Data Center des DWD": "https://www.dwd.de/DE/leistungen/cdc/climate-data-center.html",
            "Datensatz Lufttemperatur": "https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/hourly/air_temperature/historical/",
            "Datensatz Bodentemperatur": "https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/hourly/soil_temperature/historical/"
        }
    }
    ```

------------------------------
## Regionalstatistik (GENESIS)

Enthält folgende Datensätze der statistischen Ämter des Bundes und der Länder:

### Energieverwendung der Betriebe im Verarbeitenden Gewerbe (43531-01-02-4)

Jahreserhebung ü. die Energieverwendung der Betriebe im verarbeitendem Gewerbe.

Der Datensatz umfasst:
- Betriebe des Verarbeitenden Gewerbes sowie des Bergbaus
und der Gewinnung von Steinen und Erden von Unternehmen des
Produzierenden Gewerbes mit im Allgemeinen 20 und mehr
Beschäftigten.
- Betriebe des Verarbeitenden Gewerbes sowie des Bergbaus
und  der Gewinnung von Steinen und Erden mit im Allgemeinen
20 und mehr Beschäftigten von Unternehmen der übrigen
Wirtschaftsbereiche.
Die Berichterstattung schließt Verarbeitende Betriebe des
Handwerks ein.
Bei 7 Wirtschaftszweigen gilt eine Abschneidegrenze von 10
Beschäftigten. Die Merkmalswerte beziehen sich auf den
gesamten Betrieb, schließen damit die nicht produzierenden
Betriebsteile mit ein.
Maßgebend für die Zuordnung ist ab 2008 die „Klassifikation
der Wirtschaftszweige, Ausgabe 2008 (WZ 2008)“, und zwar
die Abschnitte B und C.

- Datei: `43531-01-02-4.xlsx`
- Stand: 2021

### Betriebe, tätige Personen, Bruttoentgelte (42111-01-04-5)

Jahreserhebung ü. Betriebe, tätige Personen und Bruttoentgelte der Betriebe im
verarbeitendem Gewerbe.

Der Datensatz umfasst:
- Sämtliche Betriebe des Wirtschaftsbereiches Verarbeitendes
Gewerbe sowie Bergbau und Gewinnung von Steinen und Erden,
wenn diese Betriebe zu Unternehmen des Bereiches
Verarbeitendes Gewerbe sowie Bergbau und Gewinnung von
Steinen und Erden gehören und in diesen Unternehmen
mindestens 20 Personen tätig sind;
- die Betriebe des Wirtschaftsbereiches Verarbeitendes
Gewerbe sowie Bergbau und Gewinnung von Steinen und Erden
mit mindestens 20 tätigen Personen, sofern diese Betriebe
zu Unternehmen gehören, deren wirtschaftlicher Schwerpunkt
außerhalb des Bereiches Verarbeitendes Gewerbe sowie
Bergbau und Gewinnung von Steinen und Erden liegt.
Bei 7 kleinbetrieblich strukturierten Branchen gilt eine
untere Erfassungsgrenze von 10 tätigen Personen.
Die Auswahl erfolgt jeweils nach dem Beschäftigtenstand Ende
September des Vorjahres. Die ausgewiesene Beschäftigtenzahl
betrifft dagegen die von Ende September des Berichtsjahres.
Die Merkmalswerte beziehen sich auf den gesamten Betrieb,
schließen damit die nicht produzierenden Betriebsteile mit
ein.
Maßgebend für die Zuordnung ist ab 2009 die „Klassifikation
der Wirtschaftszweige, Ausgabe 2008 (WZ 2008)“, und zwar
die Abschnitte B und C.

- Datei: `42111-01-04-5.xlsx`
- Stand: 30.09.2021

### Gebäude mit Wohnraum nach Heizungsart (31211-04-01-5-B)

Zensus 2011: Gebäude mit Wohnraum nach Heizungsart

- Datei: `31211-04-01-5-B.xlsx`
- Stand: 09.05.2011

### Gebäude mit Wohnraum nach Heizungsart (31231-02-01-5)

Bestand an Wohngebäuden und Wohnungen in Wohn- und Nichtwohngebäuden -
Fortschreibung auf Basis der endgültigen Ergebnisse der Gebäude- und
Wohnungszählung 2011 (Zensus 2011).

- Datei: `31231-02-01-5.xlsx`
- Stand: 31.12.2021

**Dataset: `raw/regiostat`**

??? metadata "Metadata"
    ```json
    {
        "Datenquellen": {
            "43531-01-02-4": "https://www.regionalstatistik.de/genesis//online?operation=table&code=43531-01-02-4",
            "42111-01-04-5": "https://www.regionalstatistik.de/genesis//online?operation=table&code=42111-01-04-5"
        }
    }
    ```

------------------------------
## Lokale Verwaltungseinheiten

Lokale Verwaltungseinheiten (LAUs) von Eurostat, mit NUTS kompatibel. Diese LAU
sind die Bausteine der NUTS und umfassen die Gemeinden und Kommunen der
Europäischen Union.

**Dataset: `raw/eurostat_lau`**

??? metadata "Metadata"
    ```json
    {
        "Datenquellen": {
            "Hauptseite": "https://ec.europa.eu/eurostat/de/web/nuts/local-administrative-units",
            "Daten": "https://ec.europa.eu/eurostat/documents/345175/501971/EU-27-LAU-2022-NUTS-2021.xlsx"
        }
    }
    ```

------------------------------
## Bevölkerungsprognose Sachsen-Anhalt

Bevölkerungsprognose je Gemeinde bis 2035 des Statistischen Landesamtes
Sachsen-Anhalt. Stand: 2021

**Dataset: `raw/stala_st_pop_prog`**

??? metadata "Metadata"
    ```json
    {
        "name": "stala_st_pop_prog",
        "title": "Regionalisierte Bev\u00f6lkerungsprognose",
        "id": "stala_st_pop_prog",
        "description": "Prognostizierter Bev\u00f6lkerungsstand in den Gemeinden, kreisfreien St\u00e4dten und Landkreisen nach Prognosejahr und Geschlecht",
        "language": [
            "de-DE"
        ],
        "subject": [],
        "keywords": [
            "Bev\u00f6lkerungsprognose",
            "population"
        ],
        "publicationDate": null,
        "context": {
            "homepage": "https://abw.rl-institut.de",
            "documentation": "https://digiplan.readthedocs.io",
            "sourceCode": "https://github.com/rl-institut/digipipe/",
            "contact": "https://reiner-lemoine-institut.de/ueber-uns/kontakt/",
            "grantNo": "None",
            "fundingAgency": "https://www.region-gestalten.bund.de",
            "fundingAgencyLogo": "https://www.region-gestalten.bund.de/Region/SiteGlobals/Frontend/Images/logo.svg",
            "publisherLogo": "https://reiner-lemoine-institut.de//wp-content/uploads/2015/09/rlilogo.png"
        },
        "spatial": {
            "location": "Sachsen-Anhalt",
            "extent": "Sachsen-Anhalt",
            "resolution": ""
        },
        "temporal": {
            "referenceDate": null,
            "timeseries": [
                {
                    "start": "2019",
                    "end": "2035",
                    "resolution": "1 year",
                    "alignment": null,
                    "aggregationType": "sum"
                }
            ]
        },
        "sources": [
            {
                "title": "1_Internettabelle_7RBP_nach_Prognosejahr_Geschlecht_alle_Ebenen",
                "description": "Prognostizierter Bev\u00f6lkerungsstand in den Gemeinden, kreisfreien St\u00e4dten und Landkreisen nach Prognosejahr und Geschlecht",
                "path": "statistik.sachsen-anhalt.de/themen/bevoelkerung-mikrozensus-freiwillige-haushaltserhebungen/bevoelkerung/bevoelkerungsprognose-und-haushalteprognose/#c312231",
                "licenses": [
                    {
                        "name": "DL-DE-BY-2.0",
                        "title": "Data licence Germany \u2013 attribution \u2013 Version 2.0",
                        "path": "https://www.govdata.de/dl-de/by-2-0",
                        "instruction": "The data and meta-data provided may, for commercial and non-commercial use, in particular be copied, printed, presented, altered, processed and transmitted to third parties; be merged with own data and with the data of others and be combined to form new and independent datasets; be integrated in internal and external business processes, products and applications in public and non-public electronic networks.",
                        "attribution": "\u00a9 2023 Landesportal Sachsen-Anhalt "
                    }
                ]
            }
        ],
        "contributors": [
            {
                "title": "hedwiglieselotte",
                "email": "hedwig.bartels@rl-institut.de",
                "date": "2023-03-28",
                "object": "metadata",
                "comment": "create metadata"
            }
        ],
        "resources": [
            {
                "profile": "tabular-data-resource",
                "name": "1_Internettabelle_7RBP_nach_Prognosejahr_Geschlecht_alle_Ebenen",
                "path": "https://statistik.sachsen-anhalt.de/fileadmin/Bibliothek/Landesaemter/StaLa/startseite/Themen/Bevoelkerung/Tabellen/Bevoelkerungsprognose/1_Internettabelle_7RBP_nach_Prognosejahr_Geschlecht_alle_Ebenen.xlsx",
                "format": "xlxs",
                "encoding": "",
                "schema": {
                    "fields": [],
                    "primaryKey": [],
                    "foreignKeys": []
                },
                "dialect": {
                    "delimiter": "",
                    "decimalSeparator": "."
                }
            }
        ],
        "@id": [],
        "@context": "https://raw.githubusercontent.com/OpenEnergyPlatform/oemetadata/develop/metadata/latest/context.json",
        "review": {
            "path": "",
            "badge": ""
        },
        "metaMetadata": {
            "metadataVersion": "OEP-1.5.2",
            "metadataLicense": {
                "name": "CC0-1.0",
                "title": "Creative Commons Zero v1.0 Universal",
                "path": "https://creativecommons.org/publicdomain/zero/1.0/"
            }
        },
        "_comment": {
            "metadata": "Metadata documentation and explanation (https://github.com/OpenEnergyPlatform/oemetadata)",
            "dates": "Dates and time must follow the ISO8601 including time zone (YYYY-MM-DD or YYYY-MM-DDThh:mm:ss\u00b1hh)",
            "units": "Use a space between numbers and units (100 m)",
            "languages": "Languages must follow the IETF (BCP47) format (en-GB, en-US, de-DE)",
            "licenses": "License name must follow the SPDX License List (https://spdx.org/licenses/)",
            "review": "Following the OEP Data Review (https://github.com/OpenEnergyPlatform/data-preprocessing/blob/master/data-review/manual/review_manual.md)",
            "null": "If not applicable use: null",
            "todo": "If a value is not yet available, use: todo"
        }
    }
    ```

------------------------------
## Municipalities population

Municipalities with population from the Federal Statistical Office.

**Dataset: `raw/destatis_gv`**

??? metadata "Metadata"
    ```json
    {
        "name": "destatis_gv",
        "title": "Adminstratives Gemeinndeverzeichnis",
        "id": "destatis_gv",
        "description": "Alle politisch selbst\u00e4ndigen Gemeinden mit ausgew\u00e4hlten Merkmalen am 31.12.2022 ",
        "language": [
            "de-DE"
        ],
        "subject": null,
        "keywords": [
            "destatis",
            "gemeindeverzeichnis"
        ],
        "publicationDate": "2023-01-12",
        "context": {
            "homepage": "https://abw.rl-institut.de",
            "documentation": "https://digiplan.readthedocs.io",
            "sourceCode": "https://github.com/rl-institut/digipipe/",
            "contact": "https://reiner-lemoine-institut.de/ueber-uns/kontakt/",
            "grantNo": null,
            "fundingAgency": "https://www.region-gestalten.bund.de",
            "fundingAgencyLogo": "https://www.region-gestalten.bund.de/Region/SiteGlobals/Frontend/Images/logo.svg",
            "publisherLogo": "https://reiner-lemoine-institut.de//wp-content/uploads/2015/09/rlilogo.png"
        },
        "spatial": {
            "location": "Germany",
            "extent": "Germany",
            "resolution": ""
        },
        "temporal": {
            "referenceDate": "2022-02-14",
            "timeseries": null
        },
        "sources": [
            {
                "title": "Statistisches Bundesamt",
                "description": "Alle politisch selbst\u00e4ndigen Gemeineden mit ausgew\u00e4hlten Merkmalen am 31.12.2022 (4.Quartal)",
                "path": "https://www.destatis.de/DE/Themen/Laender-Regionen/Regionales/Gemeindeverzeichnis/Administrativ/Archiv/GVAuszugQ/AuszugGV4QAktuell.html",
                "licenses": [
                    {
                        "name": "DL-DE-BY-2.0",
                        "title": "Data licence Germany \u2013 attribution \u2013 version 2.0",
                        "path": "https://www.govdata.de/dl-de/by-2-0",
                        "instruction": "The data and meta-data provided may, for commercial and non-commercial use, in particular be copied, printed, presented, altered, processed and transmitted to third parties; be merged with own data and with the data of others and be combined to form new and independent datasets; be integrated in internal and external business processes, products and applications in public and non-public electronic networks.",
                        "attribution": "\u00a9 Statistisches Bundesamt (Destatis), 2023"
                    }
                ]
            }
        ],
        "licenses": [
            {
                "name": "DL-DE-BY-2.0",
                "title": "Data licence Germany \u2013 attribution \u2013 version 2.0",
                "path": "https://www.govdata.de/dl-de/by-2-0",
                "instruction": "The data and meta-data provided may, for commercial and non-commercial use, in particular be copied, printed, presented, altered, processed and transmitted to third parties; be merged with own data and with the data of others and be combined to form new and independent datasets; be integrated in internal and external business processes, products and applications in public and non-public electronic networks.",
                "attribution": "\u00a9 Statistisches Bundesamt (Destatis), 2023"
            }
        ],
        "contributors": [
            {
                "title": "hedwiglieselotte",
                "email": "hedwig.bartels@rl-institut.de",
                "date": "2023-03-28",
                "object": "metadata",
                "comment": "create metadata"
            }
        ],
        "resources": [
            {
                "profile": null,
                "name": null,
                "path": null,
                "format": null,
                "encoding": null,
                "schema": {
                    "fields": [],
                    "primaryKey": null,
                    "foreignKeys": []
                },
                "dialect": {
                    "delimiter": "",
                    "decimalSeparator": "."
                }
            }
        ],
        "@id": null,
        "@context": "https://raw.githubusercontent.com/OpenEnergyPlatform/oemetadata/develop/metadata/latest/context.json",
        "review": {
            "path": "",
            "badge": ""
        },
        "metaMetadata": {
            "metadataVersion": "OEP-1.5.2",
            "metadataLicense": {
                "name": "CC0-1.0",
                "title": "Creative Commons Zero v1.0 Universal",
                "path": "https://creativecommons.org/publicdomain/zero/1.0/"
            }
        },
        "_comment": {
            "metadata": "Metadata documentation and explanation (https://github.com/OpenEnergyPlatform/oemetadata)",
            "dates": "Dates and time must follow the ISO8601 including time zone (YYYY-MM-DD or YYYY-MM-DDThh:mm:ss\u00b1hh)",
            "units": "Use a space between numbers and units (100 m)",
            "languages": "Languages must follow the IETF (BCP47) format (en-GB, en-US, de-DE)",
            "licenses": "License name must follow the SPDX License List (https://spdx.org/licenses/)",
            "review": "Following the OEP Data Review (https://github.com/OpenEnergyPlatform/data-preprocessing/blob/master/data-review/manual/review_manual.md)",
            "null": "If not applicable use: null",
            "todo": "If a value is not yet available, use: todo"
        }
    }
    ```

------------------------------
## Administative areas of Germany

Administative areas of Germany (Verwaltungsgebiete 1:250 000).

**Dataset: `raw/bkg_vg250`**

??? metadata "Metadata"
    ```json
    {
        "name": "bkg_vg250",
        "title": "Adminstrative areas of Germany",
        "id": "bkg_vb250",
        "description": "Geopackage with administative areas of Germany - Verwaltungsgebiete 1:250 000",
        "language": [
            "en-GB",
            "de-DE"
        ],
        "subject": null,
        "keywords": [
            "adminstrative areas",
            "Verwaltungsgebiete"
        ],
        "publicationDate": "2022-01-01",
        "context": {
            "homepage": "https://abw.rl-institut.de",
            "documentation": "https://digiplan.readthedocs.io",
            "sourceCode": "https://github.com/rl-institut/digipipe/",
            "contact": "https://reiner-lemoine-institut.de/ueber-uns/kontakt/",
            "grantNo": "None",
            "fundingAgency": "https://www.region-gestalten.bund.de",
            "fundingAgencyLogo": "https://www.region-gestalten.bund.de/Region/SiteGlobals/Frontend/Images/logo.svg",
            "publisherLogo": "https://reiner-lemoine-institut.de//wp-content/uploads/2015/09/rlilogo.png"
        },
        "spatial": {
            "location": "Germany",
            "extent": "Germany",
            "resolution": "1:250 000"
        },
        "temporal": {
            "referenceDate": "2022-01-01",
            "timeseries": null
        },
        "sources": [
            {
                "title": "Bundesamt f\u00fcr Kartographie und Geod\u00e4sie - Verwaltungsgebiete 1:250 000 VG250 (Ebenen)",
                "description": "Dieser Datensatz stellt die Verwaltungsgebiete 1:250 000 (VG250) mit Stand 01.01. f\u00fcr das Gebiet der Bundesrepublik Deutschland bereit.",
                "path": "https://gdz.bkg.bund.de/index.php/default/digitale-geodaten/verwaltungsgebiete/verwaltungsgebiete-1-250-000-stand-01-01-vg250-01-01.html",
                "licenses": [
                    {
                        "name": "DL-DE-BY-2.0",
                        "title": "Open Data Datenlizenz Deutschland \u2013 Namensnennung \u2013 Version 2.0",
                        "path": "http://www.govdata.de/dl-de/by-2-0",
                        "instruction": "The data and meta-data provided may, for commercial and non-commercial use, in particular be copied, printed, presented, altered, processed and transmitted to third parties; be merged with own data and with the data of others and be combined to form new and independent datasets;be integrated in internal and external business processes, products and applications in public and non-public electronic networks.",
                        "attribution": " \u00a9 GeoBasis-DE / BKG - 2022"
                    }
                ]
            }
        ],
        "contributors": [
            {
                "title": "hedwiglieselotte",
                "email": "hedwig.bartels@rl-institut.de",
                "date": "2023-03-23",
                "object": "metadata",
                "comment": "create metadata"
            }
        ],
        "resources": [
            {
                "profile": null,
                "name": null,
                "path": null,
                "format": null,
                "encoding": null,
                "schema": {
                    "fields": [],
                    "primaryKey": [],
                    "foreignKeys": []
                },
                "dialect": {
                    "delimiter": "",
                    "decimalSeparator": "."
                }
            }
        ],
        "@id": null,
        "@context": "https://raw.githubusercontent.com/OpenEnergyPlatform/oemetadata/develop/metadata/latest/context.json",
        "review": {
            "path": "",
            "badge": ""
        },
        "metaMetadata": {
            "metadataVersion": "OEP-1.5.2",
            "metadataLicense": {
                "name": "CC0-1.0",
                "title": "Creative Commons Zero v1.0 Universal",
                "path": "https://creativecommons.org/publicdomain/zero/1.0/"
            }
        },
        "_comment": {
            "metadata": "Metadata documentation and explanation (https://github.com/OpenEnergyPlatform/oemetadata)",
            "dates": "Dates and time must follow the ISO8601 including time zone (YYYY-MM-DD or YYYY-MM-DDThh:mm:ss\u00b1hh)",
            "units": "Use a space between numbers and units (100 m)",
            "languages": "Languages must follow the IETF (BCP47) format (en-GB, en-US, de-DE)",
            "licenses": "License name must follow the SPDX License List (https://spdx.org/licenses/)",
            "review": "Following the OEP Data Review (https://github.com/OpenEnergyPlatform/data-preprocessing/blob/master/data-review/manual/review_manual.md)",
            "null": "If not applicable use: null",
            "todo": "If a value is not yet available, use: todo"
        }
    }
    ```
