# 'Raw' Datasets 

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
## Erzeugungsanlagen aus Marktstammdatenregister

Ereugungsanlagen aus dem Markstammdatenregister, das mit dem Tool
[open-mastr](https://github.com/OpenEnergyPlatform/open-MaStR) erstellt und
abgelegt wurde. Die Daten wurden folgendermaßen erstellt:
```
from open_mastr import Mastr
db = Mastr()
db.download("bulk")
db.to_csv(None)  # (None for all data)
```

Die abgelegten CSV-Dateien (alle Tabellen) wurden um einen benutzerdefinierten
Export von Speichereinheiten mit
`sqlite3 -header -csv -separator "," open-mastr.db "select * from storage_units;" > bnetza_mastr_storage_unit_raw.csv`
erweitert. Anschließend wurden alle Dateien komprimiert.

Das Marktstammdatenregister (MaStR) ist ein deutsches Register, welches von der
Bundesnetzagentur (BNetza) bereitgestellt wird und alle in Deutschland
befindlichen Strom- und Gasanlagen erfasst.

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

| Datensatz                                      | Quelle                                                                                                                    | Datei                                                     |
|------------------------------------------------|---------------------------------------------------------------------------------------------------------------------------|-----------------------------------------------------------|
| Gebäude: Haushalte und GHD Energiebedarf       | [Link](https://enertile-explorer.isi.fraunhofer.de:8443/open-view/51944/21559a9532131c061668bf0751e519e3)                 | `T45-Strom_buildings_heating_demand_by_carrier.csv`       |
| Gebäude: Anzahl der Heizungen nach Technologie | [Link](https://enertile-explorer.isi.fraunhofer.de:8443/open-view/51944/21559a9532131c061668bf0751e519e3)                 | `T45-Strom_buildings_heating_structure_by_technology.csv` |
| GHD Energieträger                              | [Link](https://enertile-explorer.isi.fraunhofer.de:8443/open-view/52700/c6980ea467bb26a922d34617b4fd4798)                 | `T45-Strom_cts_demand.csv`                                |
| Haushalte Energieträger                        | [Link](https://enertile-explorer.isi.fraunhofer.de:8443/open-view/52700/c6980ea467bb26a922d34617b4fd4798)                 | `T45-Strom_hh_demand.csv`                                 |
| Industrie Energiebedarf                        | [Link](https://enertile-explorer.isi.fraunhofer.de:8443/open-view/52612/9de48084ac2d54c418daaf02a6ee26e0)                 | `T45-Strom_ind_demand.csv`                                |
| Stromsystem Deutschland Leistung               | [Link](https://enertile-explorer.isi.fraunhofer.de:8443/open-view/48766/5c11999a03c547e04e73d61e4b5fc633)                 | `T45-Strom_electricity_installed_power.csv`               |

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
## Marktstammdatenregister Datenkorrektur PV

Überprüfung und manuelle Datenkorrektur der Photovoltaikanlagen aus dem
prozessierten Marktstammdatenregister (Datensatz:
[bnetza_mastr](../bnetza_mastr/dataset.md)).

### Plausibiltätsprüfung

Um grobe Fehler herauszufiltern wird überprüft, ob
- Anlage in Betrieb ist (status = "In Betrieb"),
- Anlage Strom produziert,
- Brutto- und Nettokapazität plausibel sind und
- die Kategorisierung, d.h. Zuordnung eine PV-Anlage zu Freifläche oder Dach,
  plausibel ist (manuelle, visuelle Prüfung von geolokalisierten
  PV-Aufdachanlagen anhand von
  [Orthofotos](https://www.geodatenportal.sachsen-anhalt.de/wss/service/ST_LVermGeo_DOP_WMS_OpenData/guest))

### Dateien

- Korrektur Freiflächenanlagen `bnetza_mastr_pv_ground_region_correction.ods`
- Korrektur Aufdachanlagen `bnetza_mastr_pv_roof_region_correction.ods`

mit Spalten
- _mastr_id_: ID aus dem MaStR
- _reason_: Fehler (wrong_type, wrong_position)
- _wrong_attr_: Fehlerhaftes Attribut
- _correction_: Korrigierter Attributwert (None, wenn Korrektur nicht möglich).
  Korrigierte Geometrien liegen in EPSG:3035 vor.

**Dataset: `raw/bnetza_mastr_correction_region`**

??? metadata "Metadata"
    ```json
    {
        "name": "bnetza_mastr_correction",
        "title": "Marktstammdatenregisterdaten - Manuelle Korrektur",
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
## Geodaten PV- und Windflächenrechner

Geodaten aus dem
[PV- und Windflächenrechner](https://www.agora-energiewende.de/service/pv-und-windflaechenrechner/).

Mehr Informationen:
- [Begleitdokument](https://zenodo.org/record/6794558)
- [Geodaten Potenzialflächen](https://zenodo.org/record/6728382)

Enthält
- Geodaten
- Metadaten
- App-Datapackage

**Dataset: `raw/rli_pv_wfr`**

??? metadata "Metadata"
    ```json
    {
        "Quellen": {
            "Begleitdokument": "https://zenodo.org/record/6794558",
            "Geodaten": "https://zenodo.org/record/6728382"
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
## Verwaltungsgebiete Deutschlands

Verwaltungsgebiete Deutschlands (Verwaltungsgebiete 1:250 000).

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

------------------------------
## Regionalplan Anhalt-Bitterfeld-Wittenberg

Geodatensätze aus Teilplänen Wind 2018 und 2027 der Regionalen
Planungsgemeinschaft Anhalt-Bitterfeld-Wittenberg.

### Sachlicher Teilplan Wind 2018

Geodaten aus rechtskräftigem
[Sachlichen Teilplan Wind 2018](https://www.planungsregion-abw.de/regionalplanung/teilplan-windenergie/teilplan-2018/).

> Im Sachlichen Teilplan "Nutzung der Windenergie in der Planungsregion
> Anhalt-Bitterfeld-Wittenberg" vom 30.05.2018 werden 22 Vorranggebiete für die
> Nutzung der Windenergie mit der Wirkung von Eignungsgebieten festgelegt. Sie
> dienen der raumordnerischen Steuerung der Errichtung von raumbedeutsamen
> Windenergieanlagen in Konzentrationszonen.
>
> Die oberste Landesentwicklungsbehörde hat am 01.08.2018 die Genehmigung
> erteilt. Mit Bekanntmachung der Genehmigung tritt der Sachliche Teilplan in
> Kraft.

Dateien
- Vorrang-/Eignungsgebiete: `stp_2018_vreg.gpkg`
  ([Quelle](https://gis.planungsregion-abw.de/geoserver/stp_wind2018/ows?SERVICE=WFS&REQUEST=GetCapabilities))

### Sachlicher Teilplan Wind 2027

Geodaten aus Planentwurf des
[Sachlichen Teilplan Wind 2027](https://www.planungsregion-abw.de/regionalplanung/teilplan-windenergie/teilplan-2027/).

> Die Regionalversammlung hat am 03.03.2023 beschlossen, den Sachlichen
> Teilplan "Windenergie 2027 in der Planungsregion Anhalt-Bitterfeld-Wittenberg"
> aufzustellen und mit der Bekanntgabe der Allgemeinen Planungsabsicht die
> beabsichtigten Auswahlkriterien und mögliche Gebietskulisse der Vorranggebiete
> für die Nutzung der Windenergie bzw. für Repowering von Windenergieanlagen
> vorzustellen.

Dateien
- Suchräume: `stp_2027_suchraum.gpkg` (Quelle: RPG ABW)
- Planabsicht Vorranggebiete: `stp_2027_ideen_vr.gpkg` (Quelle: RPG ABW)
- Planabsicht Repoweringgebiete: `stp_2027_ideen_repower.gpkg` (Quelle: RPG ABW)

**Dataset: `raw/rpg_abw_regional_plan`**

??? metadata "Metadata"
    ```json
    {
        "Quellen": {
            "Sachlicher Teilplan Wind 2018": "https://www.planungsregion-abw.de/regionalplanung/teilplan-windenergie/teilplan-2018/",
            "Geodaten Sachlicher Teilplan Wind 2018": "https://gis.planungsregion-abw.de/geoserver/stp_wind2018/ows?SERVICE=WFS&REQUEST=GetCapabilities",
            "Sachlicher Teilplan Wind 2027": "https://www.planungsregion-abw.de/regionalplanung/teilplan-windenergie/teilplan-2027/"
        }
    }
    ```

------------------------------
## Emissionen

Emissionen für die Jahre 1990 und 2019 für Sachsen-Anhalt (aus
[THG-Bericht 2021](https://lau.sachsen-anhalt.de/fileadmin/Bibliothek/Politik_und_Verwaltung/MLU/LAU/Wir_ueber_uns/Publikationen/Fachberichte/Dateien/221014_THG-Bericht.pdf))
und disaggregiert für die Region ABW.

Datei: `emissions.csv`, Felder:
- `sector`: Sektor
- `cat`: Kategorie ("*" = alle)
- `subcat`: Unterkategorie ("*" = alle)
- `name`: Bezeichner
- `st`: Emissionen Sachsen-Anhalt in kt CO2-Äquivalent
- `abw`: Emissionen Region ABW in kt CO2-Äquivalent

`sector`, `cat` und `subcat` folgen der Nomenklatur des Common Reporting Formats
(CRF) nach [KSG Anlage 1](https://www.gesetze-im-internet.de/ksg/anlage_1.html).
[Grafik hierzu](https://expertenrat-klima.de/content/uploads/2023/05/ERK2023_Pruefbericht-Emissionsdaten-des-Jahres-2022.pdf)
(Abb. 2 auf S. 30).

### Disaggregation

Hier Beschreibungstext

#### Sektor Energiewirtschaft (CRF 1.A.1 + 1.B)

####### CRF 1.A.1

EnbG: Emissionen aus europäischem Emissionshandel

####### CRF 1.B

EnbG: Emissionen aus europäischem Emissionshandel

#### Sektor Industrie (CRF 1.A.2 + 2)

####### CRF 1.A.2

EnbG: Energienutzung nach Ennergieträgern

####### CRF 2 Prozessemissionen

EnbG: in Industrie beschäftigte Personen

#### Sektor Verkehr (CRF 1.A.3)

EnbG:

* Zugelassene Kraftfahrzeuge
* gewichtet mit durchschn. Fahrleistung und spez. CO2 Emission pro km und Fahrzeugklasse

#### Sektor Sonstige Energie (insbes. Gebäude) (CRF 1.A.4 + 1.A.5)

EnbG: Wärmebedarf aus Energiesystem

#### Sektor Landwirtschaft (CRF 3)

####### CRF 3.A - Landwirtschaft – Fermentation

EnbG: Viehbestände

####### CRF 3.B-J:

EnbG: landwirtschaftlich genutzte Fläche

#### Sektor Abfall und Abwasser (CRF 5)

EnbG: Bevölkerung ABW

**Dataset: `raw/emissions`**

??? metadata "Metadata"
    ```json
    {
        "Daten Sachsen-Anhalt": "https://lau.sachsen-anhalt.de/fileadmin/Bibliothek/Politik_und_Verwaltung/MLU/LAU/Wir_ueber_uns/Publikationen/Fachberichte/Dateien/221014_THG-Bericht.pdf",
        "Datens\u00e4tze Desaggregation": {
            "Industrie": ""
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
## Installierte Leistungen von Biomasse-Konversionstechnologien

Die installierten Leistungen in MW wird im Szenario 80 % Transformationspfad
und 2,6 Mio. ha Anbauflächen im Jahr 2020 und 2050 der Tabelle 13 im
Dokument ["Technoökonomische Analyse und Transformationspfade des energetischen Biomassepotentials (TATBIO)"](../dbfz_biomass_heat_capacities/metadata.json)
für die folgenden Konversionsanlagen von Biomasse entnommen:

- Biomethan-Blockheizkraftwerk
- Holzhackschnitzelkessel Sektor Industrie
- Pelletkessel Sektor GHD
- Holzhackschnitzelkessel Sektor GHD
- Scheitholzvergaserkessel
- Pelletkessel Sektor Gebäude
- Biogasanlage + Blockheizkraftwerk
- Biomethan Gas- und Dampfkombikraftwerk
- Klärschlammfaulung + Blockheizkraftwerk
- Papier-Zellstoff-KWK
- Holzvergaser + Blockheizkraftwerk
- Mikro-Holzgas-Blockheizkraftwerk

Die Konversionstechnologien sind in der Spalte "technology" gelistet, während
sich ihre installierten Leistungen für die beiden Projektionsjahre in den
Spalten "capacity_[MW]_2020" und "capacity_[MW]_2050" befinden.

In den Spalten "decentral" und "central" wird mit "x" angegeben, ob jeweils ein
dezentraler und zentraler Einsatz der Konversionsanlage Stand der Technik ist.

In der Spalte "carrier" wird analog zur Konvention der Namensgebung im
Energiesystem (siehe [esys.md](../../../../docs/sections/esys.md)) der
jeweilige in die Konversionsanlage eintretende Energieträger notiert.
Diese werden Abbildung 3 des Dokuments entommen. Der Energieträger Schwarzlauge
wird vereinfachend dem Energieträger feste Biomasse bzw. Holz zugeordnet.
Klärgas und Holzgas werden vereinfachend Biogas zugeordnet.

In der Spalte "tech" findet die Zuordnung zu der Technologie anhand der im
Energiesystem verwendeten Komponenten (siehe
[esys.md](../../../../docs/sections/esys.md)) statt.

**Dataset: `raw/dbfz_biomass_heat_capacities`**

??? metadata "Metadata"
    ```json
    {
        "Quellen": {
            "Titel": "Techno\u00f6konomische Analyse und Transformationspfade des energetischen Biomassepotentials (TATBIO)",
            "Datei": "https://www.ufz.de/export/data/2/231891_technooekonomische-analyse-und-transformationspfade-des-energetischen-biomassepotentials(1).pdf",
            "Datum": "08.05.2019",
            "Autor": "DBFZ Deutsches Biomasseforschungszentrum gemeinn\u00fctzige GmbH",
            "Seiten": "6 und 54"
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
## Bevölkerung

Einwohnerzahl nach Gemeinden des Statistischen Bundesamts.

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
## Technologiedaten

### Jahresvolllaststunden

Anhand typischer heutiger und prognostizierter Werte für Sachsen-Anhalt werden
folgende Jahresvolllaststunden angenommen:

| Technologie     | Jahr | Volllaststunden | Quelle(n) für Annahme                                                                                                                                                                                                                                                                                       | Anmerkung                                                      |
|-----------------|------|----------------:|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----------------------------------------------------------------|
| Windenergie     | 2022 |            1800 | [foederal-erneuerbar](https://www.foederal-erneuerbar.de/landesinfo/bundesland/ST/kategorie/wind/auswahl/811-durchschnittliche_ja/#goto_811)                                                                                                                                                                |                                                                |
|                 | 2045 |            2300 | [PV- und Windflächenrechner](https://zenodo.org/record/6794558)                                                                                                                                                                                                                                             |                                                                |
| Freiflächen-PV  | 2022 |             980 | [foederal-erneuerbar](https://www.foederal-erneuerbar.de/landesinfo/bundesland/ST/kategorie/solar/auswahl/813-durchschnittliche_ja/#goto_813), [ISE](https://www.ise.fraunhofer.de/content/dam/ise/de/documents/publications/studies/aktuelle-fakten-zur-photovoltaik-in-deutschland.pdf)                   |                                                                |
|                 | 2045 |             980 | [PV- und Windflächenrechner](https://zenodo.org/record/6794558), [Ariadne Szenarienreport](https://ariadneprojekt.de/media/2022/02/Ariadne_Szenarienreport_Oktober2021_corr0222_lowres.pdf)                                                                                                                 |                                                                |
| Aufdach-PV      | 2022 |             910 | [foederal-erneuerbar](https://www.foederal-erneuerbar.de/landesinfo/bundesland/ST/kategorie/solar/auswahl/813-durchschnittliche_ja/#goto_813), [ISE](https://www.ise.fraunhofer.de/content/dam/ise/de/documents/publications/studies/aktuelle-fakten-zur-photovoltaik-in-deutschland.pdf)                   |                                                                |
|                 | 2045 |             910 | [Ariadne Szenarienreport](https://ariadneprojekt.de/media/2022/02/Ariadne_Szenarienreport_Oktober2021_corr0222_lowres.pdf)                                                                                                                                                                                  |                                                                |
| Laufwasserkraft | 2022 |            3800 | [foederal-erneuerbar](https://www.foederal-erneuerbar.de/landesinfo/bundesland/ST/kategorie/wasser/auswahl/840-durchschnittliche_ja/#goto_840)                                                                                                                                                              |                                                                |
|                 | 2045 |            3800 | [foederal-erneuerbar](https://www.foederal-erneuerbar.de/landesinfo/bundesland/ST/kategorie/wasser/auswahl/840-durchschnittliche_ja/#goto_840)                                                                                                                                                              |                                                                |
| Bioenergie      | 2022 |            6000 | [foederal-erneuerbar](https://www.foederal-erneuerbar.de/landesinfo/bundesland/ST/kategorie/bioenergie/auswahl/814-durchschnittliche_ja/#goto_814), [ISE](https://www.ise.fraunhofer.de/content/dam/ise/de/documents/publications/studies/DE2018_ISE_Studie_Stromgestehungskosten_Erneuerbare_Energien.pdf) | Bioenergie-Stromerzeugung (ohne<br/>biogenen Teil des Abfalls) |
|                 |      |                 |                                                                                                                                                                                                                                                                                                             |                                                                |

Datei: `technology_data.json` -> `full_load_hours`

TBD: Generalisieren - automatische Generierung anhand von Global Wind Atlas /
Global Solar Atlas.

### Leistungsdichte

Installierbare Leistung pro Fläche / spezifischer Flächenbedarf:
- Windenergie: 21 MW/km²
- PV-Freiflächenanlagen: 100 MW/km²
- PV-Aufdachanlagen: 140 MW/km²
- Solarthermie: ? MW/km²

Quelle: [PV- und Windflächenrechner](https://zenodo.org/record/6794558)

Datei: `technology_data.json` -> `power_density`

### Kosten und Wirkungsgrade

Datei: `raw_costs_efficiencies.csv`

**Dataset: `raw/technology_data`**

??? metadata "Metadata"
    ```json
    {
        "Datenquellen": {
            "FLH": "",
            "spec_area": "",
            "emissions": "https://ens.dk/en/our-services/projections-and-models/technology-data"
        }
    }
    ```

------------------------------
## OpenStreetMap

OpenStreetMap Datenauszug Deutschland.

Quelle: https://download.geofabrik.de/europe/germany-230101.osm.pbf

Ist nicht Teil des Eingangsdaten-Packages - manueller Download erforderlich.

**Dataset: `raw/osm`**

??? metadata "Metadata"
    ```json
    {
        "name": "openstreetmap",
        "title": "",
        "id": "openstreetmap",
        "description": "OpenStreetMap extract",
        "language": [
            "de-DE",
            "en-GB"
        ],
        "subject": [],
        "keywords": [
            "openstreetmap",
            "osm"
        ],
        "publicationDate": "2023-06-30",
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
            "resolution": ""
        },
        "temporal": {
            "referenceDate": "2023-06-30",
            "timeseries": []
        },
        "sources": [
            {
                "title": "OpenStreetMap Data Extracts (Geofabrik)",
                "description": "Full data extract of OpenStreetMap data",
                "path": "https://download.geofabrik.de/europe/germany-230101.osm.pbf",
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
                "date": "2023-06-30",
                "object": "metadata",
                "comment": "Create metadata"
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
## EE-Einspeisezeitreihen

Einspeisezeitreihen für Erneuerbare Energien, normiert auf 1 MW bzw. 1 p.u.
Als Wetterjahr wird 2011 verwendet, siehe
[Szenarien](../../../../docs/sections/scenarios.md).

### Windenergie

Stündlich aufgelöste Zeitreihe der Windenergie Einspeisung über 1 Jahr auf Basis
von [MaStR](../bnetza_mastr/dataset.md) und
[renewables.ninja](http://renewables.ninja).
Auf einen Auflösung auf Gemeindeebene wird verzichtet, da die Differenz der
Produktion der Gemeinden nach renewables.ninja <5 % beträgt.

#### Windenergieanlage (2022)

Für renewables.ninja sind Position (lat, lon), Nennleistung (capacity),
Nabenhöhe und Turbinentyp erforderlich.

##### Position

Hierfür wird aus den Zentroiden der Gemeinden ein räumlicher Mittelwert
anhand des Datensatzes
[bkg_vg250_muns_region](../../datasets/bkg_vg250_muns_region/dataset.md)
(`bkg_vg250_muns_region.gpkg`) gebildet:

```
import geopandas as gpd
import os.path

def get_position(gdf):
    df = gpd.read_file(gdf)
    points_of_muns = df["geometry"].centroid
    points_of_muns_crs = points_of_muns.to_crs(4326)
    point_df = [
        points_of_muns_crs.y.sum()/len(points_of_muns),
        points_of_muns_crs.x.sum()/len(points_of_muns)
    ]
    return point_df

data_folder = os.path.join("your_data_folder")
muns_gpkg = os.path.join(data_folder, "bkg_vg250_muns_region.gpkg")
center_position = get_position(muns_gpkg)
```

##### Nennleistung

Wird auf 1 MW gesetzt/normiert.

##### Nabenhöhe

Aus dem Datensatz
[bnetza_mastr_wind_region](../../datasets/bnetza_mastr_wind_region/dataset.md)
(`bnetza_mastr_wind_agg_abw.gpkg`) wird ein Mittelwer von 100 m abgeleitet.

```
import geopandas as gpd

df = gpd.read_file("bnetza_mastr_wind_agg_abw.gpkg")
height = df[["hub_height"]].mean()
```

##### Turbinentyp

Annahme: Innerhalb eines Herstellers sind Leistungskurven sehr ähnlich.
Daher werden zwei größten Hersteller mit jeweiligen häufigsten Turbinentyp
ausgewählt - diese sind Enercon und Vestas mit ca. 70 % und ca. 30%.

```
import geopandas as gpd

df = gpd.read_file("bnetza_mastr_wind_agg_abw.gpkg")
manufacturers = df[
    ["manufacturer_name", "status"]
].groupby("manufacturer_name").count().sort_values(
    by="status", ascending=False
)
```

Häufigste Turbinentypen sind *Enercon E-70* und *Vestas V80*. Daher werden
*Enercon E70 2000* und *Vestas V80 2000* in renewables.ninja ausgewählt.

```
man_1 = manufacturers.index[0]
man_2 = manufacturers.index[1]

type_1 = df[
    ["manufacturer_name", "type_name", "status"]
].where(df["manufacturer_name"] == man_1).groupby(
    "type_name").count().sort_values(by="status", ascending=False)

type_2 = df[
    ["manufacturer_name", "type_name", "status"]
].where(df["manufacturer_name"] == man_2).groupby(
    "type_name").count().sort_values(by="status", ascending=False)
```

#### Raw Data von [renewables.ninja](http://renewables.ninja) API

Es werden zwei Zeitreihen für oben beschriebenen Vergleichsanlagen berechnet:

```
import json
import requests
import pandas as pd
import geopandas as gpd

def change_wpt(position, capacity, height, turbine):
    args = {
        'lat': 51.8000,  # 51.5000-52.0000
        'lon': 12.2000,  # 11.8000-13.1500
        'date_from': '2011-01-01',
        'date_to': '2011-12-31',
        'capacity': 1000.0,
        'height': 100,
        'turbine': 'Vestas V164 7000',
        'format': 'json',
        'local_time': 'true',
        'raw': 'false',
    }

    args['capacity'] = capacity
    args['height'] = height
    args['lat'] = position[0]
    args['lon'] = position[1]
    args['turbine'] = turbine

    return args

def get_df(args):
    token = 'Please get your own'
    api_base = 'https://www.renewables.ninja/api/'

    s = requests.session()
    # Send token header with each request
    s.headers = {'Authorization': 'Token ' + token}

    url = api_base + 'data/wind'

    r = s.get(url, params=args)

    parsed_response = json.loads(r.text)
    df = pd.read_json(
    json.dumps(parsed_response['data']),orient='index')
    metadata = parsed_response['metadata']
    return df

enercon_production = get_df(change_wpt(
    position,
    capacity=1,
    height=df[["hub_height"]].mean(),
    turbine=enercon)
)

vestas_production = get_df(change_wpt(
    position,
    capacity=1000,
    height=df[["hub_height"]].mean(),
    turbine=vestas)
)
```

#### Gewichtung und Skalierung der Zeitreihen

Um die Charakteristika der beiden o.g. Anlagentypen zu berücksichtigen, erfolgt
eine gewichtete Summierung der Zeitreihen anhand der berechneten Häufigkeit.

#### Zukunftsszenarien

Analog zu dem oben beschriebenen Vorgehen wird eine separate Zeitreihe für
zukünftige WEA berechnet. Hierbei wird eine Enercon E126 6500 mit einer
Nabenhöhe von 159 m angenommen
([PV- und Windflächenrechner](https://zenodo.org/record/6794558)).

Da die Zeitreihe sich nur marginal von der obigen Status-quo-Zeitreihe
unterscheidet, wird letztere sowohl für den Status quo als auch die
Zukunftsszenarien verwendet.

- Einspeisezeitreihe: `wind_feedin_timeseries.csv`

### Freiflächen-Photovoltaik

#### PV-Anlage (2022)

Stündlich aufgelöste Zeitreihe der Photovoltaikeinspeisung über 1 Jahr auf Basis
von [MaStR](../bnetza_mastr/dataset.md) und
[renewables.ninja](http://renewables.ninja).
Wie bei der Windeinspeisung wird auf eine Auflsöung auf Gemeindeebene aufgrund
geringer regionaler Abweichungen verzichtet.

Für die Generierung der Zeitreihe über
[renewables.ninja](http://renewables.ninja)
wird eine Position(lat, lon), Nennleistung (capacity), Verluste (system_loss)
Nachführung (tracking), Neigung (tilt) und der Azimutwinkel (azim) benötigt.

Als Position wird analog zur Windenergieanlage der räumlicher Mittelwert
verwendet. Laut MaStR werden lediglich 13 Anlagen nachgeführt (0,01 % der
Kapazität), die Nachführung wird daher vernachlässigt. Die Neigung ist aus MaStR
nicht bekannt, es dominieren jedoch Anlagen auf Freiflächen sowie Flachdächern
im landwirtschaftlichen Kontext. Nach
[Ariadne Szenarienreport](https://ariadneprojekt.de/media/2022/02/Ariadne_Szenarienreport_Oktober2021_corr0222_lowres.pdf)
wird diese mit 30° angenommen.
Die Nennleistung Wird auf 1 MW gesetzt/normiert.

#### Zukunftsszenarien

Die Status-quo-Zeitreihe wird sowohl für den Status quo als auch die
Zukunftsszenarien verwendet.

- Einspeisezeitreihe: `pv_feedin_timeseries.csv`

### Solarthermie

- Einspeisezeitreihe: `st_feedin_timeseries.csv` (Kopie von
  PV-Einspeisezeitreihe)

### Laufwasserkraft

Hier wird eine konstante Einspeisung angenommen.

- Einspeisezeitreihe: `ror_feedin_timeseries.csv`

**Dataset: `raw/renewables.ninja_feedin`**

??? metadata "Metadata"
    ```json
    {
        "Quellen": {
            "renewables.ninja": "https://www.renewables.ninja/about",
            "Marktstammdatenregister": "siehe dataset bnetza_mastr"
        }
    }
    ```

------------------------------
## Dachflächenpotenzial PV-Aufdachanlagen in ABW

Abschätzung der installierten Leistung und des Ertrags von PV-Aufdachanlagen in
Anhalt-Bitterfeld-Wittenberg der Regionalen Planungsgemeinschaft.

Dafür wurden auf Basis des
[Digitalen Oberflächenmodells (DOM2)](https://www.lvermgeo.sachsen-anhalt.de/de/dom2-landesweit.html)
Schattenberechnungen durchgeführt. Anhand des
[LoD2 3D-Gebäudemodells](https://www.lvermgeo.sachsen-anhalt.de/de/download_lod2.html)
wurden für verschiedene Dachausrichtungen (nord, ost, süd, west, flach) die
installierbare Leistung bestimmt und mittels der Globalstrahlung und typischer
technischer Parameter für jedes Gebäude und jede Dachflächenorientierung
potenzielle Erträge berechnet.

Quellen
- [Hauptseite](https://www.planungsregion-abw.de/geodaten/)
- [Geodaten](https://gis-entwicklung2.planungsregion-abw.de/geoserver/wfs?SERVICE=WFS&REQUEST=GetCapabilities)
- [Anwendung](https://ris.planungsregion-abw.de/mapbender/application/pv_dachflaechenpot_rpg_abw)

**Dataset: `raw/rpg_abw_pv_roof_potential`**

??? metadata "Metadata"
    ```json
    {
        "Quellen": {
            "Hauptseite": "https://www.planungsregion-abw.de/geodaten/",
            "Geodaten": "https://gis-entwicklung2.planungsregion-abw.de/geoserver/wfs?SERVICE=WFS&REQUEST=GetCapabilities",
            "Tool": "https://ris.planungsregion-abw.de/mapbender/application/pv_dachflaechenpot_rpg_abw"
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
