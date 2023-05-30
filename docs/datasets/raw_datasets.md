# 'Raw' Datasets 

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
## BMWK Langfristszenarien

Langfristszenarien des Bundesministerium für Wirtschaft und Klimaschutz, Daten
auf Landesebene.

Die Daten wurden über den
[Szenario Explorer](https://langfristszenarien.de/enertile-explorer-de/szenario-explorer/)
abgerufen.

Daten:
- Strom-, Wärme und Gasbedarf auf Deutschlandebene


### Haushalte: Strom

### Haushalte: Wärme/Gas


### GHD und Industrie: Strom


### GHD und Industrie: Wärme/Gas

Bedarfe und SLP-Zeitreihen je NUTS3

TBD

**Dataset: `raw/bmwk_long_term_scenarios`**

??? metadata "Metadata"
    ```json
    {
        "Datenquellen": {
            "Hauptseite": "https://langfristszenarien.de/enertile-explorer-de/szenario-explorer/",
            "Daten": "https://enertile-explorer.isi.fraunhofer.de:8443/open-view/52700/c6980ea467bb26a922d34617b4fd4798"
        }
    }
    ```

------------------------------
## DemandRegio

Regionalisierte Bevölkerungsprognose sowie Strom-, Wärme und Gasbedarf auf
Landkreisebene.

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

**Installation (digipipe venv aktiviert):**

```commandline
pip install openpyxl==3.1.0
pip install disaggregator@git+https://github.com/nesnoj/disaggregator.git#egg=disaggregator
```

Annahmen und Parameter
- Wetterjahr: Einheitlich 2011

## Details zum Datenabruf

### Bevölkerung

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

# Population
dr_hh_population = pd.DataFrame()
for year in [2010, 2015, 2017, 2020, 2021, 2022, 2025, 2030, 2035, 2040, 2045]:
    dr_hh_population[year] = round(data.population(year=year)).astype(int)

dr_hh_population.to_csv("dr_hh_population.csv")

# Households
data.households_per_size().to_csv("dr_hh_households_2011.csv")
```

### Haushalte: Strom

Bedarfe und SLP-Zeitreihen je NUTS3 mit Bottom-Up-Methode nach Haushaltsgröße.

Jahre
- 2017: Letzte verfügbare Daten
- 2022: Status quo, Fortschreibung mit Berücksichtigung Demografie und
  Wanderung
- 2035: Fortschreibungsjahr mit Berücksichtigung Demografie und Wanderung
- 2045: Fortschreibungsjahr

```python
from disaggregator import spatial, temporal

for year in [2017, 2022, 2035, 2045]:
  # Consumption
  spatial.disagg_households_power(
      by="households",
      weight_by_income=True,
      year=year,
      scale_by_pop=True,
  ).to_csv(f"dr_hh_power_demand_{year}.csv")

  # Timeseries
  temporal.disagg_temporal_power_housholds_slp(
      use_nuts3code=True,
      by="households",
      weight_by_income=True,
      year=year,
      scale_by_pop=True,
  ).to_csv(f"dr_hh_power_timeseries_{year}.csv")
```

### Haushalte: Wärme/Gas

Bedarfe und SLP-Zeitreihen je NUTS3

TBD

### GHD und Industrie: Strom

Bedarfe und Zeitreihen je NUTS3
- Bedarfe: je Wirtschaftszweig (WZ), abzüglich Eigenerzeugung
- Zeitreihen: für alle WZ aggregiert, Einzelprofile basieren je nach WZ
  auf gemessenen oder SLP inkl. Wanderung

Jahre
- 2017: Letzte verfügbare Daten
- 2022: Status quo, Fortschreibung mit Berücksichtigung Beschäftigte und
  Effizienzgewinne
- 2035: Max. Fortschreibungsjahr mit Berücksichtigung Beschäftigte und
  Effizienzgewinne

```python
from disaggregator import spatial, temporal

# CTS
for year in [2017, 2022, 2035]:
  # Consumption
  spatial.disagg_CTS_industry(
      sector='CTS',
      source='power',
      use_nuts3code=True,
      year=year,
  ).to_csv(f"dr_cts_power_demand_{year}.csv")
  # Timeseries
  temporal.disagg_temporal_power_CTS(
      detailed=False,
      use_nuts3code=True,
      year=year,
  ).to_csv(f"dr_cts_power_timeseries_{year}.csv")

# Industry
for year in [2017, 2022, 2035]:
  # Consumption
  spatial.disagg_CTS_industry(
      sector='industry',
      source='power',
      use_nuts3code=True,
      year=year,
  ).to_csv(f"dr_ind_power_demand_{year}.csv")
  # Timeseries
  temporal.disagg_temporal_industry(
      source="power",
      detailed=False,
      use_nuts3code=True,
      no_self_gen=False,
      year=year,
  ).to_csv(f"dr_ind_power_timeseries_{year}.csv")
```

### GHD und Industrie: Wärme/Gas

Bedarfe und SLP-Zeitreihen je NUTS3

TBD

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

Gemeindedaten der sozialversicherungspflichtig Beschäftigten nach Wohn- und
Arbeitsort - Deutschland, Länder, Kreise und Gemeinden (Jahreszahlen) der
Bundesagentur für Arbeit.

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
## Population prognosis

Population prognosis per municipality until 2035 by the Statistisches
Landesamt Sachsen-Anhalt.

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
