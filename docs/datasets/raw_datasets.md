# 'Raw' Datasets 

------------------------------
## OpenStreetMap

OpenStreetMap data extract for Sachsen-Anhalt.

**Dataset: `raw/osm_sachsen-anhalt`**

??? metadata "Metadata"
    ```json
    {
        "TODO": "!! PLEASE CHECK IF THE JSON IS CORRECT !!",
        "name": "openstreetmap",
        "title": "",
        "id": "openstreetmap",
        "description": "OpenStreetMap extract for federal state of Sachsen-Anhalt",
        "language": [
            "en-GB",
            "de-DE"
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
        "@id": "",
        "@context": "",
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
## Administative areas of Germany

Administative areas of Germany (Verwaltungsgebiete 1:250 000).

**Dataset: `raw/bkg_vg250`**

??? metadata "Metadata"
    ```json
    {
        "Originale Datenquelle": "https://gdz.bkg.bund.de/index.php/default/digitale-geodaten/verwaltungsgebiete/verwaltungsgebiete-1-250-000-stand-01-01-vg250-01-01.html"
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
        "Originale Datenquelle": "https://www.marktstammdatenregister.de/MaStR/Datendownload",
        "Prozessiert mit": "https://github.com/OpenEnergyPlatform/open-MaStR/}"
    }
    ```

------------------------------
## Municipalities population

Municipalities with population from the Federal Statistical Office.

**Dataset: `raw/destatis_gv`**

??? metadata "Metadata"
    ```json
    {
        "Originale Datenquelle": "https://www.destatis.de/DE/Themen/Laender-Regionen/Regionales/Gemeindeverzeichnis/Administrativ/Archiv/GVAuszugQ/AuszugGV4QAktuell.html"
    }
    ```

------------------------------
## Population prognosis

Population prognosis until 2035 by the Statistisches Landesamt Sachsen-Anhalt.

**Dataset: `raw/stala_st_pop_prog`**

??? metadata "Metadata"
    ```json
    {
        "Originale Datenquelle": [
            "https://statistik.sachsen-anhalt.de/themen/bevoelkerung-mikrozensus-freiwillige-haushaltserhebungen/bevoelkerung/bevoelkerungsprognose-und-haushalteprognose/#c312231",
            "https://statistik.sachsen-anhalt.de/fileadmin/Bibliothek/Landesaemter/StaLa/startseite/Themen/Bevoelkerung/Tabellen/Bevoelkerungsprognose/1_Internettabelle_7RBP_nach_Prognosejahr_Geschlecht_alle_Ebenen.xlsx"
        ]
    }
    ```
