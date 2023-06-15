# 'Preprocessed' Datasets 

------------------------------
## Energiedaten Sachsen-Anhalt

Datensätze zur Energie- und Wasserversorgung des Statistischen Landesamtes
Sachsen-Anhalt, extrahiert und konvertiert.

### Daten

Stromverbrauch der Industriebetriebe nach Kreisen 2003-2021 in MWh
- Datei: `power_demand_industry_st_districts.csv`

Raw dataset:
[stala_st_energy](../../raw/stala_st_energy/dataset.md)

**Dataset: `preprocessed/stala_st_energy`**

??? metadata "Metadata"
    ```json
    {
        "Datenquellen": {
            "Stromverbrauch der Industriebetriebe nach Kreisen": "https://statistik.sachsen-anhalt.de/themen/wirtschaftsbereiche/energie-und-wasserversorgung/tabellen-energieverwendung#c206986"
        }
    }
    ```

------------------------------
## sEEnergies Pan-European Thermal Atlas 5.2 (Peta5)

Wärmebedarf (extrahiert) für Europa 2015 in GJ (1ha Auflösung) für
- Haushalte: Raumwärme und Warmwasser
- GHD: Raumwärme, Warmwasser und Prozesswärme

**Dataset: `preprocessed/seenergies_peta5`**


------------------------------
## BMWK Langfristszenarien

Langfristszenarien des Bundesministerium für Wirtschaft und Klimaschutz, Daten
auf Landesebene, extrahiert.

Raw dataset:
[bmwk_long_term_scenarios](../../raw/bmwk_long_term_scenarios/dataset.md)

**Dataset: `preprocessed/bmwk_long_term_scenarios`**


------------------------------
## AGEB – Anwendungsbilanzen für die Endenergiesektoren 2011 bis 2021

Detaillierte Anwendungsbilanzen der Endenergiesektoren für 2020 und 2021 sowie
zusammenfassende Zeitreihen zum Endenergieverbrauch nach Energieträgern und
Anwendungszwecken für Jahre von 2011 bis 2021 der AG Energiebilanzen.

Aus PDF extrahierte Tabellenwerte für Haushalte, GHD und Industrie.

**Dataset: `preprocessed/ageb_energy_balance`**


------------------------------
## DemandRegio

Regionalisierte Bevölkerungsprognose sowie Strom-, Wärme und Gasbedarf auf
Landkreisebene, extrahiert.

Enthält Jahresverbräuche und Zeitreihen für die Sektoren Haushalte, Gewerbe,
Handel, Dienstleistungen (GHD) und Industrie für mehrere Zieljahre.

**Dataset: `preprocessed/demandregio`**


------------------------------
## OpenStreetMap gefiltert

OSM data nach bestimmten tags gefiltert, zu LAEA Europe (EPSG:3035) umprojiziert
und in ein Geopackage konvertiert.

**Dataset: `preprocessed/osm_filtered`**


------------------------------
## Sozialversicherungspflichtig Beschäftigte und Betriebe

Gemeindedaten der sozialversicherungspflichtig Beschäftigten am 30.06.2022 nach
Wohn- und Arbeitsort - Deutschland, Länder, Kreise und Gemeinden (Jahreszahlen)
der Bundesagentur für Arbeit.
Anzahl Beschäftigte und Betriebe extrahiert und in CSV konvertiert.

**Dataset: `preprocessed/ba_employment`**


------------------------------
## Erzeugungsanlagen aus Marktstammdatenregister

Erzeugungsanlagen aus dem MaStR für ausgewählte Technologien.

**Dataset: `preprocessed/bnetza_mastr`**


------------------------------
## Temperatur

Stündliche Mittelwerte der Luft- und Erdbodentemperatur für die Region ABW,
Mittelwert für alle Gemeinden.

Verwendet: [dwd_temperature](../../raw/dwd_temperature/dataset.md)

**Dataset: `preprocessed/dwd_temperature`**


------------------------------
## Regionalstatistik (GENESIS)

Enthält Datensätze der statistischen Ämter des Bundes und der Länder aus
[regiostat](../../raw/regiostat/dataset.md).

### Energieverwendung der Betriebe im Verarbeitenden Gewerbe (43531-01-02-4)

Jahreserhebung ü. die Energieverwendung der Betriebe im verarbeitendem Gewerbe.

Änderungen:
- Dateiformat konvertiert
- Bundesland-, Kreis und Gemeindewerte extrahiert
- Energie in TWh konvertiert

### Betriebe, tätige Personen, Bruttoentgelte (42111-01-04-5)

Jahreserhebung ü. Betriebe, tätige Personen und Bruttoentgelte der Betriebe im
verarbeitendem Gewerbe.

Änderungen:
- Dateiformat konvertiert
- Bundesland-, Kreis und Gemeindewerte extrahiert

**Dataset: `preprocessed/regiostat`**


------------------------------
## Lokale Verwaltungseinheiten

Lokale Verwaltungseinheiten (LAUs) von Eurostat, mit NUTS kompatibel. Diese
LAUs sind die Bausteine der NUTS und umfassen die Gemeinden und Kommunen der
Europäischen Union.

Daten aus Excel extrahiert und in CSV exportiert.

**Dataset: `preprocessed/eurostat_lau`**


------------------------------
## Bevölkerungsprognose Sachsen-Anhalt

Bevölkerungsprognose je Gemeinde bis 2035 des Statistischen Landesamtes
Sachsen-Anhalt, extrahiert und konvertiert.

Raw dataset:
[stala_st_pop_prog](../../raw/stala_st_pop_prog/dataset.md)

**Dataset: `preprocessed/stala_st_pop_prog`**


------------------------------
## Bevölkerung

Einwohnerzahl nach Gemeinden des Statistischen Bundesamts für die Jahre
2010, 2015, 2020, 2021, 2022.

**Dataset: `preprocessed/destatis_gv`**


------------------------------
## Administative areas of Germany

Geodata of administrative areas (Verwaltungsgebiete 1:250 000) extracted,
reprojected to LAEA Europe(EPSG:3035) and converted to Geopackage.

**Dataset: `preprocessed/bkg_vg250`**

