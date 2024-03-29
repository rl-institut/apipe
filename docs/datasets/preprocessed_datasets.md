# 'Preprocessed' Datasets 

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
[stala_st_pop_prog](../../apipe/store/raw/stala_st_pop_prog/dataset.md)

**Dataset: `preprocessed/stala_st_pop_prog`**


------------------------------

**Dataset: `preprocessed/bgr_sqr`**


------------------------------
## Dachflächenpotenzial PV-Aufdachanlagen in ABW

Abschätzung der installierten Leistung und des Ertrags von PV-Aufdachanlagen in
Anhalt-Bitterfeld-Wittenberg der Regionalen Planungsgemeinschaft, reprojizert.

Raw dataset:
[rpg_abw_pv_roof_potential](../../apipe/store/raw/rpg_abw_pv_roof_potential/dataset.md)

**Dataset: `preprocessed/rpg_abw_pv_roof_potential`**


------------------------------
## Temperatur

Stündliche Mittelwerte der Luft- und Erdbodentemperatur für die Region ABW,
Mittelwert für alle Gemeinden.

Verwendet: [dwd_temperature](../../apipe/store/raw/dwd_temperature/dataset.md)

**Dataset: `preprocessed/dwd_temperature`**


------------------------------

**Dataset: `preprocessed/oei_agri_pv`**


------------------------------
## BMWK Langfristszenarien

Langfristszenarien des Bundesministerium für Wirtschaft und Klimaschutz, Daten
auf Landesebene, extrahiert.

Raw dataset:
[bmwk_long_term_scenarios](../../apipe/store/raw/bmwk_long_term_scenarios/dataset.md)

**Dataset: `preprocessed/bmwk_long_term_scenarios`**


------------------------------
## AGEB – Anwendungsbilanzen für die Endenergiesektoren 2012 bis 2022

Detaillierte Anwendungsbilanzen der Endenergiesektoren für 2021 und 2022 sowie
zusammenfassende Zeitreihen zum Endenergieverbrauch nach Energieträgern und
Anwendungszwecken für Jahre von 2012 bis 2022 der AG Energiebilanzen.

Aus PDF extrahierte Tabellenwerte für Haushalte, GHD und Industrie.

**Dataset: `preprocessed/ageb_energy_balance`**


------------------------------
## sEEnergies Pan-European Thermal Atlas 5.2 (Peta5)

Wärmebedarf (extrahiert) für Europa 2015 in GJ (1ha Auflösung) für

- Haushalte: Raumwärme und Warmwasser
- GHD: Raumwärme, Warmwasser und Prozesswärme

**Dataset: `preprocessed/seenergies_peta5`**


------------------------------
## OpenStreetMap gefiltert

OSM data nach bestimmten Tags (s. [config.yml](../../apipe/store/preprocessed/osm_filtered/config.yml) --> `tags`) gefiltert,
zu LAEA Europe (EPSG:3035) umprojiziert und in ein Geopackage konvertiert.

**Achtung:** Konvertierungs- und Extraktionsprozess benötigt ~50 GB
Speicherplatz und kann viel Zeit in Anspruch nehmen.

**Dataset: `preprocessed/osm_filtered`**


------------------------------
## Sozialversicherungspflichtig Beschäftigte und Betriebe

Gemeindedaten der sozialversicherungspflichtig Beschäftigten am 30.06.2023 nach
Wohn- und Arbeitsort - Deutschland, Länder, Kreise und Gemeinden (Jahreszahlen)
der Bundesagentur für Arbeit.
Anzahl Beschäftigte und Betriebe extrahiert und in CSV konvertiert.

**Dataset: `preprocessed/ba_employment`**


------------------------------
## Bevölkerung

Einwohnerzahl nach Gemeinden des Statistischen Bundesamts für die Jahre
2010, 2015, 2020, 2021, 2022.

**Dataset: `preprocessed/destatis_gv`**


------------------------------
## Regionalstatistik (GENESIS)

Enthält Datensätze der statistischen Ämter des Bundes und der Länder aus
[regiostat](../../apipe/store/raw/regiostat/dataset.md).

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
## Geodaten PV- und Windflächenrechner

Geodaten aus dem [PV- und Windflächenrechner](https://www.agora-energiewende.de/service/pv-und-windflaechenrechner/),
extrahiert.

Raw dataset:
[rli_pv_windflaechenrechner](../../apipe/store/raw/rli_pv_wfr/dataset.md)

**Dataset: `preprocessed/rli_pv_wfr`**


------------------------------
## Regionalplan Anhalt-Bitterfeld-Wittenberg

Vorverarbeitete Datensätze aus Teilplänen Wind der Regionalen
Planungsgemeinschaft Anhalt-Bitterfeld-Wittenberg aus
[rpg_abw_regional_plan](../../apipe/store/raw/rpg_abw_regional_plan/dataset.md).

In der [config.yml](../../apipe/store/preprocessed/rpg_abw_regional_plan/config.yml) können Einstellungen vorgenommen werden.

**Dataset: `preprocessed/rpg_abw_regional_plan`**


------------------------------
## Erzeugungsanlagen aus Marktstammdatenregister

Erzeugungsanlagen aus dem MaStR für ausgewählte Technologien.

**Dataset: `preprocessed/bnetza_mastr`**


------------------------------
## Energiedaten Sachsen-Anhalt

Datensätze zur Energie- und Wasserversorgung des Statistischen Landesamtes
Sachsen-Anhalt, extrahiert und konvertiert.

### Daten

Stromverbrauch der Industriebetriebe nach Kreisen 2003-2021 in MWh

- Datei: `power_demand_industry_st_districts.csv`

Raw dataset:
[stala_st_energy](../../apipe/store/raw/stala_st_energy/dataset.md)

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
## DemandRegio

Regionalisierte Bevölkerungsprognose sowie Strom-, Wärme und Gasbedarf auf
Landkreisebene, extrahiert.

Enthält Jahresverbräuche und Zeitreihen für die Sektoren Haushalte, Gewerbe,
Handel, Dienstleistungen (GHD) und Industrie für mehrere Zieljahre.

**Dataset: `preprocessed/demandregio`**


------------------------------
## Regionalplan Oderland-Spree

Vorverarbeitete Datensätze aus Teilplänen Wind der Regionalen
Planungsgemeinschaft Oderland-Spree aus
[rpg_ols_regional_plan](../../apipe/store/raw/rpg_ols_regional_plan/dataset.md).

In der [config.yml](../../apipe/store/preprocessed/rpg_ols_regional_plan/config.yml) können Einstellungen vorgenommen werden.

**Dataset: `preprocessed/rpg_ols_regional_plan`**


------------------------------
## Solaratlas Brandenburg - Photovoltaikanlagen auf Dachflächen

Abschätzung der installierten Leistung und des Ertrags von PV-Aufdachanlagen in
Brandenburg der Wirtschaftsförderung Berlin-Brandenburg, reprojizert und
Attribute gefiltert und umbenannt.

Raw dataset:
[wfbb_pv_roof_potential](../../apipe/store/raw/wfbb_pv_roof_potential/dataset.md)

**Dataset: `preprocessed/wfbb_pv_roof_potential`**


------------------------------
## Anteile von Biomasse-Konversionsanlagen anhand installierter Leistung

Berechnung der Anteile der installierten Leistung an der gesamten installierten
Leistung der Biomasse-Konversionsanlagen.

Die installierten Leistungen werden
[dbfz_biomass_heat_capacities](../../apipe/store/raw/dbfz_biomass_heat_capacities/dataset.md)
entnommen. Sie werden nach Energieträger (Biogas, Methan oder Holz) und
Technologie (BHKW (bpchp), Turbine mit Kondensationsentnahme (extchp) oder
Ofen (oven)) zusammengefasst. Anschließend wird der Anteil der installierten
Leistung an der gesamten installierten Leistung der Biomasse-Konversionsanlagen
berechnet. Der Einfachheit halber werden die Projektionen für 2050 dem Jahr
2045 und die für 2020 dem Jahr 2022 zugeordnet. Der Energieträger und die
Technologie (vgl. [dbfz_biomass_heat_capacities](../../apipe/store/raw/dbfz_biomass_heat_capacities/dataset.md))
werden in einer Spalte zusammengefasst.

**Dataset: `preprocessed/dbfz_biomass_capacity_rel`**


------------------------------
## Administative areas of Germany

Geodata of administrative areas (Verwaltungsgebiete 1:250 000) extracted,
reprojected to LAEA Europe(EPSG:3035) and converted to Geopackage.

**Dataset: `preprocessed/bkg_vg250`**

