# 'Datasets' Datasets 

------------------------------
## Technologiedaten

Allgemeine Technologiedaten.

Raw dataset: [technology_data](../../apipe/store/raw/technology_data/dataset.md)

**Dataset: `datasets/technology_data`**


------------------------------
## Potenzialgebiete Windenergie

Potenzialgebiete für die Errichtung von Windenergieanlagen, basierend auf den
Teilplänen Wind der Regionalen Planungsgemeinschaft Anhalt-Bitterfeld-Wittenberg
aus
[rpg_abw_regional_plan](../../apipe/store/preprocessed/rpg_abw_regional_plan/dataset.md).

Dateien:

- STP Wind 2018 - Vorrang-/Eignungsgebiete:
  `potentialarea_wind_stp_2018_vreg.gpkg`
- STP Wind 2027 - Planabsicht Vorranggebiete:
  `potentialarea_wind_stp_2027_vr.gpkg`
- STP Wind 2027 - Planabsicht Repoweringgebiete:
  `potentialarea_wind_stp_2027_repowering.gpkg`
- STP Wind 2027 - Suchraum Wald:
  `potentialarea_wind_stp_2027_search_area_forest_area.gpkg`
- STP Wind 2027 - Suchraum Offenland:
  `potentialarea_wind_stp_2027_search_area_open_area.gpkg`

Die darin verwendeten Attributtexte werden in die Datei
`potentialarea_wind_attribute_captions.json` exportiert.

Die Flächen werden mit den Gemeindegrenzen verschnitten und den Gemeinden
zugeordnet. Je Gemeinde und obigem Flächentyp/Datei wird eine Flächensumme (in
km²) berechnet, siehe `potentialarea_wind_area_stats_muns.csv`. Die Gemeinden
werden über den Schlüssel `municipality_id` (vgl.
[bkg_vg250_muns_region](../../apipe/store/datasets/bkg_vg250_muns_region/dataset.md))
identifiziert.

**Dataset: `datasets/potentialarea_wind_region`**


------------------------------
## Speicheranlagen

Speicheranlagen in der Region aus MaStR-Registerdaten als Geopackage.
Es werden alle Anlagen berücksichtigt, die in Betrieb sind oder sich in
Planung befinden. Anlagen mit Geokoordinaten werden georeferenziert
übernommen, für Anlagen die keine Koordinaten aufweisen (üblicherweise <=30
kW Nennleistung) erfolgt ein Geocoding anhand von PLZ und Ort, um eine
ungefähre Position bereit zu stellen.

Es wird weiterhin geprüft, ob dem Speicher eine oder mehrere PV-Aufdachanlagen
zugeordnet sind, es wird die Anzahl und Summe der Nettonennleistung berechnet.

Neben einem anlagenscharfen Datensatz wird ein weiterer Datensatz erzeugt,
der alle Anlagen mit approximierter Position je Position zusammenfasst und
jeweils typische Kennwerte enthält (u.a. Anzahl Anlagen, Gesamtleistung).

Jede Anlage wird anhand ihrer Lokation einer Gemeinde (Attribut
`municipality_id`, vgl.
[bkg_vg250_muns_region](../../apipe/store/datasets/bkg_vg250_muns_region/dataset.md)) und
einem Landkreis (Attribut `district_id`, vgl.
[bkg_vg250_muns_region](../../apipe/store/datasets/bkg_vg250_districts_region/dataset.md))
zugeordnet.

Weiterhin erfolgt eine Auswertung der installierten Gesamtleistung je Gemeinde:

- Alle Speicher: `bnetza_mastr_storage_stats_muns.csv`
- Großspeicher (>=100 kWh): `bnetza_mastr_storage_large_stats_muns.csv`
- Kleinspeicher (<100 kWh): `bnetza_mastr_storage_small_stats_muns.csv`

`bnetza_mastr_storage_pv_roof.json` enthält die spezifische Speicherkapazität
sowie spezifische Nennleistung der Speicher (bezogen auf die installierte
Leistung von PV-Aufdachanlagen), aggregiert für gesamte Region, für folgende
Randbedingungen:

- Alle PV-Anlagen: `all_storages`
- PV-Anlagen mit 2..20 kWp sowie Batteriespeicher <20 kWh und <20 kW (kann in
  [config.yml](../../apipe/store/datasets/bnetza_mastr_storage_region/config.yml) unter `home_storages` konfiguriert werden):
  `home_storages`

**Dataset: `datasets/bnetza_mastr_storage_region`**


------------------------------
## EE-Einspeisezeitreihen

Einspeisezeitreihen für Erneuerbare Energien. Als Wetterjahr wird 2011
verwendet, siehe [Szenarien](../../apipe/store/../../docs/sections/scenarios.md).

Raw dataset mit methodischer Beschreibung:
[renewables.ninja_feedin](../../apipe/store/raw/renewables.ninja_feedin/dataset.md)

### Einspeisezeitreihen

Zeitreihe normiert auf Summe=1 für

- Windenergie: `wind_feedin_timeseries.csv`
- Photovoltaik: `pv_feedin_timeseries.csv`
- Solarthermie: `st_feedin_timeseries.csv`
- Laufwasserkraft: `ror_feedin_timeseries.csv`

**Dataset: `datasets/renewable_feedin`**


------------------------------
## Emissionen

Emissionen für Sachsen-Anhalt und die Region, aggregiert nach Sektoren der
CRF-Nomenklatur.

Datei `emissions.json` enthält Chartdaten.

Raw dataset: [emissions](../../apipe/store/raw/emissions/dataset.md)

**Dataset: `datasets/emissions_region`**


------------------------------
## Geo- oder Solarthermie-, Grubengas- und Klärschlamm-Anlagen

Anlagen der Geo- oder Solarthermie, Grubengas und Klärschlamm in der Region
aus MaStR-Registerdaten als Geopackage.
Es werden alle Anlagen berücksichtigt, die in Betrieb sind oder sich in
Planung befinden. Anlagen mit Geokoordinaten werden georeferenziert
übernommen, für Anlagen die keine Koordinaten aufweisen (üblicherweise <=30
kW Nennleistung) erfolgt ein Geocoding anhand von PLZ und Ort, um eine
ungefähre Position bereit zu stellen.

Neben einem anlagenscharfen Datensatz wird ein weiterer Datensatz erzeugt,
der alle Anlagen mit approximierter Position je Position zusammenfasst und
jeweils typische Kennwerte enthält (u.a. Anzahl Anlagen, Gesamtleistung).

Jede Anlage wird anhand ihrer Lokation einer Gemeinde (Attribut
`municipality_id`, vgl.
[bkg_vg250_muns_region](../../apipe/store/datasets/bkg_vg250_muns_region/dataset.md)) und
einem Landkreis (Attribut `district_id`, vgl.
[bkg_vg250_muns_region](../../apipe/store/datasets/bkg_vg250_districts_region/dataset.md))
zugeordnet.

Zusätzlich erfolgt eine statistische Auswertung der installierten Leistung in
`bnetza_mastr_gsgk_stats_muns.csv`.

**Dataset: `datasets/bnetza_mastr_gsgk_region`**


------------------------------
## Strombedarf

Nettostrombedarfe und -zeitreihen für Haushalte, GHD und Industrie je Gemeinde.

Die Berechnung der regionalen Prognosewerte je Verbrauchssektor erfolgt anhand
landesweiter Prognosen aus den
[BMWK Langfristszenarien](../../apipe/store/preprocessed/bmwk_long_term_scenarios/dataset.md).

### Haushalte

- Jährlicher Strombedarf je Gemeinde in MWh aus
  [DemandRegio](../../apipe/store/preprocessed/demandregio/dataset.md), von Landkreis- auf
  Gemeindeebene disaggregiert anhand von Bevölkerungsprognosen
  ([STALA ST](../../apipe/store/preprocessed/stala_st_pop_prog/dataset.md)).
- Prognosewerte für 2045 werden durch lineare Skalierung mittels Reduktion des
  Strombedarfs (ohne Wärmegewinnung) aus
  [BMWK Langfristszenarien](../../apipe/store/preprocessed/bmwk_long_term_scenarios/dataset.md)
  berechnet. Hierbei wird das Szenario "TN-Strom" als Grundlage für den Status
  quo verwendet und Werte für 2022 interpoliert. Die Zielwerte werden dem
  Szenario "T45-Strom" entnommen.
- Gemittelte, normierte Strombedarfszeitreihe (auf 1 MWh) aus
  [DemandRegio](../../apipe/store/preprocessed/demandregio/dataset.md)-Daten von 2022, die
  für alle Zielszenarien und Aggregationsebenen verwendet wird, da die Basis
  SLP-Profile sind und Differenzen zwischen verschiedenen Jahren nur aufgrund
  der Lage von Wochenenden und Feiertagen bestehen. Diese werden daher
  vernachlässigt.

### GHD

- Jährlicher Strombedarf je Gemeinde in MWh aus
  [DemandRegio](../../apipe/store/preprocessed/demandregio/dataset.md), von Landkreis- auf
  Gemeindeebene disaggregiert anhand von sozialversicherungspflichtig
  Beschäftigten im Jahr 2022
  ([BA für Arbeit](../../apipe/store/preprocessed/ba_employment/dataset.md)).
- Prognosewerte für 2045 werden durch lineare Skalierung mittels Reduktion des
  Strombedarfs (ohne Wärmegewinnung) aus
  [BMWK Langfristszenarien](../../apipe/store/preprocessed/bmwk_long_term_scenarios/dataset.md)
  berechnet. Hierbei wird das Szenario "TN-Strom" als Grundlage für den Status
  quo verwendet und Werte für 2022 interpoliert. Die Zielwerte werden dem
  Szenario "T45-Strom" entnommen.
- Gemittelte, normierte Strombedarfszeitreihe (auf 1 MWh) aus
  [DemandRegio](../../apipe/store/preprocessed/demandregio/dataset.md)-Daten von 2022, die
  für alle Zielszenarien und Aggregationsebenen verwendet wird. Basis bilden
  sowohl SLP- als auch branchenspezifische Profile. Aufgrund der geringen
  Differenzen zwischen den Landkreisen werden diese gemittelt. Differenzen
  zwischen verschiedenen Jahren bestehen nur aufgrund der Lage von Wochenenden
  und Feiertagen und werden daher vernachlässigt.

### Industrie

- Jährlicher Strombedarf je Gemeinde in MWh. Hierfür stehen 2 Datensätze zur
  Verfügung - welcher verwendet wird, kann in der [Konfiguration](../../apipe/store/datasets/demand_electricity_region/config.yml)
  via `ind_electricity_demand_source` eingestellt werden:
    - [DemandRegio](../../apipe/store/preprocessed/demandregio/dataset.md): Werte für alle
    Landkreise in Deutschland.
    - [STALA ST](../../apipe/store/preprocessed/stala_st_energy/dataset.md) (Standard):
    Genauere Werte, jedoch nur für Sachsen-Anhalt verfügbar.
- Die Desaggregation von Landkreis- auf Gemeindeebene erfolgt anhand der
  Beschäftigten im verarbeitenden Gewerbe im Jahr 2022
  ([Regionalstatistik](../../apipe/store/preprocessed/regiostat/dataset.md)).
- Prognosewerte für 2045 werden durch lineare Skalierung mittels Reduktion des
  industriellen Gesamtenergiebedarfs aus
  [BMWK Langfristszenarien](../../apipe/store/preprocessed/bmwk_long_term_scenarios/dataset.md)
  berechnet. Im Unterschied zu Haushalten und GHD liegen die Daten für den
  Wärme- und Stromanteil nicht getrennt vor, sodass auf den
  Gesamtenergiebedarf zurückgegriffen wird.
  Es wird das Szenario "TN-Strom" als Grundlage für den Status quo verwendet und
  Werte für 2022 interpoliert. Die Zielwerte werden dem Szenario "T45-Strom"
  entnommen.
- Gemittelte, normierte Strombedarfszeitreihe (auf 1 MWh) aus
  [DemandRegio](../../apipe/store/preprocessed/demandregio/dataset.md)-Daten von 2022, die
  für alle Zielszenarien und Aggregationsebenen verwendet wird. Basis bilden
  sowohl SLP- als auch branchenspezifische Profile. Aufgrund der geringen
  Differenzen zwischen den Landkreisen werden diese gemittelt. Differenzen
  zwischen verschiedenen Jahren bestehen nur aufgrund der Lage von Wochenenden
  und Feiertagen und werden daher vernachlässigt.

**Dataset: `datasets/demand_electricity_region`**


------------------------------
## Dachflächenpotenzial PV-Aufdachanlagen in ABW

Abschätzung der installierten Leistung und des Ertrags von PV-Aufdachanlagen in
Anhalt-Bitterfeld-Wittenberg der Regionalen Planungsgemeinschaft aus Datensatz
[rpg_abw_pv_roof_potential](../../apipe/store/raw/rpg_abw_pv_roof_potential/dataset.md).

Die Gebäudezentroide werden mit den Gemeindegrenzen verschnitten und den
Gemeinden zugeordnet.
Ergebnisdaten:

- Alle Gebäude: `potentialarea_pv_roof_area_stats_muns.csv`
- Alle nicht denkmalgeschützten Gebäude:
  `potentialarea_pv_roof_wo_historic_area_stats_muns.csv`

Des Weiteren wird je Gemeinde der relative Anteil der bereits installierten
Anlagenleistung an der theoretisch installierbaren Leistung (bei
100% Dachnutzung) berechnet.
Ergebnisdaten:

- Alle Gebäude: `potentialarea_pv_roof_deployment_stats_muns.csv`
- Alle nicht denkmalgeschützten Gebäude:
  `potentialarea_pv_roof_wo_historic_deployment_stats_muns.csv`

Die Gemeinden werden über den Schlüssel `municipality_id` (vgl.
[bkg_vg250_muns_region](../../apipe/store/datasets/bkg_vg250_muns_region/dataset.md))
identifiziert.

### Ausbauziele

Es werden PV-Ausbauziele für die Region berechnet, indem die Bundesziele aus den
[BMWK Langfristszenarien](../../apipe/store/preprocessed/bmwk_long_term_scenarios/dataset.md)
i.H.v. 428 GW
([§4 EEG 2023](https://www.gesetze-im-internet.de/eeg_2014/__4.html): 400 GW)
anhand der Gebäudegrundflächen disaggregiert werden. Hierzu wird der Anteil der
Gebäudegrundflächen in der Region an der bundesweiten Gebäudegrundflächen
berechnet (s. Datensatz [osm_buildings](../../apipe/store/datasets/osm_buildings/dataset.md)) und die
Ziele linear skaliert. Da in den o.g. Ausbauzielen nicht zwischen Freiflächen- und
Aufdach-PV unterschieden wird, wird folgende Aufteilung angenommen (änderbar in
[config.yml](../../apipe/store/datasets/potentialarea_pv_roof_region/config.yml)):

- Aufdach-PV: 52 % (221 GW)
- Freiflächen-PV (niedrig aufgeständert): 44 % (190 GW), vgl.
  [potentialarea_pv_ground_region2](../../apipe/store/datasets/potentialarea_pv_ground_region2/dataset.md)
- Agri-PV (hoch aufgeständert und vertikal bifazial): 4 % (17 GW)

File: `potentialarea_pv_roof_regionalized_targets.json`

**Dataset: `datasets/potentialarea_pv_roof_region`**


------------------------------
## Potenzialgebiete PV-Freiflächen (Deutschland gesamt)

Diese Analyse widmet sich der Identifizierung von Photovoltaik-Potenzialen auf
Agrarflächen, basierend auf einer Auswertung rasterbasierter Landnutzungsdaten.
Im Fokus stehen dabei die Bodenqualität und die Eignung landwirtschaftlicher
Dauerkulturen für die Integration von Agri-PV-Systemen.

### Datengrundlage und Methodik

Datengrundlage:

- SQR-Daten (Soil Quality Rating), s. Datensatz
  [bgr_sqr](../../apipe/store/raw/bgr_sqr/dataset.md).
  - **(I)** SQR Originaldaten
- Potenzialflächen Agri-PV , s. Datensatz
  [oei_agri_pv](../../apipe/store/raw/oei_agri_pv/dataset.md)
  - **(II)** SQR Gesamtpotenzial (`Agri-PV-Potenziale_Gesamt_100x100_EPSG3035`)
  - **(III)** SQR 50-70
    (`Agri-PV-Potenziale_SQR_50-70_NOT_Dauerkulturen_100x100_EPSG3035`)
  - **(IV)** Dauerkulturen (`Agri-PV-Potenziale_Dauerkulturen_100x100_EPSG3035`)

(III) und (IV) sind hierbei disjunkt.

Es werden folgende Flächen verwendet, auf welchen jeweils eine andere
technologische Umsetzung als Grundlage angenommen wird:

- **(A)** Auf Acker- und Grünlandflächen mit sehr geringer..geringer Bodengüte
  (SQR 0..50): Klassische, niedrig aufgeständerte FF-PV.
- **(B)** Auf Acker- und Grünlandflächen mit geringer..mittlerer Bodengüte (SQR
  50..70): Agri-PV - bifaziale, vertikal aufgeständerte PV
- **(C)** Dauerkulturen (Kernobst, Wein, Gemüse und Beeren): Agri-PV - hoch
   aufgeständerte PV-Systeme

Berechnung/Verschneidung

- **(A)** = **(II)**-(**(I)**) mit Wert>50)-(**(IV)**) mit Wert>0)
- **(B)** = **(III)**
- **(C)** = **(IV)**

### Ergebnisdaten

- **(A)**: `potentialarea_pv_ground_soil_quality_low.tif`
- **(B)**: `potentialarea_pv_ground_soil_quality_medium.tif`
- **(C)**: `potentialarea_pv_ground_permanent_crops.tif`

**Dataset: `datasets/potentialarea_pv_ground`**


------------------------------
## Biomasse-/Biogasanlagen

Biomasse-/Biogasanlagen in der Region aus MaStR-Registerdaten als Geopackage.
Es werden alle Anlagen berücksichtigt, die in Betrieb sind oder sich in
Planung befinden. Anlagen mit Geokoordinaten werden georeferenziert
übernommen, für Anlagen die keine Koordinaten aufweisen (üblicherweise <=30
kW Nennleistung) erfolgt ein Geocoding anhand von PLZ und Ort, um eine
ungefähre Position bereit zu stellen.

Neben einem anlagenscharfen Datensatz wird ein weiterer Datensatz erzeugt,
der alle Anlagen mit approximierter Position je Position zusammenfasst und
jeweils typische Kennwerte enthält (u.a. Anzahl Anlagen, Gesamtleistung).

Jede Anlage wird anhand ihrer Lokation einer Gemeinde (Attribut
`municipality_id`, vgl.
[bkg_vg250_muns_region](../../apipe/store/datasets/bkg_vg250_muns_region/dataset.md)) und
einem Landkreis (Attribut `district_id`, vgl.
[bkg_vg250_muns_region](../../apipe/store/datasets/bkg_vg250_districts_region/dataset.md))
zugeordnet.

Zusätzlich erfolgt eine statistische Auswertung der installierten Leistung in
`bnetza_mastr_biomass_stats_muns.csv`.

**Dataset: `datasets/bnetza_mastr_biomass_region`**


------------------------------
## Staat

Staatsgrenze aus Geodaten der Verwaltungsgebiete extrahiert und nach Landmasse
gefiltert (Geofaktor 4 = "mit Struktur Land").

**Dataset: `datasets/bkg_vg250_state`**


------------------------------
## Wärmepumpen COP

Zeitreihe für die Leistungszahl / Coefficient of performance (COP) für
Wärmepumpen. Berücksichtigt werden Luftwärmepumpen (ASHP) und Erdwärmepumpen
(GSHP). Der COP wird mit Hilfe von Zeitreihen der Umgebungstemperatur (ASHP)
bzw. der Bodentemperatur (GSHP) für jeden Zeitschritt berechnet.

Details zur Berechnungsmethodik können der Dokumentation von
[oemof.thermal](https://oemof-thermal.readthedocs.io/en/latest/compression_heat_pumps_and_chillers.html)
entnommen werden.

Annahmen:

- Vorlauftemperatur: 40 °C
- Gütegrad / Quality grade: 0.4 (nach
  [VDE](https://www.energiedialog2050.info/wp-content/uploads/simple-file-list/VDE_ST_ETG_Warmemarkt_RZ-web.pdf))
- Vereisungsverluste bei ASHP: 20 % bei <2 °C

Daraus ergibt sich eine mittlere Jahresarbeitszahl (JAZ) von 3,3 für ASHP und
4,3 für GSHP, die mit typischen Werten für 2019
([AEW](https://static.agora-energiewende.de/fileadmin/Projekte/2022/2022-04_DE_Scaling_up_heat_pumps/A-EW_273_Waermepumpen_WEB.pdf))
übereinstimmen. Für das Zukunftsszenario wird ferner ein Effizienzgewinn durch
technische Weiterentwicklung von 25 % angenommen
[ewi](https://www.ewi.uni-koeln.de/cms/wp-content/uploads/2015/12/2014_06_24_ENDBER_P7570_Energiereferenzprognose-GESAMT-FIN-IA.pdf).

Beide separat erstelle Zeitreihen werden anhand der heutigen Marktdurchdringung
gewichtet und in eine mittlere Zeitreihe für Wärmepumpen überführt. Im Jahr
2022 wurden 87 % ASHP und 13 % GSHP abgesetzt nach
[BWP](https://www.waermepumpe.de/fileadmin/user_upload/waermepumpe/05_Presse/01_Pressemitteilungen/BWP_Branchenstudie_2023_DRUCK.pdf),
über die letzten 10 Jahre beträgt das Verhältnis ca. 80:20.
Für 2045 wird daher ein Anteil von 80 % ASHP und 20 % GSHP angenommen.

Verwendet Datensätze:

- [dwd_temperature](../../apipe/store/preprocessed/dwd_temperature/dataset.md)

**Dataset: `datasets/heatpump_cop`**


------------------------------
## Photovoltaik-Freiflächenanlagen

Photovoltaik-Freiflächenanlagen in der Region aus MaStR-Registerdaten als
Geopackage.
Es werden alle Anlagen berücksichtigt, die in Betrieb sind oder sich in
Planung befinden. Anlagen mit Geokoordinaten werden georeferenziert
übernommen, für Anlagen die keine Koordinaten aufweisen (üblicherweise <=30
kW Nennleistung) erfolgt ein Geocoding anhand von PLZ und Ort, um eine
ungefähre Position bereit zu stellen.

Neben einem anlagenscharfen Datensatz wird ein weiterer Datensatz erzeugt,
der alle Anlagen mit approximierter Position je Position zusammenfasst und
jeweils typische Kennwerte enthält (u.a. Anzahl Anlagen, Gesamtleistung).

Jede Anlage wird anhand ihrer Lokation einer Gemeinde (Attribut
`municipality_id`, vgl.
[bkg_vg250_muns_region](../../apipe/store/datasets/bkg_vg250_muns_region/dataset.md)) und
einem Landkreis (Attribut `district_id`, vgl.
[bkg_vg250_muns_region](../../apipe/store/datasets/bkg_vg250_districts_region/dataset.md))
zugeordnet.

Zusätzlich erfolgt eine statistische Auswertung der installierten Leistung in
`bnetza_mastr_pv_ground_stats_muns.csv`.

### Datenkorrektur

Einige Anlagen sind hinsichtlich Ihrer geografischen Lage oder Typs fehlerhaft.
Anhand des Datensatzes
[bnetza_mastr_correction_region](../../apipe/store/raw/bnetza_mastr_correction_region/dataset.md)
wird für diese Anlagen eine Datenkorrektur vorgenommen.

**Dataset: `datasets/bnetza_mastr_pv_ground_region`**


------------------------------
## Landkreise

Landkreise der Region aus Geodaten der Verwaltungsgebiete extrahiert und
nach Landmasse gefiltert (Geofaktor 4 = "mit Struktur Land").

**Dataset: `datasets/bkg_vg250_districts_region`**


------------------------------
## Potenzialgebiete PV-Freiflächen (Region)

### Potenzialflächen

Potenzialgebiete für die Errichtung von PV-Freiflächenanlagen aus dem
[PV- und Windflächenrechner](https://www.agora-energiewende.de/service/pv-und-windflaechenrechner/)
(s. Datensatz [rli_pv_wfr](../../apipe/store/raw/rli_pv_wfr/dataset.md)).

Die Potenzialflächen bilden jene Flächen ab, die für die Nutzung durch
Freiflächen-Photovoltaikanlagen grundsätzlich zur Verfügung stehen. Sie
orientieren sich an der aktuellen Förderkulisse und wurden anhand des
Flächenumfangs sowie den verfügbaren Geodaten ausgewählt: Von den in §37 EEG
2021 definierten Flächen werden Flächen nach §37 Absatz 1 Nummer 2 Buchstaben c,
h und i berücksichtigt (für Details zur Methodik siehe
[methodisches Begleitdokument](https://zenodo.org/record/6794558) zum PV- und
Windflächenrechner).

Dateien:

- Freiflächen-PV auf Acker- und Grünlandflächen mit geringer Bodengüte (Soil
  Quality Rating (SQR) < 40): `potentialarea_pv_agriculture_lfa-off_region.gpkg`
- Potenzialflächen für Freiflächen-PV entlang von Bundesautobahnen und
  Schienenwegen (500m-Streifen): `potentialarea_pv_road_railway_region.gpkg`

### Statistische Auswertung

Die Flächen werden mit den Gemeindegrenzen verschnitten und den Gemeinden
zugeordnet. Je Gemeinde und obigem Flächentyp/Datei wird eine Flächensumme (in
km²) berechnet, siehe `potentialarea_pv_ground_area_stats_muns.csv`. Die
Gemeinden werden über den Schlüssel `municipality_id` (vgl.
[bkg_vg250_muns_region](../../apipe/store/datasets/bkg_vg250_muns_region/dataset.md))
identifiziert.

Des Weiteren werden die Flächenanteile der verfügbaren Potenzialgebiete - deren
Nutzung nur eingeschränkt möglich ist (z.B. durch Naturschutzgebieten etc.) -
gegenüber den gesamten Potenzialgebiete (für die Parametrierung der Regler) nach
`potentialarea_pv_ground_area_shares.json` exportiert.

### Ausbauziele

Es werden PV-Ausbauziele für die Region berechnet, indem die Bundesziele aus den
[BMWK Langfristszenarien](../../apipe/store/preprocessed/bmwk_long_term_scenarios/dataset.md)
i.H.v. 428 GW
([§4 EEG 2023](https://www.gesetze-im-internet.de/eeg_2014/__4.html): 400 GW)
anhand der regional verfügbaren Potenzialflächen disaggregiert werden. Hierzu
wird der Anteil der Flächensumme der beiden o.g. Flächentypen an den bundesweit
verfügbaren Flächen (Datensatz [rli_pv_wfr](../../apipe/store/raw/rli_pv_wfr/dataset.md))
berechnet. Da in den o.g. Ausbauzielen nicht zwischen Freiflächen- und
Aufdach-PV unterschieden wird, wird ein Verhältnis von 50:50 angenommen, d.h.
bundesweit 214 GW auf Freiflächen-PV entfallen.

Es ergeben sich folgende Flächen- und Leistungsanteile:

Gesamt: 0.38 % (819 MW)

- Entlang von BAB und Schienenwegen: 0.13 % (278 MW)
- Acker- und Grünlandflächen mit geringer Bodengüte: 0.25 % (541 MW)

Ergebnisse in `potentialarea_pv_ground_regionalized_targets.json`

### Verwandte Datensätze

- [potentialarea_pv_ground](../../apipe/store/datasets/potentialarea_pv_ground/create.smk))
- [potentialarea_pv_ground_region2](../../apipe/store/datasets/potentialarea_pv_ground_region2/dataset.md))

**Dataset: `datasets/potentialarea_pv_ground_region`**


------------------------------
## Sozialversicherungspflichtig Beschäftigte und Betriebe

Gesamtanzahl sozialversicherungspflichtig Beschäftigte und Betriebsstätten
je Gemeinde für die Region.

Raw datasets:
[ba_employment](../../apipe/store/raw/ba_employment/dataset.md),
[regiostat](../../apipe/store/raw/regiostat/dataset.md)

**Dataset: `datasets/employment_region`**


------------------------------
## Bevölkerungsentwicklung

EinwohnerInnen je Gemeinde: Historische Daten und Prognosen

### Historische Daten bis 2022

Statistisches Bundesamt (Raw dataset:
[destatis_gv](../../apipe/store/raw/destatis_gv/dataset.md))

### Prognosen bis 2035

Statistisches Landesamt Sachsen-Anhalt (Raw dataset:
[stala_st_pop_prog](../../apipe/store/raw/stala_st_pop_prog/dataset.md)). Deaktivieren
mittels entfernen der Zieljahre in [config.yml](../../apipe/store/datasets/population_region/config.yml) im Abschnitt
`prognosis_fstate_munlevel`.

Kann für andere Regionen auch durch DemandRegio (s.u.) ersetzt werden, die
tatsächliche regionale Auflösung wird dadurch reduziert.

### Prognosen bis 2045

DemandRegio (Raw dataset: [demandregio](../../apipe/store/raw/demandregio/dataset.md))
basierend auf der
[14. koordinierten Bevölkerungsvorausberechnung](https://www.destatis.de/DE/Themen/Gesellschaft-Umwelt/Bevoelkerung/Bevoelkerungsvorausberechnung/aktualisierung-bevoelkerungsvorausberechnung.html)
der Statistischen Ämter von Bund und Ländern. Diese Daten liegen auf
Landkreisebene vor, daher erfolgt eine gleichmäßige Skalierung der
dazugehörigen Gemeinden auf den jeweiligen Prognosewert.

Deaktivieren mittels entfernen der Zieljahre in [config.yml](../../apipe/store/datasets/population_region/config.yml) im
Abschnitt `prognosis_germany_districtlevel`.

### Extrapolation

Über 2045 hinaus wird lineare Extrapolation auf Basis der letzten beiden
Prognosejahre unterstützt. Um diese zu aktivieren, müssen lediglich Zieljahre
in die [config.yml](../../apipe/store/datasets/population_region/config.yml) im Abschnitt `extrapolation` eingetragen werden.

**Dataset: `datasets/population_region`**


------------------------------
## Windenergieanlagen

Windenergieanlagen in der Region aus MaStR-Registerdaten als Geopackage.
Es werden alle Anlagen berücksichtigt, die in Betrieb sind oder sich in
Planung befinden. Anlagen mit Geokoordinaten werden georeferenziert
übernommen, für Anlagen die keine Koordinaten aufweisen (üblicherweise <=30
kW Nennleistung) erfolgt ein Geocoding anhand von PLZ und Ort, um eine
ungefähre Position bereit zu stellen.

Neben einem anlagenscharfen Datensatz wird ein weiterer Datensatz erzeugt,
der alle Anlagen mit approximierter Position je Position zusammenfasst und
jeweils typische Kennwerte enthält (u.a. Anzahl Anlagen, Gesamtleistung).

Jede Anlage wird anhand ihrer Lokation einer Gemeinde (Attribut
`municipality_id`, vgl.
[bkg_vg250_muns_region](../../apipe/store/datasets/bkg_vg250_muns_region/dataset.md)) und
einem Landkreis (Attribut `district_id`, vgl.
[bkg_vg250_muns_region](../../apipe/store/datasets/bkg_vg250_districts_region/dataset.md))
zugeordnet.

Zusätzlich erfolgt eine statistische Auswertung der installierten Leistung in
`bnetza_mastr_wind_stats_muns.csv`.

**Dataset: `datasets/bnetza_mastr_wind_region`**


------------------------------
## Bundesländer

Bundesländergrenzen aus Geodaten der Verwaltungsgebiete extrahiert und nach
Landmasse gefiltert (Geofaktor 4 = "mit Struktur Land").

**Dataset: `datasets/bkg_vg250_federal_states`**


------------------------------
## Photovoltaik-Aufdachanlagen

Photovoltaik-Aufdachanlagen in der Region aus MaStR-Registerdaten als
Geopackage.
Es werden alle Anlagen berücksichtigt, die in Betrieb sind oder sich in
Planung befinden. Anlagen mit Geokoordinaten werden georeferenziert
übernommen, für Anlagen die keine Koordinaten aufweisen (üblicherweise <=30
kW Nennleistung) erfolgt ein Geocoding anhand von PLZ und Ort, um eine
ungefähre Position bereit zu stellen.

Neben einem anlagenscharfen Datensatz wird ein weiterer Datensatz erzeugt,
der alle Anlagen mit approximierter Position je Position zusammenfasst und
jeweils typische Kennwerte enthält (u.a. Anzahl Anlagen, Gesamtleistung).

Jede Anlage wird anhand ihrer Lokation einer Gemeinde (Attribut
`municipality_id`, vgl.
[bkg_vg250_muns_region](../../apipe/store/datasets/bkg_vg250_muns_region/dataset.md)) und
einem Landkreis (Attribut `district_id`, vgl.
[bkg_vg250_muns_region](../../apipe/store/datasets/bkg_vg250_districts_region/dataset.md))
zugeordnet.

Zusätzlich erfolgt eine statistische Auswertung der installierten Leistung in
`bnetza_mastr_pv_roof_stats_muns.csv`.

### Datenkorrektur

Einige Anlagen sind hinsichtlich Ihrer geografischen Lage oder Typs fehlerhaft.
Anhand des Datensatzes
[bnetza_mastr_correction_region](../../apipe/store/raw/bnetza_mastr_correction_region/dataset.md)
wird für diese Anlagen eine Datenkorrektur vorgenommen.

**Dataset: `datasets/bnetza_mastr_pv_roof_region`**


------------------------------
## Wasserkraftanlagen

Wasserkraftanlagen in der Region aus MaStR-Registerdaten als Geopackage.
Es werden alle Anlagen berücksichtigt, die in Betrieb sind oder sich in
Planung befinden. Anlagen mit Geokoordinaten werden georeferenziert
übernommen, für Anlagen die keine Koordinaten aufweisen (üblicherweise <=30
kW Nennleistung) erfolgt ein Geocoding anhand von PLZ und Ort, um eine
ungefähre Position bereit zu stellen.

Neben einem anlagenscharfen Datensatz wird ein weiterer Datensatz erzeugt,
der alle Anlagen mit approximierter Position je Position zusammenfasst und
jeweils typische Kennwerte enthält (u.a. Anzahl Anlagen, Gesamtleistung).

Jede Anlage wird anhand ihrer Lokation einer Gemeinde (Attribut
`municipality_id`, vgl.
[bkg_vg250_muns_region](../../apipe/store/datasets/bkg_vg250_muns_region/dataset.md)) und
einem Landkreis (Attribut `district_id`, vgl.
[bkg_vg250_muns_region](../../apipe/store/datasets/bkg_vg250_districts_region/dataset.md))
zugeordnet.

Zusätzlich erfolgt eine statistische Auswertung der installierten Leistung in
`bnetza_mastr_hydro_stats_muns.csv`.

**Dataset: `datasets/bnetza_mastr_hydro_region`**


------------------------------
## Wärmebedarf

Wärmebedarfe (Endenergie) Fernwärme und dezentrale Wärme sowie Wärmezeitreihen
für Haushalte, GHD und Industrie je Gemeinde.

### Gesamtwärmebedarf

Die Berechnung der regionalen Prognosewerte je Verbrauchssektor erfolgt anhand
landesweiter Prognosen aus den
[BMWK Langfristszenarien](../../apipe/store/preprocessed/bmwk_long_term_scenarios/dataset.md).

#### Haushalte

- Jährlicher Wärmebedarf je Gemeinde in MWh: Bundeswert aus
  [AG Energiebilanzen](../../apipe/store/preprocessed/ageb_energy_balance/dataset.md)
  2021 für Raumwärme, Warmwasser und Prozesswärme, desaggregiert auf Gemeinden
  mittels Wärmebedarfs-Rasterdaten aus 2015 (Wärmebedarfsdichte 1ha) aus
  [Peta5](../../apipe/store/raw/seenergies_peta5/dataset.md).
  Anm.: Die Desaggregation könnte alternativ über Zensus "Gebäude mit Wohnraum
  nach Heizungsart" (31231-02-01-5, s.
  [regiostat](../../apipe/store/raw/regiostat/dataset.md) erfolgen)
- Prognosewerte für 2045 werden durch lineare Skalierung mittels Reduktion der
  Gebäudewärmebedarfe aus
  [BMWK Langfristszenarien](../../apipe/store/preprocessed/bmwk_long_term_scenarios/dataset.md)
  berechnet. Hierbei wird das Szenario "TN-Strom" als Grundlage für den Status
  quo verwendet und Werte für 2022 interpoliert. Die Zielwerte werden dem
  Szenario "T45-Strom" entnommen.
- Gemittelte, normierte Gasbedarfszeitreihe (auf 1 MWh) aus
  [DemandRegio](../../apipe/store/preprocessed/demandregio/dataset.md)-Daten von 2022 die
  für alle Zielszenarien und Aggregationsebenen verwendet wird, da die Basis
  SLP-Profile sind und Differenzen zwischen verschiedenen Jahren nur aufgrund
  der Lage von Wochenenden und Feiertagen bestehen. Diese werden daher
  vernachlässigt.

#### GHD

- Jährlicher Wärmebedarf je Gemeinde in MWh: Bundeswert aus
  [AG Energiebilanzen](../../apipe/store/preprocessed/ageb_energy_balance/dataset.md)
  2021 für Raumwärme, Warmwasser und Prozesswärme, desaggregiert auf Gemeinden
  mittels Wärmebedarfs-Rasterdaten aus 2015 (Wärmebedarfsdichte 1ha) aus
  [Peta5](../../apipe/store/raw/seenergies_peta5/dataset.md)
- Prognosewerte für 2045 werden durch lineare Skalierung mittels Reduktion der
  Gebäudewärmebedarfe aus
  [BMWK Langfristszenarien](../../apipe/store/preprocessed/bmwk_long_term_scenarios/dataset.md)
  berechnet. Hierbei wird das Szenario "TN-Strom" als Grundlage für den Status
  quo verwendet und Werte für 2022 interpoliert. Die Zielwerte werden dem
  Szenario "T45-Strom" entnommen.
- Gemittelte, normierte Gasbedarfszeitreihe (auf 1 MWh) aus
  [DemandRegio](../../apipe/store/preprocessed/demandregio/dataset.md)-Daten von 2022 die
  für alle Zielszenarien und Aggregationsebenen verwendet wird, da die Basis
  SLP-Profile sind und Differenzen zwischen verschiedenen Jahren nur aufgrund
  der Lage von Wochenenden und Feiertagen bestehen. Diese werden daher
  vernachlässigt.

#### Industrie

- Jährlicher Wärmebedarf je Gemeinde in MWh: Bundeswert aus
  [AG Energiebilanzen](../../apipe/store/preprocessed/ageb_energy_balance/dataset.md)
  2021 für Raumwärme, Warmwasser und Prozesswärme. Die Desaggregation auf
  Landkreisebene erfolgt anhand des Gesamtenergiebedarfs im verarbeitenden
  Gewerbe aus [Regionalstatistik](../../apipe/store/preprocessed/regiostat/dataset.md).
  Die anschließende Desaggregation auf Gemeindeebene wird mittels
  Beschäftigtenzahlen im verarbeitenden Gewerbe in 2022 aus
  [Regionalstatistik](../../apipe/store/preprocessed/regiostat/dataset.md) vorgenommen.
- Prognosewerte für 2045 werden durch lineare Skalierung mittels Reduktion des
  industriellen Gesamtenergiebedarfs aus
  [BMWK Langfristszenarien](../../apipe/store/preprocessed/bmwk_long_term_scenarios/dataset.md)
  berechnet. Im Unterschied zu Haushalten und GHD liegen die Daten für den
  Wärme- und Stromanteil nicht getrennt vor, sodass auf den
  Gesamtenergiebedarf zurückgegriffen wird.
  Es wird das Szenario "TN-Strom" als Grundlage für den Status quo verwendet und
  Werte für 2022 interpoliert. Die Zielwerte werden dem Szenario "T45-Strom"
  entnommen.
- Gemittelte, normierte Gasbedarfszeitreihe (auf 1 MWh) aus
  [DemandRegio](../../apipe/store/preprocessed/demandregio/dataset.md)-Daten von 2022 die
  für alle Zielszenarien und Aggregationsebenen verwendet wird, da die Basis
  SLP-Profile sind und Differenzen zwischen verschiedenen Jahren nur aufgrund
  der Lage von Wochenenden und Feiertagen bestehen. Diese werden daher
  vernachlässigt.
- Es erfolgt keine Aufteilung des Wärmebedarfs auf unterschiedliche
  Temperaturniveaus.

### Dezentrale Wärme und Fernwärme

Der Gesamtwärmebedarf wird auf dezentrale Heizsysteme und Fernwärme aufgeteilt.
Fernwärmenetze existieren in Dessau-Roßlau, Bitterfeld-Wolfen, Köthen und
Wittenberg.

Da keine Daten zum tatsächlichen Fernwärmebedarf vorliegen, werden Annahmen auf
Basis folgender Quellen getroffen:

- [Zensus 2011: Gebäude nach Heizungsart](https://www.regionalstatistik.de/genesis//online?operation=table&code=31211-04-01-5-B)
- [BMWK Langfristszenarien: Wärmenachfrage in Wärmenetzen (HH&GHD) (2025)](https://enertile-explorer.isi.fraunhofer.de:8443/open-view/54022/62a2667df6f8c176ff129f7ede944837)
- [STALA ST: Wärmebilanz der Industriebetriebe (2021)](https://statistik.sachsen-anhalt.de/themen/wirtschaftsbereiche/energie-und-wasserversorgung/tabellen-energieverwendung#c256237)
- [STALA ST: Energie- und Wasserversorgung](https://statistik.sachsen-anhalt.de/fileadmin/Bibliothek/Landesaemter/StaLa/startseite/Themen/Energie/Berichte/6E403_2020-A.pdf)
- [WindNODE](https://windnode-abw.readthedocs.io/en/latest/energy_system_model.html#district-heating)
- [Peta5: D5 1 District Heating Areas (2020)](https://s-eenergies-open-data-euf.hub.arcgis.com/datasets/b62b8ad79f0e4ae38f032ad6aadb91a0_0/)

Annahmen zu Fernwärmeanteilen (Anteil der Endenergie aus Fernwärme an gesamter
Wärme-Endenergie) je Bedarfssektor:

| Fernwärmenetz     | Haushalte |  GHD | Industrie |
|-------------------|----------:|-----:|----------:|
| Dessau-Roßlau     |      0,36 | 0,36 |      0,19 |
| Bitterfeld-Wolfen |      0,11 | 0,11 |      0,21 |
| Köthen            |      0,07 | 0,07 |      0,21 |
| Wittenberg        |      0,15 | 0,15 |      0,01 |

Die Fernwärmeanteile können in der [config.yml](../../apipe/store/datasets/demand_heat_region/config.yml) im Abschnitt
`district_heating_share` für jeden Sektor separat angepasst werden. Es wird
vereinfachend angenommen, dass der Anteil an Fernwärme für alle
Szenarien/Zieljahre gleich bleibt.

### Beheizungsstruktur

Die Beheizungsstruktur für 2020 und 2045 wird den
[BMWK Langfristszenarien](../../apipe/store/preprocessed/bmwk_long_term_scenarios/dataset.md)
entnommen (Gebäude: Haushalte und GHD Energiebedarf) und für 2022 interpoliert.
Hierbei wird nach Technologien für dezentrale sowie Fernwärme unterschieden.
Für die Biomasse wird der relative Energiebedarf mit Hilfe von Anteilen der
installierten Leistung von spezifischen Biomasse-Konversionsanlagen
[dbfz_biomasss_capacity_rel](../../apipe/store/preprocessed/dbfz_biomass_capacity_rel/dataset.md)
je Technologie aufgelöst. Der Vereinfachung halber wird angenommen, dass die
relative installierte Leistung der relativen Energiemenge entspricht.

### Ergebnisdaten

- Haushalte: Wärmebedarf gesamt: `demand_hh_heat_demand.csv`
- Haushalte: Wärmebedarf Fernwärme: `demand_hh_heat_demand_cen.csv`
- Haushalte: Wärmebedarf dezentrale Wärme: `demand_hh_heat_demand_dec.csv`
- Haushalte: Zeitreihen: `demand_hh_heat_timeseries.csv`

- GHD: Wärmebedarf gesamt: `demand_cts_heat_demand.csv`
- GHD: Wärmebedarf Fernwärme: `demand_cts_heat_demand_cen.csv`
- GHD: Wärmebedarf dezentrale Wärme: `demand_cts_heat_demand_dec.csv`
- GHD: Zeitreihen: `demand_cts_heat_timeseries.csv`

- Industrie: Wärmebedarf gesamt: `demand_ind_heat_demand.csv`
- Industrie: Wärmebedarf Fernwärme: `demand_ind_heat_demand_cen.csv`
- Industrie: Wärmebedarf dezentrale Wärme: `demand_ind_heat_demand_dec.csv`
- GHD: Zeitreihen: `demand_ind_heat_timeseries.csv`

- Beheizungsstruktur dezentral (informativ): `demand_heat_structure_dec.csv`
- Beheizungsstruktur zentral (informativ): `demand_heat_structure_cen.csv`
- Beheizungsstruktur dezentral für Weiterverwendung im Energiesystem:
  `demand_heat_structure_esys_dec.csv`
- Beheizungsstruktur Fernwärme für Weiterverwendung im Energiesystem:
  `demand_heat_structure_esys_cen.csv`

**Dataset: `datasets/demand_heat_region`**


------------------------------
## OpenStreetMap Gebäude

OSM Gebäude aus [osm_filtered](../../apipe/store/preprocessed/osm_filtered/dataset.md)
mittels OGR extrahieren und nach Tags (s. [config.yml](../../apipe/store/datasets/osm_buildings/config.yml)) filtern.

Ziel ist die Ermittlung des regionalen Anteils Gebäudegrundflächen an der
gesamten Gebäudegrundfläche in Deutschland.

Schritte:

- Extraktion aller Gebäude in Deutschland --> `osm_buildings.gpkg`
- Zentroide und Fläche je Gebäude erstellen --> `osm_buildings_centroids.gpkg`
- Mit Region verschneiden --> `osm_buildings_centroids_region.gpkg`
- Flächensumme berechnen --> `osm_buildings_ground_area_region.gpkg`,
  `osm_buildings_ground_area_country.gpkg`
- Regionalen Anteil berechnen --> `osm_buildings_ground_area_share_region.json`

**Achtung:** Konvertierungs- und Extraktionsprozess benötigt ~15 GB
Speicherplatz und kann viel Zeit in Anspruch nehmen.

**Dataset: `datasets/osm_buildings`**


------------------------------
## Verbrennungskraftwerke

Verbrennungskraftwerke in der Region aus MaStR-Registerdaten als Geopackage.
Es werden alle Anlagen berücksichtigt, die in Betrieb sind oder sich in
Planung befinden. Anlagen mit Geokoordinaten werden georeferenziert
übernommen, für Anlagen die keine Koordinaten aufweisen (üblicherweise <=30
kW Nennleistung) erfolgt ein Geocoding anhand von PLZ und Ort, um eine
ungefähre Position bereit zu stellen.

Neben einem anlagenscharfen Datensatz wird ein weiterer Datensatz erzeugt,
der alle Anlagen mit approximierter Position je Position zusammenfasst und
jeweils typische Kennwerte enthält (u.a. Anzahl Anlagen, Gesamtleistung).

Jede Anlage wird anhand ihrer Lokation einer Gemeinde (Attribut
`municipality_id`, vgl.
[bkg_vg250_muns_region](../../apipe/store/datasets/bkg_vg250_muns_region/dataset.md)) und
einem Landkreis (Attribut `district_id`, vgl.
[bkg_vg250_muns_region](../../apipe/store/datasets/bkg_vg250_districts_region/dataset.md))
zugeordnet.

Zusätzlich erfolgt eine statistische Auswertung der installierten Leistung in
`bnetza_mastr_combustion_stats_muns.csv`.

**Dataset: `datasets/bnetza_mastr_combustion_region`**


------------------------------
## Gemeinden

Gemeinden der Region aus Geodaten der Verwaltungsgebiete extrahiert und
nach Landmasse gefiltert (Geofaktor 4 = "mit Struktur Land").

**Dataset: `datasets/bkg_vg250_muns_region`**


------------------------------
## Bezeichner und Namen aus MaStR

Bezeichner und Namen aus MaStR als Mapping <NAME_IN_GEODATEN> ->
<NAME_IN_MASTR> wobei CamelCase aus <NAME_IN_MASTR> in Leerzeichen konvertiert
werden.

**Dataset: `datasets/bnetza_mastr_captions`**


------------------------------
## Regionale Analyse von Freiflächen-PV-Potenzialen

Die Agri-PV-Potenzialflächen für Deutschland aus Datensatz
[potentialarea_pv_ground](../../apipe/store/datasets/potentialarea_pv_ground/dataset.md)
werden hier regionalisiert:

- Zuschnitt der Rasterdaten auf die Region
- Vektorisierung der Einzelflächen aller 3 Nutzungskategorien
  - Entfernung von Gebieten unterhalb einer Mindestflächengröße, einstellbar in
    [config.yml](../../apipe/store/datasets/potentialarea_pv_ground_region2/config.yml) -> `area_threshold`
  - Entfernung von Gebieten unterhalb eines Rasterwertes, einstellbar in
    [config.yml](../../apipe/store/datasets/potentialarea_pv_ground_region2/config.yml) -> `raster_value_threshold`
- Bildung eines Mittelwerts (nutzbare Fläche) je Polygon
- Zuweisung der Gemeinde-ID (`municipality_id`, vgl.
  [bkg_vg250_muns_region](../../apipe/store/apipe/store/datasets/bkg_vg250_muns_region/dataset.md))

### Ergebnisse

#### Geodaten

Die Ergebnisse bieten Einblicke in das regionale PV-Potenzial und umfassen:

- **GeoPackages**: Vektorisierte Darstellung der Agri-PV-Potenzialflächen,
  differenziert nach Bodenqualität und Kulturtyp:
    - `potentialarea_pv_ground_soil_quality_low_region.gpkg` **(A)**
    - `potentialarea_pv_ground_soil_quality_medium_region.gpkg` **(B)**
    - `potentialarea_pv_ground_permanent_crops_region.gpkg` **(C)**

#### Statistische Auswertung

Die Flächen werden mit den Gemeindegrenzen verschnitten und den Gemeinden
zugeordnet. Je Gemeinde und obigem Flächentyp/Datei wird eine Flächensumme (in
km²) berechnet. Die Gemeinden werden über den Schlüssel `municipality_id` (vgl.
[bkg_vg250_muns_region](../../apipe/store/datasets/bkg_vg250_muns_region/dataset.md))
identifiziert.

File: `potentialarea_pv_ground_area_stats_muns.csv` (Einheit: km²)

#### Regionalisierte Ausbauziele

Es werden regionalisierte PV-Ausbauziele für die Region berechnet, indem die
Bundesziele aus den
[BMWK Langfristszenarien](../../apipe/store/preprocessed/bmwk_long_term_scenarios/dataset.md)
i.H.v. 428 GW
([§4 EEG 2023](https://www.gesetze-im-internet.de/eeg_2014/__4.html): 400 GW)
anhand der regional verfügbaren Potenzialflächen disaggregiert werden. Hierzu
wird der Anteil der Flächensumme der drei o.g. Flächentypen an den bundesweit
verfügbaren Flächen (Datensatz [oei_agri_pv](../../apipe/store/raw/oei_agri_pv/dataset.md))
berechnet. Da in den o.g. Ausbauzielen nicht zwischen Freiflächen- und
Aufdach-PV unterschieden wird, wird folgende Aufteilung angenommen (änderbar in
[config.yml](../../apipe/store/datasets/potentialarea_pv_ground_region2/config.yml)), basierend auf dem
[Projektionsbericht 2024]():

TODO: Update ÖI Link Projektionsbericht 2024

- Aufdach-PV: 52 % (221 GW), vgl.
  [potentialarea_pv_roof_region](../../apipe/store/datasets/potentialarea_pv_roof_region/dataset.md)
- Freiflächen-PV (niedrig aufgeständert): 44 % (190 GW)
- Agri-PV (hoch aufgeständert und vertikal bifazial): 4 % (17 GW),
  Die Aufteilung zwischen hoch aufgeständert und vertikal bifazial erfolgt
  flächengewichtet, d.h. ein Flächenverhältnis von 1:9 führt zu 9-facher
  Nutzung von Flächen, auf denen vertikale Anlagen angenommen werden. Durch die
  spezifische Leistungsdichte (Werte s.
  [technology_data](../../apipe/store/raw/technology_data/dataset.md)) können sich andere
  Leistungspotenzial-Verhältnisse ergeben.

Ergebnisfile: `potentialarea_pv_ground_regionalized_targets.json`

- Leistungsziele: `target_power_*` (Einheit: MW)
- Flächenziele: `target_area_*` (Einheit: km²)

### Verwandte Datensätze

- [potentialarea_pv_ground](../../apipe/store/datasets/potentialarea_pv_ground/create.smk))
- [potentialarea_pv_ground_region](../../apipe/store/datasets/potentialarea_pv_ground_region/dataset.md))

**Dataset: `datasets/potentialarea_pv_ground_region2`**


------------------------------
## Geodaten PV- und Windflächenrechner

Geodaten aus dem [PV- und Windflächenrechner](https://www.agora-energiewende.de/service/pv-und-windflaechenrechner/),
extrahiert, zu LAEA Europe (EPSG:3035) umprojiziert und auf die Regionsgrenzen
zugeschnitten.

Preprocessed dataset:
[rli_pv_windflaechenrechner](../../apipe/store/preprocessed/rli_pv_wfr/dataset.md)

**Dataset: `datasets/rli_pv_wfr_region`**


------------------------------
## Region

Region aus Geodaten der Landkreise zusammengeführt.

**Dataset: `datasets/bkg_vg250_region`**


------------------------------
## Settings für App

Einstellungen für die App.

### Layerliste (rechtes Panel)

- Konfiguration: [config.yml](https://github.com/rl-institut/apipe/blob/main/apipe/store/datasets/app_settings/config.yml) --> `map_panel_layer_list`
- Ergebnisfile: `map_panel_layer_list.json`
- Wird manuell in die App eingepflegt (s.
  [map_config.py](https://github.com/rl-institut/digiplan/blob/main/digiplan/map/map_config.py))

### Settings panels

Die im linken Panel aufgeführten Einstellelemente (Slider und Schalter) werden
hier parametriert.

- Konfiguration des Templates:
  [config.yml](https://github.com/rl-institut/apipe/blob/main/apipe/store/datasets/app_settings/config.yml) --> `panel_settings_templates`
- Parametrierung der Slider und Schalter:
  [panels.py](https://github.com/rl-institut/apipe/blob/main/apipe/store/datasets/app_settings/scripts/panels.py)
- Ergebnisfiles:
    - `energy_settings_panel.json`
    - `heat_settings_panel.json`
    - `traffic_settings_panel.json`
- Werden in die App eingelesen

#### Parametrierung der Einstellelemente

Für die Slider werden folgende Attribute gesetzt:
Minimum, Maximum, Schrittweite, Startwert, Status-quo-Wert, Zielwert 2045.
Diese werden wie folgt bestimmt (vgl. auch (i)-Tooltips an den Elementen):

| **Technologie**         | **Element id** | **Maximum**                         | **Startwert**       | **Status-quo-Wert**                             | **Zielwert 2045**                              |
|-------------------------|----------------|-------------------------------------|---------------------|-------------------------------------------------|------------------------------------------------|
| Windenergie             | `s_w_1`        | Inst. Leistung in bestehenden VR/EG | Wie Status-quo-Wert | Inst. Leistung 2022                             | Aus Flächenziel Sachsen-Anhalt (2,2 % in 2032) |
|                         | `s_w_3`        | -                                   | Wie Status-quo-Wert | On                                              | -                                              |
|                         | `s_w_4`        | -                                   | Wie Status-quo-Wert | Off                                             | -                                              |
|                         | `s_w_4_1`      | -                                   | Wie Status-quo-Wert | On                                              | -                                              |
|                         | `s_w_4_2`      | -                                   | Wie Status-quo-Wert | Off                                             | -                                              |
|                         | `s_w_5`        | -                                   | Wie Status-quo-Wert | Off                                             | -                                              |
|                         | `s_w_5_1`      | 100 %                               | Wie Status-quo-Wert | Theoret. Wert berechnet aus inst. Leistung 2022 | -                                              |
|                         | `s_w_5_2`      | 100 %                               | Wie Status-quo-Wert | Theoret. Wert berechnet aus inst. Leistung 2022 | -                                              |
| Freiflächen-PV          | `s_pv_ff_1`    |                                     | Wie Status-quo-Wert | Inst. Leistung 2022                             | Aus EEG 2023 und regionalen Potenzialen        |
|                         | `s_pv_ff_3`    | 100 %                               | Wie Status-quo-Wert | Theoret. Wert berechnet aus inst. Leistung 2022 | -                                              |
|                         | `s_pv_ff_4`    | 100 %                               | Wie Status-quo-Wert | Theoret. Wert berechnet aus inst. Leistung 2022 | -                                              |
| Aufdach-PV              | `s_pv_d_1`     |                                     | Wie Status-quo-Wert | Inst. Leistung 2022                             | Aus EEG 2023 und regionalen Potenzialen        |
|                         | `s_pv_d_3`     | 100 %                               | Wie Status-quo-Wert | Theoret. Wert berechnet aus inst. Leistung 2022 | -                                              |
|                         | `s_pv_d_4`     | 100 %                               | Wie Status-quo-Wert | Aus MaStR                                       | -                                              |
| Wasserkraft             | `s_h_1`        | Inst. Leistung 2022                 | Wie Status-quo-Wert | Inst. Leistung 2022                             | Inst. Leistung 2022                            |
| Stromverbrauch          | `s_v_1`        | 200 %                               | Wie Status-quo-Wert | Verbrauch 2022 (100 %)                          | Wert 2045 aus BMWK Langfristszenarien          |
|                         | `s_v_3`        | 200 %                               | Wie Status-quo-Wert | Verbrauch 2022 (100 %)                          | Wert 2045 aus BMWK Langfristszenarien          |
|                         | `s_v_4`        | 200 %                               | Wie Status-quo-Wert | Verbrauch 2022 (100 %)                          | Wert 2045 aus BMWK Langfristszenarien          |
|                         | `s_v_5`        | 200 %                               | Wie Status-quo-Wert | Verbrauch 2022 (100 %)                          | Wert 2045 aus BMWK Langfristszenarien          |
| Batterie-Großspeicher   | `s_s_g_1`      | 50 %                                | Wie Status-quo-Wert | Aus inst. Kapazität und Einspeisung 2022        | -                                              |
|                         |                |                                     | Wie Status-quo-Wert |                                                 |                                                |
| WP dezentral            | `w_d_wp_1`     | 95 %                                | 50 %                | Inst. Leistung 2022 aus BMWK Langfristszenarien | Wert 2045 aus BMWK Langfristszenarien          |
|                         | `w_d_wp_3`     | 95 %                                | 50 %                | -                                               | -                                              |
|                         | `w_d_wp_4`     | 95 %                                | 50 %                | -                                               | -                                              |
|                         | `w_d_wp_5`     | 95 %                                | 50 %                | -                                               | -                                              |
| WP zentral              | `w_z_wp_1`     | 95 %                                | 50 %                | Inst. Leistung 2022 aus BMWK Langfristszenarien | Wert 2045 aus BMWK Langfristszenarien          |
| Wärmeverbrauch          | `w_v_1`        | 200 %                               | Wie Status-quo-Wert | Verbrauch 2022 (100 %)                          | Wert 2045 aus BMWK Langfristszenarien          |
|                         | `w_v_3`        | 200 %                               | Wie Status-quo-Wert | Verbrauch 2022 (100 %)                          | Wert 2045 aus BMWK Langfristszenarien          |
|                         | `w_v_4`        | 200 %                               | Wie Status-quo-Wert | Verbrauch 2022 (100 %)                          | Wert 2045 aus BMWK Langfristszenarien          |
|                         | `w_v_5`        | 200 %                               | Wie Status-quo-Wert | Verbrauch 2022 (100 %)                          | Wert 2045 aus BMWK Langfristszenarien          |
| Wärmespeicher dezentral | `w_d_s_1`      | 200 %                               | 100 %               | -                                               | -                                              |
| Wärmespeicher zentral   | `w_z_s_1`      | 200 %                               | 100 %               | -                                               | -                                              |

Die Maxima der Regler im Hauptpanel (`s_w_1`, `s_pv_ff_1` usw.) werden in der
App dynamisch aus den durch die UserInnen vorgenommenen Detaileinstellungen
(`s_w_3`, `s_pv_ff_1` usw.) berechnet.

**Dataset: `datasets/app_settings`**


------------------------------
## Captions

Beschriftungen für WebApp.

Dateien:

- Felder: `captions_fields.json`

**Dataset: `datasets/app_captions`**

