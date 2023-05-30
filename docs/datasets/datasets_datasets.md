# 'Datasets' Datasets 

------------------------------
## Photovoltaik-Aufdachanlagen

Photovoltaik-Aufdachanlagen in der Region aus MaStR-Registerdaten als
Geopackage.
Es werden alle Anlagen berücksichtigt, die in Betrieb sind oder sich in
Planung befinden. Anlagen mit Geokoordinaten werden gereferenziert
übernommen, für Anlagen die keine Koordinaten aufweisen (üblicherweise <=30
kW Nennleistung) erfolgt ein Geocoding anhand von PLZ und Ort, um eine
ungefähre Position bereit zu stellen.

Neben einem anlagenscharfen Datensatz wird ein weiterer Datensatz erzeugt,
der alle Anlagen mit approximierter Position je Position zusammenfasst und
jeweils typische Kennwerte enthält (u.a. Anzahl Anlagen, Gesamtleistung).

**Dataset: `datasets/bnetza_mastr_pv_roof_region`**


------------------------------
## Biomasse-/Biogasanlagen

Biomasse-/Biogasanlagen in der Region aus MaStR-Registerdaten als Geopackage.
Es werden alle Anlagen berücksichtigt, die in Betrieb sind oder sich in
Planung befinden. Anlagen mit Geokoordinaten werden gereferenziert
übernommen, für Anlagen die keine Koordinaten aufweisen (üblicherweise <=30
kW Nennleistung) erfolgt ein Geocoding anhand von PLZ und Ort, um eine
ungefähre Position bereit zu stellen.

Neben einem anlagenscharfen Datensatz wird ein weiterer Datensatz erzeugt,
der alle Anlagen mit approximierter Position je Position zusammenfasst und
jeweils typische Kennwerte enthält (u.a. Anzahl Anlagen, Gesamtleistung).

**Dataset: `datasets/bnetza_mastr_biomass_region`**


------------------------------
## Landkreise

Landkreise der Region aus Geodaten der Verwaltungsgebiete extrahiert und
gefiltert.

**Dataset: `datasets/bkg_vg250_districts_region`**


------------------------------
## Bezeichner und Namen aus MaStR

Bezeichner und Namen aus MaStR als Mapping <NAME_IN_GEODATEN> ->
<NAME_IN_MASTR> wobei CamelCase aus <NAME_IN_MASTR> in Leerzeichen konvertiert
werden.

**Dataset: `datasets/bnetza_mastr_captions`**


------------------------------
## Sozialversicherungspflichtig Beschäftigte und Betriebe

Gesamtanzahl sozialversicherungspflichtig Beschäftigte und Betriebsstätten
je Gemeinde für die Region.

**Dataset: `datasets/employment_region`**


------------------------------
## Region

Region aus Geodaten der Landkreise zusammengeführt.

**Dataset: `datasets/bkg_vg250_region`**


------------------------------
## Strombedarf

Strombedarf für Haushalte, GHD und Industrie auf Gemeindeebene.

Datengrundlage
- Strombedarf 2022: [DemandRegio](../../preprocessed/demandregio/dataset.md)
- Strombedarfsprognosen 2045:
  [BMWK Langfristszenarien](../../preprocessed/bmwk_long_term_scenarios/dataset.md)

Die Berechnung der regionalen Prognosewerte je Verbrauchssektor Haushalte, GHD
und Industrie erfolgt anhand landesweiter Prognosen. Dafür wird der anteilige
Energiebedarf der Region in 2022 am Gesamtbedarf berechnet und dieser unter der
Annahme eines gleichbleibenden regionale Anteils anschließend linear skaliert.
Die Ergebnisse liegen auf NUTS 3-Ebene vor und werden anschließend auf Basis
sektorspezifischer Parameter auf Gemeindeebene desaggregiert (s.u.)

## Haushalte

- Jährlicher Strombedarf je Gemeinde in MWh, von Landkreis- auf Gemeindeebene
  disaggregiert anhand von Bevölkerungsprognosen
  ([STALA ST](../../preprocessed/stala_st_pop_prog/dataset.md)).
- Gemittelte, normierte Bedarfszeitreihe (auf 1 MWh) aus Daten von 2022 die für
  alle Zielszenarien und Aggregationsebenen verwendet wird, da die Basis
  SLP-Profile sind und Differenzen zwischen verschiedenen Jahren nur aufgrund
  der Lage von Wochenenden und Feiertagen bestehen. Diese werden daher
  vernachlässigt.

## GHD

- Jährlicher Strombedarf je Gemeinde in MWh, von Landkreis- auf Gemeindeebene
  disaggregiert anhand von sozialversicherungspflichtig Beschäftigten im Jahr
  2022 ([BA für Arbeit](../../preprocessed/ba_employment/dataset.md)).
- Gemittelte, normierte Bedarfszeitreihe (auf 1 MWh) aus Daten von 2022 die für
  alle Zielszenarien und Aggregationsebenen verwendet wird. Basis bilden sowohl
  SLP- als auch branchenspezifische Profile. Aufgrund der geringen Differenzen
  zwischen den Landkreisen werden diese gemittelt. Differenzen zwischen
  verschiedenen Jahren bestehen nur aufgrund der Lage von Wochenenden und
  Feiertagen und werden daher vernachlässigt.

## Industrie

- Jährlicher Strombedarf je Gemeinde in MWh, von Landkreis- auf Gemeindeebene
  disaggregiert anhand der sozialversicherungspflichtig Beschäftigten im Jahr
  2022.
- Gemittelte, normierte Bedarfszeitreihe (auf 1 MWh) aus Daten von 2022 die für
  alle Zielszenarien und Aggregationsebenen verwendet wird. Basis bilden sowohl
  SLP- als auch branchenspezifische Profile. Aufgrund der geringen Differenzen
  zwischen den Landkreisen werden diese gemittelt. Differenzen zwischen
  verschiedenen Jahren bestehen nur aufgrund der Lage von Wochenenden und
  Feiertagen und werden daher vernachlässigt.

**Dataset: `datasets/demand_electricity_region`**


------------------------------
## Windenergieanlagen

Windenergieanlagen in der Region aus MaStR-Registerdaten als Geopackage.
Es werden alle Anlagen berücksichtigt, die in Betrieb sind oder sich in
Planung befinden. Anlagen mit Geokoordinaten werden gereferenziert
übernommen, für Anlagen die keine Koordinaten aufweisen (üblicherweise <=30
kW Nennleistung) erfolgt ein Geocoding anhand von PLZ und Ort, um eine
ungefähre Position bereit zu stellen.

Neben einem anlagenscharfen Datensatz wird ein weiterer Datensatz erzeugt,
der alle Anlagen mit approximierter Position je Position zusammenfasst und
jeweils typische Kennwerte enthält (u.a. Anzahl Anlagen, Gesamtleistung).

**Dataset: `datasets/bnetza_mastr_wind_region`**


------------------------------
## Verbrennungskraftwerke

Verbrennungskraftwerke in der Region aus MaStR-Registerdaten als Geopackage.
Es werden alle Anlagen berücksichtigt, die in Betrieb sind oder sich in
Planung befinden. Anlagen mit Geokoordinaten werden gereferenziert
übernommen, für Anlagen die keine Koordinaten aufweisen (üblicherweise <=30
kW Nennleistung) erfolgt ein Geocoding anhand von PLZ und Ort, um eine
ungefähre Position bereit zu stellen.

Neben einem anlagenscharfen Datensatz wird ein weiterer Datensatz erzeugt,
der alle Anlagen mit approximierter Position je Position zusammenfasst und
jeweils typische Kennwerte enthält (u.a. Anzahl Anlagen, Gesamtleistung).

**Dataset: `datasets/bnetza_mastr_combustion_region`**


------------------------------
## OpenStreetMap - Wälder

Waldflächen aus OpenStreetMap, Daten extrahiert anhand von spezifischen Tags.

**Dataset: `datasets/osm_forest`**


------------------------------
## Photovoltaik-Freiflächenanlagen

Photovoltaik-Freiflächenanlagen in der Region aus MaStR-Registerdaten als
Geopackage.
Es werden alle Anlagen berücksichtigt, die in Betrieb sind oder sich in
Planung befinden. Anlagen mit Geokoordinaten werden gereferenziert
übernommen, für Anlagen die keine Koordinaten aufweisen (üblicherweise <=30
kW Nennleistung) erfolgt ein Geocoding anhand von PLZ und Ort, um eine
ungefähre Position bereit zu stellen.

Neben einem anlagenscharfen Datensatz wird ein weiterer Datensatz erzeugt,
der alle Anlagen mit approximierter Position je Position zusammenfasst und
jeweils typische Kennwerte enthält (u.a. Anzahl Anlagen, Gesamtleistung).

**Dataset: `datasets/bnetza_mastr_pv_ground_region`**


------------------------------
## Wasserkraftanlagen

Wasserkraftanlagen in der Region aus MaStR-Registerdaten als Geopackage.
Es werden alle Anlagen berücksichtigt, die in Betrieb sind oder sich in
Planung befinden. Anlagen mit Geokoordinaten werden gereferenziert
übernommen, für Anlagen die keine Koordinaten aufweisen (üblicherweise <=30
kW Nennleistung) erfolgt ein Geocoding anhand von PLZ und Ort, um eine
ungefähre Position bereit zu stellen.

Neben einem anlagenscharfen Datensatz wird ein weiterer Datensatz erzeugt,
der alle Anlagen mit approximierter Position je Position zusammenfasst und
jeweils typische Kennwerte enthält (u.a. Anzahl Anlagen, Gesamtleistung).

**Dataset: `datasets/bnetza_mastr_hydro_region`**


------------------------------
## Geo- oder Solarthermie-, Grubengas- und Klärschlamm-Anlagen

Anlagen der Geo- oder Solarthermie, Grubengas und Klärschlamm in der Region
aus MaStR-Registerdaten als Geopackage.
Es werden alle Anlagen berücksichtigt, die in Betrieb sind oder sich in
Planung befinden. Anlagen mit Geokoordinaten werden gereferenziert
übernommen, für Anlagen die keine Koordinaten aufweisen (üblicherweise <=30
kW Nennleistung) erfolgt ein Geocoding anhand von PLZ und Ort, um eine
ungefähre Position bereit zu stellen.

Neben einem anlagenscharfen Datensatz wird ein weiterer Datensatz erzeugt,
der alle Anlagen mit approximierter Position je Position zusammenfasst und
jeweils typische Kennwerte enthält (u.a. Anzahl Anlagen, Gesamtleistung).

**Dataset: `datasets/bnetza_mastr_gsgk_region`**


------------------------------
## Bevölkerungsentwicklung

EinwohnerInnen je Gemeinde: Historische Daten und Prognosen

## Historische Daten bis 2022

Statistisches Bundesamt (Raw dataset:
[destatis_gv](../../raw/destatis_gv/dataset.md))

## Prognosen bis 2035

Statistisches Landesamt Sachsen-Anhalt (Raw dataset:
[stala_st_pop_prog](../../raw/stala_st_pop_prog/dataset.md)). Deaktivieren
mittels entfernen der Zieljahre in `config.yml` im Abschnitt
`prognosis_fstate_munlevel`.

Kann für andere Regionen auch durch DemandRegio (s.u.) ersetzt werden, die
tatsächliche regionale Auflösung wird dadurch reduziert.

## Prognosen bis 2045

DemandRegio (Raw dataset: [demandregio](../../raw/demandregio/dataset.md))
basierend auf der
[14. koordinierten Bevölkerungsvorausberechnung](https://www.destatis.de/DE/Themen/Gesellschaft-Umwelt/Bevoelkerung/Bevoelkerungsvorausberechnung/aktualisierung-bevoelkerungsvorausberechnung.html)
der Statistischen Ämter von Bund und Ländern. Diese Daten liegen auf
Landkreisebene vor, daher erfolgt eine gleichmäßige Skalierung der
dazugehörigen Gemeinden auf den jeweiligen Prognosewert.

Deaktivieren mittels entfernen der Zieljahre in `config.yml` im Abschnitt
`prognosis_germany_districtlevel`.

## Extrapolation

Über 2045 hinaus wird lineare Extrapolation auf Basis der letzten beiden
Prognosejahre unterstützt. Um diese zu aktivieren, müssen lediglich Zieljahre
in die `config.yml` im Abschnitt `extrapolation` eingetragen werden.

**Dataset: `datasets/population_region`**


------------------------------
## Speicheranlagen

Speicheranlagen in der Region aus MaStR-Registerdaten als Geopackage.
Es werden alle Anlagen berücksichtigt, die in Betrieb sind oder sich in
Planung befinden. Anlagen mit Geokoordinaten werden gereferenziert
übernommen, für Anlagen die keine Koordinaten aufweisen (üblicherweise <=30
kW Nennleistung) erfolgt ein Geocoding anhand von PLZ und Ort, um eine
ungefähre Position bereit zu stellen.

Neben einem anlagenscharfen Datensatz wird ein weiterer Datensatz erzeugt,
der alle Anlagen mit approximierter Position je Position zusammenfasst und
jeweils typische Kennwerte enthält (u.a. Anzahl Anlagen, Gesamtleistung).

**Dataset: `datasets/bnetza_mastr_storage_region`**


------------------------------
## Gemeinden

Gemeinden der Region aus Geodaten der Verwaltungsgebiete extrahiert und
gefiltert.

**Dataset: `datasets/bkg_vg250_muns_region`**

