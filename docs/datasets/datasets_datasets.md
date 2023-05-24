# 'Datasets' Datasets 

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
## OpenStreetMap - Wälder

Waldflächen aus OpenStreetMap, Daten extrahiert anhand von spezifischen Tags.

**Dataset: `datasets/osm_forest`**


------------------------------
## Gemeinden

Gemeinden der Region aus Geodaten der Verwaltungsgebiete extrahiert und
gefiltert.

**Dataset: `datasets/bkg_vg250_muns_region`**


------------------------------
## Bezeichner und Namen aus MaStR

Bezeichner und Namen aus MaStR als Mapping <NAME_IN_GEODATEN> ->
<NAME_IN_MASTR> wobei CamelCase aus <NAME_IN_MASTR> in Leerzeichen konvertiert
werden.

**Dataset: `datasets/bnetza_mastr_captions`**


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

EinwohnerInnen je Gemeinde: Historische Daten, Prognosen und darüber hinaus
linear extrapolierte Werte.

**Dataset: `datasets/population`**


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
## Landkreise

Landkreise der Region aus Geodaten der Verwaltungsgebiete extrahiert und
gefiltert.

**Dataset: `datasets/bkg_vg250_districts_region`**


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
## Region

Region aus Geodaten der Landkreise zusammengeführt.

**Dataset: `datasets/bkg_vg250_region`**

