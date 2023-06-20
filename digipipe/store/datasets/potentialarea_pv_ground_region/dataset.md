# Potenzialgebiete PV-Freiflächen

Potenzialgebiete für die Errichtung von PV-Freiflächenanlagen aus dem
[PV- und Windflächenrechner](https://www.agora-energiewende.de/service/pv-und-windflaechenrechner/)
(s. Datensatz [rli_pv_wfr](../../raw/rli_pv_wfr/dataset.md)).

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

Die darin verwendeten Attributtexte werden in die Datei
`potentialarea_pv_ground_attribute_captions.json` exportiert.

Die Flächen werden mit den Gemeindegrenzen verschnitten und den Gemeinden
zugeordnet. Je Gemeinde und obigem Flächentyp/Datei wird eine Flächensumme (in
ha) berechnet, siehe `potentialarea_pv_ground_area_stats_muns.json`. Die
Gemeinden werden über den Schlüssel `municipality_id` (vgl.
[bkg_vg250_muns_region](../../datasets/bkg_vg250_muns_region/dataset.md))
identifiziert.

Des Weiteren werden die Flächenanteile der verfügbaren Potenzialgebiete - die
zusätzlichen Einschränkungen wie Naturschutzgebieten etc. unterworfen sind -
gegenüber den gesamten Potenzialgebiete (für die Parametrierung der Regler) nach
`potentialarea_pv_ground_area_shares.json` exportiert.
