# Regionale Analyse von Freiflächen-PV-Potenzialen

Die Agri-PV-Potenzialflächen für Deutschland aus Datensatz
[potentialarea_pv_ground](../../datasets/potentialarea_pv_ground/dataset.md)
werden hier regionalisiert:

- Zuschnitt der Rasterdaten auf die Region
- Vektorisierung der Einzelflächen aller 3 Nutzungskategorien
- Entfernung von Gebieten unterhalb einer Mindestflächengröße
- Bildung eines Mittelwerts (nutzbare Fläche) je Polygon
- Zuweisung der Gemeinde-ID (`municipality_id`, vgl.
  [bkg_vg250_muns_region](../../apipe/store/datasets/bkg_vg250_muns_region/dataset.md))

Die Mindestflächengröße für die Vektorisierung kann mittels `area_threshold` in
der [config.yml](config.yml) festgelegt werden. Flächen kleiner als dieser
Grenzwert werden entfernt.

## Ergebnisse

### Geodaten

Die Ergebnisse bieten Einblicke in das regionale PV-Potenzial und umfassen:

- **GeoPackages**: Vektorisierte Darstellung der Agri-PV-Potenzialflächen,
  differenziert nach Bodenqualität und Kulturtyp:
    - `potentialarea_pv_ground_soil_quality_low_region.gpkg` **(A)**
    - `potentialarea_pv_ground_soil_quality_medium_region.gpkg` **(B)**
    - `potentialarea_pv_ground_permanent_crops_region.gpkg` **(C)**

### Statistische Auswertung

Die Flächen werden mit den Gemeindegrenzen verschnitten und den Gemeinden
zugeordnet. Je Gemeinde und obigem Flächentyp/Datei wird eine Flächensumme (in
km²) berechnet. Die Gemeinden werden über den Schlüssel `municipality_id` (vgl.
[bkg_vg250_muns_region](../../datasets/bkg_vg250_muns_region/dataset.md))
identifiziert.

File: `potentialarea_pv_ground_area_stats_muns.csv` (Einheit: km²)

### Regionalisierte Ausbauziele

Es werden regionalisierte PV-Ausbauziele für die Region berechnet, indem die
Bundesziele aus den
[BMWK Langfristszenarien](../../preprocessed/bmwk_long_term_scenarios/dataset.md)
i.H.v. 428 GW
([§4 EEG 2023](https://www.gesetze-im-internet.de/eeg_2014/__4.html): 400 GW)
anhand der regional verfügbaren Potenzialflächen disaggregiert werden. Hierzu
wird der Anteil der Flächensumme der drei o.g. Flächentypen an den bundesweit
verfügbaren Flächen (Datensatz [oei_agri_pv](../../raw/oei_agri_pv/dataset.md))
berechnet. Da in den o.g. Ausbauzielen nicht zwischen Freiflächen- und
Aufdach-PV unterschieden wird, wird folgende Aufteilung angenommen (änderbar in
[config.yml](config.yml)), basierend auf dem
[Projektionsbericht 2024]():

TODO: Update ÖI Link Projektionsbericht 2024

- Aufdach-PV: 52 % (221 GW), vgl.
  [potentialarea_pv_roof_region](../../datasets/potentialarea_pv_roof_region/dataset.md)
- Freiflächen-PV (niedrig aufgeständert): 44 % (190 GW)
- Agri-PV (hoch aufgeständert und vertikal bifazial): 4 % (17 GW),
  Die Aufteilung zwischen hoch aufgeständert und vertikal bifazial erfolgt
  flächengewichtet, d.h. ein Flächenverhältnis von 1:9 führt zu 9-facher
  Nutzung von Flächen, auf denen vertikale Anlagen angenommen werden. Durch die
  spezifische Leistungsdichte (Werte s.
  [technology_data](../../raw/technology_data/dataset.md)) können sich andere
  Leistungspotenzial-Verhältnisse ergeben.

Ergebnisfile: `potentialarea_pv_ground_regionalized_targets.json`

- Leistungsziele: `target_power_*` (Einheit: MW)
- Flächenziele: `target_area_*` (Einheit: km²)
