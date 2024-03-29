# Regionale Analyse von Freiflächen-PV-Potenzialen

Die Agri-PV-Potenzialflächen für Deutschland aus Datensatz
[potentialarea_pv_ground](../../datasets/potentialarea_pv_ground/dataset.md)
werden hier regionalisiert:

- Zuschnitt der Rasterdaten auf die Region
- Vektorisierung der Einzelflächen aller 3 Nutzungskategorien
  - Entfernung von Gebieten unterhalb einer Mindestflächengröße, einstellbar in
    [config.yml](config.yml) -> `area_threshold`
  - Entfernung von Gebieten unterhalb eines Rasterwertes, einstellbar in
    [config.yml](config.yml) -> `raster_value_threshold`
- Bildung eines Mittelwerts (nutzbare Fläche) je Polygon
- Zuweisung der Gemeinde-ID (`municipality_id`, vgl.
  [bkg_vg250_muns_region](../../apipe/store/datasets/bkg_vg250_muns_region/dataset.md))

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
km²) berechnet und in `potentialarea_pv_ground_area_stats_muns.csv` geschrieben.
Die Gemeinden werden über den Schlüssel `municipality_id` (vgl.
[bkg_vg250_muns_region](../../datasets/bkg_vg250_muns_region/dataset.md))
identifiziert.

Des Weiteren werden die Flächenanteile (in %) gegenüber der gesamten
Regionsfläche (für die Parametrierung der Regler) berechnet und nach
`potentialarea_pv_ground_area_shares.json` exportiert.

### Regionalisierte Ausbauziele

Es werden anhand überregionaler Ziele und Szenarien PV-Ausbauziele für die
Region berechnet:

File: `potentialarea_pv_ground_regionalized_targets.json`

- Leistungsziele: `target_power_*` (Einheit: MW)
- Flächenziele: `target_area_*` (Einheit: km²)

Da in den Ausbauzielen nicht zwischen Freiflächen- und Aufdach-PV unterschieden
wird, wird folgende Aufteilung angenommen, änderbar in [config.yml](config.yml),
basierend auf dem
[Projektionsbericht 2024](https://todo):

TODO: Update ÖI Link Projektionsbericht 2024

- Aufdach-PV: 52 %
- Freiflächen-PV (niedrig aufgeständert): 44 %, vgl.
  [potentialarea_pv_roof_region2](../../datasets/potentialarea_pv_roof_region2/dataset.md)
- Agri-PV (hoch aufgeständert und vertikal bifazial): 4 %
  Die Aufteilung zwischen hoch aufgeständert und vertikal bifazial erfolgt
  flächengewichtet, d.h. ein Flächenverhältnis von 1:9 führt zu 9-facher
  Nutzung von Flächen, auf denen vertikale Anlagen angenommen werden. Durch die
  spezifische Leistungsdichte (Werte s.
  [technology_data](../../raw/technology_data/dataset.md)) können sich andere
  Leistungspotenzial-Verhältnisse ergeben.

### Aus BMWK Langfristszenarien

Bundesziele aus den
[BMWK Langfristszenarien](../../preprocessed/bmwk_long_term_scenarios/dataset.md)
i.H.v. 428 GW
([§4 EEG 2023](https://www.gesetze-im-internet.de/eeg_2014/__4.html): 400 GW)
werden anhand der regional verfügbaren Potenzialflächen disaggregiert. Hierzu
wird der Anteil der Flächensumme der drei o.g. Flächentypen an den bundesweit
verfügbaren Flächen (Datensatz [oei_agri_pv](../../raw/oei_agri_pv/dataset.md))
berechnet und die Ziele linear skaliert.

Key: `bmwk_de`

### Aus Energiestrategie Brandenburg 2040

Die Brandenburger Ziele für 2030 und 2040 (vgl. Datensatz
[mwae_bb_energy_strategy_region](../../datasets/mwae_bb_energy_strategy_region/dataset.md))
werden anhand der Regionsfläche (15,48 %) linear skaliert.

Key: `mwae_bb`

## Verwandte Datensätze

- [potentialarea_pv_ground](../../datasets/potentialarea_pv_ground/create.smk))
- [potentialarea_pv_ground_region](../../datasets/potentialarea_pv_ground_region/dataset.md))
