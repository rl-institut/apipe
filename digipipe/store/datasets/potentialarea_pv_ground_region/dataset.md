# Potenzialgebiete PV-Freiflächen

## Potenzialflächen

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

## Statistische Auswertung

Die Flächen werden mit den Gemeindegrenzen verschnitten und den Gemeinden
zugeordnet. Je Gemeinde und obigem Flächentyp/Datei wird eine Flächensumme (in
km²) berechnet, siehe `potentialarea_pv_ground_area_stats_muns.csv`. Die
Gemeinden werden über den Schlüssel `municipality_id` (vgl.
[bkg_vg250_muns_region](../../datasets/bkg_vg250_muns_region/dataset.md))
identifiziert.

Des Weiteren werden die Flächenanteile der verfügbaren Potenzialgebiete - deren
Nutzung nur eingeschränkt möglich ist (z.B. durch Naturschutzgebieten etc.) -
gegenüber den gesamten Potenzialgebiete (für die Parametrierung der Regler) nach
`potentialarea_pv_ground_area_shares.json` exportiert.

## Ausbauziele

Es werden PV-Ausbauziele für die Region berechnet, indem die Bundesziele aus den
[BMWK Langfristszenarien](../../preprocessed/bmwk_long_term_scenarios/dataset.md)
i.H.v. 428 GW
([§4 EEG 2023](https://www.gesetze-im-internet.de/eeg_2014/__4.html): 400 GW)
anhand der regional verfügbaren Potenzialflächen disaggregiert werden. Hierzu
wird der Anteil der Flächensumme der beiden o.g. Flächentypen an den bundesweit
verfügbaren Flächen (Datensatz [rli_pv_wfr](../../raw/rli_pv_wfr/dataset.md))
berechnet. Da in den o.g. Ausbauzielen nicht zwischen Freiflächen- und
Aufdach-PV unterschieden wird, wird ein Verhältnis von 50:50 angenommen, d.h.
bundesweit 214 GW auf Freiflächen-PV entfallen.

Es ergeben sich folgende Flächen- und Leistungsanteile:

Gesamt: 0.38 % (819 MW)

- Entlang von BAB und Schienenwegen: 0.13 % (278 MW)
- Acker- und Grünlandflächen mit geringer Bodengüte: 0.25 % (541 MW)

Ergebnisse in `potentialarea_pv_ground_regionalized_targets.json`
