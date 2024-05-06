# Dachflächenpotenzial PV-Aufdachanlagen in der Region

Berechnung der installierten Leistung und des Ertrags von PV-Aufdachanlagen in
der Region aus Datensatz
[wfbb_pv_roof_potential](../../preprocessed/wfbb_pv_roof_potential/dataset.md).

Es werden nur Dächer verwendet, deren Eignung über 60 % beträgt, d.h. geeignet
oder gut geeignet sind (Klassifikation s.
[wfbb_pv_roof_potential](../../preprocessed/wfbb_pv_roof_potential/dataset.md)).
Der Grenzwert `roof_suitability_threshold` ist in [config.yml](config.yml)
änderbar.

## Ergebnisse

### Geodaten

- `potentialarea_pv_roof_region.gpkg`

### Statistische Auswertung

Die Gebäudezentroide werden mit den Gemeindegrenzen verschnitten und den
Gemeinden zugeordnet. Je Gemeinde und obigem Flächentyp/Datei wird eine
Flächensumme (in km²) berechnet und in
`potentialarea_pv_roof_area_stats_muns.csv` geschrieben.

Des Weiteren wird je Gemeinde der relative Anteil der bereits installierten
Anlagenleistung an der theoretisch installierbaren Leistung (bei
100% Dachnutzung) berechnet.
Ergebnisdaten: `potentialarea_pv_roof_deployment_stats_muns.csv`

Die Gemeinden werden über den Schlüssel `municipality_id` (vgl.
[bkg_vg250_muns_region](../../datasets/bkg_vg250_muns_region/dataset.md))
identifiziert.

### Regionalisierte Ausbauziele

Es werden anhand überregionaler Ziele und Szenarien PV-Ausbauziele für die
Region berechnet:

File: `potentialarea_pv_roof_regionalized_targets.json`

Da in den Ausbauzielen nicht zwischen Freiflächen- und Aufdach-PV unterschieden
wird, wird folgende Aufteilung angenommen (Parameter`pv_roof_share` in
[config.yml](config.yml)), basierend auf dem
[Projektionsbericht 2024](https://todo):

TODO: Update ÖI Link Projektionsbericht 2024

- Aufdach-PV: 52 %
- Freiflächen-PV (niedrig aufgeständert): 44 %, vgl.
  [potentialarea_pv_ground_region2](../../datasets/potentialarea_pv_ground_region2/dataset.md)
- Agri-PV (hoch aufgeständert und vertikal bifazial): 4 %

### Aus BMWK Langfristszenarien

Bundesziele aus den
[BMWK Langfristszenarien](../../preprocessed/bmwk_long_term_scenarios/dataset.md)
i.H.v. 428 GW
([§4 EEG 2023](https://www.gesetze-im-internet.de/eeg_2014/__4.html): 400 GW)
werden anhand der Gebäudegrundflächen disaggregiert. Hierzu wird der Anteil der
Gebäudegrundflächen in der Region an der bundesweiten Gebäudegrundflächen
berechnet (s. Datensatz [osm_buildings](../osm_buildings/dataset.md)) und die
Ziele linear skaliert.

Key: `bmwk_de`

### Aus Energiestrategie Brandenburg 2040

Die Brandenburger Ziele für 2030 und 2040 (vgl. Datensatz
[mwae_bb_energy_strategy_region](../../datasets/mwae_bb_energy_strategy_region/dataset.md))
werden anhand der Regionsfläche (15,48 %) linear skaliert.

Key: `mwae_bb`
