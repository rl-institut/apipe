# Dachflächenpotenzial PV-Aufdachanlagen in der Region

Berechnung der installierten Leistung und des Ertrags von PV-Aufdachanlagen in
der Region aus Datensatz
[wfbb_pv_roof_potential](../../preprocessed/wfbb_pv_roof_potential/dataset.md).

Es werden nur Dächer verwendet, deren Eignung über 60 % beträgt, d.h. geeignet
oder gut geeignet sind (Klassifikation s.
[wfbb_pv_roof_potential](../../preprocessed/wfbb_pv_roof_potential/dataset.md)).
Der Grenzwert `roof_suitability_threshold` ist in [config.yml](config.yml)
änderbar.

Es werden Statistiken je Gemeinde erstellt, hierfür werden die Gebäudezentroide
mit den Gemeindegrenzen verschnitten und den Gemeinden zugeordnet.
Ergebnisdaten: `potentialarea_pv_roof_area_stats_muns.csv`

Des Weiteren wird je Gemeinde der relative Anteil der bereits installierten
Anlagenleistung an der theoretisch installierbaren Leistung (bei
100% Dachnutzung) berechnet.
Ergebnisdaten: `potentialarea_pv_roof_deployment_stats_muns.csv`

Die Gemeinden werden über den Schlüssel `municipality_id` (vgl.
[bkg_vg250_muns_region](../../datasets/bkg_vg250_muns_region/dataset.md))
identifiziert.

## Ausbauziele

Es werden PV-Ausbauziele für die Region berechnet, indem die Bundesziele aus den
[BMWK Langfristszenarien](../../preprocessed/bmwk_long_term_scenarios/dataset.md)
i.H.v. 428 GW
([§4 EEG 2023](https://www.gesetze-im-internet.de/eeg_2014/__4.html): 400 GW)
anhand der Gebäudegrundflächen disaggregiert werden. Hierzu wird der Anteil der
Gebäudegrundflächen in der Region an der bundesweiten Gebäudegrundflächen
berechnet (s. Datensatz [osm_buildings](../osm_buildings/dataset.md)) und die
Ziele linear skaliert. Da in den o.g. Ausbauzielen nicht zwischen Freiflächen-
und Aufdach-PV unterschieden wird, wird folgende Aufteilung angenommen
(Parameter`pv_roof_share`, änderbar in [config.yml](config.yml)):

- Aufdach-PV: 52 % (221 GW)
- Freiflächen-PV (niedrig aufgeständert): 44 % (190 GW), vgl.
  [potentialarea_pv_ground_region2](../../datasets/potentialarea_pv_ground_region2/dataset.md)
- Agri-PV (hoch aufgeständert und vertikal bifazial): 4 % (17 GW)

File: `potentialarea_pv_roof_regionalized_targets.json`
